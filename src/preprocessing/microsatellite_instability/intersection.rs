//! intersection.rs
//!
//! Streaming intersection alongside with dummy indel
//! injection and variant annotation for MSI detection.
//!
//! This module provides:
//! 1. Streaming intersection of BED regions with VCF variants;
//! 2. Perfect microsatellite repeat detection;
//! 3. Variant filtering;
//! 4. Dummy indel injection for MS regions without variants;
//! 5. MS region annotation (INFO/MS_REGION).
//!
//! The streaming approach processes BED regions sequentially while maintaining
//! a sliding window of VCF variants, enabling memory-efficient analysis of
//! large datasets.
//!
//! Note: This module assumes that both the BED and VCF files are sorted by chromosome
//! (lexicographically) and position, which is a common requirement for genomic analyses.

use std::collections::VecDeque;
use std::path::Path;

use anyhow::{Context, Result};
use bio::io::bed;
use log::debug;
use rust_htslib::bcf::{self, Read};

use crate::errors::Error;
use crate::utils::bcf_utils::get_chrom;
use crate::utils::genomics::calculate_indel_position;
use crate::utils::ms_bed::{parse_bed_record, BedRegion};

use super::variant_analysis::should_include_variant;
use super::writer::{inject_dummy_deletion, write_variant, VariantInWindow};

/* ============ Data Structures =================== */

/// Statistics from MSI preprocessing.
#[derive(Debug)]
pub(super) struct PreprocessingStats {
    /// Total number of BED regions processed (including invalid motifs).
    pub total_regions: usize,
    /// Number of valid BED regions with motif length 1-6 bp.
    pub valid_regions: usize,
    /// Number of indels that were annotated with MS region information.
    pub annotated_indels: usize,
    /// Number of dummy indels injected for MS regions without variants.
    pub dummy_indels: usize,
}

/* ================================================ */

/* ============ Streaming Intersection ============ */

/// Check if a variant's indel overlaps with a BED region.
///
/// Accounts for VCF anchor bases by calculating where the indel actually
/// occurs (after the anchor), then checking if that position falls within
/// the pure repeat region. Also filters out complex variants and SNVs.
///
/// # Algorithm
/// 1. Get VCF position (0-based, points to anchor base)
/// 2. Extract REF and ALT alleles for this specific alt_idx
/// 3. Calculate actual indel position using genomics::calculate_indel_position
///    - Returns None if: complex variant (both tails non-empty), SNV, or identical sequences
///    - Returns Some(pos) if clean indel
/// 4. Check if indel position is within BED region [start, end)
///
/// # Note:
/// We perform point-based overlap, which means we only consider the variant position.
///
/// # Arguments
/// * `record` - BCF record representing the variant
/// * `region` - BED region to check overlap against
/// * `alt_idx` - Index of the alternate allele to analyze (0-based into ALT array)
///
/// # Returns
/// * `true` if clean indel position is within region
/// * `false` if complex variant, SNV, or outside region
///
/// # Example
/// VCF: POS=18630802 (0-based), REF=GCCT, ALT=G
/// Anchor = G (1 base)
/// Indel position = 18630802 + 1 = 18630803
/// Region: [18630803, 18630833) (pure CCT repeat)
/// Result: 18630803 >= 18630803 && 18630803 < 18630833 , TRUE
#[inline]
fn variant_overlaps_region(record: &bcf::Record, region: &BedRegion, alt_idx: usize) -> bool {
    let vcf_pos = record.pos() as u64;
    let alleles = record.alleles();
    let ref_allele = alleles[0];
    let alt_allele = alleles[alt_idx + 1]; // +1 because alleles[0] is REF

    // Calculate where the indel occurs (None if complex variant/SNV)
    match calculate_indel_position(vcf_pos, ref_allele, alt_allele) {
        Some(indel_pos) => indel_pos >= region.start && indel_pos < region.end,
        None => false,
    }
}

/// Streams through sorted VCF variants and BED regions to identify perfect microsatellite
/// indels. Annotates matching variants with `INFO/REGION_ID` field containing
/// comma-separated region identifiers for all overlapping regions. Injects dummy
/// deletions for regions with no observed perfect indels.
///
/// # Algorithm
/// 1. Pre-scan BED file to collect all chromosomes, add missing contigs to output VCF header
/// 2. Stream through BED regions sequentially (bounded memory usage)
/// 3. Maintain sliding window of VCF variants around current region
/// 4. For each BED region:
///    - Remove and write variants before region (past window boundary)
///    - Load variants until past region end
///    - Accumulate REGION_ID annotations for overlapping perfect MS indels
///    - Inject dummy deletion if no perfect indel found
/// 5. Write all remaining variants after BED exhausted
///
/// # Output Behavior
/// - All input variants are written to output (preserves complete VCF)
/// - Perfect MS indels receive `INFO/REGION_ID` annotation with overlapping region IDs
/// - Non-MS variants are written unchanged (no annotation)
/// - Dummy indels are injected for regions without observed perfect indels
///
/// # Perfect MS Indel Criteria
/// - Must be a clean indel (insertion or deletion, not complex)
/// - Changed sequence must be exact tandem repeats of the region's motif
/// - Must overlap the BED region's coordinates (anchor-aware positioning)
///
/// # Arguments
/// * `vcf_path` - Path to candidate VCF file (must contain indels)
/// * `bed_path` - Path to microsatellite regions BED file (format: chrom, start, end, NxMOTIF)
/// * `output` - Output path (None = stdout)
///
/// # Returns
/// * `Ok(PreprocessingStats)` on success
/// * `Err` if fails
///
/// # Example
/// assert!(process_and_annotate(Path::new("candidates.vcf"), Path::new("ms_regions.bed"), Some(Path::new("output.vcf"))).is_ok());
pub(super) fn process_and_annotate(
    input_vcf: &mut bcf::Reader,
    bed_path: &Path,
    mut writer: &mut bcf::Writer,
) -> Result<PreprocessingStats> {
    /* ========== Setup ========== */
    let header_view = input_vcf.header().clone();

    /* ===== Main processing loop ===== */
    let mut bed_reader = bed::Reader::from_file(bed_path).context("Failed to open BED file")?;

    let mut total_regions = 0;
    let mut skipped_invalid_regions = 0;
    let mut total_annotated_indels = 0;
    let mut total_dummy_indels = 0;
    let mut variant_window: VecDeque<VariantInWindow> = VecDeque::new();
    let mut seen_any_chrom_overlap = false;

    /* ========== Main Loop: Process each BED region ========== */
    for (line_num, bed_result) in bed_reader.records().enumerate() {
        let bed_record = bed_result.map_err(|e| Error::BedRecordReadFailed {
            line: line_num + 1,
            details: e.to_string(),
        })?;
        let region = parse_bed_record(&bed_record)?;

        total_regions += 1;

        if !region.is_valid_motif() {
            skipped_invalid_regions += 1;
            debug!(
                "Skipping region {} with invalid motif length {}",
                region.region_id(),
                region.motif_length(),
            );
            continue;
        }

        /* ===== STEP 1: Remove and write variants before this region ===== */
        while let Some(variant_info) = variant_window.front() {
            if variant_info.chrom < region.chrom {
                let variant_info = variant_window.pop_front().unwrap();
                write_variant(&mut writer, variant_info, &mut total_annotated_indels)?;
            } else if variant_info.chrom == region.chrom {
                let pos = variant_info.record.pos() as u64;
                let alleles = variant_info.record.alleles();
                let ref_allele = alleles[0];

                let max_indel_pos = (1..alleles.len())
                    .filter_map(|alt_idx| {
                        let alt_allele = alleles[alt_idx];
                        calculate_indel_position(pos, ref_allele, alt_allele)
                    })
                    .max()
                    .unwrap_or(pos);

                if max_indel_pos < region.start {
                    let variant_info = variant_window.pop_front().unwrap();
                    write_variant(&mut writer, variant_info, &mut total_annotated_indels)?;
                } else {
                    break;
                }
            } else {
                break;
            }
        }

        /* ===== STEP 2: Load variants until we pass region ===== */
        loop {
            if let Some(variant_info) = variant_window.back() {
                if variant_info.chrom > region.chrom {
                    break;
                } else if variant_info.chrom == region.chrom {
                    let pos = variant_info.record.pos() as u64;
                    let alleles = variant_info.record.alleles();
                    let ref_allele = alleles[0];

                    let min_indel_pos = (1..alleles.len())
                        .filter_map(|alt_idx| {
                            let alt_allele = alleles[alt_idx];
                            calculate_indel_position(pos, ref_allele, alt_allele)
                        })
                        .min()
                        .unwrap_or(pos);

                    if min_indel_pos >= region.end {
                        break;
                    }
                }
            }

            let mut next_record = input_vcf.empty_record();
            match input_vcf.read(&mut next_record) {
                None => break,
                Some(Err(e)) => {
                    return Err(Error::VcfRecordReadFailed {
                        details: e.to_string(),
                    }
                    .into());
                }
                Some(Ok(())) => {
                    let chrom = get_chrom(&next_record, &header_view)?;
                    let pos = next_record.pos();

                    if pos < 0 {
                        debug!(
                            "Skipping malformed VCF record with invalid position: {}:{}",
                            chrom,
                            pos + 1
                        );
                        continue;
                    }

                    variant_window.push_back(VariantInWindow {
                        record: next_record,
                        chrom,
                        matching_regions: Vec::new(),
                    });
                }
            }
        }

        /* ===== STEP 3: Accumulate region IDs for overlapping variants ===== */
        let mut found_perfect_indel_in_region = false;

        for variant_info in &mut variant_window {
            if variant_info.chrom != region.chrom {
                continue;
            }

            if !seen_any_chrom_overlap {
                seen_any_chrom_overlap = true;
            }

            let allele_count = variant_info.record.allele_count() as usize;
            for alt_idx in 0..(allele_count - 1) {
                if !variant_overlaps_region(&variant_info.record, &region, alt_idx) {
                    continue;
                }

                if should_include_variant(&variant_info.record, &header_view, alt_idx, &region)? {
                    let region_id = region.region_id();

                    if !variant_info.matching_regions.contains(&region_id) {
                        variant_info.matching_regions.push(region_id.clone());
                    }

                    found_perfect_indel_in_region = true;
                    break;
                }
            }
        }

        /* ===== Inject dummy indel if no perfect indel found ===== */
        if !found_perfect_indel_in_region {
            inject_dummy_deletion(&mut writer, &region)?;
            total_dummy_indels += 1;
        }
    }

    /* ========== Write remaining variants in window ========== */
    while let Some(variant_info) = variant_window.pop_front() {
        write_variant(&mut writer, variant_info, &mut total_annotated_indels)?;
    }

    /* ========== Finalization ========== */
    if !seen_any_chrom_overlap {
        return Err(Error::MsiVcfChromMismatch.into());
    }

    Ok(PreprocessingStats {
        total_regions,
        valid_regions: total_regions - skipped_invalid_regions,
        annotated_indels: total_annotated_indels,
        dummy_indels: total_dummy_indels,
    })
}

/* ================================================ */

#[cfg(test)]
mod tests {
    use super::*;

    use crate::utils::bcf_utils::tests::{
        create_test_vcf, read_first_record_simple, TestVcfConfig,
    };

    /* ============ Bed Helper  ====================== */

    // fn create_multi_region_bed() -> NamedTempFile {
    //     let tmp_bed = NamedTempFile::new().unwrap();
    //     writeln!(tmp_bed.as_file(), "chr1\t97\t104\t7xT").unwrap(); // pos 100 inside [97, 104)
    //     writeln!(tmp_bed.as_file(), "chr1\t197\t204\t7xT").unwrap(); // pos 200 inside [197, 204)
    //     writeln!(tmp_bed.as_file(), "chr2\t147\t154\t7xT").unwrap(); // pos 150 inside [147, 154)
    //     writeln!(tmp_bed.as_file(), "chrX\t172\t179\t7xT").unwrap(); // pos 175 inside [172, 179)
    //     tmp_bed
    // }

    /* ====== variant_overlaps_region tests ========== */

    #[test]
    fn test_variant_overlaps_region() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"AT"],
            ..Default::default()
        });

        let (_reader, record) = read_first_record_simple(tmp_vcf.path());

        // VCF POS=99 (0-based), REF=A, ALT=AT
        // Anchor = A (1 base)
        // Indel position = 99 + 1 = 100

        // At start (inclusive) - region [100, 200), indel at 100
        assert!(variant_overlaps_region(
            &record,
            &BedRegion {
                chrom: "chr1".to_string(),
                start: 100,
                end: 200,
                motif: "T".to_string(),
            },
            0
        ));

        // Before region - region [101, 200), indel at 100
        assert!(!variant_overlaps_region(
            &record,
            &BedRegion {
                chrom: "chr1".to_string(),
                start: 101,
                end: 200,
                motif: "T".to_string(),
            },
            0
        ));

        // At end (exclusive) - region [99, 100), indel at 100
        assert!(!variant_overlaps_region(
            &record,
            &BedRegion {
                chrom: "chr1".to_string(),
                start: 99,
                end: 100,
                motif: "T".to_string(),
            },
            0
        ));
    }

    #[test]
    fn test_variant_overlaps_region_multiple_anchors() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"TGCCT",
            alt_alleles: vec![b"TG"],
            ..Default::default()
        });

        let (_reader, record) = read_first_record_simple(tmp_vcf.path());

        // VCF POS=99, REF=TGCCT, ALT=TG
        // Anchor = TG (2 bases)
        // Indel position = 99 + 2 = 101

        // Region [101, 105): indel at 101: INSIDE
        assert!(variant_overlaps_region(
            &record,
            &BedRegion {
                chrom: "chr1".to_string(),
                start: 101,
                end: 105,
                motif: "CCT".to_string(),
            },
            0
        ));

        // Region [102, 105): indel at 101: OUTSIDE
        assert!(!variant_overlaps_region(
            &record,
            &BedRegion {
                chrom: "chr1".to_string(),
                start: 102,
                end: 105,
                motif: "CCT".to_string(),
            },
            0
        ));
    }

    #[test]
    fn test_variant_overlaps_region_complex_variant_tails_non_empty() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"ATT",
            alt_alleles: vec![b"AG"],
            ..Default::default()
        });

        let (_reader, record) = read_first_record_simple(tmp_vcf.path());

        // Complex variant returns None: should NOT overlap
        assert!(!variant_overlaps_region(
            &record,
            &BedRegion {
                chrom: "chr1".to_string(),
                start: 99,
                end: 200,
                motif: "T".to_string(),
            },
            0
        ));
    }

    #[test]
    fn test_variant_overlaps_region_snv() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"T"],
            ..Default::default()
        });

        let (_reader, record) = read_first_record_simple(tmp_vcf.path());

        // SNV returns None: should NOT overlap
        assert!(!variant_overlaps_region(
            &record,
            &BedRegion {
                chrom: "chr1".to_string(),
                start: 99,
                end: 100,
                motif: "A".to_string(),
            },
            0
        ));
    }

    #[test]
    fn test_variant_overlaps_region_multi_alt() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"GCCT",
            alt_alleles: vec![b"G", b"GCCTCCT"],
            ..Default::default()
        });

        let (_reader, record) = read_first_record_simple(tmp_vcf.path());

        // VCF POS=99
        // ALT1: GCCT -> G, anchor=1, indel_pos=100
        // ALT2: GCCT -> GCCTCCT, anchor=4, indel_pos=103

        let region1 = BedRegion {
            chrom: "chr1".to_string(),
            start: 100,
            end: 102,
            motif: "CCT".to_string(),
        };

        let region2 = BedRegion {
            chrom: "chr1".to_string(),
            start: 103,
            end: 110,
            motif: "CCT".to_string(),
        };

        // ALT1 overlaps region1 (indel at 100)
        assert!(variant_overlaps_region(&record, &region1, 0));

        // ALT1 does NOT overlap region2
        assert!(!variant_overlaps_region(&record, &region2, 0));

        // ALT2 does NOT overlap region1
        assert!(!variant_overlaps_region(&record, &region1, 1));

        // ALT2 overlaps region2 (indel at 103)
        assert!(variant_overlaps_region(&record, &region2, 1));
    }
}
