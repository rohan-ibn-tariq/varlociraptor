//! intersection.rs
//!
//! Streaming intersection of BED regions with VCF variants.
//!
//! This module provides:
//! 1. `variant_overlaps_region` function to determine if a variant's indel overlaps a BED region,
//!     accounting for VCF anchor bases and filtering out complex variants and SNVs.
//! 2. `process_and_annotate` function as a core streaming algorithm for MSI preprocessing:
//!     - Maintains a sliding window of VCF variants for memory efficiency
//!     - Coordinates variant analysis (via variant_analysis module)
//!     - Coordinates VCF output (via writer module)
//!     - Returns processing statistics
//!
//!     The streaming approach processes BED regions sequentially while maintaining
//!     a sliding window of VCF variants, enabling memory-efficient analysis of
//!     large datasets.
//!
//! Note: This module assumes that both the BED and VCF files are sorted by chromosome
//! (lexicographically) and position, which is a common requirement for genomic analyses.
//!

use std::collections::VecDeque;
use std::path::Path;

use anyhow::{Context, Result};
use bio::io::bed;
use log::debug;
use rust_htslib::bcf::{self, Read};

use crate::errors::Error;
use crate::utils::aux_info::AuxInfoCollector;
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

impl PreprocessingStats {
    /// Log a summary of the preprocessing statistics.
    /// - Total BED regions processed
    /// - Valid regions (1-6 bp motif)
    /// - Annotated MS indels
    /// - Dummy indels injected
    pub fn log_stats(&self) {
        info!("  Total BED regions: {}", self.total_regions);
        info!("  Valid regions (1-6bp motif): {}", self.valid_regions);
        info!("  Annotated MS indels: {}", self.annotated_indels);
        info!("  Dummy indels injected: {}", self.dummy_indels);
    }
}

/* ============ Functions =========================== */

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
/// assert!(variant_overlaps_region(&record, &region, 0));
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
/// 1. Stream through BED regions sequentially (bounded memory usage)
/// 2. Maintain sliding window of VCF variants around current region
/// 3. For each BED region:
///    - Remove and write variants before region (past window boundary)
///    - Load variants until past region end
///    - Check each variant for perfect repeat status (via variant_analysis::should_include_variant)
///    - Accumulate REGION_ID annotations for overlapping perfect MS indels
///    - Inject dummy deletion if no perfect indel found (via writer::inject_dummy_deletion)
/// 4. Write all remaining variants after BED exhausted
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
/// # Delegation
/// - Variant filtering: `variant_analysis::should_include_variant`
/// - VCF writing: `writer::write_variant`, `writer::inject_dummy_deletion`
///
/// # Arguments
/// * `input_vcf` - Opened BCF reader for input VCF
/// * `bed_path` - Path to microsatellite regions BED file (format: chrom, start, end, NxMOTIF)
/// * `writer` - Opened BCF writer for output
/// * `aux_info_collector` - Collector for auxiliary INFO fields to propagate (used for keeping
///    track of which INFO fields to copy from input to output)
///
/// # Returns
/// * `Ok(PreprocessingStats)` with processing statistics
/// * `Err` if processing fails
///
/// # Example
/// assert!(process_and_annotate(&mut input_vcf, bed_path, &mut writer, &aux_info_collector).is_ok());
pub(super) fn process_and_annotate(
    input_vcf: &mut bcf::Reader,
    bed_path: &Path,
    mut writer: &mut bcf::Writer,
    aux_info_collector: &AuxInfoCollector,
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
                        matching_region: None,
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

                    if let Some(existing) = &variant_info.matching_region {
                        return Err(Error::MsiBedRegionsOverlapping {
                            pos: variant_info.record.pos(),
                            existing_region: existing.clone(),
                            new_region: region_id,
                        }
                        .into());
                    }

                    variant_info.matching_region = Some(region_id);
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

/* =============== Tests ========================== */

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    use tempfile::NamedTempFile;

    use crate::preprocessing::microsatellite_instability::header::prepare_header;
    use crate::utils::aux_info::tests::make_aux_collector;
    use crate::utils::bcf_utils::tests::{
        create_test_vcf, read_first_record_simple, TestVcfConfig,
    };

    // Helper: Create BED file
    fn create_bed_file(regions: &[(&str, u64, u64, &str)]) -> NamedTempFile {
        let mut tmp = NamedTempFile::new().unwrap();
        for (chrom, start, end, motif) in regions {
            writeln!(tmp, "{}\t{}\t{}\t{}", chrom, start, end, motif).unwrap();
        }
        tmp.flush().unwrap();
        tmp
    }

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

    /* ====== process_and_annotate tests ============= */

    #[test]
    fn test_process_and_annotate_basic_annotation() {
        // VCF: chr1:99 ACAG to ACAGCAG (perfect CAG insertion)
        // BED: chr1:100-121 7xCAG
        // Expected: 1 annotated, 0 dummy

        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"ACAG",
            alt_alleles: vec![b"ACAGCAG"],
            ..Default::default()
        });
        let aux = make_aux_collector(tmp_vcf.path(), &[]);

        let tmp_bed = create_bed_file(&[("chr1", 100, 121, "7xCAG")]);
        let tmp_output = NamedTempFile::new().unwrap();

        let mut input_vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = prepare_header(input_vcf.header(), tmp_bed.path(), &aux).unwrap();
        let mut writer =
            bcf::Writer::from_path(tmp_output.path(), &header, true, bcf::Format::Vcf).unwrap();

        let stats =
            process_and_annotate(&mut input_vcf, tmp_bed.path(), &mut writer, &aux).unwrap();

        assert_eq!(stats.annotated_indels, 1);
        assert_eq!(stats.dummy_indels, 0);
        assert_eq!(stats.total_regions, 1);
        assert_eq!(stats.valid_regions, 1);
    }

    #[test]
    fn test_process_and_annotate_dummy_injection() {
        // VCF: chr1:200 A to T (SNV, not MS indel)
        // BED: chr1:100-121 7xCAG
        // Expected: 0 annotated, 1 dummy

        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"T"],
            ..Default::default()
        });
        let aux = make_aux_collector(tmp_vcf.path(), &[]);

        let tmp_bed = create_bed_file(&[("chr1", 100, 121, "7xCAG")]);
        let tmp_output = NamedTempFile::new().unwrap();

        let mut input_vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = prepare_header(input_vcf.header(), tmp_bed.path(), &aux).unwrap();
        let mut writer =
            bcf::Writer::from_path(tmp_output.path(), &header, true, bcf::Format::Vcf).unwrap();

        let stats =
            process_and_annotate(&mut input_vcf, tmp_bed.path(), &mut writer, &aux).unwrap();

        assert_eq!(stats.annotated_indels, 0);
        assert_eq!(stats.dummy_indels, 1);
        assert_eq!(stats.total_regions, 1);
        assert_eq!(stats.valid_regions, 1);
    }

    #[test]
    fn test_process_and_annotate_invalid_motif_skipped() {
        // BED: chr1:100-160 20xCAGCAGCAG (invalid: motif >6bp)
        // BED: chr1:200-221 7xCAG (valid, but no overlap with variant at pos 99)
        // Expected: 2 total regions, 1 valid processed, 1 dummy (invalid skipped, valid gets dummy indel)

        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"AATAT"],
            ..Default::default()
        });
        let aux = make_aux_collector(tmp_vcf.path(), &[]);

        let tmp_bed = create_bed_file(&[
            ("chr1", 100, 160, "20xCAGCAGCAG"), // Invalid (motif >6bp)
            ("chr1", 200, 221, "7xCAG"),        // Valid but no variant overlap
        ]);
        let tmp_output = NamedTempFile::new().unwrap();

        let mut input_vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = prepare_header(input_vcf.header(), tmp_bed.path(), &aux).unwrap();
        let mut writer =
            bcf::Writer::from_path(tmp_output.path(), &header, true, bcf::Format::Vcf).unwrap();

        let stats =
            process_and_annotate(&mut input_vcf, tmp_bed.path(), &mut writer, &aux).unwrap();

        assert_eq!(stats.total_regions, 2);
        assert_eq!(stats.valid_regions, 1);
        assert_eq!(stats.dummy_indels, 1);
    }

    #[test]
    fn test_process_and_annotate_imperfect_indel_gets_dummy() {
        // VCF: chr1:99 ACAG to ACAGCAT (imperfect: not pure CAG repeat)
        // BED: chr1:100-121 7xCAG
        // Expected: 0 annotated, 1 dummy (imperfect doesn't count)

        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"ACAG",
            alt_alleles: vec![b"ACAGCAT"],
            ..Default::default()
        });
        let aux = make_aux_collector(tmp_vcf.path(), &[]);

        let tmp_bed = create_bed_file(&[("chr1", 100, 121, "7xCAG")]);
        let tmp_output = NamedTempFile::new().unwrap();

        let mut input_vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = prepare_header(input_vcf.header(), tmp_bed.path(), &aux).unwrap();
        let mut writer =
            bcf::Writer::from_path(tmp_output.path(), &header, true, bcf::Format::Vcf).unwrap();

        let stats =
            process_and_annotate(&mut input_vcf, tmp_bed.path(), &mut writer, &aux).unwrap();

        assert_eq!(stats.annotated_indels, 0);
        assert_eq!(stats.dummy_indels, 1);
    }

    #[test]
    fn test_process_and_annotate_overlapping_bed_regions_errors() {
        // VCF: chr1:99 ACAG to ACAGCAG (perfect CAG insertion at pos 100)
        // BED: chr1:94-106 4xCAG    (overlaps)
        //      chr1:100-121 7xCAG   (overlaps)
        //      chr1:150-158 8xA     (no overlap)
        // Expected: 1 variant annotated with 2 region IDs, 1 dummy for chr1:150-158

        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"ACAG",
            alt_alleles: vec![b"ACAGCAG"],
            ..Default::default()
        });
        let aux = make_aux_collector(tmp_vcf.path(), &[]);

        let tmp_bed = create_bed_file(&[
            ("chr1", 94, 106, "4xCAG"),
            ("chr1", 100, 121, "7xCAG"),
            ("chr1", 150, 158, "8xA"),
        ]);
        let tmp_output = NamedTempFile::new().unwrap();

        let mut input_vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = prepare_header(input_vcf.header(), tmp_bed.path(), &aux).unwrap();
        let mut writer =
            bcf::Writer::from_path(tmp_output.path(), &header, true, bcf::Format::Vcf).unwrap();

        let result = process_and_annotate(&mut input_vcf, tmp_bed.path(), &mut writer, &aux);

        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("overlaps multiple BED regions"));
    }
}
