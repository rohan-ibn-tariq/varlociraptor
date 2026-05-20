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
use log::{debug, info};
use rust_htslib::bcf::{self, header::HeaderView, Read, Writer};

use crate::errors::Error;
use crate::utils::bcf_utils::{
    get_chrom, get_svlen, is_breakend, is_reference_allele, is_spanning_deletion, is_symbolic,
};
use crate::utils::genomics::{
    calculate_anchor_length, calculate_indel_position, is_clean_indel, is_indel,
};
use crate::utils::ms_bed::{parse_bed_record, BedRegion};

use super::variant_analysis::should_include_variant;

/* ============ Data Structures =================== */

/// Variant in sliding window with accumulated region annotations.
struct VariantInWindow {
    /// The VCF record
    record: bcf::Record,
    /// Chromosome name
    chrom: String,
    /// Accumulated region IDs for MS indels (empty if not MS indel)
    matching_regions: Vec<String>,
}

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

/* =========== Writer Helper Functions ============ */

/// Inject a dummy deletion for a region with no observed perfect indels.
///
/// Creates a hypothetical deletion of one motif unit positioned after the
/// first repeat. Uses the last base of the first motif as anchor, avoiding
/// the need for flanking sequence outside the region.
///
/// # Arguments
/// * `writer` - VCF writer
/// * `region` - BED region requiring dummy indel
/// * `header` - VCF header
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if chromosome not found or motif is empty
///
/// # Example
/// assert!(inject_dummy_deletion(&mut writer, &region, &header).is_ok());
fn inject_dummy_deletion(writer: &mut Writer, region: &BedRegion) -> Result<()> {
    let mut record = writer.empty_record();

    let header = writer.header();

    let rid =
        header
            .name2rid(region.chrom.as_bytes())
            .map_err(|_| Error::MsiChromosomeNotFound {
                chrom: region.chrom.clone(),
            })?;
    record.set_rid(Some(rid));

    let motif_bytes = region.motif.as_bytes();

    if motif_bytes.is_empty() {
        return Err(Error::MsiBedMotifInvalid {
            motif: "(empty)".to_string(),
        }
        .into());
    }

    let deletion_pos = region.start + (motif_bytes.len() as u64) - 1;
    record.set_pos(deletion_pos as i64);

    let anchor = vec![motif_bytes[motif_bytes.len() - 1]];

    let mut ref_allele = anchor.clone();
    ref_allele.extend_from_slice(motif_bytes);

    let alt_allele = anchor;

    record.set_alleles(&[&ref_allele, &alt_allele])?;

    let region_id = region.region_id();
    record.push_info_string(b"REGION_ID", &[region_id.as_bytes()])?;

    writer.write(&record)?;

    Ok(())
}

/// Write variant to output, with region annotations if present.
///
/// All variants are written to maintain complete VCF. Perfect MS indels
/// receive `INFO/REGION_ID` annotation with comma-separated region IDs.
///
/// # Arguments
/// * `writer` - VCF writer
/// * `variant_info` - Variant with accumulated region annotations
/// * `counter` - Counter for annotated MS indels (incremented if annotations present)
fn write_variant(
    writer: &mut Writer,
    variant_info: VariantInWindow,
    counter: &mut usize,
) -> Result<()> {
    let mut output_record = variant_info.record.clone();

    // Clear any existing REGION_ID annotations from previous preprocessing.
    output_record.clear_info_string(b"REGION_ID")?;

    if !variant_info.matching_regions.is_empty() {
        let region_id_bytes: Vec<&[u8]> = variant_info
            .matching_regions
            .iter()
            .map(|s| s.as_bytes())
            .collect();

        output_record.push_info_string(b"REGION_ID", &region_id_bytes)?;
        *counter += 1;
    }

    writer.write(&output_record)?;

    Ok(())
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
    use std::io::Write;

    use rust_htslib::bcf::{self, Read};
    use tempfile::NamedTempFile;

    use crate::utils::bcf_utils::tests::{
        create_multi_chromosome_vcf, create_test_record, create_test_vcf, read_first_record,
        read_first_record_simple, TestVcfConfig,
    };
    use crate::utils::stats::TEST_EPSILON;

    /* ============ Bed Helper  ====================== */

    // fn create_multi_region_bed() -> NamedTempFile {
    //     let tmp_bed = NamedTempFile::new().unwrap();
    //     writeln!(tmp_bed.as_file(), "chr1\t97\t104\t7xT").unwrap(); // pos 100 inside [97, 104)
    //     writeln!(tmp_bed.as_file(), "chr1\t197\t204\t7xT").unwrap(); // pos 200 inside [197, 204)
    //     writeln!(tmp_bed.as_file(), "chr2\t147\t154\t7xT").unwrap(); // pos 150 inside [147, 154)
    //     writeln!(tmp_bed.as_file(), "chrX\t172\t179\t7xT").unwrap(); // pos 175 inside [172, 179)
    //     tmp_bed
    // }

    /* ============ VCF Helpers  ======================= */

    /// Create minimal VCF header for dummy indel tests
    fn create_minimal_vcf_header() -> bcf::Header {
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(
            br##"##INFO=<ID=REGION_ID,Number=.,Type=String,Description="BED region ID">"##,
        );
        header
    }

    /* ====== write_variant tests ==================== */

    #[test]
    fn test_write_variant_with_annotation() {
        let tmp: NamedTempFile = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let record = create_test_record(&writer, 0, 1000, b"ACAG", b"ACAGCAG");

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions: vec!["chr1:1000-1020".to_string()],
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        assert_eq!(counter, 1);

        let (_reader, record) = read_first_record_simple(tmp.path());

        let region_id = record.info(b"REGION_ID").string().unwrap();
        assert!(region_id.is_some());
        assert_eq!(region_id.unwrap()[0], b"chr1:1000-1020");
    }

    #[test]
    fn test_write_variant_with_multiple_regions() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let record = create_test_record(&writer, 0, 1015, b"ACAG", b"ACAGCAG");

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions: vec!["chr1:1000-1020".to_string(), "chr1:1010-1030".to_string()],
        };

        let mut counter = 1;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        assert_eq!(counter, 2);

        let (_reader, record) = read_first_record_simple(tmp.path());

        let region_ids = record.info(b"REGION_ID").string().unwrap().unwrap();
        assert_eq!(region_ids.len(), 2);
        assert_eq!(region_ids[0], b"chr1:1000-1020");
        assert_eq!(region_ids[1], b"chr1:1010-1030");
    }

    #[test]
    fn test_write_variant_without_annotation() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let record = create_test_record(&writer, 0, 2000, b"A", b"T");

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions: Vec::new(),
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        assert_eq!(counter, 0);

        let (_reader, record) = read_first_record_simple(tmp.path());

        let region_id = record.info(b"REGION_ID").string().unwrap();
        assert!(region_id.is_none());

        assert_eq!(record.pos(), 2000);
        let alleles = record.alleles();
        assert_eq!(alleles[0], b"A");
        assert_eq!(alleles[1], b"T");
    }

    #[test]
    fn test_write_variant_preserves_variant_data() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let record = create_test_record(&writer, 0, 5000, b"GCAG", b"G");

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions: vec!["chr1:5000-5020".to_string()],
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        let (_reader, record) = read_first_record_simple(tmp.path());

        assert_eq!(record.pos(), 5000);

        let alleles = record.alleles();
        assert_eq!(alleles[0], b"GCAG");
        assert_eq!(alleles[1], b"G");

        let region_id = record.info(b"REGION_ID").string().unwrap().unwrap();
        assert_eq!(region_id[0], b"chr1:5000-5020");
    }

    #[test]
    fn test_write_variant_removes_existing_region_id_non_ms() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let mut record = create_test_record(&writer, 0, 2000, b"A", b"T");
        record
            .push_info_string(b"REGION_ID", &[b"old:value"])
            .unwrap();

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions: Vec::new(),
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        assert_eq!(counter, 0);

        let (_reader, record) = read_first_record_simple(tmp.path());

        let region_id = record.info(b"REGION_ID").string().unwrap();
        assert!(region_id.is_none(), "REGION_ID should be removed");
    }

    #[test]
    fn test_write_variant_replaces_existing_region_id_ms_indel() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let mut record = create_test_record(&writer, 0, 1000, b"ACAG", b"ACAGCAG");
        record
            .push_info_string(b"REGION_ID", &[b"old:value"])
            .unwrap();

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_regions: vec!["chr1:1000-1020".to_string()],
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        assert_eq!(counter, 1);

        let (_reader, record) = read_first_record_simple(tmp.path());

        let region_id = record.info(b"REGION_ID").string().unwrap().unwrap();
        assert_eq!(region_id.len(), 1);
        assert_eq!(region_id[0], b"chr1:1000-1020");
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

    /* ====== inject_dummy_deletion tests ============ */

    #[test]
    fn test_inject_dummy_deletion_simple_motif() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 1000,
            end: 1021,
            motif: "CAG".to_string(),
        };

        inject_dummy_deletion(&mut writer, &region).unwrap();
        drop(writer);

        let mut reader = bcf::Reader::from_path(tmp.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        let alleles = record.alleles();

        assert_eq!(record.pos(), 1002, "Position: last base of 1st CAG");
        assert_eq!(alleles[0], b"GCAG", "REF: anchor G + motif CAG");
        assert_eq!(alleles[1], b"G", "ALT: just anchor G");

        let region_id = record.info(b"REGION_ID").string().unwrap();
        assert_eq!(region_id.as_ref().unwrap()[0], b"chr1:1000-1021");
    }

    #[test]
    fn test_inject_dummy_deletion_different_motifs() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let cases = vec![
            ("A", 1000, 1020, 1000, b"AA" as &[u8], b"A" as &[u8]),
            ("AT", 1100, 1120, 1101, b"TAT" as &[u8], b"T" as &[u8]),
            ("CAG", 1200, 1221, 1202, b"GCAG" as &[u8], b"G" as &[u8]),
            ("AAAG", 1300, 1320, 1303, b"GAAAG" as &[u8], b"G" as &[u8]),
        ];

        for (motif, start, end, _, _, _) in &cases {
            let region = BedRegion {
                chrom: "chr1".to_string(),
                start: *start,
                end: *end,
                motif: motif.to_string(),
            };
            inject_dummy_deletion(&mut writer, &region).unwrap();
        }
        drop(writer);

        let mut reader = bcf::Reader::from_path(tmp.path()).unwrap();
        for (motif, _, _, expected_pos, expected_ref, expected_alt) in &cases {
            let record = reader.records().next().unwrap().unwrap();
            let alleles = record.alleles();

            assert_eq!(
                record.pos(),
                *expected_pos as i64,
                "Position mismatch for motif {}",
                motif
            );
            assert_eq!(
                alleles[0], *expected_ref,
                "REF mismatch for motif {}",
                motif
            );
            assert_eq!(
                alleles[1], *expected_alt,
                "ALT mismatch for motif {}",
                motif
            );
        }
    }

    #[test]
    fn test_inject_dummy_deletion_chromosome_not_found() {
        let tmp = NamedTempFile::new().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");

        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let region = BedRegion {
            chrom: "chr2".to_string(), // Not in header!
            start: 1000,
            end: 1021,
            motif: "CAG".to_string(),
        };

        let result = inject_dummy_deletion(&mut writer, &region);

        assert!(result.is_err());
        assert!(format!("{}", result.unwrap_err()).contains("chr2"));
    }

    #[test]
    fn test_inject_dummy_deletion_empty_motif() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 1000,
            end: 1020,
            motif: "".to_string(), // Empty!
        };

        let result = inject_dummy_deletion(&mut writer, &region);

        assert!(result.is_err());
        assert!(format!("{}", result.unwrap_err()).contains("motif"));
    }
}
