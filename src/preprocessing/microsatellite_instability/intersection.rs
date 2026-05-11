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

/* ============ Data Structures =================== */

/// Classification of indel repeat pattern relative to microsatellite motif.
#[derive(Debug, Clone, PartialEq)]
enum RepeatStatus {
    /// Variant indel is a perfect tandem repeat of the motif.
    Perfect,
    /// Variant indel does not match motif pattern or fails validation.
    NA,
}

/* ================================================ */

/* ======== Variant Analysis Functions ============ */

/// Check if an indel is a perfect tandem repeat of a microsatellite motif.
///
/// Determines if the changed sequence (insertion or deletion) consists
/// entirely of complete motif units.
///
/// # Algorithm
/// 1. Find anchor (common prefix between REF and ALT)
/// 2. Extract changed sequence after anchor
/// 3. Verify changed sequence length matches |SVLEN|
/// 4. Check if changed sequence is exact motif repeats
///
/// # Arguments
/// * `alt_seq` - Alternate allele sequence
/// * `svlen` - Structural variant length (positive=insertion, negative=deletion)
/// * `motif` - Microsatellite motif to match
/// * `ref_seq` - Reference allele sequence
///
/// # Returns
/// * `RepeatStatus::Perfect` - Changed sequence is exact motif repeats
/// * `RepeatStatus::NA` - Does not match or fails validation
///
/// # Example
/// assert_eq!(is_perfect_repeat(b"ACAGCAG", 3, "CAG", b"ACAG"), RepeatStatus::Perfect);
fn is_perfect_repeat(alt_seq: &[u8], svlen: i32, motif: &str, ref_seq: &[u8]) -> RepeatStatus {
    // 0. Handling Edge Cases
    if ref_seq.is_empty()
        || alt_seq.is_empty()
        || !ref_seq[0].eq_ignore_ascii_case(&alt_seq[0])
        || svlen == 0
    {
        return RepeatStatus::NA;
    }

    // 1. Use genomics utility to check if clean indel
    if !is_clean_indel(ref_seq, alt_seq) {
        return RepeatStatus::NA;
    }

    // 2. Finding the anchor length and absolute SVLEN
    // Note: Anchor length 0 is not errored as a valid indel
    let abs_svlen = svlen.unsigned_abs() as usize;
    let anchor_len = calculate_anchor_length(ref_seq, alt_seq);

    // 3. Extracting the changed sequence
    let changed_seq = if svlen > 0 {
        if anchor_len < alt_seq.len() {
            &alt_seq[anchor_len..]
        } else {
            return RepeatStatus::NA;
        }
    } else if anchor_len < ref_seq.len() {
        &ref_seq[anchor_len..]
    } else {
        return RepeatStatus::NA;
    };

    // Validate: changed sequence length should match SVLEN
    if changed_seq.len() != abs_svlen || changed_seq.is_empty() {
        return RepeatStatus::NA;
    }

    // 4. Check if changed sequence is a perfect repeat of the motif
    let motif_bytes: Vec<u8> = motif.bytes().map(|b| b.to_ascii_uppercase()).collect();
    let motif_len = motif_bytes.len();

    /* NOTE: In case we toggle the Error on motif.len() in bed parsing off, turn this on.*/
    // if motif_len == 0 {
    //     return RepeatStatus::NA
    // }

    if changed_seq.len() % motif_len != 0 {
        return RepeatStatus::NA;
    }

    for (i, &base) in changed_seq.iter().enumerate() {
        let expected_base = motif_bytes[i % motif_len];
        if base.to_ascii_uppercase() != expected_base {
            return RepeatStatus::NA;
        }
    }

    RepeatStatus::Perfect
}

/// Analyze a variant if it should be included in preprocessed output.
///
/// This function processes a single variant-allele pair to determine if it's
/// relevant for preprocessed output. It performs the following steps:
///
/// # Algorithm
/// 1. **Filtering**: Skip variants that aren't relevant for MSI:
///    - Reference alleles (ALT=<REF>)
///    - Symbolic alleles (<DEL>, <INS>)
///    - Breakends (complex structural variants)
///    - Spanning deletions (*)
///    - Non-indel variants
/// 2. **Check Perfect Repeat Status**:
///    - Calculate SVLEN (indel length)
///    - Verify indel is perfect tandem repeat of motif
///
/// # Arguments
/// * `record` - BCF record representing the variant
/// * `header` - BCF header for metadata access
/// * `alt_idx` - Index of the alternate allele to analyze
/// * `region` - BED region context for motif information
///
/// # Returns
/// * `Ok(true)` - Variant is a perfect MS indel at this locus
/// * `Ok(false)` - Variant should be skipped (not relevant for MSI quantification)
/// * `Err` - Error reading variant data
///
/// # Example
/// assert!(should_include_variant(&record, &header, 0, &region).unwrap());
fn should_include_variant(
    record: &bcf::Record,
    header: &HeaderView,
    alt_idx: usize,
    region: &BedRegion,
) -> Result<bool> {
    let alleles = record.alleles();
    let ref_allele = alleles[0];
    let alt_allele = alleles[alt_idx + 1]; // +1 because alleles[0] is REF

    /* 1. Filter non indel variants */
    if is_reference_allele(alt_allele)
        || is_symbolic(alt_allele)
        || is_breakend(alt_allele)
        || is_spanning_deletion(alt_allele)
        || !is_indel(ref_allele, alt_allele)
    {
        debug!(
            "Filtering non-indel variant at {}:{}",
            get_chrom(record, header)?,
            record.pos() + 1
        );
        return Ok(false);
    }

    /* 2. Check Perfect Repeat Status */
    let svlen = get_svlen(record, alt_idx, ref_allele, alt_allele)?;
    let repeat_status = is_perfect_repeat(alt_allele, svlen, &region.motif, ref_allele);

    Ok(repeat_status == RepeatStatus::Perfect)
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
fn inject_dummy_deletion(
    writer: &mut Writer,
    region: &BedRegion,
    header: &HeaderView,
) -> Result<()> {
    let mut record = writer.empty_record();

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

/* ================================================ */

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    use rust_htslib::bcf::{self, Read};
    use tempfile::NamedTempFile;

    use crate::utils::bcf_utils::tests::{
        create_multi_chromosome_vcf, create_test_vcf, TestVcfConfig,
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
            br##"##INFO=<ID=REGION_ID,Number=1,Type=String,Description="BED region ID">"##,
        );
        header
    }

    /* ========== is_perfect_repeat tests ============ */

    #[test]
    fn test_is_perfect_repeat_insertion() {
        assert_eq!(
            is_perfect_repeat(b"ACAGCAG", 3, "CAG", b"ACAG"),
            RepeatStatus::Perfect
        );
        assert_eq!(
            is_perfect_repeat(b"ACAGCAGCAG", 6, "CAG", b"ACAG"),
            RepeatStatus::Perfect
        );
        assert_eq!(
            is_perfect_repeat(b"ATT", 1, "T", b"AT"),
            RepeatStatus::Perfect
        );
        /* We are considering even low indels like in the last case. TODO: Check if discarding is better behaviour. */
    }

    #[test]
    fn test_is_perfect_repeat_deletion() {
        assert_eq!(
            is_perfect_repeat(b"ACAG", -3, "CAG", b"ACAGCAG"),
            RepeatStatus::Perfect
        );
    }

    #[test]
    fn test_is_perfect_repeat_case_insensitive() {
        assert_eq!(
            is_perfect_repeat(b"acagCAG", 3, "cag", b"acag"),
            RepeatStatus::Perfect
        );
    }

    #[test]
    fn test_is_perfect_repeat_case_insensitive_first_byte() {
        assert_eq!(
            is_perfect_repeat(b"AcagCAG", 3, "cag", b"acag"),
            RepeatStatus::Perfect
        );
    }

    #[test]
    fn test_is_perfect_repeat_not_perfect() {
        assert_eq!(
            is_perfect_repeat(b"ACAGCAT", 3, "CAG", b"ACAG"),
            RepeatStatus::NA
        );
        assert_eq!(
            is_perfect_repeat(b"ACAGCA", 2, "CAG", b"ACAG"),
            RepeatStatus::NA
        );
        assert_eq!(
            is_perfect_repeat(b"AAAGAGAGAGA", 7, "GA", b"AAAT"),
            RepeatStatus::NA
        );
        /* Last test: Special Case when tail remains in ref/alt apart from anchor.
            Here we consider it NA as it's not a clean indel.
        */
    }

    #[test]
    fn test_is_perfect_repeat_edge_cases() {
        assert_eq!(is_perfect_repeat(b"CAG", 3, "CAG", b""), RepeatStatus::NA);
        assert_eq!(is_perfect_repeat(b"", 0, "CAG", b"CAG"), RepeatStatus::NA);
        assert_eq!(
            is_perfect_repeat(b"ACAT", 0, "CAG", b"ACAG"),
            RepeatStatus::NA
        );
        assert_eq!(
            is_perfect_repeat(b"TCAG", 3, "TCAG", b"A"),
            RepeatStatus::NA
        );
    }

    /* ========== analyze_variant tests ============== */

    #[test]
    fn test_should_include_variant_filters_snv() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"T"],
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "A".to_string(),
        };

        let result = should_include_variant(&record, &header, 0, &region).unwrap();
        assert!(!result); // SNV should be filtered out
    }

    #[test]
    fn test_should_include_variant_perfect_indel() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"ACAG",
            alt_alleles: vec![b"ACAGCAG"],
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 201,
            motif: "CAG".to_string(),
        };

        let result = should_include_variant(&record, &header, 0, &region).unwrap();

        assert!(result); // Perfect indel should be included
    }

    #[test]
    fn test_should_include_variant_multi_allelic() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"T", b"ATG"],
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "TG".to_string(),
        };

        // alt_idx=0 is SNV (T) - should filter out
        assert!(!should_include_variant(&record, &header, 0, &region).unwrap());

        // alt_idx=1 is indel (ATG) - should include
        assert!(should_include_variant(&record, &header, 1, &region).unwrap());
    }

    #[test]
    fn test_should_include_variant_filters_symbolic() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"<DEL>"],
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 200,
            motif: "A".to_string(),
        };

        let result = should_include_variant(&record, &header, 0, &region).unwrap();
        assert!(!result); // Symbolic allele should be filtered out
    }

    #[test]
    fn test_should_include_variant_imperfect_repeat() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"ACAG",
            alt_alleles: vec![b"ACAGCAT"], // Not perfect CAG repeat
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 201,
            motif: "CAG".to_string(),
        };

        let result = should_include_variant(&record, &header, 0, &region).unwrap();
        assert!(!result); // Imperfect repeat should be filtered
    }

    /* ====== variant_overlaps_region tests ========== */

    #[test]
    fn test_variant_overlaps_region() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"AT"],
            ..Default::default()
        });
        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();

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

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();

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

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();

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

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();

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

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let record = reader.records().next().unwrap().unwrap();

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

        let header_view = writer.header().clone();
        inject_dummy_deletion(&mut writer, &region, &header_view).unwrap();
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
        let header_view = writer.header().clone();

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
            inject_dummy_deletion(&mut writer, &region, &header_view).unwrap();
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

        let header_view = writer.header().clone();
        let result = inject_dummy_deletion(&mut writer, &region, &header_view);

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

        let header_view = writer.header().clone();
        let result = inject_dummy_deletion(&mut writer, &region, &header_view);

        assert!(result.is_err());
        assert!(format!("{}", result.unwrap_err()).contains("motif"));
    }
}
