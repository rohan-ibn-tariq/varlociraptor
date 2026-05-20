//! variant_analysis.rs
//!
//! Variant Analysis Utilities for MSI preprocessing.
//!
//! This module provides:
//! 1. `is_perfect_repeat` function to determine if an indel is a perfect tandem repeat of a microsatellite motif.
//! 2. `should_include_variant` function to analyze a variant-allele pair and determine if it should be included
//!     in the preprocessed output based on MSI relevance and perfect repeat status.
//!

use anyhow::Result;
use log::debug;
use rust_htslib::bcf::{self, header::HeaderView};

use crate::utils::bcf_utils::{
    get_chrom, get_svlen, is_breakend, is_reference_allele, is_spanning_deletion, is_symbolic,
};
use crate::utils::genomics::{calculate_anchor_length, is_clean_indel, is_indel};
use crate::utils::ms_bed::BedRegion;


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
pub(super) fn should_include_variant(
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


#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::bcf_utils::tests::{create_test_vcf, read_first_record, TestVcfConfig};
    use crate::utils::ms_bed::BedRegion;


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

    /* ======= should_include_variant tests ========== */

    #[test]
    fn test_should_include_variant_filters_snv() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"T"],
            ..Default::default()
        });

        let (_, header, record) = read_first_record(tmp_vcf.path());

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

        let (_, header, record) = read_first_record(tmp_vcf.path());

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

        let (_, header, record) = read_first_record(tmp_vcf.path());

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

        let (_reader, header, record) = read_first_record(tmp_vcf.path());

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

        let (_, header, record) = read_first_record(tmp_vcf.path());

        let region = BedRegion {
            chrom: "chr1".to_string(),
            start: 0,
            end: 201,
            motif: "CAG".to_string(),
        };

        let result = should_include_variant(&record, &header, 0, &region).unwrap();
        assert!(!result); // Imperfect repeat should be filtered
    }
}
