//! writer.rs
//!
//! VCF writing operations for MSI preprocessing.
//!
//! This module provides:
//! 1. `VariantInWindow` - Structure for variants with accumulated region annotations
//! 2. `inject_dummy_deletion` - Creates dummy indels for regions without perfect repeats
//! 3. `write_variant` - Writes variants with REGION_ID annotations if applicable
//! 4.  Unit tests for both (2,3) functions, covering typical and edge cases
//!

use anyhow::Result;
use rust_htslib::bcf::{self, Writer};

use crate::errors::Error;
use crate::utils::ms_bed::BedRegion;

/* ============ Data Structures =================== */

/// Variant in sliding window with accumulated region annotations.
pub(super) struct VariantInWindow {
    /// The VCF record
    pub record: bcf::Record,
    /// Chromosome name
    pub chrom: String,
    /// Accumulated region IDs for MS indels (empty if not MS indel)
    pub matching_regions: Vec<String>,
}

/* ============ Functions ========================= */

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
pub(super) fn inject_dummy_deletion(writer: &mut Writer, region: &BedRegion) -> Result<()> {
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
pub(super) fn write_variant(
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

/* =============== Tests ========================== */

#[cfg(test)]
mod tests {
    use super::*;

    use rust_htslib::bcf::{self, Read};
    use tempfile::NamedTempFile;

    use crate::utils::bcf_utils::tests::{create_test_record, read_first_record_simple};

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
}
