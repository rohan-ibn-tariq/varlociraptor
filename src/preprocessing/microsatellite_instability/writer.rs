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

use crate::constants::{
    MSI_DUMMY_HEADER, MSI_DUMMY_TAG, MSI_REGION_ID_HEADER, MSI_REGION_ID_TAG,
    PREPROCESS_MSI_COPY_FIELDS, PREPROCESS_MSI_OMIT_AUX,
};
use crate::errors::Error;
use crate::utils::aux_info::AuxInfo;
use crate::utils::bcf_utils::copy_info_fields;
use crate::utils::ms_bed::BedRegion;

/* ============ Data Structures =================== */

/// Variant in sliding window with accumulated region annotations.
pub(super) struct VariantInWindow {
    /// The VCF record
    pub record: bcf::Record,
    /// Chromosome name
    pub chrom: String,
    /// Region ID for MS indels (None if not MS indel)
    pub matching_region: Option<String>,
    /// User-specified INFO fields collected from input record.
    /// Written to output via aux_info.write() during write_variant.
    /// Fields in PREPROCESS_MSI_OMIT_AUX are excluded to prevent
    /// double-writing with copy_info_fields().
    pub aux_info: AuxInfo,
}

/* ============ Functions ========================= */

/// Inject a dummy deletion for a region with no observed perfect indels.
///
/// Creates a hypothetical deletion of one motif unit positioned after the
/// first repeat. Uses the last base of the first motif as anchor, avoiding
/// the need for flanking sequence outside the region.
///
/// Output record contains only REGION_ID and MSI_DUMMY INFO fields.
/// No FORMAT or sample data is written.
///
/// # Arguments
/// * `writer` - VCF writer
/// * `region` - BED region requiring dummy indel
///
/// # Returns
/// * `Ok(())` on success
/// * `Err` if chromosome not found or motif is empty
///
/// # Example
/// assert!(inject_dummy_deletion(&mut writer, &region).is_ok());
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
    record.push_info_string(MSI_REGION_ID_TAG, &[region_id.as_bytes()])?;

    record.push_info_flag(MSI_DUMMY_TAG)?;

    writer.write(&record)?;

    Ok(())
}

/// Write variant to output, with region annotations if present.
///
/// Creates a fresh output record containing:
/// 1. Position and alleles from input record
/// 2. Standard INFO fields copied via copy_info_fields()
/// 3. User-specified INFO fields via aux_info.write()
/// 4. REGION_ID annotation if variant overlaps a perfect MS region
///
/// FORMAT and sample data are intentionally omitted.
///
/// # Arguments
/// * `writer` - VCF writer
/// * `variant_info` - Variant with accumulated region annotation and aux info
/// * `counter` - Counter for annotated MS indels (incremented if annotations present)
pub(super) fn write_variant(
    writer: &mut Writer,
    variant_info: VariantInWindow,
    counter: &mut usize,
) -> Result<()> {
    let mut output_record = writer.empty_record();
    output_record.set_rid(variant_info.record.rid());
    output_record.set_pos(variant_info.record.pos());
    output_record.set_alleles(&variant_info.record.alleles())?;

    copy_info_fields(
        &variant_info.record,
        &mut output_record,
        &PREPROCESS_MSI_COPY_FIELDS,
    )?;

    variant_info
        .aux_info
        .write(&mut output_record, &PREPROCESS_MSI_OMIT_AUX)?;

    if let Some(region_id) = &variant_info.matching_region {
        output_record.push_info_string(MSI_REGION_ID_TAG, &[region_id.as_bytes()])?;
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

    use crate::utils::aux_info::tests::make_aux_collector;
    use crate::utils::bcf_utils::tests::{
        create_test_record, create_test_vcf, read_first_record_simple, TestVcfConfig,
    };

    /* ============ VCF Helpers  ======================= */

    /// Create minimal VCF header for dummy indel tests
    fn create_minimal_vcf_header() -> bcf::Header {
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(
            br##"##INFO=<ID=REGION_ID,Number=1,Type=String,Description="BED region ID">"##,
        );
        header.push_record(
            br##"##INFO=<ID=MSI_DUMMY,Number=0,Type=Flag,Description="Dummy deletion injected for MS region with no observed indel">"##
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
            matching_region: Some("chr1:1000-1020".to_string()),
            aux_info: AuxInfo::default(),
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
    fn test_write_variant_without_annotation() {
        let tmp = NamedTempFile::new().unwrap();
        let header = create_minimal_vcf_header();
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let record = create_test_record(&writer, 0, 2000, b"A", b"T");

        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_region: None,
            aux_info: AuxInfo::default(),
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
            matching_region: Some("chr1:5000-5020".to_string()),
            aux_info: AuxInfo::default(),
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
            matching_region: None,
            aux_info: AuxInfo::default(),
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
            matching_region: Some("chr1:1000-1020".to_string()),
            aux_info: AuxInfo::default(),
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

    #[test]
    fn test_write_variant_propagates_aux_info() {
        let (tmp_src, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"ACAG",
            alt_alleles: vec![b"ACAGCAG"],
            extra_info_fields: vec![(b"COSMIC_ID", b"1", b"String", b"COSMIC ID", b"COSM123")],
            ..Default::default()
        });

        let aux = make_aux_collector(tmp_src.path(), &["COSMIC_ID"]);
        let mut src_reader = bcf::Reader::from_path(tmp_src.path()).unwrap();
        let src_record = src_reader.records().next().unwrap().unwrap();
        let aux_info = aux.collect(&src_record).unwrap();

        let tmp = NamedTempFile::new().unwrap();
        let mut header = create_minimal_vcf_header();
        header.push_record(
            br##"##INFO=<ID=COSMIC_ID,Number=1,Type=String,Description="COSMIC ID">"##,
        );
        let mut writer = Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();

        let record = create_test_record(&writer, 0, 1000, b"ACAG", b"ACAGCAG");
        let variant_info = VariantInWindow {
            record,
            chrom: "chr1".to_string(),
            matching_region: Some("chr1:1000-1020".to_string()),
            aux_info,
        };

        let mut counter = 0;
        write_variant(&mut writer, variant_info, &mut counter).unwrap();
        drop(writer);

        let (_reader, record) = read_first_record_simple(tmp.path());
        let cosmic = record.info(b"COSMIC_ID").string().unwrap();
        assert!(cosmic.is_some(), "COSMIC_ID propagated via aux_info");
        assert_eq!(cosmic.unwrap()[0], b"COSM123");
    }
}
