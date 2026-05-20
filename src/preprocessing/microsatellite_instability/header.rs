//! header.rs
//!
//! Output VCF header preparation for MSI preprocessing.
//!
//! This module provides:
//! 1. `prepare_header` function to create a modified VCF header for MSI preprocessing output.
//!

use std::path::Path;

use anyhow::Result;
use log::debug;
use rust_htslib::bcf;

use crate::errors::Error;
use crate::utils::ms_bed::collect_bed_chromosomes;

/// Prepare output VCF header for MSI preprocessing.
///
/// Creates a modified header from the input VCF by:
/// 1. Copying all existing fields from input header
/// 2. Adding INFO/REGION_ID field definition
/// 3. Adding contig definitions for chromosomes present in BED but missing from input VCF
///
/// # Arguments
/// * `input_header` - HeaderView from input VCF reader
/// * `bed_path` - Path to BED file (for collecting chromosome names)
///
/// # Returns
/// * `Ok(Header)` - Prepared header ready for writer creation
/// * `Err` if BED file cannot be read or has no valid regions
///
/// # Example
/// assert!(prepare_header(&input_vcf.header(), Path::new("ms_regions.bed")).is_ok());
pub(super) fn prepare_header(
    input_header: &bcf::header::HeaderView,
    bed_path: &Path,
) -> Result<bcf::Header> {
    let mut header = bcf::Header::from_template(input_header);

    header.push_record(
        br##"##INFO=<ID=REGION_ID,Number=.,Type=String,Description="Comma-separated BED region IDs for MS indels">"##
    );

    let bed_chroms = collect_bed_chromosomes(bed_path)?;

    if bed_chroms.is_empty() {
        return Err(Error::BedFileNoValidRegions.into());
    }

    for chrom in &bed_chroms {
        if input_header.name2rid(chrom.as_bytes()).is_err() {
            let contig_line = format!("##contig=<ID={}>", chrom);
            header.push_record(contig_line.as_bytes());
            debug!("Added missing contig to VCF header: {}", chrom);
        }
    }

    Ok(header)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    use rust_htslib::bcf::{self, Format, Read};
    use tempfile::NamedTempFile;

    /// Helper: Create minimal test VCF file with given contigs
    fn create_test_vcf(contigs: &[&str]) -> NamedTempFile {
        let mut tmp = NamedTempFile::new().unwrap();
        write!(tmp, "##fileformat=VCFv4.2\n").unwrap();
        for contig in contigs {
            writeln!(tmp, "##contig=<ID={}>", contig).unwrap();
        }
        writeln!(tmp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        tmp.flush().unwrap();
        tmp
    }

    /// Helper: Create test BED file with given chromosomes
    /// NOTE: Each region is 100-121 with a 7xCAG motif
    /// (inappropriate for MSI but good for testing header prep)
    fn create_test_bed(chroms: &[&str]) -> NamedTempFile {
        let mut tmp = NamedTempFile::new().unwrap();
        for chrom in chroms {
            writeln!(tmp, "{}\t100\t121\t7xCAG", chrom).unwrap();
        }
        tmp.flush().unwrap();
        tmp
    }

    /// Helper: Write header to file and read back as string for assertions
    fn header_to_string(header: &bcf::Header) -> String {
        let tmp = NamedTempFile::new().unwrap();
        let writer = bcf::Writer::from_path(tmp.path(), header, true, Format::Vcf).unwrap();
        drop(writer);
        std::fs::read_to_string(tmp.path()).unwrap()
    }

    #[test]
    fn test_prepare_header_adds_region_id_field() {
        let vcf_file = create_test_vcf(&["chr1"]);
        let bed_file = create_test_bed(&["chr1"]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), bed_file.path()).unwrap();

        let content = header_to_string(&result);
        assert!(content.contains("ID=REGION_ID"));
    }

    #[test]
    fn test_prepare_header_adds_missing_contigs() {
        let vcf_file = create_test_vcf(&["chr1"]);
        let bed_file = create_test_bed(&["chr1", "chr2"]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), bed_file.path()).unwrap();

        let content = header_to_string(&result);
        assert!(content.contains("ID=REGION_ID"));
        assert!(content.contains("contig=<ID=chr2>"));
    }

    #[test]
    fn test_prepare_header_errors_on_empty_bed() {
        let vcf_file = create_test_vcf(&["chr1"]);
        let empty_bed = NamedTempFile::new().unwrap();

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), empty_bed.path());

        assert!(result.is_err());
    }

    #[test]
    fn test_prepare_header_preserves_existing_contigs() {
        let vcf_file = create_test_vcf(&["chr1", "chr2"]);
        let bed_file = create_test_bed(&["chr1"]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), bed_file.path()).unwrap();

        let content = header_to_string(&result);
        // chr2 should still be in header (not removed)
        assert!(content.contains("contig=<ID=chr1>"));
        assert!(content.contains("contig=<ID=chr2>"));
    }

    #[test]
    fn test_prepare_header_no_duplicate_contigs() {
        let vcf_file = create_test_vcf(&["chr1"]);
        let bed_file = create_test_bed(&["chr1"]);

        let reader = bcf::Reader::from_path(vcf_file.path()).unwrap();
        let result = prepare_header(reader.header(), bed_file.path()).unwrap();

        let content = header_to_string(&result);

        let count = content.matches("contig=<ID=chr1>").count();
        assert_eq!(count, 1, "chr1 should appear exactly once, not duplicated");
    }
}
