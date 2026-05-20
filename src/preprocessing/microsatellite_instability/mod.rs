//! mod.rs
//!
//! Microsatellite instability (MSI) preprocessing module.
//! Enriches candidate VCF with:
//! - Dummy indel records for MS regions without variants
//! - MS region annotations (INFO/MS_REGION)
//!
//! **Note: This feature is experimental.**
//!
//! This module provides:
//! 1. @TODO: A brief overview of the MSI preprocessing workflow, e.g.:
//!

mod header;
mod intersection;
mod variant_analysis;
mod writer;

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use log::info;
use rust_htslib::bcf::{self, Format, Read};

use crate::cli::{validate_bed_file, validate_vcf_file};
use crate::errors::Error;
use crate::utils::bcf_utils::validate_vcf_file as validate_ms_vcf_file;
use crate::utils::ms_bed::validate_bed_file as validate_ms_bed_file;

use header::prepare_header;
use intersection::process_and_annotate;

/* ======== CLI CONFIGURATION ===================== */

/// Configuration for MS candidate preprocessing pipeline.
///
/// Contains all parameters needed for the MSI preprocessing workflow,
/// populated from CLI arguments.
#[derive(Debug)]
pub struct PreprocessMSIConfig {
    /// Path to BED file (sorted) with microsatellite regions.
    pub microsatellite_bed: PathBuf,
    /// Path to candidate VCF/BCF file (sorted) with variant calls (gnomAD-annotated with INFO/HETEROZYGOSITY).
    pub candidate_vcf: PathBuf,
    /// Output file path (VCF, BCF, or VCF.GZ; if omitted, writes BCF to STDOUT).
    pub output: Option<PathBuf>,
}

impl PreprocessMSIConfig {
    /// Validate PreprocessMSIConfig
    pub fn validate(&self) -> Result<()> {
        // Use existing validators from cli.rs to validate file formats.
        validate_bed_file(&self.microsatellite_bed)?;
        validate_vcf_file(&self.candidate_vcf)?;

        let mut vcf =
            bcf::Reader::from_path(&self.candidate_vcf).context("Failed to open candidate VCF")?;

        // Additional file specifc validators:
        validate_ms_bed_file(&self.microsatellite_bed)?;
        info!("Validating VCF file: {}", self.candidate_vcf.display());
        validate_ms_vcf_file(&mut vcf)?;

        Ok(())
    }

    /// Determine output format and compression from file extension.
    ///
    /// # Format Detection
    /// - `.bcf` to BCF binary (always BGZF compressed)
    /// - `.vcf.gz` or `.vcf.bgz` to VCF text, BGZF compressed
    /// - `.vcf` to VCF text, uncompressed
    /// - Anything else to VCF text, uncompressed (safest default)
    /// - None (stdout) to VCF text, uncompressed
    ///
    /// # Returns
    /// Tuple of (Format, uncompressed_flag) for writer creation
    pub fn determine_output_format(&self) -> (Format, bool) {
        match self.output.as_deref() {
            None => (Format::Vcf, true),
            Some(path) => {
                let path_str = path.to_string_lossy().to_lowercase();

                if path_str.ends_with(".bcf") {
                    (Format::Bcf, false)
                } else if path_str.ends_with(".vcf.gz") || path_str.ends_with(".vcf.bgz") {
                    (Format::Vcf, false)
                } else {
                    // .vcf and any other extensions default to uncompressed VCF
                    (Format::Vcf, true)
                }
            }
        }
    }
}

/* ================================================ */

/* ======== Main Workflow Orchestration ============ */

/// Main Orchestrator for MSI candidate preprocessing workflow.
pub fn preprocess_ms_candidates(config: PreprocessMSIConfig) -> Result<()> {
    info!("----------------------------------------------");
    info!("Step 1: Config Stats");
    info!("----------------------------------------------");

    info!("Input files:");
    info!("BED file: {}", config.microsatellite_bed.display());
    info!("VCF/BCF file: {}", config.candidate_vcf.display());

    if let Some(ref output) = config.output {
        info!("Output: {}", output.display());
    } else {
        info!("Output: STDOUT (VCF format)");
    }

    info!("----------------------------------------------");
    info!("Step 2: Streaming Intersection & Output Generation");
    info!("----------------------------------------------");
    info!("Starting step 2.");

    info!("Opening input VCF: {}", config.candidate_vcf.display());
    let mut input_vcf =
        bcf::Reader::from_path(&config.candidate_vcf).context("Failed to open input VCF")?;

    info!("Preparing output header...");
    let output_header = prepare_header(input_vcf.header(), &config.microsatellite_bed)?;

    info!("Determining output format based on output path...");
    let (output_format, output_uncompressed) = config.determine_output_format();

    info!("Creating output writer...");
    let mut writer = match config.output {
        Some(ref path) => {
            bcf::Writer::from_path(path, &output_header, output_uncompressed, output_format)?
        }
        None => bcf::Writer::from_stdout(&output_header, output_uncompressed, output_format)?,
    };

    info!("\nStarting streaming intersection...");
    let stats = process_and_annotate(&mut input_vcf, &config.microsatellite_bed, &mut writer)?;
    info!("Intersection finished and output generated.");

    info!("----------------------------------------------");
    info!("Step 3(Final): Logging Final Stats");
    info!("----------------------------------------------");
    info!("  Total BED regions: {}", stats.total_regions);
    info!("  Valid regions (1-6bp motif): {}", stats.valid_regions);
    info!("  Annotated MS indels: {}", stats.annotated_indels);
    info!("  Dummy indels injected: {}", stats.dummy_indels);

    info!("==============================================");
    info!("Preprocessing MSI complete");
    info!("==============================================");

    Ok(())
}

/* ================================================ */

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    use rust_htslib::bcf::Format;

    /// Helper function to create a test PreprocessMSIConfig with optional output path.
    fn create_test_config(output: Option<PathBuf>) -> PreprocessMSIConfig {
        PreprocessMSIConfig {
            microsatellite_bed: PathBuf::from("test.bed"),
            candidate_vcf: PathBuf::from("test.vcf"),
            output,
        }
    }

    #[test]
    fn test_determine_output_format_stdout() {
        let config = create_test_config(None);
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(uncompressed, "Stdout should be uncompressed");
    }

    #[test]
    fn test_determine_output_format_bcf() {
        let config = create_test_config(Some(PathBuf::from("output.bcf")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Bcf);
        assert!(!uncompressed, "BCF should be compressed");
    }

    #[test]
    fn test_determine_output_format_vcf_uncompressed() {
        let config = create_test_config(Some(PathBuf::from("output.vcf")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(uncompressed, "Plain .vcf should be uncompressed");
    }

    #[test]
    fn test_determine_output_format_vcf_gz() {
        let config = create_test_config(Some(PathBuf::from("output.vcf.gz")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(!uncompressed, ".vcf.gz should be BGZF compressed");
    }

    #[test]
    fn test_determine_output_format_vcf_bgz() {
        let config = create_test_config(Some(PathBuf::from("output.vcf.bgz")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(!uncompressed, ".vcf.bgz should be BGZF compressed");
    }

    #[test]
    fn test_determine_output_format_case_insensitive() {
        let config = create_test_config(Some(PathBuf::from("OUTPUT.BCF")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Bcf);
        assert!(!uncompressed);

        let config = create_test_config(Some(PathBuf::from("output.VCF.GZ")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(!uncompressed);
    }

    #[test]
    fn test_determine_output_format_unknown_extension() {
        let config = create_test_config(Some(PathBuf::from("weird.txt")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(
            uncompressed,
            "Unknown extension should default to uncompressed"
        );

        let config = create_test_config(Some(PathBuf::from("output")));
        let (format, uncompressed) = config.determine_output_format();
        assert_eq!(format, Format::Vcf);
        assert!(uncompressed, "No extension should default to uncompressed");
    }
}
