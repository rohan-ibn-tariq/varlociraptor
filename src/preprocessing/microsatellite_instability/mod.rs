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

mod intersection;

use std::path::PathBuf;

use anyhow::{Context, Result};
use log::info;
use rust_htslib::bcf;

use crate::cli::{validate_bed_file, validate_vcf_file};
use crate::errors::Error;
use crate::utils::bcf_utils::validate_vcf_file as validate_ms_vcf_file;
use crate::utils::ms_bed::validate_bed_file as validate_ms_bed_file;

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
}

/* ================================================ */

/// Main Orchestrator for MSI candidate preprocessing workflow.
pub fn preprocess_ms_candidates(config: PreprocessMSIConfig) -> Result<()> {
    Ok(())
}
