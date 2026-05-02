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

mod annotation;
mod dummy_indels;
mod intersection;

use std::path::PathBuf;

use anyhow::Result;
use log::info;

use crate::errors::Error;

/* ======== CLI CONFIGURATION ===================== */

/// Configuration for MS candidate preprocessing pipeline.
///
/// Contains all parameters needed for the MSI preprocessing subcommand,
/// populated from CLI arguments.
#[derive(Debug)]
pub struct PreprocessMSIConfig {
    /// Path to BED file(should be sorted) with microsatellite regions.
    pub microsatellite_bed: PathBuf,
    /// Candidate VCF/BCF file with variant calls (gnomAD-annotated with INFO/HETEROZYGOSITY).
    pub calls: PathBuf,
    /// Output VCF file (if omitted, writes to STDOUT; logs go to STDERR).
    pub output_vcf: Option<PathBuf>,
}

impl PreprocessMSIConfig {
    /// Validate PreprocessMSIConfig
    pub fn validate(&self) -> Result<()> {
        // Validate BED file
        if !self.microsatellite_bed.exists() {
            return Err(anyhow::anyhow!(
                "BED file does not exist: {}",
                self.microsatellite_bed.display()
            ));
        }

        // Validate VCF file
        if !self.calls.exists() {
            return Err(anyhow::anyhow!(
                "VCF file does not exist: {}",
                self.calls.display()
            ));
        }

        // Use existing validators from cli.rs to validate file formats.
        crate::cli::validate_bed_file(&self.microsatellite_bed)?;
        crate::cli::validate_vcf_file(&self.calls)?;

        Ok(())
    }
}
/* ================================================ */

/// Main Orchestrator for MSI candidate preprocessing workflow.
pub fn preprocess_ms_candidates(config: PreprocessMSIConfig) -> Result<()> {
    Ok(())
}
