//! mod.rs
//!
//! Microsatellite instability (MSI) calling module.
//! Calculates MSI scores from enriched, called VCF/BCF for a single sample.
//!
//! **Note: This feature is experimental.**
//!
//! This module provides:
//! 1. @TODO: A brief overview of the MSI calling workflow, e.g.:
//! 2. Calculates MSI scores from enriched, called VCF/BCF for a single sample.

mod dp_analysis;
mod extraction;
mod output;

use std::path::PathBuf;

use anyhow::Result;
use log::info;

use crate::constants::MIN_MSI_THRESHOLD;
use crate::errors::Error;

/* ======== CLI CONFIGURATION ===================== */

/// Configuration for MSI calling pipeline.
///
/// Contains all parameters needed for the MSI calling workflow,
/// populated from CLI arguments.
#[derive(Debug)]
pub(crate) struct MSIConfig {
    /// Path to VCF/BCF file(should be sorted) with variant calls.
    pub calls: PathBuf,
    /// Number of threads (None = use rayon default).
    pub threads: Option<usize>,
    /// MSI-High classification threshold (percentage), default: 3.5.
    pub msi_threshold: f64,
    /// Sample name to process from VCF/BCF.
    pub sample: String,
    /// Event types to consider for MSI calling (e.g. "Somatic", "Germline").
    pub events: Vec<String>,
    /// Whether the probabilities in the VCF are PHRED-scaled.
    pub is_phred: bool,
    /// Allele frequency thresholds to consider for AF evolution analysis
    /// when generating pseudotime outputs. (default: [1.0,0.8,0.6,0.4,0.2,0.0])
    /// If no pseudotime outputs are requested, this will be set to [0.0] to optimize computation.
    /// Note: This field is populated during CLI parsing and currently not validated for non [0-1] values
    /// as the this field is hidden constant set at CLI level. So future changes to expose this field
    /// to users should include validation for this field.
    pub af_thresholds: Vec<f64>,
    /// Sliding window size (in base pairs) for MSI score calculation for heatmap. If None, default is used.
    pub sliding_window: Option<u64>,
    /// Output path for distribution plot (Vega-Lite JSON).
    pub plot_distribution: Option<PathBuf>,
    /// Output path for pseudotime plot (Vega-Lite JSON).
    pub plot_pseudotime: Option<PathBuf>,
    /// Output path for heatmap plot (Vega-Lite JSON).
    pub plot_heatmap: Option<PathBuf>,
    /// Output path for distribution data (TSV).
    pub data_distribution: Option<PathBuf>,
    /// Output path for pseudotime data (TSV).
    pub data_pseudotime: Option<PathBuf>,
    /// Output path for heatmap data (TSV).
    pub data_heatmap: Option<PathBuf>,
}

impl MSIConfig {
    /// Validate the MSI configuration.
    /// Responsibilities:
    /// 1. Checks for valid MSI threshold;
    /// 2. Checks at least one output specified;
    /// 3. Validates the calls file format using existing validators from cli.rs.
    pub fn validate(&self) -> Result<()> {
        if self.msi_threshold <= MIN_MSI_THRESHOLD {
            return Err(Error::MsiConfigThresholdInvalid {
                threshold: self.msi_threshold,
            }
            .into());
        }

        if self.plot_distribution.is_none()
            && self.plot_pseudotime.is_none()
            && self.plot_heatmap.is_none()
            && self.data_distribution.is_none()
            && self.data_pseudotime.is_none()
            && self.data_heatmap.is_none()
        {
            return Err(Error::MsiConfigOutputMissing.into());
        }

        // Use existing validators from cli.rs to validate file format.
        crate::cli::validate_vcf_file(&self.calls)?;

        Ok(())
    }
}
/* ================================================ */

/// Orchestrates the MSI calling workflow based on the provided configuration.
pub fn call_msi(config: MSIConfig) -> Result<()> {
    Ok(())
}
