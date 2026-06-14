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

use anyhow::{Context, Result};
use log::info;
use rust_htslib::bcf::{self, header::HeaderView, Read};

use crate::cli::validate_vcf_file;
use crate::constants::{DEFAULT_SLIDING_WINDOW_SIZE, MIN_MSI_THRESHOLD};
use crate::errors::Error;
use crate::utils::bcf_utils::{
    is_phred_scaled_from_path, validate_events_exist, validate_samples_exist,
};

/* ======== CLI CONFIGURATION ===================== */

/// Configuration for MSI calling pipeline.
///
/// Contains all parameters needed for per-sample MSI analysis,
/// populated from CLI arguments.
#[derive(Debug)]
pub(crate) struct MSIConfig {
    /// Path to VCF/BCF file with variant calls.
    pub calls: PathBuf,
    /// Sample name to analyze from VCF/BCF.
    pub sample: String,
    /// Event names to combine for MSI probability (e.g., ["somatic_tumor", "high_vaf"].
    pub events: Vec<String>,
    /// MSI-High classification threshold (percentage), default: 3.5.
    pub msi_threshold: f64,
    /// Allele frequency thresholds to consider for AF evolution analysis
    /// when generating pseudotime outputs. (default: [1.0,0.8,0.6,0.4,0.2,0.0])
    /// If no pseudotime outputs are requested, this will be set to [0.0] to optimize computation.
    /// Note: This field is populated during CLI parsing and currently not validated for non [0-1] values
    /// as the this field is hidden constant set at CLI level. So future changes to expose this field
    /// to users should include validation for this field.
    pub af_thresholds: Vec<f64>,
    /// Sliding window size (bp) for regional MSI heatmap analysis (default: 1,000,000).
    /// Only used if --plot-heatmap or --data-heatmap specified.
    pub sliding_window: Option<u64>,
    /// Number of threads (None = use rayon default).
    pub threads: Option<usize>,
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
    /// Whether the probabilities in the VCF are PHRED-scaled.
    pub is_phred: bool,
}

impl MSIConfig {
    /// Set default values for MSIConfig fields based on functional logic.
    /// Defaults Set:
    /// 1. Set sliding window to default if not specified and heatmap output(s) requested.
    /// 2. Set is_phred based on events sepcification.
    pub fn set_defaults(&mut self) -> Result<()> {
        if (self.plot_heatmap.is_some() || self.data_heatmap.is_some())
            && self.sliding_window.is_none()
        {
            self.sliding_window = Some(DEFAULT_SLIDING_WINDOW_SIZE);
        }

        self.is_phred = is_phred_scaled_from_path(&self.calls)?;
        info!(
            "  - Probabilities are {} scaled",
            if self.is_phred { "PHRED" } else { "linear" }
        );

        Ok(())
    }

    /// Validate the MSI configuration.
    /// Responsibilities:
    /// 1. Checks for valid MSI threshold;
    /// 2. Checks at least one output specified;
    /// 3. Validates the calls file format using existing validators from cli.rs;
    /// 4. Validates the calls file contains the specified sample.
    /// 5. Validates the calls file contains the specified events in the header.
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

        let vcf =
            bcf::Reader::from_path(&self.calls.as_path()).context("Failed to open VCF file")?;
        let header: HeaderView = vcf.header().clone();

        validate_vcf_file(&self.calls)?;
        validate_samples_exist(&header, &[self.sample.clone()])?;
        validate_events_exist(&header, &self.events)?;

        Ok(())
    }
}
/* ================================================ */

/// Orchestrates the MSI calling workflow based on the provided configuration.
pub fn call_msi(config: MSIConfig) -> Result<()> {
    // Extract
    //let (regions, stats) = extract_regions(...)?;
    // stats.log_stats(&regions);

    // NOTE: Move regions so it gets dropped before the next step, as it may be large
    // and we want to keep space efficient.

    // let (global, windows) = {
    //     let r = regions; // move into block
    //     run_af_evolution_analysis(&r, ...)?
    //     // r dropped
    // };

    Ok(())
}
