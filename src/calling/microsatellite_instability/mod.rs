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
use crate::constants::{EPSILON, MIN_MSI_THRESHOLD};
use crate::errors::Error;
use crate::utils::bcf_utils::{
    is_phred_scaled_from_path, validate_events_exist, validate_samples_exist,
};

use dp_analysis::{AnalysisConfig, OutputRequirements};

/* ======== CLI CONFIGURATION ===================== */

/// Configuration for MSI calling pipeline.
///
/// Contains all parameters needed for per-sample MSI analysis,
/// populated from CLI arguments.
#[derive(Debug)]
pub struct MSIConfig {
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
    /// Fixed AF at which the full P(K=k) distribution is reported.
    /// Must be present in `af_thresholds` - enforced in `validate()`.
    pub distribution_af: f64,
    /// Fixed AF used for windowed heatmap analysis. Independent of `af_thresholds`.
    pub windowed_af: f64,
    /// Sliding window size (bp) for regional MSI heatmap analysis (default: 1,000,000).
    /// Only used if --plot-heatmap or --data-heatmap specified.
    pub sliding_window: u64,
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
    /// 1. Set is_phred based on events sepcification.
    pub fn set_defaults(&mut self) -> Result<()> {
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
    /// 2. Every entry of `af_thresholds` in [0.0, 1.0]. Required even though
    ///    `--af-thresholds` is a `hidden` CLI flag - hidden only suppresses it
    ///    from `--help`, it does not prevent a caller (e.g. an evaluation pipeline)
    ///    from overriding it with an out-of-range value.
    /// 3. `distribution_af` / `windowed_af` individually in [0.0, 1.0].
    /// 4. `distribution_af` present in `af_thresholds` (epsilon-compared) -
    ///    required because `calculate_msi_metrics` only populates `distribution`
    ///    when the `af-thresholds` contains `distribution_af` exactly.
    /// 5. Checks at least one output specified;
    /// 6. Validates the calls file format using existing validators from cli.rs;
    /// 7. Validates the calls file contains the specified sample.
    /// 8. Validates the calls file contains the specified events in the header.
    pub fn validate(&self) -> Result<()> {
        if self.msi_threshold <= MIN_MSI_THRESHOLD {
            return Err(Error::MsiConfigThresholdInvalid {
                threshold: self.msi_threshold,
            }
            .into());
        }

        for &af in &self.af_thresholds {
            if !(0.0..=1.0).contains(&af) {
                return Err(Error::MsiConfigAfThresholdInvalid { threshold: af }.into());
            }
        }

        if !(0.0..=1.0).contains(&self.distribution_af) {
            return Err(Error::MsiConfigAfThresholdInvalid {
                threshold: self.distribution_af,
            }
            .into());
        }

        if !(0.0..=1.0).contains(&self.windowed_af) {
            return Err(Error::MsiConfigAfThresholdInvalid {
                threshold: self.windowed_af,
            }
            .into());
        }

        let is_distribution_af = self
            .af_thresholds
            .iter()
            .any(|af| (af - self.distribution_af).abs() < EPSILON);
        if !is_distribution_af {
            return Err(Error::MsiConfigDistributionAfMissing {
                distribution_af: self.distribution_af,
                af_thresholds: self.af_thresholds.clone(),
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
            bcf::Reader::from_path(self.calls.as_path()).context("Failed to open VCF file")?;
        let header: HeaderView = vcf.header().clone();

        validate_vcf_file(&self.calls)?;
        validate_samples_exist(&header, &[self.sample.clone()])?;
        validate_events_exist(&header, &self.events)?;

        Ok(())
    }
}
/* ================================================ */

/// Orchestrates the MSI calling workflow: extract regions from a preprocessed +
/// called VCF/BCF, run AF-evolution analysis, then write requested outputs.
pub fn call_msi(config: MSIConfig) -> Result<()> {
    info!("----------------------------------------------");
    info!("Step 1: Data Extraction");
    info!("----------------------------------------------");

    let (global_results, window_results) = {
        let mut vcf =
            bcf::Reader::from_path(&config.calls).context("Failed to open calls VCF/BCF")?;
        let header = vcf.header().clone();

        let sample_idx = header.sample_id(config.sample.as_bytes()).ok_or_else(|| {
            Error::VcfSamplesNotFound {
                sample: config.sample.clone(),
            }
        })?;

        let (regions, stats) = extraction::extract_regions(
            &mut vcf,
            &config.sample,
            sample_idx,
            &config.events,
            config.is_phred,
        )?;
        stats.log_stats(&regions);

        info!("----------------------------------------------");
        info!("Step 2: AF Evolution Analysis");
        info!("----------------------------------------------");

        let output_req = OutputRequirements {
            needs_pseudotime: config.plot_pseudotime.is_some() || config.data_pseudotime.is_some(),
            needs_distribution: config.plot_distribution.is_some()
                || config.data_distribution.is_some(),
            needs_heatmap: config.plot_heatmap.is_some() || config.data_heatmap.is_some(),
        };

        let analysis_config = AnalysisConfig {
            total_regions: stats.total_ms_regions,
            sample: &config.sample,
            msi_high_threshold: config.msi_threshold,
            af_thresholds: &config.af_thresholds,
            num_threads: config.threads,
            window_size: config.sliding_window,
            distribution_af: config.distribution_af,
            windowed_af: config.windowed_af,
        };

        dp_analysis::run_af_evolution_analysis(&regions, analysis_config, output_req)?
    };

    info!("----------------------------------------------");
    info!("Step 3: Output Generation");
    info!("----------------------------------------------");

    if let Some(ref path) = config.data_distribution {
        output::write_distribution_data(
            &global_results,
            &config.sample,
            path,
            config.msi_threshold,
            config.distribution_af,
        )?;
    }

    if let Some(ref path) = config.plot_distribution {
        output::generate_distribution_plot_spec(
            &global_results,
            &config.sample,
            path,
            config.msi_threshold,
            config.distribution_af,
        )?;
    }

    if let Some(ref path) = config.data_pseudotime {
        output::write_pseudotime_data(&global_results, &config.sample, path, config.msi_threshold)?;
    }

    if let Some(ref path) = config.plot_pseudotime {
        output::generate_pseudotime_plot_spec(
            &global_results,
            &config.sample,
            path,
            config.msi_threshold,
        )?;
    }

    if let Some(ref path) = config.data_heatmap {
        output::write_heatmap_data(&window_results, &config.sample, path, config.msi_threshold)?;
    }

    if let Some(ref path) = config.plot_heatmap {
        output::generate_heatmap_plot_spec(
            &window_results,
            &config.sample,
            path,
            config.msi_threshold,
        )?;
    }

    info!("==============================================");
    info!("MSI calling complete");
    info!("==============================================");

    Ok(())
}
