//! constants.rs
//!
//! This module defines global constants used across varlociraptor.
//! Centralizing constants here promotes maintainability and consistency across the codebase.
//!
//! Note: Constants should be defined under proper sections (e.g. MSI: CONSTANTS) with clear documentation on their purpose and usage.
//!
//! Sections of Constants included are:
//! 1. Generic: Constants - General constants used across the codebase (e.g. default values, thresholds).
//! 2. MSI: Constants - Microsatellite Instability related constants.
//!

use std::collections::HashSet;

/* =============== GENERIC: CONSTANTS ============== */

// Standard INFO fields always propagated in varlociraptor preprocessing.
// These fields are explicitly handled and must be omitted from aux_info.write()
// to prevent double-writing.
//
// Use standard_omit_aux_info!() macro instead of a constant for this string.

/* =============== MSI: CONSTANTS ================== */

/// Default MSI-High threshold (percentage).
///
/// Variants with MSI score >= this threshold are classified as MSI-High.
/// This is the standard clinical cutoff used in MSI analysis.
pub(crate) const DEFAULT_MSI_THRESHOLD: &str = "3.5";

/// Minimum allowed MSI threshold (percentage).
///
/// MSI thresholds must be positive values. Zero or negative thresholds
/// are not meaningful for MSI classification.
pub(crate) const MIN_MSI_THRESHOLD: f64 = 0.0;

/// INFO fields injected by MSI preprocessing output.
/// Added to omit set to prevent double-writing via --propagate-info-fields.
pub(crate) const MSI_OUTPUT_INFO_FIELDS: &[&str] = &["REGION_ID", "MSI_DUMMY"];

/// Sliding window size (bp) for regional MSI heatmap analysis.
pub(crate) const DEFAULT_SLIDING_WINDOW_SIZE: u64 = 1_000_000;

/// INFO field tag for MSI region annotation.
/// Added to variants overlapping a microsatellite region.
pub(crate) const MSI_REGION_ID_TAG: &[u8] = b"REGION_ID";

/// INFO field tag for dummy deletion flag.
/// Set on synthetic records injected for MS regions without observed indels.
pub(crate) const MSI_DUMMY_TAG: &[u8] = b"MSI_DUMMY";

/// Full VCF header declaration for REGION_ID INFO field.
pub(crate) const MSI_REGION_ID_HEADER: &[u8] =
    br##"##INFO=<ID=REGION_ID,Number=A,Type=String,Description="BED region ID for the overlapping microsatellite locus">"##;

/// Full VCF header declaration for MSI_DUMMY INFO field.
pub(crate) const MSI_DUMMY_HEADER: &[u8] =
    br##"##INFO=<ID=MSI_DUMMY,Number=0,Type=Flag,Description="Dummy deletion injected for MS region with no observed indel">"##;

lazy_static! {
    /// INFO fields copied explicitly from input to output in MSI preprocessing.
    /// Derived from preprocess_msi_omit_aux_info!() macro.
    pub(crate) static ref PREPROCESS_MSI_COPY_FIELDS: Vec<&'static str> =
        preprocess_msi_omit_aux_info!()
            .split(", ")
            .collect();

    /// INFO fields omitted from aux_info.write() in MSI preprocessing.
    /// Combines standard propagated fields with MSI-specific output fields
    /// to prevent double-writing.
    pub(crate) static ref PREPROCESS_MSI_OMIT_AUX: HashSet<Vec<u8>> = {
        preprocess_msi_omit_aux_info!()
            .split(", ")
            .chain(MSI_OUTPUT_INFO_FIELDS.iter().copied())
            .map(|s| s.as_bytes().to_vec())
            .collect()
    };
}

/* ================================================ */
