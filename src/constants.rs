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
pub const DEFAULT_MSI_THRESHOLD: &str = "3.5";

/// Minimum allowed MSI threshold (percentage).
///
/// MSI thresholds must be positive values. Zero or negative thresholds
/// are not meaningful for MSI classification.
pub const MIN_MSI_THRESHOLD: f64 = 0.0;

/// INFO fields injected by MSI preprocessing output.
/// Added to omit set to prevent double-writing via --propagate-info-fields.
pub const MSI_OUTPUT_INFO_FIELDS: &[&str] = &[
    "REGION_ID",
    "MSI_DUMMY",
];

/// Sliding window size (bp) for regional MSI heatmap analysis.
pub const DEFAULT_SLIDING_WINDOW_SIZE: u64 = 1_000_000;

/* ================================================ */
