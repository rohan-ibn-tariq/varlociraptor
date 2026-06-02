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

/// Standard INFO fields always propagated in varlociraptor preprocessing.
/// These fields are explicitly handled and must be omitted from aux_info.write()
/// to prevent double-writing.
///
/// TODO: Replace OMIT_AUX_INFO in calling/variants/mod.rs with this constant.
pub const STANDARD_OMIT_AUX_INFO: &[&[u8]] = &[
    b"MATEID",
    b"EVENT",
    b"SVLEN",
    b"SVTYPE",
    b"END",
];


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

/// MSI-specific fields to omit from aux_info.write() in addition to
/// STANDARD_OMIT_AUX_INFO. These are fields explicitly handled by
/// MSI preprocessing.
pub const MSI_OMIT_AUX_EXTRA: &[&[u8]] = &[
    b"REGION_ID",
    b"MSI_DUMMY",
];

/// Sliding window size (bp) for regional MSI heatmap analysis.
pub const DEFAULT_SLIDING_WINDOW_SIZE: u64 = 1_000_000;

/* ================================================ */
