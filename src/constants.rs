//! constants.rs
//!
//! This module defines global constants used across varlociraptor.
//! Centralizing constants here promotes maintainability and consistency across the codebase.
//!
//! Note: Constants should be defined under proper sections (e.g. MSI: CONSTANTS) with clear documentation on their purpose and usage.
//!
//! Sections of Constants included are:
//! 1. MSI: Constants - Microsatellite Instability related constants.
//!

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

/* ================================================ */
