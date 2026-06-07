//! extraction.rs
//!

use std::collections::HashMap;

use anyhow::Result;
use log::{debug, info, warn};
use rust_htslib::bcf::{self, header::HeaderView, Read};

use crate::errors::Error;
use crate::utils::bcf_utils::{
    get_chrom, get_events_probability, get_info_strings,
    /* get_sample_af,*/ record_has_info_flag, record_has_info_string,
};

/* ============ Data Structures =================== */

/// A single variant contributing to MSI analysis for one region.
///
/// Note: In multi-allelic records each ALT allele yields a separate Variant.
#[derive(Debug, Clone)]
pub(super) struct Variant {
    /// P(variant absent) = 1.0 - P(at least one specified event).
    pub prob_absent: f64,
    /// Allele frequency for the sample (FORMAT:AF).
    pub af: f64,
}

/// All variants observed within one microsatellite region.
///
/// Region is identified by its `REGION_ID` string ("chrom:start-end").
#[derive(Debug)]
pub(super) struct RegionSummary {
    /// Region identifier, format: "chrom:start-end"
    pub region_id: String,
    /// Chromosome — used for window assignment in heatmap analysis.
    pub chrom: String,
    /// Region start (0-based) — used for window assignment in heatmap analysis.
    pub start: u64,
    /// Variants in this region (real or dummy, ordered by encounter).
    pub variants: Vec<Variant>,
    /// True if at least one real (non-dummy) indel was observed in this region.
    /// False means the region only has a dummy indel (no indel observed in reads).
    pub has_real_indel: bool,
}
