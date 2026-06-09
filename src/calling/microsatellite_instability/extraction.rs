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
    /// True if at least one non-dummy indel record was encountered,
    /// regardless of whether prob or AF extraction succeeded.
    /// Used to report the "real indel only" MSI score alongside
    /// the full (real + dummy) score.
    pub has_real_indel: bool,
}

/// Statistics collected during extraction.
#[derive(Debug, Default)]
pub(super) struct ExtractionStats {
    /// Total VCF records read (MS and non-MS).
    pub total_records: usize,
    /// Records skipped — no REGION_ID, not MS-relevant.
    pub skipped_non_ms: usize,
    /// Unique MS regions encountered.
    pub total_ms_regions: usize,
    /// Dummy records processed (MSI_DUMMY flag set).
    pub dummy_records: usize,
    /// Alleles skipped — FORMAT:AF absent or missing for this sample.
    pub skipped_missing_af: usize,
    /// Alleles skipped — event probability field absent or missing.
    pub skipped_missing_prob: usize,
}
