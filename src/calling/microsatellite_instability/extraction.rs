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

impl RegionSummary {
    /// Parse a REGION_ID string ("chrom:start-end") into a RegionSummary.
    ///
    /// Returns Err if the format is invalid so the caller can warn and skip.
    fn from_region_id(region_id: &str) -> Result<Self> {
        let (chrom, coords) =
            region_id
                .split_once(':')
                .ok_or_else(|| Error::MsiRegionIdMalformed {
                    region_id: region_id.to_string(),
                    details: "missing ':' separator".to_string(),
                })?;

        let (start_str, _end_str) =
            coords
                .split_once('-')
                .ok_or_else(|| Error::MsiRegionIdMalformed {
                    region_id: region_id.to_string(),
                    details: "missing '-' in coordinate part".to_string(),
                })?;

        let start: u64 = start_str.parse().map_err(|_| Error::MsiRegionIdMalformed {
            region_id: region_id.to_string(),
            details: format!("start '{}' is not a valid integer", start_str),
        })?;

        Ok(RegionSummary {
            region_id: region_id.to_string(),
            chrom: chrom.to_string(),
            start,
            variants: Vec::new(),
            has_real_indel: false,
        })
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf;
    use tempfile::NamedTempFile;

    use crate::utils::stats::test_constants::{TEST_EPSILON, TEST_EPSILON_LOOSE};

    /* ====== RegionSummary::from_region_id ======= */

    #[test]
    fn test_from_region_id_valid() {
        let r = RegionSummary::from_region_id("chr1:1000-2000").unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 1000);
        assert_eq!(r.region_id, "chr1:1000-2000");
        assert!(!r.has_real_indel);
        assert!(r.variants.is_empty());
    }

    #[test]
    fn test_from_region_id_missing_colon() {
        let err = RegionSummary::from_region_id("chr1_1000-2000").unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("chr1_1000-2000"));
        assert!(msg.contains("missing ':'"));
    }

    #[test]
    fn test_from_region_id_missing_dash() {
        let err = RegionSummary::from_region_id("chr1:10002000").unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("chr1:10002000"));
        assert!(msg.contains("missing '-'"));
    }

    #[test]
    fn test_from_region_id_non_numeric_start() {
        let err = RegionSummary::from_region_id("chr1:abc-2000").unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("chr1:abc-2000"));
        assert!(msg.contains("not a valid integer"));
    }
}
