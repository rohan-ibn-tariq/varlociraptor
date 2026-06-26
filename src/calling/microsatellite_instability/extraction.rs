//! extraction.rs
//!

use std::collections::HashMap;

use anyhow::Result;
use log::{debug, info, warn};
use rust_htslib::bcf::{self, header::HeaderView, Read};

use crate::constants::{MSI_REGION_ID_TAG, MSI_DUMMY_TAG}
use crate::errors::Error;
use crate::utils::bcf_utils::{
    get_chrom, get_events_probability, get_info_strings, get_sample_af, record_has_info_flag,
    record_has_info_string,
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
    /// Includes SNVs and non-perfect indels that passed through preprocessing unchanged.
    pub skipped_non_ms: usize,
    /// Unique MS regions encountered — used as the MSI score denominator.
    pub total_ms_regions: usize,
    /// Dummy indel records processed (MSI_DUMMY flag set).
    /// Each corresponds to a region where no perfect indel was observed in reads.
    pub dummy_records: usize,
    /// Alleles skipped — FORMAT:AF absent or missing for this sample.
    pub skipped_missing_af: usize,
    /// Alleles skipped — event probability field absent or missing.
    pub skipped_missing_prob: usize,
}

impl ExtractionStats {
    /// Log extraction statistics alongside region-level counts.
    pub fn log_stats(&self, regions: &[RegionSummary]) {
        let with_real = regions.iter().filter(|r| r.has_real_indel).count();

        info!("Extraction statistics:");
        info!("  Total records read:          {}", self.total_records);
        info!("  Non-MS records skipped:      {}", self.skipped_non_ms);
        info!("  Total MS regions:            {}", self.total_ms_regions);
        info!("  Regions with real indel:     {}", with_real);
        info!("  Regions with dummy indels:   {}", self.dummy_records);
        info!(
            "  Alleles skipped (missing AF):   {}",
            self.skipped_missing_af
        );
        info!(
            "  Alleles skipped (missing prob): {}",
            self.skipped_missing_prob
        );
    }
}

/// Stream a preprocessed+called VCF and extract per-region variant summaries.
///
/// Processes only records with `INFO/REGION_ID`, grouping variants by region
/// via a HashMap index into a Vec.
///
/// `has_real_indel` is set for regions where at least one non-dummy record
/// was encountered, regardless of whether prob or AF extraction succeeded.
/// This allows reporting a real-indel-only MSI score alongside the full score.
///
/// # Arguments
/// * `vcf`        - Reader for the called VCF
/// * `sample`     - Sample name (used in log messages)
/// * `sample_idx` - Index of this sample in FORMAT columns
/// * `events`     - Event names to combine (e.g., `["somatic"]`)
/// * `is_phred`   - Whether `INFO/PROB_*` values are PHRED-scaled
///
/// # Returns
/// `(Vec<RegionSummary>, ExtractionStats)` where Vec length equals
/// `stats.total_ms_regions`, which is the MSI score denominator.
///
/// # Example
/// assert!(extract_regions(&mut vcf, "sample1", 0, &["somatic".to_string()], true).is_ok());
pub(super) fn extract_regions(
    vcf: &mut bcf::Reader,
    sample: &str,
    sample_idx: usize,
    events: &[String],
    is_phred: bool,
) -> Result<(Vec<RegionSummary>, ExtractionStats)> {
    let header: HeaderView = vcf.header().clone();
    let mut regions: Vec<RegionSummary> = Vec::new();
    let mut region_index: HashMap<String, usize> = HashMap::new();
    let mut stats = ExtractionStats::default();
    let mut record = vcf.empty_record();

    loop {
        match vcf.read(&mut record) {
            None => break,
            Some(Err(e)) => {
                return Err(Error::VcfRecordReadFailed {
                    details: e.to_string(),
                }
                .into())
            }
            Some(Ok(())) => {}
        }

        // Preprocessing guarantees REGION_ID is only written when at least one
        // ALT matched a BED region - all-dot records are never written.
        // So presence of the field means at least one non-dot value exists.
        if !record_has_info_string(&record, MSI_REGION_ID_TAG) {
            stats.skipped_non_ms += 1;
            continue;
        }

        let is_dummy = record_has_info_flag(&record, MSI_DUMMY_TAG);
        if is_dummy {
            stats.dummy_records += 1;
        }

        let region_id =
            match get_info_strings(&record, MSI_REGION_ID_TAG).and_then(|v| v.into_iter().next()) {
                Some(id) => id,
                None => {
                    warn!(
                        "REGION_ID missing or empty at {}:{} — skipping",
                        get_chrom(&record, &header).unwrap_or_default(),
                        record.pos() + 1
                    );
                    continue;
                }
            };

        let region_idx = match region_index.get(&region_id) {
            Some(&idx) => idx,
            None => match RegionSummary::from_region_id(&region_id) {
                Ok(summary) => {
                    let idx = regions.len();
                    regions.push(summary);
                    region_index.insert(region_id, idx);
                    stats.total_ms_regions += 1;
                    idx
                }
                Err(e) => {
                    warn!("Skipping record with malformed REGION_ID: {}", e);
                    continue;
                }
            },
        };

        // Mark region as having a real (non-dummy) indel if applicable.
        // This is set regardless of whether prob/AF extraction succeeds -
        // the indel structurally existed in the reads.
        if !is_dummy {
            regions[region_idx].has_real_indel = true;
        }

        // Extract one Variant per ALT allele.
        let allele_count = record.allele_count() as usize;
        for alt_idx in 0..allele_count.saturating_sub(1) {
            // Combine event probabilities, P(at least one event)
            let prob_events =
                match get_events_probability(&record, &header, alt_idx, events, is_phred)? {
                    Some(p) => p,
                    None => {
                        stats.skipped_missing_prob += 1;
                        debug!(
                            "Event prob missing at {}:{} alt={} — skipping allele",
                            get_chrom(&record, &header).unwrap_or_default(),
                            record.pos() + 1,
                            alt_idx
                        );
                        continue;
                    }
                };

            // FORMAT:AF for this sample and ALT allele
            let af = match get_sample_af(&record, &header, sample_idx, alt_idx)? {
                Some(a) => a,
                None => {
                    stats.skipped_missing_af += 1;
                    debug!(
                        "FORMAT:AF missing for '{}' at {}:{} alt={} — skipping allele",
                        sample,
                        get_chrom(&record, &header).unwrap_or_default(),
                        record.pos() + 1,
                        alt_idx
                    );
                    continue;
                }
            };

            // prob_absent = 1 - P(events).
            regions[region_idx].variants.push(Variant {
                prob_absent: (1.0 - prob_events).max(0.0),
                af,
            });
        }
    }

    regions.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));

    Ok((regions, stats))
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
