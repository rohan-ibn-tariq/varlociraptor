//! dp_analysis.rs
//!
//! Dynamic programming analysis for MSI estimation
//! Implements AF evolution analysis using dynamic programming.
//! Computes MSI probability distributions across allele frequency thresholds.
//!
//! # Overview
//! 1. Core Dp Algorithm
//! 2. Filtering Variants Per Region by Sample and AF
//! 3. Calculating MSI Metrics
//! 4. Af Evolution Function: Generates results for all samples and AF thresholds based on user input.

use std::collections::HashMap;
use std::sync::Mutex;

use anyhow::Result;
use log::info;
use rayon::prelude::*;
use rust_decimal::prelude::*;
use rust_decimal::Decimal;
use serde::Serialize;

use super::extraction::{RegionSummary, Variant};
use crate::utils::genomics::classify_msi_status;
use crate::utils::stats::calculate_percentage_exact;

/* ============ Data Structures =================== */

/// Region probability for DP computation:
///
/// Represents the probability that a microsatellite region is stable,
/// calculated as the product of all variant absence probabilities.
///
/// p_stable = P(all variants absent) = Π(prob_absent)
/// p_unstable = 1 - p_stable (computed in DP algorithm)
#[derive(Debug, Clone)]
struct RegionProbability {
    p_stable: f64,
}

/// DP result for one k value:
///
/// Represents the probability that exactly k regions are unstable,
/// along with the corresponding MSI score.
#[derive(Debug, Clone, Serialize)]
pub(super) struct DpResult {
    pub k: usize,
    pub msi_score: f64,
    pub probability: f64,
}

/// Result for one sample at one AF threshold:
///
/// Contains the MSI analysis results including MAP estimate,
/// uncertainty bounds, and optionally the full probability distribution.
#[derive(Debug, Clone, Serialize)]
pub(super) struct AfEvolutionResult {
    pub sample: String,
    pub af_threshold: f64,
    // Only computed if pseudotime output requested:
    #[serde(skip_serializing_if = "Option::is_none")]
    pub k_map: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub msi_score_map: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub regions_with_variants: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub msi_status: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub uncertainty_lower: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub uncertainty_upper: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub map_std_dev: Option<f64>,
    // Only computed if distribution output requested AND AF=0.0:
    #[serde(skip_serializing_if = "Option::is_none")]
    pub distribution: Option<Vec<DpResult>>,
}

/// MSI score for one genomic window at one AF threshold.
///
/// Forms one cell in the heatmap:
/// - X-axis: genomic position (`window_start`)
/// - Y-axis: AF threshold (`af_threshold`)
/// - Color:  MSI score (`msi_score`)
///
/// Produced by `run_windowed_analysis` when `--sliding-window` is specified.
#[derive(Debug)]
pub(super) struct WindowResult {
    /// Chromosome name.
    pub chrom: String,
    /// Window start position (0-based, inclusive).
    pub window_start: u64,
    /// Window end position (0-based, exclusive).
    pub window_end: u64,
    /// AF threshold (Y-axis in heatmap).
    pub af_threshold: f64,
    /// MSI score for this window at this AF threshold.
    pub msi_score: f64,
    /// Total MS regions in this window - denominator for MSI score.
    pub regions_in_window: usize,
}

/// Filtered regions view
///
/// Uses a flat vector of variant references with region boundary indices
/// to avoid allocating thousands of small Vecs (~1-3 per region).
#[derive(Debug)]
struct FilteredRegions<'a> {
    /// Flat vector of all variant references
    /// Example: [&V1, &V2, &V3, &V4, &V5]
    variants: Vec<&'a Variant>,

    /// Region boundaries (start indices)
    /// Example: [0, 2, 3] means:
    ///   Region 0: variants[0..2] = [&V1, &V2]
    ///   Region 1: variants[2..3] = [&V3]
    ///   Region 2: variants[3..5] = [&V4, &V5]
    region_starts: Vec<usize>,
}

impl<'a> FilteredRegions<'a> {
    /// Get variants for a specific region
    fn get_region(&self, region_idx: usize) -> &[&'a Variant] {
        let start = self.region_starts[region_idx];
        let end = self
            .region_starts
            .get(region_idx + 1)
            .copied()
            .unwrap_or(self.variants.len());
        &self.variants[start..end]
    }

    /// Number of regions
    fn len(&self) -> usize {
        self.region_starts.len()
    }
}

/// Output requirements for optimizing calculations(conditional generation of metrics):
///
/// Determines which expensive computations are needed based on
/// which output files the user requested.
/// - needs_pseudotime: For pseudotime data
/// - needs_distribution: For distribution data
/// - needs_heatmap: For heatmap data
#[derive(Debug, Clone, Copy)]
pub(super) struct OutputRequirements {
    /// Whether to compute uncertainty bounds (std dev, lower/upper)
    pub needs_pseudotime: bool,
    /// Whether to compute full probability distribution
    pub needs_distribution: bool,
    /// Whether to compute windowed heatmap analysis
    pub needs_heatmap: bool,
}

/* ================================================ */

/* =========== Core DP Algorithm ================== */
/// Execute dynamic programming algorithm for MSI variant probability distribution.
///
/// Implements matrix-based approach where:
/// - Matrix has n columns (regions 0 to n-1) and n+1 rows (unstable counts 0 to n)
/// - Our implementation stores row vectors of length n+1
/// - Rows represent possible unstable region counts (0 to n)
/// - Each cell contains P(exactly i unstable regions using first k regions)
///
/// # Algorithm
///
/// - **Initialize column 0**: `prev_col[0] = p_stable_0`, `prev_col[1] = p_unstable_0`
/// - **Recurrence**: `curr_col[i] = prev_col[i] × p_stable_k + prev_col[i-1] × p_unstable_k`
///
/// # Arguments
///
/// * `region_probs` - Regions with instability probabilities where:
///   - `p_stable` = P(all variants absent) = Π(prob_absent)
///   - `p_unstable` = 1 - `p_stable`
///
/// # Returns
///
/// Probability distribution `[P(0 unstable), P(1 unstable), ..., P(n unstable)]`
/// where sum equals 1.0 and each `P(i)` represents probability of exactly
/// i regions being unstable.
///
/// # Note
///
/// - Uses space-optimized approach storing only previous column instead of full matrix
/// - `p_unstable + p_stable` should equal 1.0
/// - Floating-point precision: PHRED probability conversion has ~0.005 tolerance for
///   floating-point precision
///
/// # Example
///
/// ```text
/// Input: 2 regions with p_unstable = [0.3, 0.4]
///
/// Step 1 (Region 0):
///   prev_col[0] = 0.7  // P(0 unstable) = p_stable_0
///   prev_col[1] = 0.3  // P(1 unstable) = p_unstable_0
///
/// Step 2 (Region 1):
///   p_stable_1 = 0.6, p_unstable_1 = 0.4
///   
///   curr_col[0] = prev_col[0] × 0.6 = 0.7 × 0.6 = 0.42
///   curr_col[1] = prev_col[1] × 0.6 + prev_col[0] × 0.4 = 0.3 × 0.6 + 0.7 × 0.4 = 0.46
///   curr_col[2] = prev_col[1] × 0.4 = 0.3 × 0.4 = 0.12
///
/// Result: [0.42, 0.46, 0.12]
/// Meaning: P(0 unstable)=42%, P(1 unstable)=46%, P(2 unstable)=12%
/// ```
fn run_msi_dp(region_probs: &[RegionProbability]) -> Vec<f64> {
    let n = region_probs.len();

    // Base case: no regions means P(0 unstable) = 100%
    if n == 0 {
        return vec![1.0];
    }

    // Initialize probability vector: [P(0), P(1), P(2), ..., P(n)]
    let mut prev_col = vec![0.0; n + 1];

    // Column 0: Initialize with first region probabilities
    let p_stable_0 = region_probs[0].p_stable;
    let p_unstable_0 = 1.0 - p_stable_0;

    prev_col[0] = p_stable_0; // P(0 unstable) = first region stable
    prev_col[1] = p_unstable_0; // P(1 unstable) = first region unstable

    // Process remaining regions using recurrence relation
    #[allow(clippy::needless_range_loop)]
    for k in 1..n {
        let p_stable_k = region_probs[k].p_stable;
        let p_unstable_k = 1.0 - p_stable_k;

        let mut curr_col = vec![0.0; n + 1];

        // Base case: P(0 unstable) = previous P(0 unstable) × current region stable
        curr_col[0] = prev_col[0] * p_stable_k;

        // Recurrence: P(exactly i unstable) has two paths:
        // Path 1: Had i unstable, current region stable → still i unstable
        // Path 2: Had i-1 unstable, current region unstable → now i unstable
        for i in 1..=k + 1 {
            // Up to k+1 total unstable regions possible
            curr_col[i] = prev_col[i] * p_stable_k      // Path 1: Stay at i
                        + prev_col[i - 1] * p_unstable_k; // Path 2: Increment to i
        }

        prev_col = curr_col; // Update for next iteration
    }

    prev_col // Final distribution: [P(0), P(1), ..., P(n)]
}

/* ================================================ */

/* =========== Parallelization Setup ============== */

/// Configure the rayon global thread pool.
///
/// Called once at the start of `run_af_evolution_analysis` before
/// any parallel work begins. Both `run_global_analysis` and
/// `run_windowed_analysis` share this pool automatically.
///
/// Silently falls back to the existing pool if already configured.
fn setup_thread_pool(num_threads: Option<usize>) {
    if let Some(threads) = num_threads {
        // Note: In CLI usage, build_global() should always succeed on first call.
        // Error handling included for future library usage or testing scenarios.
        match rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
        {
            Ok(_) => {
                info!("Using {} threads (CLI specified)", threads);
            }
            Err(_) => {
                info!(
                    "Using {} threads (global pool already configured, requested {} ignored)",
                    rayon::current_num_threads(),
                    threads
                );
            }
        }
    } else {
        info!(
            "Using {} threads (rayon default)",
            rayon::current_num_threads()
        );
    }
}

/* ================================================ */

/* === FILTERING REGIONS VARIANTS BY AF============ */

/// Filter regions by AF threshold for a specific sample
///
/// Creates a flat filtered view avoiding Vec-per-region allocation overhead.
/// Optimized for MSI analysis where most regions have ~ 1-3 variants.
///
/// # Arguments
/// * `regions` - All regions with variants (lifetime 'a)
/// * `af_threshold` - Minimum AF
///
/// # Returns
/// Flat filtered view with region boundary tracking.
fn filter_regions_by_af<'a>(
    regions: &'a [RegionSummary],
    af_threshold: f64,
) -> FilteredRegions<'a> {
    let mut all_variants = Vec::new();
    let mut region_starts = Vec::new();

    for region in regions {
        let region_start = all_variants.len();
        let mut found_any = false;

        for variant in &region.variants {
            if variant.af >= af_threshold {
                all_variants.push(variant);
                found_any = true;
            }
        }

        if found_any {
            region_starts.push(region_start);
        }
    }

    FilteredRegions {
        variants: all_variants,
        region_starts,
    }
}

/* ================================================ */

/* ============= MSI Metrics ====================== */

/// Compute DP probability distribution from pre-filtered regions.
///
/// Converts filtered regions into region stability probabilities
/// then runs the DP algorithm to get P(k regions unstable).
///
/// # Arguments
/// * `filtered` - Pre-filtered region view from `filter_regions_by_af`
///
/// # Returns
/// Probability distribution `[P(0 unstable), ..., P(n unstable)]`
/// summing to 1.0. Returns `[1.0]` if no regions passed filtering.
fn run_dp_for_regions(filtered: &FilteredRegions) -> Vec<f64> {
    let region_probs: Vec<RegionProbability> = (0..filtered.len())
        .map(|i| {
            let p_all_absent: f64 = filtered
                .get_region(i)
                .iter()
                .map(|v| v.prob_absent)
                .product();
            RegionProbability {
                p_stable: p_all_absent,
            }
        })
        .collect();

    if region_probs.is_empty() {
        vec![1.0]
    } else {
        run_msi_dp(&region_probs)
    }
}

/// Find MAP estimate and MSI score from a DP distribution.
///
/// Finds k that maximizes P(K=k) - the maximum a posteriori estimate
/// of unstable region count. Ties broken by taking the highest k.
///
/// # Arguments
/// * `dist`          - DP probability distribution from `run_dp_for_regions`
/// * `total_regions` - Denominator for MSI score (total MS regions)
///
/// # Returns
/// `(k_map, msi_score)` where `msi_score = k_map / total_regions × 100`
fn find_map_estimate(dist: &[f64], total_regions: usize) -> (usize, f64) {
    let k_map = dist
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(k, _)| k)
        .unwrap_or(0);

    let msi_score = calculate_percentage_exact(k_map, total_regions);
    (k_map, msi_score)
}

/// Calculate MSI metrics for filtered regions
///
/// Computes region instability probabilities, runs DP, calculates MAP estimate,
/// and computes uncertainty bounds if requested.
///
/// # Steps
/// 1. **Region Probabilities**: P(region unstable) = 1 - Π(p_absent)
/// 2. **DP**: Compute probability distribution P(K = k)
/// 3. **MAP**: Find k that maximizes P(K = k)
/// 4. **MSI Score**: (k_map / total_regions) × 100
/// 5. **Uncertainty** (if requested): Standard deviation bounds
/// 6. **Distribution** (if requested and AF=0.0): Full P(K = k) array
///
/// # Arguments
/// * `filtered` - Filtered view of variants passing AF threshold
/// * `total_regions` - Total number of MS regions
/// * `msi_high_threshold` - Threshold for MSI-High classification (e.g., 3.5%)
/// * `af_threshold` - AF threshold used for filtering
/// * `sample` - The sample being processed
/// * `output_req` - Which outputs are requested (controls what to compute)
///
/// # Returns
/// Complete MSI analysis result for this sample/AF combination
fn calculate_msi_metrics(
    filtered: &FilteredRegions,
    total_regions: usize,
    msi_high_threshold: f64,
    af_threshold: f64,
    sample: String,
    output_req: OutputRequirements,
) -> AfEvolutionResult {
    // Step 1: Compute region instability probabilities
    // For each region, calculate P(at least one variant present)
    let region_probs: Vec<RegionProbability> = (0..filtered.len())
        .map(|i| {
            let region_variants = filtered.get_region(i);

            // P(all variants absent) = Π(p_absent)
            let p_all_absent: f64 = region_variants.iter().map(|v| v.prob_absent).product();

            RegionProbability {
                p_stable: p_all_absent,
            }
        })
        .collect();

    // Step 2: Run DP to get probability distribution
    let distribution_raw = if !region_probs.is_empty() {
        run_msi_dp(&region_probs)
    } else {
        vec![1.0] // No regions → P(0 unstable) = 1.0
    };

    // Step 3: Compute k_map, msi_score_map, regions_with_variants, msi_status ONLY if pseudotime needed
    let (k_map, msi_score_map, msi_status, regions_with_variants) =
        if output_req.needs_pseudotime && total_regions > 0 {
            // Note: partial_cmp is safe here because upstream validation ensures no NaN values.
            let k_map = distribution_raw
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .map(|(k, _)| k)
                .unwrap_or(0);
            let regions_with_variants = region_probs.len();
            let msi_score_map = calculate_percentage_exact(k_map, total_regions);
            let msi_status = classify_msi_status(msi_score_map, msi_high_threshold).to_string();

            (
                Some(k_map),
                Some(msi_score_map),
                Some(msi_status),
                Some(regions_with_variants),
            )
        } else {
            (None, None, None, None)
        };

    // Step 4: Calculate uncertainty using exact Decimal arithmetic
    let (uncertainty_lower, uncertainty_upper, map_std_dev) =
        if output_req.needs_pseudotime && total_regions > 0 {
            let k_map_decimal = Decimal::from(k_map.unwrap());

            // Calculate variance: Var(K) = Σ[(k - k_map)² × P(k)]
            let variance_decimal: Decimal = distribution_raw
                .iter()
                .enumerate()
                .map(|(k, &prob)| {
                    let k_decimal = Decimal::from(k);
                    let diff = k_decimal - k_map_decimal;
                    let diff_squared = diff * diff;
                    let prob_decimal = Decimal::from_f64_retain(prob).unwrap_or(Decimal::from(0));
                    diff_squared * prob_decimal
                })
                .sum();

            // Calculate standard deviation
            let std_dev_decimal = variance_decimal.sqrt().unwrap_or(Decimal::from(0));

            // Get confidence bounds
            let total_decimal = Decimal::from(total_regions);
            let hundred = Decimal::from(100);
            let zero = Decimal::from(0);

            let lower_k_decimal = (k_map_decimal - std_dev_decimal).max(zero);
            let upper_k_decimal = (k_map_decimal + std_dev_decimal).min(total_decimal);
            let lower_percentage_decimal = (lower_k_decimal / total_decimal) * hundred;
            let upper_percentage_decimal = (upper_k_decimal / total_decimal) * hundred;
            let lower = lower_percentage_decimal
                .max(zero)
                .min(hundred)
                .to_f64()
                .unwrap_or(0.0);
            let upper = upper_percentage_decimal
                .max(zero)
                .min(hundred)
                .to_f64()
                .unwrap_or(0.0);
            let std_dev_f64 = std_dev_decimal.to_f64().unwrap_or(0.0);

            (Some(lower), Some(upper), Some(std_dev_f64))
        } else {
            (None, None, None)
        };

    // Step 5: Create full distribution (only if distribution output requested AND AF=0.0)
    let distribution = if output_req.needs_distribution && af_threshold == 0.0 {
        Some(
            distribution_raw
                .iter()
                .enumerate()
                .map(|(k, &probability)| {
                    let msi_score = calculate_percentage_exact(k, total_regions);
                    DpResult {
                        k,
                        msi_score,
                        probability,
                    }
                })
                .collect(),
        )
    } else {
        None
    };

    AfEvolutionResult {
        sample,
        af_threshold,
        k_map,
        msi_score_map,
        regions_with_variants,
        msi_status,
        uncertainty_lower,
        uncertainty_upper,
        map_std_dev,
        distribution,
    }
}

/* ================================================ */

/* ============= AF ANALYSIS UTILITIES ============ */

/// Get a zero-copy sub-slice of regions falling within a genomic window.
///
/// # Precondition
/// `regions` must be sorted by `(chrom, start)`. This should be guaranteed by
/// `extract_regions`, otherwise results would be incorrect.
///
/// Regions must be sorted by (chrom, start) - guaranteed by extraction order.
///
/// # Arguments
/// * `regions`      - All extracted regions sorted by `(chrom, start)`
/// * `chrom`        - Chromosome to filter for
/// * `window_start` - Window start position (0-based, inclusive)
/// * `window_end`   - Window end position (0-based, exclusive)
///
/// # Returns
/// Sub-slice of `regions` whose chrom matches and start falls in
/// `[window_start, window_end)`. Returns empty slice if no match.
///
/// # Example
/// assert_eq!(get_window_slice(&regions, "chr1", 1000, 2000), &[RegionSummary { chrom: "chr1", start: 1500, ... }]);
fn get_window_slice<'a>(
    regions: &'a [RegionSummary],
    chrom: &str,
    window_start: u64,
    window_end: u64,
) -> &'a [RegionSummary] {
    // Find first index where chrom matches and start >= window_start
    let start_idx = regions
        .iter()
        .position(|r| r.chrom == chrom && r.start >= window_start);

    let start_idx = match start_idx {
        Some(i) => i,
        None => return &[],
    };

    // Find end index where chrom still matches and start < window_end
    let end_idx = regions[start_idx..]
        .iter()
        .position(|r| r.chrom != chrom || r.start >= window_end)
        .map(|i| start_idx + i)
        .unwrap_or(regions.len());

    &regions[start_idx..end_idx]
}

/* == DP ANALYSIS CORE: RUN AF EVOLUTION ========== */

/// Run AF evolution analysis across AF thresholds for one sample.
///
/// Performs parallel MSI analysis for each AF threshold,
/// computing probability distributions and MSI scores.
///
/// # Arguments
/// * `regions` - Regions with variants from extraction.
/// * `total_regions` - Total MS valid regions in BED (for MSI score denominator)
/// * `sample`        - Sample name (for result labeling)
/// * `msi_high_threshold` - Threshold for MSI-High classification (percentage)
/// * `af_thresholds` - AF thresholds to analyze (typically [1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.02, 0.0] or [0.0])
/// * `output_req` - Which outputs requested (controls computation)
/// * `num_threads` - Thread count (None = use rayon default)
///
/// # Returns
/// HashMap keyed by AF threshold string (e.g. "0.00", "0.50") -> MSI result.
pub(super) fn run_af_evolution_analysis(
    regions: &[RegionSummary],
    total_regions: usize,
    sample: &str,
    msi_high_threshold: f64,
    af_thresholds: &[f64],
    output_req: OutputRequirements,
    num_threads: Option<usize>,
) -> Result<HashMap<String, AfEvolutionResult>> {
    info!("Sample: {:?}", sample);
    info!("AF thresholds: {:?}", af_thresholds);
    info!("MSI-High threshold: {}%", msi_high_threshold);
    info!("Total regions (BED): {}", total_regions);
    info!("Regions with variants: {}", regions.len());
    info!("Output requirements:");
    info!("    - Uncertainty metrics: {}", output_req.needs_pseudotime);
    info!("    - Full distribution: {}", output_req.needs_distribution);

    // Step 1: Configure rayon thread pool
    setup_thread_pool(num_threads); // once shared by both analysis functions

    // Create all (AF) based work items for parallel processing
    let work_items: Vec<f64> = af_thresholds.to_vec();

    info!("Total parallel tasks: {}", work_items.len());

    let results = Mutex::new(HashMap::new());

    work_items.par_iter().for_each(|af_threshold| {
        let filtered = filter_regions_by_af(regions, *af_threshold);

        let result = calculate_msi_metrics(
            &filtered,
            total_regions,
            msi_high_threshold,
            *af_threshold,
            sample.to_string(),
            output_req,
        );

        results
            .lock()
            .unwrap()
            .insert(format!("{:.2}", af_threshold), result);
    });

    let all_results = results.into_inner().unwrap();

    Ok(all_results)
}

/* ================================================ */

#[cfg(test)]
mod tests {
    use super::*;

    use crate::utils::stats::test_constants::TEST_EPSILON;

    /// Create a test Variant with given absence probability and allele frequency.
    fn make_variant(prob: f64, af: f64) -> Variant {
        Variant {
            prob_absent: prob,
            af,
        }
    }

    /// Create a test RegionSummary with full control over all fields.
    fn make_region(
        region_id: &str,
        chrom: &str,
        start: u64,
        variants: Vec<Variant>,
        has_real_indel: bool,
    ) -> RegionSummary {
        RegionSummary {
            region_id: region_id.to_string(),
            chrom: chrom.to_string(),
            start,
            variants,
            has_real_indel,
        }
    }

    /// Create a test RegionSummary with default location fields.
    /// Use when the test only cares about variants, not genomic position.
    fn make_region_simple(variants: Vec<Variant>) -> RegionSummary {
        make_region("chr1:0-100", "chr1", 0, variants, true)
    }

    /* ============ DP Core Tests ==================== */

    #[test]
    fn test_run_msi_dp_empty() {
        let probs = vec![];
        let dist = run_msi_dp(&probs);

        assert_eq!(dist.len(), 1);
        assert!((dist[0] - 1.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_msi_dp_single_region() {
        let probs = vec![RegionProbability { p_stable: 0.7 }];
        let dist = run_msi_dp(&probs);

        assert_eq!(dist.len(), 2);
        assert!((dist[0] - 0.7).abs() < TEST_EPSILON); // P(0 unstable) = 0.7
        assert!((dist[1] - 0.3).abs() < TEST_EPSILON); // P(1 unstable) = 0.3
    }

    #[test]
    fn test_run_msi_dp_two_regions() {
        let probs = vec![
            RegionProbability { p_stable: 0.7 },
            RegionProbability { p_stable: 0.6 },
        ];
        let dist = run_msi_dp(&probs);

        assert_eq!(dist.len(), 3);
        // P(0) = 0.7 × 0.6 = 0.42
        // P(1) = 0.7 × 0.4 + 0.3 × 0.6 = 0.46
        // P(2) = 0.3 × 0.4 = 0.12
        assert!((dist[0] - 0.42).abs() < TEST_EPSILON);
        assert!((dist[1] - 0.46).abs() < TEST_EPSILON);
        assert!((dist[2] - 0.12).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_msi_dp_distribution_sums_to_one() {
        let probs = vec![
            RegionProbability { p_stable: 0.2 },
            RegionProbability { p_stable: 0.5 },
            RegionProbability { p_stable: 0.8 },
        ];
        let dist = run_msi_dp(&probs);

        let sum: f64 = dist.iter().sum();
        assert!((sum - 1.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_msi_dp_all_stable() {
        let probs = vec![
            RegionProbability { p_stable: 1.0 },
            RegionProbability { p_stable: 1.0 },
        ];
        let dist = run_msi_dp(&probs);
        assert!((dist[0] - 1.0).abs() < TEST_EPSILON);
        assert!((dist[1] - 0.0).abs() < TEST_EPSILON);
        assert!((dist[2] - 0.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_msi_dp_all_unstable() {
        let probs = vec![
            RegionProbability { p_stable: 0.0 },
            RegionProbability { p_stable: 0.0 },
        ];
        let dist = run_msi_dp(&probs);
        assert!((dist[0] - 0.0).abs() < TEST_EPSILON);
        assert!((dist[1] - 0.0).abs() < TEST_EPSILON);
        assert!((dist[2] - 1.0).abs() < TEST_EPSILON);
    }

    /* ============ FilteredRegions Tests ============ */

    #[test]
    fn test_filtered_regions_get_region() {
        let v1 = make_variant(0.01, 0.8);
        let v2 = make_variant(0.02, 0.9);
        let v3 = make_variant(0.03, 0.7);

        let filtered = FilteredRegions {
            variants: vec![&v1, &v2, &v3],
            region_starts: vec![0, 2],
        };

        // Region 0: variants[0..2]
        let r0 = filtered.get_region(0);
        assert_eq!(r0.len(), 2);
        assert!((r0[0].prob_absent - 0.01).abs() < TEST_EPSILON);
        assert!((r0[1].prob_absent - 0.02).abs() < TEST_EPSILON);

        // Region 1: variants[2..3]
        let r1 = filtered.get_region(1);
        assert_eq!(r1.len(), 1);
        assert!((r1[0].prob_absent - 0.03).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_filtered_regions_len() {
        let v1 = make_variant(0.01, 0.8);

        let filtered = FilteredRegions {
            variants: vec![&v1],
            region_starts: vec![0],
        };

        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered.get_region(0).len(), 1);
    }

    #[test]
    fn test_filtered_regions_empty() {
        let filtered = FilteredRegions {
            variants: vec![],
            region_starts: vec![],
        };
        assert_eq!(filtered.len(), 0);
    }

    /* ============ AF Filtering Tests =============== */

    #[test]
    fn test_filter_by_af_basic() {
        let regions = vec![make_region_simple(vec![
            make_variant(0.01, 0.8), // Pass
            make_variant(0.02, 0.4), // Fail
            make_variant(0.03, 0.7), // Pass
        ])];

        let filtered = filter_regions_by_af(&regions, 0.6);

        assert_eq!(filtered.len(), 1);
        let region_1 = filtered.get_region(0);
        assert_eq!(region_1.len(), 2); // V1 and V3
        assert!((region_1[0].prob_absent - 0.01).abs() < TEST_EPSILON);
        assert!((region_1[1].prob_absent - 0.03).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_filter_af_zero_includes_all() {
        let regions = vec![make_region_simple(vec![
            make_variant(0.01, 0.0),
            make_variant(0.02, 0.5),
            make_variant(0.03, 1.0),
        ])];

        let filtered = filter_regions_by_af(&regions, 0.0);

        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered.get_region(0).len(), 3);
    }

    #[test]
    fn test_filter_multiple_regions_some_empty() {
        let regions = vec![
            make_region(
                "chr1:100-130",
                "chr1",
                100,
                vec![make_variant(0.01, 0.8)],
                true,
            ),
            make_region(
                "chr1:200-230",
                "chr1",
                200,
                vec![make_variant(0.02, 0.3)],
                true,
            ), // Filtered out
            make_region(
                "chr1:300-330",
                "chr1",
                300,
                vec![make_variant(0.03, 0.9)],
                true,
            ),
        ];

        let filtered = filter_regions_by_af(&regions, 0.5);
        assert_eq!(filtered.len(), 2); // Only first and third
        assert!(filtered.get_region(0)[0].af - 0.8 < TEST_EPSILON);
        assert!(filtered.get_region(1)[0].af - 0.9 < TEST_EPSILON);
    }

    /* ======= MSI Metrics tests ===================== */
    /* ======= run_dp_for_regions tests ============== */

    #[test]
    fn test_run_dp_for_regions_empty() {
        let filtered = FilteredRegions {
            variants: vec![],
            region_starts: vec![],
        };
        let dist = run_dp_for_regions(&filtered);
        assert_eq!(dist, vec![1.0]);
    }

    #[test]
    fn test_run_dp_for_regions_single_region() {
        // prob_absent=0.7, p_stable=0.7
        let v = make_variant(0.7, 0.8);
        let filtered = FilteredRegions {
            variants: vec![&v],
            region_starts: vec![0],
        };
        let dist = run_dp_for_regions(&filtered);
        assert_eq!(dist.len(), 2);
        assert!((dist[0] - 0.7).abs() < TEST_EPSILON); // P(0 unstable)
        assert!((dist[1] - 0.3).abs() < TEST_EPSILON); // P(1 unstable)
    }

    #[test]
    fn test_run_dp_for_regions_sums_to_one() {
        let v1 = make_variant(0.3, 0.8);
        let v2 = make_variant(0.5, 0.9);
        let v3 = make_variant(0.7, 0.7);
        let filtered = FilteredRegions {
            variants: vec![&v1, &v2, &v3],
            region_starts: vec![0, 1, 2],
        };
        let dist = run_dp_for_regions(&filtered);
        let sum: f64 = dist.iter().sum();
        assert!((sum - 1.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_run_dp_for_regions_multiple_variants_per_region() {
        // Region with two variants: p_stable = 0.3 × 0.5 = 0.15
        let v1 = make_variant(0.3, 0.8);
        let v2 = make_variant(0.5, 0.9);
        let filtered = FilteredRegions {
            variants: vec![&v1, &v2],
            region_starts: vec![0], // one region with two variants
        };
        let dist = run_dp_for_regions(&filtered);
        assert_eq!(dist.len(), 2);
        assert!((dist[0] - 0.15).abs() < TEST_EPSILON); // P(0) = 0.15
        assert!((dist[1] - 0.85).abs() < TEST_EPSILON); // P(1) = 0.85
    }

    /* ======= find_map_estimate tests =============== */

    #[test]
    fn test_find_map_estimate_basic() {
        // P(0)=0.42, P(1)=0.46, P(2)=0.12, k_map=1
        let dist = vec![0.42, 0.46, 0.12];
        let (k_map, score) = find_map_estimate(&dist, 2);
        assert_eq!(k_map, 1);
        assert!((score - 50.0).abs() < TEST_EPSILON); // 1/2 × 100
    }

    #[test]
    fn test_find_map_estimate_tie_takes_higher_k() {
        // Equal probabilities, higher k wins
        let dist = vec![0.5, 0.5];
        let (k_map, _) = find_map_estimate(&dist, 10);
        assert_eq!(k_map, 1); // higher k wins on tie
    }

    #[test]
    fn test_find_map_estimate_all_stable() {
        // P(0)=1.0, k_map=0, score=0%
        let dist = vec![1.0, 0.0];
        let (k_map, score) = find_map_estimate(&dist, 2);
        assert_eq!(k_map, 0);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_find_map_estimate_all_unstable() {
        // P(2)=1.0, k_map=2, score=100%
        let dist = vec![0.0, 0.0, 1.0];
        let (k_map, score) = find_map_estimate(&dist, 2);
        assert_eq!(k_map, 2);
        assert!((score - 100.0).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_find_map_estimate_denominator() {
        // k_map=2, total=10, 2/10 × 100 = 20%
        let dist = vec![0.1, 0.2, 0.7];
        let (k_map, score) = find_map_estimate(&dist, 10);
        assert_eq!(k_map, 2);
        assert!((score - 20.0).abs() < TEST_EPSILON);
    }

    /* ======= calculate_msi_metrics tests =========== */

    #[test]
    fn test_calculate_msi_metrics_empty() {
        let filtered = FilteredRegions {
            variants: vec![],
            region_starts: vec![],
        };

        /* In theory should never come to this
            still keeping the test for safety and future
            flexibility.
        */
        let output_req = OutputRequirements {
            needs_pseudotime: false,
            needs_distribution: false,
            needs_heatmap: false,
        };

        let result =
            calculate_msi_metrics(&filtered, 100, 3.5, 0.0, "sample1".to_string(), output_req);

        assert!(result.k_map.is_none());
        assert!(result.msi_score_map.is_none());
        assert!(result.regions_with_variants.is_none());
        assert!(result.msi_status.is_none());
    }

    #[test]
    fn test_calculate_msi_metrics_with_uncertainty() {
        let v1 = make_variant(0.1, 0.8);
        let v2 = make_variant(0.2, 0.9);

        let filtered = FilteredRegions {
            variants: vec![&v1, &v2],
            region_starts: vec![0, 1],
        };

        let output_req = OutputRequirements {
            needs_pseudotime: true,
            needs_distribution: false,
            needs_heatmap: false,
        };

        let result =
            calculate_msi_metrics(&filtered, 100, 3.5, 0.5, "sample1".to_string(), output_req);
        let lower = result.uncertainty_lower.unwrap();
        let upper = result.uncertainty_upper.unwrap();
        let std_dev: f64 = result.map_std_dev.unwrap();

        assert!(result.distribution.is_none());
        assert!(lower <= result.msi_score_map.unwrap());
        assert!(result.msi_score_map.unwrap() <= upper);
        assert!(std_dev >= 0.0);
        assert_eq!(result.k_map.unwrap(), 2);
        assert!(result.uncertainty_lower.is_some());
        assert!(result.uncertainty_upper.is_some());
        assert!(result.map_std_dev.is_some());
        assert!(result.distribution.is_none());
    }

    #[test]
    fn test_calculate_msi_metrics_distribution_only_at_af_zero() {
        let v1 = make_variant(0.1, 0.8);

        let filtered = FilteredRegions {
            variants: vec![&v1],
            region_starts: vec![0],
        };

        let output_req = OutputRequirements {
            needs_pseudotime: false,
            needs_distribution: true,
            needs_heatmap: false,
        };

        // AF=0.0 should include distribution
        let result_af0 =
            calculate_msi_metrics(&filtered, 100, 3.5, 0.0, "sample1".to_string(), output_req);
        assert!(result_af0.distribution.is_some());
        let dist = result_af0.distribution.unwrap();
        assert_eq!(dist.len(), 2);
        let prob_sum: f64 = dist.iter().map(|d| d.probability).sum();
        assert!((prob_sum - 1.0).abs() < TEST_EPSILON);
        assert_eq!(dist[0].k, 0);
        assert_eq!(dist[1].k, 1);

        // AF=0.5 should NOT include distribution
        let result_af5 =
            calculate_msi_metrics(&filtered, 100, 3.5, 0.5, "sample1".to_string(), output_req);
        assert!(result_af5.distribution.is_none());
    }

    /* ==== AF Analysis Utilities ==================== */

    #[test]
    fn test_get_window_slice_basic() {
        let regions = vec![
            make_region("chr1:100-130", "chr1", 100, vec![], true),
            make_region("chr1:500-530", "chr1", 500, vec![], true),
            make_region("chr1:1100-1130", "chr1", 1100, vec![], true),
            make_region("chr2:100-130", "chr2", 100, vec![], true),
        ];
        let slice = get_window_slice(&regions, "chr1", 0, 1000);
        assert_eq!(slice.len(), 2);
        assert_eq!(slice[0].start, 100);
        assert_eq!(slice[1].start, 500);

        let slice = get_window_slice(&regions, "chr1", 1000, 2000);
        assert_eq!(slice.len(), 1);
        assert_eq!(slice[0].start, 1100);

        let slice = get_window_slice(&regions, "chr2", 0, 1000);
        assert_eq!(slice.len(), 1);
        assert_eq!(slice[0].start, 100);
    }

    #[test]
    fn test_get_window_slice_empty_window() {
        let regions = vec![make_region("chr1:100-130", "chr1", 100, vec![], true)];
        let slice = get_window_slice(&regions, "chr1", 5000, 6000);
        assert_eq!(slice.len(), 0);
    }

    #[test]
    fn test_get_window_slice_chrom_not_found() {
        let regions = vec![make_region("chr1:100-130", "chr1", 100, vec![], true)];
        let slice = get_window_slice(&regions, "chr99", 0, 1000);
        assert_eq!(slice.len(), 0);
    }

    #[test]
    fn test_get_window_slice_all_regions_in_window() {
        let regions = vec![
            make_region("chr1:100-130", "chr1", 100, vec![], true),
            make_region("chr1:200-230", "chr1", 200, vec![], true),
        ];
        let slice = get_window_slice(&regions, "chr1", 0, 1_000_000);
        assert_eq!(slice.len(), 2);
    }

    #[test]
    fn test_get_window_slice_boundary_inclusive_exclusive() {
        let regions = vec![
            make_region("chr1:0-30", "chr1", 0, vec![], true),
            make_region("chr1:1000-1030", "chr1", 1000, vec![], true),
            make_region("chr1:2000-2030", "chr1", 2000, vec![], true),
        ];

        // window_start is inclusive: region at exactly window_start=1000 included
        let slice = get_window_slice(&regions, "chr1", 1000, 1500);
        assert_eq!(slice.len(), 1);
        assert_eq!(slice[0].start, 1000);

        // window_end is exclusive: region at exactly window_end=2000 excluded
        let slice = get_window_slice(&regions, "chr1", 0, 1000);
        assert_eq!(slice.len(), 1);
        assert_eq!(slice[0].start, 0);
    }

    /* ==== run_af_evolution_analysis tests ========== */

    #[test]
    fn test_run_af_evolution_analysis_basic() {
        let regions = vec![make_region_simple(vec![make_variant(0.1, 0.5)])];

        let output_req = OutputRequirements {
            needs_pseudotime: false,
            needs_distribution: true,
            needs_heatmap: false,
        };

        let results = run_af_evolution_analysis(
            &regions,
            100,
            "sample1",
            3.5,
            &vec![0.0],
            output_req,
            Some(1),
        )
        .unwrap();

        let result = &results["0.00"];

        assert!(result.k_map.is_none());
        assert!(result.msi_score_map.is_none());
        assert!(result.regions_with_variants.is_none());
        assert!(result.msi_status.is_none());
        assert!(result.uncertainty_lower.is_none());
        assert!(result.uncertainty_upper.is_none());
        assert!(result.map_std_dev.is_none());

        // Distribution SHOULD exist (needs_distribution: true AND af=0.0)
        assert!(result.distribution.is_some());
        let dist = result.distribution.as_ref().unwrap();
        assert_eq!(dist.len(), 2);
        let prob_sum: f64 = dist.iter().map(|d| d.probability).sum();
        assert!((prob_sum - 1.0).abs() < TEST_EPSILON);
        // With prob_absent=0.1, expect P(k=1) > P(k=0)
        assert!(dist[1].probability > dist[0].probability);
    }

    #[test]
    fn test_run_af_evolution_analysis_multiple_af_thresholds() {
        let regions = vec![
            make_region(
                "chr1:100-130",
                "chr1",
                100,
                vec![make_variant(0.1, 0.9)],
                true,
            ),
            make_region(
                "chr1:200-230",
                "chr1",
                200,
                vec![make_variant(0.05, 0.5)],
                true,
            ),
        ];

        let output_req = OutputRequirements {
            needs_pseudotime: true,
            needs_distribution: false,
            needs_heatmap: false,
        };

        let results = run_af_evolution_analysis(
            &regions,
            100,
            "sample1",
            3.5,
            &vec![0.0, 0.5, 1.0],
            output_req,
            Some(1),
        )
        .unwrap();

        assert_eq!(results.len(), 3, "Should have 3 AF threshold results");

        // AF=1.0: no variants pass (0.9 < 1.0, 0.5 < 1.0)
        let af_1_0 = &results["1.00"];
        assert_eq!(
            af_1_0.regions_with_variants.unwrap(),
            0,
            "No regions should pass AF=1.0 threshold"
        );
        assert_eq!(
            af_1_0.k_map.unwrap(),
            0,
            "k_map should be 0 with no regions"
        );
        assert_eq!(af_1_0.msi_score_map.unwrap(), 0.0, "MSI score should be 0%");
        assert_eq!(
            af_1_0.msi_status.as_deref(),
            Some("MSS"),
            "Status should be MSS"
        );
        assert!(
            af_1_0.uncertainty_lower.is_some(),
            "Uncertainty should exist with pseudotime=true"
        );

        // AF=0.5: both variants pass (0.9 ≥ 0.5, 0.5 ≥ 0.5)
        let af_0_5 = &results["0.50"];
        assert_eq!(
            af_0_5.regions_with_variants.unwrap(),
            2,
            "Both regions should pass AF=0.5 threshold"
        );
        assert_eq!(
            af_0_5.k_map.unwrap(),
            2,
            "k_map should be 2 with high instability"
        );
        assert!(
            (af_0_5.msi_score_map.unwrap() - 2.0).abs() < TEST_EPSILON,
            "MSI score should be 2.0%, got {}",
            af_0_5.msi_score_map.unwrap()
        );
        assert_eq!(af_0_5.msi_status.as_deref(), Some("MSS"), "2% < 3.5% : MSS");
        assert!(af_0_5.uncertainty_lower.is_some());
        assert!(af_0_5.uncertainty_upper.is_some());

        // AF=0.0: both variants pass (all pass)
        let af_0_0 = &results["0.00"];
        assert_eq!(
            af_0_0.regions_with_variants.unwrap(),
            2,
            "Both regions should pass AF=0.0 threshold"
        );
        assert_eq!(af_0_0.k_map.unwrap(), 2);
        assert!((af_0_0.msi_score_map.unwrap() - 2.0).abs() < TEST_EPSILON);
        assert!(af_0_0.uncertainty_lower.is_some());

        // Verify no distribution (needs_distribution=false)
        assert!(af_1_0.distribution.is_none());
        assert!(af_0_5.distribution.is_none());
        assert!(af_0_0.distribution.is_none());

        assert_eq!(
            af_0_5.msi_score_map, af_0_0.msi_score_map,
            "MSI scores should match when same regions pass"
        );
    }
}
