//! bcf_utils.rs
//!
//! Utilities for VCF/BCF file handling and parsing.
//!
//! This module provides:
//! 1. Sample information extraction
//! 2. Record field extraction (chromosome, SVLEN, probabilities, allele frequencies)
//! 3. Allele type classification (indel, symbolic, breakend, reference, spanning deletion)
//! 4. VCF/BCF fields validation
//! 5. VCF file validation

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use log::{debug, info, warn};
use rust_htslib::bcf::header::{HeaderView, TagLength, TagType};
use rust_htslib::bcf::{self, record::Numeric, Read};

use crate::errors::Error;
use crate::utils::genomics::calculate_dynamic_svlen;
use crate::utils::is_phred_scaled;
use crate::utils::stats::phred_to_prob;

const EPSILON: f64 = 1e-6;

/* ============ Data Structures =================== */

/// Sample information extracted from VCF header.
///
/// Contains sample names and their corresponding indices
/// in the VCF header for efficient lookup during processing.
#[derive(Debug, Clone)]
pub(crate) struct SampleInfo {
    /// Sample names to process
    pub samples: Vec<String>,
    /// Map of sample name to VCF header index
    pub samples_index_map: HashMap<String, usize>,
}

/* ================================================ */

/* ========= BCF Extraction Functions ============= */

/// Extract sample names from VCF header
///
/// # Arguments
/// * `vcf` - VCF reader
///
/// # Returns
/// Vector of sample names as Strings
///
/// # Example
/// assert_eq!(extract_sample_names(&vcf), vec!["sample1", "sample2"]);
pub(crate) fn extract_sample_names(vcf: &bcf::Reader) -> Vec<String> {
    let header = vcf.header();
    header
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect()
}

/// Get chromosome name from a VCF record
///
/// # Arguments
/// * `record` - VCF record
/// * `header` - VCF header (for resolving RID to name)
///
/// # Returns
/// Chromosome name as String
///
/// # Errors
/// Returns error if:
/// - RID is missing in record
/// - Chromosome name resolution fails
///
/// # Example
/// assert_eq!(get_chrom(&record, &header).unwrap(), "chr1");
pub(crate) fn get_chrom(record: &bcf::Record, header: &HeaderView) -> Result<String> {
    let rid = record.rid().ok_or_else(|| Error::VcfRecordChromMissing {
        pos: record.pos() + 1,
    })?;

    let chrom_bytes = header
        .rid2name(rid)
        .map_err(|_| Error::VcfRecordChromResolveFailed {
            pos: record.pos() + 1,
            rid,
            details: "Failed to resolve chromosome name".to_string(),
        })?;

    let chrom = String::from_utf8_lossy(chrom_bytes).to_string();

    Ok(chrom)
}

/// Get SVLEN from INFO field or calculate dynamically.
///
/// Returns length difference between ALT and REF alleles.
/// (above referes to SVLEN, not simple length difference)
/// Positive = insertion, Negative = deletion.
///
/// # Arguments
/// * `record` - VCF record
/// * `alt_idx` - Index into ALT alleles (0 = first ALT)
/// * `ref_seq` - Reference allele sequence
/// * `alt_seq` - Alternate allele sequence
///
/// # Example
/// REF=ACAG, ALT=ACAGCAG : SVLEN=+3
/// assert_eq!(get_svlen(&record, 0, b"ACAG", b"ACAGCAG").unwrap(), 3);
pub(crate) fn get_svlen(
    record: &bcf::Record,
    alt_idx: usize,
    ref_seq: &[u8],
    alt_seq: &[u8],
) -> Result<i32> {
    // Try to get SVLEN from INFO field first
    if let Ok(Some(svlens)) = record.info(b"SVLEN").integer() {
        if let Some(&svlen) = svlens.get(alt_idx) {
            if !svlen.is_missing() {
                return Ok(svlen);
            }
        }
    }

    // Fallback: calculate dynamically using anchor detection
    Ok(calculate_dynamic_svlen(ref_seq, alt_seq))
}

/// Get combined probability that variant is absent (artifact)
///
/// Extracts and combines two sources of artifact evidence:
/// 1. PROBABILITY_ABSENT tag: Primary absent probability
/// 2. PROBABILITY_ARTIFACT tag: Additional artifact probability
///
/// NOTE: Both can be either PHRED-scaled or direct probabilities.
///
/// # Formula
/// ```text
/// P(absent) = P(artifact_PA) + P(artifact_ART)
/// ```
/// Notes: Currently caters Phred and Linear Prob.
/// Check if LogProb is also possible and extend accordingly if required.
///
/// # Arguments
/// * `record` - VCF record
/// * `header` - VCF header (for error messages)
/// * `alt_idx` - ALT allele index (0 = first ALT)
/// * `is_phred` - Whether probabilities are PHRED-scaled
///
/// # Returns
/// - `Some(prob)` if both PROB_ABSENT and PROB_ARTIFACT exist for this ALT
/// - `None` if either tag is missing or missing at this ALT index
///
/// # Errors
/// Returns error if:
/// - Probability is NaN
/// - Final combined probability is outside [0, 1]
pub(crate) fn get_prob_absent(
    record: &bcf::Record,
    header: &HeaderView,
    alt_idx: usize,
    is_phred: bool,
) -> Result<Option<f64>> {
    let prob_absent = match record.info(b"PROB_ABSENT").float()? {
        Some(p) if alt_idx < p.len() => p[alt_idx],
        _ => return Ok(None),
    };

    let prob_artifact = match record.info(b"PROB_ARTIFACT").float()? {
        Some(p) if alt_idx < p.len() => p[alt_idx],
        _ => return Ok(None),
    };

    if prob_absent.is_missing() || prob_artifact.is_missing() {
        return Ok(None);
    }

    if prob_absent.is_nan() {
        return Err(Error::VcfProbabilityValueInvalid {
            field: "PROB_ABSENT".to_string(),
            value: prob_absent,
            chrom: get_chrom(record, header)?,
            pos: record.pos() + 1,
        }
        .into());
    }

    if prob_artifact.is_nan() {
        return Err(Error::VcfProbabilityValueInvalid {
            field: "PROB_ARTIFACT".to_string(),
            value: prob_artifact,
            chrom: get_chrom(record, header)?,
            pos: record.pos() + 1,
        }
        .into());
    }

    let probability_absent = if is_phred {
        phred_to_prob(prob_absent as f64)
    } else {
        prob_absent as f64
    };

    let probability_artifact = if is_phred {
        phred_to_prob(prob_artifact as f64)
    } else {
        prob_artifact as f64
    };

    let probability = probability_absent + probability_artifact;

    // Optional: Adjust with gnomAD population frequency
    let gnomad_af = match record.info(b"POPULATION_AF").float() {
        Ok(Some(p)) if alt_idx < p.len() && !p[alt_idx].is_missing() => p[alt_idx] as f64,
        Ok(Some(_)) | Ok(None) | Err(_) => 0.0, // Missing/invalid: use 0.0
    };

    // Efficient calculation:
    // prob_absent_final = 1 - ((1 - prob_absent_base) × (1 - gnomad_af))
    // Expanded: prob_absent_base + gnomad_af × (1 - prob_absent_base)
    let probability_final = probability + gnomad_af * (1.0 - probability);

    let valid_range = -EPSILON..=1.0 + EPSILON;
    if !valid_range.contains(&probability_final) {
        return Err(Error::VcfProbabilityValueInvalid {
            field: "ADJUSTED PROBABILITY ABSENT (with gnomAD)".to_string(),
            value: probability_final as f32,
            chrom: get_chrom(record, header)?,
            pos: record.pos() + 1,
        }
        .into());
    }

    let probability = probability_final.clamp(0.0, 1.0);

    Ok(Some(probability))
}

/// Extract per-sample allele frequencies for a specific ALT allele.
///
/// Reads FORMAT:AF field and returns AF values for each sample
/// in the samples_index_map.
///
/// # Arguments
/// * `record` - VCF record
/// * `header` - VCF header
/// * `samples_index_map` - Map of sample names to VCF indices
/// * `alt_idx` - Index into ALT alleles (0 = first ALT)
///
/// # Returns
/// HashMap of sample name to AF value. Empty if AF field missing.
///
/// # Errors
/// Returns error if AF value is NaN or outside [0.0, 1.0].
pub(crate) fn get_sample_afs(
    record: &bcf::Record,
    header: &HeaderView,
    samples_index_map: &HashMap<String, usize>,
    alt_idx: usize,
) -> Result<HashMap<String, f64>> {
    let mut sample_afs = HashMap::new();

    let afs = match record.format(b"AF").float() {
        Ok(a) => a,
        Err(_) => {
            warn!(
                "AF field missing at {}:{} - variant will have no AF data",
                get_chrom(record, header)?,
                record.pos() + 1
            );
            return Ok(sample_afs);
        }
    };

    for (sample_name, &vcf_header_idx) in samples_index_map {
        let Some(sample_af_values) = afs.get(vcf_header_idx) else {
            continue;
        };

        let Some(&af) = sample_af_values.get(alt_idx) else {
            continue;
        };

        if af.is_missing() {
            continue;
        }

        if af.is_nan() || !(0.0..=1.0).contains(&af) {
            return Err(Error::VcfAlleleFrequencyInvalid {
                sample: sample_name.clone(),
                af,
                chrom: get_chrom(record, header)?,
                pos: record.pos() + 1,
            }
            .into());
        }

        sample_afs.insert(sample_name.clone(), af as f64);
    }

    if sample_afs.is_empty() {
        debug!(
            "No valid AF values for any sample at {}:{} - variant will be skipped in analysis",
            get_chrom(record, header)?,
            record.pos() + 1
        );
    }

    Ok(sample_afs)
}

/* ================================================ */

/* ===== BCF Specification Check Functions ======== */

/// Check if VCF uses PHRED-scaled probabilities from file path.
///
/// Opens the VCF file and checks the description of PROB_{event}s
/// to determine if probabilities are PHRED-scaled or linear.
///
/// # Arguments
/// * `vcf_path` - Path to VCF/BCF file
///
/// # Returns
/// * `Ok(true)` if probabilities are PHRED-scaled
/// * `Ok(false)` if probabilities are linear
/// * `Err` if file cannot be opened
///
/// # Example
/// assert!(is_phred_scaled_from_path(Path::new("variants.vcf")).is_ok());
pub(crate) fn is_phred_scaled_from_path(vcf_path: &Path) -> Result<bool> {
    let vcf = bcf::Reader::from_path(vcf_path)?;
    Ok(is_phred_scaled(&vcf))
}

/* ================================================ */

/* ======== BCF Allele Type Check Functions ======= */

// /// Check if all alleles are SNVs (all same length as ref)
// /// Example: REF=A, ALT=T,G -> all length 1
// NOTE: Not required, can be toggled on if required in other utilities.
// pub fn all_alleles_snv(ref_allele: &[u8], alt_alleles: &[&[u8]]) -> bool {
//     let ref_allele_len = ref_allele.len();
//     alt_alleles.iter().all(|alt| alt.len() == ref_allele_len)
// }

/// Check if allele represents no variant (reference)
///
/// Matches alleles that indicate "no alternative allele":
/// - `.` - Missing/no ALT allele (VCF spec)
/// - `<REF>` - Explicit reference allele (rare)
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is reference, `false` otherwise
///
/// # Example
/// assert!(is_reference_allele(b"."));
pub(crate) fn is_reference_allele(allele: &[u8]) -> bool {
    allele == b"." || allele == b"<REF>"
}

/// Check if allele is symbolic (starts with <)
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is symbolic, `false` otherwise
///
/// Examples
///  assert!(is_symbolic(b"<DEL>"));
pub(crate) fn is_symbolic(allele: &[u8]) -> bool {
    allele.len() >= 3 && allele.starts_with(b"<") && allele.ends_with(b">")
}

/// Check if allele is a breakend (contains [ or ])
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is a breakend, `false` otherwise
///
/// # Examples
/// assert!(is_breakend(b"A[chr2:100["));
pub(crate) fn is_breakend(allele: &[u8]) -> bool {
    allele.iter().any(|&c| c == b'[' || c == b']')
}

/// Check if allele is a spanning deletion (*)
///
/// # Arguments
/// * `allele` - Allele sequence as byte slice
///
/// # Returns
/// `true` if allele is a spanning deletion, `false` otherwise
///
/// # Example
/// assert!(is_spanning_deletion(b"*"));
pub(crate) fn is_spanning_deletion(allele: &[u8]) -> bool {
    allele == b"*"
}

/* ================================================ */

/* ========= BCF Validation Functions ============= */

/// Validate a single VCF header field has correct type and number.
///
/// # Arguments
/// * `field_result` - Result from `header.info_type()` or `header.format_type()`
/// * `location` - "INFO" or "FORMAT" (for error messages)
/// * `field_name` - Name of the field being validated
/// * `expected_type` - Expected TagType (e.g., `TagType::Float`)
/// * `expected_length` - Expected TagLength (e.g., `TagLength::AltAlleles`)
///
/// # Returns
/// * `Ok(())` if field exists with correct type and number
/// * `Err` if field missing or has wrong type/number
fn validate_vcf_header_field(
    field_result: std::result::Result<(TagType, TagLength), rust_htslib::errors::Error>,
    location: &str,
    field_name: &str,
    expected_type: TagType,
    expected_length: TagLength,
) -> Result<()> {
    match field_result {
        Ok((actual_type, actual_length)) => {
            if actual_type != expected_type || actual_length != expected_length {
                Err(Error::VcfHeaderFieldTypeInvalid {
                    location: location.to_string(),
                    field: field_name.to_string(),
                    expected: format!("Type={:?}, Number={:?}", expected_type, expected_length),
                    found: format!("Type={:?}, Number={:?}", actual_type, actual_length),
                }
                .into())
            } else {
                Ok(())
            }
        }
        Err(_) => Err(Error::VcfHeaderFieldMissing {
            field: field_name.to_string(),
            location: location.to_string(),
        }
        .into()),
    }
}

/// Validate that required samples exist in VCF header.
///
/// Checks VCF header for presence of all requested sample names.
///
/// # Arguments
/// * `header` - VCF header to validate
/// * `required_samples` - Slice of sample names to validate
///
/// # Returns
/// * `Ok(())` if all samples exist
/// * `Err` if any sample is missing
///
/// # Errors
/// Returns error if any sample name is not found in VCF header.
/// Error message includes comma-separated list of missing samples.
///
/// # Example
/// assert!(validate_samples_exist(&header, &vec!["tumor".to_string()]).is_ok());
pub(crate) fn validate_samples_exist(
    header: &HeaderView,
    required_samples: &[String],
) -> Result<()> {
    let mut missing = Vec::new();

    for sample_name in required_samples {
        if header.sample_id(sample_name.as_bytes()).is_none() {
            missing.push(sample_name.clone());
        }
    }

    if !missing.is_empty() {
        return Err(Error::VcfSamplesNotFound {
            sample: missing.join(", "),
        }
        .into());
    }

    info!("  - Samples validated: {:?}", required_samples);

    Ok(())
}

/// Validate that required events exist in VCF header.
///
/// Checks VCF header for presence of INFO/PROB_{EVENT} fields
/// for each requested event name. Event names are automatically
/// converted to uppercase and prefixed with "PROB_".
///
/// # Arguments
/// * `header` - VCF header to check for event fields
/// * `event_names` - Slice of event names to validate (e.g., ["somatic_tumor"])
///
/// # Returns
/// * `Ok(())` if all event fields exist
/// * `Err` if any event field is missing
///
/// # Errors
/// Returns error if any INFO/PROB_{EVENT} field is not found in VCF header.
/// Error message includes comma-separated list of missing events.
///
/// # Event Field Mapping
/// Event name          -> INFO field checked
/// "somatic_tumor"     -> INFO/PROB_SOMATIC_TUMOR
/// "germline_normal"   -> INFO/PROB_GERMLINE_NORMAL
///
/// # Example
/// assert!(validate_events_exist(&header, &vec!["somatic_tumor".to_string()]).is_ok());
pub(crate) fn validate_events_exist(header: &HeaderView, event_names: &[String]) -> Result<()> {
    let mut missing_events = Vec::new();

    for event_name in event_names {
        let field_name = format!("PROB_{}", event_name.to_uppercase());

        if header.info_type(field_name.as_bytes()).is_err() {
            missing_events.push(event_name.clone());
        }
    }

    if !missing_events.is_empty() {
        return Err(Error::VcfEventsMissing {
            events: missing_events.join(", "),
        }
        .into());
    }

    info!("  - Events validated: {:?}", event_names);
    Ok(())
}

/// Validate VCF header contains required fields for MSI analysis.
///
/// Validates that mandatory INFO or FORMAT fields exist with correct types:
/// - FORMAT:AF (Type=Float, Number=A)
///
/// # Arguments
/// * `header` - VCF header to validate
///
/// # Returns
/// * `Ok(())` if all required fields present with correct types
/// * `Err` if any field is missing or has incorrect type
pub(crate) fn validate_required_vcf_fields_msi(header: &HeaderView) -> Result<()> {
    validate_vcf_header_field(
        header.format_type(b"AF"),
        "FORMAT",
        "AF",
        TagType::Float,
        TagLength::AltAlleles,
    )?;

    Ok(())
}

/// Validate VCF file has at least one variant.
///
/// Simple sanity check - ensures file is not empty and logs first variant position.
///
/// # Arguments
/// * `vcf` - VCF reader (must be at start of file)
///
/// # Returns
/// * `Ok(())` if file has at least one variant
/// * `Err` if file is empty or read fails
pub(crate) fn validate_vcf_file(vcf: &mut bcf::Reader) -> Result<()> {
    let header = vcf.header().clone();

    match vcf.records().next() {
        None => Err(Error::VcfFileEmpty.into()),
        Some(Err(e)) => Err(Error::VcfRecordReadFailed {
            details: e.to_string(),
        }
        .into()),
        Some(Ok(record)) => {
            let chrom = get_chrom(&record, &header)?;
            let pos = record.pos();
            info!("  - First variant: {}:{}", chrom, pos + 1);
            info!("  - VCF file validated successfully");
            Ok(())
        }
    }
}

/// Validate VCF file and extract sample information.
///
/// Performs validation checks:
/// 1. File can be opened
/// 2. Contains required fields for MSI analysis
/// 3. Contains at least one sample
/// 4. Excluded samples exist in file
/// 5. At least one sample remains after exclusion
/// 6. Contains at least one variant record
///
/// # Arguments
/// * `vcf_path` - Path to VCF/BCF file
/// * `samples_exclusion` - Sample names to exclude from processing
///
/// # Returns
/// * `Ok((SampleInfo, bool))` - Sample info and whether probabilities are PHRED-scaled
/// * `Err` - Validation failed
///
/// # Example
/// assert!(validate_vcf_file(&vcf_path, &vec!["sample1".to_string()]).is_ok());
pub(crate) fn validate_vcf_file_ex(
    vcf_path: &Path,
    samples_exclusion: &[String],
) -> Result<(SampleInfo, bool)> {
    info!("Validating VCF file format: {}", vcf_path.display());

    let mut vcf = bcf::Reader::from_path(vcf_path).context("Failed to open VCF/BCF file")?;
    let header = vcf.header().clone();

    // Validate required fields for MSI analysis
    validate_required_vcf_fields_msi(&header)?;
    info!("  - Required header fields validated: PROB_ABSENT, PROB_ARTIFACT, AF");

    // Check if probabilities are PHRED-scaled
    let is_phred = is_phred_scaled(&vcf);
    info!(
        "  - Probabilities are {} scaled",
        if is_phred { "PHRED" } else { "linear" }
    );

    // Extract sample names
    let sample_names = extract_sample_names(&vcf);

    if sample_names.is_empty() {
        return Err(Error::VcfSamplesMissing.into());
    }

    let mut invalid_exclusions: Vec<String> = Vec::new();

    for excluded_sample in samples_exclusion {
        if !sample_names.contains(excluded_sample) {
            invalid_exclusions.push(excluded_sample.clone());
        }
    }

    if !invalid_exclusions.is_empty() {
        return Err(Error::VcfSampleExclusionInvalid {
            samples: invalid_exclusions.join(", "),
        }
        .into());
    }

    let mut remaining_samples: Vec<String> = vec![];
    let mut samples_index_map: HashMap<String, usize> = HashMap::new();

    for (i, s) in sample_names.iter().enumerate() {
        if !samples_exclusion.contains(s) {
            remaining_samples.push(s.to_string());
            samples_index_map.insert(s.to_string(), i);
        }
    }

    if remaining_samples.is_empty() {
        return Err(Error::VcfSamplesEmptyAfterExclusion.into());
    }

    info!("  - Samples to process: {}", remaining_samples.len());
    info!("  - Sample names: {:?}", remaining_samples);

    // Check for at least one variant record
    match vcf.records().next() {
        None => {
            return Err(Error::VcfFileEmpty.into());
        }
        Some(Err(e)) => {
            return Err(Error::VcfRecordReadFailed {
                details: e.to_string(),
            }
            .into());
        }
        Some(Ok(record)) => {
            let chrom = get_chrom(&record, &header)?;
            let pos = record.pos();
            info!("  - First variant: {}:{}", chrom, pos + 1);
        }
    }

    info!("  - VCF file format validated successfully");

    Ok((
        SampleInfo {
            samples: remaining_samples,
            samples_index_map,
        },
        is_phred,
    ))
}

/* ================================================ */

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use rust_htslib::bcf::record::Numeric;
    use tempfile::NamedTempFile;

    use crate::utils::genomics::is_indel;
    use crate::utils::stats::TEST_EPSILON;

    /// Encodes a single allele for the GT field in BCF format.
    ///
    /// BCF stores each allele as: (allele_index + 1) << 1 | phased_bit
    /// - allele_index: -1 for missing, 0 for REF, 1 for ALT1, etc.
    /// - phased_bit: 1 for phased ('|'), 0 for unphased ('/')
    ///
    /// Example:
    ///   encode_genotype_allele(0, false) -> 2   // allele 0, unphased
    ///
    /// A full genotype like "0/0" would be encoded as: [2, 2]
    ///
    /// Reference: https://samtools.github.io/hts-specs/BCFv2_qref.pdf
    pub(crate) fn encode_genotype_allele(allele_index: i32, phased: bool) -> i32 {
        let phased_flag = if phased { 1 } else { 0 };
        (allele_index + 1) * 2 | phased_flag
    }

    /// Configuration for test VCF creation
    pub(crate) struct TestVcfConfig<'a> {
        pub ref_allele: &'a [u8],
        pub alt_alleles: Vec<&'a [u8]>,
        pub af_values: Option<Vec<f32>>,
        pub prob_absent: Option<Vec<f32>>,
        pub prob_artifact: Option<Vec<f32>>,
        pub prob_somatic: Option<Vec<f32>>,
        pub prob_high_vaf: Option<Vec<f32>>,
        pub num_samples: usize,
        pub use_phred: bool,
    }

    impl<'a> Default for TestVcfConfig<'a> {
        fn default() -> Self {
            Self {
                ref_allele: b"A",
                alt_alleles: vec![b"AT"],
                af_values: None,
                prob_absent: None,
                prob_artifact: None,
                prob_somatic: None,
                prob_high_vaf: None,
                num_samples: 2,
                use_phred: false,
            }
        }
    }

    pub(crate) fn create_test_vcf(config: TestVcfConfig) -> (NamedTempFile, Vec<String>) {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        let mut header = rust_htslib::bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##);

        // Conditionally set PROB_ABSENT, PROB_ARTIFACT and PROB_SOMATIC headers based on use_phred
        if config.use_phred {
            header.push_record(br##"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,Description="Probability absent (PHRED)">"##);
            header.push_record(br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Probability artifact (PHRED)">"##);
            header.push_record(br##"##INFO=<ID=PROB_SOMATIC,Number=A,Type=Float,Description="Probability somatic (PHRED)">"##);
            header.push_record(br##"##INFO=<ID=PROB_HIGH_VAF,Number=A,Type=Float,Description="Probability high VAF (PHRED)">"##);
        } else {
            header.push_record(br##"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,Description="Probability absent (linear)">"##);
            header.push_record(br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Probability artifact (linear)">"##);
            header.push_record(br##"##INFO=<ID=PROB_SOMATIC,Number=A,Type=Float,Description="Probability somatic (linear)">"##);
            header.push_record(br##"##INFO=<ID=PROB_HIGH_VAF,Number=A,Type=Float,Description="Probability high VAF (linear)">"##);
        }

        header.push_record(br##"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"##);
        header.push_record(
            br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">"##,
        );

        let mut sample_names = Vec::new();
        for i in 0..config.num_samples {
            let name = format!("sample{}", i + 1);
            header.push_sample(name.as_bytes());
            sample_names.push(name);
        }

        let mut wtr = rust_htslib::bcf::Writer::from_path(
            path,
            &header,
            false,
            rust_htslib::bcf::Format::Vcf,
        )
        .unwrap();

        let mut rec = wtr.empty_record();
        rec.set_rid(Some(0));
        rec.set_pos(99);

        // Build alleles: [REF, ALT1, ALT2, ...]
        let mut alleles_vec = vec![config.ref_allele];
        alleles_vec.extend(config.alt_alleles.iter());
        rec.set_alleles(&alleles_vec).unwrap();

        // SVLEN for each ALT
        let svlen_values: Vec<i32> = config
            .alt_alleles
            .iter()
            .map(|alt| (alt.len() as i32) - (config.ref_allele.len() as i32))
            .collect();
        rec.push_info_integer(b"SVLEN", &svlen_values).unwrap();

        // PROB_ABSENT
        let prob_absent = match config.prob_absent {
            Some(pa) => pa,
            None => {
                if config.use_phred {
                    vec![20.0; config.alt_alleles.len()]
                } else {
                    vec![0.005; config.alt_alleles.len()]
                }
            }
        };
        rec.push_info_float(b"PROB_ABSENT", &prob_absent).unwrap();

        // PROB_ARTIFACT
        let prob_artifact = match config.prob_artifact {
            Some(pa) => pa,
            None => {
                if config.use_phred {
                    vec![10.0; config.alt_alleles.len()]
                } else {
                    vec![0.005; config.alt_alleles.len()]
                }
            }
        };
        rec.push_info_float(b"PROB_ARTIFACT", &prob_artifact)
            .unwrap();

        // PROB_SOMATIC
        let prob_somatic = match config.prob_somatic {
            Some(ps) => ps,
            None => {
                if config.use_phred {
                    vec![30.0; config.alt_alleles.len()]
                } else {
                    vec![0.001; config.alt_alleles.len()]
                }
            }
        };
        rec.push_info_float(b"PROB_SOMATIC", &prob_somatic).unwrap();

        // PROB_HIGH_VAF
        let prob_high_vaf = match config.prob_high_vaf {
            Some(phv) => phv,
            None => {
                if config.use_phred {
                    vec![40.0; config.alt_alleles.len()]
                } else {
                    vec![0.0001; config.alt_alleles.len()]
                }
            }
        };
        rec.push_info_float(b"PROB_HIGH_VAF", &prob_high_vaf)
            .unwrap();

        // GT
        let mut genotypes = Vec::new();
        for _ in 0..config.num_samples {
            genotypes.push(encode_genotype_allele(0, false));
            genotypes.push(encode_genotype_allele(1, false));
        }
        rec.push_format_integer(b"GT", &genotypes).unwrap();

        // AF values
        if let Some(afs) = config.af_values {
            rec.push_format_float(b"AF", &afs).unwrap();
        }

        wtr.write(&rec).unwrap();

        (tmp, sample_names)
    }

    pub(crate) fn create_multi_chromosome_vcf() -> NamedTempFile {
        let tmp = NamedTempFile::new().unwrap();
        let path = tmp.path();

        let mut header = rust_htslib::bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(br"##contig=<ID=chr2,length=1000000>");
        header.push_record(br"##contig=<ID=chrX,length=1000000>");
        header.push_record(br##"##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="SV length">"##);
        header.push_record(br##"##INFO=<ID=PROB_ABSENT,Number=A,Type=Float,Description="Probability absent (linear)">"##);
        header.push_record(br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Probability artifact (linear)">"##);
        header.push_record(
            br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">"##,
        );
        header.push_sample(b"sample1");

        let mut wtr = rust_htslib::bcf::Writer::from_path(
            path,
            &header,
            false,
            rust_htslib::bcf::Format::Vcf,
        )
        .unwrap();

        for (rid, pos) in &[(0, 100), (0, 200), (1, 150), (2, 175)] {
            let mut rec = wtr.empty_record();
            rec.set_rid(Some(*rid));
            rec.set_pos(*pos);
            rec.set_alleles(&[b"A", b"AT"]).unwrap();
            rec.push_info_integer(b"SVLEN", &[1]).unwrap();
            rec.push_info_float(b"PROB_ABSENT", &[0.01]).unwrap();
            rec.push_info_float(b"PROB_ARTIFACT", &[0.005]).unwrap();
            rec.push_format_float(b"AF", &[0.5]).unwrap();
            wtr.write(&rec).unwrap();
        }

        tmp
    }

    fn create_header_view(header_lines: &[&[u8]]) -> HeaderView {
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");

        for line in header_lines {
            header.push_record(line);
        }

        let tmp = NamedTempFile::new().unwrap();
        let writer = bcf::Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();
        drop(writer);

        let reader = bcf::Reader::from_path(tmp.path()).unwrap();
        reader.header().clone()
    }

    /// Read the first record from a VCF/BCF file (without header).
    ///
    /// # Returns
    /// Tuple of (reader, record) for testing
    pub fn read_first_record_simple(path: &Path) -> (bcf::Reader, bcf::Record) {
        let mut reader = bcf::Reader::from_path(path).unwrap();
        let record = reader.records().next().unwrap().unwrap();
        (reader, record)
    }

    /// Read the first record from a VCF/BCF file with header.
    ///
    /// Internally calls `read_first_record_simple` and clones the header.
    ///
    /// # Returns
    /// Tuple of (reader, header, record) for testing
    pub fn read_first_record(path: &Path) -> (bcf::Reader, HeaderView, bcf::Record) {
        let (reader, record) = read_first_record_simple(path);
        let header = reader.header().clone();
        (reader, header, record)
    }

    /// Create a test VCF record with specified parameters.
    ///
    /// # Arguments
    /// * `writer` - VCF writer
    /// * `rid` - Reference ID (chromosome index in header)
    /// * `pos` - Position (0-based)
    /// * `ref_allele` - Reference allele (e.g., b"A")
    /// * `alt_allele` - Alternate allele (e.g., b"T")
    ///
    /// # Returns
    /// Configured VCF record ready for testing
    pub fn create_test_record(
        writer: &bcf::Writer,
        rid: u32,
        pos: i64,
        ref_allele: &[u8],
        alt_allele: &[u8],
    ) -> bcf::Record {
        let mut record = writer.empty_record();
        record.set_rid(Some(rid));
        record.set_pos(pos);
        record.set_alleles(&[ref_allele, alt_allele]).unwrap();
        record
    }

    /* ==== BCF Extraction Function(s) tests ========= */

    #[test]
    fn test_get_chrom_and_svlen() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig::default());

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let chrom = get_chrom(&record, &header).unwrap();
        assert_eq!(chrom, "chr1");

        let svlen = get_svlen(&record, 0, b"A", b"AT").unwrap();
        assert_eq!(svlen, 1); // A -> AT is an insertion, svlen = 1
    }

    #[test]
    fn test_get_prob_absent() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            prob_absent: Some(vec![0.01]),
            prob_artifact: Some(vec![0.005]),
            num_samples: 1,
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let prob = get_prob_absent(&record, &header, 0, false)
            .unwrap()
            .unwrap();
        assert!((prob - 0.015).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_get_prob_absent_phred() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            use_phred: true,
            num_samples: 1,
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let prob = get_prob_absent(&record, &header, 0, true).unwrap().unwrap();
        // PHRED 20 ≈ 0.01, PHRED 10 ≈ 0.1, sum ≈ 0.11
        assert!((prob - 0.11).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_get_prob_absent_invalid_sum() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            prob_absent: Some(vec![0.7]),
            prob_artifact: Some(vec![0.6]), // Sum > 1.0
            num_samples: 1,
            ..Default::default()
        });

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let result = get_prob_absent(&record, &header, 0, false);
        assert!(result.is_err());
    }

    #[test]
    fn test_get_sample_afs() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![0.45, 0.98]),
            ..Default::default()
        });

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);
        samples_index_map.insert("sample2".to_string(), 1);

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let sample_afs = get_sample_afs(&record, &header, &samples_index_map, 0).unwrap();

        assert_eq!(sample_afs.len(), 2);
        assert!((sample_afs["sample1"] - 0.45).abs() < TEST_EPSILON);
        assert!((sample_afs["sample2"] - 0.98).abs() < TEST_EPSILON);
    }

    #[test]
    fn test_get_sample_afs_invalid() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![1.5]), // Invalid: > 1.0
            num_samples: 1,
            ..Default::default()
        });

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let result = get_sample_afs(&record, &header, &samples_index_map, 0);
        assert!(result.is_err());
    }

    #[test]
    fn test_get_sample_afs_missing() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            af_values: Some(vec![f32::missing()]),
            num_samples: 1,
            ..Default::default()
        });

        let mut samples_index_map = HashMap::new();
        samples_index_map.insert("sample1".to_string(), 0);

        let mut reader = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = reader.header().clone();
        let record = reader.records().next().unwrap().unwrap();

        let sample_afs = get_sample_afs(&record, &header, &samples_index_map, 0).unwrap();
        assert!(sample_afs.is_empty());
    }

    /* ====== BCF Specification check tests ===== ==== */

    #[test]
    fn test_is_phred_scaled_from_path_linear() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            use_phred: false,
            num_samples: 1,
            ..Default::default()
        });

        let is_phred = is_phred_scaled_from_path(tmp_vcf.path()).unwrap();
        assert!(!is_phred, "Should detect linear probabilities");
    }

    #[test]
    fn test_is_phred_scaled_from_path_phred() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            use_phred: true,
            num_samples: 1,
            ..Default::default()
        });

        let is_phred = is_phred_scaled_from_path(tmp_vcf.path()).unwrap();
        assert!(is_phred, "Should detect PHRED-scaled probabilities");
    }

    #[test]
    fn test_is_phred_scaled_from_path_invalid_file() {
        let result = is_phred_scaled_from_path(Path::new("/nonexistent/file.vcf"));
        assert!(result.is_err(), "Should fail for nonexistent file");
    }

    /* ====== BCF ALLELE Type tests ================== */

    #[test]
    fn test_allele_type_checks() {
        // Reference allele
        assert!(is_reference_allele(b"."));
        assert!(is_reference_allele(b"<REF>"));

        // Indel
        assert!(is_indel(b"AC", b"ACG")); // insertion
        assert!(is_indel(b"ACG", b"AC")); // deletion
        assert!(!is_indel(b"AC", b"AG"));

        // Symbolic
        assert!(is_symbolic(b"<DEL>"));
        assert!(is_symbolic(b"<DUP:TANDEM>"));
        assert!(!is_symbolic(b"ACG"));

        // Breakend
        assert!(is_breakend(b"A[chr2:100["));
        assert!(!is_breakend(b"ACG"));

        // Spanning deletion
        assert!(is_spanning_deletion(b"*"));
        assert!(!is_spanning_deletion(b"AC"));
    }

    /* ======== validate_vcf_file tests ============== */

    #[test]
    fn test_validate_samples_exist_single_sample() {
        let (tmp_vcf, sample_names) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result = validate_samples_exist(header, &[sample_names[0].clone()]);

        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_samples_exist_multiple_samples() {
        let (tmp_vcf, sample_names) = create_test_vcf(TestVcfConfig {
            num_samples: 3,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let to_validate = vec![sample_names[0].clone(), sample_names[2].clone()];

        let result = validate_samples_exist(header, &to_validate);

        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_samples_exist_missing_sample() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 2,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result = validate_samples_exist(header, &["nonexistent_sample".to_string()]);

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("nonexistent_sample"));
    }

    #[test]
    fn test_validate_samples_exist_multiple_missing() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result =
            validate_samples_exist(header, &["missing1".to_string(), "missing2".to_string()]);

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("missing1"));
        assert!(err_msg.contains("missing2"));
    }

    #[test]
    fn test_validate_samples_exist_case_sensitive() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        // Sample is "sample1" (lowercase)
        let result = validate_samples_exist(
            header,
            &["Sample1".to_string()], // Uppercase S
        );

        assert!(result.is_err(), "Sample names are case-sensitive");
    }

    #[test]
    fn test_validate_events_exist_single_event() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result = validate_events_exist(header, &["somatic".to_string()]);

        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_events_exist_multiple_events() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result =
            validate_events_exist(header, &["somatic".to_string(), "high_vaf".to_string()]);

        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_events_exist_missing_event() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result = validate_events_exist(header, &["nonexistent_event".to_string()]);

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("nonexistent_event"));
        assert!(err_msg.contains("PROB_"));
    }

    #[test]
    fn test_validate_events_exist_case_insensitive() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        // Event field is PROB_SOMATIC (uppercase)
        // User provides lowercase
        let result = validate_events_exist(header, &["somatic".to_string()]);

        assert!(result.is_ok(), "Should handle case conversion");
    }

    #[test]
    fn test_validate_events_exist_multiple_missing() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        let result =
            validate_events_exist(header, &["missing1".to_string(), "missing2".to_string()]);

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("missing1"));
        assert!(err_msg.contains("missing2"));
    }

    #[test]
    fn test_validate_events_exist_partial_match() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            num_samples: 1,
            ..Default::default()
        });

        let vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();
        let header = vcf.header();

        // One valid, one missing
        let result = validate_events_exist(
            header,
            &["somatic".to_string(), "missing_event".to_string()],
        );

        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("missing_event"));
        assert!(!err_msg.contains("somatic")); // Valid one not in error
    }

    #[test]
    fn test_validate_required_vcf_fields_all_valid() {
        let header =
            create_header_view(&[br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Test">"##]);

        assert!(validate_required_vcf_fields_msi(&header).is_ok());
    }

    #[test]
    fn test_validate_required_vcf_fields_missing_field() {
        let header = create_header_view(&[
            br##"##INFO=<ID=PROB_ARTIFACT,Number=A,Type=Float,Description="Test">"##,
        ]);

        let result = validate_required_vcf_fields_msi(&header);
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("AF") && err_msg.contains("missing"));
    }

    #[test]
    fn test_validate_required_vcf_fields_wrong_type() {
        let header = create_header_view(&[
            br##"##FORMAT=<ID=AF,Number=A,Type=Integer,Description="Test">"##,
        ]);

        let result = validate_required_vcf_fields_msi(&header);
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("incorrect type"));
    }

    #[test]
    fn test_validate_vcf_file_with_variant() {
        let (tmp_vcf, _) = create_test_vcf(TestVcfConfig {
            ref_allele: b"A",
            alt_alleles: vec![b"T"],
            af_values: Some(vec![0.5]),
            num_samples: 1,
            ..Default::default()
        });

        let mut vcf = bcf::Reader::from_path(tmp_vcf.path()).unwrap();

        let result = validate_vcf_file(&mut vcf);
        assert!(result.is_ok());
    }

    #[test]
    fn test_validate_vcf_file_empty() {
        let tmp = NamedTempFile::new().unwrap();
        let mut header = bcf::Header::new();
        header.push_record(br"##fileformat=VCFv4.2");
        header.push_record(br"##contig=<ID=chr1,length=1000000>");
        header.push_record(br##"##FORMAT=<ID=AF,Number=A,Type=Float,Description="AF">"##);
        header.push_sample(b"sample1");

        let writer = bcf::Writer::from_path(tmp.path(), &header, false, bcf::Format::Vcf).unwrap();
        drop(writer); // No records written

        let mut vcf = bcf::Reader::from_path(tmp.path()).unwrap();

        let result = validate_vcf_file(&mut vcf);
        assert!(result.is_err());

        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("empty") || err_msg.contains("no variant"));
    }

    #[test]
    fn test_validate_vcf_file_corrupted() {
        let tmp = NamedTempFile::new().unwrap();

        // Write valid header but invalid record line
        std::fs::write(
            tmp.path(),
            b"##fileformat=VCFv4.2\n\
            ##contig=<ID=chr1,length=1000>\n\
            ##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n\
            chr1\tINVALID_POS\t.\tA\tT\t.\t.\t.\tAF\t0.5\n",
        )
        .unwrap();

        let mut vcf = bcf::Reader::from_path(tmp.path()).unwrap();

        let result = validate_vcf_file(&mut vcf);
        assert!(result.is_err());

        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("read") || err_msg.contains("parse"));
    }
}
