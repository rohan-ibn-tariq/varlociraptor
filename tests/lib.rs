#[macro_use]
extern crate paste;
#[macro_use]
extern crate lazy_static;

use anyhow::Context;
use anyhow::Result;
use std::sync::Mutex;
use std::{fs, path::Path, path::PathBuf};

use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use rust_htslib::bcf::{self, Read, Reader};
use varlociraptor::utils::bcf_utils;
use varlociraptor::{testcase, testcase_should_panic};

testcase!(test01, exact, fast);
testcase!(test02, exact, fast);
testcase!(test03, exact, fast);
testcase!(test04, exact, fast);
testcase!(test05, exact, fast);
testcase!(test06, exact, fast);
testcase!(test07, exact, fast);
testcase!(test08, exact, fast);
testcase!(test09, exact, fast);
testcase!(test10, exact, fast);
testcase!(test11, exact, fast);
testcase!(test12, exact, fast);
testcase!(test13, exact, fast);
testcase!(test14, exact, fast);
testcase!(test15, exact, fast);
testcase!(test16, exact, fast);
testcase!(test17, exact, fast);
testcase!(test18, exact, fast);
testcase!(test19, exact, fast);
testcase!(test20, exact, fast);
// skip the next test because this insertion cannot currently be resolved properly
// TODO find a way to fix this.
// testcase!(test21, exact, fast);
testcase!(test22, exact, fast);
testcase!(test23, exact, fast);
testcase!(test24, exact, fast);
testcase!(test25, exact, fast);
testcase!(test26, exact, fast);
testcase!(test27, exact, fast);
testcase!(test28, exact, fast);
testcase!(test29, exact, fast);
testcase!(test30, exact, fast);
testcase!(test31, exact, fast);
testcase!(test32, exact, fast);
testcase!(test33, exact, fast);
testcase!(test34, exact, fast);
testcase!(test36, exact, fast);
testcase!(test37, exact, fast);
// Skip this test. It is most likely a strand bias artifact, which is correctly recognized.
// However, there are also very few reads with nonstandard orientation, which are on the other
// strand.
//testcase!(test38, exact, fast);
testcase!(test39, exact, fast);
testcase!(test40, exact, fast);
testcase!(test41, exact, fast);
testcase!(test42, exact, fast);
testcase!(test43, exact, fast);
// Fast mode fails here, because there is a read with two insertions against
// the alt allele. This is very unlikely to happen, but it happens here.
// In the exact mode, there are various paths around this alignment which rescue
// the alt allele probability. With fast mode, these are missed, making the
// probability artificially small. That leads to Varlociraptor evaluating the
// locus to be heterozygous although it is homozygous in reality.
testcase!(test44, exact);
testcase!(test45, exact, fast);

testcase!(test47, exact, fast);
testcase!(test48, exact, fast);
testcase!(test49, exact);
testcase!(test50, exact);
testcase!(test51, exact);
testcase!(test52, exact, fast);
testcase!(test53, exact, fast);
testcase!(test54, exact, fast);
testcase!(test55, exact, fast);
testcase!(test57, exact);
testcase!(test58, exact);
testcase!(test59, exact);
testcase!(test60, exact);
testcase!(test61, exact);
testcase!(test62, exact);
testcase!(test63, exact);
testcase!(test64, exact);
testcase!(test65, exact);
testcase!(test69, exact);
testcase!(test70, exact);
testcase!(test71, exact);
testcase!(test72, exact);
testcase!(test74, exact);
testcase!(test75, exact);
testcase!(test77, exact);
testcase!(test78, exact);
testcase!(test79, exact);

// Skip test_giab_01 for now as the call seems to be correct.
// TODO try to find out what is wrong in the GIAB callset at that location.
testcase!(test_giab_01, exact);
testcase!(test_giab_02, exact);
testcase!(test_giab_03, exact);
testcase!(test_giab_04, exact);

testcase!(test_giab_05, exact);

testcase!(test_giab_06, exact);
testcase!(test_giab_07, exact);
testcase!(test_giab_08, exact);
testcase!(test_giab_09, exact);
testcase!(test_giab_10, exact);
testcase!(test_giab_11, exact);
testcase!(test_giab_12, exact);
testcase!(test_giab_13, exact);
// Skip test_giab_14. It is just bad luck that the reads here look like a homopolymer artifact although the variant is real.
// See testcase.yaml for details.
//testcase!(test_giab_14, exact);
testcase!(test_giab_15, exact);
testcase!(test_giab_16, exact);
testcase!(test_giab_17, exact);
testcase!(test_giab_18, exact);
testcase!(test_giab_19, exact);
testcase!(test_giab_20, exact);
testcase!(test_giab_21, exact);
testcase!(test_giab_22, exact);
testcase!(test_giab_23, exact);
// Skip test_giab_24. It is simply an unlucky combination of homology artifacts that looks strong.
// At least our probability is weaker than Freebayes's already.
// testcase!(test_giab_24, exact);
testcase!(test_giab_25, exact);
testcase!(test_giab_26, exact);
testcase!(test_giab_27, exact);
// Skip test_giab_28. It is simply an unlucky combination of homology artifacts that looks strongly like an artifact.
//testcase!(test_giab_28, exact);
testcase!(test_giab_29, exact);
testcase!(test_giab_30, exact);
testcase!(test_giab_31, exact);
testcase!(test_giab_32, exact);
testcase!(test_giab_33, exact);
testcase!(test_giab_34, exact);
testcase!(test_giab_35, exact);
testcase!(test_mapq_meth, exact);

testcase!(test_pcr_homopolymer_error1, exact);
testcase!(test_pcr_homopolymer_error2, exact);
testcase!(test_pcr_homopolymer_error3, exact);

testcase!(test_mendelian_prior, exact);
testcase!(pattern_too_long, exact, fast);
testcase!(test_long_pattern, exact, fast);
testcase!(test_contig_universe, exact, fast);
testcase!(test_expressions, exact);
testcase!(omit_sb, exact);
testcase!(test_panel_overlap, exact);
testcase!(test_panel_unknown_orientation_bias, exact);
testcase!(issue_154, exact, fast);
testcase!(test_low_cov_vaf, exact);
testcase_should_panic!(test_overlapping_events, exact);

testcase!(test_l2fc, exact, fast);
testcase!(test_cmp, exact, fast);

testcase!(test_nanopore_01, homopolymer);
testcase!(test_nanopore_02, homopolymer);
testcase!(test_nanopore_03, homopolymer);
testcase!(test_nanopore_04, homopolymer);
testcase!(test_nanopore_05, homopolymer);

testcase!(test_haplotype_absent, exact);
testcase!(test_haplotype_present, exact);
testcase!(test_haplotype_singleton, exact);

testcase!(test_alt_locus_bias_01, exact);
testcase!(test_alt_locus_bias_02, exact);
testcase!(test_uzuner_fn_mnv, exact);
testcase!(test_uzuner_fp_mnv1, exact);

testcase!(test_prinz_af_scan, exact);
testcase!(test_prinz_call_meth_1, exact);
testcase!(test_prinz_call_meth_2, exact);
testcase!(test_prinz_pacbio_zero, exact);

testcase!(test_imprecise_fusion, exact);
testcase!(test_imprecise_fusion_absent, exact);

testcase!(test_uzuner_clonal_1, exact);
testcase!(test_uzuner_clonal_2, exact);
testcase!(test_uzuner_clonal_3, exact);
testcase!(test_uzuner_fp_snv_on_ins, exact);
testcase!(test_false_negative_indel_call, exact);
testcase!(test_hiv_vaf_higher_than_expected, exact);
testcase!(test_uzuner_only_N, exact);
testcase!(test_moelder_floatisnan, exact);
testcase!(test_alt_locus_mapq_only, exact);
testcase!(test_single_value_afd, exact);

fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}

fn cleanup_file(f: &str) {
    if Path::new(f).exists() {
        fs::remove_file(f).unwrap();
    }
}

fn control_fdr(
    test: &str,
    events: &[&str],
    alpha: f64,
    local: bool,
    smart: bool,
    smart_retain_artifacts: bool,
    vartype: Option<&varlociraptor::variants::model::VariantType>,
) {
    let basedir = basedir(test);
    let output = format!("{}/calls.filtered.bcf", basedir);
    cleanup_file(&output);
    let event_strs: Vec<varlociraptor::SimpleEvent> = events
        .iter()
        .map(|&event_str| varlociraptor::SimpleEvent {
            name: event_str.to_owned(),
        })
        .collect();
    varlociraptor::filtration::fdr::control_fdr(
        &format!("{}/calls.matched.bcf", basedir),
        Some(&output),
        &event_strs,
        vartype,
        LogProb::from(Prob(alpha)),
        local,
        smart,
        smart_retain_artifacts,
    )
    .unwrap();
}

fn assert_call_number(test: &str, expected_calls: usize) {
    let basedir = basedir(test);

    let mut reader = bcf::Reader::from_path(format!("{}/calls.filtered.bcf", basedir)).unwrap();

    let calls = reader.records().map(|r| r.unwrap()).collect_vec();

    let ok = if expected_calls > 50 {
        // allow one more or less, in order to be robust to numeric fluctuations
        (calls.len() as i32 - expected_calls as i32).abs() <= 1
    } else {
        calls.len() == expected_calls
    };

    assert!(
        ok,
        "unexpected number of calls ({} vs {})",
        calls.len(),
        expected_calls
    );
}

#[test]
fn test_fdr_control1() {
    control_fdr(
        "test_fdr_ev_1",
        &["SOMATIC"],
        0.05,
        false,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    //assert_call_number("test_fdr_ev_1", 974);
}

#[test]
fn test_fdr_control2() {
    control_fdr(
        "test_fdr_ev_2",
        &["SOMATIC"],
        0.05,
        false,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_ev_2", 985);
}

/// same test, but low alpha
#[test]
fn test_fdr_control3() {
    control_fdr(
        "test_fdr_ev_3",
        &["ABSENT"],
        0.001,
        false,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_ev_3", 0);
}

#[test]
fn test_fdr_control4() {
    control_fdr(
        "test_fdr_ev_4",
        &["SOMATIC_TUMOR"],
        0.05,
        false,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_ev_4", 0);
}

#[test]
fn test_fdr_control_local1() {
    control_fdr(
        "test_fdr_local1",
        &["SOMATIC"],
        0.05,
        true,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_local1", 0);
}

#[test]
fn test_fdr_control_local2() {
    control_fdr(
        "test_fdr_local2",
        &["SOMATIC"],
        0.25,
        true,
        false,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_local2", 1);
}

#[test]
fn test_fdr_control_local2_smart() {
    control_fdr(
        "test_fdr_local2_smart",
        &["SOMATIC"],
        0.08,
        true,
        true,
        false,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_local2_smart", 1);
}

#[test]
fn test_fdr_control_local2_smart_retain_artifacts() {
    control_fdr(
        "test_fdr_local2_smart",
        &["SOMATIC"],
        0.08,
        true,
        true,
        true,
        Some(&varlociraptor::variants::model::VariantType::Deletion(
            Some(1..30),
        )),
    );
    assert_call_number("test_fdr_local2_smart", 1);
}

#[test]
fn test_fdr_control_local3() {
    control_fdr(
        "test_fdr_local3",
        &["GERMLINE", "SOMATIC_TUMOR_LOW"],
        0.05,
        true,
        false,
        false,
        None,
    );
    assert_call_number("test_fdr_local3", 0);
}

// TODO enable this test again once https://github.com/samtools/bcftools/issues/874 is truly fixed upstream
// Then, also encode SVLEN as negative again for deletions.
//#[test]
// fn test_fdr_control5() {
//     control_fdr(
//         "test_fdr_control_out_of_bounds",
//         &["PRESENT"],
//         0.05,
//         false,
//         Some(&varlociraptor::variants::model::VariantType::Deletion(
//             Some(1..30),
//         )),
//     );
// }

//####################################################################################################################################################
// Tests for methylation candidates
//####################################################################################################################################################

fn control_meth_candidates(test: &str) -> Result<()> {
    let basedir = basedir(test);
    let output = format!("{}/candidates.bcf", basedir);
    cleanup_file(&output);
    varlociraptor::candidates::methylation::find_candidates(
        PathBuf::from(format!("{}/genome.fasta", basedir)),
        vec![varlociraptor::candidates::methylation::MethylationMotif::CG],
        Some(PathBuf::from(output)),
    )
    .with_context(|| "error computing methylation candidates".to_string())?;
    Ok(())
}

fn assert_candidates_number(test: &str, expected_calls: usize) -> Result<()> {
    let basedir = basedir(test);

    let mut reader = Reader::from_path(format!("{}/candidates.bcf", basedir))
        .with_context(|| "error reading BCF file".to_string())?;
    let calls = reader.records().map(|record| record.unwrap()).collect_vec();

    let ok = calls.len() == expected_calls;

    assert!(
        ok,
        "unexpected number of calls ({} vs {})",
        calls.len(),
        expected_calls
    );
    Ok(())
}

#[test]
fn test_meth_candidates1() -> Result<()> {
    control_meth_candidates("test_meth_ev_1")?;
    assert_candidates_number("test_meth_ev_1", 6)?;
    Ok(())
}

//####################################################################################################################################################
// Tests for Preprocessing: Microsatellite Instability
//####################################################################################################################################################

fn run_msi_preprocess_with_fields(
    test: &str,
    propagate_info_fields: Vec<String>,
) -> Result<PathBuf> {
    let basedir = basedir(test);
    let output = format!("{}/output.vcf", basedir);
    cleanup_file(&output);

    varlociraptor::preprocess_ms_candidates(varlociraptor::PreprocessMSIConfig {
        microsatellite_bed: PathBuf::from(format!("{}/regions.bed", basedir)),
        candidate_vcf: PathBuf::from(format!("{}/input.vcf", basedir)),
        output: Some(PathBuf::from(&output)),
        propagate_info_fields,
    })?;

    Ok(PathBuf::from(output))
}

fn run_msi_preprocess(test: &str) -> Result<PathBuf> {
    run_msi_preprocess_with_fields(test, vec![])
}

fn run_msi_preprocess_expect_err(test: &str) -> anyhow::Error {
    run_msi_preprocess_with_fields(test, vec![]).unwrap_err()
}

#[test]
fn test_preprocess_msi_basic_annotation() -> Result<()> {
    // VCF: 1 perfect CAG insertion at chr1:95
    // BED: chr1:94-124 10xCAG
    // Expected: 1 annotated record, 0 dummies
    let output = run_msi_preprocess("test_preprocess_msi_basic")?;
    let records = bcf_utils::read_bcf_records(&output)?;

    assert_eq!(records.len(), 1, "Expected 1 record total");
    assert!(
        bcf_utils::record_has_info_string(&records[0], b"REGION_ID"),
        "Real variant should have REGION_ID"
    );
    assert!(
        !bcf_utils::record_has_info_flag(&records[0], b"MSI_DUMMY"),
        "Should NOT have MSI_DUMMY flag"
    );
    assert_eq!(
        bcf_utils::get_info_strings(&records[0], b"REGION_ID")
            .and_then(|v| v.into_iter().next())
            .unwrap(),
        "chr1:94-124",
        "REGION_ID should match the BED region"
    );

    Ok(())
}

#[test]
fn test_preprocess_msi_dummy_injection() -> Result<()> {
    let output = run_msi_preprocess("test_preprocess_msi_dummy")?;
    let records = bcf_utils::read_bcf_records(&output)?;

    assert_eq!(records.len(), 2, "Expected SNV + dummy = 2 records");

    let dummies: Vec<_> = records
        .iter()
        .filter(|r| bcf_utils::record_has_info_flag(r, b"MSI_DUMMY"))
        .collect();
    let non_dummies: Vec<_> = records
        .iter()
        .filter(|r| !bcf_utils::record_has_info_flag(r, b"MSI_DUMMY"))
        .collect();

    assert_eq!(dummies.len(), 1, "Expected exactly 1 dummy");
    assert_eq!(non_dummies.len(), 1, "Expected exactly 1 original SNV");

    assert!(
        bcf_utils::record_has_info_string(dummies[0], b"REGION_ID"),
        "Dummy should have REGION_ID"
    );
    assert_eq!(
        bcf_utils::get_info_strings(dummies[0], b"REGION_ID")
            .and_then(|v| v.into_iter().next())
            .unwrap(),
        "chr1:94-124"
    );
    assert_eq!(dummies[0].pos(), 96, "Dummy pos = 94 + 3 - 1 = 96");

    assert!(
        !bcf_utils::record_has_info_string(non_dummies[0], b"REGION_ID"),
        "SNV should NOT have REGION_ID"
    );
    assert!(
        !bcf_utils::record_has_info_flag(non_dummies[0], b"MSI_DUMMY"),
        "SNV should NOT have MSI_DUMMY flag"
    );

    Ok(())
}

#[test]
fn test_preprocess_msi_imperfect_indel_gets_dummy() -> Result<()> {
    let output = run_msi_preprocess("test_preprocess_msi_imperfect")?;
    let records = bcf_utils::read_bcf_records(&output)?;

    assert_eq!(
        records.len(),
        2,
        "Expected imperfect indel + dummy = 2 records"
    );

    let dummies: Vec<_> = records
        .iter()
        .filter(|r| bcf_utils::record_has_info_flag(r, b"MSI_DUMMY"))
        .collect();
    let non_dummies: Vec<_> = records
        .iter()
        .filter(|r| !bcf_utils::record_has_info_flag(r, b"MSI_DUMMY"))
        .collect();

    assert_eq!(dummies.len(), 1, "Expected 1 dummy");
    assert_eq!(non_dummies.len(), 1, "Expected 1 imperfect indel");

    assert!(
        bcf_utils::record_has_info_string(dummies[0], b"REGION_ID"),
        "Dummy should have REGION_ID"
    );
    assert_eq!(
        bcf_utils::get_info_strings(dummies[0], b"REGION_ID")
            .and_then(|v| v.into_iter().next())
            .unwrap(),
        "chr1:94-124"
    );
    assert_eq!(dummies[0].pos(), 96, "Dummy pos = 94 + 3 - 1 = 96");

    assert!(
        !bcf_utils::record_has_info_string(non_dummies[0], b"REGION_ID"),
        "Imperfect indel should NOT have REGION_ID"
    );
    assert!(
        !bcf_utils::record_has_info_flag(non_dummies[0], b"MSI_DUMMY"),
        "Imperfect indel should NOT have MSI_DUMMY flag"
    );

    Ok(())
}

#[test]
fn test_preprocess_msi_multi_region() -> Result<()> {
    let output = run_msi_preprocess("test_preprocess_msi_multi")?;
    let records = bcf_utils::read_bcf_records(&output)?;

    assert_eq!(
        records.len(),
        2,
        "Expected annotated variant + dummy = 2 records"
    );

    let annotated: Vec<_> = records
        .iter()
        .filter(|r| {
            bcf_utils::record_has_info_string(r, b"REGION_ID")
                && !bcf_utils::record_has_info_flag(r, b"MSI_DUMMY")
        })
        .collect();
    let dummies: Vec<_> = records
        .iter()
        .filter(|r| bcf_utils::record_has_info_flag(r, b"MSI_DUMMY"))
        .collect();

    assert_eq!(annotated.len(), 1, "Expected 1 annotated real variant");
    assert_eq!(dummies.len(), 1, "Expected 1 dummy for second region");

    assert_eq!(
        bcf_utils::get_info_strings(annotated[0], b"REGION_ID")
            .and_then(|v| v.into_iter().next())
            .unwrap(),
        "chr1:94-124"
    );
    assert_eq!(
        bcf_utils::get_info_strings(dummies[0], b"REGION_ID")
            .and_then(|v| v.into_iter().next())
            .unwrap(),
        "chr1:200-230"
    );
    assert_eq!(dummies[0].pos(), 202, "Dummy pos = 200 + 3 - 1 = 202");

    Ok(())
}

#[test]
fn test_preprocess_msi_overlapping_bed_errors() -> Result<()> {
    let err = run_msi_preprocess_expect_err("test_preprocess_msi_overlap_error");
    assert!(
        err.to_string().contains("overlaps multiple BED regions"),
        "Error message should mention overlapping regions, got: {}",
        err
    );
    Ok(())
}

#[test]
fn test_preprocess_msi_propagates_info_fields() -> Result<()> {
    let output = run_msi_preprocess_with_fields(
        "test_preprocess_msi_propagate",
        vec!["COSMIC_ID".to_string()],
    )?;
    let records = bcf_utils::read_bcf_records(&output)?;

    assert_eq!(records.len(), 1, "Expected 1 annotated record");
    assert!(
        bcf_utils::record_has_info_string(&records[0], b"COSMIC_ID"),
        "COSMIC_ID should be propagated"
    );
    assert_eq!(
        bcf_utils::get_info_strings(&records[0], b"COSMIC_ID")
            .and_then(|v| v.into_iter().next())
            .unwrap(),
        "COSM123"
    );
    Ok(())
}

//####################################################################################################################################################
// Tests for Call: Microsatellite Instability
//####################################################################################################################################################

/// Allowed error when checking that several probabilities add up to 1.0.
/// Adding many floats might accumulate small rounding errors, so this is looser
/// than EXACT_MATCH_EPSILON below.
const PROB_SUM_EPSILON: f64 = 1e-6;

/// Allowed rounding error when comparing the same number as written to two
/// different output files (e.g. the TSV and the JSON plot). This tolerance is
/// small because it only needs to cover that rounding, not any real computation
/// difference.
const EXACT_MATCH_EPSILON: f64 = 1e-9;

fn default_msi_call_config(basedir: &str) -> varlociraptor::MSIConfig {
    varlociraptor::MSIConfig {
        calls: PathBuf::from(format!("{}/calls.vcf", basedir)),
        sample: "tumor".to_string(),
        events: vec!["somatic".to_string()],
        msi_threshold: 3.5,
        af_thresholds: vec![1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.02, 0.0],
        distribution_af: 0.05,
        windowed_af: 0.05,
        sliding_window: 1_000_000,
        threads: Some(1),
        plot_distribution: None,
        plot_pseudotime: None,
        plot_heatmap: None,
        data_distribution: None,
        data_pseudotime: None,
        data_heatmap: None,
        is_phred: false,
    }
}

fn run_msi_call(
    basedir: &str,
    configure: impl FnOnce(&mut varlociraptor::MSIConfig),
) -> Result<()> {
    let mut config = default_msi_call_config(&basedir);
    configure(&mut config);
    config.set_defaults()?;
    config.validate()?;
    varlociraptor::call_msi(config)
}

/* ====== 1. Distribution ====== */

#[test]
fn test_call_msi_distribution_basic() -> Result<()> {
    let basedir = basedir("test_call_msi_distribution_basic");
    let data = format!("{}/dist.tsv", basedir);
    let plot = format!("{}/dist.vl.json", basedir);
    cleanup_file(&data);
    cleanup_file(&plot);

    run_msi_call(&basedir, |c| {
        c.data_distribution = Some(PathBuf::from(&data));
        c.plot_distribution = Some(PathBuf::from(&plot));
    })?;

    /****** TSV Checks ************/
    let content = fs::read_to_string(&data)?;
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 4, "header + k=0,1,2 rows for 2 regions");

    let tsv_probs: Vec<f64> = lines[1..]
        .iter()
        .map(|l| l.split('\t').last().unwrap().parse::<f64>().unwrap())
        .collect();
    let tsv_sum: f64 = tsv_probs.iter().sum();
    assert!(
        (tsv_sum - 1.0).abs() < PROB_SUM_EPSILON,
        "TSV distribution must sum to 1.0"
    );

    /****** Plot Checks ************/
    let plot_value: serde_json::Value = serde_json::from_str(&fs::read_to_string(&plot)?)?;
    let data_values = plot_value["data"]["values"]
        .as_array()
        .expect("plot should have a data.values array");
    assert_eq!(data_values.len(), 3, "one plot point per k=0,1,2");

    let plot_probs: Vec<f64> = data_values
        .iter()
        .map(|v| {
            v["probability"]
                .as_f64()
                .expect("probability should be numeric")
        })
        .collect();
    let plot_sum: f64 = plot_probs.iter().sum();
    assert!(
        (plot_sum - 1.0).abs() < PROB_SUM_EPSILON,
        "plot distribution must sum to 1.0"
    );

    /****** Cross-check: TSV and plot must agree exactly, for every k ************/
    for k in 0..3 {
        assert!(
            (tsv_probs[k] - plot_probs[k]).abs() < EXACT_MATCH_EPSILON,
            "TSV/plot mismatch at k={}: tsv={}, plot={}",
            k,
            tsv_probs[k],
            plot_probs[k]
        );
    }

    Ok(())
}

/* ====== 2. Pseudotime ====== */

#[test]
fn test_call_msi_pseudotime_basic() -> Result<()> {
    let basedir = basedir("test_call_msi_pseudotime_basic");
    let data = format!("{}/pseudo.tsv", basedir);
    let plot = format!("{}/pseudo.vl.json", basedir);
    cleanup_file(&data);
    cleanup_file(&plot);

    run_msi_call(&basedir, |c| {
        c.data_pseudotime = Some(PathBuf::from(&data));
        c.plot_pseudotime = Some(PathBuf::from(&plot));
    })?;

    /****** TSV Checks ************/
    let content = fs::read_to_string(&data)?;
    let rows: Vec<Vec<&str>> = content.lines().skip(1).map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows.len(), 9, "one row per AF threshold in the thresholds list");

    // af=1.0: no variant (AF 0.9, 0.6, 0.3) passes -> empty DP -> k_map=0
    let row_af_1_0 = rows.iter().find(|r| r[1] == "1.00").expect("af=1.0 row");
    assert_eq!(row_af_1_0[3], "0", "k_map=0 when no regions pass the filter");
    assert_eq!(row_af_1_0[2], "0.00");
    assert_eq!(row_af_1_0[5], "MSS");

    // af=0.8: only the AF=0.9 region passes -> single region, p_stable=0.4 -> k_map=1
    let row_af_0_8 = rows.iter().find(|r| r[1] == "0.80").expect("af=0.8 row");
    assert_eq!(row_af_0_8[3], "1", "single region with p_stable=0.4 -> k_map=1");
    assert_eq!(row_af_0_8[2], "33.33");

    // af=0.6: AF=0.9 and AF=0.6 regions pass -> p_stable=[0.4,0.45] -> DP favors k=1 (0.49 vs 0.18, 0.33)
    let row_af_0_6 = rows.iter().find(|r| r[1] == "0.60").expect("af=0.6 row");
    assert_eq!(row_af_0_6[3], "1", "two-region DP still favors k_map=1 here");
    assert_eq!(row_af_0_6[2], "33.33");

    // af=0.0: all three regions pass -> p_stable=[0.4,0.45,0.6] -> DP favors k=2 (0.394 vs 0.366/0.132/0.108)
    let row_af_0_0 = rows.iter().find(|r| r[1] == "0.00").expect("af=0.0 row");
    assert_eq!(row_af_0_0[3], "2", "three-region DP shifts the mode to k_map=2");
    assert_eq!(row_af_0_0[2], "66.67");
    assert_eq!(row_af_0_0[5], "MSI-High");

    /****** Plot Checks ************/
    let plot_value: serde_json::Value = serde_json::from_str(&fs::read_to_string(&plot)?)?;
    let data_values = plot_value["data"]["values"]
        .as_array()
        .expect("plot should have a data.values array");
    assert_eq!(data_values.len(), 9, "one plot point per AF threshold");

    let plot_af_0_0 = data_values
        .iter()
        .find(|v| (v["af_threshold"].as_f64().unwrap() - 0.0).abs() < EXACT_MATCH_EPSILON)
        .expect("plot should have an af_threshold=0.0 point");
    let plot_msi_score = plot_af_0_0["msi_score"].as_f64().expect("msi_score should be numeric");
    assert!((plot_msi_score - 66.67).abs() < 0.01, "plot msi_score should match TSV at af=0.0");

    Ok(())
}