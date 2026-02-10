
use super::*;

fn dummy_inputs() -> Stage5Inputs<'static> {
    let axes = Axes {
        tbi: vec![0.5],
        rci: vec![0.2],
        pds: vec![0.3],
        trs: vec![0.4],
        nsai: vec![0.1],
        iaa: vec![0.0],
        dfa: vec![0.0],
        cea: vec![0.0],
    };
    let drivers = vec![AxisDrivers {
        expressed_genes: 50,
        gene_entropy: 0.0,
        panel_entropy: 0.0,
        max_program_share: 0.0,
        tf_entropy: 0.0,
        stress_ratio: 0.0,
        dev_ratio: 0.0,
        iaa_raw: 0.0,
        dfa_raw: 0.0,
        cea_raw: 0.0,
        axis_variance: 0.0,
    }];
    let thresholds = ThresholdProfile::default_v1();
    Stage5Inputs {
        axes: Box::leak(Box::new(axes)),
        drivers: Box::leak(Box::new(drivers)),
        thresholds: Box::leak(Box::new(thresholds)),
        n_genes_mappable: Some(100),
        key_panel_coverage_median: Some(Box::leak(Box::new(vec![0.8]))),
        ambient_rna_risk: Some(Box::leak(Box::new(vec![false]))),
        key_panels_missing: Some(Box::leak(Box::new(vec![false]))),
        panel_nonzero_fraction: Some(Box::leak(Box::new(vec![0.5]))),
        axis_p90: Some([0.9, 0.1, 0.1]),
        scoring_mode: NuclearScoringMode::ImmuneAware,
    }
}

#[test]
fn test_composite_formula() {
    let inputs = dummy_inputs();
    let out = run_stage5(&inputs);
    let nps = clip01(0.45 * 0.5 + 0.35 * 0.2 - 0.20 * 0.3 - 0.20 * 0.4);
    let ci = clip01(0.55 * 0.4 + 0.45 * 0.3 - 0.15 * 0.5);
    let rls_base = clip01(0.45 * 0.5 + 0.35 * 0.2 - 0.25 * 0.3 - 0.15 * 0.1);
    assert!((out.scores.nps[0] - nps).abs() < 1e-6);
    assert!((out.scores.ci[0] - ci).abs() < 1e-6);
    assert!(out.scores.rls[0] <= rls_base + 1e-6);
    assert_eq!(out.scores.confidence_breakdown[0].len(), 4);
}

#[test]
fn test_confidence_degradation() {
    let axes = Axes {
        tbi: vec![0.5],
        rci: vec![0.2],
        pds: vec![0.3],
        trs: vec![0.4],
        nsai: vec![0.1],
        iaa: vec![0.0],
        dfa: vec![0.0],
        cea: vec![0.0],
    };
    let drivers = vec![AxisDrivers::default()];
    let thresholds = ThresholdProfile::default_v1();
    let inputs = Stage5Inputs {
        axes: &axes,
        drivers: &drivers,
        thresholds: &thresholds,
        n_genes_mappable: None,
        key_panel_coverage_median: None,
        ambient_rna_risk: None,
        key_panels_missing: None,
        panel_nonzero_fraction: None,
        axis_p90: None,
        scoring_mode: NuclearScoringMode::ImmuneAware,
    };
    let out = run_stage5(&inputs);
    assert_eq!(out.scores.confidence[0], 0.0);
}

#[test]
fn test_rls_floor_with_high_p90_iaa() {
    let inputs = dummy_inputs();
    let out = run_stage5(&inputs);
    assert!(out.scores.rls[0] >= 0.1);
}

#[test]
fn test_confidence_not_low_when_structure_high() {
    let mut inputs = dummy_inputs();
    let mut drivers = (*inputs.drivers).to_vec();
    drivers[0].axis_variance = 0.1;
    let axes = (*inputs.axes).clone();
    let thresholds = (*inputs.thresholds).clone();
    let inputs = Stage5Inputs {
        axes: &axes,
        drivers: &drivers,
        thresholds: &thresholds,
        n_genes_mappable: Some(100),
        key_panel_coverage_median: Some(Box::leak(Box::new(vec![0.9]))),
        ambient_rna_risk: Some(Box::leak(Box::new(vec![false]))),
        key_panels_missing: Some(Box::leak(Box::new(vec![false]))),
        panel_nonzero_fraction: Some(Box::leak(Box::new(vec![0.5]))),
        axis_p90: Some([0.9, 0.2, 0.2]),
        scoring_mode: NuclearScoringMode::ImmuneAware,
    };
    let out = run_stage5(&inputs);
    assert!(out.scores.confidence[0] >= 0.2);
}

#[test]
fn test_driver_ordering() {
    let inputs = dummy_inputs();
    let out = run_stage5(&inputs);
    let drivers = &out.drivers.nps[0];
    for w in drivers.windows(2) {
        let a = w[0].1.abs();
        let b = w[1].1.abs();
        assert!(a >= b || (a - b).abs() < 1e-12);
    }
}

#[test]
fn test_determinism_bits() {
    let inputs = dummy_inputs();
    let out_a = run_stage5(&inputs);
    let out_b = run_stage5(&inputs);
    assert_eq!(out_a.scores.nps[0].to_bits(), out_b.scores.nps[0].to_bits());
    assert_eq!(out_a.scores.ci[0].to_bits(), out_b.scores.ci[0].to_bits());
    assert_eq!(out_a.scores.rls[0].to_bits(), out_b.scores.rls[0].to_bits());
}
