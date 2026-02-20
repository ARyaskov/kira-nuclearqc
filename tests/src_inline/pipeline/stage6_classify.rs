use super::*;

struct TestInputs {
    tbi: Vec<f32>,
    rci: Vec<f32>,
    pds: Vec<f32>,
    trs: Vec<f32>,
    nsai: Vec<f32>,
    iaa: Vec<f32>,
    dfa: Vec<f32>,
    cea: Vec<f32>,
    rss: Vec<f32>,
    drbi: Vec<f32>,
    cci: Vec<f32>,
    trci: Vec<f32>,
    scores: CompositeScores,
    drivers: Vec<AxisDrivers>,
    thresholds: ThresholdProfile,
    key_panel_coverage_median: Option<Vec<f32>>,
    key_panels_missing: Option<Vec<bool>>,
    sum_tf_panels: Option<Vec<f32>>,
    ambient_rna_risk: Option<Vec<bool>>,
    proliferation_program_share: Option<Vec<f32>>,
    program_sum: Option<Vec<f32>>,
}

impl TestInputs {
    fn as_inputs(&self) -> Stage6Inputs<'_> {
        Stage6Inputs {
            tbi: &self.tbi,
            rci: &self.rci,
            pds: &self.pds,
            trs: &self.trs,
            nsai: &self.nsai,
            iaa: &self.iaa,
            dfa: &self.dfa,
            cea: &self.cea,
            rss: &self.rss,
            drbi: &self.drbi,
            cci: &self.cci,
            trci: &self.trci,
            scores: &self.scores,
            drivers: &self.drivers,
            thresholds: &self.thresholds,
            scoring_mode: NuclearScoringMode::ImmuneAware,
            key_panel_coverage_median: self.key_panel_coverage_median.as_deref(),
            key_panels_missing: self.key_panels_missing.as_deref(),
            sum_tf_panels: self.sum_tf_panels.as_deref(),
            ambient_rna_risk: self.ambient_rna_risk.as_deref(),
            proliferation_program_share: self.proliferation_program_share.as_deref(),
            program_sum: self.program_sum.as_deref(),
        }
    }
}

fn base_inputs() -> TestInputs {
    TestInputs {
        tbi: vec![0.2],
        rci: vec![0.2],
        pds: vec![0.2],
        trs: vec![0.2],
        nsai: vec![0.2],
        iaa: vec![0.2],
        dfa: vec![0.2],
        cea: vec![0.2],
        rss: vec![0.2],
        drbi: vec![0.2],
        cci: vec![0.2],
        trci: vec![0.2],
        scores: CompositeScores {
            nps: vec![0.2],
            ci: vec![0.2],
            rls: vec![0.2],
            confidence: vec![0.5],
            confidence_breakdown: vec![[0.0, 0.0, 0.0, 0.0]],
        },
        drivers: vec![AxisDrivers {
            expressed_genes: 50,
            gene_entropy: 0.2,
            panel_entropy: 0.2,
            max_program_share: 0.2,
            tf_entropy: 0.2,
            stress_ratio: 0.2,
            dev_ratio: 0.2,
            iaa_raw: 0.0,
            dfa_raw: 0.0,
            cea_raw: 0.0,
            axis_variance: 0.0,
        }],
        thresholds: ThresholdProfile::default_v1(),
        key_panel_coverage_median: None,
        key_panels_missing: None,
        sum_tf_panels: None,
        ambient_rna_risk: None,
        proliferation_program_share: None,
        program_sum: None,
    }
}

#[test]
fn test_first_match_collapsed() {
    let mut inputs = base_inputs();
    inputs.drivers[0].expressed_genes = 0;
    let out = run_stage6(&inputs.as_inputs());
    assert_eq!(out[0].regime, NuclearRegime::TranscriptionallyCollapsed);
}

#[test]
fn test_rigid_deg() {
    let mut inputs = base_inputs();
    inputs.trs[0] = 0.8;
    inputs.nsai[0] = 0.6;
    inputs.rci[0] = 0.3;
    let out = run_stage6(&inputs.as_inputs());
    assert_eq!(out[0].regime, NuclearRegime::RigidDegenerative);
}

#[test]
fn test_committed_state() {
    let mut inputs = base_inputs();
    inputs.trs[0] = 0.72;
    inputs.pds[0] = 0.62;
    inputs.tbi[0] = 0.40;
    inputs.nsai[0] = 0.30;
    let out = run_stage6(&inputs.as_inputs());
    assert_eq!(out[0].regime, NuclearRegime::CommittedState);
}

#[test]
fn test_stress_adaptive() {
    let mut inputs = base_inputs();
    inputs.nsai[0] = 0.70;
    inputs.rci[0] = 0.40;
    inputs.tbi[0] = 0.36;
    let out = run_stage6(&inputs.as_inputs());
    assert_eq!(out[0].regime, NuclearRegime::StressAdaptive);
}

#[test]
fn test_plastic_adaptive() {
    let mut inputs = base_inputs();
    inputs.scores.nps[0] = 0.65;
    inputs.trs[0] = 0.40;
    inputs.pds[0] = 0.40;
    let out = run_stage6(&inputs.as_inputs());
    assert_eq!(out[0].regime, NuclearRegime::PlasticAdaptive);
}

#[test]
fn test_transient_adaptive() {
    let mut inputs = base_inputs();
    inputs.scores.nps[0] = 0.50;
    inputs.trs[0] = 0.50;
    inputs.pds[0] = 0.60;
    inputs.iaa[0] = 0.40;
    let out = run_stage6(&inputs.as_inputs());
    assert_eq!(out[0].regime, NuclearRegime::TransientAdaptive);
}

#[test]
fn test_unclassified() {
    let inputs = base_inputs();
    let out = run_stage6(&inputs.as_inputs());
    assert_eq!(out[0].regime, NuclearRegime::Unclassified);
}

#[test]
fn test_flags() {
    let mut inputs = base_inputs();
    inputs.drivers[0].expressed_genes = 1;
    inputs.scores.confidence[0] = 0.2;
    inputs.key_panel_coverage_median = Some(vec![0.2]);
    inputs.key_panels_missing = Some(vec![true]);
    inputs.sum_tf_panels = Some(vec![0.0]);
    inputs.ambient_rna_risk = Some(vec![true]);
    inputs.proliferation_program_share = Some(vec![0.9]);
    inputs.pds[0] = 0.9;
    inputs.nsai[0] = 0.9;
    inputs.rss[0] = 0.8;
    inputs.drbi[0] = 0.9;
    inputs.cci[0] = 0.8;
    inputs.trci[0] = 0.8;

    let out = run_stage6(&inputs.as_inputs());
    let flags = &out[0].flags;
    assert!(flags.contains(&Flag::LowExprGenes));
    assert!(flags.contains(&Flag::LowPanelCoverage));
    assert!(flags.contains(&Flag::MissingKeyPanels));
    assert!(flags.contains(&Flag::HighProgramDominance));
    assert!(flags.contains(&Flag::HighStressBias));
    assert!(flags.contains(&Flag::LowTfSignal));
    assert!(flags.contains(&Flag::AmbientRnaRisk));
    assert!(flags.contains(&Flag::CellCycleConfounder));
    assert!(flags.contains(&Flag::LowConfidence));
    assert!(flags.contains(&Flag::HighReplicationStress));
    assert!(flags.contains(&Flag::HrDominantRepair));
    assert!(flags.contains(&Flag::ChromatinHypercompact));
    assert!(flags.contains(&Flag::HighTrConflict));
}

#[test]
fn test_ddr_repair_bias_flags() {
    let mut inputs = base_inputs();
    inputs.drbi[0] = 0.2;
    let out = run_stage6(&inputs.as_inputs());
    assert!(out[0].flags.contains(&Flag::NhejDominantRepair));

    let mut inputs = base_inputs();
    inputs.drbi[0] = 0.9;
    let out = run_stage6(&inputs.as_inputs());
    assert!(out[0].flags.contains(&Flag::HrDominantRepair));
}

#[test]
fn test_determinism() {
    let inputs = base_inputs();
    let a = run_stage6(&inputs.as_inputs());
    let b = run_stage6(&inputs.as_inputs());
    assert_eq!(a[0].regime, b[0].regime);
    assert_eq!(a[0].flags, b[0].flags);
}
