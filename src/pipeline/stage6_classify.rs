use crate::model::axes::AxisDrivers;
use crate::model::flags::{Flag, flag_order};
use crate::model::regimes::NuclearRegime;
use crate::model::scores::CompositeScores;
use crate::model::thresholds::{AxisActivationMode, NuclearScoringMode, ThresholdProfile};

#[derive(Debug, Clone)]
pub struct Classification {
    pub regime: NuclearRegime,
    pub flags: Vec<Flag>,
}

#[derive(Debug, Clone)]
pub struct Stage6Inputs<'a> {
    pub tbi: &'a [f32],
    pub rci: &'a [f32],
    pub pds: &'a [f32],
    pub trs: &'a [f32],
    pub nsai: &'a [f32],
    pub iaa: &'a [f32],
    pub dfa: &'a [f32],
    pub cea: &'a [f32],
    pub scores: &'a CompositeScores,
    pub drivers: &'a [AxisDrivers],
    pub thresholds: &'a ThresholdProfile,
    pub scoring_mode: NuclearScoringMode,
    pub key_panel_coverage_median: Option<&'a [f32]>,
    pub key_panels_missing: Option<&'a [bool]>,
    pub sum_tf_panels: Option<&'a [f32]>,
    pub ambient_rna_risk: Option<&'a [bool]>,
    pub proliferation_program_share: Option<&'a [f32]>,
    pub program_sum: Option<&'a [f32]>,
}

pub fn run_stage6(inputs: &Stage6Inputs<'_>) -> Vec<Classification> {
    let n_cells = inputs.tbi.len();
    let mut out = Vec::with_capacity(n_cells);

    for cell in 0..n_cells {
        let regime = classify_cell(inputs, cell);
        let flags = collect_flags(inputs, cell);
        out.push(Classification { regime, flags });
    }

    out
}

fn classify_cell(inputs: &Stage6Inputs<'_>, cell: usize) -> NuclearRegime {
    let expressed_genes = inputs.drivers[cell].expressed_genes;
    let gene_entropy = inputs.drivers[cell].gene_entropy;
    let program_sum = inputs
        .program_sum
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(0.0);

    let tbi = inputs.tbi[cell];
    let rci = inputs.rci[cell];
    let pds = inputs.pds[cell];
    let trs = inputs.trs[cell];
    let nsai = inputs.nsai[cell];

    if expressed_genes < inputs.thresholds.min_expr_genes
        || (tbi < 0.15 && gene_entropy < 0.10 && program_sum < inputs.thresholds.program_min_sum)
    {
        return NuclearRegime::TranscriptionallyCollapsed;
    }

    if trs >= 0.75 && nsai >= 0.55 && rci <= 0.35 {
        return NuclearRegime::RigidDegenerative;
    }

    if trs >= 0.70 && pds >= 0.60 && tbi <= 0.45 && nsai < 0.55 {
        return NuclearRegime::CommittedState;
    }

    if nsai >= 0.65 && rci >= 0.35 && (tbi >= 0.35 || pds <= 0.60) {
        return NuclearRegime::StressAdaptive;
    }

    if inputs.scores.nps[cell] >= 0.60 && trs <= 0.45 && pds <= 0.50 {
        return NuclearRegime::PlasticAdaptive;
    }

    if inputs.scoring_mode == NuclearScoringMode::ImmuneAware {
        if (inputs.scores.nps[cell] >= 0.45 || inputs.iaa[cell] >= 0.35 || inputs.dfa[cell] >= 0.35)
            && trs <= 0.55
            && pds <= 0.65
        {
            return NuclearRegime::TransientAdaptive;
        }
    }

    NuclearRegime::Unclassified
}

fn collect_flags(inputs: &Stage6Inputs<'_>, cell: usize) -> Vec<Flag> {
    let mut flags = Vec::new();

    let expressed_genes = inputs.drivers[cell].expressed_genes;
    let key_cov = inputs
        .key_panel_coverage_median
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(0.0);
    let missing_key = inputs
        .key_panels_missing
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(false);
    let sum_tf = inputs
        .sum_tf_panels
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(0.0);
    let ambient = inputs
        .ambient_rna_risk
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(false);
    let proliferation_share = inputs
        .proliferation_program_share
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(0.0);

    let _tbi = inputs.tbi[cell];
    let pds = inputs.pds[cell];
    let nsai = inputs.nsai[cell];
    let confidence = inputs.scores.confidence[cell];
    let axis_var = inputs.drivers[cell].axis_variance;

    if expressed_genes < inputs.thresholds.min_expr_genes {
        flags.push(Flag::LowExprGenes);
    }
    if key_cov < 0.4 {
        flags.push(Flag::LowPanelCoverage);
    }
    if missing_key {
        flags.push(Flag::MissingKeyPanels);
    }
    if pds > 0.75 {
        flags.push(Flag::HighProgramDominance);
    }
    if nsai > 0.75 {
        flags.push(Flag::HighStressBias);
    }
    if sum_tf < inputs.thresholds.tf_min_sum {
        flags.push(Flag::LowTfSignal);
    }
    if ambient {
        flags.push(Flag::AmbientRnaRisk);
    }
    if proliferation_share > 0.5 {
        flags.push(Flag::CellCycleConfounder);
    }
    if confidence < inputs.thresholds.confidence_low
        && (inputs.scoring_mode == NuclearScoringMode::StrictBulk || axis_var < 0.01)
    {
        flags.push(Flag::LowConfidence);
    }

    let model_limitation = inputs.thresholds.activation_mode != AxisActivationMode::Absolute
        || inputs.iaa[cell] > 0.0
        || inputs.dfa[cell] > 0.0
        || inputs.cea[cell] > 0.0;
    if model_limitation {
        flags.push(Flag::ModelLimitation);
    } else if confidence >= inputs.thresholds.confidence_low {
        flags.push(Flag::BiologicalSilence);
    }

    // stable ordering
    let mut ordered = Vec::new();
    for flag in flag_order() {
        if flags.contains(flag) {
            ordered.push(*flag);
        }
    }
    ordered
}

#[cfg(test)]
mod tests {
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
    }

    #[test]
    fn test_determinism() {
        let inputs = base_inputs();
        let a = run_stage6(&inputs.as_inputs());
        let b = run_stage6(&inputs.as_inputs());
        assert_eq!(a[0].regime, b[0].regime);
        assert_eq!(a[0].flags, b[0].flags);
    }
}
