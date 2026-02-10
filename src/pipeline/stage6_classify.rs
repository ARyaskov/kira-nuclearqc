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
#[path = "../../tests/src_inline/pipeline/stage6_classify.rs"]
mod tests;
