use crate::model::axes::{Axes, AxisDrivers, clip01};
use crate::model::drivers::ScoreDrivers;
use crate::model::scores::CompositeScores;
use crate::model::thresholds::{NuclearScoringMode, ThresholdProfile};

#[derive(Debug)]
pub struct Stage5Output {
    pub scores: CompositeScores,
    pub drivers: ScoreDrivers,
}

#[derive(Debug, Clone)]
pub struct Stage5Inputs<'a> {
    pub axes: &'a Axes,
    pub drivers: &'a [AxisDrivers],
    pub thresholds: &'a ThresholdProfile,
    pub n_genes_mappable: Option<u32>,
    pub key_panel_coverage_median: Option<&'a [f32]>,
    pub ambient_rna_risk: Option<&'a [bool]>,
    pub key_panels_missing: Option<&'a [bool]>,
    pub panel_nonzero_fraction: Option<&'a [f32]>,
    pub axis_p90: Option<[f32; 3]>,
    pub scoring_mode: NuclearScoringMode,
    pub include_ddr: bool,
}

pub fn run_stage5(inputs: &Stage5Inputs<'_>) -> Stage5Output {
    let n_cells = inputs.axes.tbi.len();
    let mut scores = CompositeScores {
        nps: vec![0.0; n_cells],
        ci: vec![0.0; n_cells],
        rls: vec![0.0; n_cells],
        confidence: vec![0.0; n_cells],
        confidence_breakdown: vec![[0.0, 0.0, 0.0, 0.0]; n_cells],
    };

    let mut drivers_out = ScoreDrivers {
        nps: vec![Vec::new(); n_cells],
        ci: vec![Vec::new(); n_cells],
        rls: vec![Vec::new(); n_cells],
    };

    for cell in 0..n_cells {
        let tbi = inputs.axes.tbi[cell];
        let rci = inputs.axes.rci[cell];
        let pds = inputs.axes.pds[cell];
        let trs = inputs.axes.trs[cell];
        let nsai = inputs.axes.nsai[cell];

        let nps = clip01(0.45 * tbi + 0.35 * rci - 0.20 * pds - 0.20 * trs);
        let mut ci = clip01(0.55 * trs + 0.45 * pds - 0.15 * tbi);
        let (confidence, breakdown) = match inputs.scoring_mode {
            NuclearScoringMode::StrictBulk => compute_confidence_legacy(inputs, cell),
            NuclearScoringMode::ImmuneAware => compute_confidence(inputs, cell),
        };
        let mut rls = match inputs.scoring_mode {
            NuclearScoringMode::StrictBulk => compute_rls_legacy(inputs, cell, confidence),
            NuclearScoringMode::ImmuneAware => compute_rls(inputs, cell, confidence),
        };

        if inputs.include_ddr {
            rls = clip01(rls - 0.25 * inputs.axes.rss[cell] - 0.20 * inputs.axes.trci[cell]);
            ci = clip01(ci + 0.15 * inputs.axes.cci[cell]);
        }

        scores.nps[cell] = nps;
        scores.ci[cell] = ci;
        scores.rls[cell] = rls;
        scores.confidence[cell] = confidence;
        scores.confidence_breakdown[cell] = breakdown;

        drivers_out.nps[cell] = top_k_drivers(vec![
            ("high_tbi", 0.45 * tbi),
            ("high_rci", 0.35 * rci),
            ("high_pds", -0.20 * pds),
            ("high_trs", -0.20 * trs),
        ]);

        let mut ci_drivers = vec![
            ("high_trs", 0.55 * trs),
            ("high_pds", 0.45 * pds),
            ("high_tbi", -0.15 * tbi),
        ];
        if inputs.include_ddr {
            ci_drivers.push(("high_cci", 0.15 * inputs.axes.cci[cell]));
        }
        drivers_out.ci[cell] = top_k_drivers(ci_drivers);

        let mut rls_drivers = vec![
            ("high_tbi", 0.45 * tbi),
            ("high_rci", 0.35 * rci),
            ("high_pds", -0.25 * pds),
            ("high_nsai", -0.15 * nsai),
        ];
        if inputs.include_ddr {
            rls_drivers.push(("high_rss", -0.25 * inputs.axes.rss[cell]));
            rls_drivers.push(("high_trci", -0.20 * inputs.axes.trci[cell]));
        }
        drivers_out.rls[cell] = top_k_drivers(rls_drivers);
    }

    Stage5Output {
        scores,
        drivers: drivers_out,
    }
}

fn compute_confidence(inputs: &Stage5Inputs<'_>, cell: usize) -> (f32, [f32; 4]) {
    let key_cov = inputs
        .key_panel_coverage_median
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(0.0);
    let missing_key = inputs
        .key_panels_missing
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(false);

    let n_genes_mappable = inputs
        .n_genes_mappable
        .unwrap_or(inputs.thresholds.min_expr_genes.max(1));
    let expressed = inputs
        .drivers
        .get(cell)
        .map(|d| d.expressed_genes)
        .unwrap_or(0);
    let expr_frac = clip01(expressed as f32 / n_genes_mappable as f32);
    let panel_nonzero = inputs
        .panel_nonzero_fraction
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(expr_frac);

    let _ambient_risk = inputs
        .ambient_rna_risk
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(true);

    let panel_coverage_score = if missing_key {
        0.0
    } else {
        clip01(key_cov / 0.6)
    };
    let expression_support_score = clip01(panel_nonzero.sqrt());
    let axis_structure_score = clip01(inputs.drivers[cell].axis_variance / 0.05);
    let consistency_score = consistency_score(inputs, cell);

    let conf = clip01(
        0.30 * panel_coverage_score
            + 0.25 * expression_support_score
            + 0.25 * axis_structure_score
            + 0.20 * consistency_score,
    );

    let conf = if inputs.key_panel_coverage_median.is_none()
        && inputs.panel_nonzero_fraction.is_none()
        && axis_structure_score == 0.0
    {
        0.0
    } else if !missing_key && axis_structure_score >= 0.2 {
        conf.max(0.2)
    } else {
        conf
    };

    (
        conf,
        [
            panel_coverage_score,
            expression_support_score,
            axis_structure_score,
            consistency_score,
        ],
    )
}

fn compute_confidence_legacy(inputs: &Stage5Inputs<'_>, cell: usize) -> (f32, [f32; 4]) {
    let key_cov = inputs
        .key_panel_coverage_median
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(0.0);
    let n_genes_mappable = inputs
        .n_genes_mappable
        .unwrap_or(inputs.thresholds.min_expr_genes.max(1));
    let expressed = inputs
        .drivers
        .get(cell)
        .map(|d| d.expressed_genes)
        .unwrap_or(0);
    let expr_frac = clip01(expressed as f32 / n_genes_mappable as f32);
    let ambient_risk = inputs
        .ambient_rna_risk
        .and_then(|v| v.get(cell).copied())
        .unwrap_or(true);
    let ambient = if ambient_risk { 1.0 } else { 0.0 };

    let a = 0.5 * key_cov;
    let b = 0.3 * expr_frac;
    let c = 0.2 * (1.0 - ambient);
    let conf = clip01(a * b * c);
    (conf, [a, b, c, 0.0])
}

fn consistency_score(inputs: &Stage5Inputs<'_>, cell: usize) -> f32 {
    let tbi = inputs.axes.tbi[cell];
    let trs = inputs.axes.trs[cell];
    let pds = inputs.axes.pds[cell];
    let rci = inputs.axes.rci[cell];

    let mut penalty = 0.0;
    let c1 = (trs + tbi - 1.2).max(0.0);
    let c2 = (pds + tbi - 1.2).max(0.0);
    let c3 = (trs + rci - 1.3).max(0.0);
    penalty += c1 + c2 + c3;

    clip01(1.0 - penalty)
}

fn compute_rls(inputs: &Stage5Inputs<'_>, cell: usize, confidence: f32) -> f32 {
    let tbi = inputs.axes.tbi[cell];
    let dfa = inputs.axes.dfa[cell];
    let iaa = inputs.axes.iaa[cell];
    let nsai = inputs.axes.nsai[cell];
    let trs = inputs.axes.trs[cell];
    let pds = inputs.axes.pds[cell];
    let axis_var = clip01(inputs.drivers[cell].axis_variance / 0.05);

    let rigid_commit = trs.max(pds);
    let mut rls = clip01(
        0.35 * tbi + 0.20 * dfa + 0.20 * iaa + 0.15 * nsai + 0.10 * axis_var - 0.30 * rigid_commit,
    );

    let allow_zero =
        tbi < 0.2 && dfa < 0.2 && iaa < 0.2 && nsai < 0.2 && axis_var < 0.05 && confidence >= 0.6;

    if !allow_zero {
        if let Some(p90) = inputs.axis_p90 {
            if p90[0] >= 0.8 || p90[1] >= 0.8 || p90[2] >= 0.8 {
                rls = rls.max(0.1);
            }
        }
    }

    rls
}

fn compute_rls_legacy(inputs: &Stage5Inputs<'_>, cell: usize, confidence: f32) -> f32 {
    let tbi = inputs.axes.tbi[cell];
    let rci = inputs.axes.rci[cell];
    let pds = inputs.axes.pds[cell];
    let nsai = inputs.axes.nsai[cell];
    let rls_base = clip01(0.45 * tbi + 0.35 * rci - 0.25 * pds - 0.15 * nsai);
    rls_base * confidence
}

fn top_k_drivers(items: Vec<(&'static str, f32)>) -> Vec<(String, f32)> {
    let mut v: Vec<(String, f32)> = items
        .into_iter()
        .map(|(name, value)| (name.to_string(), value))
        .collect();

    v.sort_by(|a, b| {
        let am = a.1.abs();
        let bm = b.1.abs();
        match bm.partial_cmp(&am).unwrap_or(std::cmp::Ordering::Equal) {
            std::cmp::Ordering::Equal => a.0.cmp(&b.0),
            other => other,
        }
    });

    if v.len() > 5 {
        v.truncate(5);
    }
    v
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage5_scores.rs"]
mod tests;
