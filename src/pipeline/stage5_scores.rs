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
        let ci = clip01(0.55 * trs + 0.45 * pds - 0.15 * tbi);
        let (confidence, breakdown) = match inputs.scoring_mode {
            NuclearScoringMode::StrictBulk => compute_confidence_legacy(inputs, cell),
            NuclearScoringMode::ImmuneAware => compute_confidence(inputs, cell),
        };
        let rls = match inputs.scoring_mode {
            NuclearScoringMode::StrictBulk => compute_rls_legacy(inputs, cell, confidence),
            NuclearScoringMode::ImmuneAware => compute_rls(inputs, cell, confidence),
        };

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

        drivers_out.ci[cell] = top_k_drivers(vec![
            ("high_trs", 0.55 * trs),
            ("high_pds", 0.45 * pds),
            ("high_tbi", -0.15 * tbi),
        ]);

        drivers_out.rls[cell] = top_k_drivers(vec![
            ("high_tbi", 0.45 * tbi),
            ("high_rci", 0.35 * rci),
            ("high_pds", -0.25 * pds),
            ("high_nsai", -0.15 * nsai),
        ]);
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
mod tests {
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
}
