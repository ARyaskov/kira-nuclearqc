use crate::model::axes::{Axes, AxisDrivers, AxisFlags, clip01};
use crate::model::ddr::{DdrMetrics, compute_ddr_metrics};
use crate::model::thresholds::{AxisActivationMode, ThresholdProfile};
use crate::panels::defs::PanelGroup;
use crate::panels::{PanelScores, PanelSet};
use crate::pipeline::stage2_normalize::ExprAccessor;
use crate::simd;

#[derive(Debug)]
pub struct Stage4Output {
    pub axes: Axes,
    pub drivers: Vec<AxisDrivers>,
    pub flags: Vec<AxisFlags>,
    pub ddr: DdrMetrics,
}

pub fn run_stage4(
    accessor: &dyn ExprAccessor,
    panel_set: &PanelSet,
    panel_scores: &PanelScores,
    thresholds: &ThresholdProfile,
) -> Stage4Output {
    let n_cells = accessor.n_cells();
    let n_panels = panel_set.panels.len();

    let mut program_panels = Vec::new();
    let mut tf_panels = Vec::new();
    let mut chromatin_panels = Vec::new();
    let mut stress_panels = Vec::new();
    let mut dev_panels = Vec::new();

    for (idx, panel) in panel_set.panels.iter().enumerate() {
        match panel.group {
            PanelGroup::Program => program_panels.push(idx),
            PanelGroup::Tf => tf_panels.push(idx),
            PanelGroup::Chromatin => chromatin_panels.push(idx),
            PanelGroup::Stress => stress_panels.push(idx),
            PanelGroup::Developmental => dev_panels.push(idx),
            _ => {}
        }
    }

    let mut axes = Axes {
        tbi: vec![0.0; n_cells],
        rci: vec![0.0; n_cells],
        pds: vec![0.0; n_cells],
        trs: vec![0.0; n_cells],
        nsai: vec![0.0; n_cells],
        iaa: vec![0.0; n_cells],
        dfa: vec![0.0; n_cells],
        cea: vec![0.0; n_cells],
        rss: vec![0.0; n_cells],
        drbi: vec![0.0; n_cells],
        cci: vec![0.0; n_cells],
        trci: vec![0.0; n_cells],
    };
    let mut drivers = vec![AxisDrivers::default(); n_cells];
    let mut flags = vec![AxisFlags::default(); n_cells];

    let mut value_buf: Vec<f32> = Vec::new();
    let mut program_buf: Vec<f32> = Vec::with_capacity(program_panels.len());
    let mut tf_buf: Vec<f32> = Vec::with_capacity(tf_panels.len() + chromatin_panels.len());
    let mut iaa_raw = vec![0.0f32; n_cells];
    let mut dfa_raw = vec![0.0f32; n_cells];
    let mut cea_raw = vec![0.0f32; n_cells];

    let replication_stress_panel = find_panel(panel_set, "replication_stress_genes");
    let checkpoint_activation_panel = find_panel(panel_set, "checkpoint_activation");
    let replication_fork_stability_panel = find_panel(panel_set, "replication_fork_stability");
    let dna_repair_hr_panel = find_panel(panel_set, "dna_repair_hr");
    let dna_repair_nhej_panel = find_panel(panel_set, "dna_repair_nhej");
    let chromatin_compaction_panel = find_panel(panel_set, "chromatin_compaction");
    let chromatin_open_state_panel = find_panel(panel_set, "chromatin_open_state");

    let iaa_panel = find_panel(panel_set, "immune_activation");
    let dfa_panel = find_panel(panel_set, "differentiation_flux");
    let cea_panel = find_panel(panel_set, "clonal_engagement");

    let mut replication_stress_raw = vec![0.0f32; n_cells];
    let mut checkpoint_activation_raw = vec![0.0f32; n_cells];
    let mut replication_fork_stability_raw = vec![0.0f32; n_cells];
    let mut hr_raw = vec![0.0f32; n_cells];
    let mut nhej_raw = vec![0.0f32; n_cells];
    let mut chromatin_compaction_raw = vec![0.0f32; n_cells];
    let mut chromatin_open_raw = vec![0.0f32; n_cells];

    for cell in 0..n_cells {
        if let Some(p) = iaa_panel {
            iaa_raw[cell] = panel_scores.panel_sum[cell][p];
        }
        if let Some(p) = dfa_panel {
            dfa_raw[cell] = panel_scores.panel_sum[cell][p];
        }
        if let Some(p) = cea_panel {
            cea_raw[cell] = panel_scores.panel_sum[cell][p];
        }
        if let Some(p) = replication_stress_panel {
            replication_stress_raw[cell] = panel_scores.panel_sum[cell][p];
        }
        if let Some(p) = checkpoint_activation_panel {
            checkpoint_activation_raw[cell] = panel_scores.panel_sum[cell][p];
        }
        if let Some(p) = replication_fork_stability_panel {
            replication_fork_stability_raw[cell] = panel_scores.panel_sum[cell][p];
        }
        if let Some(p) = dna_repair_hr_panel {
            hr_raw[cell] = panel_scores.panel_sum[cell][p];
        }
        if let Some(p) = dna_repair_nhej_panel {
            nhej_raw[cell] = panel_scores.panel_sum[cell][p];
        }
        if let Some(p) = chromatin_compaction_panel {
            chromatin_compaction_raw[cell] = panel_scores.panel_sum[cell][p];
        }
        if let Some(p) = chromatin_open_state_panel {
            chromatin_open_raw[cell] = panel_scores.panel_sum[cell][p];
        }
    }

    let iaa_rel = compute_relative_scores(&iaa_raw, thresholds);
    let dfa_rel = compute_relative_scores(&dfa_raw, thresholds);
    let cea_rel = compute_relative_scores(&cea_raw, thresholds);
    let replication_stress_norm = compute_relative_scores(&replication_stress_raw, thresholds);
    let checkpoint_activation_norm =
        compute_relative_scores(&checkpoint_activation_raw, thresholds);
    let replication_fork_stability_norm =
        compute_relative_scores(&replication_fork_stability_raw, thresholds);
    let hr_norm = compute_relative_scores(&hr_raw, thresholds);
    let nhej_norm = compute_relative_scores(&nhej_raw, thresholds);
    let chromatin_compaction_norm = compute_relative_scores(&chromatin_compaction_raw, thresholds);
    let chromatin_open_norm = compute_relative_scores(&chromatin_open_raw, thresholds);

    for cell in 0..n_cells {
        value_buf.clear();
        let nnz = accessor.nnz(cell) as usize;
        if value_buf.capacity() < nnz {
            value_buf.reserve(nnz - value_buf.capacity());
        }

        let mut expressed_genes = 0u32;
        accessor.for_cell(cell, &mut |_gene_id, value| {
            if value > 0.0 {
                value_buf.push(value);
            }
            if value > thresholds.expr_min {
                expressed_genes += 1;
            }
        });

        let n_genes_mappable = accessor.n_genes() as f32;
        let frac = if n_genes_mappable > 0.0 {
            expressed_genes as f32 / n_genes_mappable
        } else {
            0.0
        };
        let frac_norm = rescale01(
            frac,
            thresholds.frac_rescale_min,
            thresholds.frac_rescale_max,
        );

        let (gene_entropy, gene_entropy_norm) = entropy_norm_from_values(&value_buf);

        program_buf.clear();
        for &idx in &program_panels {
            program_buf.push(panel_scores.panel_sum[cell][idx]);
        }
        let (panel_entropy_norm, panel_entropy) = panel_entropy_program(&program_buf);

        let tbi = thresholds.tbi_w1 * frac_norm
            + thresholds.tbi_w2 * gene_entropy_norm
            + thresholds.tbi_w3 * panel_entropy_norm;

        tf_buf.clear();
        for &idx in tf_panels.iter().chain(chromatin_panels.iter()) {
            tf_buf.push(panel_scores.panel_sum[cell][idx]);
        }
        let (rci, tf_entropy, low_tf) = rci_score(&tf_buf, thresholds.tf_min_sum);

        let (pds, max_share) = pds_score(&program_buf, thresholds.program_min_sum);

        let trs = clip01(
            thresholds.trs_a * (1.0 - tbi)
                + thresholds.trs_b * (1.0 - rci)
                + thresholds.trs_c * pds,
        );

        let (nsai, stress_ratio, dev_ratio) = nsai_score(
            cell,
            panel_scores,
            &stress_panels,
            &dev_panels,
            &program_panels,
            thresholds.program_min_sum,
            thresholds.stress_boost,
        );

        let iaa = activate_axis(iaa_raw[cell], iaa_rel[cell], thresholds);
        let dfa = activate_axis(dfa_raw[cell], dfa_rel[cell], thresholds);
        let cea = activate_axis(cea_raw[cell], cea_rel[cell], thresholds);

        axes.tbi[cell] = clip01(tbi);
        axes.rci[cell] = clip01(rci);
        axes.pds[cell] = clip01(pds);
        axes.trs[cell] = trs;
        axes.nsai[cell] = nsai;
        axes.iaa[cell] = iaa;
        axes.dfa[cell] = dfa;
        axes.cea[cell] = cea;

        drivers[cell] = AxisDrivers {
            expressed_genes,
            gene_entropy,
            panel_entropy,
            max_program_share: max_share,
            tf_entropy,
            stress_ratio,
            dev_ratio,
            iaa_raw: iaa_raw[cell],
            dfa_raw: dfa_raw[cell],
            cea_raw: cea_raw[cell],
            axis_variance: 0.0,
        };
        flags[cell] = AxisFlags {
            low_tf_signal: low_tf,
        };

        let _ = n_panels; // silence unused variable warning if panels unused in build.
    }

    let ddr = compute_ddr_metrics(
        &replication_stress_norm,
        &checkpoint_activation_norm,
        &replication_fork_stability_norm,
        &hr_norm,
        &nhej_norm,
        &chromatin_compaction_norm,
        &chromatin_open_norm,
        &axes.tbi,
    );

    for cell in 0..n_cells {
        axes.rss[cell] = ddr.rss[cell];
        axes.drbi[cell] = ddr.drbi[cell];
        axes.cci[cell] = ddr.cci[cell];
        axes.trci[cell] = ddr.trci[cell];

        let axis_variance = axis_variance(
            axes.tbi[cell],
            axes.rci[cell],
            axes.pds[cell],
            axes.trs[cell],
            axes.nsai[cell],
            axes.iaa[cell],
            axes.dfa[cell],
            axes.cea[cell],
            axes.rss[cell],
            axes.drbi[cell],
            axes.cci[cell],
            axes.trci[cell],
        );
        drivers[cell].axis_variance = axis_variance;
    }

    Stage4Output {
        axes,
        drivers,
        flags,
        ddr,
    }
}

fn find_panel(panel_set: &PanelSet, id: &str) -> Option<usize> {
    panel_set.panels.iter().position(|p| p.id == id)
}

fn compute_relative_scores(values: &[f32], thresholds: &ThresholdProfile) -> Vec<f32> {
    if values.is_empty() {
        return Vec::new();
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    let p70 = sorted[((n - 1) as f32 * thresholds.rel_p70).ceil() as usize];
    let p85 = sorted[((n - 1) as f32 * thresholds.rel_p85).ceil() as usize];
    let mut out = Vec::with_capacity(values.len());
    for &v in values {
        if p85 <= p70 {
            out.push(0.0);
        } else {
            let x = (v - p70) / (p85 - p70);
            out.push(clip01(x));
        }
    }
    out
}

fn activate_axis(raw: f32, rel: f32, thresholds: &ThresholdProfile) -> f32 {
    match thresholds.activation_mode {
        AxisActivationMode::Absolute => clip01(raw),
        AxisActivationMode::Relative => rel,
        AxisActivationMode::Hybrid => clip01(0.5 * clip01(raw) + 0.5 * rel),
    }
}

fn axis_variance(
    tbi: f32,
    rci: f32,
    pds: f32,
    trs: f32,
    nsai: f32,
    iaa: f32,
    dfa: f32,
    cea: f32,
    rss: f32,
    drbi: f32,
    cci: f32,
    trci: f32,
) -> f32 {
    let vals = [
        tbi, rci, pds, trs, nsai, iaa, dfa, cea, rss, drbi, cci, trci,
    ];
    let mut mean = 0f64;
    for v in vals {
        mean += v as f64;
    }
    mean /= vals.len() as f64;
    let mut var = 0f64;
    for v in vals {
        let d = v as f64 - mean;
        var += d * d;
    }
    (var / vals.len() as f64) as f32
}
fn rescale01(x: f32, min: f32, max: f32) -> f32 {
    if max <= min {
        return 0.0;
    }
    let v = (x - min) / (max - min);
    clip01(v)
}

fn entropy_norm_from_values(values: &[f32]) -> (f32, f32) {
    if values.is_empty() {
        return (0.0, 0.0);
    }
    let h = simd::entropy_f32(values) as f64;
    let n = values.len();
    let h_norm = if n >= 2 {
        let denom = (n as f64).ln();
        if denom > 0.0 { h / denom } else { 0.0 }
    } else {
        0.0
    };
    (h as f32, h_norm as f32)
}

fn panel_entropy_program(values: &[f32]) -> (f32, f32) {
    if values.is_empty() {
        return (0.0, 0.0);
    }
    let sum = simd::sum_f32_f64(values);
    let mut nonzero = 0usize;
    for &v in values {
        if v > 0.0 {
            nonzero += 1;
        }
    }
    if sum <= 0.0 || nonzero < 2 {
        return (0.0, 0.0);
    }
    let mut h = 0f64;
    for &v in values {
        if v > 0.0 {
            let p = (v as f64) / sum;
            h -= p * p.ln();
        }
    }
    let h_norm = h / (nonzero as f64).ln();
    (h_norm as f32, h as f32)
}

fn rci_score(values: &[f32], tf_min_sum: f32) -> (f32, f32, bool) {
    let sum = simd::sum_f32_f64(values);
    let max = simd::max_f32(values) as f64;
    let mut nonzero = 0usize;
    for &v in values {
        if v > 0.0 {
            nonzero += 1;
        }
    }

    if sum < tf_min_sum as f64 {
        return (0.0, 0.0, true);
    }

    let mut h = 0f64;
    for &v in values {
        if v > 0.0 {
            let p = (v as f64) / sum;
            h -= p * p.ln();
        }
    }

    let entropy_norm = if nonzero >= 2 {
        let denom = (nonzero as f64).ln();
        if denom > 0.0 { h / denom } else { 0.0 }
    } else {
        0.0
    };

    let anti_dom = if sum > 0.0 { 1.0 - (max / sum) } else { 0.0 };
    let rci = 0.5 * entropy_norm + 0.5 * anti_dom;
    (rci as f32, entropy_norm as f32, false)
}

fn pds_score(values: &[f32], program_min_sum: f32) -> (f32, f32) {
    let sum = simd::sum_f32_f64(values);
    let mut top1 = 0f64;
    let mut top2 = 0f64;
    let mut top3 = 0f64;

    for &v in values {
        if v > 0.0 {
            let v = v as f64;
            if v >= top1 {
                top3 = top2;
                top2 = top1;
                top1 = v;
            } else if v >= top2 {
                top3 = top2;
                top2 = v;
            } else if v > top3 {
                top3 = v;
            }
        }
    }

    if sum < program_min_sum as f64 || sum == 0.0 {
        return (0.0, 0.0);
    }

    let max_share = top1 / sum;
    let top3_share = (top1 + top2 + top3) / sum;
    let pds = 0.7 * max_share + 0.3 * top3_share;
    (pds as f32, max_share as f32)
}

fn nsai_score(
    cell: usize,
    panel_scores: &PanelScores,
    stress_panels: &[usize],
    dev_panels: &[usize],
    program_panels: &[usize],
    program_min_sum: f32,
    stress_boost: f32,
) -> (f32, f32, f32) {
    let mut program_sum = 0f64;
    for &idx in program_panels {
        let v = panel_scores.panel_sum[cell][idx] as f64;
        if v > 0.0 {
            program_sum += v;
        }
    }
    if program_sum < program_min_sum as f64 || program_sum == 0.0 {
        return (0.0, 0.0, 0.0);
    }

    let mut stress_sum = 0f64;
    for &idx in stress_panels {
        let v = panel_scores.panel_sum[cell][idx] as f64;
        if v > 0.0 {
            stress_sum += v;
        }
    }

    let mut dev_sum = 0f64;
    for &idx in dev_panels {
        let v = panel_scores.panel_sum[cell][idx] as f64;
        if v > 0.0 {
            dev_sum += v;
        }
    }

    let stress_ratio = (stress_sum / program_sum) as f32;
    let dev_ratio = (dev_sum / program_sum) as f32;
    let nsai = clip01(stress_ratio - dev_ratio + stress_boost);
    (nsai, stress_ratio, dev_ratio)
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage4_axes.rs"]
mod tests;
