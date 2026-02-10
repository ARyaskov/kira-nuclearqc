use std::cell::Cell;
use std::collections::BTreeMap;
use std::fmt::from_fn;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::model::drivers::ScoreDrivers;
use crate::model::flags::{Flag, flag_order};
use crate::model::regimes::NuclearRegime;
use crate::model::scores::CompositeScores;
use crate::panels::{PanelAudit, PanelScores, PanelSet};
use crate::report::json::render_summary_json;
use crate::report::text::render_report_text;
use crate::report::{
    NamedStats, RegimeStat, ReportContext, SummaryData, bool_fraction, format_f32_6, median, p10,
    p90, p99,
};

#[derive(Debug, Clone, Copy)]
pub enum ReportMode {
    Cell,
    Sample,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunMode {
    Standalone,
    Pipeline,
}

#[derive(Debug, Clone)]
pub struct PipelineContext {
    pub input_dir: String,
    pub input_source: String,
    pub shared_bin: Option<String>,
    pub run_mode: String,
}

#[derive(Debug, Clone)]
pub struct Stage7Input<'a> {
    pub barcodes: &'a [String],
    pub sample: Option<&'a [String]>,
    pub condition: Option<&'a [String]>,
    pub species_per_cell: Option<&'a [String]>,
    pub species_global: String,

    pub libsize: &'a [f32],
    pub nnz: &'a [u32],
    pub expressed_genes: &'a [u32],

    pub axes_tbi: &'a [f32],
    pub axes_rci: &'a [f32],
    pub axes_pds: &'a [f32],
    pub axes_trs: &'a [f32],
    pub axes_nsai: &'a [f32],
    pub axes_iaa: &'a [f32],
    pub axes_dfa: &'a [f32],
    pub axes_cea: &'a [f32],

    pub scores: &'a CompositeScores,
    pub drivers: &'a ScoreDrivers,
    pub activation_mode: String,
    pub scoring_mode: String,
    pub pipeline_context: Option<PipelineContext>,

    pub classifications: &'a [crate::pipeline::stage6_classify::Classification],

    pub panel_set: &'a PanelSet,
    pub panel_audits: &'a [PanelAudit],
    pub panel_scores: &'a PanelScores,

    pub tool_name: String,
    pub tool_version: String,
    pub git_hash: Option<String>,
    pub simd_backend: String,

    pub n_genes_raw: usize,
    pub n_genes_mappable: usize,

    pub normalize: bool,
    pub scale: f32,
    pub log1p: bool,
    pub confidence_breakdown: Option<&'a [[f32; 4]]>,
}

pub fn write_reports(
    input: &Stage7Input<'_>,
    out_dir: &Path,
    mode: ReportMode,
) -> std::io::Result<()> {
    fs::create_dir_all(out_dir)?;

    let nuclearqc_path = out_dir.join("nuclearqc.tsv");
    match mode {
        ReportMode::Cell => write_cell_tsv(input, &nuclearqc_path)?,
        ReportMode::Sample => write_sample_tsv(input, &nuclearqc_path)?,
    }

    let summary_path = out_dir.join("summary.json");
    let summary = build_summary(input, mode);
    let json = render_summary_json(&summary);
    write_text(&summary_path, &json)?;

    let report_path = out_dir.join("report.txt");
    let report_ctx = build_report_context(input, &summary);
    let report = render_report_text(&report_ctx);
    write_text(&report_path, &report)?;

    let panels_path = out_dir.join("panels_report.tsv");
    write_panels_report(input, &panels_path)?;

    if let Some(ctx) = &input.pipeline_context {
        if ctx.run_mode != "pipeline" {
            return Ok(());
        }
        let pipeline_path = out_dir.join("pipeline_step.json");
        let json = render_pipeline_step_json(&summary);
        write_text(&pipeline_path, &json)?;
    }

    Ok(())
}

fn write_cell_tsv(input: &Stage7Input<'_>, path: &Path) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    let header = [
        "barcode",
        "sample",
        "condition",
        "species",
        "libsize",
        "nnz",
        "expressed_genes",
        "confidence",
        "a1_tbi",
        "a2_rci",
        "a3_pds",
        "a4_trs",
        "a5_nsai",
        "a6_iaa",
        "a7_dfa",
        "a8_cea",
        "c1_nps",
        "c2_ci",
        "c3_rls",
        "regime",
        "flags",
        "drivers_nps",
        "drivers_ci",
        "drivers_rls",
        "top_program_panel",
        "top_program_share",
        "activation_mode",
    ]
    .join("\t");
    writeln!(w, "{}", header)?;

    let program_panels = program_panel_indices(input.panel_set);

    let n_cells = input.barcodes.len();
    let mut row_order = (0..n_cells).collect::<Vec<_>>();
    row_order.sort_by(|&a, &b| match input.barcodes[a].cmp(&input.barcodes[b]) {
        std::cmp::Ordering::Equal => a.cmp(&b),
        other => other,
    });

    for cell in row_order {
        let barcode = &input.barcodes[cell];
        let sample = input
            .sample
            .and_then(|v| v.get(cell))
            .cloned()
            .unwrap_or_default();
        let condition = input
            .condition
            .and_then(|v| v.get(cell))
            .cloned()
            .unwrap_or_default();
        let species = input
            .species_per_cell
            .and_then(|v| v.get(cell))
            .cloned()
            .unwrap_or_else(|| input.species_global.clone());

        let flags = format_flags(&input.classifications[cell].flags);
        let regime = regime_name(input.classifications[cell].regime);

        let drivers_nps = format_drivers(&input.drivers.nps[cell]);
        let drivers_ci = format_drivers(&input.drivers.ci[cell]);
        let drivers_rls = format_drivers(&input.drivers.rls[cell]);

        let (top_panel, top_share) =
            top_program_panel(cell, &program_panels, input.panel_set, input.panel_scores);

        let row = vec![
            barcode.to_string(),
            sample,
            condition,
            species,
            format_f32_6(input.libsize[cell]),
            input.nnz[cell].to_string(),
            input.expressed_genes[cell].to_string(),
            format_f32_6(input.scores.confidence[cell]),
            format_f32_6(input.axes_tbi[cell]),
            format_f32_6(input.axes_rci[cell]),
            format_f32_6(input.axes_pds[cell]),
            format_f32_6(input.axes_trs[cell]),
            format_f32_6(input.axes_nsai[cell]),
            format_f32_6(input.axes_iaa[cell]),
            format_f32_6(input.axes_dfa[cell]),
            format_f32_6(input.axes_cea[cell]),
            format_f32_6(input.scores.nps[cell]),
            format_f32_6(input.scores.ci[cell]),
            format_f32_6(input.scores.rls[cell]),
            regime.to_string(),
            flags,
            drivers_nps,
            drivers_ci,
            drivers_rls,
            top_panel,
            format_f32_6(top_share),
            input.activation_mode.clone(),
        ]
        .join("\t");
        writeln!(w, "{}", row)?;
    }

    Ok(())
}

fn write_sample_tsv(input: &Stage7Input<'_>, path: &Path) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);

    let regime_names = regime_names();

    let mut header = String::new();
    header.push_str("sample\tn_cells\t");
    for name in [
        "a1_tbi", "a2_rci", "a3_pds", "a4_trs", "a5_nsai", "a6_iaa", "a7_dfa", "a8_cea", "c1_nps",
        "c2_ci", "c3_rls",
    ] {
        header.push_str(name);
        header.push_str("_median\t");
        header.push_str(name);
        header.push_str("_p90\t");
        header.push_str(name);
        header.push_str("_p99\t");
    }
    header.push_str("regime_majority\t");
    for name in regime_names {
        header.push_str("regime_frac_");
        header.push_str(name);
        header.push('\t');
    }
    header.push_str("trs_ge_0_75\tnps_ge_0_60\trls_le_0_35");

    writeln!(w, "{}", header)?;

    let mut sample_map: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    let n_cells = input.barcodes.len();
    for cell in 0..n_cells {
        let key = input
            .sample
            .and_then(|v| v.get(cell))
            .cloned()
            .unwrap_or_default();
        sample_map.entry(key).or_default().push(cell);
    }

    for (sample, idxs) in sample_map {
        let n = idxs.len();
        let mut a1 = Vec::with_capacity(n);
        let mut a2 = Vec::with_capacity(n);
        let mut a3 = Vec::with_capacity(n);
        let mut a4 = Vec::with_capacity(n);
        let mut a5 = Vec::with_capacity(n);
        let mut a6 = Vec::with_capacity(n);
        let mut a7 = Vec::with_capacity(n);
        let mut a8 = Vec::with_capacity(n);
        let mut c1 = Vec::with_capacity(n);
        let mut c2 = Vec::with_capacity(n);
        let mut c3 = Vec::with_capacity(n);

        let mut trs_tail = 0usize;
        let mut nps_tail = 0usize;
        let mut rls_tail = 0usize;

        let mut regime_counts: BTreeMap<&'static str, usize> = BTreeMap::new();

        for &cell in &idxs {
            a1.push(input.axes_tbi[cell]);
            a2.push(input.axes_rci[cell]);
            a3.push(input.axes_pds[cell]);
            a4.push(input.axes_trs[cell]);
            a5.push(input.axes_nsai[cell]);
            a6.push(input.axes_iaa[cell]);
            a7.push(input.axes_dfa[cell]);
            a8.push(input.axes_cea[cell]);
            c1.push(input.scores.nps[cell]);
            c2.push(input.scores.ci[cell]);
            c3.push(input.scores.rls[cell]);

            if input.axes_trs[cell] >= 0.75 {
                trs_tail += 1;
            }
            if input.scores.nps[cell] >= 0.60 {
                nps_tail += 1;
            }
            if input.scores.rls[cell] <= 0.35 {
                rls_tail += 1;
            }

            let r = regime_name(input.classifications[cell].regime);
            *regime_counts.entry(r).or_insert(0) += 1;
        }

        let majority = majority_regime(&regime_counts, regime_names);

        let mut line = String::new();
        line.push_str(&sample);
        line.push('\t');
        line.push_str(&n.to_string());
        line.push('\t');

        for v in [
            stats(&a1),
            stats(&a2),
            stats(&a3),
            stats(&a4),
            stats(&a5),
            stats(&a6),
            stats(&a7),
            stats(&a8),
            stats(&c1),
            stats(&c2),
            stats(&c3),
        ] {
            line.push_str(&format_f32_6(v.0));
            line.push('\t');
            line.push_str(&format_f32_6(v.1));
            line.push('\t');
            line.push_str(&format_f32_6(v.2));
            line.push('\t');
        }

        line.push_str(majority);
        line.push('\t');
        for name in regime_names {
            let count = *regime_counts.get(name).unwrap_or(&0) as f32;
            let frac = if n > 0 { count / n as f32 } else { 0.0 };
            line.push_str(&format_f32_6(frac));
            line.push('\t');
        }

        line.push_str(&format_f32_6(trs_tail as f32 / n as f32));
        line.push('\t');
        line.push_str(&format_f32_6(nps_tail as f32 / n as f32));
        line.push('\t');
        line.push_str(&format_f32_6(rls_tail as f32 / n as f32));

        writeln!(w, "{}", line)?;
    }

    Ok(())
}

fn write_panels_report(input: &Stage7Input<'_>, path: &Path) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "panel_id\tpanel_name\tpanel_group\tpanel_size_defined\tpanel_size_mappable\tmissing_genes\tcoverage_median\tcoverage_p10\tsum_median\tsum_p90\tsum_p99"
    )?;

    let n_cells = input.barcodes.len();
    let n_panels = input.panel_set.panels.len();

    for panel_idx in 0..n_panels {
        let panel = &input.panel_set.panels[panel_idx];
        let audit = input
            .panel_audits
            .iter()
            .find(|a| a.panel_id == panel.id)
            .cloned();

        let mut coverage = Vec::with_capacity(n_cells);
        let mut sums = Vec::with_capacity(n_cells);
        for cell in 0..n_cells {
            coverage.push(input.panel_scores.panel_coverage[cell][panel_idx]);
            sums.push(input.panel_scores.panel_sum[cell][panel_idx]);
        }

        let missing = audit
            .as_ref()
            .map(|a| a.missing_genes.join(","))
            .unwrap_or_default();
        let size_defined = audit.as_ref().map(|a| a.panel_size_defined).unwrap_or(0);
        let size_mappable = audit.as_ref().map(|a| a.panel_size_mappable).unwrap_or(0);

        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            panel.id,
            panel.name,
            panel_group_name(panel.group),
            size_defined,
            size_mappable,
            missing,
            format_f32_6(median(&coverage)),
            format_f32_6(p10(&coverage)),
            format_f32_6(median(&sums)),
            format_f32_6(p90(&sums)),
            format_f32_6(p99(&sums)),
        )?;
    }

    Ok(())
}

fn build_summary(input: &Stage7Input<'_>, mode: ReportMode) -> SummaryData {
    let n_cells = input.barcodes.len();

    let confidence = input.scores.confidence.to_vec();
    let confidence_median = median(&confidence);
    let confidence_p10 = p10(&confidence);

    let low_conf = input
        .classifications
        .iter()
        .map(|c| c.flags.contains(&Flag::LowConfidence))
        .collect::<Vec<_>>();
    let low_expr = input
        .classifications
        .iter()
        .map(|c| c.flags.contains(&Flag::LowExprGenes))
        .collect::<Vec<_>>();

    let axes = vec![
        named_stats("a1_tbi", input.axes_tbi),
        named_stats("a2_rci", input.axes_rci),
        named_stats("a3_pds", input.axes_pds),
        named_stats("a4_trs", input.axes_trs),
        named_stats("a5_nsai", input.axes_nsai),
        named_stats("a6_iaa", input.axes_iaa),
        named_stats("a7_dfa", input.axes_dfa),
        named_stats("a8_cea", input.axes_cea),
    ];
    let composites = vec![
        named_stats("c1_nps", &input.scores.nps),
        named_stats("c2_ci", &input.scores.ci),
        named_stats("c3_rls", &input.scores.rls),
    ];

    let regimes = regime_stats(input.classifications, n_cells);

    let trs_ge_0_75 = fraction_threshold(input.axes_trs, |v| v >= 0.75);
    let nps_ge_0_60 = fraction_threshold(&input.scores.nps, |v| v >= 0.60);
    let rls_le_0_35 = fraction_threshold(&input.scores.rls, |v| v <= 0.35);

    let missing_genes_by_panel = input
        .panel_audits
        .iter()
        .map(|a| (a.panel_id.clone(), a.missing_genes.clone()))
        .collect::<Vec<_>>();
    let rls_contributors_top = top_rls_contributors(input);

    SummaryData {
        tool_name: input.tool_name.clone(),
        tool_version: input.tool_version.clone(),
        git_hash: input.git_hash.clone(),
        simd_backend: input.simd_backend.clone(),
        run_mode: input
            .pipeline_context
            .as_ref()
            .map(|ctx| ctx.run_mode.clone())
            .unwrap_or_else(|| "standalone".to_string()),
        resolution: match mode {
            ReportMode::Cell => "cell".to_string(),
            ReportMode::Sample => "sample".to_string(),
        },

        n_cells,
        n_genes_raw: input.n_genes_raw,
        n_genes_mappable: input.n_genes_mappable,
        species: input.species_global.clone(),

        normalize: input.normalize,
        scale: input.scale,
        log1p: input.log1p,
        axis_activation_mode: input.activation_mode.clone(),
        confidence_breakdown: input
            .confidence_breakdown
            .map(|v| confidence_breakdown_median(v)),
        scoring_mode: input.scoring_mode.clone(),

        confidence_median,
        confidence_p10,
        low_confidence_fraction: bool_fraction(&low_conf),
        low_expr_fraction: bool_fraction(&low_expr),

        axes,
        composites,
        regimes,

        trs_ge_0_75,
        nps_ge_0_60,
        rls_le_0_35,

        missing_genes_by_panel,
        rls_contributors_top,
    }
}

fn top_rls_contributors(input: &Stage7Input<'_>) -> Vec<String> {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for drivers in &input.drivers.rls {
        for (name, _val) in drivers {
            *counts.entry(name.clone()).or_insert(0) += 1;
        }
    }
    let mut items = counts.into_iter().collect::<Vec<_>>();
    items.sort_by(|a, b| match b.1.cmp(&a.1) {
        std::cmp::Ordering::Equal => a.0.cmp(&b.0),
        other => other,
    });
    items.into_iter().take(3).map(|(n, _)| n).collect()
}

fn confidence_breakdown_median(values: &[[f32; 4]]) -> [f32; 4] {
    if values.is_empty() {
        return [0.0, 0.0, 0.0, 0.0];
    }
    let a = values.iter().map(|v| v[0]).collect::<Vec<_>>();
    let b = values.iter().map(|v| v[1]).collect::<Vec<_>>();
    let c = values.iter().map(|v| v[2]).collect::<Vec<_>>();
    let d = values.iter().map(|v| v[3]).collect::<Vec<_>>();
    [median(&a), median(&b), median(&c), median(&d)]
}

fn build_report_context(input: &Stage7Input<'_>, summary: &SummaryData) -> ReportContext {
    let ambient = input
        .classifications
        .iter()
        .map(|c| c.flags.contains(&Flag::AmbientRnaRisk))
        .collect::<Vec<_>>();
    let cell_cycle = input
        .classifications
        .iter()
        .map(|c| c.flags.contains(&Flag::CellCycleConfounder))
        .collect::<Vec<_>>();

    ReportContext {
        n_cells: input.barcodes.len(),
        regimes: summary.regimes.clone(),
        nps_median: median(&input.scores.nps),
        ci_median: median(&input.scores.ci),
        nsai_median: median(input.axes_nsai),
        rls_median: median(&input.scores.rls),
        low_confidence_fraction: summary.low_confidence_fraction,
        low_expr_fraction: summary.low_expr_fraction,
        ambient_rna_fraction: bool_fraction(&ambient),
        cell_cycle_fraction: bool_fraction(&cell_cycle),
        immune_note: input.activation_mode != "Absolute",
        confidence_breakdown: summary.confidence_breakdown,
        rls_contributors_top: summary.rls_contributors_top.clone(),
        rls_tail_fraction: summary.rls_le_0_35,
        immune_tail_note: immune_tail_note(input),
        scoring_mode: summary.scoring_mode.clone(),
        axis_activation_mode: summary.axis_activation_mode.clone(),
        confidence_model: if summary.scoring_mode.contains("strict") {
            "legacy multiplicative".to_string()
        } else {
            "immune-calibrated additive".to_string()
        },
    }
}

fn render_pipeline_step_json(summary: &SummaryData) -> String {
    let mut out = String::new();
    out.push('{');
    push_kv_str(&mut out, "tool", "kira-nuclearqc");
    out.push(',');
    push_kv_str(&mut out, "mode", "pipeline");
    out.push(',');

    out.push_str("\"artifacts\":{");
    push_kv_str(&mut out, "summary", "summary.json");
    out.push(',');
    push_kv_str(&mut out, "primary_metrics", "nuclearqc.tsv");
    out.push_str("},");

    out.push_str("\"cell_metrics\":{");
    push_kv_str(&mut out, "file", "nuclearqc.tsv");
    out.push(',');
    push_kv_str(&mut out, "regime_column", "regime");
    out.push(',');
    push_kv_str(&mut out, "confidence_column", "confidence");
    out.push_str("},");

    let nps_median = stat_median(&summary.composites, "c1_nps");
    let ci_median = stat_median(&summary.composites, "c2_ci");
    let rls_median = stat_median(&summary.composites, "c3_rls");
    let trs_p90 = stat_p90(&summary.axes, "a4_trs");

    out.push_str("\"key_metrics\":{");
    push_kv_num(&mut out, "nps_median", nps_median as f64);
    out.push(',');
    push_kv_num(&mut out, "ci_median", ci_median as f64);
    out.push(',');
    push_kv_num(&mut out, "rls_median", rls_median as f64);
    out.push(',');
    push_kv_num(&mut out, "trs_p90", trs_p90 as f64);
    out.push(',');
    push_kv_num(
        &mut out,
        "low_confidence_fraction",
        summary.low_confidence_fraction as f64,
    );
    out.push_str("}");

    out.push('}');
    out
}

fn stat_median(stats: &[NamedStats], name: &str) -> f32 {
    for s in stats {
        if s.name == name {
            return s.median;
        }
    }
    0.0
}

fn stat_p90(stats: &[NamedStats], name: &str) -> f32 {
    for s in stats {
        if s.name == name {
            return s.p90;
        }
    }
    0.0
}

fn push_kv_str(out: &mut String, key: &str, value: &str) {
    out.push('"');
    out.push_str(key);
    out.push('"');
    out.push(':');
    out.push('"');
    out.push_str(&escape_json(value));
    out.push('"');
}

fn push_kv_num(out: &mut String, key: &str, value: f64) {
    out.push('"');
    out.push_str(key);
    out.push('"');
    out.push(':');
    let _ = std::fmt::Write::write_fmt(out, format_args!("{}", format_f32_6(value as f32)));
}

fn push_str_val(out: &mut String, value: &str) {
    out.push('"');
    out.push_str(&escape_json(value));
    out.push('"');
}

fn escape_json(input: &str) -> String {
    let mut out = String::new();
    for ch in input.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            _ => out.push(ch),
        }
    }
    out
}

fn immune_tail_note(input: &Stage7Input<'_>) -> bool {
    let p90_iaa = p90(input.axes_iaa);
    let p90_dfa = p90(input.axes_dfa);
    let p90_cea = p90(input.axes_cea);
    p90_iaa >= 0.8 || p90_dfa >= 0.8 || p90_cea >= 0.8
}

fn write_text(path: &Path, contents: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    w.write_all(contents.as_bytes())?;
    Ok(())
}

fn named_stats(name: &'static str, values: &[f32]) -> NamedStats {
    NamedStats {
        name,
        median: median(values),
        p90: p90(values),
        p99: p99(values),
    }
}

fn fraction_threshold(values: &[f32], predicate: impl Fn(f32) -> bool) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mut count = 0usize;
    for &v in values {
        if predicate(v) {
            count += 1;
        }
    }
    count as f32 / values.len() as f32
}

fn format_flags(flags: &[Flag]) -> String {
    let order = flag_order();
    let idx = Cell::new(0usize);
    let wrote_any = Cell::new(false);
    from_fn(|f| {
        while idx.get() < order.len() {
            let flag = order[idx.get()];
            idx.set(idx.get() + 1);
            if !flags.contains(&flag) {
                continue;
            }
            if wrote_any.get() {
                f.write_str(",")?;
            } else {
                wrote_any.set(true);
            }
            f.write_str(flag_name(flag))?;
        }
        Ok(())
    })
    .to_string()
}

fn format_drivers(drivers: &[(String, f32)]) -> String {
    let idx = Cell::new(0usize);
    from_fn(|f| {
        while idx.get() < drivers.len() {
            let (name, value) = &drivers[idx.get()];
            idx.set(idx.get() + 1);
            if idx.get() > 1 {
                f.write_str(",")?;
            }
            f.write_str(name)?;
            f.write_str(":")?;
            f.write_fmt(format_args!("{:.6}", *value))?;
        }
        Ok(())
    })
    .to_string()
}

fn top_program_panel(
    cell: usize,
    program_panels: &[usize],
    panel_set: &PanelSet,
    panel_scores: &PanelScores,
) -> (String, f32) {
    if program_panels.is_empty() {
        return (String::new(), 0.0);
    }
    let mut sum = 0f64;
    let mut max = -1f64;
    let mut max_idx = None;
    for &idx in program_panels {
        let v = panel_scores.panel_sum[cell][idx] as f64;
        sum += v;
        if v > max {
            max = v;
            max_idx = Some(idx);
        }
    }
    if sum <= 0.0 {
        return (String::new(), 0.0);
    }
    let idx = max_idx.unwrap();
    let share = (panel_scores.panel_sum[cell][idx] as f64 / sum) as f32;
    (panel_set.panels[idx].id.to_string(), share)
}

fn program_panel_indices(panel_set: &PanelSet) -> Vec<usize> {
    let mut out = Vec::new();
    for (idx, panel) in panel_set.panels.iter().enumerate() {
        if panel.group == crate::panels::defs::PanelGroup::Program {
            out.push(idx);
        }
    }
    out
}

fn stats(values: &[f32]) -> (f32, f32, f32) {
    (median(values), p90(values), p99(values))
}

fn majority_regime<'a>(counts: &BTreeMap<&'a str, usize>, order: &'a [&'a str]) -> &'a str {
    let mut best = ("Unclassified", 0usize);
    for name in order {
        let count = *counts.get(*name).unwrap_or(&0);
        if count > best.1 {
            best = (*name, count);
        }
    }
    best.0
}

fn regime_stats(
    classifications: &[crate::pipeline::stage6_classify::Classification],
    n_cells: usize,
) -> Vec<RegimeStat> {
    let names = regime_names();
    let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
    for c in classifications {
        let name = regime_name(c.regime);
        *counts.entry(name).or_insert(0) += 1;
    }
    let mut out = Vec::new();
    for name in names {
        let count = *counts.get(name).unwrap_or(&0);
        let fraction = if n_cells > 0 {
            count as f32 / n_cells as f32
        } else {
            0.0
        };
        out.push(RegimeStat {
            name,
            count,
            fraction,
        });
    }
    out
}

fn regime_names() -> &'static [&'static str] {
    &[
        "PlasticAdaptive",
        "StressAdaptive",
        "CommittedState",
        "RigidDegenerative",
        "TranscriptionallyCollapsed",
        "TransientAdaptive",
        "Unclassified",
    ]
}

fn regime_name(r: NuclearRegime) -> &'static str {
    match r {
        NuclearRegime::PlasticAdaptive => "PlasticAdaptive",
        NuclearRegime::StressAdaptive => "StressAdaptive",
        NuclearRegime::CommittedState => "CommittedState",
        NuclearRegime::RigidDegenerative => "RigidDegenerative",
        NuclearRegime::TranscriptionallyCollapsed => "TranscriptionallyCollapsed",
        NuclearRegime::TransientAdaptive => "TransientAdaptive",
        NuclearRegime::Unclassified => "Unclassified",
    }
}

fn flag_name(flag: Flag) -> &'static str {
    match flag {
        Flag::LowExprGenes => "LOW_EXPR_GENES",
        Flag::LowPanelCoverage => "LOW_PANEL_COVERAGE",
        Flag::MissingKeyPanels => "MISSING_KEY_PANELS",
        Flag::HighProgramDominance => "HIGH_PROGRAM_DOMINANCE",
        Flag::HighStressBias => "HIGH_STRESS_BIAS",
        Flag::LowTfSignal => "LOW_TF_SIGNAL",
        Flag::AmbientRnaRisk => "AMBIENT_RNA_RISK",
        Flag::CellCycleConfounder => "CELL_CYCLE_CONFOUNDER",
        Flag::LowConfidence => "LOW_CONFIDENCE",
        Flag::ModelLimitation => "MODEL_LIMITATION",
        Flag::BiologicalSilence => "BIOLOGICAL_SILENCE",
    }
}

fn panel_group_name(group: crate::panels::defs::PanelGroup) -> &'static str {
    use crate::panels::defs::PanelGroup;
    match group {
        PanelGroup::Housekeeping => "housekeeping",
        PanelGroup::Tf => "tf",
        PanelGroup::Chromatin => "chromatin",
        PanelGroup::Stress => "stress",
        PanelGroup::Developmental => "developmental",
        PanelGroup::Proliferation => "proliferation",
        PanelGroup::Program => "program",
        PanelGroup::Confounder => "confounder",
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/pipeline/stage7_report.rs"]
mod tests;
