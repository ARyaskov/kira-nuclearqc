use std::fmt::Write;

use crate::report::{SummaryData, format_f32_6};

pub fn render_summary_json(data: &SummaryData) -> String {
    let mut out = String::new();
    out.push('{');

    // Aggregator contract keys.
    push_kv_str(&mut out, "tool", "kira-nuclearqc");
    out.push(',');
    out.push_str("\"input\":{");
    push_kv_str(&mut out, "mode", &data.run_mode);
    out.push(',');
    push_kv_str(&mut out, "resolution", &data.resolution);
    out.push(',');
    push_kv_num(&mut out, "n_cells", data.n_cells as f64);
    out.push(',');
    push_kv_num(&mut out, "n_genes_raw", data.n_genes_raw as f64);
    out.push(',');
    push_kv_num(&mut out, "n_genes_mappable", data.n_genes_mappable as f64);
    out.push(',');
    push_kv_str(&mut out, "species", &data.species);
    out.push(',');
    push_kv_str(&mut out, "scoring_mode", &data.scoring_mode);
    out.push(',');
    out.push_str("\"normalization\":{");
    push_kv_bool(&mut out, "normalize", data.normalize);
    out.push(',');
    push_kv_num(&mut out, "scale", data.scale as f64);
    out.push(',');
    push_kv_bool(&mut out, "log1p", data.log1p);
    out.push(',');
    push_kv_str(&mut out, "axis_activation_mode", &data.axis_activation_mode);
    if let Some(breakdowns) = &data.confidence_breakdown {
        out.push(',');
        out.push_str("\"confidence_breakdown_median\":{");
        push_kv_num(&mut out, "panel_coverage", breakdowns[0] as f64);
        out.push(',');
        push_kv_num(&mut out, "expr_fraction", breakdowns[1] as f64);
        out.push(',');
        push_kv_num(&mut out, "ambient_inverse", breakdowns[2] as f64);
        out.push(',');
        push_kv_num(&mut out, "consistency", breakdowns[3] as f64);
        out.push('}');
    }
    out.push_str("}");
    out.push_str("},");

    out.push_str("\"regimes\":{");
    for (i, r) in data.regimes.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        push_kv_num(&mut out, r.name, r.fraction as f64);
    }
    out.push_str("},");

    out.push_str("\"composites\":{");
    push_kv_num(
        &mut out,
        "nps_median",
        stat_median(&data.composites, "c1_nps") as f64,
    );
    out.push(',');
    push_kv_num(
        &mut out,
        "ci_median",
        stat_median(&data.composites, "c2_ci") as f64,
    );
    out.push(',');
    push_kv_num(
        &mut out,
        "rls_median",
        stat_median(&data.composites, "c3_rls") as f64,
    );
    out.push_str("},");

    out.push_str("\"tails\":{");
    push_kv_num(&mut out, "trs_p90", stat_p90(&data.axes, "a4_trs") as f64);
    out.push(',');
    push_kv_num(&mut out, "trs_ge_0_75", data.trs_ge_0_75 as f64);
    out.push(',');
    push_kv_num(&mut out, "nps_ge_0_60", data.nps_ge_0_60 as f64);
    out.push(',');
    push_kv_num(&mut out, "rls_le_0_35", data.rls_le_0_35 as f64);
    out.push_str("},");

    out.push_str("\"qc\":{");
    push_kv_num(
        &mut out,
        "low_confidence_fraction",
        data.low_confidence_fraction as f64,
    );
    out.push(',');
    push_kv_num(&mut out, "confidence_median", data.confidence_median as f64);
    out.push(',');
    push_kv_num(&mut out, "confidence_p10", data.confidence_p10 as f64);
    out.push(',');
    push_kv_num(
        &mut out,
        "low_expr_genes_fraction",
        data.low_expr_fraction as f64,
    );
    out.push_str("},");

    // Existing extended metadata and distributions.
    out.push_str("\"tool_meta\":{");
    push_kv_str(&mut out, "name", &data.tool_name);
    out.push(',');
    push_kv_str(&mut out, "version", &data.tool_version);
    out.push(',');
    out.push_str("\"git_hash\":");
    match &data.git_hash {
        Some(h) => push_str_val(&mut out, h),
        None => out.push_str("null"),
    }
    out.push(',');
    push_kv_str(&mut out, "simd_backend", &data.simd_backend);
    out.push_str("},");

    out.push_str("\"distributions\":{");
    out.push_str("\"axes\":[");
    for (i, s) in data.axes.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push('{');
        push_kv_str(&mut out, "name", s.name);
        out.push(',');
        push_kv_num(&mut out, "median", s.median as f64);
        out.push(',');
        push_kv_num(&mut out, "p90", s.p90 as f64);
        out.push(',');
        push_kv_num(&mut out, "p99", s.p99 as f64);
        out.push('}');
    }
    out.push_str("],");
    out.push_str("\"composites\":[");
    for (i, s) in data.composites.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push('{');
        push_kv_str(&mut out, "name", s.name);
        out.push(',');
        push_kv_num(&mut out, "median", s.median as f64);
        out.push(',');
        push_kv_num(&mut out, "p90", s.p90 as f64);
        out.push(',');
        push_kv_num(&mut out, "p99", s.p99 as f64);
        out.push('}');
    }
    out.push_str("]},");

    out.push_str("\"regime_stats\":[");
    for (i, r) in data.regimes.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push('{');
        push_kv_str(&mut out, "name", r.name);
        out.push(',');
        push_kv_num(&mut out, "count", r.count as f64);
        out.push(',');
        push_kv_num(&mut out, "fraction", r.fraction as f64);
        out.push('}');
    }
    out.push_str("],");

    out.push_str("\"regime_counts\":{");
    for (i, r) in data.regimes.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        push_kv_num(&mut out, r.name, r.count as f64);
    }
    out.push_str("},");

    out.push_str("\"panels\":{");
    out.push_str("\"missing_genes_by_panel\":{");
    for (i, (panel, genes)) in data.missing_genes_by_panel.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        push_str_key(&mut out, panel);
        out.push(':');
        out.push('[');
        for (j, g) in genes.iter().enumerate() {
            if j > 0 {
                out.push(',');
            }
            push_str_val(&mut out, g);
        }
        out.push(']');
    }
    out.push_str("},");
    out.push_str("\"rls_contributors_top\":[");
    for (i, name) in data.rls_contributors_top.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        push_str_val(&mut out, name);
    }
    out.push_str("]}");

    out.push('}');
    out
}

fn stat_median(stats: &[crate::report::NamedStats], name: &str) -> f32 {
    for s in stats {
        if s.name == name {
            return s.median;
        }
    }
    0.0
}

fn stat_p90(stats: &[crate::report::NamedStats], name: &str) -> f32 {
    for s in stats {
        if s.name == name {
            return s.p90;
        }
    }
    0.0
}

fn push_kv_str(out: &mut String, key: &str, value: &str) {
    push_str_key(out, key);
    out.push(':');
    push_str_val(out, value);
}

fn push_kv_num(out: &mut String, key: &str, value: f64) {
    push_str_key(out, key);
    out.push(':');
    let _ = write!(out, "{}", format_f32_6(value as f32));
}

fn push_kv_bool(out: &mut String, key: &str, value: bool) {
    push_str_key(out, key);
    out.push(':');
    out.push_str(if value { "true" } else { "false" });
}

fn push_str_key(out: &mut String, key: &str) {
    out.push('"');
    out.push_str(&escape_json(key));
    out.push('"');
}

fn push_str_val(out: &mut String, value: &str) {
    out.push('"');
    out.push_str(&escape_json(value));
    out.push('"');
}

fn escape_json(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 4);
    for ch in s.chars() {
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
