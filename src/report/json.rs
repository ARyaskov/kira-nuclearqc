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
    out.push_str("\"ddr_metrics\":{");
    push_kv_obj_stats(&mut out, "rss", &data.ddr_metrics);
    out.push(',');
    push_kv_obj_stats(&mut out, "drbi", &data.ddr_metrics);
    out.push(',');
    push_kv_obj_stats(&mut out, "cci", &data.ddr_metrics);
    out.push(',');
    push_kv_obj_stats(&mut out, "trci", &data.ddr_metrics);
    out.push_str("},");
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
    out.push(',');
    out.push_str("\"genome_stability\":{");
    push_kv_str(
        &mut out,
        "panel_version",
        data.genome_stability.panel_version,
    );
    out.push(',');
    out.push_str("\"thresholds\":{");
    push_kv_num(
        &mut out,
        "replication_stress_high_rss",
        data.genome_stability.thresholds.replication_stress_high_rss as f64,
    );
    out.push(',');
    push_kv_num(
        &mut out,
        "ddr_high",
        data.genome_stability.thresholds.ddr_high as f64,
    );
    out.push(',');
    push_kv_num(
        &mut out,
        "checkpoint_addicted_cds",
        data.genome_stability.thresholds.checkpoint_addicted_cds as f64,
    );
    out.push(',');
    push_kv_num(
        &mut out,
        "checkpoint_addicted_sas_max",
        data.genome_stability.thresholds.checkpoint_addicted_sas_max as f64,
    );
    out.push(',');
    push_kv_num(
        &mut out,
        "senescent_like_sas",
        data.genome_stability.thresholds.senescent_like_sas as f64,
    );
    out.push(',');
    push_kv_num(
        &mut out,
        "senescent_like_sphase_z_max",
        data.genome_stability.thresholds.senescent_like_sphase_z_max as f64,
    );
    out.push(',');
    push_kv_num(
        &mut out,
        "gir_rb_nhej_skew_max",
        data.genome_stability.thresholds.gir_rb_nhej_skew_max as f64,
    );
    out.push(',');
    push_kv_num(
        &mut out,
        "panel_min_genes",
        data.genome_stability.thresholds.panel_min_genes as f64,
    );
    out.push(',');
    push_kv_num(
        &mut out,
        "trimmed_mean_fraction",
        data.genome_stability.thresholds.trimmed_mean_fraction as f64,
    );
    out.push(',');
    push_kv_num(
        &mut out,
        "robust_z_eps",
        data.genome_stability.thresholds.robust_z_eps as f64,
    );
    out.push_str("},");
    out.push_str("\"panel_audits\":[");
    for (idx, audit) in data.genome_stability.panel_audits.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push('{');
        push_kv_str(&mut out, "panel_id", &audit.panel_id);
        out.push(',');
        push_kv_num(
            &mut out,
            "panel_size_defined",
            audit.panel_size_defined as f64,
        );
        out.push(',');
        push_kv_num(
            &mut out,
            "panel_size_mappable",
            audit.panel_size_mappable as f64,
        );
        out.push(',');
        out.push_str("\"missing_genes\":[");
        for (gene_idx, gene) in audit.missing_genes.iter().enumerate() {
            if gene_idx > 0 {
                out.push(',');
            }
            push_str_val(&mut out, gene);
        }
        out.push(']');
        out.push('}');
    }
    out.push_str("],");
    out.push_str("\"global_stats\":{");
    out.push_str("\"robust_norm\":[");
    for (idx, stat) in data
        .genome_stability
        .global_stats
        .robust_norm
        .iter()
        .enumerate()
    {
        if idx > 0 {
            out.push(',');
        }
        out.push('{');
        push_kv_str(&mut out, "name", stat.name);
        out.push(',');
        push_kv_num(&mut out, "median", stat.median as f64);
        out.push(',');
        push_kv_num(&mut out, "mad", stat.mad as f64);
        out.push('}');
    }
    out.push_str("],");
    out.push_str("\"distributions\":");
    push_genome_metric_distributions(&mut out, &data.genome_stability.global_stats.distributions);
    out.push(',');
    out.push_str("\"flag_fractions\":");
    push_genome_flag_fractions(&mut out, &data.genome_stability.global_stats.flag_fractions);
    out.push(',');
    out.push_str("\"missing_fraction\":");
    push_genome_flag_fractions(
        &mut out,
        &data.genome_stability.global_stats.missing_fraction,
    );
    out.push(',');
    out.push_str("\"top_clusters_by_rss\":");
    push_top_clusters(
        &mut out,
        &data.genome_stability.global_stats.top_clusters_by_rss,
    );
    out.push(',');
    out.push_str("\"top_clusters_by_checkpoint_addicted\":");
    push_top_clusters(
        &mut out,
        &data
            .genome_stability
            .global_stats
            .top_clusters_by_checkpoint_addicted,
    );
    out.push_str("},");
    out.push_str("\"cluster_stats\":[");
    for (idx, cluster) in data.genome_stability.cluster_stats.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push('{');
        push_kv_str(&mut out, "cluster_id", &cluster.cluster_id);
        out.push(',');
        push_kv_num(&mut out, "n_cells", cluster.n_cells as f64);
        out.push(',');
        out.push_str("\"metrics\":");
        push_genome_metric_distributions(&mut out, &cluster.metrics);
        out.push(',');
        out.push_str("\"flag_fractions\":");
        push_genome_flag_fractions(&mut out, &cluster.flag_fractions);
        out.push('}');
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

fn push_kv_obj_stats(out: &mut String, key: &str, stats: &[crate::report::NamedStats]) {
    out.push('"');
    out.push_str(key);
    out.push_str("\":{");
    push_kv_num(out, "median", stat_median(stats, key) as f64);
    out.push(',');
    push_kv_num(out, "p90", stat_p90(stats, key) as f64);
    out.push(',');
    push_kv_num(out, "p99", stat_p99(stats, key) as f64);
    out.push('}');
}

fn stat_p99(stats: &[crate::report::NamedStats], name: &str) -> f32 {
    for s in stats {
        if s.name == name {
            return s.p99;
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

fn push_genome_metric_distributions(
    out: &mut String,
    values: &[crate::metrics::genome_stability::aggregate::GenomeMetricDistribution],
) {
    out.push('[');
    for (idx, metric) in values.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push('{');
        push_kv_str(out, "name", metric.name);
        out.push(',');
        push_kv_num(out, "median", metric.median as f64);
        out.push(',');
        push_kv_num(out, "p10", metric.p10 as f64);
        out.push(',');
        push_kv_num(out, "p90", metric.p90 as f64);
        out.push('}');
    }
    out.push(']');
}

fn push_genome_flag_fractions(
    out: &mut String,
    values: &[crate::metrics::genome_stability::aggregate::GenomeFlagFraction],
) {
    out.push('[');
    for (idx, flag) in values.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push('{');
        push_kv_str(out, "name", flag.name);
        out.push(',');
        push_kv_num(out, "fraction", flag.fraction as f64);
        out.push('}');
    }
    out.push(']');
}

fn push_top_clusters(
    out: &mut String,
    values: &[crate::metrics::genome_stability::aggregate::TopClusterStat],
) {
    out.push('[');
    for (idx, cluster) in values.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        out.push('{');
        push_kv_str(out, "cluster_id", &cluster.cluster_id);
        out.push(',');
        push_kv_num(out, "value", cluster.value as f64);
        out.push('}');
    }
    out.push(']');
}
