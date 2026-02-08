use crate::report::{RegimeStat, ReportContext, format_f32_6};

pub fn render_report_text(ctx: &ReportContext) -> String {
    let mut out = String::new();

    out.push_str("Nuclear State & Transcriptional Plasticity Report\n");
    out.push_str("=============================================\n\n");

    out.push_str("1. Overall nuclear state\n");
    out.push_str(&format!("Nuclear scoring mode: {}\n", ctx.scoring_mode));
    out.push_str(&format!(
        "Axis activation mode: {}\n",
        ctx.axis_activation_mode
    ));
    out.push_str(&format!("Confidence model: {}\n", ctx.confidence_model));
    let dominant = dominant_regimes(&ctx.regimes);
    out.push_str(&format!("Dominant regimes: {}\n", dominant));
    out.push_str(&format!(
        "Overall state: {}\n\n",
        overall_state_label(&ctx.regimes)
    ));

    out.push_str("2. Plasticity vs commitment\n");
    out.push_str(&format!(
        "NPS median: {}\nCI median: {}\n",
        format_f32_6(ctx.nps_median),
        format_f32_6(ctx.ci_median)
    ));
    out.push_str(&format!(
        "{}\n\n",
        plasticity_statement(ctx.nps_median, ctx.ci_median)
    ));

    out.push_str("3. Stress adaptation\n");
    out.push_str(&format!("NSAI median: {}\n", format_f32_6(ctx.nsai_median)));
    out.push_str(&format!("{}\n\n", stress_statement(ctx.nsai_median)));

    out.push_str("4. Reversibility outlook\n");
    out.push_str(&format!("RLS median: {}\n", format_f32_6(ctx.rls_median)));
    if !ctx.rls_contributors_top.is_empty() {
        out.push_str(&format!(
            "RLS contributors: {}\n",
            ctx.rls_contributors_top.join(", ")
        ));
    }
    out.push_str(&format!(
        "Conclusion: {}\n\n",
        reversibility_statement(ctx.rls_median, ctx.rls_tail_fraction)
    ));

    out.push_str("5. Quality and caveats\n");
    out.push_str(&format!(
        "LOW_CONFIDENCE fraction: {}\n",
        format_f32_6(ctx.low_confidence_fraction)
    ));
    out.push_str(&format!(
        "LOW_EXPR_GENES fraction: {}\n",
        format_f32_6(ctx.low_expr_fraction)
    ));
    out.push_str(&format!(
        "AMBIENT_RNA_RISK fraction: {}\n",
        format_f32_6(ctx.ambient_rna_fraction)
    ));
    out.push_str(&format!(
        "CELL_CYCLE_CONFOUNDER fraction: {}\n",
        format_f32_6(ctx.cell_cycle_fraction)
    ));
    if ctx.immune_note {
        out.push_str("Note: Immune-like scRNA detected; using relative nuclear scoring.\n");
    }
    if ctx.immune_tail_note {
        out.push_str("High IAA/DFA/CEA tails indicate immune activation subpopulation; consider cell-type gating for GC B cells\n");
    }
    if let Some([pc, es, ax, cs]) = ctx.confidence_breakdown {
        out.push_str(&format!(
            "Confidence breakdown (median): panel_coverage={}, expr_support={}, axis_structure={}, consistency={}\n",
            format_f32_6(pc),
            format_f32_6(es),
            format_f32_6(ax),
            format_f32_6(cs)
        ));
    }

    out
}

fn dominant_regimes(regimes: &[RegimeStat]) -> String {
    let mut sorted = regimes.to_vec();
    sorted.sort_by(|a, b| {
        match b
            .fraction
            .partial_cmp(&a.fraction)
            .unwrap_or(std::cmp::Ordering::Equal)
        {
            std::cmp::Ordering::Equal => a.name.cmp(b.name),
            other => other,
        }
    });
    let top = sorted.iter().take(2);
    let mut parts = Vec::new();
    for r in top {
        parts.push(format!("{} ({})", r.name, format_f32_6(r.fraction)));
    }
    parts.join(", ")
}

fn overall_state_label(regimes: &[RegimeStat]) -> &'static str {
    let mut sorted = regimes.to_vec();
    sorted.sort_by(|a, b| {
        match b
            .fraction
            .partial_cmp(&a.fraction)
            .unwrap_or(std::cmp::Ordering::Equal)
        {
            std::cmp::Ordering::Equal => a.name.cmp(b.name),
            other => other,
        }
    });
    let top = sorted.first().map(|r| r.name).unwrap_or("Unclassified");
    match top {
        "PlasticAdaptive" => "plastic",
        "StressAdaptive" => "adaptive",
        "CommittedState" => "committed",
        "RigidDegenerative" => "rigid",
        "TranscriptionallyCollapsed" => "rigid",
        _ => "mixed",
    }
}

fn plasticity_statement(nps: f32, ci: f32) -> &'static str {
    if nps >= 0.60 && ci <= 0.40 {
        "Plasticity signal is high with low commitment."
    } else if ci >= 0.60 {
        "Commitment signal is high with reduced plasticity."
    } else {
        "Plasticity and commitment are balanced."
    }
}

fn stress_statement(nsai: f32) -> &'static str {
    if nsai >= 0.60 {
        "Stress adaptation signal is high."
    } else if nsai <= 0.40 {
        "Stress adaptation signal is low."
    } else {
        "Stress adaptation signal is moderate."
    }
}

fn reversibility_statement(rls: f32, tail_frac: f32) -> &'static str {
    if rls >= 0.60 {
        "likely reversible"
    } else if rls >= 0.40 {
        "partially reversible"
    } else if rls >= 0.20 && tail_frac > 0.10 {
        "majority low-RLS with adaptive tails"
    } else {
        "low reversibility signal"
    }
}
