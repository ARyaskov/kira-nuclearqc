mod input;
mod model;
mod panels;
mod pipeline;
mod report;
mod simd;
mod tracing;

use std::path::{Path, PathBuf};

use crate::input::{load_input_organelle, load_input_tenx, resolve_shared_bin};
use crate::model::thresholds::{NuclearScoringMode, ThresholdProfile};
use crate::pipeline::stage2_normalize::{Stage2Params, build_expr_accessor};
use crate::pipeline::stage3_panels::run_stage3;
use crate::pipeline::stage4_axes::run_stage4;
use crate::pipeline::stage5_scores::{Stage5Inputs, run_stage5};
use crate::pipeline::stage6_classify::{Stage6Inputs, run_stage6};
use crate::pipeline::stage7_report::{
    PipelineContext, ReportMode, RunMode, Stage7Input, write_reports,
};
use crate::report::p90;

fn main() {
    println!("SIMD backend: {}", simd::backend_name());
    if let Err(err) = run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let config = parse_args(&args)?;

    let out_dir = resolve_output_dir(&config.out_dir, config.run_mode);

    let (bundle, input_source, shared_bin) = match config.run_mode {
        RunMode::Standalone => (
            load_input_tenx(&config.input_dir, config.meta_path.as_deref())
                .map_err(|e| e.to_string())?,
            "10x".to_string(),
            None,
        ),
        RunMode::Pipeline => {
            let resolution = resolve_shared_bin(&config.input_dir).map_err(|e| e.to_string())?;
            if resolution.exists {
                (
                    load_input_organelle(
                        &config.input_dir,
                        config.meta_path.as_deref(),
                        &resolution.path,
                    )
                    .map_err(|e| e.to_string())?,
                    "kira-organelle.bin".to_string(),
                    Some(resolution.name),
                )
            } else {
                crate::warn!(
                    "--run-mode pipeline requested but shared cache file {} was not found; falling back to 10x MTX reading (slower).",
                    resolution.name
                );
                (
                    load_input_tenx(&config.input_dir, config.meta_path.as_deref())
                        .map_err(|e| e.to_string())?,
                    "10x".to_string(),
                    None,
                )
            }
        }
    };

    let stage2 = Stage2Params {
        normalize: config.normalize,
        cache_normalized: config.cache_normalized,
        cache_path: None,
    };
    let accessor = build_expr_accessor(&bundle, &stage2).map_err(|e| e.to_string())?;

    let stage3 = run_stage3(&bundle, accessor.as_ref()).map_err(|e| e.to_string())?;
    let thresholds = match config.scoring_mode {
        NuclearScoringMode::ImmuneAware => ThresholdProfile::immune_v1(),
        NuclearScoringMode::StrictBulk => ThresholdProfile::default_v1(),
    };
    let stage4 = run_stage4(
        accessor.as_ref(),
        &stage3.panels,
        &stage3.scores,
        &thresholds,
    );
    log_scoring_mode(config.scoring_mode, &stage3, &stage4);

    let key_panel_coverage_median = compute_key_panel_coverage(&stage3.panels, &stage3.scores);
    let ambient_rna_risk = vec![false; bundle.n_cells];
    let axis_p90 = [
        p90(&stage4.axes.iaa),
        p90(&stage4.axes.dfa),
        p90(&stage4.axes.nsai),
    ];

    let (program_sum, sum_tf, proliferation_share, key_panels_missing, panel_nonzero_fraction) =
        compute_panel_signals(&stage3.panels, &stage3.scores, &stage3.audits);

    let stage5 = run_stage5(&Stage5Inputs {
        axes: &stage4.axes,
        drivers: &stage4.drivers,
        thresholds: &thresholds,
        n_genes_mappable: Some(bundle.n_genes_indexed as u32),
        key_panel_coverage_median: Some(&key_panel_coverage_median),
        ambient_rna_risk: Some(&ambient_rna_risk),
        key_panels_missing: Some(&key_panels_missing),
        panel_nonzero_fraction: Some(&panel_nonzero_fraction),
        axis_p90: Some(axis_p90),
        scoring_mode: config.scoring_mode,
    });

    let stage6 = run_stage6(&Stage6Inputs {
        tbi: &stage4.axes.tbi,
        rci: &stage4.axes.rci,
        pds: &stage4.axes.pds,
        trs: &stage4.axes.trs,
        nsai: &stage4.axes.nsai,
        iaa: &stage4.axes.iaa,
        dfa: &stage4.axes.dfa,
        cea: &stage4.axes.cea,
        scores: &stage5.scores,
        drivers: &stage4.drivers,
        thresholds: &thresholds,
        scoring_mode: config.scoring_mode,
        key_panel_coverage_median: Some(&key_panel_coverage_median),
        key_panels_missing: Some(&key_panels_missing),
        sum_tf_panels: Some(&sum_tf),
        ambient_rna_risk: Some(&ambient_rna_risk),
        proliferation_program_share: Some(&proliferation_share),
        program_sum: Some(&program_sum),
    });

    let (sample, condition, species_per_cell) = extract_meta(&bundle);

    let mut libsize_vec = Vec::with_capacity(bundle.n_cells);
    let mut nnz_vec = Vec::with_capacity(bundle.n_cells);
    for cell in 0..bundle.n_cells {
        libsize_vec.push(accessor.libsize(cell));
        nnz_vec.push(accessor.nnz(cell));
    }
    let expressed_vec = stage4
        .drivers
        .iter()
        .map(|d| d.expressed_genes)
        .collect::<Vec<_>>();

    let input = Stage7Input {
        barcodes: &bundle.barcodes,
        sample: sample.as_deref(),
        condition: condition.as_deref(),
        species_per_cell: species_per_cell.as_deref(),
        species_global: format!("{:?}", bundle.species),

        libsize: &libsize_vec,
        nnz: &nnz_vec,
        expressed_genes: &expressed_vec,

        axes_tbi: &stage4.axes.tbi,
        axes_rci: &stage4.axes.rci,
        axes_pds: &stage4.axes.pds,
        axes_trs: &stage4.axes.trs,
        axes_nsai: &stage4.axes.nsai,
        axes_iaa: &stage4.axes.iaa,
        axes_dfa: &stage4.axes.dfa,
        axes_cea: &stage4.axes.cea,

        scores: &stage5.scores,
        drivers: &stage5.drivers,
        activation_mode: format!("{:?}", thresholds.activation_mode),

        classifications: &stage6,

        panel_set: &stage3.panels,
        panel_audits: &stage3.audits,
        panel_scores: &stage3.scores,

        tool_name: "kira-nuclearqc".to_string(),
        tool_version: env!("CARGO_PKG_VERSION").to_string(),
        git_hash: read_git_hash(&PathBuf::from(".")),
        simd_backend: simd::backend_name().to_string(),

        n_genes_raw: bundle.n_features_raw,
        n_genes_mappable: bundle.n_genes_indexed,

        normalize: config.normalize,
        scale: 10_000.0,
        log1p: config.normalize,
        confidence_breakdown: Some(&stage5.scores.confidence_breakdown),
        scoring_mode: match config.scoring_mode {
            NuclearScoringMode::ImmuneAware => "immune-aware (default)".to_string(),
            NuclearScoringMode::StrictBulk => "strict (bulk-oriented)".to_string(),
        },
        pipeline_context: if config.run_mode == RunMode::Pipeline {
            Some(PipelineContext {
                input_dir: config.input_dir.display().to_string(),
                input_source,
                shared_bin,
                run_mode: "pipeline".to_string(),
            })
        } else {
            None
        },
    };

    write_reports(&input, &out_dir, config.report_mode).map_err(|e| e.to_string())?;

    Ok(())
}

#[derive(Debug, Clone)]
struct RunConfig {
    input_dir: PathBuf,
    out_dir: PathBuf,
    report_mode: ReportMode,
    meta_path: Option<PathBuf>,
    normalize: bool,
    cache_normalized: bool,
    scoring_mode: NuclearScoringMode,
    run_mode: RunMode,
}

fn parse_args(args: &[String]) -> Result<RunConfig, String> {
    if args.is_empty() {
        return Err("missing command".to_string());
    }
    let mut args = args.to_vec();
    let cmd = args.remove(0);
    if cmd != "run" {
        return Err("unsupported command".to_string());
    }

    let mut input_dir: Option<PathBuf> = None;
    let mut out_dir: Option<PathBuf> = None;
    let mut report_mode = ReportMode::Cell;
    let mut meta_path: Option<PathBuf> = None;
    let mut normalize = false;
    let mut cache_normalized = false;
    let mut scoring_mode = NuclearScoringMode::ImmuneAware;
    let mut run_mode = RunMode::Standalone;

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "--input" => {
                i += 1;
                if i >= args.len() {
                    return Err("missing value for --input".to_string());
                }
                input_dir = Some(PathBuf::from(&args[i]));
            }
            "--out" => {
                i += 1;
                if i >= args.len() {
                    return Err("missing value for --out".to_string());
                }
                out_dir = Some(PathBuf::from(&args[i]));
            }
            "--mode" => {
                i += 1;
                if i >= args.len() {
                    return Err("missing value for --mode".to_string());
                }
                report_mode = match args[i].as_str() {
                    "cell" => ReportMode::Cell,
                    "sample" => ReportMode::Sample,
                    _ => return Err("invalid --mode (use cell|sample)".to_string()),
                };
            }
            "--meta" => {
                i += 1;
                if i >= args.len() {
                    return Err("missing value for --meta".to_string());
                }
                meta_path = Some(PathBuf::from(&args[i]));
            }
            "--normalize" => {
                normalize = true;
            }
            "--cache-normalized" => {
                cache_normalized = true;
            }
            "--strict-nuclear" => {
                scoring_mode = NuclearScoringMode::StrictBulk;
            }
            "--run-mode" => {
                i += 1;
                if i >= args.len() {
                    return Err("missing value for --run-mode".to_string());
                }
                run_mode = match args[i].as_str() {
                    "standalone" => RunMode::Standalone,
                    "pipeline" => RunMode::Pipeline,
                    _ => return Err("invalid --run-mode (use standalone|pipeline)".to_string()),
                };
            }
            other => {
                return Err(format!("unknown argument: {}", other));
            }
        }
        i += 1;
    }

    Ok(RunConfig {
        input_dir: input_dir.ok_or_else(|| "missing --input".to_string())?,
        out_dir: out_dir.ok_or_else(|| "missing --out".to_string())?,
        report_mode,
        meta_path,
        normalize,
        cache_normalized,
        scoring_mode,
        run_mode,
    })
}

fn resolve_output_dir(base: &Path, run_mode: RunMode) -> PathBuf {
    match run_mode {
        RunMode::Standalone => base.to_path_buf(),
        RunMode::Pipeline => base.join("kira-nuclearqc"),
    }
}

fn extract_meta(
    bundle: &input::InputBundle,
) -> (
    Option<Vec<String>>,
    Option<Vec<String>>,
    Option<Vec<String>>,
) {
    let mut sample: Option<Vec<String>> = None;
    let mut condition: Option<Vec<String>> = None;
    let mut species: Option<Vec<String>> = None;
    if let Some(meta) = &bundle.meta {
        let mut sample_idx = None;
        let mut condition_idx = None;
        let mut species_idx = None;
        for (i, name) in meta.columns.iter().enumerate() {
            let lower = name.to_ascii_lowercase();
            if lower == "sample" {
                sample_idx = Some(i);
            } else if lower == "condition" {
                condition_idx = Some(i);
            } else if lower == "species" {
                species_idx = Some(i);
            }
        }
        if let Some(idx) = sample_idx {
            sample = Some(
                meta.rows
                    .iter()
                    .map(|r| r.get(idx).cloned().unwrap_or_default())
                    .collect(),
            );
        }
        if let Some(idx) = condition_idx {
            condition = Some(
                meta.rows
                    .iter()
                    .map(|r| r.get(idx).cloned().unwrap_or_default())
                    .collect(),
            );
        }
        if let Some(idx) = species_idx {
            species = Some(
                meta.rows
                    .iter()
                    .map(|r| r.get(idx).cloned().unwrap_or_default())
                    .collect(),
            );
        }
    }
    (sample, condition, species)
}

fn compute_key_panel_coverage(
    panel_set: &panels::PanelSet,
    scores: &panels::PanelScores,
) -> Vec<f32> {
    let n_cells = scores.panel_coverage.len();
    let n_panels = panel_set.panels.len();
    let mut out = Vec::with_capacity(n_cells);
    for cell in 0..n_cells {
        if n_panels == 0 {
            out.push(0.0);
            continue;
        }
        let mut values = Vec::with_capacity(n_panels);
        for p in 0..n_panels {
            values.push(scores.panel_coverage[cell][p]);
        }
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let idx = values.len() / 2;
        out.push(values[idx]);
    }
    out
}

fn compute_panel_signals(
    panel_set: &panels::PanelSet,
    scores: &panels::PanelScores,
    audits: &[panels::PanelAudit],
) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<bool>, Vec<f32>) {
    let n_cells = scores.panel_sum.len();
    let mut program_sum = vec![0.0f32; n_cells];
    let mut tf_sum = vec![0.0f32; n_cells];
    let mut proliferation_sum = vec![0.0f32; n_cells];
    let mut nonzero_frac = vec![0.0f32; n_cells];

    let mut key_panels_missing = false;
    for audit in audits {
        if audit.panel_size_mappable == 0 {
            key_panels_missing = true;
        }
    }

    for (idx, panel) in panel_set.panels.iter().enumerate() {
        for cell in 0..n_cells {
            let v = scores.panel_sum[cell][idx];
            match panel.group {
                panels::defs::PanelGroup::Program => program_sum[cell] += v,
                panels::defs::PanelGroup::Tf | panels::defs::PanelGroup::Chromatin => {
                    tf_sum[cell] += v
                }
                panels::defs::PanelGroup::Proliferation => proliferation_sum[cell] += v,
                _ => {}
            }
        }
    }

    for cell in 0..n_cells {
        let mut detected = 0u32;
        let mut total = 0u32;
        for (idx, panel) in panel_set.panels.iter().enumerate() {
            detected += scores.panel_detected[cell][idx];
            total += panel.genes.len() as u32;
        }
        if total > 0 {
            nonzero_frac[cell] = detected as f32 / total as f32;
        }
    }

    let mut proliferation_share = vec![0.0f32; n_cells];
    for cell in 0..n_cells {
        let denom = program_sum[cell];
        if denom > 0.0 {
            proliferation_share[cell] = proliferation_sum[cell] / denom;
        }
    }

    (
        program_sum,
        tf_sum,
        proliferation_share,
        vec![key_panels_missing; n_cells],
        nonzero_frac,
    )
}

fn read_git_hash(repo_root: &Path) -> Option<String> {
    let head = repo_root.join(".git/HEAD");
    let content = std::fs::read_to_string(head).ok()?;
    if let Some(ref_line) = content.strip_prefix("ref: ") {
        let ref_path = repo_root.join(".git").join(ref_line.trim());
        return std::fs::read_to_string(ref_path)
            .ok()
            .map(|s| s.trim().to_string());
    }
    Some(content.trim().to_string())
}

fn log_scoring_mode(
    mode: NuclearScoringMode,
    stage3: &pipeline::stage3_panels::Stage3Output,
    stage4: &pipeline::stage4_axes::Stage4Output,
) {
    match mode {
        NuclearScoringMode::ImmuneAware => {
            eprintln!("INFO  Immune-aware nuclear scoring enabled (default)");
            if immune_like_detected(stage3, stage4) {
                eprintln!("INFO  Immune-like scRNA detected; relative nuclear scoring in effect");
            }
        }
        NuclearScoringMode::StrictBulk => {
            eprintln!(
                "WARN  Strict nuclear mode enabled (--strict-nuclear); immune dynamics may be underdetected"
            );
        }
    }
}

fn immune_like_detected(
    stage3: &pipeline::stage3_panels::Stage3Output,
    stage4: &pipeline::stage4_axes::Stage4Output,
) -> bool {
    let has_immune_panels = stage3
        .panels
        .panels
        .iter()
        .any(|p| p.id == "immune_activation" || p.id == "clonal_engagement");
    let p90_iaa = p90(&stage4.axes.iaa);
    let p90_dfa = p90(&stage4.axes.dfa);
    let p90_cea = p90(&stage4.axes.cea);
    has_immune_panels && (p90_iaa > 0.5 || p90_dfa > 0.5 || p90_cea > 0.5)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_args_default_run_mode_standalone() {
        let args = vec![
            "run".to_string(),
            "--input".to_string(),
            "data".to_string(),
            "--out".to_string(),
            "out".to_string(),
        ];
        let parsed = parse_args(&args).unwrap();
        assert_eq!(parsed.run_mode, RunMode::Standalone);
    }

    #[test]
    fn test_parse_args_pipeline_run_mode() {
        let args = vec![
            "run".to_string(),
            "--input".to_string(),
            "data".to_string(),
            "--out".to_string(),
            "out".to_string(),
            "--run-mode".to_string(),
            "pipeline".to_string(),
        ];
        let parsed = parse_args(&args).unwrap();
        assert_eq!(parsed.run_mode, RunMode::Pipeline);
    }

    #[test]
    fn test_resolve_output_dir_pipeline() {
        let out = resolve_output_dir(Path::new("/tmp/out"), RunMode::Pipeline);
        assert_eq!(out, PathBuf::from("/tmp/out/kira-nuclearqc"));
    }

    #[test]
    fn test_resolve_output_dir_standalone() {
        let out = resolve_output_dir(Path::new("/tmp/out"), RunMode::Standalone);
        assert_eq!(out, PathBuf::from("/tmp/out"));
    }
}
