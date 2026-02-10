
use super::*;
use crate::model::drivers::ScoreDrivers;
use crate::model::flags::Flag;
use crate::model::regimes::NuclearRegime;
use crate::model::scores::CompositeScores;
use crate::panels::{Panel, PanelAudit, PanelScores, PanelSet};
use std::sync::atomic::{AtomicUsize, Ordering};

static DIR_COUNTER: AtomicUsize = AtomicUsize::new(0);

fn make_temp_dir() -> std::path::PathBuf {
    let mut dir = std::env::temp_dir();
    let id = DIR_COUNTER.fetch_add(1, Ordering::SeqCst);
    dir.push(format!("kira_report_test_{}_{}", std::process::id(), id));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn build_input() -> Stage7Input<'static> {
    let barcodes = vec!["c1".to_string(), "c2".to_string()];
    let sample = vec!["s1".to_string(), "s1".to_string()];
    let condition = vec!["a".to_string(), "a".to_string()];
    let species = vec!["Human".to_string(), "Human".to_string()];

    let libsize = vec![10.0, 20.0];
    let nnz = vec![1u32, 2u32];
    let expr = vec![5u32, 6u32];

    let axes_tbi = vec![0.1, 0.2];
    let axes_rci = vec![0.2, 0.3];
    let axes_pds = vec![0.3, 0.4];
    let axes_trs = vec![0.4, 0.5];
    let axes_nsai = vec![0.1, 0.2];
    let axes_iaa = vec![0.1, 0.2];
    let axes_dfa = vec![0.1, 0.2];
    let axes_cea = vec![0.1, 0.2];

    let scores = CompositeScores {
        nps: vec![0.1, 0.2],
        ci: vec![0.2, 0.3],
        rls: vec![0.3, 0.4],
        confidence: vec![0.9, 0.8],
        confidence_breakdown: vec![[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]],
    };
    let drivers = ScoreDrivers {
        nps: vec![
            vec![("high_tbi".to_string(), 0.1)],
            vec![("high_tbi".to_string(), 0.2)],
        ],
        ci: vec![
            vec![("high_trs".to_string(), 0.2)],
            vec![("high_trs".to_string(), 0.3)],
        ],
        rls: vec![
            vec![("high_rci".to_string(), 0.3)],
            vec![("high_rci".to_string(), 0.4)],
        ],
    };

    let classifications = vec![
        crate::pipeline::stage6_classify::Classification {
            regime: NuclearRegime::PlasticAdaptive,
            flags: vec![Flag::LowConfidence],
        },
        crate::pipeline::stage6_classify::Classification {
            regime: NuclearRegime::Unclassified,
            flags: vec![],
        },
    ];

    let panels = PanelSet {
        panels: vec![Panel {
            id: "p1",
            name: "P1",
            group: crate::panels::defs::PanelGroup::Program,
            genes: vec![0],
            missing: vec![],
        }],
    };
    let panel_audits = vec![PanelAudit {
        panel_id: "p1".to_string(),
        panel_size_defined: 1,
        panel_size_mappable: 1,
        missing_genes: vec![],
    }];
    let panel_scores = PanelScores {
        panel_sum: vec![vec![1.0], vec![2.0]],
        panel_detected: vec![vec![1], vec![1]],
        panel_coverage: vec![vec![1.0], vec![1.0]],
    };

    Stage7Input {
        barcodes: Box::leak(Box::new(barcodes)),
        sample: Some(Box::leak(Box::new(sample))),
        condition: Some(Box::leak(Box::new(condition))),
        species_per_cell: Some(Box::leak(Box::new(species))),
        species_global: "Human".to_string(),

        libsize: Box::leak(Box::new(libsize)),
        nnz: Box::leak(Box::new(nnz)),
        expressed_genes: Box::leak(Box::new(expr)),

        axes_tbi: Box::leak(Box::new(axes_tbi)),
        axes_rci: Box::leak(Box::new(axes_rci)),
        axes_pds: Box::leak(Box::new(axes_pds)),
        axes_trs: Box::leak(Box::new(axes_trs)),
        axes_nsai: Box::leak(Box::new(axes_nsai)),
        axes_iaa: Box::leak(Box::new(axes_iaa)),
        axes_dfa: Box::leak(Box::new(axes_dfa)),
        axes_cea: Box::leak(Box::new(axes_cea)),

        scores: Box::leak(Box::new(scores)),
        drivers: Box::leak(Box::new(drivers)),

        classifications: Box::leak(Box::new(classifications)),

        panel_set: Box::leak(Box::new(panels)),
        panel_audits: Box::leak(Box::new(panel_audits)),
        panel_scores: Box::leak(Box::new(panel_scores)),

        tool_name: "kira-nuclearqc".to_string(),
        tool_version: "0.1.0".to_string(),
        git_hash: None,
        simd_backend: "scalar".to_string(),

        n_genes_raw: 10,
        n_genes_mappable: 8,

        normalize: true,
        scale: 10000.0,
        log1p: true,
        activation_mode: "Hybrid".to_string(),
        confidence_breakdown: None,
        scoring_mode: "immune-aware (default)".to_string(),
        pipeline_context: None,
    }
}

#[test]
fn test_cell_tsv_header() {
    let input = build_input();
    let dir = make_temp_dir();
    write_reports(&input, &dir, ReportMode::Cell).unwrap();
    let text = std::fs::read_to_string(dir.join("nuclearqc.tsv")).unwrap();
    let header = text.lines().next().unwrap();
    assert!(header.starts_with("barcode\tsample\tcondition\tspecies\tlibsize"));
}

#[test]
fn test_json_schema() {
    let input = build_input();
    let dir = make_temp_dir();
    write_reports(&input, &dir, ReportMode::Cell).unwrap();
    let text = std::fs::read_to_string(dir.join("summary.json")).unwrap();
    assert!(text.contains("\"tool\""));
    assert!(text.contains("\"input\""));
    assert!(text.contains("\"qc\""));
    assert!(text.contains("\"regimes\""));
}

#[test]
fn test_deterministic_output() {
    let input = build_input();
    let dir = make_temp_dir();
    write_reports(&input, &dir, ReportMode::Cell).unwrap();
    let a = std::fs::read_to_string(dir.join("summary.json")).unwrap();
    write_reports(&input, &dir, ReportMode::Cell).unwrap();
    let b = std::fs::read_to_string(dir.join("summary.json")).unwrap();
    assert_eq!(a, b);
}

#[test]
fn test_pipeline_step_json_schema_and_determinism() {
    let mut input = build_input();
    input.pipeline_context = Some(PipelineContext {
        input_dir: "/tmp/input".to_string(),
        input_source: "10x".to_string(),
        shared_bin: None,
        run_mode: "pipeline".to_string(),
    });

    let dir = make_temp_dir();
    write_reports(&input, &dir, ReportMode::Cell).unwrap();
    let first = std::fs::read_to_string(dir.join("pipeline_step.json")).unwrap();
    assert!(first.contains("\"tool\""));
    assert!(first.contains("\"artifacts\""));
    assert!(first.contains("\"cell_metrics\""));
    assert!(first.contains("\"key_metrics\""));
    assert!(first.contains("\"mode\":\"pipeline\""));

    write_reports(&input, &dir, ReportMode::Cell).unwrap();
    let second = std::fs::read_to_string(dir.join("pipeline_step.json")).unwrap();
    assert_eq!(first, second);
}
