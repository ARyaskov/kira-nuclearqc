use super::*;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::input::load_input;
use crate::pipeline::stage2_normalize::{Stage2Params, build_expr_accessor};

static DIR_COUNTER: AtomicUsize = AtomicUsize::new(0);

fn make_temp_dir() -> PathBuf {
    let mut dir = std::env::temp_dir();
    let id = DIR_COUNTER.fetch_add(1, Ordering::SeqCst);
    dir.push(format!(
        "kira_nuclearqc_stage3_{}_{}",
        std::process::id(),
        id
    ));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_file(path: &Path, contents: &str) {
    let mut f = BufWriter::new(File::create(path).unwrap());
    f.write_all(contents.as_bytes()).unwrap();
}

fn write_mtx(path: &Path, rows: usize, cols: usize, entries: &[(usize, usize, i64)]) {
    let mut out = String::new();
    out.push_str("%%MatrixMarket matrix coordinate integer general\n");
    out.push_str("% generated\n");
    out.push_str(&format!("{} {} {}\n", rows, cols, entries.len()));
    for (r, c, v) in entries {
        out.push_str(&format!("{} {} {}\n", r, c, v));
    }
    write_file(path, &out);
}

fn setup_bundle(
    dir: &Path,
    rows: usize,
    cols: usize,
    entries: &[(usize, usize, i64)],
) -> InputBundle {
    let mtx_path = dir.join("matrix.mtx");
    let features_path = dir.join("features.tsv");
    let barcodes_path = dir.join("barcodes.tsv");

    let mut feats = String::new();
    let symbols = ["ACTB", "GAPDH", "SOX2", "FOS", "MKI67"];
    for i in 0..rows {
        feats.push_str(&format!("G{}\t{}\tGene Expression\n", i + 1, symbols[i]));
    }
    write_file(&features_path, &feats);

    let mut bcs = String::new();
    for i in 0..cols {
        bcs.push_str(&format!("CELL-{}\n", i + 1));
    }
    write_file(&barcodes_path, &bcs);

    write_mtx(&mtx_path, rows, cols, entries);

    load_input(dir, None).unwrap()
}

#[test]
fn test_panel_scoring_simple() {
    let dir = make_temp_dir();
    let bundle = setup_bundle(&dir, 5, 2, &[(1, 1, 2), (2, 1, 1), (3, 2, 3), (4, 2, 4)]);

    let accessor = build_expr_accessor(
        &bundle,
        &Stage2Params {
            normalize: false,
            cache_normalized: false,
            cache_path: None,
        },
    )
    .unwrap();

    let output = run_stage3(&bundle, accessor.as_ref()).unwrap();
    let panels = &output.panels.panels;

    let hk_idx = panels
        .iter()
        .position(|p| p.id == "housekeeping_core")
        .unwrap();
    let tf_idx = panels.iter().position(|p| p.id == "tf_basic").unwrap();
    let stress_idx = panels
        .iter()
        .position(|p| p.id == "stress_response")
        .unwrap();

    assert_eq!(output.scores.panel_sum[0][hk_idx], 3.0);
    assert_eq!(output.scores.panel_sum[1][tf_idx], 3.0);
    assert_eq!(output.scores.panel_sum[1][stress_idx], 4.0);
}

#[test]
fn test_determinism() {
    let dir = make_temp_dir();
    let bundle = setup_bundle(&dir, 5, 2, &[(1, 1, 2), (2, 1, 1), (3, 2, 3), (4, 2, 4)]);

    let accessor = build_expr_accessor(
        &bundle,
        &Stage2Params {
            normalize: false,
            cache_normalized: false,
            cache_path: None,
        },
    )
    .unwrap();

    let a = run_stage3(&bundle, accessor.as_ref()).unwrap();
    let b = run_stage3(&bundle, accessor.as_ref()).unwrap();

    assert_eq!(a.scores.panel_sum, b.scores.panel_sum);
    assert_eq!(a.scores.panel_detected, b.scores.panel_detected);
    assert_eq!(a.scores.panel_coverage, b.scores.panel_coverage);
}
