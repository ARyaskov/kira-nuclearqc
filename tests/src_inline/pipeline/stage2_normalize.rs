
use super::*;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::input::load_input;

static DIR_COUNTER: AtomicUsize = AtomicUsize::new(0);

fn make_temp_dir() -> PathBuf {
    let mut dir = std::env::temp_dir();
    let id = DIR_COUNTER.fetch_add(1, Ordering::SeqCst);
    dir.push(format!(
        "kira_nuclearqc_stage2_{}_{}",
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
    for i in 0..rows {
        feats.push_str(&format!("G{}\tGene{}\tGene Expression\n", i + 1, i + 1));
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
fn test_small_mtx_stats_and_values() {
    let dir = make_temp_dir();
    let bundle = setup_bundle(&dir, 2, 2, &[(1, 1, 1), (2, 1, 2), (2, 2, 3)]);

    let params = Stage2Params {
        normalize: false,
        cache_normalized: false,
        cache_path: None,
    };
    let accessor = build_expr_accessor(&bundle, &params).unwrap();

    assert_eq!(accessor.n_cells(), 2);
    assert_eq!(accessor.n_genes(), 2);
    assert_eq!(accessor.libsize(0), 3.0);
    assert_eq!(accessor.libsize(1), 3.0);
    assert_eq!(accessor.nnz(0), 2);
    assert_eq!(accessor.nnz(1), 1);

    let mut values = Vec::new();
    accessor.for_cell(0, &mut |g, v| values.push((g, v)));
    assert_eq!(values, vec![(0, 1.0), (1, 2.0)]);

    values.clear();
    accessor.for_cell(1, &mut |g, v| values.push((g, v)));
    assert_eq!(values, vec![(1, 3.0)]);
}

#[test]
fn test_normalized_vs_raw() {
    let dir = make_temp_dir();
    let bundle = setup_bundle(&dir, 2, 2, &[(1, 1, 1), (2, 1, 3)]);

    let raw = build_expr_accessor(
        &bundle,
        &Stage2Params {
            normalize: false,
            cache_normalized: false,
            cache_path: None,
        },
    )
    .unwrap();

    let norm = build_expr_accessor(
        &bundle,
        &Stage2Params {
            normalize: true,
            cache_normalized: false,
            cache_path: None,
        },
    )
    .unwrap();

    let mut raw_vals = Vec::new();
    raw.for_cell(0, &mut |g, v| raw_vals.push((g, v)));

    let mut norm_vals = Vec::new();
    norm.for_cell(0, &mut |g, v| norm_vals.push((g, v)));

    assert_eq!(raw_vals.len(), norm_vals.len());
    assert!(norm_vals[0].1 > 0.0);
}

#[test]
fn test_cache_round_trip() {
    let dir = make_temp_dir();
    let bundle = setup_bundle(&dir, 2, 2, &[(1, 1, 1), (2, 1, 2), (2, 2, 3)]);

    let cache_path = dir.join("cache.bin");
    let params = Stage2Params {
        normalize: true,
        cache_normalized: true,
        cache_path: Some(cache_path.clone()),
    };
    let accessor_a = build_expr_accessor(&bundle, &params).unwrap();
    let accessor_b = build_expr_accessor(&bundle, &params).unwrap();

    let mut a_vals = Vec::new();
    accessor_a.for_cell(0, &mut |g, v| a_vals.push((g, v)));

    let mut b_vals = Vec::new();
    accessor_b.for_cell(0, &mut |g, v| b_vals.push((g, v)));

    assert_eq!(a_vals, b_vals);
}

#[test]
fn test_determinism_bitwise() {
    let dir = make_temp_dir();
    let bundle = setup_bundle(&dir, 3, 3, &[(1, 1, 1), (2, 1, 2), (3, 2, 3), (1, 3, 4)]);

    let params = Stage2Params {
        normalize: true,
        cache_normalized: false,
        cache_path: None,
    };
    let a = build_expr_accessor(&bundle, &params).unwrap();
    let b = build_expr_accessor(&bundle, &params).unwrap();

    for cell in 0..3 {
        let mut av = Vec::new();
        let mut bv = Vec::new();
        a.for_cell(cell, &mut |g, v| av.push((g, v.to_bits())));
        b.for_cell(cell, &mut |g, v| bv.push((g, v.to_bits())));
        assert_eq!(av, bv);
    }
}
