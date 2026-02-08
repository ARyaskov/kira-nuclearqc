use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use super::barcodes::parse_barcodes;
use super::features::{Feature, normalize_symbol, parse_features};
use super::meta::load_meta;
use super::{Species, build_gene_index, detect_prefix, detect_species, resolve_shared_bin};

static DIR_COUNTER: AtomicUsize = AtomicUsize::new(0);

fn make_temp_dir() -> PathBuf {
    let mut dir = std::env::temp_dir();
    let id = DIR_COUNTER.fetch_add(1, Ordering::SeqCst);
    dir.push(format!("kira_nuclearqc_test_{}_{}", std::process::id(), id));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_file(path: &Path, contents: &str) {
    let mut f = BufWriter::new(File::create(path).unwrap());
    f.write_all(contents.as_bytes()).unwrap();
}

fn write_gz(path: &Path, contents: &str) {
    let tmp_dir = make_temp_dir();
    let src = tmp_dir.join("src.txt");
    write_file(&src, contents);

    let output = Command::new("gzip")
        .arg("-c")
        .arg("-n")
        .arg(&src)
        .output()
        .unwrap();
    assert!(output.status.success());
    fs::write(path, output.stdout).unwrap();
}

#[test]
fn test_feature_parsing_v2_v3() {
    let dir = make_temp_dir();
    let v2_path = dir.join("genes.tsv");
    let v3_path = dir.join("features.tsv.gz");

    write_file(&v2_path, "GENE1\tActb\nGENE2\tGapdh\n");
    write_gz(
        &v3_path,
        "ENSG0001\tHLA-A\tGene Expression\nENSG0002\tHLA-B\tGene Expression\n",
    );

    let v2 = parse_features(&v2_path).unwrap();
    assert_eq!(v2.len(), 2);
    assert_eq!(v2[0].feature_type, None);

    let v3 = parse_features(&v3_path).unwrap();
    assert_eq!(v3.len(), 2);
    assert_eq!(v3[0].feature_type.as_deref(), Some("Gene Expression"));
}

#[test]
fn test_gene_symbol_normalization() {
    assert_eq!(normalize_symbol("  ensG000001.12 "), "ENSG000001");
    assert_eq!(normalize_symbol(" HLA-a "), "HLA-A");
    assert_eq!(normalize_symbol(""), "");
}

#[test]
fn test_duplicate_symbol_handling() {
    let features = vec![
        Feature {
            id: "1".to_string(),
            symbol_raw: "GeneA".to_string(),
            symbol_norm: "GENEA".to_string(),
            feature_type: None,
        },
        Feature {
            id: "2".to_string(),
            symbol_raw: "GeneA".to_string(),
            symbol_norm: "GENEA".to_string(),
            feature_type: None,
        },
    ];

    let index = build_gene_index(&features);
    assert_eq!(index.symbols_by_gene_id, vec!["GENEA".to_string()]);
    assert_eq!(index.gene_id_by_feature, vec![Some(0), Some(0)]);
}

#[test]
fn test_species_detection() {
    let human_features = vec![
        Feature {
            id: "1".to_string(),
            symbol_raw: "HLA-A".to_string(),
            symbol_norm: "HLA-A".to_string(),
            feature_type: None,
        },
        Feature {
            id: "2".to_string(),
            symbol_raw: "HLA-B".to_string(),
            symbol_norm: "HLA-B".to_string(),
            feature_type: None,
        },
        Feature {
            id: "3".to_string(),
            symbol_raw: "HLA-C".to_string(),
            symbol_norm: "HLA-C".to_string(),
            feature_type: None,
        },
        Feature {
            id: "4".to_string(),
            symbol_raw: "HLA-DRA".to_string(),
            symbol_norm: "HLA-DRA".to_string(),
            feature_type: None,
        },
        Feature {
            id: "5".to_string(),
            symbol_raw: "HLA-DRB1".to_string(),
            symbol_norm: "HLA-DRB1".to_string(),
            feature_type: None,
        },
    ];

    let mouse_features = vec![
        Feature {
            id: "1".to_string(),
            symbol_raw: "H2-K1".to_string(),
            symbol_norm: "H2-K1".to_string(),
            feature_type: None,
        },
        Feature {
            id: "2".to_string(),
            symbol_raw: "H2-D1".to_string(),
            symbol_norm: "H2-D1".to_string(),
            feature_type: None,
        },
        Feature {
            id: "3".to_string(),
            symbol_raw: "H2-AB1".to_string(),
            symbol_norm: "H2-AB1".to_string(),
            feature_type: None,
        },
        Feature {
            id: "4".to_string(),
            symbol_raw: "H2-AA".to_string(),
            symbol_norm: "H2-AA".to_string(),
            feature_type: None,
        },
        Feature {
            id: "5".to_string(),
            symbol_raw: "H2-EB1".to_string(),
            symbol_norm: "H2-EB1".to_string(),
            feature_type: None,
        },
    ];

    let unknown_features = vec![
        Feature {
            id: "1".to_string(),
            symbol_raw: "GENE1".to_string(),
            symbol_norm: "GENE1".to_string(),
            feature_type: None,
        },
        Feature {
            id: "2".to_string(),
            symbol_raw: "GENE2".to_string(),
            symbol_norm: "GENE2".to_string(),
            feature_type: None,
        },
    ];

    assert_eq!(detect_species(&human_features), Species::Human);
    assert_eq!(detect_species(&mouse_features), Species::Mouse);
    assert_eq!(detect_species(&unknown_features), Species::Unknown);
}

#[test]
fn test_metadata_join() {
    let dir = make_temp_dir();
    let meta_path = dir.join("meta.tsv");

    write_file(
        &meta_path,
        "barcode\tsample\tcondition\nAA-1\tS1\tC1\nCC-1\tS2\tC2\n",
    );

    let barcodes = vec!["AA-1".to_string(), "BB-1".to_string(), "CC-1".to_string()];
    let meta = load_meta(&meta_path, &barcodes).unwrap();

    assert_eq!(
        meta.columns,
        vec!["sample".to_string(), "condition".to_string()]
    );
    assert_eq!(meta.rows.len(), 3);
    assert_eq!(meta.rows[0], vec!["S1".to_string(), "C1".to_string()]);
    assert_eq!(meta.rows[1], vec!["".to_string(), "".to_string()]);
    assert_eq!(meta.rows[2], vec!["S2".to_string(), "C2".to_string()]);
}

#[test]
fn test_barcodes_parse_order() {
    let dir = make_temp_dir();
    let path = dir.join("barcodes.tsv");
    write_file(&path, "AA-1\nBB-1\nCC-1\n");

    let barcodes = parse_barcodes(&path).unwrap();
    assert_eq!(barcodes, vec!["AA-1", "BB-1", "CC-1"]);
}

#[test]
fn test_detect_prefix_present() {
    let dir = make_temp_dir();
    write_file(&dir.join("GSM123_matrix.mtx"), "x");
    let prefix = detect_prefix(&dir).unwrap();
    assert_eq!(prefix.as_deref(), Some("GSM123"));
}

#[test]
fn test_detect_prefix_absent() {
    let dir = make_temp_dir();
    write_file(&dir.join("matrix.mtx"), "x");
    write_file(&dir.join("features.tsv"), "x");
    write_file(&dir.join("barcodes.tsv"), "x");
    let prefix = detect_prefix(&dir).unwrap();
    assert_eq!(prefix, None);
}

#[test]
fn test_resolve_shared_bin_filename_prefixed() {
    let dir = make_temp_dir();
    write_file(&dir.join("GSM1_matrix.mtx"), "x");
    let res = resolve_shared_bin(&dir).unwrap();
    assert_eq!(res.name, "GSM1.kira-organelle.bin");
}

#[test]
fn test_resolve_shared_bin_filename_default() {
    let dir = make_temp_dir();
    write_file(&dir.join("matrix.mtx"), "x");
    let res = resolve_shared_bin(&dir).unwrap();
    assert_eq!(res.name, "kira-organelle.bin");
}
