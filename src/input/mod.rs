use std::collections::HashMap;
use std::path::{Path, PathBuf};

pub mod barcodes;
pub mod cache;
pub mod features;
pub mod meta;
pub mod mtx;
pub mod organelle_bin;

use barcodes::parse_barcodes;
use features::{Feature, parse_features};
use meta::{CellMeta, load_meta};
use mtx::find_matrix_path;
use organelle_bin::{OrganelleBin, read_organelle_bin};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Species {
    Human,
    Mouse,
    Unknown,
}

#[derive(Debug, Clone)]
pub struct GeneIndex {
    pub gene_id_by_feature: Vec<Option<usize>>,
    pub symbols_by_gene_id: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct InputBundle {
    pub mtx_path: PathBuf,
    pub features_path: PathBuf,
    pub barcodes_path: PathBuf,
    pub n_cells: usize,
    pub n_features_raw: usize,
    pub n_genes_indexed: usize,
    pub species: Species,
    pub gene_index: GeneIndex,
    pub barcodes: Vec<String>,
    pub meta: Option<CellMeta>,
    pub source: InputSourceKind,
    pub organelle: Option<OrganelleBin>,
    pub shared_bin_path: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputSourceKind {
    TenX,
    OrganelleBin,
}

#[derive(Debug)]
pub enum InputError {
    Io(std::io::Error),
    MissingInput(String),
    InvalidInput(String),
    Parse(String),
}

impl std::fmt::Display for InputError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InputError::Io(e) => write!(f, "IO error: {e}"),
            InputError::MissingInput(msg) => write!(f, "missing input: {msg}"),
            InputError::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
            InputError::Parse(msg) => write!(f, "parse error: {msg}"),
        }
    }
}

impl std::error::Error for InputError {}

impl From<std::io::Error> for InputError {
    fn from(value: std::io::Error) -> Self {
        InputError::Io(value)
    }
}

pub fn load_input(input_dir: &Path, meta_path: Option<&Path>) -> Result<InputBundle, InputError> {
    load_input_tenx(input_dir, meta_path)
}

pub fn load_input_tenx(
    input_dir: &Path,
    meta_path: Option<&Path>,
) -> Result<InputBundle, InputError> {
    let mtx_path = find_matrix_path(input_dir)?;
    let features_path = find_features_path(input_dir)?;
    let barcodes_path = find_barcodes_path(input_dir)?;

    crate::info!(
        "discovered input files: mtx={}, features={}, barcodes={}",
        mtx_path.display(),
        features_path.display(),
        barcodes_path.display()
    );

    let features = parse_features(&features_path)?;
    let n_features_raw = features.len();

    let gene_index = build_gene_index(&features);
    let n_genes_indexed = gene_index.symbols_by_gene_id.len();

    let species = detect_species(&features);

    let barcodes = parse_barcodes(&barcodes_path)?;
    let n_cells = barcodes.len();

    let meta = if let Some(path) = meta_path {
        Some(load_meta(path, &barcodes)?)
    } else {
        None
    };

    Ok(InputBundle {
        mtx_path,
        features_path,
        barcodes_path,
        n_cells,
        n_features_raw,
        n_genes_indexed,
        species,
        gene_index,
        barcodes,
        meta,
        source: InputSourceKind::TenX,
        organelle: None,
        shared_bin_path: None,
    })
}

pub fn load_input_organelle(
    input_dir: &Path,
    meta_path: Option<&Path>,
    bin_path: &Path,
) -> Result<InputBundle, InputError> {
    let bin = read_organelle_bin(bin_path)?;
    let gene_symbols = bin.genes.clone();
    let barcodes = bin.barcodes.clone();

    let features = build_features_from_symbols(&gene_symbols);
    let n_features_raw = features.len();
    let gene_index = build_gene_index(&features);
    let n_genes_indexed = gene_index.symbols_by_gene_id.len();
    let species = detect_species(&features);
    let n_cells = barcodes.len();

    let meta = if let Some(path) = meta_path {
        Some(load_meta(path, &barcodes)?)
    } else {
        None
    };

    Ok(InputBundle {
        mtx_path: input_dir.join("kira-organelle.bin"),
        features_path: input_dir.join("kira-organelle.bin"),
        barcodes_path: input_dir.join("kira-organelle.bin"),
        n_cells,
        n_features_raw,
        n_genes_indexed,
        species,
        gene_index,
        barcodes,
        meta,
        source: InputSourceKind::OrganelleBin,
        organelle: Some(bin),
        shared_bin_path: Some(bin_path.to_path_buf()),
    })
}

pub fn build_gene_index(features: &[Feature]) -> GeneIndex {
    let mut symbols_by_gene_id: Vec<String> = Vec::new();
    let mut symbol_to_gene_id: HashMap<String, usize> = HashMap::new();
    let mut gene_id_by_feature: Vec<Option<usize>> = Vec::with_capacity(features.len());

    for (idx, feature) in features.iter().enumerate() {
        if feature.symbol_norm.is_empty() {
            gene_id_by_feature.push(None);
            continue;
        }
        if let Some(existing) = symbol_to_gene_id.get(feature.symbol_norm.as_str()) {
            crate::warn!(
                "duplicate gene symbol; mapping to existing gene id: feature_index={}, symbol={}",
                idx,
                feature.symbol_norm
            );
            gene_id_by_feature.push(Some(*existing));
            continue;
        }
        let gene_id = symbols_by_gene_id.len();
        symbols_by_gene_id.push(feature.symbol_norm.clone());
        symbol_to_gene_id.insert(feature.symbol_norm.clone(), gene_id);
        gene_id_by_feature.push(Some(gene_id));
    }

    GeneIndex {
        gene_id_by_feature,
        symbols_by_gene_id,
    }
}

pub fn detect_species(features: &[Feature]) -> Species {
    const HUMAN_SYMBOLS: &[&str] = &[
        "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-E", "HLA-F",
        "HLA-G",
    ];
    const MOUSE_SYMBOLS: &[&str] = &[
        "H2-K1", "H2-D1", "H2-AB1", "H2-AA", "H2-EB1", "H2-EA", "H2-Q7", "H2-Q10", "H2-T23",
        "H2-M2",
    ];

    let mut human = 0usize;
    let mut mouse = 0usize;
    for feature in features {
        let s = feature.symbol_norm.as_str();
        if s.is_empty() {
            continue;
        }
        if HUMAN_SYMBOLS.iter().any(|&x| x == s) {
            human += 1;
        }
        if MOUSE_SYMBOLS.iter().any(|&x| x == s) {
            mouse += 1;
        }
    }

    const MIN_MATCHES: usize = 3;
    const MIN_DELTA: usize = 2;

    if human >= MIN_MATCHES && human >= mouse + MIN_DELTA {
        Species::Human
    } else if mouse >= MIN_MATCHES && mouse >= human + MIN_DELTA {
        Species::Mouse
    } else {
        Species::Unknown
    }
}

fn find_features_path(input_dir: &Path) -> Result<PathBuf, InputError> {
    let candidates = [
        "features.tsv",
        "features.tsv.gz",
        "genes.tsv",
        "genes.tsv.gz",
    ];
    for name in candidates {
        let path = input_dir.join(name);
        if path.exists() {
            return Ok(path);
        }
    }
    Err(InputError::MissingInput(
        "missing features.tsv(.gz) or genes.tsv".to_string(),
    ))
}

fn find_barcodes_path(input_dir: &Path) -> Result<PathBuf, InputError> {
    let candidates = ["barcodes.tsv", "barcodes.tsv.gz"];
    for name in candidates {
        let path = input_dir.join(name);
        if path.exists() {
            return Ok(path);
        }
    }
    Err(InputError::MissingInput(
        "missing barcodes.tsv or barcodes.tsv.gz".to_string(),
    ))
}

#[derive(Debug, Clone)]
pub struct SharedBinResolution {
    pub name: String,
    pub path: PathBuf,
    pub exists: bool,
    pub prefix: Option<String>,
}

pub fn resolve_shared_bin(input_dir: &Path) -> Result<SharedBinResolution, InputError> {
    let prefix = detect_prefix(input_dir)?;
    let name = match &prefix {
        Some(p) => format!("{}.kira-organelle.bin", p),
        None => "kira-organelle.bin".to_string(),
    };
    let path = input_dir.join(&name);
    let exists = path.exists();
    Ok(SharedBinResolution {
        name,
        path,
        exists,
        prefix,
    })
}

pub fn detect_prefix(input_dir: &Path) -> Result<Option<String>, InputError> {
    let mut prefixes = std::collections::BTreeSet::new();
    for entry in std::fs::read_dir(input_dir)? {
        let entry = entry?;
        let name = entry.file_name();
        let name = name.to_string_lossy();
        for suffix in [
            "_matrix.mtx",
            "_matrix.mtx.gz",
            "_features.tsv",
            "_features.tsv.gz",
            "_barcodes.tsv",
            "_barcodes.tsv.gz",
        ] {
            if let Some(prefix) = name.strip_suffix(suffix) {
                if !prefix.is_empty() {
                    prefixes.insert(prefix.to_string());
                }
            }
        }
    }
    if let Some(prefix) = prefixes.into_iter().next() {
        return Ok(Some(prefix));
    }
    Ok(None)
}

fn build_features_from_symbols(symbols: &[String]) -> Vec<Feature> {
    let mut out = Vec::with_capacity(symbols.len());
    for (idx, symbol) in symbols.iter().enumerate() {
        let norm = crate::input::features::normalize_symbol(symbol);
        out.push(Feature {
            id: idx.to_string(),
            symbol_raw: symbol.clone(),
            symbol_norm: norm,
            feature_type: None,
        });
    }
    out
}

#[cfg(test)]
mod tests;
