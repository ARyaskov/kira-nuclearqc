use std::path::PathBuf;

use crate::input::cache::{
    CacheMeta, CachedNormalizedData, cache_path_default, hash_bytes, hash_file,
    read_normalized_cache, write_normalized_cache,
};
use crate::input::mtx::{CscMatrix, read_mtx_csc};
use crate::input::organelle_bin::OrganelleBin;
use crate::input::{GeneIndex, InputBundle, InputError, InputSourceKind};

#[derive(Debug)]
pub enum Stage2Error {
    Input(InputError),
    Cache(String),
}

impl From<InputError> for Stage2Error {
    fn from(value: InputError) -> Self {
        Stage2Error::Input(value)
    }
}

impl std::fmt::Display for Stage2Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage2Error::Input(e) => write!(f, "{e}"),
            Stage2Error::Cache(msg) => write!(f, "cache error: {msg}"),
        }
    }
}

impl std::error::Error for Stage2Error {}

pub trait ExprAccessor {
    fn n_cells(&self) -> usize;
    fn n_genes(&self) -> usize;
    fn for_cell(&self, cell: usize, f: &mut dyn FnMut(u32, f32));
    fn libsize(&self, cell: usize) -> f32;
    fn nnz(&self, cell: usize) -> u32;
}

pub struct RawCountsAccessor {
    cols: Vec<Vec<(u32, i64)>>,
    libsizes: Vec<f32>,
    nnz: Vec<u32>,
    n_genes: usize,
    normalize: bool,
    scale: f32,
}

impl ExprAccessor for RawCountsAccessor {
    fn n_cells(&self) -> usize {
        self.cols.len()
    }

    fn n_genes(&self) -> usize {
        self.n_genes
    }

    fn for_cell(&self, cell: usize, f: &mut dyn FnMut(u32, f32)) {
        let lib = self.libsizes[cell] as f64;
        for &(gene_id, count) in &self.cols[cell] {
            let value = if self.normalize {
                if lib == 0.0 {
                    0.0
                } else {
                    let scaled = (count as f64) / lib * (self.scale as f64);
                    (scaled.ln_1p()) as f32
                }
            } else {
                count as f32
            };
            f(gene_id, value);
        }
    }

    fn libsize(&self, cell: usize) -> f32 {
        self.libsizes[cell]
    }

    fn nnz(&self, cell: usize) -> u32 {
        self.nnz[cell]
    }
}

pub struct CachedNormalizedAccessor {
    cols: Vec<Vec<(u32, f32)>>,
    libsizes: Vec<f32>,
    nnz: Vec<u32>,
    n_genes: usize,
}

pub struct OrganelleCountsAccessor {
    bin: OrganelleBin,
    gene_index: GeneIndex,
    libsizes: Vec<f32>,
    nnz: Vec<u32>,
    normalize: bool,
    scale: f32,
    n_genes: usize,
}

impl ExprAccessor for OrganelleCountsAccessor {
    fn n_cells(&self) -> usize {
        self.bin.csc.n_cells
    }

    fn n_genes(&self) -> usize {
        self.n_genes
    }

    fn for_cell(&self, cell: usize, f: &mut dyn FnMut(u32, f32)) {
        let start = self.bin.csc.col_ptr[cell] as usize;
        let end = self.bin.csc.col_ptr[cell + 1] as usize;
        let lib = self.libsizes[cell] as f64;
        for idx in start..end {
            let feature = self.bin.csc.row_idx[idx] as usize;
            if let Some(gene_id) = self.gene_index.gene_id_by_feature[feature] {
                let count = self.bin.csc.values[idx] as f64;
                let value = if self.normalize {
                    if lib == 0.0 {
                        0.0
                    } else {
                        let scaled = count / lib * (self.scale as f64);
                        scaled.ln_1p() as f32
                    }
                } else {
                    count as f32
                };
                f(gene_id as u32, value);
            }
        }
    }

    fn libsize(&self, cell: usize) -> f32 {
        self.libsizes[cell]
    }

    fn nnz(&self, cell: usize) -> u32 {
        self.nnz[cell]
    }
}

impl ExprAccessor for CachedNormalizedAccessor {
    fn n_cells(&self) -> usize {
        self.cols.len()
    }

    fn n_genes(&self) -> usize {
        self.n_genes
    }

    fn for_cell(&self, cell: usize, f: &mut dyn FnMut(u32, f32)) {
        for &(gene_id, value) in &self.cols[cell] {
            f(gene_id, value);
        }
    }

    fn libsize(&self, cell: usize) -> f32 {
        self.libsizes[cell]
    }

    fn nnz(&self, cell: usize) -> u32 {
        self.nnz[cell]
    }
}

#[derive(Debug, Clone)]
pub struct Stage2Params {
    pub normalize: bool,
    pub cache_normalized: bool,
    pub cache_path: Option<PathBuf>,
}

pub fn build_expr_accessor(
    bundle: &InputBundle,
    params: &Stage2Params,
) -> Result<Box<dyn ExprAccessor>, Stage2Error> {
    let scale = 10_000f32;
    let normalize = params.normalize;

    if bundle.source == InputSourceKind::OrganelleBin {
        let bin = bundle
            .organelle
            .as_ref()
            .ok_or_else(|| InputError::InvalidInput("missing organelle bin".to_string()))?
            .clone();
        let n_genes = bundle.gene_index.symbols_by_gene_id.len();

        if normalize && params.cache_normalized {
            let meta = build_cache_meta_organelle(bundle, &bin, scale, true)?;
            let cache_path = params
                .cache_path
                .clone()
                .unwrap_or_else(|| cache_path_default(bundle.shared_bin_path.as_deref().unwrap()));

            if let Some(cached) = read_normalized_cache(&cache_path, &meta)? {
                let accessor = CachedNormalizedAccessor {
                    cols: cached.columns,
                    libsizes: cached.libsizes,
                    nnz: cached.nnz,
                    n_genes,
                };
                return Ok(Box::new(accessor));
            }

            let (libsizes, nnz, normalized_cols) =
                normalize_organelle(&bin, &bundle.gene_index, scale);
            let data = CachedNormalizedData {
                libsizes: libsizes.clone(),
                nnz: nnz.clone(),
                columns: normalized_cols.clone(),
            };
            write_normalized_cache(&cache_path, &meta, &data)?;

            let accessor = CachedNormalizedAccessor {
                cols: normalized_cols,
                libsizes,
                nnz,
                n_genes,
            };
            return Ok(Box::new(accessor));
        }

        let (libsizes, nnz) = compute_stats_organelle(&bin, &bundle.gene_index);
        let accessor = OrganelleCountsAccessor {
            bin,
            gene_index: bundle.gene_index.clone(),
            libsizes,
            nnz,
            normalize,
            scale,
            n_genes,
        };
        return Ok(Box::new(accessor));
    }

    let csc = read_mtx_csc(
        &bundle.mtx_path,
        bundle.n_features_raw,
        bundle.n_cells,
        &bundle.gene_index,
    )?;

    let n_genes = bundle.gene_index.symbols_by_gene_id.len();

    if normalize && params.cache_normalized {
        let meta = build_cache_meta(bundle, scale, true)?;
        let cache_path = params
            .cache_path
            .clone()
            .unwrap_or_else(|| cache_path_default(&bundle.mtx_path));

        if let Some(cached) = read_normalized_cache(&cache_path, &meta)? {
            let accessor = CachedNormalizedAccessor {
                cols: cached.columns,
                libsizes: cached.libsizes,
                nnz: cached.nnz,
                n_genes,
            };
            return Ok(Box::new(accessor));
        }

        let (libsizes, nnz, normalized_cols) = normalize_csc(&csc, scale);
        let data = CachedNormalizedData {
            libsizes: libsizes.clone(),
            nnz: nnz.clone(),
            columns: normalized_cols.clone(),
        };
        write_normalized_cache(&cache_path, &meta, &data)?;

        let accessor = CachedNormalizedAccessor {
            cols: normalized_cols,
            libsizes,
            nnz,
            n_genes,
        };
        return Ok(Box::new(accessor));
    }

    let (libsizes, nnz) = compute_stats(&csc);

    let accessor = RawCountsAccessor {
        cols: csc.cols,
        libsizes,
        nnz,
        n_genes,
        normalize,
        scale,
    };
    Ok(Box::new(accessor))
}

fn compute_stats(csc: &CscMatrix) -> (Vec<f32>, Vec<u32>) {
    let mut libsizes = Vec::with_capacity(csc.n_cols);
    let mut nnz = Vec::with_capacity(csc.n_cols);
    for col in &csc.cols {
        let mut sum = 0f64;
        for &(_, v) in col {
            sum += v as f64;
        }
        libsizes.push(sum as f32);
        nnz.push(col.len() as u32);
    }
    (libsizes, nnz)
}

fn normalize_csc(csc: &CscMatrix, scale: f32) -> (Vec<f32>, Vec<u32>, Vec<Vec<(u32, f32)>>) {
    let mut libsizes = Vec::with_capacity(csc.n_cols);
    let mut nnz = Vec::with_capacity(csc.n_cols);
    let mut out_cols: Vec<Vec<(u32, f32)>> = Vec::with_capacity(csc.n_cols);

    for col in &csc.cols {
        let mut sum = 0f64;
        for &(_, v) in col {
            sum += v as f64;
        }
        let lib = sum;
        libsizes.push(lib as f32);
        nnz.push(col.len() as u32);

        let mut out_col = Vec::with_capacity(col.len());
        if lib == 0.0 {
            for &(gene, _) in col {
                out_col.push((gene, 0.0));
            }
        } else {
            let denom = lib;
            for &(gene, v) in col {
                let scaled = (v as f64) / denom * (scale as f64);
                let val = scaled.ln_1p() as f32;
                out_col.push((gene, val));
            }
        }
        out_cols.push(out_col);
    }

    (libsizes, nnz, out_cols)
}

fn build_cache_meta(
    bundle: &InputBundle,
    scale: f32,
    log1p: bool,
) -> Result<CacheMeta, InputError> {
    let hash_mtx = hash_file(&bundle.mtx_path)?;
    let hash_features = hash_file(&bundle.features_path)?;
    let hash_barcodes = hash_file(&bundle.barcodes_path)?;
    let hash_gene_index = hash_gene_index(&bundle.gene_index);

    Ok(CacheMeta {
        n_cells: bundle.n_cells as u32,
        n_genes: bundle.gene_index.symbols_by_gene_id.len() as u32,
        hash_mtx,
        hash_features,
        hash_barcodes,
        hash_gene_index,
        scale,
        log1p,
    })
}

fn build_cache_meta_organelle(
    bundle: &InputBundle,
    bin: &OrganelleBin,
    scale: f32,
    log1p: bool,
) -> Result<CacheMeta, InputError> {
    let bin_path = bundle
        .shared_bin_path
        .as_ref()
        .ok_or_else(|| InputError::InvalidInput("missing shared bin path".to_string()))?;
    let hash_bin = hash_file(bin_path)?;
    let hash_gene_index = hash_bytes(bundle.gene_index.symbols_by_gene_id.join("|").as_bytes());

    Ok(CacheMeta {
        n_cells: bin.csc.n_cells as u32,
        n_genes: bundle.gene_index.symbols_by_gene_id.len() as u32,
        hash_mtx: hash_bin,
        hash_features: hash_bin,
        hash_barcodes: hash_bin,
        hash_gene_index,
        scale,
        log1p,
    })
}

fn compute_stats_organelle(bin: &OrganelleBin, gene_index: &GeneIndex) -> (Vec<f32>, Vec<u32>) {
    let n_cells = bin.csc.n_cells;
    let mut libsizes = vec![0f32; n_cells];
    let mut nnz = vec![0u32; n_cells];
    for cell in 0..n_cells {
        let start = bin.csc.col_ptr[cell] as usize;
        let end = bin.csc.col_ptr[cell + 1] as usize;
        let mut sum = 0f64;
        let mut count = 0u32;
        for idx in start..end {
            let feature = bin.csc.row_idx[idx] as usize;
            if gene_index.gene_id_by_feature[feature].is_some() {
                sum += bin.csc.values[idx] as f64;
                count += 1;
            }
        }
        libsizes[cell] = sum as f32;
        nnz[cell] = count;
    }
    (libsizes, nnz)
}

fn normalize_organelle(
    bin: &OrganelleBin,
    gene_index: &GeneIndex,
    scale: f32,
) -> (Vec<f32>, Vec<u32>, Vec<Vec<(u32, f32)>>) {
    let n_cells = bin.csc.n_cells;
    let mut libsizes = vec![0f32; n_cells];
    let mut nnz = vec![0u32; n_cells];
    let mut out_cols = Vec::with_capacity(n_cells);

    for cell in 0..n_cells {
        let start = bin.csc.col_ptr[cell] as usize;
        let end = bin.csc.col_ptr[cell + 1] as usize;
        let mut sum = 0f64;
        for idx in start..end {
            let feature = bin.csc.row_idx[idx] as usize;
            if gene_index.gene_id_by_feature[feature].is_some() {
                sum += bin.csc.values[idx] as f64;
            }
        }
        let lib = sum;
        libsizes[cell] = lib as f32;

        let mut out_col = Vec::new();
        for idx in start..end {
            let feature = bin.csc.row_idx[idx] as usize;
            if let Some(gene_id) = gene_index.gene_id_by_feature[feature] {
                let count = bin.csc.values[idx] as f64;
                let val = if lib == 0.0 {
                    0.0
                } else {
                    let scaled = count / lib * (scale as f64);
                    scaled.ln_1p() as f32
                };
                out_col.push((gene_id as u32, val));
            }
        }
        nnz[cell] = out_col.len() as u32;
        out_cols.push(out_col);
    }

    (libsizes, nnz, out_cols)
}

fn hash_gene_index(index: &GeneIndex) -> u64 {
    let mut data = Vec::new();
    for symbol in &index.symbols_by_gene_id {
        data.extend_from_slice(symbol.as_bytes());
        data.push(0);
    }
    hash_bytes(&data)
}

#[cfg(test)]
mod tests {
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
}
