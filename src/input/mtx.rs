use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::input::{GeneIndex, InputError};

pub fn find_matrix_path(input_dir: &Path) -> Result<PathBuf, InputError> {
    let ds = kira_scio::discover(input_dir).map_err(|e| InputError::InvalidInput(e.message))?;
    Ok(ds.matrix)
}

#[derive(Debug, Clone)]
pub struct CscMatrix {
    pub n_rows: usize,
    pub n_cols: usize,
    pub cols: Vec<Vec<(u32, i64)>>,
}

pub fn read_mtx_csc(
    path: &Path,
    n_features_raw: usize,
    n_cells: usize,
    gene_index: &GeneIndex,
) -> Result<CscMatrix, InputError> {
    let matrix = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_matrix()
    .map_err(|e| InputError::Parse(e.message))?;

    if matrix.n_genes != n_features_raw {
        return Err(InputError::InvalidInput(format!(
            "matrix row count {} does not match features {}",
            matrix.n_genes, n_features_raw
        )));
    }
    if matrix.n_cells != n_cells {
        return Err(InputError::InvalidInput(format!(
            "matrix column count {} does not match barcodes {}",
            matrix.n_cells, n_cells
        )));
    }

    let mut per_col: Vec<BTreeMap<u32, i64>> = vec![BTreeMap::new(); matrix.n_cells];

    for (col_idx, window) in matrix.col_ptr.windows(2).enumerate() {
        let start = window[0];
        let end = window[1];
        for idx in start..end {
            let feature_idx = matrix.row_idx[idx];
            let val_f = matrix.values[idx];
            if val_f == 0.0 {
                continue;
            }
            let val = val_f as i64;
            if let Some(gene_id) = gene_index
                .gene_id_by_feature
                .get(feature_idx)
                .and_then(|v| *v)
            {
                let entry = per_col[col_idx].entry(gene_id as u32).or_insert(0);
                *entry += val;
            }
        }
    }

    let mut cols_vec: Vec<Vec<(u32, i64)>> = Vec::with_capacity(matrix.n_cells);
    for map in per_col {
        let mut col_vec = Vec::with_capacity(map.len());
        for (gene, v) in map {
            col_vec.push((gene, v));
        }
        cols_vec.push(col_vec);
    }

    Ok(CscMatrix {
        n_rows: matrix.n_genes,
        n_cols: matrix.n_cells,
        cols: cols_vec,
    })
}
