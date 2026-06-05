use std::path::{Path, PathBuf};

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::input::{GeneIndex, InputError};

pub fn find_matrix_path(input_dir: &Path) -> Result<PathBuf, InputError> {
    let ds = kira_scio::discover(input_dir).map_err(|e| InputError::InvalidInput(e.message))?;
    Ok(ds.matrix)
}

/// Per-cell CSC view after feature→gene_id remapping. `cols[cell]` is sorted
/// by gene_id ascending with no duplicate gene_ids.
#[derive(Debug, Clone)]
pub struct CscMatrix {
    pub n_rows: usize,
    pub n_cols: usize,
    pub cols: Vec<Vec<(u32, f32)>>,
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

    let mut cols_vec: Vec<Vec<(u32, f32)>> = Vec::with_capacity(matrix.n_cells);

    // Without duplicate symbols the gene_ids are already increasing (just
    // filter + copy); otherwise dedup_sort_merge collapses equal gene_ids.
    for window in matrix.col_ptr.windows(2) {
        let start = window[0] as usize;
        let end = window[1] as usize;
        let mut col: Vec<(u32, f32)> = Vec::with_capacity(end - start);

        for idx in start..end {
            let feature_idx = matrix.row_idx[idx] as usize;
            let val = matrix.values[idx];
            if val == 0.0 {
                continue;
            }
            if let Some(gene_id) = gene_index
                .gene_id_by_feature
                .get(feature_idx)
                .and_then(|v| *v)
            {
                col.push((gene_id as u32, val));
            }
        }

        if gene_index.has_duplicates {
            dedup_sort_merge(&mut col);
        }
        cols_vec.push(col);
    }

    Ok(CscMatrix {
        n_rows: matrix.n_genes,
        n_cols: matrix.n_cells,
        cols: cols_vec,
    })
}

/// Sort by gene_id ascending, then in-place merge equal gene_ids by summing.
pub(crate) fn dedup_sort_merge(col: &mut Vec<(u32, f32)>) {
    if col.len() <= 1 {
        return;
    }
    col.sort_by_key(|&(g, _)| g);
    let mut write = 0usize;
    for read in 1..col.len() {
        if col[read].0 == col[write].0 {
            col[write].1 += col[read].1;
        } else {
            write += 1;
            col[write] = col[read];
        }
    }
    col.truncate(write + 1);
}
