use std::collections::BTreeMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use crate::input::cache::open_maybe_gz;
use crate::input::{GeneIndex, InputError};

pub fn find_matrix_path(input_dir: &Path) -> Result<PathBuf, InputError> {
    let plain = input_dir.join("matrix.mtx");
    if plain.exists() {
        return Ok(plain);
    }
    let gz = input_dir.join("matrix.mtx.gz");
    if gz.exists() {
        return Ok(gz);
    }
    Err(InputError::MissingInput(
        "missing matrix.mtx or matrix.mtx.gz".to_string(),
    ))
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
    let mut reader = open_maybe_gz(path)?;
    let mut buf = String::new();

    // Header
    buf.clear();
    let read = reader.read_line(&mut buf)?;
    if read == 0 {
        return Err(InputError::Parse("matrix.mtx is empty".to_string()));
    }
    let header = buf.trim_end();
    if !header.starts_with("%%MatrixMarket") {
        return Err(InputError::Parse("missing MatrixMarket header".to_string()));
    }

    // Skip comments to size line
    let (rows, cols, _nnz) = loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
            return Err(InputError::Parse("missing matrix size line".to_string()));
        }
        let line = buf.trim_end();
        if line.starts_with('%') || line.is_empty() {
            continue;
        }
        let mut parts = line.split_whitespace();
        let row_raw = parts.next();
        let col_raw = parts.next();
        let nnz_raw = parts.next();
        if row_raw.is_none() || col_raw.is_none() || nnz_raw.is_none() {
            return Err(InputError::Parse("invalid matrix size line".to_string()));
        }
        let rows: usize = row_raw
            .unwrap()
            .parse()
            .map_err(|_| InputError::Parse("invalid row count".to_string()))?;
        let cols: usize = col_raw
            .unwrap()
            .parse()
            .map_err(|_| InputError::Parse("invalid column count".to_string()))?;
        let nnz: usize = nnz_raw
            .unwrap()
            .parse()
            .map_err(|_| InputError::Parse("invalid nnz count".to_string()))?;
        break (rows, cols, nnz);
    };

    if rows != n_features_raw {
        return Err(InputError::InvalidInput(format!(
            "matrix row count {} does not match features {}",
            rows, n_features_raw
        )));
    }
    if cols != n_cells {
        return Err(InputError::InvalidInput(format!(
            "matrix column count {} does not match barcodes {}",
            cols, n_cells
        )));
    }

    let mut per_col: Vec<BTreeMap<u32, i64>> = vec![BTreeMap::new(); cols];

    let mut line_no = 0usize;
    loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
            break;
        }
        line_no += 1;
        let line = buf.trim_end();
        if line.is_empty() || line.starts_with('%') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let row_raw = parts.next();
        let col_raw = parts.next();
        let val_raw = parts.next();
        if row_raw.is_none() || col_raw.is_none() || val_raw.is_none() {
            return Err(InputError::Parse(format!(
                "invalid matrix entry at line {}",
                line_no
            )));
        }
        let row: usize = row_raw
            .unwrap()
            .parse()
            .map_err(|_| InputError::Parse("invalid row index".to_string()))?;
        let col: usize = col_raw
            .unwrap()
            .parse()
            .map_err(|_| InputError::Parse("invalid col index".to_string()))?;
        let val: i64 = val_raw
            .unwrap()
            .parse()
            .map_err(|_| InputError::Parse("invalid value".to_string()))?;
        if row == 0 || row > rows || col == 0 || col > cols {
            return Err(InputError::Parse(format!(
                "matrix entry out of bounds at line {}",
                line_no
            )));
        }
        if val == 0 {
            continue;
        }
        let feature_idx = row
            .checked_sub_signed(1)
            .expect("row is validated as 1-based and non-zero");
        let col_idx = col
            .checked_sub_signed(1)
            .expect("col is validated as 1-based and non-zero");
        if let Some(gene_id) = gene_index
            .gene_id_by_feature
            .get(feature_idx)
            .and_then(|v| *v)
        {
            let entry = per_col[col_idx].entry(gene_id as u32).or_insert(0);
            *entry += val;
        }
    }

    let mut cols_vec: Vec<Vec<(u32, i64)>> = Vec::with_capacity(cols);
    for map in per_col {
        let mut col_vec = Vec::with_capacity(map.len());
        for (gene, v) in map {
            col_vec.push((gene, v));
        }
        cols_vec.push(col_vec);
    }

    Ok(CscMatrix {
        n_rows: rows,
        n_cols: cols,
        cols: cols_vec,
    })
}
