use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use crate::input::InputError;
use crate::input::cache::open_maybe_gz;

#[derive(Debug, Clone)]
pub struct CellMeta {
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

pub fn load_meta(path: &Path, barcodes: &[String]) -> Result<CellMeta, InputError> {
    let mut reader = open_maybe_gz(path)?;
    let mut buf = String::new();

    buf.clear();
    let read = reader.read_line(&mut buf)?;
    if read == 0 {
        return Err(InputError::Parse("meta file is empty".to_string()));
    }
    let header_line = buf.trim_end();
    let header_cols: Vec<String> = header_line
        .split('\t')
        .map(|s| s.trim().to_string())
        .collect();
    if header_cols.is_empty() {
        return Err(InputError::Parse("meta file header is empty".to_string()));
    }

    let mut barcode_col = 0usize;
    for (idx, name) in header_cols.iter().enumerate() {
        let lower = name.to_ascii_lowercase();
        if lower == "barcode" || lower == "barcodes" {
            barcode_col = idx;
            break;
        }
    }

    let mut columns = Vec::new();
    for (idx, name) in header_cols.iter().enumerate() {
        if idx != barcode_col {
            columns.push(name.to_string());
        }
    }

    let mut map: HashMap<String, Vec<String>> = HashMap::new();
    let mut line_no = 1usize;

    loop {
        buf.clear();
        let read = reader.read_line(&mut buf)?;
        if read == 0 {
            break;
        }
        line_no += 1;
        let line = buf.trim_end();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }
        if barcode_col >= fields.len() {
            crate::warn!(
                "meta line has no barcode column; skipping (line {})",
                line_no
            );
            continue;
        }
        let barcode = fields[barcode_col].trim().to_string();
        if barcode.is_empty() {
            crate::warn!("meta line has empty barcode; skipping (line {})", line_no);
            continue;
        }
        if map.contains_key(&barcode) {
            crate::warn!(
                "duplicate barcode in metadata; keeping first (line {}, barcode {})",
                line_no,
                barcode
            );
            continue;
        }

        let mut row = Vec::with_capacity(columns.len());
        for (idx, _name) in header_cols.iter().enumerate() {
            if idx == barcode_col {
                continue;
            }
            let value = fields.get(idx).map(|s| s.trim()).unwrap_or("");
            row.push(value.to_string());
        }
        map.insert(barcode, row);
    }

    let mut rows = Vec::with_capacity(barcodes.len());
    for bc in barcodes {
        if let Some(row) = map.get(bc) {
            rows.push(row.clone());
        } else {
            rows.push(vec![String::new(); columns.len()]);
        }
    }

    Ok(CellMeta { columns, rows })
}
