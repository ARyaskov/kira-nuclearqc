use std::io::BufRead;
use std::path::Path;

use crate::input::InputError;
use crate::input::cache::open_maybe_gz;

#[derive(Debug, Clone)]
pub struct Feature {
    pub id: String,
    pub symbol_raw: String,
    pub symbol_norm: String,
    pub feature_type: Option<String>,
}

pub fn parse_features(path: &Path) -> Result<Vec<Feature>, InputError> {
    let mut reader = open_maybe_gz(path)?;
    let mut buf = String::new();
    let mut features = Vec::new();
    let mut line_no = 0usize;
    let mut format_cols: Option<usize> = None;

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
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            return Err(InputError::Parse(format!(
                "features line {} has <2 columns",
                line_no
            )));
        }
        if format_cols.is_none() {
            format_cols = Some(cols.len());
        } else if let Some(expected) = format_cols {
            if expected == 2 && cols.len() >= 3 {
                crate::warn!(
                    "features file appears to be v3 even though first line looked v2 (line {})",
                    line_no
                );
            }
        }
        let id = cols[0].trim().to_string();
        let symbol_raw = cols[1].trim().to_string();
        let symbol_norm = normalize_symbol(&symbol_raw);
        let feature_type = if cols.len() >= 3 {
            Some(cols[2].trim().to_string())
        } else {
            None
        };
        features.push(Feature {
            id,
            symbol_raw,
            symbol_norm,
            feature_type,
        });
    }

    if features.is_empty() {
        return Err(InputError::Parse("features file is empty".to_string()));
    }

    Ok(features)
}

pub fn normalize_symbol(raw: &str) -> String {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return String::new();
    }
    let upper = trimmed.to_ascii_uppercase();
    if let Some((left, right)) = upper.rsplit_once('.') {
        if left.starts_with("ENS") && right.chars().all(|c| c.is_ascii_digit()) {
            return left.to_string();
        }
    }
    upper
}
