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

/// Parse a 10x features.tsv(.gz) or genes.tsv(.gz) file standalone.
///
/// Tab-separated, accepting `id`, `id<TAB>symbol`, or
/// `id<TAB>symbol<TAB>feature_type[<TAB>...]`.
pub fn parse_features(path: &Path) -> Result<Vec<Feature>, InputError> {
    let mut reader = open_maybe_gz(path)?;
    let mut out = Vec::new();
    let mut buf = String::new();
    loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
            break;
        }
        let line = buf.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        let mut parts = line.split('\t');
        let id = parts.next().unwrap_or("").trim().to_string();
        let symbol = parts.next().unwrap_or("").trim().to_string();
        let feature_type = parts.next().map(|s| s.trim().to_string());

        // Single-column rows fall back to using id as the symbol.
        let symbol_raw = if symbol.is_empty() {
            id.clone()
        } else {
            symbol
        };
        let symbol_norm = normalize_symbol(&symbol_raw);
        out.push(Feature {
            id: if id.is_empty() {
                out.len().to_string()
            } else {
                id
            },
            symbol_raw,
            symbol_norm,
            feature_type: feature_type.filter(|s| !s.is_empty()),
        });
    }
    if out.is_empty() {
        return Err(InputError::Parse("features file is empty".to_string()));
    }
    Ok(out)
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
