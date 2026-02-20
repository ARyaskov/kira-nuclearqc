use std::path::Path;

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::input::InputError;

#[derive(Debug, Clone)]
pub struct Feature {
    pub id: String,
    pub symbol_raw: String,
    pub symbol_norm: String,
    pub feature_type: Option<String>,
}

pub fn parse_features(path: &Path) -> Result<Vec<Feature>, InputError> {
    let md = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_metadata()
    .map_err(|e| InputError::Parse(e.message))?;

    if md.gene_symbols.is_empty() {
        return Err(InputError::Parse("features file is empty".to_string()));
    }

    let mut features = Vec::with_capacity(md.gene_symbols.len());
    for (idx, symbol) in md.gene_symbols.iter().enumerate() {
        let id = md
            .gene_ids
            .get(idx)
            .cloned()
            .unwrap_or_else(|| idx.to_string());
        features.push(Feature {
            id,
            symbol_raw: symbol.clone(),
            symbol_norm: normalize_symbol(symbol),
            feature_type: None,
        });
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
