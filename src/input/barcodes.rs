use std::path::Path;

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use crate::input::InputError;

pub fn parse_barcodes(path: &Path) -> Result<Vec<String>, InputError> {
    let md = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_metadata()
    .map_err(|e| InputError::Parse(e.message))?;

    if md.barcodes.is_empty() {
        return Err(InputError::Parse("barcodes file is empty".to_string()));
    }

    Ok(md.barcodes)
}
