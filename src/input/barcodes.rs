use std::io::BufRead;
use std::path::Path;

use crate::input::InputError;
use crate::input::cache::open_maybe_gz;

/// Parse a 10x barcodes.tsv(.gz) file as a flat list, one barcode per line.
///
/// Reads directly via `open_maybe_gz` rather than `kira_scio::Reader`, which
/// requires the full bundle (matrix + features + barcodes) to be present.
pub fn parse_barcodes(path: &Path) -> Result<Vec<String>, InputError> {
    let reader = open_maybe_gz(path)?;
    let mut out = Vec::new();
    let mut buf = String::new();
    let mut reader = reader;
    loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
            break;
        }
        let trimmed = buf.trim_end_matches(['\r', '\n']);
        if trimmed.is_empty() {
            continue;
        }
        out.push(trimmed.to_string());
    }
    if out.is_empty() {
        return Err(InputError::Parse("barcodes file is empty".to_string()));
    }
    Ok(out)
}
