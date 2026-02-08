use std::io::BufRead;
use std::path::Path;

use crate::input::InputError;
use crate::input::cache::open_maybe_gz;

pub fn parse_barcodes(path: &Path) -> Result<Vec<String>, InputError> {
    let mut reader = open_maybe_gz(path)?;
    let mut buf = String::new();
    let mut barcodes = Vec::new();

    loop {
        buf.clear();
        let read = reader.read_line(&mut buf)?;
        if read == 0 {
            break;
        }
        let line = buf.trim_end();
        if line.is_empty() {
            continue;
        }
        barcodes.push(line.trim().to_string());
    }

    if barcodes.is_empty() {
        return Err(InputError::Parse("barcodes file is empty".to_string()));
    }

    Ok(barcodes)
}
