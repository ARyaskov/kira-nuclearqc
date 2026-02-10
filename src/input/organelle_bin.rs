use std::fs::File;
use std::path::Path;

use memmap2::Mmap;

use crate::input::InputError;

#[derive(Debug, Clone)]
pub struct HeaderV1 {
    pub n_genes: u64,
    pub n_cells: u64,
    pub nnz: u64,
    pub genes_table_offset: u64,
    pub genes_table_bytes: u64,
    pub barcodes_table_offset: u64,
    pub barcodes_table_bytes: u64,
    pub col_ptr_offset: u64,
    pub row_idx_offset: u64,
    pub values_u32_offset: u64,
    pub n_blocks: u64,
    pub blocks_offset: u64,
    pub file_bytes: u64,
    pub header_crc64: u64,
    pub data_crc64: u64,
}

#[derive(Debug, Clone)]
pub struct CscView {
    pub n_genes: usize,
    pub n_cells: usize,
    pub nnz: usize,
    pub col_ptr: Vec<u64>,
    pub row_idx: Vec<u32>,
    pub values: Vec<u32>,
}

#[derive(Debug, Clone)]
pub struct OrganelleBin {
    pub header: HeaderV1,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
    pub csc: CscView,
}

pub fn read_organelle_bin(path: &Path) -> Result<OrganelleBin, InputError> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let bytes = &mmap[..];
    if bytes.len() < 256 {
        return Err(InputError::InvalidInput(
            "kira-organelle.bin too small".to_string(),
        ));
    }

    let header = parse_header(&bytes)?;
    validate_header(&header, &bytes)?;

    let genes = parse_string_table(
        &bytes,
        header.genes_table_offset,
        header.genes_table_bytes,
        header.n_genes as usize,
    )?;
    let barcodes = parse_string_table(
        &bytes,
        header.barcodes_table_offset,
        header.barcodes_table_bytes,
        header.n_cells as usize,
    )?;

    let col_ptr = read_u64_vec(&bytes, header.col_ptr_offset, header.n_cells as usize + 1)?;
    let row_idx = read_u32_vec(&bytes, header.row_idx_offset, header.nnz as usize)?;
    let values = read_u32_vec(&bytes, header.values_u32_offset, header.nnz as usize)?;

    let csc = CscView {
        n_genes: header.n_genes as usize,
        n_cells: header.n_cells as usize,
        nnz: header.nnz as usize,
        col_ptr,
        row_idx,
        values,
    };
    validate_csc(&csc)?;

    Ok(OrganelleBin {
        header,
        genes,
        barcodes,
        csc,
    })
}

fn parse_header(bytes: &[u8]) -> Result<HeaderV1, InputError> {
    if &bytes[0..4] != b"KORG" {
        return Err(InputError::InvalidInput(
            "invalid magic; expected KORG".to_string(),
        ));
    }
    let version_major = read_u16(bytes, 4);
    let version_minor = read_u16(bytes, 6);
    if version_major != 1 || version_minor != 0 {
        return Err(InputError::InvalidInput(format!(
            "unsupported version: {}.{}",
            version_major, version_minor
        )));
    }
    if read_u32(bytes, 8) != 0x1234_5678 {
        return Err(InputError::InvalidInput(
            "unsupported endianness tag".to_string(),
        ));
    }
    if read_u32(bytes, 12) != 256 {
        return Err(InputError::InvalidInput(
            "invalid header_size; expected 256".to_string(),
        ));
    }

    Ok(HeaderV1 {
        n_genes: read_u64(bytes, 16),
        n_cells: read_u64(bytes, 24),
        nnz: read_u64(bytes, 32),
        genes_table_offset: read_u64(bytes, 40),
        genes_table_bytes: read_u64(bytes, 48),
        barcodes_table_offset: read_u64(bytes, 56),
        barcodes_table_bytes: read_u64(bytes, 64),
        col_ptr_offset: read_u64(bytes, 72),
        row_idx_offset: read_u64(bytes, 80),
        values_u32_offset: read_u64(bytes, 88),
        n_blocks: read_u64(bytes, 96),
        blocks_offset: read_u64(bytes, 104),
        file_bytes: read_u64(bytes, 112),
        header_crc64: read_u64(bytes, 120),
        data_crc64: read_u64(bytes, 128),
    })
}

fn validate_header(header: &HeaderV1, bytes: &[u8]) -> Result<(), InputError> {
    if header.file_bytes as usize != bytes.len() {
        return Err(InputError::InvalidInput(
            "file_bytes does not match file length".to_string(),
        ));
    }
    if header.n_blocks != 0 {
        return Err(InputError::InvalidInput(
            "unsupported optional blocks in header".to_string(),
        ));
    }
    if header.blocks_offset != 0 {
        return Err(InputError::InvalidInput(
            "blocks_offset must be zero when n_blocks=0".to_string(),
        ));
    }
    if header.data_crc64 != 0 {
        return Err(InputError::InvalidInput(
            "data_crc64 not supported in v1".to_string(),
        ));
    }
    if header.header_crc64 == 0 {
        return Err(InputError::InvalidInput(
            "header_crc64 missing or zero".to_string(),
        ));
    }
    let mut hdr = bytes[0..256].to_vec();
    hdr[120..128].fill(0);
    let crc = crc64_ecma(&hdr);
    if crc != header.header_crc64 {
        return Err(InputError::InvalidInput(
            "header_crc64 mismatch".to_string(),
        ));
    }
    Ok(())
}

fn parse_string_table(
    bytes: &[u8],
    offset: u64,
    table_bytes: u64,
    expected_count: usize,
) -> Result<Vec<String>, InputError> {
    let offset = offset as usize;
    let table_bytes = table_bytes as usize;
    if offset + table_bytes > bytes.len() || table_bytes < 8 {
        return Err(InputError::InvalidInput(
            "string table out of bounds".to_string(),
        ));
    }
    let tbl = &bytes[offset..offset + table_bytes];
    let count = read_u32(tbl, 0) as usize;
    if count != expected_count {
        return Err(InputError::InvalidInput(
            "string table count mismatch".to_string(),
        ));
    }
    let offsets_bytes = (count + 1)
        .checked_mul(4)
        .ok_or_else(|| InputError::InvalidInput("offset overflow".to_string()))?;
    if 4 + offsets_bytes > tbl.len() {
        return Err(InputError::InvalidInput(
            "string table offsets out of bounds".to_string(),
        ));
    }
    let blob = &tbl[(4 + offsets_bytes)..];
    let mut offsets = Vec::with_capacity(count + 1);
    for i in 0..=count {
        offsets.push(read_u32(tbl, 4 + i * 4) as usize);
    }
    if offsets[count] > blob.len() {
        return Err(InputError::InvalidInput(
            "string table blob length mismatch".to_string(),
        ));
    }
    for i in 0..count {
        if offsets[i] > offsets[i + 1] {
            return Err(InputError::InvalidInput(
                "string table offsets not monotonic".to_string(),
            ));
        }
    }

    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let start = offsets[i];
        let end = offsets[i + 1];
        let s = std::str::from_utf8(&blob[start..end])
            .map_err(|_| InputError::InvalidInput("invalid utf-8 in string table".to_string()))?;
        out.push(s.to_string());
    }
    Ok(out)
}

fn validate_csc(csc: &CscView) -> Result<(), InputError> {
    if csc.col_ptr.len() != csc.n_cells + 1 {
        return Err(InputError::InvalidInput(
            "col_ptr length mismatch".to_string(),
        ));
    }
    if csc.row_idx.len() != csc.nnz || csc.values.len() != csc.nnz {
        return Err(InputError::InvalidInput(
            "row_idx/values length mismatch".to_string(),
        ));
    }
    if csc.col_ptr[0] != 0 {
        return Err(InputError::InvalidInput("col_ptr[0] must be 0".to_string()));
    }
    if csc.col_ptr[csc.n_cells] != csc.nnz as u64 {
        return Err(InputError::InvalidInput(
            "col_ptr[n_cells] must equal nnz".to_string(),
        ));
    }
    for i in 1..csc.col_ptr.len() {
        if csc.col_ptr[i - 1] > csc.col_ptr[i] {
            return Err(InputError::InvalidInput(
                "col_ptr not monotonic".to_string(),
            ));
        }
    }
    for cell in 0..csc.n_cells {
        let start = csc.col_ptr[cell] as usize;
        let end = csc.col_ptr[cell + 1] as usize;
        if start > end || end > csc.nnz {
            return Err(InputError::InvalidInput(
                "col_ptr bounds invalid".to_string(),
            ));
        }
        let mut prev = None;
        for k in start..end {
            let row = csc.row_idx[k] as usize;
            if row >= csc.n_genes {
                return Err(InputError::InvalidInput(
                    "row_idx out of bounds".to_string(),
                ));
            }
            if let Some(p) = prev {
                if row <= p {
                    return Err(InputError::InvalidInput(
                        "row_idx not strictly increasing within column".to_string(),
                    ));
                }
            }
            prev = Some(row);
        }
    }
    Ok(())
}

fn read_u16(bytes: &[u8], offset: usize) -> u16 {
    let raw = bytes[offset..offset + 2]
        .as_array()
        .expect("slice length should be exactly 2");
    u16::from_le_bytes(*raw)
}

fn read_u32(bytes: &[u8], offset: usize) -> u32 {
    let raw = bytes[offset..offset + 4]
        .as_array()
        .expect("slice length should be exactly 4");
    u32::from_le_bytes(*raw)
}

fn read_u64(bytes: &[u8], offset: usize) -> u64 {
    let raw = bytes[offset..offset + 8]
        .as_array()
        .expect("slice length should be exactly 8");
    u64::from_le_bytes(*raw)
}

fn read_u32_vec(bytes: &[u8], offset: u64, len: usize) -> Result<Vec<u32>, InputError> {
    let offset = offset as usize;
    let byte_len = len
        .checked_mul(4)
        .ok_or_else(|| InputError::InvalidInput("u32 size overflow".to_string()))?;
    if offset + byte_len > bytes.len() {
        return Err(InputError::InvalidInput(
            "u32 vector out of bounds".to_string(),
        ));
    }
    Ok(bytes[offset..offset + byte_len]
        .chunks_exact(4)
        .map(|chunk| {
            u32::from_le_bytes(
                *chunk
                    .as_array()
                    .expect("chunk size from chunks_exact(4) must be 4"),
            )
        })
        .collect())
}

fn read_u64_vec(bytes: &[u8], offset: u64, len: usize) -> Result<Vec<u64>, InputError> {
    let offset = offset as usize;
    let byte_len = len
        .checked_mul(8)
        .ok_or_else(|| InputError::InvalidInput("u64 size overflow".to_string()))?;
    if offset + byte_len > bytes.len() {
        return Err(InputError::InvalidInput(
            "u64 vector out of bounds".to_string(),
        ));
    }
    Ok(bytes[offset..offset + byte_len]
        .chunks_exact(8)
        .map(|chunk| {
            u64::from_le_bytes(
                *chunk
                    .as_array()
                    .expect("chunk size from chunks_exact(8) must be 8"),
            )
        })
        .collect())
}

pub fn crc64_ecma(bytes: &[u8]) -> u64 {
    let mut crc = 0u64;
    for &b in bytes {
        crc ^= (b as u64) << 56;
        for _ in 0..8 {
            if (crc & 0x8000_0000_0000_0000) != 0 {
                crc = (crc << 1) ^ 0x42F0_E1EB_A9EA_3693;
            } else {
                crc <<= 1;
            }
        }
    }
    crc
}

#[cfg(test)]
#[path = "../../tests/src_inline/input/organelle_bin.rs"]
mod tests;
