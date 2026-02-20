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
    let shared = kira_shared_sc_cache::read_shared_cache_owned(path).map_err(map_err)?;

    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    if mmap.len() < 256 {
        return Err(InputError::InvalidInput(
            "kira-organelle.bin too small".to_string(),
        ));
    }

    let header = parse_header(&mmap[..256])?;
    let csc = CscView {
        n_genes: shared.n_genes as usize,
        n_cells: shared.n_cells as usize,
        nnz: shared.nnz as usize,
        col_ptr: shared.col_ptr,
        row_idx: shared.row_idx,
        values: shared.values_u32,
    };

    Ok(OrganelleBin {
        header,
        genes: shared.genes,
        barcodes: shared.barcodes,
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

#[inline]
fn read_u16(bytes: &[u8], offset: usize) -> u16 {
    let mut arr = [0u8; 2];
    arr.copy_from_slice(&bytes[offset..offset + 2]);
    u16::from_le_bytes(arr)
}

#[inline]
fn read_u64(bytes: &[u8], offset: usize) -> u64 {
    let mut arr = [0u8; 8];
    arr.copy_from_slice(&bytes[offset..offset + 8]);
    u64::from_le_bytes(arr)
}

fn map_err(err: kira_shared_sc_cache::SharedCacheError) -> InputError {
    match err {
        kira_shared_sc_cache::SharedCacheError::Io { source, .. } => InputError::Io(source),
        kira_shared_sc_cache::SharedCacheError::Format { message, .. } => {
            InputError::InvalidInput(message)
        }
    }
}

fn crc64_ecma(bytes: &[u8]) -> u64 {
    kira_shared_sc_cache::crc64_ecma(bytes)
}

#[cfg(test)]
#[path = "../../tests/src_inline/input/organelle_bin.rs"]
mod tests;
