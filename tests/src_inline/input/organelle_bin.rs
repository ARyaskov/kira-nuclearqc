use super::*;
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

static DIR_COUNTER: AtomicUsize = AtomicUsize::new(0);

fn make_temp_dir() -> PathBuf {
    let mut dir = std::env::temp_dir();
    let id = DIR_COUNTER.fetch_add(1, Ordering::SeqCst);
    dir.push(format!(
        "kira_organelle_bin_test_{}_{}",
        std::process::id(),
        id
    ));
    fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn test_reader_minimal_bin() {
    let dir = make_temp_dir();
    let path = dir.join("kira-organelle.bin");
    let file = build_test_bin();
    fs::write(&path, file).unwrap();

    let bin = read_organelle_bin(&path).unwrap();
    assert_eq!(bin.genes, vec!["GENEA", "GENEB", "GENEC"]);
    assert_eq!(bin.barcodes, vec!["BC1", "BC2"]);
    assert_eq!(bin.csc.col_ptr, vec![0, 2, 3]);
    assert_eq!(bin.csc.row_idx, vec![0, 2, 1]);
    assert_eq!(bin.csc.values, vec![5, 1, 7]);
}

fn build_test_bin() -> Vec<u8> {
    let genes = build_string_table(&["GENEA", "GENEB", "GENEC"]);
    let barcodes = build_string_table(&["BC1", "BC2"]);
    let col_ptr = [0u64, 2, 3];
    let row_idx = [0u32, 2, 1];
    let values = [5u32, 1, 7];

    let mut offset = 256usize;
    let genes_offset = offset;
    offset = align64(offset + genes.len());
    let barcodes_offset = offset;
    offset = align64(offset + barcodes.len());
    let col_ptr_offset = offset;
    offset = align64(offset + col_ptr.len() * 8);
    let row_idx_offset = offset;
    offset = align64(offset + row_idx.len() * 4);
    let values_offset = offset;
    offset = align64(offset + values.len() * 4);

    let mut bytes = vec![0u8; offset];
    let file_bytes = bytes.len() as u64;
    bytes[0..4].copy_from_slice(b"KORG");
    bytes[4..6].copy_from_slice(&1u16.to_le_bytes());
    bytes[6..8].copy_from_slice(&0u16.to_le_bytes());
    bytes[8..12].copy_from_slice(&0x1234_5678u32.to_le_bytes());
    bytes[12..16].copy_from_slice(&256u32.to_le_bytes());
    bytes[16..24].copy_from_slice(&3u64.to_le_bytes());
    bytes[24..32].copy_from_slice(&2u64.to_le_bytes());
    bytes[32..40].copy_from_slice(&3u64.to_le_bytes());
    bytes[40..48].copy_from_slice(&(genes_offset as u64).to_le_bytes());
    bytes[48..56].copy_from_slice(&(genes.len() as u64).to_le_bytes());
    bytes[56..64].copy_from_slice(&(barcodes_offset as u64).to_le_bytes());
    bytes[64..72].copy_from_slice(&(barcodes.len() as u64).to_le_bytes());
    bytes[72..80].copy_from_slice(&(col_ptr_offset as u64).to_le_bytes());
    bytes[80..88].copy_from_slice(&(row_idx_offset as u64).to_le_bytes());
    bytes[88..96].copy_from_slice(&(values_offset as u64).to_le_bytes());
    bytes[96..104].copy_from_slice(&0u64.to_le_bytes());
    bytes[104..112].copy_from_slice(&0u64.to_le_bytes());
    bytes[112..120].copy_from_slice(&file_bytes.to_le_bytes());
    bytes[128..136].copy_from_slice(&0u64.to_le_bytes());

    bytes[genes_offset..genes_offset + genes.len()].copy_from_slice(&genes);
    bytes[barcodes_offset..barcodes_offset + barcodes.len()].copy_from_slice(&barcodes);
    for (i, v) in col_ptr.iter().enumerate() {
        let p = col_ptr_offset + i * 8;
        bytes[p..p + 8].copy_from_slice(&v.to_le_bytes());
    }
    for (i, v) in row_idx.iter().enumerate() {
        let p = row_idx_offset + i * 4;
        bytes[p..p + 4].copy_from_slice(&v.to_le_bytes());
    }
    for (i, v) in values.iter().enumerate() {
        let p = values_offset + i * 4;
        bytes[p..p + 4].copy_from_slice(&v.to_le_bytes());
    }

    let mut hdr = bytes[0..256].to_vec();
    hdr[120..128].fill(0);
    let crc = crc64_ecma(&hdr);
    bytes[120..128].copy_from_slice(&crc.to_le_bytes());
    bytes
}

fn build_string_table(items: &[&str]) -> Vec<u8> {
    let mut blob = Vec::new();
    let mut offsets = Vec::with_capacity(items.len() + 1);
    offsets.push(0u32);
    for item in items {
        blob.extend_from_slice(item.as_bytes());
        offsets.push(blob.len() as u32);
    }
    let mut out = Vec::new();
    out.extend_from_slice(&(items.len() as u32).to_le_bytes());
    for off in offsets {
        out.extend_from_slice(&off.to_le_bytes());
    }
    out.extend_from_slice(&blob);
    out
}

fn align64(v: usize) -> usize {
    let rem = v % 64;
    if rem == 0 { v } else { v + (64 - rem) }
}
