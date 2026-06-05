use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;

use crate::input::InputError;

pub fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>, InputError> {
    let file = File::open(path)?;
    if path.extension().is_some_and(|ext| ext == "gz") {
        // Buffer the raw file so GzDecoder gets larger reads.
        let buffered_file = BufReader::with_capacity(64 * 1024, file);
        Ok(Box::new(BufReader::new(GzDecoder::new(buffered_file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// CRC64-ECMA over a file's contents, streamed in 64 KiB blocks.
pub fn hash_file(path: &Path) -> Result<u64, InputError> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(64 * 1024, file);
    let mut crc: u64 = 0;
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = reader.read(&mut buf)?;
        if n == 0 {
            break;
        }
        crc = crc64_ecma_continue(crc, &buf[..n]);
    }
    Ok(crc)
}

pub fn hash_bytes(data: &[u8]) -> u64 {
    kira_shared_sc_cache::crc64_ecma(data)
}

#[derive(Debug, Clone)]
pub struct CacheMeta {
    pub n_cells: u32,
    pub n_genes: u32,
    pub hash_mtx: u64,
    pub hash_features: u64,
    pub hash_barcodes: u64,
    pub hash_gene_index: u64,
    pub scale: f32,
    pub log1p: bool,
}

#[derive(Debug, Clone)]
pub struct CachedNormalizedData {
    pub libsizes: Vec<f32>,
    pub nnz: Vec<u32>,
    pub columns: Vec<Vec<(u32, f32)>>,
}

const CACHE_MAGIC: &[u8; 8] = b"KIRAQC2\0";
const CACHE_VERSION: u32 = 2;

pub fn cache_path_default(mtx_path: &Path) -> PathBuf {
    let dir = mtx_path.parent().unwrap_or_else(|| Path::new("."));
    dir.join("kira_nuclearqc.normcache")
}

pub fn write_normalized_cache(
    path: &Path,
    meta: &CacheMeta,
    data: &CachedNormalizedData,
) -> Result<(), InputError> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);

    w.write_all(CACHE_MAGIC)?;
    write_u32(&mut w, CACHE_VERSION)?;
    write_f32(&mut w, meta.scale)?;
    write_u8(&mut w, if meta.log1p { 1 } else { 0 })?;
    w.write_all(&[0u8; 3])?;
    write_u32(&mut w, meta.n_cells)?;
    write_u32(&mut w, meta.n_genes)?;
    write_u64(&mut w, meta.hash_mtx)?;
    write_u64(&mut w, meta.hash_features)?;
    write_u64(&mut w, meta.hash_barcodes)?;
    write_u64(&mut w, meta.hash_gene_index)?;

    // Pack libsizes and nnz as contiguous LE blobs in one write_all each.
    w.write_all(slice_as_bytes_f32(&data.libsizes))?;
    w.write_all(slice_as_bytes_u32(&data.nnz))?;

    // Columns: pack each column's (u32, f32) pairs into a scratch buffer.
    let mut scratch: Vec<u8> = Vec::new();
    for col in &data.columns {
        if col.is_empty() {
            continue;
        }
        scratch.clear();
        scratch.reserve(col.len() * 8);
        for &(gene_id, value) in col {
            scratch.extend_from_slice(&gene_id.to_le_bytes());
            scratch.extend_from_slice(&value.to_le_bytes());
        }
        w.write_all(&scratch)?;
    }

    w.flush()?;
    Ok(())
}

pub fn read_normalized_cache(
    path: &Path,
    meta: &CacheMeta,
) -> Result<Option<CachedNormalizedData>, InputError> {
    if !path.exists() {
        return Ok(None);
    }
    let file = File::open(path)?;
    let mut r = BufReader::with_capacity(1 << 20, file);

    let mut magic = [0u8; 8];
    r.read_exact(&mut magic)?;
    if &magic != CACHE_MAGIC {
        return Ok(None);
    }
    let version = read_u32(&mut r)?;
    if version != CACHE_VERSION {
        return Ok(None);
    }
    let scale = read_f32(&mut r)?;
    let log1p = read_u8(&mut r)? != 0;
    let mut _reserved = [0u8; 3];
    r.read_exact(&mut _reserved)?;
    let n_cells = read_u32(&mut r)?;
    let n_genes = read_u32(&mut r)?;
    let hash_mtx = read_u64(&mut r)?;
    let hash_features = read_u64(&mut r)?;
    let hash_barcodes = read_u64(&mut r)?;
    let hash_gene_index = read_u64(&mut r)?;

    if scale != meta.scale
        || log1p != meta.log1p
        || n_cells != meta.n_cells
        || n_genes != meta.n_genes
        || hash_mtx != meta.hash_mtx
        || hash_features != meta.hash_features
        || hash_barcodes != meta.hash_barcodes
        || hash_gene_index != meta.hash_gene_index
    {
        return Ok(None);
    }

    let mut libsizes = vec![0f32; n_cells as usize];
    r.read_exact(slice_as_bytes_mut_f32(&mut libsizes))?;
    let mut nnz = vec![0u32; n_cells as usize];
    r.read_exact(slice_as_bytes_mut_u32(&mut nnz))?;
    // On-disk format is LE; the byte views above/below assume a LE host.
    #[cfg(not(target_endian = "little"))]
    compile_error!("kira-nuclearqc requires a little-endian host");

    let mut columns = Vec::with_capacity(n_cells as usize);
    let mut byte_scratch: Vec<u8> = Vec::new();
    for &count in &nnz {
        if count == 0 {
            columns.push(Vec::new());
            continue;
        }
        let bytes = count as usize * 8;
        byte_scratch.resize(bytes, 0);
        r.read_exact(&mut byte_scratch)?;
        let mut col = Vec::with_capacity(count as usize);
        for chunk in byte_scratch.chunks_exact(8) {
            let gene_id = u32::from_le_bytes(chunk[..4].try_into().unwrap());
            let value = f32::from_le_bytes(chunk[4..].try_into().unwrap());
            col.push((gene_id, value));
        }
        columns.push(col);
    }

    Ok(Some(CachedNormalizedData {
        libsizes,
        nnz,
        columns,
    }))
}

fn write_u8<W: Write>(w: &mut W, v: u8) -> Result<(), InputError> {
    w.write_all(&[v])?;
    Ok(())
}

fn write_u32<W: Write>(w: &mut W, v: u32) -> Result<(), InputError> {
    w.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_u64<W: Write>(w: &mut W, v: u64) -> Result<(), InputError> {
    w.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn write_f32<W: Write>(w: &mut W, v: f32) -> Result<(), InputError> {
    w.write_all(&v.to_le_bytes())?;
    Ok(())
}

fn read_u8<R: Read>(r: &mut R) -> Result<u8, InputError> {
    let mut buf = [0u8; 1];
    r.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_u32<R: Read>(r: &mut R) -> Result<u32, InputError> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64<R: Read>(r: &mut R) -> Result<u64, InputError> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_f32<R: Read>(r: &mut R) -> Result<f32, InputError> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(f32::from_le_bytes(buf))
}

// Bytewise views over Vec<f32>/Vec<u32>; LE hosts only (see compile_error above).
fn slice_as_bytes_f32(slice: &[f32]) -> &[u8] {
    // SAFETY: f32 is POD and the host is LE, so bytes match the on-disk layout.
    unsafe {
        std::slice::from_raw_parts(slice.as_ptr() as *const u8, std::mem::size_of_val(slice))
    }
}

fn slice_as_bytes_u32(slice: &[u32]) -> &[u8] {
    // SAFETY: see slice_as_bytes_f32.
    unsafe {
        std::slice::from_raw_parts(slice.as_ptr() as *const u8, std::mem::size_of_val(slice))
    }
}

fn slice_as_bytes_mut_f32(slice: &mut [f32]) -> &mut [u8] {
    // SAFETY: see slice_as_bytes_f32; mutable variant.
    unsafe {
        std::slice::from_raw_parts_mut(
            slice.as_mut_ptr() as *mut u8,
            std::mem::size_of_val(slice),
        )
    }
}

fn slice_as_bytes_mut_u32(slice: &mut [u32]) -> &mut [u8] {
    // SAFETY: see slice_as_bytes_f32; mutable variant.
    unsafe {
        std::slice::from_raw_parts_mut(
            slice.as_mut_ptr() as *mut u8,
            std::mem::size_of_val(slice),
        )
    }
}

// Resumable CRC64-ECMA, same polynomial as kira_shared_sc_cache.
const CRC64_ECMA_POLY: u64 = 0x42F0_E1EB_A9EA_3693;

static CRC_TABLE: std::sync::OnceLock<[u64; 256]> = std::sync::OnceLock::new();

fn crc_table() -> &'static [u64; 256] {
    CRC_TABLE.get_or_init(|| {
        let mut table = [0u64; 256];
        let mut byte = 0usize;
        while byte < 256 {
            let mut crc = (byte as u64) << 56;
            let mut bit = 0;
            while bit < 8 {
                if crc & 0x8000_0000_0000_0000 != 0 {
                    crc = (crc << 1) ^ CRC64_ECMA_POLY;
                } else {
                    crc <<= 1;
                }
                bit += 1;
            }
            table[byte] = crc;
            byte += 1;
        }
        table
    })
}

fn crc64_ecma_continue(mut crc: u64, bytes: &[u8]) -> u64 {
    let table = crc_table();
    for &b in bytes {
        let idx = ((crc >> 56) as u8) ^ b;
        crc = (crc << 8) ^ table[idx as usize];
    }
    crc
}
