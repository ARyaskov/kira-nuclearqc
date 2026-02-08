use std::fs::{self, File};
use std::io::{BufRead, BufReader, Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::input::InputError;

pub fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>, InputError> {
    if path.extension().is_some_and(|ext| ext == "gz") {
        let output = Command::new("gzip").arg("-cd").arg(path).output()?;
        if !output.status.success() {
            return Err(InputError::InvalidInput(format!(
                "gzip failed to decompress {}",
                path.display()
            )));
        }
        Ok(Box::new(BufReader::new(Cursor::new(output.stdout))))
    } else {
        let file = File::open(path)?;
        Ok(Box::new(BufReader::new(file)))
    }
}

pub fn hash_file(path: &Path) -> Result<u64, InputError> {
    let mut file = File::open(path)?;
    let mut buf = [0u8; 8192];
    let mut hasher = Fnv64::new();
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hasher.finish())
}

pub fn hash_bytes(data: &[u8]) -> u64 {
    let mut hasher = Fnv64::new();
    hasher.update(data);
    hasher.finish()
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
const CACHE_VERSION: u32 = 1;

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
    let mut file = File::create(path)?;
    file.write_all(CACHE_MAGIC)?;
    write_u32(&mut file, CACHE_VERSION)?;
    write_f32(&mut file, meta.scale)?;
    write_u8(&mut file, if meta.log1p { 1 } else { 0 })?;
    file.write_all(&[0u8; 3])?;
    write_u32(&mut file, meta.n_cells)?;
    write_u32(&mut file, meta.n_genes)?;
    write_u64(&mut file, meta.hash_mtx)?;
    write_u64(&mut file, meta.hash_features)?;
    write_u64(&mut file, meta.hash_barcodes)?;
    write_u64(&mut file, meta.hash_gene_index)?;

    for &lib in &data.libsizes {
        write_f32(&mut file, lib)?;
    }
    for &n in &data.nnz {
        write_u32(&mut file, n)?;
    }

    for col in &data.columns {
        for (gene_id, value) in col {
            write_u32(&mut file, *gene_id)?;
            write_f32(&mut file, *value)?;
        }
    }
    Ok(())
}

pub fn read_normalized_cache(
    path: &Path,
    meta: &CacheMeta,
) -> Result<Option<CachedNormalizedData>, InputError> {
    if !path.exists() {
        return Ok(None);
    }
    let mut file = File::open(path)?;
    let mut magic = [0u8; 8];
    file.read_exact(&mut magic)?;
    if &magic != CACHE_MAGIC {
        return Ok(None);
    }
    let version = read_u32(&mut file)?;
    if version != CACHE_VERSION {
        return Ok(None);
    }
    let scale = read_f32(&mut file)?;
    let log1p = read_u8(&mut file)? != 0;
    let mut _reserved = [0u8; 3];
    file.read_exact(&mut _reserved)?;
    let n_cells = read_u32(&mut file)?;
    let n_genes = read_u32(&mut file)?;
    let hash_mtx = read_u64(&mut file)?;
    let hash_features = read_u64(&mut file)?;
    let hash_barcodes = read_u64(&mut file)?;
    let hash_gene_index = read_u64(&mut file)?;

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
    for item in &mut libsizes {
        *item = read_f32(&mut file)?;
    }
    let mut nnz = vec![0u32; n_cells as usize];
    for item in &mut nnz {
        *item = read_u32(&mut file)?;
    }

    let mut columns = Vec::with_capacity(n_cells as usize);
    for &count in &nnz {
        let mut col = Vec::with_capacity(count as usize);
        for _ in 0..count {
            let gene_id = read_u32(&mut file)?;
            let value = read_f32(&mut file)?;
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

struct Fnv64 {
    hash: u64,
}

impl Fnv64 {
    fn new() -> Self {
        Self {
            hash: 0xcbf29ce484222325,
        }
    }

    fn update(&mut self, data: &[u8]) {
        let mut h = self.hash;
        for &b in data {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
        self.hash = h;
    }

    fn finish(&self) -> u64 {
        self.hash
    }
}
