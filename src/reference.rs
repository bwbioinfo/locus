use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;

use crate::region::Region;

#[derive(Debug, Clone)]
pub struct ReferenceSlice {
    pub start: u64,
    pub bases: Vec<u8>,
}

impl ReferenceSlice {
    pub fn base_at(&self, ref_pos: u64) -> Option<u8> {
        if ref_pos < self.start {
            return None;
        }
        self.bases.get((ref_pos - self.start) as usize).copied()
    }

    pub fn end(&self) -> u64 {
        self.start + self.bases.len() as u64
    }
}

#[derive(Debug)]
pub struct ReferenceStore {
    source: ReferenceSource,
}

#[derive(Debug)]
enum ReferenceSource {
    Indexed {
        path: PathBuf,
        records: HashMap<String, FaiRecord>,
    },
    InMemory {
        contigs: HashMap<String, Vec<u8>>,
    },
}

#[derive(Debug, Clone)]
struct FaiRecord {
    len: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

impl ReferenceStore {
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let fai_path = path.with_extension(format!(
            "{}fai",
            path.extension()
                .and_then(|ext| ext.to_str())
                .map(|ext| format!("{ext}."))
                .unwrap_or_default()
        ));

        if fai_path.exists() && !is_gzip_like(path) {
            let records = load_fai(&fai_path)
                .with_context(|| format!("loading FASTA index {}", fai_path.display()))?;
            return Ok(Self {
                source: ReferenceSource::Indexed {
                    path: path.to_path_buf(),
                    records,
                },
            });
        }

        let contigs = load_fasta_into_memory(path)
            .with_context(|| format!("loading reference FASTA {}", path.display()))?;
        Ok(Self {
            source: ReferenceSource::InMemory { contigs },
        })
    }

    pub fn fetch(&self, region: &Region) -> Result<Option<ReferenceSlice>> {
        match &self.source {
            ReferenceSource::Indexed { path, records } => {
                let Some(record) = records.get(&region.contig) else {
                    return Ok(None);
                };
                let start = region.start.min(record.len);
                let end = region.end.min(record.len);
                if start >= end {
                    return Ok(Some(ReferenceSlice {
                        start,
                        bases: Vec::new(),
                    }));
                }
                let bases = fetch_indexed(path, record, start, end)?;
                Ok(Some(ReferenceSlice { start, bases }))
            }
            ReferenceSource::InMemory { contigs } => {
                let Some(seq) = contigs.get(&region.contig) else {
                    return Ok(None);
                };
                let start = region.start.min(seq.len() as u64);
                let end = region.end.min(seq.len() as u64);
                Ok(Some(ReferenceSlice {
                    start,
                    bases: seq[start as usize..end as usize].to_vec(),
                }))
            }
        }
    }
}

fn load_fai(path: &Path) -> Result<HashMap<String, FaiRecord>> {
    let reader = BufReader::new(File::open(path)?);
    let mut records = HashMap::new();
    for line_result in reader.lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        let record = FaiRecord {
            len: fields[1].parse()?,
            offset: fields[2].parse()?,
            line_bases: fields[3].parse()?,
            line_width: fields[4].parse()?,
        };
        records.insert(fields[0].to_string(), record);
    }
    Ok(records)
}

fn fetch_indexed(path: &Path, record: &FaiRecord, start: u64, end: u64) -> io::Result<Vec<u8>> {
    let mut file = File::open(path)?;
    let mut bases = Vec::with_capacity((end - start) as usize);
    let byte_offset = record.offset
        + (start / record.line_bases) * record.line_width
        + (start % record.line_bases);
    file.seek(SeekFrom::Start(byte_offset))?;

    let mut buf = [0u8; 8192];
    while bases.len() < (end - start) as usize {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        for &byte in &buf[..n] {
            if byte == b'\n' || byte == b'\r' {
                continue;
            }
            bases.push(byte.to_ascii_uppercase());
            if bases.len() == (end - start) as usize {
                break;
            }
        }
    }
    Ok(bases)
}

fn load_fasta_into_memory(path: &Path) -> io::Result<HashMap<String, Vec<u8>>> {
    let mut reader = open_fasta_reader(path)?;
    let mut text = String::new();
    reader.read_to_string(&mut text)?;

    let mut contigs = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = Vec::new();

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = current_name.replace(header_name(rest)) {
                contigs.insert(name, std::mem::take(&mut current_seq));
            }
        } else {
            current_seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }

    if let Some(name) = current_name {
        contigs.insert(name, current_seq);
    }

    Ok(contigs)
}

fn open_fasta_reader(path: &Path) -> io::Result<Box<dyn Read>> {
    let file = File::open(path)?;
    if is_gzip_like(path) {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

fn header_name(header: &str) -> String {
    header
        .split_whitespace()
        .next()
        .unwrap_or(header)
        .to_string()
}

fn is_gzip_like(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| matches!(ext.to_ascii_lowercase().as_str(), "gz" | "bgz" | "bgzip"))
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loads_plain_fasta_regions() {
        let path = std::env::temp_dir().join(format!("locus-ref-test-{}.fa", std::process::id()));
        std::fs::write(&path, ">chr1 description\nACGT\nNNAA\n").unwrap();

        let reference = ReferenceStore::load(&path).unwrap();
        let slice = reference
            .fetch(&Region::new("chr1", 2, 7))
            .unwrap()
            .unwrap();

        assert_eq!(slice.start, 2);
        assert_eq!(slice.end(), 7);
        assert_eq!(slice.bases, b"GTNNA");

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn uses_fai_for_plain_fasta_regions() {
        let path =
            std::env::temp_dir().join(format!("locus-ref-index-test-{}.fa", std::process::id()));
        let fai_path = path.with_extension("fa.fai");
        std::fs::write(&path, ">chr1\nACGT\nNNAA\n").unwrap();
        std::fs::write(&fai_path, "chr1\t8\t6\t4\t5\n").unwrap();

        let reference = ReferenceStore::load(&path).unwrap();
        let slice = reference
            .fetch(&Region::new("chr1", 3, 8))
            .unwrap()
            .unwrap();

        assert_eq!(slice.bases, b"TNNAA");

        let _ = std::fs::remove_file(path);
        let _ = std::fs::remove_file(fai_path);
    }
}
