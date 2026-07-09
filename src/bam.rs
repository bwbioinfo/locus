use std::{
    cell::RefCell,
    fs::File,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::record::data::field::{Tag, Value, value::Array},
};

use crate::cache::RenderRead;
use crate::error::LocusError;
use crate::methylation::parse_modified_bases;
use crate::region::Region;

pub struct ContigInfo {
    pub name: String,
    pub length: u64,
}

/// BAM 4-bit base encoding: index = nibble value, value = ASCII base.
const BAM_BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Decode a BAM 4-bit packed sequence into ASCII bytes.
/// `n_bases` comes from the read length (sum of CIGAR ops that consume the read).
fn decode_sequence(encoded: &[u8], n_bases: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(n_bases);
    for i in 0..n_bases {
        let Some(&byte) = encoded.get(i / 2) else {
            out.push(b'N');
            continue;
        };
        let nibble = if i % 2 == 0 {
            (byte >> 4) & 0xf
        } else {
            byte & 0xf
        };
        out.push(BAM_BASES[nibble as usize]);
    }
    out
}

pub struct BamSource {
    pub path: PathBuf,
    pub header: sam::Header,
    pub contigs: Vec<ContigInfo>,
    // Persistent reader: avoids re-opening + re-reading the index on every fetch.
    // RefCell gives interior mutability for &self fetch_reads.
    reader: RefCell<bam::io::IndexedReader<bgzf::Reader<File>>>,
}

impl BamSource {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();

        if !path.exists() {
            return Err(LocusError::BamNotFound(path.display().to_string()).into());
        }

        let mut reader = bam::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| {
                format!(
                    "failed to open BAM (check that {}.bai exists)",
                    path.display()
                )
            })?;

        let header = reader.read_header().context("failed to read BAM header")?;

        let contigs = header
            .reference_sequences()
            .iter()
            .map(|(name, map)| ContigInfo {
                name: name.to_string(),
                length: map.length().get() as u64,
            })
            .collect();

        Ok(Self {
            path: path.to_path_buf(),
            header,
            contigs,
            reader: RefCell::new(reader),
        })
    }

    pub fn contig_len(&self, name: &str) -> Option<u64> {
        self.contigs
            .iter()
            .find(|c| c.name == name)
            .map(|c| c.length)
    }

    pub fn resolve_region(&self, region: &Region) -> Result<Region> {
        let len = self
            .contig_len(&region.contig)
            .ok_or_else(|| LocusError::UnknownContig(region.contig.clone()))?;

        if region.end == u64::MAX {
            let end = len.min(region.start + 1_000);
            Ok(Region::new(region.contig.clone(), region.start, end))
        } else {
            Ok(region.clone())
        }
    }

    /// Fetch reads overlapping `region`, reusing the persistent reader (no file re-open).
    pub fn fetch_reads(&self, region: &Region) -> Result<Vec<RenderRead>> {
        let contig_len = self
            .contig_len(&region.contig)
            .ok_or_else(|| LocusError::UnknownContig(region.contig.clone()))?;

        let noodles_start = (region.start + 1).min(contig_len);
        let noodles_end = region.end.min(contig_len);

        let start_pos = Position::try_from(noodles_start as usize)
            .map_err(|_| LocusError::MalformedRegion(region.to_string()))?;
        let end_pos = Position::try_from(noodles_end as usize)
            .map_err(|_| LocusError::MalformedRegion(region.to_string()))?;

        let noodles_region = noodles_core::Region::new(region.contig.as_str(), start_pos..=end_pos);

        // borrow_mut lasts only until the query results are collected into `reads`
        let mut reader = self.reader.borrow_mut();
        let query = reader
            .query(&self.header, &noodles_region)
            .with_context(|| format!("querying region {}", region))?;

        let reads = query
            .filter_map(|r| r.ok())
            .filter_map(|rec| record_to_render(&rec))
            .collect();

        Ok(reads)
    }
}

fn record_to_render(record: &bam::Record) -> Option<RenderRead> {
    use crate::cache::{CigarOp, Strand};
    use noodles_sam::alignment::record::cigar::op::Kind;

    let flags = record.flags();
    if flags.is_unmapped() {
        return None;
    }

    let start_pos = record.alignment_start()?.ok()?;
    let start_0based = (start_pos.get() as u64).saturating_sub(1);

    // Walk CIGAR: compute ref span, read length, and collect ops.
    let mut cigar_ops: Vec<CigarOp> = Vec::new();
    let mut ref_span: u64 = 0;
    let mut read_len: usize = 0;

    for op_result in record.cigar().iter() {
        let op = op_result.ok()?;
        let len = op.len() as u64;
        let kind = op.kind();

        let render_op = match kind {
            Kind::Match | Kind::SequenceMatch => CigarOp::Match(len),
            Kind::SequenceMismatch => CigarOp::Mismatch(len),
            Kind::Insertion => CigarOp::Insertion(len),
            Kind::Deletion => CigarOp::Deletion(len),
            Kind::Skip => CigarOp::Skip(len),
            Kind::SoftClip => CigarOp::SoftClip(len),
            Kind::HardClip | Kind::Pad => continue,
        };

        if kind.consumes_reference() {
            ref_span += len;
        }
        if kind.consumes_read() {
            read_len += len as usize;
        }
        cigar_ops.push(render_op);
    }

    // Decode 4-bit packed sequence. as_ref() gives the raw encoded bytes.
    let sequence = decode_sequence(record.sequence().as_ref(), read_len);
    let methylation = parse_record_methylation(record, &sequence);

    let end_0based = start_0based + ref_span;
    let name = record
        .name()
        .map(|n| n.to_string())
        .unwrap_or_else(|| "*".to_string());
    let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
    let strand = if flags.is_reverse_complemented() {
        Strand::Reverse
    } else {
        Strand::Forward
    };

    Some(RenderRead {
        name,
        start: start_0based,
        end: end_0based,
        strand,
        mapq,
        cigar_ops,
        sequence,
        methylation,
        is_secondary: flags.is_secondary(),
        is_supplementary: flags.is_supplementary(),
        is_duplicate: flags.is_duplicate(),
    })
}

fn parse_record_methylation(
    record: &bam::Record,
    sequence: &[u8],
) -> Vec<crate::cache::ModifiedBaseCall> {
    let data = record.data();
    let mm = match data.get(&Tag::new(b'M', b'M')).and_then(Result::ok) {
        Some(Value::String(value)) => value,
        _ => return Vec::new(),
    };
    let Ok(mm) = std::str::from_utf8(mm.as_ref()) else {
        return Vec::new();
    };

    let ml = data
        .get(&Tag::new(b'M', b'L'))
        .and_then(Result::ok)
        .and_then(|value| match value {
            Value::Array(Array::UInt8(values)) => {
                Some(values.iter().filter_map(Result::ok).collect::<Vec<_>>())
            }
            _ => None,
        });

    parse_modified_bases(mm, ml.as_deref(), sequence)
}
