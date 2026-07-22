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

use crate::cache::{ReadPhase, RenderRead};
use crate::error::LocusError;
use crate::methylation::parse_modified_bases;
use crate::region::Region;

pub struct ContigInfo {
    pub name: String,
    pub length: u64,
}

const DEFAULT_VIEW_SPAN: u64 = 1_000;

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

    /// Return a default viewport centered on the first mapped alignment in file order.
    pub fn first_mapped_region(&self) -> Result<Option<Region>> {
        let mut reader = self.reader.borrow_mut();

        for result in reader.records() {
            let record = result.context("failed to read BAM while locating its first alignment")?;
            if record.flags().is_unmapped() {
                continue;
            }

            let Some(reference_sequence_id) = record
                .reference_sequence_id()
                .transpose()
                .context("invalid reference sequence ID on first mapped alignment")?
            else {
                continue;
            };
            let Some(alignment_start) = record
                .alignment_start()
                .transpose()
                .context("invalid position on first mapped alignment")?
            else {
                continue;
            };
            let contig = self.contigs.get(reference_sequence_id).with_context(|| {
                format!(
                    "first mapped alignment references unknown sequence ID {reference_sequence_id}"
                )
            })?;
            let start_0based = (alignment_start.get() as u64).saturating_sub(1);

            return Ok(Some(initial_region_for_alignment(contig, start_0based)));
        }

        Ok(None)
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
    let phase = parse_record_phase(record);

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
        phase,
        is_secondary: flags.is_secondary(),
        is_supplementary: flags.is_supplementary(),
        is_duplicate: flags.is_duplicate(),
    })
}

fn parse_record_phase(record: &bam::Record) -> ReadPhase {
    let data = record.data();
    let haplotype =
        parse_unsigned_integer_tag(&data, Tag::new(b'H', b'P')).filter(|haplotype| *haplotype > 0);
    let phase_set = parse_unsigned_integer_tag(&data, Tag::new(b'P', b'S'));

    ReadPhase {
        haplotype,
        phase_set,
    }
}

fn parse_unsigned_integer_tag(data: &bam::record::Data<'_>, tag: Tag) -> Option<u32> {
    data.get(&tag)
        .and_then(Result::ok)
        .and_then(parse_unsigned_integer_value)
}

fn parse_unsigned_integer_value(value: Value<'_>) -> Option<u32> {
    value.as_int().and_then(|value| u32::try_from(value).ok())
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

fn initial_region_for_alignment(contig: &ContigInfo, alignment_start: u64) -> Region {
    let span = DEFAULT_VIEW_SPAN.min(contig.length);
    let start = alignment_start
        .saturating_sub(span / 2)
        .min(contig.length.saturating_sub(span));

    Region::new(contig.name.clone(), start, start + span)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase_tag_values_accept_unsigned_integer_encodings() {
        assert_eq!(parse_unsigned_integer_value(Value::Int8(1)), Some(1));
        assert_eq!(
            parse_unsigned_integer_value(Value::UInt32(123_456)),
            Some(123_456)
        );
    }

    #[test]
    fn phase_tag_values_reject_negative_and_non_integer_values() {
        assert_eq!(parse_unsigned_integer_value(Value::Int32(-1)), None);
        assert_eq!(parse_unsigned_integer_value(Value::Float(1.0)), None);
    }

    #[test]
    fn initial_region_centers_alignment_and_clamps_to_contig() {
        let contig = ContigInfo {
            name: "chr1".to_string(),
            length: 1_000_000,
        };

        assert_eq!(
            initial_region_for_alignment(&contig, 500_000),
            Region::new("chr1", 499_500, 500_500)
        );
        assert_eq!(
            initial_region_for_alignment(&contig, 100),
            Region::new("chr1", 0, 1_000)
        );
        assert_eq!(
            initial_region_for_alignment(&contig, 999_900),
            Region::new("chr1", 999_000, 1_000_000)
        );
    }

    #[test]
    fn first_mapped_region_uses_first_demo_alignment_and_keeps_queries_working() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/demo/demo.sorted.bam");
        let source = BamSource::open(path).expect("open demo BAM");

        let region = source
            .first_mapped_region()
            .expect("locate first mapped alignment")
            .expect("demo BAM has mapped alignments");

        assert_eq!(region, Region::new("chrDemo", 0, 154));

        let reads = source.fetch_reads(&region).expect("query initial region");
        assert_eq!(
            reads.first().map(|read| read.name.as_str()),
            Some("read_ins_meth")
        );
    }

    #[test]
    fn demo_reads_parse_tagged_untagged_and_malformed_phase_values() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/demo/demo.sorted.bam");
        let source = BamSource::open(path).expect("open demo BAM");
        let reads = source
            .fetch_reads(&Region::new("chrDemo", 0, 154))
            .expect("fetch demo reads");

        let phase_for = |name: &str| {
            reads
                .iter()
                .find(|read| read.name == name)
                .map(|read| read.phase)
                .unwrap_or_else(|| panic!("missing demo read {name}"))
        };

        assert_eq!(
            phase_for("read_ins_meth"),
            ReadPhase {
                haplotype: Some(1),
                phase_set: Some(50),
            }
        );
        assert_eq!(
            phase_for("read_del"),
            ReadPhase {
                haplotype: Some(2),
                phase_set: Some(50),
            }
        );
        assert_eq!(phase_for("read_reverse_meth"), ReadPhase::default());
        assert_eq!(phase_for("read_bad_phase"), ReadPhase::default());
    }
}
