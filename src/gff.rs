use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, Read},
    path::Path,
};

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;

use crate::cache::Strand;
use crate::region::Region;

/// A single GFF3/GTF feature record (normalised to 0-based half-open coordinates).
#[derive(Debug, Clone)]
pub struct GffFeature {
    pub seqname: String,
    pub feature_type: String,
    /// 0-based inclusive start
    pub start: u64,
    /// 0-based exclusive end
    pub end: u64,
    pub strand: Option<Strand>,
    /// GFF3 ID= attribute or GTF gene_id/transcript_id fallback
    pub id: Option<String>,
    /// GFF3 Name= attribute or GTF gene_name fallback
    pub name: Option<String>,
    /// GFF3 Parent= attribute or GTF parent-like fallback
    #[allow(dead_code)]
    pub parent: Option<String>,
    /// gene_name= / gene= fallback
    pub gene_name: Option<String>,
}

impl GffFeature {
    /// Best display name: Name > gene_name > gene > ID > "<type>"
    pub fn display_name(&self) -> &str {
        self.name
            .as_deref()
            .or(self.gene_name.as_deref())
            .or(self.id.as_deref())
            .unwrap_or(&self.feature_type)
    }

    #[allow(dead_code)]
    pub fn to_region(&self) -> Region {
        Region::new(self.seqname.clone(), self.start, self.end)
    }
}

/// Holds all parsed GFF3/GTF features with lookup indices.
pub struct GffStore {
    pub features: Vec<GffFeature>,
    /// lowercase-name → feature indices for fast search
    name_index: HashMap<String, Vec<usize>>,
}

impl GffStore {
    /// Load a GFF3/GTF file. Plain text, gzip, and BGZF streams are supported.
    /// Comment lines and embedded FASTA sections are skipped.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let reader = open_annotation_reader(path)
            .with_context(|| format!("opening annotation file {}", path.display()))?;

        let mut features = Vec::new();
        let mut in_fasta = false;

        for line_result in reader.lines() {
            let line = line_result?;
            let trimmed = line.trim();

            if trimmed == "##FASTA" {
                in_fasta = true;
                continue;
            }
            if in_fasta || trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            if let Some(feat) = parse_feature_line(trimmed) {
                features.push(feat);
            }
        }

        let name_index = build_name_index(&features);
        Ok(Self {
            features,
            name_index,
        })
    }

    /// Return feature indices whose display name contains `query` (case-insensitive).
    /// Results are sorted: exact-match first, then prefix, then substring.
    pub fn search(&self, query: &str) -> Vec<usize> {
        let q = query.to_lowercase();
        if q.is_empty() {
            return Vec::new();
        }

        // Collect all indices where any searchable name contains q
        let mut exact: Vec<usize> = Vec::new();
        let mut prefix: Vec<usize> = Vec::new();
        let mut substr: Vec<usize> = Vec::new();
        let mut seen = std::collections::HashSet::new();

        for (key, indices) in &self.name_index {
            if key == &q {
                for &i in indices {
                    if seen.insert(i) {
                        exact.push(i);
                    }
                }
            } else if key.starts_with(&q as &str) {
                for &i in indices {
                    if seen.insert(i) {
                        prefix.push(i);
                    }
                }
            } else if key.contains(&q as &str) {
                for &i in indices {
                    if seen.insert(i) {
                        substr.push(i);
                    }
                }
            }
        }

        exact.extend(prefix);
        exact.extend(substr);
        exact
    }

    /// Return features overlapping [start, end) on the given contig.
    pub fn features_in_region<'a>(
        &'a self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> impl Iterator<Item = &'a GffFeature> {
        self.features
            .iter()
            .filter(move |f| f.seqname == contig && f.start < end && f.end > start)
    }
}

// ─── Parsing ─────────────────────────────────────────────────────────────────

fn open_annotation_reader(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    let reader: Box<dyn Read> = if is_gzip_like(path) {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    Ok(Box::new(BufReader::new(reader)))
}

fn is_gzip_like(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| matches!(ext.to_ascii_lowercase().as_str(), "gz" | "bgz" | "bgzip"))
        .unwrap_or(false)
}

fn parse_feature_line(line: &str) -> Option<GffFeature> {
    let fields: Vec<&str> = line.splitn(9, '\t').collect();
    if fields.len() < 8 {
        return None;
    }

    let seqname = fields[0].to_string();
    let feature_type = fields[2].to_string();

    // start / end are 1-based inclusive in GFF3
    let start_1: u64 = fields[3].parse().ok()?;
    let end_1: u64 = fields[4].parse().ok()?;
    if start_1 == 0 || end_1 < start_1 {
        return None;
    }
    let start = start_1 - 1; // 0-based
    let end = end_1; // 0-based exclusive

    let strand = match fields[6] {
        "+" => Some(Strand::Forward),
        "-" => Some(Strand::Reverse),
        _ => None,
    };

    let (id, name, parent, gene_name) = if fields.len() >= 9 {
        parse_attributes(fields[8])
    } else {
        (None, None, None, None)
    };

    Some(GffFeature {
        seqname,
        feature_type,
        start,
        end,
        strand,
        id,
        name,
        parent,
        gene_name,
    })
}

/// Extract ID, Name, Parent, and gene_name/gene from GFF3 or GTF attributes.
fn parse_attributes(
    attrs: &str,
) -> (
    Option<String>,
    Option<String>,
    Option<String>,
    Option<String>,
) {
    let mut id = None;
    let mut name = None;
    let mut parent = None;
    let mut gene_name = None;
    let mut gene_id = None;
    let mut transcript_id = None;

    for kv in attrs.split(';') {
        let kv = kv.trim();
        if kv.is_empty() {
            continue;
        }
        let Some((key, val)) = parse_attribute_pair(kv) else {
            continue;
        };

        // take only the first value if comma-separated
        let first_val = val.split(',').next().unwrap_or("").to_string();
        let first_val = if first_val.is_empty() {
            continue;
        } else {
            first_val
        };

        match key.as_str() {
            "ID" => id = Some(first_val),
            "Name" => name = Some(first_val),
            "Parent" => parent = Some(first_val),
            "gene_id" => {
                gene_id.get_or_insert(first_val);
            }
            "transcript_id" => {
                transcript_id.get_or_insert(first_val);
            }
            "gene_name" | "gene" => {
                gene_name.get_or_insert(first_val.clone());
                name.get_or_insert(first_val);
            }
            _ => {}
        }
    }

    if id.is_none() {
        id = transcript_id.clone().or_else(|| gene_id.clone());
    }
    if parent.is_none() {
        parent = gene_id.clone();
    }
    if gene_name.is_none() {
        gene_name = gene_id;
    }

    (id, name, parent, gene_name)
}

fn parse_attribute_pair(kv: &str) -> Option<(String, String)> {
    if let Some(eq) = kv.find('=') {
        let key = kv[..eq].trim();
        let val = clean_attribute_value(kv[eq + 1..].trim());
        return Some((key.to_string(), val));
    }

    let mut parts = kv.splitn(2, char::is_whitespace);
    let key = parts.next()?.trim();
    let val = parts.next()?.trim();
    if key.is_empty() || val.is_empty() {
        return None;
    }

    Some((key.to_string(), clean_attribute_value(val)))
}

fn clean_attribute_value(value: &str) -> String {
    let value = value.trim();
    let value = value
        .strip_prefix('"')
        .and_then(|v| v.strip_suffix('"'))
        .unwrap_or(value);
    percent_decode(value)
}

/// Minimal percent-decoding for the characters GFF3 commonly encodes.
fn percent_decode(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            if let (Some(h), Some(l)) = (hex_val(bytes[i + 1]), hex_val(bytes[i + 2])) {
                out.push(char::from(h << 4 | l));
                i += 3;
                continue;
            }
        }
        out.push(bytes[i] as char);
        i += 1;
    }
    out
}

fn hex_val(b: u8) -> Option<u8> {
    match b {
        b'0'..=b'9' => Some(b - b'0'),
        b'A'..=b'F' => Some(b - b'A' + 10),
        b'a'..=b'f' => Some(b - b'a' + 10),
        _ => None,
    }
}

/// Build a lowercase name → indices map covering ID, Name, gene_name.
fn build_name_index(features: &[GffFeature]) -> HashMap<String, Vec<usize>> {
    let mut map: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, feat) in features.iter().enumerate() {
        for key in [
            feat.id.as_deref(),
            feat.name.as_deref(),
            feat.gene_name.as_deref(),
        ]
        .into_iter()
        .flatten()
        {
            map.entry(key.to_lowercase()).or_default().push(i);
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use flate2::{Compression, write::GzEncoder};

    use super::*;

    #[test]
    fn test_parse_basic_line() {
        let line =
            "chr1\tEnsembl\tgene\t11869\t14409\t.\t+\t.\tID=ENSG001;Name=DDX11L1;gene_name=DDX11L1";
        let feat = parse_feature_line(line).unwrap();
        assert_eq!(feat.seqname, "chr1");
        assert_eq!(feat.feature_type, "gene");
        assert_eq!(feat.start, 11868); // 0-based
        assert_eq!(feat.end, 14409); // exclusive
        assert_eq!(feat.id.as_deref(), Some("ENSG001"));
        assert_eq!(feat.name.as_deref(), Some("DDX11L1"));
        assert_eq!(feat.strand, Some(Strand::Forward));
    }

    #[test]
    fn test_parse_parent() {
        let line = "chr1\t.\texon\t100\t200\t.\t-\t.\tParent=ENSG001;ID=exon1";
        let feat = parse_feature_line(line).unwrap();
        assert_eq!(feat.parent.as_deref(), Some("ENSG001"));
        assert_eq!(feat.strand, Some(Strand::Reverse));
    }

    #[test]
    fn test_skip_comment() {
        assert!(parse_feature_line("# comment").is_none());
        assert!(parse_feature_line("").is_none());
    }

    #[test]
    fn test_percent_decode() {
        assert_eq!(percent_decode("hello%20world"), "hello world");
        assert_eq!(percent_decode("no%25encoding"), "no%encoding");
    }

    #[test]
    fn test_parse_gtf_attributes() {
        let line = "chr13\tRefSeq\texon\t16451446\t16451550\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\"; gene_name \"TPTE2\";";
        let feat = parse_feature_line(line).unwrap();

        assert_eq!(feat.seqname, "chr13");
        assert_eq!(feat.feature_type, "exon");
        assert_eq!(feat.start, 16_451_445);
        assert_eq!(feat.end, 16_451_550);
        assert_eq!(feat.id.as_deref(), Some("TX1"));
        assert_eq!(feat.parent.as_deref(), Some("GENE1"));
        assert_eq!(feat.name.as_deref(), Some("TPTE2"));
        assert_eq!(feat.gene_name.as_deref(), Some("TPTE2"));
    }

    #[test]
    fn test_load_gzipped_gtf() {
        let path =
            std::env::temp_dir().join(format!("locus-gff-test-{}.gtf.gz", std::process::id()));
        let file = std::fs::File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        writeln!(
            encoder,
            "chr13\tRefSeq\tgene\t16451446\t18451446\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"TPTE2\";"
        )
        .unwrap();
        encoder.finish().unwrap();

        let store = GffStore::load(&path).unwrap();
        let _ = std::fs::remove_file(path);

        assert_eq!(store.features.len(), 1);
        assert_eq!(store.features[0].display_name(), "TPTE2");
        assert_eq!(store.search("tpte2").len(), 1);
    }

    #[test]
    fn test_search_ranking() {
        let features = vec![
            make_feature("chr1", "gene", 0, 100, Some("BRCA1"), Some("BRCA1")),
            make_feature("chr1", "gene", 0, 100, Some("BRCA2"), Some("BRCA2")),
            make_feature("chr1", "gene", 0, 100, Some("NOTBRCA"), Some("NOTBRCA")),
        ];
        let name_index = build_name_index(&features);
        let store = GffStore {
            features,
            name_index,
        };
        let results = store.search("brca1");
        // BRCA1 exact match should appear first
        assert!(!results.is_empty());
        assert_eq!(store.features[results[0]].name.as_deref(), Some("BRCA1"));
        // BRCA2 (prefix of "brca2", not exact) and NOTBRCA (substring) come after
        let names: Vec<_> = results
            .iter()
            .map(|&i| store.features[i].name.as_deref().unwrap())
            .collect();
        assert!(names.contains(&"BRCA1"));
    }

    fn make_feature(
        seq: &str,
        ty: &str,
        s: u64,
        e: u64,
        id: Option<&str>,
        name: Option<&str>,
    ) -> GffFeature {
        GffFeature {
            seqname: seq.to_string(),
            feature_type: ty.to_string(),
            start: s,
            end: e,
            strand: None,
            id: id.map(|s| s.to_string()),
            name: name.map(|s| s.to_string()),
            parent: None,
            gene_name: None,
        }
    }
}
