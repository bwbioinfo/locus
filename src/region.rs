use crate::error::LocusError;

/// A genomic region. Coordinates are 0-based half-open [start, end).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub contig: String,
    /// 0-based inclusive start
    pub start: u64,
    /// 0-based exclusive end
    pub end: u64,
}

impl Region {
    pub fn new(contig: impl Into<String>, start: u64, end: u64) -> Self {
        Self {
            contig: contig.into(),
            start,
            end,
        }
    }

    pub fn len(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    /// Return a padded region extended by `pad` bases on each side, clamped to [0, contig_len).
    pub fn padded(&self, pad: u64, contig_len: u64) -> Self {
        Self {
            contig: self.contig.clone(),
            start: self.start.saturating_sub(pad),
            end: (self.end + pad).min(contig_len),
        }
    }
}

impl std::fmt::Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Display as 1-based closed interval
        write!(f, "{}:{}-{}", self.contig, self.start + 1, self.end)
    }
}

/// Parse a region string such as:
///   chr1
///   chr1:1000-2000
///   chr1:1,000-2,000
///   chr1:1000
pub fn parse_region(s: &str) -> Result<Region, LocusError> {
    // strip any whitespace
    let s = s.trim();

    let colon = s.find(':');
    if colon.is_none() {
        // whole-contig region; we don't know the length yet, use u64::MAX as sentinel
        return Ok(Region::new(s, 0, u64::MAX));
    }

    let colon = colon.unwrap();
    let contig = &s[..colon];
    let rest = &s[colon + 1..];

    // strip commas (1,000 -> 1000)
    let rest_clean = rest.replace(',', "");

    if let Some(dash) = rest_clean.find('-') {
        let start_str = &rest_clean[..dash];
        let end_str = &rest_clean[dash + 1..];

        let start_1based: u64 = start_str
            .parse()
            .map_err(|_| LocusError::MalformedRegion(s.to_string()))?;
        let end_1based: u64 = end_str
            .parse()
            .map_err(|_| LocusError::MalformedRegion(s.to_string()))?;

        if start_1based == 0 || end_1based < start_1based {
            return Err(LocusError::MalformedRegion(s.to_string()));
        }

        // convert 1-based closed [start, end] -> 0-based half-open [start-1, end)
        Ok(Region::new(contig, start_1based - 1, end_1based))
    } else {
        // single position
        let pos: u64 = rest_clean
            .parse()
            .map_err(|_| LocusError::MalformedRegion(s.to_string()))?;
        if pos == 0 {
            return Err(LocusError::MalformedRegion(s.to_string()));
        }
        let start = pos - 1;
        Ok(Region::new(contig, start, start + 1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_whole_contig() {
        let r = parse_region("chr1").unwrap();
        assert_eq!(r.contig, "chr1");
        assert_eq!(r.start, 0);
        assert_eq!(r.end, u64::MAX);
    }

    #[test]
    fn test_basic_range() {
        let r = parse_region("chr1:1000-2000").unwrap();
        assert_eq!(r.start, 999); // 0-based
        assert_eq!(r.end, 2000);
    }

    #[test]
    fn test_commas() {
        let r = parse_region("chr1:1,000-2,000").unwrap();
        assert_eq!(r.start, 999);
        assert_eq!(r.end, 2000);
    }

    #[test]
    fn test_single_pos() {
        let r = parse_region("chr1:500").unwrap();
        assert_eq!(r.start, 499);
        assert_eq!(r.end, 500);
    }

    #[test]
    fn test_bad_region() {
        assert!(parse_region("chr1:abc-def").is_err());
        assert!(parse_region("chr1:0-100").is_err());
    }

    #[test]
    fn test_display_is_1based() {
        let r = Region::new("chr1", 999, 2000);
        assert_eq!(r.to_string(), "chr1:1000-2000");
    }

    #[test]
    fn test_padded() {
        let r = Region::new("chr1", 1000, 2000);
        let p = r.padded(500, 100_000);
        assert_eq!(p.start, 500);
        assert_eq!(p.end, 2500);
    }
}
