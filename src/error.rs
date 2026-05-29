use thiserror::Error;

#[derive(Debug, Error)]
#[allow(dead_code)]
pub enum LocusError {
    #[error("BAM file not found: {0}")]
    BamNotFound(String),

    #[error("BAM index not found: {0}")]
    IndexNotFound(String),

    #[error("Unknown contig: {0}")]
    UnknownContig(String),

    #[error("Malformed region string: {0}")]
    MalformedRegion(String),

    #[error("BAM not indexed or not coordinate-sorted")]
    NotIndexed,

    #[error("Terminal too small (need at least 40x10)")]
    TerminalTooSmall,

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}
