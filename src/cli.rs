use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(name = "locus", about = "Terminal genome browser")]
pub struct Args {
    #[command(subcommand)]
    pub command: Option<Command>,

    /// BAM file to open (must have .bai index)
    pub bam: Option<String>,

    /// Jump to region on startup, e.g. chr1:100000-101000
    #[arg(long, short)]
    pub region: Option<String>,

    /// GFF3/GTF annotation file for feature track and gene search
    #[arg(long, short = 'a')]
    pub gff: Option<String>,

    /// Reference FASTA for the reference track and mismatch coloring
    #[arg(long, short = 'f')]
    pub reference: Option<String>,

    /// Start with the light color theme
    #[arg(long)]
    pub light: bool,

    /// Hide reads with mapping quality below this threshold
    #[arg(long, default_value_t = 0)]
    pub min_mapq: u8,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Sort, BGZF-compress, and tabix-index a GFF3/GTF annotation file.
    PrepareAnnotations {
        /// Input GFF3/GTF annotation file.
        input: String,

        /// Output BGZF-compressed annotation path.
        #[arg(long, short)]
        output: String,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_min_mapq_threshold() {
        let args = Args::try_parse_from(["locus", "sample.bam", "--min-mapq", "30"])
            .expect("parse valid minimum MAPQ");

        assert_eq!(args.min_mapq, 30);
    }

    #[test]
    fn defaults_min_mapq_to_zero() {
        let args = Args::try_parse_from(["locus", "sample.bam"]).expect("parse defaults");

        assert_eq!(args.min_mapq, 0);
    }

    #[test]
    fn rejects_min_mapq_above_u8_range() {
        assert!(Args::try_parse_from(["locus", "sample.bam", "--min-mapq", "256"]).is_err());
    }
}
