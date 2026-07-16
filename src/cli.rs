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
