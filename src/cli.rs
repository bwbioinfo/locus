use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "locus", about = "Terminal genome browser")]
pub struct Args {
    /// BAM file to open (must have .bai index)
    pub bam: String,

    /// Jump to region on startup, e.g. chr1:100000-101000
    #[arg(long, short)]
    pub region: Option<String>,

    /// GFF3/GTF annotation file for feature track and gene search
    #[arg(long, short = 'a')]
    pub gff: Option<String>,

    /// Reference FASTA (optional, reserved for future use)
    #[arg(long, short = 'f')]
    pub reference: Option<String>,
}
