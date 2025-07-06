#[cfg(feature = "cli")]
use clap::Parser;
#[cfg(feature = "cli")]
use seqrush::Args;

#[cfg(feature = "cli")]
#[derive(Parser)]
#[command(name = "seqrush", about = "Build pangenome graphs")] 
struct CliArgs {
    /// Input FASTA file
    #[arg(short = 's', long)]
    sequences: String,

    /// Output GFA file
    #[arg(short = 'o', long)]
    output: String,

    /// Number of worker threads
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Minimum match length
    #[arg(short = 'k', long = "min-match-length", default_value_t = 15)]
    min_match_length: usize,
}

#[cfg(feature = "cli")]
/// Parse command line arguments into `Args` using clap.
pub fn parse() -> Args {
    let cli = CliArgs::parse();
    Args {
        sequences: cli.sequences,
        output: cli.output,
        threads: cli.threads,
        min_match_length: cli.min_match_length,
    }
}

