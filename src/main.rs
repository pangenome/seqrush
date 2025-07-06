use seqrush::{Args, run_seqrush};

#[cfg(feature = "cli")]
mod cli;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    #[cfg(feature = "cli")]
    let args = cli::parse();

    #[cfg(not(feature = "cli"))]
    let args = {
        let mut iter = std::env::args().skip(1);
        let sequences = iter.next().expect("input FASTA required");
        let output = iter.next().expect("output file required");
        Args {
            sequences,
            output,
            threads: 1,
            min_match_length: 15,
        }
    };

    run_seqrush(args)?;
    Ok(())
}

