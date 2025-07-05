use seqrush::seqrush::{Args, run_seqrush};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut iter = std::env::args().skip(1);
    let sequences = iter.next().expect("input FASTA required");
    let output = iter.next().expect("output file required");
    let args = Args {
        sequences,
        output,
        threads: 1,
        min_match_length: 15,
    };
    run_seqrush(args)?;
    Ok(())
}