use clap::Parser;
use seqrush::seqrush::{Args, run_seqrush};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    run_seqrush(args)
}