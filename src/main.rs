use clap::Parser;
use seqrush::seqrush::{run_seqrush, Args};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    run_seqrush(args)
}
