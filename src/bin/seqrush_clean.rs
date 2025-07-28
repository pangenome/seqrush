use clap::Parser;
use seqrush::seqrush_clean::{Args, run_seqrush};

fn main() {
    let args = Args::parse();
    
    if let Err(e) = run_seqrush(args) {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}