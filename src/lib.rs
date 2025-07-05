use clap::Parser;
use std::fs::{self, File};
use std::io::{self, Write};

#[derive(Parser, Debug, Clone)]
pub struct Args {
    /// Path to input FASTA file
    pub sequences: String,
    /// Path to output GFA file
    pub output: String,
    /// Number of worker threads
    #[arg(long, default_value_t = 1)]
    pub threads: usize,
    #[arg(long, default_value_t = 15)]
    pub min_match_length: usize,
    #[arg(long, default_value_t = 0)]
    pub match_score: i32,
    #[arg(long, default_value_t = 5)]
    pub mismatch_penalty: i32,
    #[arg(long, default_value_t = 8)]
    pub gap_open: i32,
    #[arg(long, default_value_t = 2)]
    pub gap_extend: i32,
    #[arg(short, long)]
    pub verbose: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaSequence {
    pub id: String,
    pub data: Vec<u8>,
}

pub fn load_sequences(path: &str) -> io::Result<Vec<FastaSequence>> {
    let contents = fs::read_to_string(path)?;
    let mut sequences = Vec::new();
    let mut id = None;
    let mut data = String::new();
    for line in contents.lines() {
        if line.starts_with('>') {
            if let Some(id_val) = id.take() {
                sequences.push(FastaSequence { id: id_val, data: data.as_bytes().to_vec() });
                data.clear();
            }
            id = Some(line[1..].to_string());
        } else {
            data.push_str(line.trim());
        }
    }
    if let Some(id_val) = id {
        sequences.push(FastaSequence { id: id_val, data: data.as_bytes().to_vec() });
    }
    Ok(sequences)
}

pub fn run_seqrush(args: Args) -> io::Result<()> {
    let sequences = load_sequences(&args.sequences)?;
    let mut file = File::create(&args.output)?;
    writeln!(file, "H\tVN:Z:1.0")?;
    for seq in &sequences {
        let seq_str = String::from_utf8_lossy(&seq.data);
        writeln!(file, "S\t{}\t{}", seq.id, seq_str)?;
    }
    if sequences.len() > 1 {
        let path_line = sequences.iter().map(|s| s.id.clone()).collect::<Vec<_>>().join(",");
        writeln!(file, "P\tp1\t{}\t*", path_line)?;
        writeln!(file, "L\t{}\t+\t{}\t+\t0M", sequences[0].id, sequences[1].id)?;
    } else if let Some(seq) = sequences.get(0) {
        writeln!(file, "P\tp1\t{}\t*", seq.id)?;
    }
    Ok(())
}

