pub mod seqrush {
    use std::fs::{self, File};
    use std::io::{self, Write};

    #[derive(Debug, Clone)]
    pub struct Args {
        /// Path to input FASTA file
        pub sequences: String,
        /// Path to output GFA file
        pub output: String,
        /// Number of worker threads (unused in this stub)
        pub threads: usize,
        /// Minimum match length (unused)
        pub min_match_length: usize,
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
                    sequences.push(FastaSequence {
                        id: id_val,
                        data: data.as_bytes().to_vec(),
                    });
                    data.clear();
                }
                id = Some(line[1..].to_string());
            } else {
                data.push_str(line.trim());
            }
        }
        if let Some(id_val) = id {
            sequences.push(FastaSequence {
                id: id_val,
                data: data.as_bytes().to_vec(),
            });
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
            let path_line = sequences
                .iter()
                .map(|s| s.id.clone())
                .collect::<Vec<_>>()
                .join(",");
            writeln!(file, "P\tp1\t{}\t*", path_line)?;
            writeln!(file, "L\t{}\t+\t{}\t+\t0M", sequences[0].id, sequences[1].id)?;
        } else if let Some(seq) = sequences.get(0) {
            writeln!(file, "P\tp1\t{}\t*", seq.id)?;
        }
        Ok(())
    }
}

pub use seqrush::{load_sequences, run_seqrush, Args, FastaSequence};

