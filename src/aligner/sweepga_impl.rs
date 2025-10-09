/// SweepGA aligner implementation
use super::{Aligner, AlignmentRecord, AlignmentSequence};
use anyhow::Result;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use fastga_rs::{Config, FastGA};
use tempfile::NamedTempFile;

pub struct SweepgaAligner {
    fastga: FastGA,
    verbose: bool,
    temp_dir: Option<PathBuf>,
}

impl SweepgaAligner {
    pub fn new(
        frequency: Option<usize>,
        threads: usize,
        verbose: bool,
    ) -> Result<Self, Box<dyn Error>> {
        let mut config_builder = Config::builder()
            .num_threads(threads)
            .verbose(verbose);

        if let Some(freq) = frequency {
            config_builder = config_builder.adaptive_seed_cutoff(freq);
        }

        let config = config_builder.build();
        let fastga = FastGA::new(config)
            .map_err(|e| format!("Failed to create FastGA: {}", e))?;

        Ok(SweepgaAligner {
            fastga,
            verbose,
            temp_dir: None,
        })
    }

    /// Write sequences to a temporary FASTA file
    fn write_sequences_to_fasta(
        &self,
        sequences: &[AlignmentSequence],
    ) -> Result<NamedTempFile, Box<dyn Error>> {
        let mut temp_fasta = if let Some(dir) = &self.temp_dir {
            NamedTempFile::new_in(dir)?
        } else {
            NamedTempFile::new()?
        };

        for seq in sequences {
            writeln!(temp_fasta, ">{}", seq.id)?;
            // Write sequence in 80-character lines (standard FASTA format)
            for chunk in seq.seq.chunks(80) {
                temp_fasta.write_all(chunk)?;
                writeln!(temp_fasta)?;
            }
        }

        temp_fasta.flush()?;

        if self.verbose {
            eprintln!(
                "[sweepga] Wrote {} sequences to {}",
                sequences.len(),
                temp_fasta.path().display()
            );
        }

        Ok(temp_fasta)
    }

    /// Parse PAF file into alignment records
    fn parse_paf(&self, paf_path: &Path) -> Result<Vec<AlignmentRecord>, Box<dyn Error>> {
        let file = File::open(paf_path)?;
        let reader = BufReader::new(file);
        let mut records = Vec::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 12 {
                continue;
            }

            // Extract CIGAR string from tags
            let mut cigar = String::new();
            for field in &fields[12..] {
                if let Some(cg) = field.strip_prefix("cg:Z:") {
                    cigar = cg.to_string();
                    break;
                }
            }

            let record = AlignmentRecord {
                query_name: fields[0].to_string(),
                query_len: fields[1].parse().unwrap_or(0),
                query_start: fields[2].parse().unwrap_or(0),
                query_end: fields[3].parse().unwrap_or(0),
                strand: fields[4].chars().next().unwrap_or('+'),
                target_name: fields[5].to_string(),
                target_len: fields[6].parse().unwrap_or(0),
                target_start: fields[7].parse().unwrap_or(0),
                target_end: fields[8].parse().unwrap_or(0),
                cigar,
            };

            records.push(record);
        }

        if self.verbose {
            eprintln!("[sweepga] Parsed {} alignment records", records.len());
        }

        Ok(records)
    }
}

impl Aligner for SweepgaAligner {
    fn align_sequences(
        &self,
        sequences: &[AlignmentSequence],
    ) -> Result<Vec<AlignmentRecord>, Box<dyn Error>> {
        if self.verbose {
            eprintln!("[sweepga] Aligning {} sequences", sequences.len());
        }

        // Write sequences to a persistent temp file (not auto-deleted)
        let temp_fasta = self.write_sequences_to_fasta(sequences)?;
        let fasta_path = temp_fasta.path().to_path_buf();

        // Keep the temp file alive by persisting it
        let (_file, fasta_path) = temp_fasta.keep()?;

        if self.verbose {
            eprintln!("[sweepga] FASTA written to: {}", fasta_path.display());
        }

        // Use fastga-rs to align and write directly to PAF
        if self.verbose {
            eprintln!("[sweepga] Running FastGA alignment...");
        }

        // Create temp PAF file
        let paf_temp = NamedTempFile::new()?;

        // Call FastGA align_to_file
        let num_alignments = self.fastga.align_to_file(&fasta_path, &fasta_path, paf_temp.path())
            .map_err(|e| anyhow::anyhow!("FastGA alignment failed: {}", e))?;

        if self.verbose {
            eprintln!("[sweepga] FastGA produced {} alignments", num_alignments);
        }

        // Parse PAF output
        let records = self.parse_paf(paf_temp.path())?;

        if self.verbose {
            eprintln!("[sweepga] Alignment complete: {} records", records.len());
        }

        // Clean up temp FASTA
        let _ = std::fs::remove_file(&fasta_path);
        let gdb_path = fasta_path.with_extension("1gdb");
        let _ = std::fs::remove_file(&gdb_path);
        let bps_path = format!(".{}.bps", fasta_path.file_name().unwrap().to_string_lossy());
        let _ = std::fs::remove_file(fasta_path.parent().unwrap().join(&bps_path));

        Ok(records)
    }
}
