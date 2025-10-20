/// SweepGA aligner implementation using modern sweepga library API
use super::{Aligner, AlignmentRecord, AlignmentSequence};
use std::error::Error;
use std::io::Write;
use std::path::Path;
use sweepga::fastga_integration::FastGAIntegration;
use sweepga::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};
use tempfile::NamedTempFile;

pub struct SweepgaAligner {
    frequency: Option<usize>,
    threads: usize,
    verbose: bool,
}

impl SweepgaAligner {
    pub fn new(
        frequency: Option<usize>,
        threads: usize,
        verbose: bool,
    ) -> Result<Self, Box<dyn Error>> {
        Ok(SweepgaAligner {
            frequency,
            threads,
            verbose,
        })
    }

    /// Write sequences to a temporary FASTA file
    fn write_sequences_to_fasta(
        &self,
        sequences: &[AlignmentSequence],
    ) -> Result<NamedTempFile, Box<dyn Error>> {
        let mut temp_fasta = NamedTempFile::new()?;

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
                "[sweepga] Wrote {} sequences to temporary FASTA",
                sequences.len()
            );
        }

        Ok(temp_fasta)
    }

    /// Parse PAF file into alignment records
    fn parse_paf(&self, paf_path: &Path) -> Result<Vec<AlignmentRecord>, Box<dyn Error>> {
        use std::fs::File;
        use std::io::{BufRead, BufReader};

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
            eprintln!("[sweepga] Parsed {} alignment records from PAF", records.len());
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
            eprintln!("[sweepga] Starting alignment of {} sequences", sequences.len());
        }

        // Step 1: Write sequences to temp FASTA
        let temp_fasta = self.write_sequences_to_fasta(sequences)?;
        let fasta_path = temp_fasta.path();

        if self.verbose {
            eprintln!("[sweepga] Running FastGA alignment (frequency: {:?}, threads: {})",
                self.frequency, self.threads);
        }

        // Step 2: Run FastGA alignment using sweepga library
        let fastga = FastGAIntegration::new(self.frequency, self.threads);
        let temp_paf = fastga.align_to_temp_paf(fasta_path, fasta_path)
            .map_err(|e| format!("FastGA alignment failed: {}", e))?;

        if self.verbose {
            // Count lines in PAF
            use std::fs;
            let paf_content = fs::read_to_string(temp_paf.path())?;
            let line_count = paf_content.lines()
                .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
                .count();
            eprintln!("[sweepga] FastGA produced {} raw alignments", line_count);
        }

        // Step 3: Apply plane sweep filtering with 1:1 mode for graph construction
        // This keeps one best mapping per query-target pair while allowing
        // multiple mappings per sequence (needed for graph construction)
        if self.verbose {
            eprintln!("[sweepga] Applying plane sweep filtering (1:1 mode)...");
        }

        let filter_config = FilterConfig {
            mapping_filter_mode: FilterMode::OneToOne,
            mapping_max_per_query: Some(1),
            mapping_max_per_target: Some(1),
            scoring_function: ScoringFunction::LogLengthIdentity,
            overlap_threshold: 0.95,
            min_block_length: 0,  // Keep all alignments (no length filter for graph construction)
            scaffold_filter_mode: FilterMode::ManyToMany,
            scaffold_max_per_query: None,
            scaffold_max_per_target: None,
            plane_sweep_secondaries: 0,
            sparsity: 1.0,
            no_merge: true,
            chain_gap: 2000,
            scaffold_gap: 10000,
            min_scaffold_length: 0,
            scaffold_overlap_threshold: 0.95,
            scaffold_max_deviation: 0,
            prefix_delimiter: '#',
            skip_prefix: false,
            min_identity: 0.0,
            min_scaffold_identity: 0.0,
        };

        let temp_filtered = NamedTempFile::new()?;
        let filter = PafFilter::new(filter_config);
        filter.filter_paf(
            temp_paf.path().to_str().unwrap(),
            temp_filtered.path().to_str().unwrap()
        ).map_err(|e| format!("PAF filtering failed: {}", e))?;

        if self.verbose {
            // Count filtered lines
            use std::fs;
            let filtered_content = fs::read_to_string(temp_filtered.path())?;
            let filtered_count = filtered_content.lines()
                .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
                .count();
            eprintln!("[sweepga] After 1:1 filtering: {} alignments", filtered_count);
        }

        // Step 4: Parse filtered PAF
        let records = self.parse_paf(temp_filtered.path())?;

        if self.verbose {
            eprintln!("[sweepga] Alignment complete: {} final records", records.len());
        }

        Ok(records)
    }
}
