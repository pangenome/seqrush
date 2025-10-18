use crate::bidirected_ops::BidirectedGraph;
use crate::bidirected_union_find::BidirectedUnionFind;
use crate::graph_ops::{Edge, Graph, Node};
use crate::pos::{make_pos, Pos};
#[cfg(feature = "use-allwave")]
use allwave::{AlignmentParams, AllPairIterator, SparsificationStrategy};
use clap::Parser;
use rayon::prelude::*;
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::Arc;
use std::sync::Mutex;

#[derive(Parser)]
#[command(
    name = "seqrush",
    version = "0.4.0",
    about = "Dynamic pangenome graph construction"
)]
pub struct Args {
    /// Input FASTA file
    #[arg(short, long)]
    pub sequences: String,

    /// Input PAF file (optional, if not provided will compute alignments internally)
    #[arg(short = 'p', long)]
    pub paf: Option<String>,

    /// Output GFA file
    #[arg(short, long, default_value = "output.gfa")]
    pub output: String,

    /// Number of threads
    #[arg(short, long, default_value = "4")]
    pub threads: usize,

    /// Minimum match length
    #[arg(short = 'k', long, default_value = "0")]
    pub min_match_length: usize,

    /// Alignment scores (match,mismatch,gap_open,gap_extend[,gap2_open,gap2_extend])
    #[arg(short = 'S', long = "scores", default_value = "0,5,8,2,24,1")]
    pub scores: String,

    /// Orientation check alignment scores for fast edit distance (match,mismatch,gap_open,gap_extend)
    #[arg(long = "orientation-scores", default_value = "0,1,1,1")]
    pub orientation_scores: String,

    /// Maximum divergence threshold (0.0-1.0, e.g., 0.1 = 10% divergence)
    #[arg(short = 'd', long = "max-divergence")]
    pub max_divergence: Option<f64>,

    /// Verbose output
    #[arg(short, long)]
    pub verbose: bool,

    /// Test mode - disable bidirectional alignment
    #[arg(long, hide = true)]
    pub test_mode: bool,

    /// Disable node compaction (compaction merges linear chains of nodes)
    #[arg(long = "no-compact", default_value = "false")]
    pub no_compact: bool,

    /// Sparsification factor (keep this fraction of alignment pairs, 1.0 = keep all, 'auto' for automatic)
    #[arg(short = 'x', long = "sparsify", default_value = "1.0")]
    pub sparsification: String,

    /// Output alignments to PAF file
    #[arg(long = "output-alignments")]
    pub output_alignments: Option<String>,

    /// Validate PAF records as they are generated (helps catch bugs immediately)
    #[arg(long = "validate-paf", default_value = "true")]
    pub validate_paf: bool,

    /// Use seqwish-style builder (experimental toggle used in tests)
    #[arg(long = "seqwish-style", default_value = "false", hide = true)]
    pub seqwish_style: bool,

    /// Skip all sorting (by default, uses full Ygs pipeline: Y=SGD, g=groom, s=topological sort)
    #[arg(long = "no-sort", default_value = "false")]
    pub no_sort: bool,

    /// Skip the path-guided SGD phase (Y) of the Ygs pipeline
    #[arg(long = "skip-sgd", default_value = "false")]
    pub skip_sgd: bool,

    /// Skip the grooming phase (g) of the Ygs pipeline
    #[arg(long = "skip-groom", default_value = "false")]
    pub skip_groom: bool,

    /// Skip the topological sort phase (s) of the Ygs pipeline
    #[arg(long = "skip-topo", default_value = "false")]
    pub skip_topo: bool,

    /// Number of SGD iterations (default: 100, matching ODGI)
    #[arg(long = "sgd-iter-max", default_value = "100")]
    pub sgd_iter_max: u64,

    /// SGD learning rate parameter eta_max (default: calculated from graph)
    #[arg(long = "sgd-eta-max")]
    pub sgd_eta_max: Option<f64>,

    /// SGD cooling/momentum parameter theta (default: 0.99)
    #[arg(long = "sgd-theta", default_value = "0.99")]
    pub sgd_theta: f64,

    /// SGD convergence parameter eps (default: 0.01)
    #[arg(long = "sgd-eps", default_value = "0.01")]
    pub sgd_eps: f64,

    /// SGD cooling start (default: 0.5)
    #[arg(long = "sgd-cooling-start", default_value = "0.5")]
    pub sgd_cooling_start: f64,

    /// [DEPRECATED] Use sort-groom-sort strategy (use --skip-* flags instead)
    #[arg(long = "sort-groom-sort", default_value = "false", hide = true)]
    pub sort_groom_sort: bool,

    /// [DEPRECATED] Apply iterative grooming (use Ygs pipeline instead)
    #[arg(long = "iterative-groom", value_name = "N", hide = true)]
    pub iterative_groom: Option<usize>,

    /// [DEPRECATED] Apply ODGI-style grooming (use Ygs pipeline instead)
    #[arg(long = "odgi-groom", default_value = "false", hide = true)]
    pub odgi_style_groom: bool,

    /// [DEPRECATED] Old flag for SGD sorting (now enabled by default, use --skip-sgd to disable)
    #[arg(long = "sgd-sort", default_value = "false", hide = true)]
    pub sgd_sort: bool,

    /// [DEPRECATED] Use --skip-groom instead
    #[arg(long = "groom", default_value = "false", hide = true)]
    pub groom: bool,

    /// Aligner to use: 'allwave', 'sweepga' (default: allwave)
    #[arg(long = "aligner", default_value = "allwave")]
    pub aligner: String,

    /// K-mer frequency threshold for SweepGA (use k-mers occurring ≤ N times)
    #[arg(long = "frequency", short = 'f')]
    pub frequency: Option<usize>,
}

#[derive(Clone, Debug)]
pub struct AlignmentScores {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap1_open: i32,
    pub gap1_extend: i32,
    pub gap2_open: Option<i32>,
    pub gap2_extend: Option<i32>,
}

impl AlignmentScores {
    pub fn parse(scores_str: &str) -> Result<Self, String> {
        let parts: Vec<&str> = scores_str.split(',').collect();

        if parts.len() < 4 {
            return Err(
                "Scores must have at least 4 values: match,mismatch,gap1_open,gap1_extend"
                    .to_string(),
            );
        }

        if parts.len() > 6 {
            return Err("Too many score values provided (max 6)".to_string());
        }

        let match_score = parts[0]
            .parse::<i32>()
            .map_err(|_| format!("Invalid match score: {}", parts[0]))?;
        let mismatch_penalty = parts[1]
            .parse::<i32>()
            .map_err(|_| format!("Invalid mismatch penalty: {}", parts[1]))?;
        let gap1_open = parts[2]
            .parse::<i32>()
            .map_err(|_| format!("Invalid gap1_open penalty: {}", parts[2]))?;
        let gap1_extend = parts[3]
            .parse::<i32>()
            .map_err(|_| format!("Invalid gap1_extend penalty: {}", parts[3]))?;

        let (gap2_open, gap2_extend) = if parts.len() >= 6 {
            (
                Some(
                    parts[4]
                        .parse::<i32>()
                        .map_err(|_| format!("Invalid gap2_open penalty: {}", parts[4]))?,
                ),
                Some(
                    parts[5]
                        .parse::<i32>()
                        .map_err(|_| format!("Invalid gap2_extend penalty: {}", parts[5]))?,
                ),
            )
        } else {
            (None, None)
        };

        Ok(AlignmentScores {
            match_score,
            mismatch_penalty,
            gap1_open,
            gap1_extend,
            gap2_open,
            gap2_extend,
        })
    }

    pub fn parse_orientation(scores_str: &str) -> Result<Self, String> {
        let parts: Vec<&str> = scores_str.split(',').collect();

        if parts.len() != 4 {
            return Err(
                "Orientation scores must have exactly 4 values: match,mismatch,gap_open,gap_extend"
                    .to_string(),
            );
        }

        let match_score = parts[0]
            .parse::<i32>()
            .map_err(|_| format!("Invalid match score: {}", parts[0]))?;
        let mismatch_penalty = parts[1]
            .parse::<i32>()
            .map_err(|_| format!("Invalid mismatch penalty: {}", parts[1]))?;
        let gap_open = parts[2]
            .parse::<i32>()
            .map_err(|_| format!("Invalid gap_open penalty: {}", parts[2]))?;
        let gap_extend = parts[3]
            .parse::<i32>()
            .map_err(|_| format!("Invalid gap_extend penalty: {}", parts[3]))?;

        Ok(AlignmentScores {
            match_score,
            mismatch_penalty,
            gap1_open: gap_open,
            gap1_extend: gap_extend,
            gap2_open: None,
            gap2_extend: None,
        })
    }

    /// Calculate maximum acceptable score for a given sequence length and divergence
    pub fn max_score_for_divergence(&self, seq_len: usize, max_divergence: f64) -> i32 {
        // Estimate: worst case is all mismatches up to divergence threshold
        let max_mismatches = (seq_len as f64 * max_divergence).ceil() as i32;
        let max_gaps = (seq_len as f64 * max_divergence * 0.5).ceil() as i32; // Assume gaps are half of divergence

        // Score calculation: matches give 0, mismatches and gaps add penalties
        let mismatch_score = max_mismatches * self.mismatch_penalty;
        let gap_score = if max_gaps > 0 {
            self.gap1_open + (max_gaps - 1) * self.gap1_extend
        } else {
            0
        };

        // Ensure minimum threshold is reasonable
        let threshold = mismatch_score + gap_score;
        threshold.max(self.mismatch_penalty * 2) // At least 2 mismatches allowed
    }
}

#[derive(Clone)]
pub struct Sequence {
    pub id: String,
    pub data: Vec<u8>,
    pub offset: usize, // Offset in concatenated sequence space
}

impl Sequence {
    /// Get reverse complement of the sequence
    pub fn reverse_complement(&self) -> Vec<u8> {
        self.data
            .iter()
            .rev()
            .map(|&base| {
                match base {
                    b'A' | b'a' => b'T',
                    b'T' | b't' => b'A',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    _ => base, // Keep N and other ambiguous bases as-is
                }
            })
            .collect()
    }
}

pub struct SeqRush {
    pub sequences: Vec<Sequence>,
    pub total_length: usize,
    pub union_find: BidirectedUnionFind,
    pub sparsity_threshold: u64,
    /// Maps each position to its orientation relative to its union representative
    /// true = reverse orientation relative to representative
    pub orientation_map: HashMap<Pos, bool>,
}

impl SeqRush {
    pub fn new(sequences: Vec<Sequence>, sparsity_threshold: u64) -> Self {
        // Validate that no sequences are empty
        for seq in &sequences {
            if seq.data.is_empty() {
                panic!(
                    "Empty sequences are not allowed: sequence '{}' has length 0",
                    seq.id
                );
            }
        }

        let total_length = sequences.iter().map(|s| s.data.len()).sum();
        let union_find = BidirectedUnionFind::new(total_length);

        // Unite forward and reverse orientations of each position by default
        // This matches test expectations that F/R of the same position share a representative
        for i in 0..total_length {
            let f = make_pos(i, false);
            let r = make_pos(i, true);
            union_find.unite(f, r);
        }

        Self {
            sequences,
            total_length,
            union_find,
            sparsity_threshold,
            orientation_map: HashMap::new(),
        }
    }

    pub fn build_graph(&mut self, args: &Args) {
        println!(
            "Building graph with {} sequences (total length: {})",
            self.sequences.len(),
            self.total_length
        );

        // Phase 1: Align all pairs and update union-find
        self.align_and_unite(args);

        // Phase 2: Write graph by walking sequences through union-find
        self.write_gfa(args).expect("Failed to write GFA");
    }

    pub fn align_and_unite(&self, args: &Args) {
        if let Some(paf_path) = &args.paf {
            // Read alignments from PAF file
            self.align_and_unite_from_paf(paf_path, args);
        } else {
            // Check sequence lengths and provide guidance
            let min_len = self.sequences.iter().map(|s| s.data.len()).min().unwrap_or(0);
            let max_len = self.sequences.iter().map(|s| s.data.len()).max().unwrap_or(0);

            // Warn about very short sequences (potential alignment issues)
            if min_len < 100 && args.verbose {
                eprintln!("Warning: Found sequences shorter than 100 bp (min: {} bp, max: {} bp)",
                         min_len, max_len);
                eprintln!("         Short sequences may produce suboptimal alignments.");
            }

            // Auto-select aligner based on sequence length if using default
            let chosen_aligner = if args.aligner == "allwave" && min_len >= 100 {
                // Default is allwave, but for longer sequences, sweepga might be better
                #[cfg(feature = "use-sweepga")]
                {
                    if args.verbose {
                        eprintln!("Info: Sequences ≥100 bp detected. Consider --aligner sweepga for faster alignment.");
                    }
                }
                "allwave"
            } else {
                args.aligner.as_str()
            };

            // Use the configured aligner backend
            match chosen_aligner {
                #[cfg(feature = "use-allwave")]
                "allwave" => self.align_and_unite_with_allwave(args),

                #[cfg(feature = "use-sweepga")]
                "sweepga" => self.align_and_unite_with_sweepga(args),

                _ => {
                    eprintln!("Error: Unknown aligner '{}'. Available aligners:", chosen_aligner);
                    #[cfg(feature = "use-allwave")]
                    eprintln!("  - allwave (works with any sequence length)");
                    #[cfg(feature = "use-sweepga")]
                    eprintln!("  - sweepga (requires sequences ≥20 bp, recommended for ≥100 bp)");
                    std::process::exit(1);
                }
            }
        }
    }

    pub fn align_and_unite_from_paf(&self, paf_path: &str, args: &Args) {
        let verbose = args.verbose;
        if verbose {
            println!("Reading alignments from PAF file: {}", paf_path);
        }

        // Build a map from sequence names to indices
        let seq_name_to_idx: HashMap<String, usize> = self
            .sequences
            .iter()
            .enumerate()
            .map(|(idx, seq)| (seq.id.clone(), idx))
            .collect();

        // Read PAF file
        let file = File::open(paf_path).expect("Failed to open PAF file");
        let reader = BufReader::new(file);
        let mut alignment_count = 0;

        for line in reader.lines() {
            let line = line.expect("Failed to read PAF line");
            if line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 12 {
                eprintln!("Warning: Invalid PAF line (not enough fields): {}", line);
                continue;
            }

            // Parse PAF fields
            let query_name = fields[0];
            let _query_len = fields[1].parse::<usize>().unwrap();
            let query_start = fields[2].parse::<usize>().unwrap();
            let query_end = fields[3].parse::<usize>().unwrap();
            let query_strand = fields[4];
            let target_name = fields[5];
            let _target_len = fields[6].parse::<usize>().unwrap();
            let target_start = fields[7].parse::<usize>().unwrap();
            let target_end = fields[8].parse::<usize>().unwrap();

            // Find CIGAR string in optional fields
            let mut cigar = None;
            for field in &fields[12..] {
                if field.starts_with("cg:Z:") {
                    cigar = Some(&field[5..]);
                    break;
                }
            }

            let cigar = cigar.unwrap_or_else(|| {
                eprintln!("Warning: No CIGAR string found in PAF line: {}", line);
                ""
            });

            // Look up sequences by name
            let query_idx = seq_name_to_idx.get(query_name);
            let target_idx = seq_name_to_idx.get(target_name);

            if query_idx.is_none() || target_idx.is_none() {
                eprintln!(
                    "Warning: Unknown sequence name(s) in PAF: {} or {}",
                    query_name, target_name
                );
                continue;
            }

            let query_idx = *query_idx.unwrap();
            let target_idx = *target_idx.unwrap();
            let query_is_rc = query_strand == "-";

            if verbose {
                println!(
                    "Processing alignment: {} vs {} (query_is_rc: {})",
                    query_name, target_name, query_is_rc
                );
            }

            // Process the alignment
            self.process_alignment(
                cigar,
                &self.sequences[query_idx],
                &self.sequences[target_idx],
                args.min_match_length,
                query_is_rc,
                query_start,
                query_end,
                target_start,
                target_end,
                verbose,
            );

            alignment_count += 1;
        }

        if verbose {
            println!("Processed {} alignments from PAF file", alignment_count);
        }
    }

    #[cfg(feature = "use-allwave")]
    pub fn align_and_unite_with_allwave(&self, args: &Args) {
        // Parse alignment scores
        let scores = match AlignmentScores::parse(&args.scores) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Error parsing scores: {}", e);
                std::process::exit(1);
            }
        };

        // Parse orientation check scores
        let orientation_scores = match AlignmentScores::parse_orientation(&args.orientation_scores)
        {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Error parsing orientation scores: {}", e);
                std::process::exit(1);
            }
        };

        if args.verbose {
            println!("Using alignment scores: {:?}", scores);
            println!("Using orientation scores: {:?}", orientation_scores);
        }

        // Convert sequences to allwave format
        let allwave_sequences: Vec<allwave::Sequence> = self
            .sequences
            .iter()
            .map(|s| allwave::Sequence {
                id: s.id.clone(),
                seq: s.data.clone(),
            })
            .collect();

        // Create alignment parameters
        let params = AlignmentParams {
            match_score: scores.match_score,
            mismatch_penalty: scores.mismatch_penalty,
            gap_open: scores.gap1_open,
            gap_extend: scores.gap1_extend,
            gap2_open: scores.gap2_open,
            gap2_extend: scores.gap2_extend,
            max_divergence: args.max_divergence,
        };

        let orientation_params = AlignmentParams {
            match_score: orientation_scores.match_score,
            mismatch_penalty: orientation_scores.mismatch_penalty,
            gap_open: orientation_scores.gap1_open,
            gap_extend: orientation_scores.gap1_extend,
            gap2_open: orientation_scores.gap2_open,
            gap2_extend: orientation_scores.gap2_extend,
            max_divergence: None,
        };

        // Parse sparsification strategy
        let sparsification = match args.sparsification.as_str() {
            "1.0" => SparsificationStrategy::None,
            "auto" => SparsificationStrategy::Auto,
            s => match s.parse::<f64>() {
                Ok(factor) if factor > 0.0 && factor <= 1.0 => {
                    SparsificationStrategy::Random(factor)
                }
                _ => {
                    eprintln!(
                        "Invalid sparsification factor: {}. Using 1.0 (no sparsification)",
                        s
                    );
                    SparsificationStrategy::None
                }
            },
        };

        // Create PAF writer if requested
        let paf_writer = if let Some(paf_path) = &args.output_alignments {
            match File::create(paf_path) {
                Ok(file) => {
                    println!("Writing alignments to {}", paf_path);
                    Some(Arc::new(Mutex::new(BufWriter::new(file))))
                }
                Err(e) => {
                    eprintln!("Warning: Failed to create PAF file {}: {}", paf_path, e);
                    None
                }
            }
        } else {
            None
        };

        // Handle PAF output if requested
        if let Some(writer) = &paf_writer {
            let paf_aligner = AllPairIterator::with_options(
                &allwave_sequences,
                params.clone(),
                false, // include self-alignments
                false, // use_mash_orientation
                SparsificationStrategy::None, // No sparsification for PAF
            )
            .with_orientation_params(orientation_params.clone());

            // Convert to parallel iterator and process
            paf_aligner.into_par_iter().for_each(|alignment| {
                let paf_record = allwave::alignment_to_paf(&alignment, &allwave_sequences);
                if let Ok(mut w) = writer.lock() {
                    let _ = writeln!(w, "{}", paf_record);
                }
            });

            // Flush PAF writer
            if let Ok(mut w) = writer.lock() {
                let _ = w.flush();
            }
        }

        // For graph construction, we do all-vs-all including self
        let n = self.sequences.len();
        let total_pairs = n * n; // all-vs-all including self

        println!(
            "Total sequence pairs: {} (sparsification: {:?})",
            total_pairs, sparsification
        );

        // Create aligner for graph construction
        let graph_aligner = AllPairIterator::with_options(
            &allwave_sequences,
            params,
            false, // include self-alignments
            false, // use_mash_orientation
            sparsification,
        )
        .with_orientation_params(orientation_params);

        // Process alignments for graph construction
        graph_aligner.into_par_iter().for_each(|alignment| {
            // Process all alignments (both directions) for graph construction
            let cigar = allwave::cigar_bytes_to_string(&alignment.cigar_bytes);
            let query_seq = &self.sequences[alignment.query_idx];
            let target_seq = &self.sequences[alignment.target_idx];
            // Allwave aligns full sequences, so coordinates are 0 to length
            self.process_alignment(
                &cigar,
                query_seq,
                target_seq,
                args.min_match_length,
                alignment.is_reverse,
                0,
                query_seq.data.len(),
                0,
                target_seq.data.len(),
                args.verbose,
            );
        });
    }

    #[cfg(feature = "use-sweepga")]
    pub fn align_and_unite_with_sweepga(&self, args: &Args) {
        use std::path::Path;
        use sweepga::fastga_integration::FastGAIntegration;
        use sweepga::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};
        use tempfile::NamedTempFile;

        // Guard against very short sequences that cause FastGA to segfault
        let min_len = self.sequences.iter().map(|s| s.data.len()).min().unwrap_or(0);
        const MIN_SEQ_LENGTH: usize = 20;

        if min_len < MIN_SEQ_LENGTH {
            eprintln!(
                "Error: SweepGA aligner requires sequences of at least {} bp.",
                MIN_SEQ_LENGTH
            );
            eprintln!("       Found sequence with length {} bp.", min_len);
            eprintln!("       Use --aligner allwave for shorter sequences.");
            std::process::exit(1);
        }

        if args.verbose {
            println!("Using SweepGA aligner (FastGA + plane sweep filtering)");
            if let Some(freq) = args.frequency {
                println!("K-mer frequency threshold: {}", freq);
            }
        }

        // Get input FASTA path
        let fasta_path = Path::new(&args.sequences);

        // Run FastGA alignment using modern sweepga API
        if args.verbose {
            println!("Running FastGA alignment on {}...", fasta_path.display());
        }

        let fastga = FastGAIntegration::new(args.frequency, args.threads);
        let temp_paf = match fastga.align_to_temp_paf(fasta_path, fasta_path) {
            Ok(paf) => paf,
            Err(e) => {
                eprintln!("FastGA alignment failed: {}", e);
                std::process::exit(1);
            }
        };

        if args.verbose {
            // Count raw alignments
            use std::fs;
            let paf_content = fs::read_to_string(temp_paf.path()).expect("Failed to read PAF");
            let line_count = paf_content.lines()
                .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
                .count();
            println!("FastGA produced {} raw alignments", line_count);
        }

        // Apply 1:1 plane sweep filtering for clean graph construction
        if args.verbose {
            println!("Applying 1:1 plane sweep filtering...");
        }

        let filter_config = FilterConfig {
            mapping_filter_mode: FilterMode::OneToOne,
            mapping_max_per_query: Some(1),
            mapping_max_per_target: Some(1),
            scoring_function: ScoringFunction::LogLengthIdentity,
            overlap_threshold: 0.95,
            min_block_length: 100,  // Filter out very short alignments
            scaffold_filter_mode: FilterMode::ManyToMany,  // No scaffolding
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

        let temp_filtered = NamedTempFile::new().expect("Failed to create temp file");
        let filter = PafFilter::new(filter_config);
        if let Err(e) = filter.filter_paf(
            temp_paf.path().to_str().unwrap(),
            temp_filtered.path().to_str().unwrap()
        ) {
            eprintln!("PAF filtering failed: {}", e);
            std::process::exit(1);
        }

        if args.verbose {
            // Count filtered alignments
            use std::fs;
            let filtered_content = fs::read_to_string(temp_filtered.path()).expect("Failed to read filtered PAF");
            let filtered_count = filtered_content.lines()
                .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
                .count();
            println!("After 1:1 filtering: {} alignments", filtered_count);
        }

        // Now process the filtered PAF using existing align_and_unite_from_paf
        self.align_and_unite_from_paf(temp_filtered.path().to_str().unwrap(), args);
    }

    fn process_alignment(
        &self,
        cigar: &str,
        seq1: &Sequence,
        seq2: &Sequence,
        min_match_len: usize,
        query_is_rc: bool,
        query_start: usize,
        query_end: usize,
        target_start: usize,
        target_end: usize,
        verbose: bool,
    ) {
        if verbose && (seq1.id == seq2.id || seq1.id.contains("seq1") || seq2.id.contains("seq1")) {
            eprintln!(
                "Processing alignment: {} vs {} (query_is_rc: {})",
                seq1.id, seq2.id, query_is_rc
            );
        }
        // IMPORTANT: query_is_rc means the QUERY (seq1) was reverse complemented for alignment
        // Debug specific alignment
        let debug_this = query_is_rc && verbose; // Debug RC alignments

        // CRITICAL FIX: Don't RC the entire sequence!
        // Instead, we'll access it backwards and RC bases on-the-fly
        // This matches seqwish's approach and preserves correct offsets

        // Helper to get base at position (with RC if needed)
        let get_query_base = |local_pos: usize| -> u8 {
            if query_is_rc {
                // Access backwards and RC the base
                let base = seq1.data[seq1.data.len() - 1 - local_pos];
                match base {
                    b'A' | b'a' => b'T',
                    b'T' | b't' => b'A',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    _ => base,
                }
            } else {
                seq1.data[local_pos]
            }
        };

        // Validation function - CRITICAL! (user: "oh shit. you do need validation btw. be careful.")
        let validate_match = |start1: usize, start2: usize, len: usize| {
            for i in 0..len {
                let base1 = get_query_base(start1 + i);
                let base2 = seq2.data[start2 + i];

                if base1 != base2 {
                    panic!(
                        "VALIDATION ERROR: Attempting to unite non-matching bases!\n\
                         Alignment: {} -> {}\n\
                         Position {}: query_base[{}]='{}' != target[{}]='{}'\n\
                         Match region: query[{}..{}] -> target[{}..{}] len={}\n\
                         query_is_rc: {}",
                        seq1.id,
                        seq2.id,
                        i,
                        start1 + i,
                        base1 as char,
                        start2 + i,
                        base2 as char,
                        start1,
                        start1 + len,
                        start2,
                        start2 + len,
                        len,
                        query_is_rc
                    );
                }
            }
        };

        // CRITICAL: Start from PAF coordinates, not from 0!
        let mut pos1 = query_start;
        let mut pos2 = target_start;
        let mut count = 0;

        // Track consecutive matches across multiple M operations
        let mut in_match_run = false;
        let mut match_run_start1 = 0;
        let mut match_run_start2 = 0;
        let mut match_run_len = 0;

        if debug_this {
            eprintln!("\n=== DEBUG ALIGNMENT ===");
            eprintln!(
                "Seq1: {} (len: {}, offset: {})",
                seq1.id,
                seq1.data.len(),
                seq1.offset
            );
            eprintln!(
                "Seq2: {} (len: {}, offset: {}) QUERY_RC={}",
                seq2.id,
                seq2.data.len(),
                seq2.offset,
                query_is_rc
            );
            eprintln!("CIGAR: {}", cigar);
            eprintln!(
                "Seq1 start: {}",
                String::from_utf8_lossy(&seq1.data[0..20.min(seq1.data.len())])
            );
            // Debug output removed - seq1_data no longer exists
            eprintln!(
                "Seq2 start: {}",
                String::from_utf8_lossy(&seq2.data[0..20.min(seq2.data.len())])
            );
        }

        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                count = count * 10 + (ch as usize - '0' as usize);
            } else {
                if count == 0 {
                    count = 1;
                }

                match ch {
                    'M' | '=' => {
                        if debug_this && pos1 < 100 {
                            eprintln!("  CIGAR {}{}  pos1={} pos2={}", count, ch, pos1, pos2);
                        }
                        // if verbose {
                        //     println!("    Processing {}M operation at pos1={}, pos2={}", count, pos1, pos2);
                        // }
                        // For 'M' operations, we need to check if bases actually match
                        // since 'M' can represent either match or mismatch

                        // Check all positions in this M operation
                        for k in 0..count {
                            if pos1 + k < seq1.data.len() && pos2 + k < seq2.data.len() {
                                let base1 = get_query_base(pos1 + k);
                                let base2 = seq2.data[pos2 + k];
                                if debug_this && ch == '=' && k < 3 {
                                    eprintln!(
                                        "    '=' op: comparing seq1_data[{}]='{}' vs seq2[{}]='{}'",
                                        pos1 + k,
                                        base1 as char,
                                        pos2 + k,
                                        base2 as char
                                    );
                                }
                                // if verbose && k < 3 {
                                //     println!("      pos {}+{}: {} vs {} = {}", pos1, k, base1 as char, base2 as char, base1 == base2);
                                // }
                                if base1 == base2 {
                                    // Bases match - extend or start match run
                                    if debug_this && pos1 + k < 20 {
                                        eprintln!(
                                            "    Match at pos1={} pos2={}: '{}' == '{}'",
                                            pos1 + k,
                                            pos2 + k,
                                            base1 as char,
                                            base2 as char
                                        );
                                    }
                                    if !in_match_run {
                                        in_match_run = true;
                                        match_run_start1 = pos1 + k;
                                        match_run_start2 = pos2 + k;
                                        match_run_len = 1;
                                        // if verbose {
                                        //     println!("      Starting match run at pos1={}, pos2={}", match_run_start1, match_run_start2);
                                        // }
                                    } else {
                                        match_run_len += 1;
                                    }
                                } else {
                                    // Mismatch - process any accumulated match run
                                    if debug_this && in_match_run {
                                        eprintln!("    Breaking match run at mismatch: len={} >= min_len={}? {}", 
                                            match_run_len, min_match_len, match_run_len >= min_match_len);
                                    }
                                    if in_match_run && match_run_len >= min_match_len {
                                        // if verbose {
                                        //     println!("    Uniting match region: seq1[{}..{}] <-> seq2[{}..{}] (rc={})",
                                        //             match_run_start1, match_run_start1 + match_run_len,
                                        //             match_run_start2, match_run_start2 + match_run_len,
                                        //             query_is_rc);
                                        // }
                                        // Validate before uniting
                                        validate_match(
                                            match_run_start1,
                                            match_run_start2,
                                            match_run_len,
                                        );
                                        if debug_this {
                                            eprintln!("UNITING: seq1[{}..{}] <-> seq2[{}..{}] len={} (RC={})",
                                                match_run_start1, match_run_start1 + match_run_len,
                                                match_run_start2, match_run_start2 + match_run_len,
                                                match_run_len, query_is_rc);
                                        }
                                        if query_is_rc {
                                            // Query (seq1) was RC'd for alignment
                                            self.union_find.unite_matching_region(
                                                seq1.offset,
                                                seq2.offset,
                                                match_run_start1,
                                                match_run_start2,
                                                match_run_len,
                                                true,            // seq1 is RC
                                                seq1.data.len(), // seq1 length for RC transformation
                                            );
                                        } else {
                                            // Normal forward-forward alignment
                                            self.union_find.unite_matching_region(
                                                seq1.offset,
                                                seq2.offset,
                                                match_run_start1,
                                                match_run_start2,
                                                match_run_len,
                                                false,
                                                0, // Not used for forward alignments
                                            );
                                        }
                                    }
                                    in_match_run = false;
                                    match_run_len = 0;
                                }
                            }
                        }

                        pos1 += count;
                        pos2 += count;
                    }
                    _ => {
                        // Not a match operation - process any accumulated match run
                        if in_match_run && match_run_len >= min_match_len {
                            // if verbose {
                            //     println!("    Uniting match region: seq1[{}..{}] <-> seq2[{}..{}] (rc={})",
                            //             match_run_start1, match_run_start1 + match_run_len,
                            //             match_run_start2, match_run_start2 + match_run_len,
                            //             query_is_rc);
                            // }
                            // Validate before uniting
                            validate_match(
                                match_run_start1,
                                match_run_start2,
                                match_run_len,
                            );
                            if debug_this {
                                eprintln!("UNITING (end of op): seq1[{}..{}] <-> seq2[{}..{}] len={} (RC={})",
                                    match_run_start1, match_run_start1 + match_run_len,
                                    match_run_start2, match_run_start2 + match_run_len,
                                    match_run_len, query_is_rc);
                            }
                            if query_is_rc {
                                // Query (seq1) was RC'd for alignment
                                self.union_find.unite_matching_region(
                                    seq1.offset,
                                    seq2.offset,
                                    match_run_start1,
                                    match_run_start2,
                                    match_run_len,
                                    true,            // seq1 is RC
                                    seq1.data.len(), // seq1 length for RC transformation
                                );
                            } else {
                                // Normal forward-forward alignment
                                self.union_find.unite_matching_region(
                                    seq1.offset,
                                    seq2.offset,
                                    match_run_start1,
                                    match_run_start2,
                                    match_run_len,
                                    false,
                                    0, // Not used for forward alignments
                                );
                            }
                        }
                        in_match_run = false;
                        match_run_len = 0;

                        match ch {
                            'X' => {
                                // Mismatch
                                if debug_this && pos1 < 100 {
                                    eprintln!("  CIGAR {}X  pos1={} pos2={}", count, pos1, pos2);
                                }
                                pos1 += count;
                                pos2 += count;
                            }
                            'I' => {
                                // Insertion in query (seq1) relative to target (seq2)
                                if debug_this && pos1 < 100 {
                                    eprintln!("  CIGAR {}I  pos1={} pos2={}", count, pos1, pos2);
                                }
                                pos1 += count;
                            }
                            'D' => {
                                // Deletion in query (seq1) relative to target (seq2)
                                if debug_this && pos1 < 100 {
                                    eprintln!("  CIGAR {}D  pos1={} pos2={}", count, pos1, pos2);
                                }
                                pos2 += count;
                            }
                            _ => {}
                        }
                    }
                }
                count = 0;
            }
        }

        // Process final match run if alignment ends with matches
        // if verbose {
        //     println!("    End of CIGAR: in_match_run={}, match_run_len={}, min_match_len={}",
        //             in_match_run, match_run_len, min_match_len);
        // }
        if in_match_run && match_run_len >= min_match_len {
            // if verbose {
            //     println!("    Final uniting match region: seq1[{}..{}] <-> seq2[{}..{}] (rc={})",
            //             match_run_start1, match_run_start1 + match_run_len,
            //             match_run_start2, match_run_start2 + match_run_len,
            //             query_is_rc);
            // }
            // Validate before uniting
            validate_match(
                match_run_start1,
                match_run_start2,
                match_run_len,
            );
            if debug_this {
                eprintln!(
                    "UNITING (final): seq1[{}..{}] <-> seq2[{}..{}] len={} (RC={})",
                    match_run_start1,
                    match_run_start1 + match_run_len,
                    match_run_start2,
                    match_run_start2 + match_run_len,
                    match_run_len,
                    query_is_rc
                );
            }
            self.union_find.unite_matching_region(
                seq1.offset,
                seq2.offset,
                match_run_start1,
                match_run_start2,
                match_run_len,
                query_is_rc,
                seq1.data.len(), // Fix: should be seq1 length for RC transformation
            );
        }
    }

    fn write_gfa(&self, args: &Args) -> Result<(), Box<dyn std::error::Error>> {
        let output_path = &args.output;
        let verbose = args.verbose;
        let _test_mode = args.test_mode;

        // Use direct BidirectedGraph writing (no HashGraph conversion)
        self.write_bidirected_gfa(
            output_path,
            args.no_compact,
            args.no_sort,
            args.skip_sgd,
            args.skip_groom,
            args.skip_topo,
            args.sgd_iter_max,
            args.sgd_eta_max,
            args.sgd_theta,
            args.sgd_eps,
            args.sgd_cooling_start,
            args.threads,
            verbose,
        )
    }

    fn build_initial_graph(&self, verbose: bool) -> Result<Graph, Box<dyn std::error::Error>> {
        // Build bidirected graph first, then convert to simple graph for compatibility
        let bi_graph = self.build_bidirected_graph(verbose)?;
        // For now, convert back to simple graph until we update the rest of the pipeline
        Ok(self.bidirected_to_simple_graph(bi_graph))
    }

    fn build_old_graph(&self, verbose: bool) -> Result<Graph, Box<dyn std::error::Error>> {
        let mut graph = Graph::new();

        // Track which union representatives we've seen and their node IDs
        let mut union_to_node: HashMap<Pos, usize> = HashMap::new();
        let mut next_node_id = 1;

        // Build paths and discover nodes
        let mut all_union_reps = std::collections::HashSet::new();
        for seq in &self.sequences {
            let mut path = Vec::new();

            for i in 0..seq.data.len() {
                // Check both forward and reverse orientations to find the union representative
                let pos_fwd = make_pos(seq.offset + i, false);
                let pos_rev = make_pos(seq.offset + i, true);

                // Find union representatives for both orientations
                let union_fwd = self.union_find.find(pos_fwd);
                let union_rev = self.union_find.find(pos_rev);

                // Use the one that has been united with other sequences
                // (has a different representative than itself)
                if verbose && i < 5 {
                    eprintln!("  [SEQRUSH_DEBUG] Path building: {} pos {} (offset {}) - fwd: {} -> {}, rev: {} -> {}", 
                        seq.id, i, seq.offset + i, pos_fwd, union_fwd, pos_rev, union_rev);

                    // Also check what seq1's reverse positions unite to
                    if i < seq.data.len() {
                        let seq1_offset = 0; // seq1 starts at offset 0
                        let seq1_pos_rev = make_pos(seq1_offset + (11 - i), true); // RC mapping
                        let seq1_union = self.union_find.find(seq1_pos_rev);
                        eprintln!(
                            "    [SEQRUSH_DEBUG] seq1_rev pos {} -> union {}",
                            seq1_pos_rev, seq1_union
                        );
                    }
                }

                // Choose the union representative that has the most connections
                // If a position was involved in an RC alignment, one of its orientations
                // will have a different union representative than itself
                let union_rep = if union_fwd != pos_fwd && union_rev != pos_rev {
                    // Both orientations were united - choose the one with smaller ID
                    // to ensure consistency across sequences
                    std::cmp::min(union_fwd, union_rev)
                } else if union_fwd != pos_fwd {
                    union_fwd
                } else if union_rev != pos_rev {
                    union_rev
                } else {
                    // Neither orientation was united, use forward by default
                    union_fwd
                };

                if verbose && i < 5 {
                    eprintln!("    -> chosen union_rep {}", union_rep);
                }
                all_union_reps.insert(union_rep);

                // if verbose {
                //     eprintln!("DEBUG: seq {} pos {} -> union_rep {} (offset={}, rev={})",
                //              seq.id, seq.offset + i, union_rep, union_rep >> 1, union_rep & 1);
                // }

                // Get or create node ID for this union representative
                let node_id = match union_to_node.get(&union_rep) {
                    Some(&id) => {
                        if verbose && i < 5 {
                            eprintln!(
                                "    [SEQRUSH_DEBUG] Found existing node {} for union_rep {}",
                                id, union_rep
                            );
                        }
                        // Node already exists - verify consistency
                        if let Some(node) = graph.nodes.get_mut(&id) {
                            // Check if the base matches what we expect
                            let current_base = seq.data[i];
                            if !node.sequence.is_empty() && node.sequence[0] != current_base {
                                // This can happen when forward and reverse complement positions
                                // are united - we need to choose a canonical representation
                                if verbose {
                                    eprintln!("INFO: Union representative {} contains positions with different bases: {} and {} (choosing first)", 
                                             union_rep, node.sequence[0] as char, current_base as char);
                                }
                                // Keep the first base we saw for this union
                            } else if node.sequence.is_empty() {
                                // Initialize with this base
                                node.sequence = vec![current_base];
                            }
                        }
                        id
                    }
                    None => {
                        // First time seeing this union representative - create node
                        let id = next_node_id;
                        next_node_id += 1;
                        if verbose && i < 5 {
                            eprintln!(
                                "    [SEQRUSH_DEBUG] Creating new node {} for union_rep {}",
                                id, union_rep
                            );
                        }
                        union_to_node.insert(union_rep, id);

                        // Create node with single character sequence
                        let node = Node {
                            id,
                            sequence: vec![seq.data[i]],
                            rank: id as f64, // Initial rank based on creation order
                        };
                        graph.nodes.insert(id, node);
                        id
                    }
                };

                // Always add every position to the path - this preserves sequence length
                path.push(node_id);
            }

            graph.paths.push((seq.id.clone(), path));
        }

        // Build edges from paths (including self-loops for consecutive duplicates)
        for (_, path) in &graph.paths {
            for window in path.windows(2) {
                if let [from, to] = window {
                    graph.edges.insert(Edge {
                        from: *from,
                        to: *to,
                    });
                }
            }
        }

        // Basic path integrity check during initial graph construction
        if verbose {
            for (seq, (_, path)) in self.sequences.iter().zip(&graph.paths) {
                if let Err(error) = graph.verify_path_integrity(&seq.id, path, &seq.data) {
                    eprintln!("WARNING: {}", error);
                } else {
                    println!("✓ Path integrity verified for sequence {}", seq.id);
                }
            }

            eprintln!(
                "Total unique union representatives: {}",
                all_union_reps.len()
            );
            if all_union_reps.len() <= 20 {
                let mut reps: Vec<_> = all_union_reps.iter().collect();
                reps.sort();
                eprintln!("Union representatives: {:?}", reps);
            }
        }

        Ok(graph)
    }

    fn write_graph_to_gfa(
        &self,
        graph: &Graph,
        output_path: &str,
        _verbose: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);

        // GFA header
        writeln!(writer, "H\tVN:Z:1.0")?;

        // Write segments (nodes) in sorted order
        let mut sorted_nodes: Vec<_> = graph.nodes.values().collect();
        sorted_nodes.sort_by_key(|n| n.id);

        for node in sorted_nodes {
            let seq = String::from_utf8_lossy(&node.sequence);
            writeln!(writer, "S\t{}\t{}", node.id, seq)?;
        }

        // Write edges in sorted order (L records must come before P records)
        let mut sorted_edges: Vec<_> = graph.edges.iter().collect();
        sorted_edges.sort_by_key(|e| (e.from, e.to));

        for edge in sorted_edges {
            writeln!(writer, "L\t{}\t+\t{}\t+\t0M", edge.from, edge.to)?;
        }

        // Write paths (P records come last)
        for (seq_id, path) in &graph.paths {
            if !path.is_empty() {
                let path_str = path
                    .iter()
                    .map(|&node| format!("{}+", node))
                    .collect::<Vec<_>>()
                    .join(",");
                writeln!(writer, "P\t{}\t{}\t*", seq_id, path_str)?;
            }
        }

        println!(
            "Graph written to {}: {} nodes, {} edges",
            output_path,
            graph.nodes.len(),
            graph.edges.len()
        );

        Ok(())
    }

    fn compute_path_hashes(&self, graph: &BidirectedGraph) -> HashMap<String, Vec<u8>> {
        let mut hashes = HashMap::new();

        for path in &graph.paths {
            let sequence = path.get_sequence(|id| graph.nodes.get(id).and_then(|n| n.as_ref()));
            let mut hasher = Sha256::new();
            hasher.update(&sequence);
            let hash = hasher.finalize().to_vec();
            hashes.insert(path.name.clone(), hash);
        }

        hashes
    }

    fn compute_graph_size(&self, graph: &BidirectedGraph) -> usize {
        graph.nodes.iter().filter_map(|n| n.as_ref()).map(|node| node.sequence.len()).sum()
    }

    pub fn validate_paths_match_sequences(
        &self,
        graph: &BidirectedGraph,
    ) -> Result<(), Box<dyn std::error::Error>> {
        for seq in &self.sequences {
            // Find corresponding path in graph
            let path = graph
                .paths
                .iter()
                .find(|p| p.name == seq.id)
                .ok_or_else(|| format!("Path '{}' not found in graph", seq.id))?;

            // Extract sequence from path
            let path_sequence = path.get_sequence(|id| graph.nodes.get(id).and_then(|n| n.as_ref()));

            // Compare with original
            if path_sequence != seq.data {
                // Compute hashes for better error reporting
                let mut orig_hasher = Sha256::new();
                orig_hasher.update(&seq.data);
                let orig_hash = orig_hasher.finalize();

                let mut path_hasher = Sha256::new();
                path_hasher.update(&path_sequence);
                let path_hash = path_hasher.finalize();

                eprintln!("ERROR: Path '{}' does not match original sequence!", seq.id);
                eprintln!("  Original length: {} bp", seq.data.len());
                eprintln!("  Path length:     {} bp", path_sequence.len());

                // Show first difference and context
                for (i, (&orig, &path)) in seq.data.iter().zip(path_sequence.iter()).enumerate() {
                    if orig != path {
                        eprintln!(
                            "  First difference at position {}: '{}' (expected) vs '{}' (got)",
                            i, orig as char, path as char
                        );

                        // Show context around the difference
                        let start = i.saturating_sub(5);
                        let end = (i + 6).min(seq.data.len()).min(path_sequence.len());
                        eprintln!("  Context:");
                        eprintln!("    Original[{}..{}]: '{}'",
                            start, end,
                            String::from_utf8_lossy(&seq.data[start..end]));
                        eprintln!("    Path[{}..{}]:     '{}'",
                            start, end,
                            String::from_utf8_lossy(&path_sequence[start..end]));
                        break;
                    }
                }

                return Err(format!("Path '{}' does not match original sequence", seq.id).into());
            }
        }

        Ok(())
    }
}

pub fn load_sequences(file_path: &str) -> Result<Vec<Sequence>, Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut sequences = Vec::new();
    let mut current_id = String::new();
    let mut current_data = Vec::new();
    let mut offset = 0;

    for line in reader.lines() {
        let line = line?;
        if let Some(stripped) = line.strip_prefix('>') {
            if !current_id.is_empty() {
                sequences.push(Sequence {
                    id: current_id.clone(),
                    data: current_data.clone(),
                    offset,
                });
                offset += current_data.len();
                current_data.clear();
            }
            // Only take the first word (before any whitespace) as the sequence ID
            current_id = stripped.split_whitespace().next().unwrap_or("").to_string();
        } else {
            current_data.extend(line.trim().bytes());
        }
    }

    if !current_id.is_empty() {
        sequences.push(Sequence {
            id: current_id,
            data: current_data,
            offset,
        });
    }

    Ok(sequences)
}

pub fn run_seqrush(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    // Only initialize thread pool if not already initialized
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global();

    let sequences = load_sequences(&args.sequences)?;
    println!("Loaded {} sequences", sequences.len());

    // Calculate sparsification threshold
    let sparsity_threshold = if args.sparsification == "auto" {
        // Use Erdős-Rényi model: keep 10 * log(n) / n fraction of alignment pairs
        let n = sequences.len() as f64;
        let fraction = (10.0 * n.ln() / n).min(1.0);
        println!(
            "Auto sparsification: keeping {:.2}% of alignment pairs (n={})",
            fraction * 100.0,
            sequences.len()
        );
        (fraction * u64::MAX as f64) as u64
    } else {
        match args.sparsification.parse::<f64>() {
            Ok(frac) if (0.0..=1.0).contains(&frac) => {
                if frac == 1.0 {
                    u64::MAX
                } else {
                    println!(
                        "Sparsification: keeping {:.2}% of alignment pairs",
                        frac * 100.0
                    );
                    (frac * u64::MAX as f64) as u64
                }
            }
            _ => {
                eprintln!(
                    "Invalid sparsification value: {}. Using 1.0 (no sparsification)",
                    args.sparsification
                );
                u64::MAX
            }
        }
    };

    let mut seqrush = SeqRush::new(sequences, sparsity_threshold);
    seqrush.build_graph(&args);

    println!("Graph written to {}", args.output);
    Ok(())
}
