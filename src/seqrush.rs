use clap::Parser;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use std::fs::File;
use std::io::{BufWriter, Write, BufRead, BufReader};
use uf_rush::UFRush;
use lib_wfa2::affine_wavefront::{AffineWavefronts, MemoryMode, AlignmentStatus};

#[derive(Parser)]
#[command(name = "seqrush", version = "0.4.0", about = "Dynamic pangenome graph construction")]
pub struct Args {
    /// Input FASTA file
    #[arg(short, long)]
    pub sequences: String,
    
    /// Output GFA file
    #[arg(short, long, default_value = "output.gfa")]
    pub output: String,
    
    /// Number of threads
    #[arg(short, long, default_value = "4")]
    pub threads: usize,
    
    /// Minimum match length
    #[arg(short = 'k', long, default_value = "15")]
    pub min_match_length: usize,
    
    /// Alignment scores (match,mismatch,gap_open,gap_extend)
    #[arg(short = 'S', long = "scores", default_value = "0,5,8,2")]
    pub scores: String,
    
    /// Maximum divergence threshold (0.0-1.0, e.g., 0.1 = 10% divergence)
    #[arg(short = 'd', long = "max-divergence")]
    pub max_divergence: Option<f64>,
    
    /// Verbose output
    #[arg(short, long)]
    pub verbose: bool,
    
    /// Test mode - disable bidirectional alignment
    #[arg(long, hide = true)]
    pub test_mode: bool,
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
            return Err("Scores must have at least 4 values: match,mismatch,gap1_open,gap1_extend".to_string());
        }
        
        if parts.len() > 6 {
            return Err("Too many score values provided (max 6)".to_string());
        }
        
        let match_score = parts[0].parse::<i32>()
            .map_err(|_| format!("Invalid match score: {}", parts[0]))?;
        let mismatch_penalty = parts[1].parse::<i32>()
            .map_err(|_| format!("Invalid mismatch penalty: {}", parts[1]))?;
        let gap1_open = parts[2].parse::<i32>()
            .map_err(|_| format!("Invalid gap1_open penalty: {}", parts[2]))?;
        let gap1_extend = parts[3].parse::<i32>()
            .map_err(|_| format!("Invalid gap1_extend penalty: {}", parts[3]))?;
        
        let (gap2_open, gap2_extend) = if parts.len() >= 6 {
            (
                Some(parts[4].parse::<i32>()
                    .map_err(|_| format!("Invalid gap2_open penalty: {}", parts[4]))?),
                Some(parts[5].parse::<i32>()
                    .map_err(|_| format!("Invalid gap2_extend penalty: {}", parts[5]))?)
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
    pub offset: usize,  // Offset in concatenated sequence space
}

impl Sequence {
    /// Get reverse complement of the sequence
    pub fn reverse_complement(&self) -> Vec<u8> {
        self.data.iter().rev().map(|&base| {
            match base {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                _ => base,  // Keep N and other ambiguous bases as-is
            }
        }).collect()
    }
}

pub struct SeqRush {
    pub sequences: Vec<Sequence>,
    pub total_length: usize,
    pub union_find: Arc<UFRush>,
}

impl SeqRush {
    pub fn new(sequences: Vec<Sequence>) -> Self {
        let total_length = sequences.iter().map(|s| s.data.len()).sum();
        let union_find = Arc::new(UFRush::new(total_length));
        
        Self {
            sequences,
            total_length,
            union_find,
        }
    }
    
    pub fn build_graph(&mut self, args: &Args) {
        println!("Building graph with {} sequences (total length: {})", 
                 self.sequences.len(), self.total_length);
        
        // Phase 1: Align all pairs and update union-find
        self.align_and_unite(args);
        
        // Phase 2: Write graph by walking sequences through union-find
        self.write_gfa(&args.output, args.verbose).expect("Failed to write GFA");
    }
    
    pub fn align_and_unite(&self, args: &Args) {
        // Parse alignment scores
        let scores = match AlignmentScores::parse(&args.scores) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Error parsing scores: {}", e);
                std::process::exit(1);
            }
        };
        
        if args.verbose {
            println!("Using alignment scores: {:?}", scores);
        }
        
        let n = self.sequences.len();
        let pairs: Vec<(usize, usize)> = (0..n)
            .flat_map(|i| (i+1..n).map(move |j| (i, j)))
            .collect();
        
        println!("Aligning {} sequence pairs", pairs.len());
        
        // Process pairs in parallel
        pairs.par_iter().for_each(|&(i, j)| {
            self.align_pair(i, j, args, &scores);
        });
    }
    
    fn align_pair(&self, idx1: usize, idx2: usize, args: &Args, scores: &AlignmentScores) {
        let seq1 = &self.sequences[idx1];
        let seq2 = &self.sequences[idx2];
        
        if args.verbose {
            println!("Aligning {} vs {}", seq1.id, seq2.id);
        }
        
        // Calculate max score threshold if divergence is specified
        let max_score = args.max_divergence.map(|div| {
            let avg_len = (seq1.data.len() + seq2.data.len()) / 2;
            scores.max_score_for_divergence(avg_len, div)
        });
        
        // Create WFA aligner
        let mut wf = if scores.gap2_open.is_some() && scores.gap2_extend.is_some() {
            // Two-piece affine gap model
            // TODO: Update when WFA2 bindings support two-piece affine
            // For now, fall back to single affine
            AffineWavefronts::with_penalties(
                scores.match_score,
                scores.mismatch_penalty,
                scores.gap1_open,
                scores.gap1_extend
            )
        } else {
            // Single affine gap model
            AffineWavefronts::with_penalties(
                scores.match_score,
                scores.mismatch_penalty,
                scores.gap1_open,
                scores.gap1_extend
            )
        };
        wf.set_memory_mode(MemoryMode::Ultralow);
        
        // Try forward-forward alignment
        let status_ff = wf.align(&seq1.data, &seq2.data);
        let score_ff = if matches!(status_ff, AlignmentStatus::Completed) {
            wf.score().abs()  // Convert to positive for easier comparison
        } else {
            i32::MAX
        };
        
        // Try reverse complement alignment unless in test mode
        let (_status_fr, score_fr) = if !args.test_mode {
            let seq2_rc = seq2.reverse_complement();
            let status = wf.align(&seq1.data, &seq2_rc);
            let score = if matches!(status, AlignmentStatus::Completed) {
                wf.score().abs()  // Convert to positive for easier comparison
            } else {
                i32::MAX
            };
            (status, score)
        } else {
            (AlignmentStatus::Undefined, i32::MAX)
        };
        
        if args.verbose {
            println!("  Forward-Forward score: {}", score_ff);
            println!("  Forward-Reverse score: {}", score_fr);
        }
        
        // Check if either alignment meets threshold
        if let Some(threshold) = max_score {
            if score_ff > threshold && score_fr > threshold {
                if args.verbose {
                    println!("  Both alignments exceed divergence threshold ({})", threshold);
                }
                return;
            }
        }
        
        // Use the better alignment
        if score_ff <= score_fr && score_ff != i32::MAX {
            // Re-run forward-forward to get CIGAR
            let status = wf.align(&seq1.data, &seq2.data);
            if matches!(status, AlignmentStatus::Completed) {
                let cigar_bytes = wf.cigar();
                let cigar = String::from_utf8_lossy(cigar_bytes);
                if args.verbose {
                    println!("  Using Forward-Forward alignment");
                    println!("  CIGAR: {}", cigar);
                }
                self.process_alignment(&cigar, seq1, seq2, args.min_match_length, false);
            }
        } else if score_fr != i32::MAX && !args.test_mode {
            // Re-run forward-reverse to get CIGAR
            let seq2_rc = seq2.reverse_complement();
            let status = wf.align(&seq1.data, &seq2_rc);
            if matches!(status, AlignmentStatus::Completed) {
                let cigar_bytes = wf.cigar();
                let cigar = String::from_utf8_lossy(cigar_bytes);
                if args.verbose {
                    println!("  Using Forward-Reverse alignment");
                    println!("  CIGAR: {}", cigar);
                }
                // Note: seq2 is reverse complemented in this alignment
                self.process_alignment(&cigar, seq1, seq2, args.min_match_length, true);
            }
        }
    }
    
    fn process_alignment(&self, cigar: &str, seq1: &Sequence, seq2: &Sequence, min_match_len: usize, seq2_is_rc: bool) {
        let mut pos1 = 0;
        let mut pos2 = 0;
        let mut count = 0;
        
        // Track consecutive matches across multiple M operations
        let mut in_match_run = false;
        let mut match_run_start1 = 0;
        let mut match_run_start2 = 0;
        let mut match_run_len = 0;
        
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                count = count * 10 + (ch as usize - '0' as usize);
            } else {
                if count == 0 { count = 1; }
                
                match ch {
                    'M' | '=' => {
                        // For 'M' operations, we need to check if bases actually match
                        // since 'M' can represent either match or mismatch
                        
                        // Check all positions in this M operation
                        for k in 0..count {
                            if pos1 + k < seq1.data.len() && pos2 + k < seq2.data.len() {
                                if seq1.data[pos1 + k] == seq2.data[pos2 + k] {
                                    // Bases match - extend or start match run
                                    if !in_match_run {
                                        in_match_run = true;
                                        match_run_start1 = pos1 + k;
                                        match_run_start2 = pos2 + k;
                                        match_run_len = 1;
                                    } else {
                                        match_run_len += 1;
                                    }
                                } else {
                                    // Mismatch - process any accumulated match run
                                    if in_match_run && match_run_len >= min_match_len {
                                        for j in 0..match_run_len {
                                            let global_pos1 = seq1.offset + match_run_start1 + j;
                                            let global_pos2 = if seq2_is_rc {
                                                // Map reverse complement position back to forward strand
                                                seq2.offset + (seq2.data.len() - 1 - (match_run_start2 + j))
                                            } else {
                                                seq2.offset + match_run_start2 + j
                                            };
                                            self.union_find.unite(global_pos1, global_pos2);
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
                            for j in 0..match_run_len {
                                let global_pos1 = seq1.offset + match_run_start1 + j;
                                let global_pos2 = if seq2_is_rc {
                                    // Map reverse complement position back to forward strand
                                    seq2.offset + (seq2.data.len() - 1 - (match_run_start2 + j))
                                } else {
                                    seq2.offset + match_run_start2 + j
                                };
                                self.union_find.unite(global_pos1, global_pos2);
                            }
                        }
                        in_match_run = false;
                        match_run_len = 0;
                        
                        match ch {
                            'X' => { 
                                // Mismatch
                                pos1 += count; 
                                pos2 += count; 
                            }
                            'I' => { 
                                // Insertion in seq2
                                pos2 += count; 
                            }
                            'D' => { 
                                // Deletion in seq2
                                pos1 += count; 
                            }
                            _ => {}
                        }
                    }
                }
                count = 0;
            }
        }
        
        // Process final match run if alignment ends with matches
        if in_match_run && match_run_len >= min_match_len {
            for j in 0..match_run_len {
                let global_pos1 = seq1.offset + match_run_start1 + j;
                let global_pos2 = if seq2_is_rc {
                    // Map reverse complement position back to forward strand
                    seq2.offset + (seq2.data.len() - 1 - (match_run_start2 + j))
                } else {
                    seq2.offset + match_run_start2 + j
                };
                self.union_find.unite(global_pos1, global_pos2);
            }
        }
    }
    
    fn write_gfa(&self, output_path: &str, verbose: bool) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);
        
        // GFA header
        writeln!(writer, "H\tVN:Z:1.0")?;
        
        // Track which unions we've seen and their node IDs
        let mut union_to_node: HashMap<usize, usize> = HashMap::new();
        let mut node_sequences: HashMap<usize, Vec<u8>> = HashMap::new();
        let mut next_node_id = 1;
        
        // Build paths and discover nodes
        let mut paths = Vec::new();
        
        for seq in &self.sequences {
            let mut path = Vec::new();
            let mut last_node = None;
            
            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                let union = self.union_find.find(global_pos);
                
                // Get or create node ID for this union
                let node_id = match union_to_node.get(&union) {
                    Some(&id) => {
                        // Verify this character matches the node's sequence
                        if let Some(node_seq) = node_sequences.get(&id) {
                            if !node_seq.is_empty() && node_seq[0] != seq.data[i] {
                                eprintln!("WARNING: Character mismatch in union {}: {} vs {}", 
                                         union, node_seq[0] as char, seq.data[i] as char);
                            }
                        }
                        id
                    }
                    None => {
                        // First time seeing this union - create node
                        let id = next_node_id;
                        next_node_id += 1;
                        union_to_node.insert(union, id);
                        
                        // The character at this position becomes the node's sequence
                        node_sequences.insert(id, vec![seq.data[i]]);
                        id
                    }
                };
                
                // Add to path if different from previous
                if last_node != Some(node_id) {
                    path.push(node_id);
                    last_node = Some(node_id);
                }
            }
            
            paths.push((seq.id.clone(), path));
        }
        
        // Write segments (nodes)
        for node_id in 1..next_node_id {
            let seq = node_sequences.get(&node_id)
                .map(|s| String::from_utf8_lossy(s))
                .unwrap_or_else(|| "N".into());
            writeln!(writer, "S\t{}\t{}", node_id, seq)?;
        }
        
        // Write paths and verify integrity
        for (seq, path) in self.sequences.iter().zip(&paths) {
            let (seq_id, path) = (&seq.id, &path.1);
            if !path.is_empty() {
                // Verify path integrity - reconstruct sequence from graph
                if verbose {
                    let mut reconstructed = Vec::new();
                    for &node_id in path {
                        if let Some(node_seq) = node_sequences.get(&node_id) {
                            reconstructed.extend_from_slice(node_seq);
                        }
                    }
                    
                    if reconstructed != seq.data {
                        eprintln!("WARNING: Path integrity check failed for sequence {}", seq_id);
                        eprintln!("  Original length: {}, Reconstructed length: {}", 
                                 seq.data.len(), reconstructed.len());
                        
                        // Find first mismatch
                        for (i, (&orig, &recon)) in seq.data.iter().zip(&reconstructed).enumerate() {
                            if orig != recon {
                                eprintln!("  First mismatch at position {}: {} vs {}", 
                                         i, orig as char, recon as char);
                                break;
                            }
                        }
                    } else if verbose {
                        println!("Path integrity verified for sequence {}", seq_id);
                    }
                }
                
                let path_str = path.iter()
                    .map(|&node| format!("{}+", node))
                    .collect::<Vec<_>>()
                    .join(",");
                writeln!(writer, "P\t{}\t{}\t*", seq_id, path_str)?;
            }
        }
        
        // Write edges
        let mut edges = HashSet::new();
        for (_, path) in &paths {
            for window in path.windows(2) {
                if let [from, to] = window {
                    edges.insert((*from, *to));
                }
            }
        }
        
        for (from, to) in edges {
            writeln!(writer, "L\t{}\t+\t{}\t+\t0M", from, to)?;
        }
        
        println!("Created {} nodes", next_node_id - 1);
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
        if line.starts_with('>') {
            if !current_id.is_empty() {
                sequences.push(Sequence {
                    id: current_id.clone(),
                    data: current_data.clone(),
                    offset,
                });
                offset += current_data.len();
                current_data.clear();
            }
            current_id = line[1..].trim().to_string();
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
    
    let mut seqrush = SeqRush::new(sequences);
    seqrush.build_graph(&args);
    
    println!("Graph written to {}", args.output);
    Ok(())
}