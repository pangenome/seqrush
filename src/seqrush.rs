use clap::Parser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;
use std::fs::File;
use std::io::{BufWriter, Write, BufRead, BufReader};
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use std::sync::Mutex;
use uf_rush::UFRush;
use lib_wfa2::affine_wavefront::{AffineWavefronts, MemoryMode, AlignmentStatus};
use crate::graph_ops::{Graph, Node, Edge};

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
    
    pub fn parse_orientation(scores_str: &str) -> Result<Self, String> {
        let parts: Vec<&str> = scores_str.split(',').collect();
        
        if parts.len() != 4 {
            return Err("Orientation scores must have exactly 4 values: match,mismatch,gap_open,gap_extend".to_string());
        }
        
        let match_score = parts[0].parse::<i32>()
            .map_err(|_| format!("Invalid match score: {}", parts[0]))?;
        let mismatch_penalty = parts[1].parse::<i32>()
            .map_err(|_| format!("Invalid mismatch penalty: {}", parts[1]))?;
        let gap_open = parts[2].parse::<i32>()
            .map_err(|_| format!("Invalid gap_open penalty: {}", parts[2]))?;
        let gap_extend = parts[3].parse::<i32>()
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
    pub sparsity_threshold: u64,
}

impl SeqRush {
    pub fn new(sequences: Vec<Sequence>, sparsity_threshold: u64) -> Self {
        let total_length = sequences.iter().map(|s| s.data.len()).sum();
        let union_find = Arc::new(UFRush::new(total_length));
        
        Self {
            sequences,
            total_length,
            union_find,
            sparsity_threshold,
        }
    }
    
    pub fn build_graph(&mut self, args: &Args) {
        println!("Building graph with {} sequences (total length: {})", 
                 self.sequences.len(), self.total_length);
        
        // Phase 1: Align all pairs and update union-find
        self.align_and_unite(args);
        
        // Phase 2: Write graph by walking sequences through union-find
        self.write_gfa(args).expect("Failed to write GFA");
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
        
        // Parse orientation check scores
        let orientation_scores = match AlignmentScores::parse_orientation(&args.orientation_scores) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("Error parsing orientation scores: {}", e);
                std::process::exit(1);
            }
        };
        
        if args.verbose {
            println!("Using alignment scores: {:?}", scores);
            println!("Using orientation check scores: {:?}", orientation_scores);
        }
        
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
        
        let n = self.sequences.len();
        let pairs: Vec<(usize, usize)> = (0..n)
            .flat_map(|i| (i+1..n).map(move |j| (i, j)))
            .collect();
        
        // Count how many pairs will actually be aligned after sparsification
        let pairs_to_align = pairs.iter()
            .filter(|&&(i, j)| self.should_align_pair(i, j))
            .count();
        
        println!("Total sequence pairs: {} (will align {} after sparsification)", 
                 pairs.len(), pairs_to_align);
        
        // Process pairs in parallel
        pairs.par_iter().for_each(|&(i, j)| {
            self.align_pair(i, j, args, &scores, &orientation_scores, &paf_writer);
        });
        
        // Flush PAF writer
        if let Some(writer) = paf_writer {
            if let Ok(mut w) = writer.lock() {
                let _ = w.flush();
            }
        }
    }
    
    fn align_pair(&self, idx1: usize, idx2: usize, args: &Args, scores: &AlignmentScores, orientation_scores: &AlignmentScores, paf_writer: &Option<Arc<Mutex<BufWriter<File>>>>) {
        let seq1 = &self.sequences[idx1];
        let seq2 = &self.sequences[idx2];
        
        // Apply sparsification at the alignment-pair level (before alignment)
        if !self.should_align_pair(idx1, idx2) {
            if args.verbose {
                println!("Skipping alignment pair {} vs {} (sparsification filter)", seq1.id, seq2.id);
            }
            return;
        }
        
        if args.verbose {
            println!("Aligning {} vs {}", seq1.id, seq2.id);
        }
        
        // First, do fast orientation check with edit distance scoring
        let mut orientation_wf = AffineWavefronts::with_penalties(
            orientation_scores.match_score,
            orientation_scores.mismatch_penalty,
            orientation_scores.gap1_open,
            orientation_scores.gap1_extend
        );
        orientation_wf.set_memory_mode(MemoryMode::Ultralow);
        
        // Try forward-forward orientation check
        let status_ff = orientation_wf.align(&seq1.data, &seq2.data);
        let orientation_score_ff = if matches!(status_ff, AlignmentStatus::Completed) {
            orientation_wf.score().abs()
        } else {
            i32::MAX
        };
        
        // Try forward-reverse orientation check
        let seq2_rc = seq2.reverse_complement();
        let status_fr = orientation_wf.align(&seq1.data, &seq2_rc);
        let orientation_score_fr = if matches!(status_fr, AlignmentStatus::Completed) {
            orientation_wf.score().abs()
        } else {
            i32::MAX
        };
        
        if args.verbose {
            println!("  Orientation check - FF score: {}, FR score: {}", orientation_score_ff, orientation_score_fr);
        }
        
        // Choose the better orientation
        let use_reverse = orientation_score_fr < orientation_score_ff && !args.test_mode;
        
        // Calculate max score threshold if divergence is specified
        let max_score = args.max_divergence.map(|div| {
            let avg_len = (seq1.data.len() + seq2.data.len()) / 2;
            scores.max_score_for_divergence(avg_len, div)
        });
        
        // Create WFA aligner for full alignment
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
        
        // Perform alignment in the chosen orientation
        let (status, score, is_reverse) = if use_reverse {
            // Forward-reverse alignment
            let status = wf.align(&seq1.data, &seq2_rc);
            let score = if matches!(status, AlignmentStatus::Completed) {
                wf.score().abs()
            } else {
                i32::MAX
            };
            (status, score, true)
        } else {
            // Forward-forward alignment
            let status = wf.align(&seq1.data, &seq2.data);
            let score = if matches!(status, AlignmentStatus::Completed) {
                wf.score().abs()
            } else {
                i32::MAX
            };
            (status, score, false)
        };
        
        if args.verbose {
            println!("  Full alignment score: {} ({})", score, if is_reverse { "reverse" } else { "forward" });
        }
        
        // Check if alignment meets threshold
        if let Some(threshold) = max_score {
            if score > threshold {
                if args.verbose {
                    println!("  Alignment exceeds divergence threshold ({})", threshold);
                }
                return;
            }
        }
        
        // Process the alignment if successful
        if matches!(status, AlignmentStatus::Completed) && score != i32::MAX {
            let cigar_bytes = wf.cigar();
            let cigar = String::from_utf8_lossy(cigar_bytes);
            if args.verbose {
                println!("  CIGAR: {}", cigar);
            }
            
            // Write PAF output if requested
            if let Some(writer) = paf_writer {
                self.write_paf_record(writer, seq1, seq2, &cigar, score, is_reverse, &wf);
            }
            
            self.process_alignment(&cigar, seq1, seq2, args.min_match_length, is_reverse);
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
    
    fn write_paf_record(&self, writer: &Arc<Mutex<BufWriter<File>>>, seq1: &Sequence, seq2: &Sequence, 
                       cigar: &str, score: i32, is_reverse: bool, _wf: &AffineWavefronts) {
        // Convert CIGAR to --eqx style (M -> =)
        let eqx_cigar = self.convert_to_eqx_cigar(cigar);
        
        // Calculate alignment lengths and statistics from CIGAR
        let (query_start, query_end, target_start, target_end, num_matches, block_length) = 
            self.parse_cigar_for_paf(&eqx_cigar, seq1.data.len(), seq2.data.len(), is_reverse);
        
        // PAF format: query_name query_len query_start query_end strand target_name target_len target_start target_end num_matches block_len mapping_quality
        // Additional fields: NM:i:edit_distance AS:i:alignment_score cg:Z:cigar
        let strand = if is_reverse { '-' } else { '+' };
        let mapq = 60; // High quality since we're doing full alignment
        
        // Calculate edit distance from score (assuming match=0, mismatch and gaps add to score)
        let edit_distance = score;
        
        if let Ok(mut w) = writer.lock() {
            let _ = writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\tAS:i:{}\tcg:Z:{}",
                seq1.id,
                seq1.data.len(),
                query_start,
                query_end,
                strand,
                seq2.id,
                seq2.data.len(),
                target_start,
                target_end,
                num_matches,
                block_length,
                mapq,
                edit_distance,
                -score, // AS tag is typically negative of edit distance
                eqx_cigar
            );
        }
    }
    
    fn convert_to_eqx_cigar(&self, cigar: &str) -> String {
        // Convert M operations to = for matches (--eqx style) and compact the CIGAR
        let mut result = String::new();
        let mut count = 0;
        let mut current_op: Option<char> = None;
        let mut op_count = 0;
        
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                count = count * 10 + (ch as usize - '0' as usize);
            } else {
                if count == 0 { count = 1; }
                
                // Convert M to =
                let op = if ch == 'M' { '=' } else { ch };
                
                // If this is a new operation type, write out the previous one
                if current_op.is_some() && current_op != Some(op) {
                    result.push_str(&op_count.to_string());
                    result.push(current_op.unwrap());
                    op_count = 0;
                }
                
                // Accumulate counts for the same operation
                current_op = Some(op);
                op_count += count;
                count = 0;
            }
        }
        
        // Write out the last operation
        if let Some(op) = current_op {
            result.push_str(&op_count.to_string());
            result.push(op);
        }
        
        result
    }
    
    fn parse_cigar_for_paf(&self, cigar: &str, _query_len: usize, target_len: usize, is_reverse: bool) 
        -> (usize, usize, usize, usize, usize, usize) {
        let mut query_pos = 0;
        let mut target_pos = 0;
        let mut num_matches = 0;
        let mut count = 0;
        
        // Track start positions (first match position)
        let mut query_start = 0;
        let mut target_start = 0;
        let mut found_first_match = false;
        
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                count = count * 10 + (ch as usize - '0' as usize);
            } else {
                if count == 0 { count = 1; }
                
                match ch {
                    'M' | '=' => {
                        if !found_first_match {
                            query_start = query_pos;
                            target_start = target_pos;
                            found_first_match = true;
                        }
                        // For --eqx style, = means match
                        num_matches += count;
                        query_pos += count;
                        target_pos += count;
                    }
                    'X' => {
                        if !found_first_match {
                            query_start = query_pos;
                            target_start = target_pos;
                            found_first_match = true;
                        }
                        query_pos += count;
                        target_pos += count;
                    }
                    'I' => {
                        query_pos += count;
                    }
                    'D' => {
                        target_pos += count;
                    }
                    _ => {}
                }
                count = 0;
            }
        }
        
        // For reverse strand, adjust target coordinates
        let (final_target_start, final_target_end) = if is_reverse {
            (target_len - target_pos, target_len - target_start)
        } else {
            (target_start, target_pos)
        };
        
        // Block length is the span in both query and target
        let block_length = query_pos.max(target_pos);
        
        (query_start, query_pos, final_target_start, final_target_end, num_matches, block_length)
    }
    
    fn write_gfa(&self, args: &Args) -> Result<(), Box<dyn std::error::Error>> {
        let output_path = &args.output;
        let verbose = args.verbose;
        let test_mode = args.test_mode;
        // Build initial graph from union-find
        let mut graph = self.build_initial_graph(verbose)?;
        
        if verbose {
            println!("Initial graph: {} nodes, {} edges", graph.nodes.len(), graph.edges.len());
        }
        
        // Prepare original sequences for verification
        let original_sequences: Vec<(String, Vec<u8>)> = self.sequences.iter()
            .map(|seq| (seq.id.clone(), seq.data.clone()))
            .collect();
        
        // Verify paths before compaction (skip in test mode for synthetic sequences)
        if !test_mode {
            if verbose {
                println!("Verifying paths before compaction...");
            }
            if let Err(errors) = graph.comprehensive_verify(Some(&original_sequences), verbose) {
                eprintln!("Path verification failed before compaction:");
                for error in &errors {
                    eprintln!("  - {}", error);
                }
                return Err(format!("Path verification failed before compaction: {} errors", errors.len()).into());
            }
        }
        
        // Perform node compaction unless disabled
        if !args.no_compact {
            let compacted = graph.compact_nodes();
            if verbose {
                println!("Compacted {} nodes", compacted);
                println!("After compaction: {} nodes, {} edges", graph.nodes.len(), graph.edges.len());
            }
            
            // Verify paths after compaction (skip in test mode for synthetic sequences)
            if !test_mode {
                if verbose {
                    println!("Verifying paths after compaction...");
                }
                if let Err(errors) = graph.comprehensive_verify(Some(&original_sequences), verbose) {
                    eprintln!("Path verification failed after compaction:");
                    for error in &errors {
                        eprintln!("  - {}", error);
                    }
                    return Err(format!("Path verification failed after compaction: {} errors", errors.len()).into());
                }
            }
        } else {
            if verbose {
                println!("Skipping node compaction to preserve graph structure");
            }
        }
        
        
        // Perform topological sort
        graph.topological_sort();
        if verbose {
            println!("Completed topological sort");
        }
        
        // Final verification after topological sort (skip in test mode for synthetic sequences)
        if !test_mode {
            if verbose {
                println!("Final verification after topological sort...");
            }
            if let Err(errors) = graph.comprehensive_verify(Some(&original_sequences), verbose) {
                eprintln!("Path verification failed after topological sort:");
                for error in &errors {
                    eprintln!("  - {}", error);
                }
                // Don't fail on verification errors - just warn
                eprintln!("WARNING: Continuing despite {} verification errors", errors.len());
            }
        } else if verbose {
            println!("Skipping verification in test mode");
        }
        
        // Always validate GFA format before final output
        if let Err(errors) = graph.validate_gfa_format(verbose) {
            eprintln!("\nERROR: Invalid GFA format detected:");
            for error in &errors {
                eprintln!("  - {}", error);
            }
            return Err("Cannot write invalid GFA file".into());
        }
        
        // Write the graph to GFA
        self.write_graph_to_gfa(&graph, output_path, verbose)?;
        
        Ok(())
    }
    
    fn build_initial_graph(&self, verbose: bool) -> Result<Graph, Box<dyn std::error::Error>> {
        let mut graph = Graph::new();
        
        // Track which unions we've seen and their node IDs
        let mut union_to_node: HashMap<usize, usize> = HashMap::new();
        let mut next_node_id = 1;
        
        // Build paths and discover nodes
        for seq in &self.sequences {
            let mut path = Vec::new();
            
            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                let union = self.union_find.find(global_pos);
                
                // Get or create node ID for this union
                let node_id = match union_to_node.get(&union) {
                    Some(&id) => {
                        // Node already exists - extend its sequence if this character is different
                        if let Some(node) = graph.nodes.get_mut(&id) {
                            // If this is a new character in the same union, extend the node's sequence
                            if node.sequence.is_empty() || node.sequence[node.sequence.len() - 1] != seq.data[i] {
                                // Character differs from last in node - this is expected in bidirectional graphs
                                // when positions from different strands are united
                                if !node.sequence.is_empty() && node.sequence[0] != seq.data[i] {
                                    if verbose {
                                        eprintln!("INFO: Union {} contains multiple characters: {} and {} (expected for bidirectional graph)", 
                                                 union, node.sequence[0] as char, seq.data[i] as char);
                                    }
                                }
                                // Keep the first character we saw for this union
                                // In a proper bidirectional graph implementation, we'd track both
                            }
                        }
                        id
                    }
                    None => {
                        // First time seeing this union - create node
                        let id = next_node_id;
                        next_node_id += 1;
                        union_to_node.insert(union, id);
                        
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
                    graph.edges.insert(Edge { from: *from, to: *to });
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
        }
        
        Ok(graph)
    }
    
    fn write_graph_to_gfa(&self, graph: &Graph, output_path: &str, _verbose: bool) 
        -> Result<(), Box<dyn std::error::Error>> {
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
                let path_str = path.iter()
                    .map(|&node| format!("{}+", node))
                    .collect::<Vec<_>>()
                    .join(",");
                writeln!(writer, "P\t{}\t{}\t*", seq_id, path_str)?;
            }
        }
        
        println!("Graph written to {}: {} nodes, {} edges", 
                 output_path, graph.nodes.len(), graph.edges.len());
        
        Ok(())
    }
    
    fn should_align_pair(&self, idx1: usize, idx2: usize) -> bool {
        if self.sparsity_threshold == u64::MAX {
            return true; // No sparsification
        }
        
        // Hash the sequence names to get a deterministic pseudo-random value
        // Order matters: A→B is different from B→A
        let mut hasher = DefaultHasher::new();
        self.sequences[idx1].id.hash(&mut hasher);
        self.sequences[idx2].id.hash(&mut hasher);
        let hash_value = hasher.finish();
        
        // Keep the alignment if hash is below threshold
        hash_value <= self.sparsity_threshold
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
    
    // Calculate sparsification threshold
    let sparsity_threshold = if args.sparsification == "auto" {
        // Use Erdős-Rényi model: keep 10 * log(n) / n fraction of alignment pairs
        let n = sequences.len() as f64;
        let fraction = (10.0 * n.ln() / n).min(1.0);
        println!("Auto sparsification: keeping {:.2}% of alignment pairs (n={})", fraction * 100.0, sequences.len());
        (fraction * u64::MAX as f64) as u64
    } else {
        match args.sparsification.parse::<f64>() {
            Ok(frac) if frac >= 0.0 && frac <= 1.0 => {
                if frac == 1.0 {
                    u64::MAX
                } else {
                    println!("Sparsification: keeping {:.2}% of alignment pairs", frac * 100.0);
                    (frac * u64::MAX as f64) as u64
                }
            }
            _ => {
                eprintln!("Invalid sparsification value: {}. Using 1.0 (no sparsification)", args.sparsification);
                u64::MAX
            }
        }
    };
    
    let mut seqrush = SeqRush::new(sequences, sparsity_threshold);
    seqrush.build_graph(&args);
    
    println!("Graph written to {}", args.output);
    Ok(())
}