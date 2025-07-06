use crate::pos::*;
use crate::bidirected_graph::{Handle, BiPath, reverse_complement};
use crate::bidirected_ops::BidirectedGraph;
use crate::seqrush::{Args, AlignmentScores, Sequence};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, BufRead, BufReader};
use uf_rush::UFRush;
use lib_wfa2::affine_wavefront::{AffineWavefronts, MemoryMode, AlignmentStatus};

pub struct SeqRushBidirected {
    sequences: Vec<Sequence>,
    union_find: UFRush,
    #[allow(dead_code)]
    total_positions: usize,
}

impl SeqRushBidirected {
    pub fn new(sequences: Vec<Sequence>) -> Self {
        // Calculate total positions (each position can be forward or reverse)
        let total_positions: usize = sequences.iter().map(|s| s.data.len()).sum::<usize>() * 2;
        
        let union_find = UFRush::new(total_positions);
        
        SeqRushBidirected {
            sequences,
            union_find,
            total_positions,
        }
    }
    
    /// Convert an oriented position to a union-find index
    fn pos_to_uf_index(&self, pos: Pos) -> usize {
        // Each sequence position has two indices: forward and reverse
        offset(pos) * 2 + (is_rev(pos) as usize)
    }
    
    /// Convert a union-find index back to an oriented position
    #[allow(dead_code)]
    fn uf_index_to_pos(&self, index: usize) -> Pos {
        let offset = index / 2;
        let is_reverse = (index % 2) == 1;
        make_pos(offset, is_reverse)
    }
    
    pub fn build_graph(&mut self, args: &Args) -> Result<(), Box<dyn std::error::Error>> {
        if args.verbose {
            println!("Building bidirectional pangenome graph...");
        }
        
        // Parse alignment scores
        let scores = AlignmentScores::parse(&args.scores)?;
        
        // Set number of threads
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
        
        // Process all pairs of sequences
        self.process_alignments(&scores, args.min_match_length, args.max_divergence, args.verbose)?;
        
        // Write the graph
        self.write_gfa(&args.output, args.verbose)?;
        
        Ok(())
    }
    
    fn process_alignments(
        &mut self, 
        scores: &AlignmentScores,
        min_match_length: usize,
        max_divergence: Option<f64>,
        verbose: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let n = self.sequences.len();
        
        // Generate all pairs
        let pairs: Vec<(usize, usize)> = (0..n)
            .flat_map(|i| (i..n).map(move |j| (i, j)))
            .collect();
        
        if verbose {
            println!("Processing {} sequence pairs", pairs.len());
        }
        
        // Clone sequences for parallel processing
        let sequences = self.sequences.clone();
        
        // Process pairs in parallel
        let results: Vec<_> = pairs.par_iter()
            .flat_map(|&(i, j)| {
                let mut pair_results = Vec::new();
                
                // Always do forward-forward alignment
                if let Ok(result) = Self::align_sequences_static(&sequences[i], &sequences[j], false, false, scores, max_divergence) {
                    if result.score >= 0 {
                        pair_results.push((i, j, false, false, result.cigar.clone(), result.score));
                    }
                }
                
                // For different sequences, also try forward-reverse
                if i != j {
                    if let Ok(result) = Self::align_sequences_static(&sequences[i], &sequences[j], false, true, scores, max_divergence) {
                        if result.score >= 0 {
                            pair_results.push((i, j, false, true, result.cigar.clone(), result.score));
                        }
                    }
                }
                
                pair_results
            })
            .collect();
        
        if verbose {
            println!("Generated {} alignments", results.len());
        }
        
        // Process alignments sequentially to update union-find
        for (i, j, _seq1_is_rc, seq2_is_rc, cigar, score) in results {
            if verbose && i != j {
                let orientation = if seq2_is_rc { "forward-reverse" } else { "forward-forward" };
                println!("  Processing alignment {} to {} ({}): score={}, cigar={}", 
                         self.sequences[i].id, self.sequences[j].id, orientation, score, cigar);
            }
            
            // Clone sequences to avoid borrowing issues
            let seq1 = self.sequences[i].clone();
            let seq2 = self.sequences[j].clone();
            
            self.process_cigar_bidirected(
                &seq1, 
                &seq2, 
                &cigar, 
                seq2_is_rc,
                min_match_length,
                verbose
            );
        }
        
        Ok(())
    }
    
    fn align_sequences_static(
        seq1: &Sequence,
        seq2: &Sequence,
        seq1_is_rc: bool,
        seq2_is_rc: bool,
        scores: &AlignmentScores,
        max_divergence: Option<f64>,
    ) -> Result<AlignmentResult, Box<dyn std::error::Error>> {
        let pattern = if seq1_is_rc {
            reverse_complement(&seq1.data)
        } else {
            seq1.data.clone()
        };
        
        let text = if seq2_is_rc {
            reverse_complement(&seq2.data)
        } else {
            seq2.data.clone()
        };
        
        // Create wavefront aligner
        let mut aligner = if scores.gap2_open.is_some() && scores.gap2_extend.is_some() {
            AffineWavefronts::with_penalties_affine2p(
                scores.match_score,
                scores.mismatch_penalty,
                scores.gap1_open,
                scores.gap1_extend,
                scores.gap2_open.unwrap(),
                scores.gap2_extend.unwrap()
            )
        } else {
            AffineWavefronts::with_penalties(
                scores.match_score,
                scores.mismatch_penalty,
                scores.gap1_open,
                scores.gap1_extend
            )
        };
        aligner.set_memory_mode(MemoryMode::Ultralow);
        
        // Perform alignment
        let status = aligner.align(&pattern, &text);
        
        match status {
            AlignmentStatus::Completed => {
                let score = aligner.score();
                
                // Check divergence if specified
                if let Some(max_div) = max_divergence {
                    let max_allowed_score = scores.max_score_for_divergence(pattern.len().min(text.len()), max_div);
                    if score > max_allowed_score {
                        return Ok(AlignmentResult {
                            score: -1,
                            cigar: String::new(),
                        });
                    }
                }
                
                let cigar = String::from_utf8_lossy(&aligner.cigar()).to_string();
                Ok(AlignmentResult { score, cigar })
            }
            _ => Err(format!("Alignment failed with status: {:?}", status).into()),
        }
    }
    
    fn process_cigar_bidirected(
        &mut self,
        seq1: &Sequence,
        seq2: &Sequence,
        cigar: &str,
        seq2_is_rc: bool,
        min_match_length: usize,
        verbose: bool,
    ) {
        if verbose {
            println!("    CIGAR: {}", cigar);
        }
        
        let mut pos1 = 0;
        let mut pos2 = 0;
        let mut count = 0;
        
        // Track match runs
        let mut in_match_run = false;
        let mut match_run_start1 = 0;
        let mut match_run_start2 = 0;
        let mut match_run_len = 0;
        
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                count = count * 10 + (ch as usize - '0' as usize);
            } else {
                // If no number precedes the operation, count is 1
                if count == 0 {
                    count = 1;
                }
                match ch {
                    'M' | '=' => {
                        // Process matches in runs
                        if !in_match_run {
                            in_match_run = true;
                            match_run_start1 = pos1;
                            match_run_start2 = pos2;
                            match_run_len = count;
                        } else {
                            match_run_len += count;
                        }
                        
                        // Check if we should process this run
                        if match_run_len >= min_match_length {
                            if verbose {
                                println!("    Processing match run of length {} at seq1[{}] seq2[{}]", 
                                         match_run_len, match_run_start1, match_run_start2);
                            }
                            // Process the match run
                            for j in 0..match_run_len {
                                // Create oriented positions
                                let pos1_fwd = make_pos(seq1.offset + match_run_start1 + j, false);
                                
                                let pos2_oriented = if seq2_is_rc {
                                    // For reverse complement, we need to map back
                                    let rev_offset = seq2.offset + (seq2.data.len() - 1 - (match_run_start2 + j));
                                    make_pos(rev_offset, true)
                                } else {
                                    make_pos(seq2.offset + match_run_start2 + j, false)
                                };
                                
                                if verbose && j < 2 {
                                    println!("      Unite: pos1[{}]={}{} <-> pos2[{}]={}{}", 
                                             offset(pos1_fwd), seq1.data[match_run_start1 + j] as char, orientation_char(pos1_fwd),
                                             offset(pos2_oriented), seq2.data[match_run_start2 + j] as char, orientation_char(pos2_oriented));
                                }
                                
                                // Unite the oriented positions
                                let idx1 = self.pos_to_uf_index(pos1_fwd);
                                let idx2 = self.pos_to_uf_index(pos2_oriented);
                                self.union_find.unite(idx1, idx2);
                                
                                // Also create the implicit reverse complement relationships
                                // If A+ matches B+, then A- matches B-
                                // If A+ matches B-, then A- matches B+
                                let pos1_rev = flip_orientation(pos1_fwd);
                                let pos2_flipped = flip_orientation(pos2_oriented);
                                
                                let idx1_rev = self.pos_to_uf_index(pos1_rev);
                                let idx2_flipped = self.pos_to_uf_index(pos2_flipped);
                                self.union_find.unite(idx1_rev, idx2_flipped);
                            }
                        }
                        
                        pos1 += count;
                        pos2 += count;
                    }
                    _ => {
                        // Not a match - reset match run tracking
                        in_match_run = false;
                        match_run_len = 0;
                        
                        match ch {
                            'X' => { 
                                pos1 += count; 
                                pos2 += count; 
                            }
                            'I' => { 
                                pos2 += count; 
                            }
                            'D' => { 
                                pos1 += count; 
                            }
                            _ => {}
                        }
                    }
                }
                count = 0;
            }
        }
    }
    
    fn write_gfa(&self, output_path: &str, verbose: bool) -> Result<(), Box<dyn std::error::Error>> {
        // Build the bidirectional graph
        let graph = self.build_bidirected_graph(verbose)?;
        
        if verbose {
            println!("Built bidirectional graph: {} nodes, {} edges, {} paths", 
                     graph.nodes.len(), graph.edges.len(), graph.paths.len());
        }
        
        // Write to file
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);
        graph.write_gfa(&mut writer)?;
        
        Ok(())
    }
    
    fn build_bidirected_graph(&self, _verbose: bool) -> Result<BidirectedGraph, Box<dyn std::error::Error>> {
        let mut graph = BidirectedGraph::new();
        
        // Track which unions we've seen and their node IDs
        let mut union_to_node: HashMap<usize, usize> = HashMap::new();
        let mut next_node_id = 1;
        
        // Build paths and discover nodes
        for seq in &self.sequences {
            let mut path = BiPath::new(seq.id.clone());
            
            for i in 0..seq.data.len() {
                // Use forward orientation for the path
                let pos = make_pos(seq.offset + i, false);
                let uf_idx = self.pos_to_uf_index(pos);
                let union = self.union_find.find(uf_idx);
                
                // Get or create node for this union
                let node_id = match union_to_node.get(&union) {
                    Some(&id) => id,
                    None => {
                        let id = next_node_id;
                        next_node_id += 1;
                        union_to_node.insert(union, id);
                        
                        // Create node with the base at this position
                        graph.add_node(id, vec![seq.data[i]]);
                        id
                    }
                };
                
                // Add step to path (always forward for now)
                path.add_step(Handle::forward(node_id));
            }
            
            graph.paths.push(path);
        }
        
        // Build edges from paths
        let paths = graph.paths.clone();
        for path in &paths {
            for window in path.steps.windows(2) {
                if let [from, to] = window {
                    graph.add_edge(*from, *to);
                }
            }
        }
        
        Ok(graph)
    }
}

#[derive(Debug)]
struct AlignmentResult {
    score: i32,
    cigar: String,
}

/// Run the bidirectional version of SeqRush
pub fn run_seqrush_bidirected(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    // Read sequences
    let sequences = read_fasta(&args.sequences)?;
    
    if args.verbose {
        println!("Loaded {} sequences", sequences.len());
        for seq in &sequences {
            println!("  {}: {} bp", seq.id, seq.data.len());
        }
    }
    
    // Create and run SeqRush
    let mut seqrush = SeqRushBidirected::new(sequences);
    seqrush.build_graph(&args)?;
    
    Ok(())
}

fn read_fasta(path: &str) -> Result<Vec<Sequence>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    
    let mut sequences = Vec::new();
    let mut current_id = String::new();
    let mut current_seq = Vec::new();
    let mut offset = 0;
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_id.is_empty() {
                sequences.push(Sequence {
                    id: current_id.clone(),
                    data: current_seq.clone(),
                    offset,
                });
                offset += current_seq.len();
                current_seq.clear();
            }
            current_id = line[1..].split_whitespace().next().unwrap_or("").to_string();
        } else {
            current_seq.extend(line.trim().as_bytes());
        }
    }
    
    if !current_id.is_empty() {
        sequences.push(Sequence {
            id: current_id,
            data: current_seq,
            offset,
        });
    }
    
    Ok(sequences)
}