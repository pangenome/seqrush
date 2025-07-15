use crate::pos::*;
use crate::bidirected_graph::{Handle, BiPath, reverse_complement};
use crate::bidirected_ops::BidirectedGraph;
use crate::seqrush::{Args, AlignmentScores, Sequence};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, BufRead, BufReader};
use uf_rush::UFRush;
use lib_wfa2::affine_wavefront::{AffineWavefronts, MemoryMode, AlignmentStatus};

pub struct SeqRushBidirectedSimple {
    sequences: Vec<Sequence>,
    union_find: UFRush,
}

impl SeqRushBidirectedSimple {
    pub fn new(sequences: Vec<Sequence>) -> Self {
        // Calculate total positions (each position can be forward or reverse)
        let total_positions: usize = sequences.iter().map(|s| s.data.len()).sum::<usize>() * 2;
        
        let union_find = UFRush::new(total_positions);
        
        SeqRushBidirectedSimple {
            sequences,
            union_find,
        }
    }
    
    /// Convert an oriented position to a union-find index
    fn pos_to_uf_index(&self, pos: Pos) -> usize {
        // Each sequence position has two indices: forward and reverse
        offset(pos) * 2 + (is_rev(pos) as usize)
    }
    
    pub fn build_graph(&mut self, args: &Args) -> Result<(), Box<dyn std::error::Error>> {
        if args.verbose {
            println!("Building bidirectional pangenome graph...");
        }
        
        // Parse alignment scores
        let scores = AlignmentScores::parse(&args.scores)?;
        
        // Set number of threads (only if not already initialized)
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global();
        
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
        
        // Process all pairs
        for i in 0..n {
            for j in i..n {
                // Forward-forward alignment
                if let Ok((score, cigar)) = self.align_sequences(i, j, false, false, scores, max_divergence) {
                    if score >= 0 {
                        if verbose && i != j {
                            println!("  Aligned {} to {} (forward-forward): score={}", 
                                     self.sequences[i].id, self.sequences[j].id, score);
                        }
                        self.process_alignment(i, j, false, false, &cigar, min_match_length, verbose);
                    }
                }
                
                // Forward-reverse alignment for different sequences
                if i != j {
                    if let Ok((score, cigar)) = self.align_sequences(i, j, false, true, scores, max_divergence) {
                        if score >= 0 {
                            if verbose {
                                println!("  Aligned {} to {} (forward-reverse): score={}", 
                                         self.sequences[i].id, self.sequences[j].id, score);
                            }
                            self.process_alignment(i, j, false, true, &cigar, min_match_length, verbose);
                        }
                    }
                }
            }
        }
        
        Ok(())
    }
    
    fn align_sequences(
        &self,
        i: usize,
        j: usize,
        seq1_is_rc: bool,
        seq2_is_rc: bool,
        scores: &AlignmentScores,
        max_divergence: Option<f64>,
    ) -> Result<(i32, String), Box<dyn std::error::Error>> {
        let pattern = if seq1_is_rc {
            reverse_complement(&self.sequences[i].data)
        } else {
            self.sequences[i].data.clone()
        };
        
        let text = if seq2_is_rc {
            reverse_complement(&self.sequences[j].data)
        } else {
            self.sequences[j].data.clone()
        };
        
        // Create wavefront aligner
        let aligner = if scores.gap2_open.is_some() && scores.gap2_extend.is_some() {
            AffineWavefronts::with_penalties_affine2p_and_memory_mode(
                scores.match_score,
                scores.mismatch_penalty,
                scores.gap1_open,
                scores.gap1_extend,
                scores.gap2_open.unwrap(),
                scores.gap2_extend.unwrap(),
                MemoryMode::Ultralow
            )
        } else {
            AffineWavefronts::with_penalties_and_memory_mode(
                scores.match_score,
                scores.mismatch_penalty,
                scores.gap1_open,
                scores.gap1_extend,
                MemoryMode::Ultralow
            )
        };
        
        // Perform alignment
        let status = aligner.align(&pattern, &text);
        
        match status {
            AlignmentStatus::Completed => {
                let score = aligner.score();
                
                // Check divergence if specified
                if let Some(max_div) = max_divergence {
                    let max_allowed_score = scores.max_score_for_divergence(pattern.len().min(text.len()), max_div);
                    if score > max_allowed_score {
                        return Ok((-1, String::new()));
                    }
                }
                
                let cigar = String::from_utf8_lossy(&aligner.cigar()).to_string();
                Ok((score, cigar))
            }
            _ => Err(format!("Alignment failed with status: {:?}", status).into()),
        }
    }
    
    fn process_alignment(
        &mut self,
        seq1_idx: usize,
        seq2_idx: usize,
        _seq1_is_rc: bool,
        seq2_is_rc: bool,
        cigar: &str,
        min_match_length: usize,
        verbose: bool,
    ) {
        
        if verbose {
            println!("    Processing CIGAR: {}", cigar);
        }
        
        let mut pos1 = 0;
        let mut pos2 = 0;
        
        // Parse CIGAR and collect match runs
        let mut operations = Vec::new();
        let mut count = 0;
        
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                count = count * 10 + (ch as usize - '0' as usize);
            } else {
                if count == 0 {
                    count = 1;
                }
                operations.push((ch, count));
                count = 0;
            }
        }
        
        // Process operations and find match runs
        let mut in_match_run = false;
        let mut match_run_start1 = 0;
        let mut match_run_start2 = 0;
        let mut match_run_len = 0;
        
        for (op, count) in operations {
            match op {
                'M' | '=' => {
                    if !in_match_run {
                        in_match_run = true;
                        match_run_start1 = pos1;
                        match_run_start2 = pos2;
                        match_run_len = count;
                    } else {
                        match_run_len += count;
                    }
                    pos1 += count;
                    pos2 += count;
                }
                _ => {
                    // End of match run - process it if long enough
                    if in_match_run && match_run_len >= min_match_length {
                        self.process_match_run(
                            seq1_idx, seq2_idx, 
                            match_run_start1, match_run_start2, match_run_len,
                            seq2_is_rc, verbose
                        );
                    }
                    
                    in_match_run = false;
                    match_run_len = 0;
                    
                    match op {
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
        }
        
        // Process final match run if any
        if in_match_run && match_run_len >= min_match_length {
            self.process_match_run(
                seq1_idx, seq2_idx, 
                match_run_start1, match_run_start2, match_run_len,
                seq2_is_rc, verbose
            );
        }
    }
    
    fn process_match_run(
        &mut self,
        seq1_idx: usize,
        seq2_idx: usize,
        start1: usize,
        start2: usize,
        length: usize,
        seq2_is_rc: bool,
        verbose: bool,
    ) {
        let seq1 = &self.sequences[seq1_idx];
        let seq2 = &self.sequences[seq2_idx];
        if verbose {
            println!("    Processing match run: seq1[{}..{}] <-> seq2[{}..{}] (len={})", 
                     start1, start1 + length, start2, start2 + length, length);
        }
        
        for j in 0..length {
            let pos1_fwd = make_pos(seq1.offset + start1 + j, false);
            
            let pos2_oriented = if seq2_is_rc {
                let rev_offset = seq2.offset + (seq2.data.len() - 1 - (start2 + j));
                make_pos(rev_offset, true)
            } else {
                make_pos(seq2.offset + start2 + j, false)
            };
            
            if verbose && j < 3 {
                let base1 = seq1.data[start1 + j] as char;
                let actual_pos2 = if seq2_is_rc {
                    seq2.data.len() - 1 - (start2 + j)
                } else {
                    start2 + j
                };
                let base2 = seq2.data[actual_pos2] as char;
                
                let idx1 = self.pos_to_uf_index(pos1_fwd);
                let idx2 = self.pos_to_uf_index(pos2_oriented);
                
                println!("      Unite: seq1[{}]={}{} <-> seq2[{}]={}{} (actual seq2[{}]={}) [uf_idx: {} <-> {}]", 
                         offset(pos1_fwd), base1, orientation_char(pos1_fwd),
                         offset(pos2_oriented), base2, orientation_char(pos2_oriented),
                         actual_pos2, base2, idx1, idx2);
            }
            
            // Unite forward positions
            let idx1 = self.pos_to_uf_index(pos1_fwd);
            let idx2 = self.pos_to_uf_index(pos2_oriented);
            self.union_find.unite(idx1, idx2);
            
            // Unite reverse complement positions
            let pos1_rev = flip_orientation(pos1_fwd);
            let pos2_flipped = flip_orientation(pos2_oriented);
            
            let idx1_rev = self.pos_to_uf_index(pos1_rev);
            let idx2_flipped = self.pos_to_uf_index(pos2_flipped);
            self.union_find.unite(idx1_rev, idx2_flipped);
        }
    }
    
    fn write_gfa(&self, output_path: &str, verbose: bool) -> Result<(), Box<dyn std::error::Error>> {
        let graph = self.build_bidirected_graph(verbose)?;
        
        if verbose {
            println!("Built bidirectional graph: {} nodes, {} edges, {} paths", 
                     graph.nodes.len(), graph.edges.len(), graph.paths.len());
        }
        
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);
        graph.write_gfa(&mut writer)?;
        
        Ok(())
    }
    
    fn build_bidirected_graph(&self, verbose: bool) -> Result<BidirectedGraph, Box<dyn std::error::Error>> {
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
                
                if verbose && i < 3 {
                    println!("    Graph building: seq {} pos {} -> uf_idx {} -> union {}", 
                             seq.id, i, uf_idx, union);
                }
                
                // Get or create node for this union
                let node_id = match union_to_node.get(&union) {
                    Some(&id) => {
                        if verbose && i < 3 {
                            println!("    Reusing node {} for seq {} pos {}", id, seq.id, i);
                        }
                        id
                    },
                    None => {
                        let id = next_node_id;
                        next_node_id += 1;
                        union_to_node.insert(union, id);
                        
                        // Create node with the base at this position
                        graph.add_node(id, vec![seq.data[i]]);
                        if verbose && i < 3 {
                            println!("    Creating node {} for seq {} pos {} (union {})", id, seq.id, i, union);
                        }
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

/// Run the simplified bidirectional version of SeqRush
pub fn run_seqrush_bidirected_simple(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    // Read sequences
    let sequences = read_fasta(&args.sequences)?;
    
    if args.verbose {
        println!("Loaded {} sequences", sequences.len());
        for seq in &sequences {
            println!("  {}: {} bp", seq.id, seq.data.len());
        }
    }
    
    // Create and run SeqRush
    let mut seqrush = SeqRushBidirectedSimple::new(sequences);
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