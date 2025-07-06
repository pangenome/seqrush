use crate::pos::*;
use crate::bidirected_graph::{Handle, BiPath, reverse_complement};
use crate::bidirected_ops::BidirectedGraph;
use crate::seqrush::{Args, AlignmentScores, Sequence};
use crate::cigar_analysis::{find_potential_inversion_sites, is_potential_inversion};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, BufRead, BufReader};
use uf_rush::UFRush;
use lib_wfa2::affine_wavefront::{AffineWavefronts, MemoryMode, AlignmentStatus};

pub struct InversionAwareSeqRush {
    sequences: Vec<Sequence>,
    union_find: UFRush,
}

#[derive(Debug)]
struct Alignment {
    seq1_idx: usize,
    seq2_idx: usize,
    seq1_start: usize,
    seq1_end: usize,
    seq2_start: usize,
    seq2_end: usize,
    seq2_is_rc: bool,
    cigar: String,
    score: i32,
}

impl InversionAwareSeqRush {
    pub fn new(sequences: Vec<Sequence>) -> Self {
        let total_positions: usize = sequences.iter().map(|s| s.data.len()).sum::<usize>() * 2;
        let union_find = UFRush::new(total_positions);
        
        InversionAwareSeqRush {
            sequences,
            union_find,
        }
    }
    
    fn pos_to_uf_index(&self, pos: Pos) -> usize {
        offset(pos) * 2 + (is_rev(pos) as usize)
    }
    
    pub fn build_graph(&mut self, args: &Args) -> Result<(), Box<dyn std::error::Error>> {
        if args.verbose {
            println!("Building inversion-aware bidirectional pangenome graph...");
        }
        
        let scores = AlignmentScores::parse(&args.scores)?;
        
        // Collect all alignments including patched inversions
        let alignments = self.find_all_alignments(&scores, args.min_match_length, args.max_divergence, args.verbose)?;
        
        if args.verbose {
            println!("Found {} alignments (including {} patched inversions)", 
                     alignments.len(),
                     alignments.iter().filter(|a| a.seq2_is_rc && (a.seq1_start > 0 || a.seq2_start > 0)).count());
        }
        
        // Process alignments to build union-find
        for alignment in &alignments {
            self.process_alignment(alignment, args.min_match_length, args.verbose);
        }
        
        // Write the graph
        self.write_gfa(&args.output, args.verbose)?;
        
        Ok(())
    }
    
    fn find_all_alignments(
        &self,
        scores: &AlignmentScores,
        min_match_length: usize,
        max_divergence: Option<f64>,
        verbose: bool,
    ) -> Result<Vec<Alignment>, Box<dyn std::error::Error>> {
        let mut all_alignments = Vec::new();
        let n = self.sequences.len();
        
        for i in 0..n {
            for j in i..n {
                if verbose && i != j {
                    println!("Aligning {} vs {}...", self.sequences[i].id, self.sequences[j].id);
                }
                // Full sequence alignments
                if let Ok(mut alignments) = self.align_sequences_with_inversions(
                    i, j, scores, min_match_length, max_divergence, verbose
                ) {
                    all_alignments.append(&mut alignments);
                }
            }
        }
        
        Ok(all_alignments)
    }
    
    fn align_sequences_with_inversions(
        &self,
        seq1_idx: usize,
        seq2_idx: usize,
        scores: &AlignmentScores,
        min_match_length: usize,
        max_divergence: Option<f64>,
        verbose: bool,
    ) -> Result<Vec<Alignment>, Box<dyn std::error::Error>> {
        let mut alignments = Vec::new();
        
        // First, try full forward-forward alignment
        match self.align_region(
            seq1_idx, 0, self.sequences[seq1_idx].data.len(),
            seq2_idx, 0, self.sequences[seq2_idx].data.len(),
            false, scores, max_divergence
        ) {
            Ok((score, cigar)) => {
                if score >= 0 {
                if verbose {
                    println!("  Forward-forward alignment {} vs {}: score={}, cigar={}",
                             self.sequences[seq1_idx].id, self.sequences[seq2_idx].id,
                             score, cigar);
                }
                
                let alignment = Alignment {
                    seq1_idx,
                    seq2_idx,
                    seq1_start: 0,
                    seq1_end: self.sequences[seq1_idx].data.len(),
                    seq2_start: 0,
                    seq2_end: self.sequences[seq2_idx].data.len(),
                    seq2_is_rc: false,
                    cigar: cigar.clone(),
                    score,
                };
                
                // Check for potential inversions in the alignment
                let gaps = find_potential_inversion_sites(&cigar, min_match_length * 2);
                if verbose && !gaps.is_empty() {
                    println!("  Found {} potential gaps for inversions", gaps.len());
                }
                let mut _has_inversions = false;
                
                for gap in &gaps {
                    if is_potential_inversion(gap, min_match_length * 2) {
                        _has_inversions = true;
                        if verbose {
                            println!("  Found potential inversion in {} vs {}: query[{}..{}] target[{}..{}]",
                                     self.sequences[seq1_idx].id, self.sequences[seq2_idx].id,
                                     gap.query_start, gap.query_end,
                                     gap.target_start, gap.target_end);
                        }
                        
                        // Try to align the gap region in reverse orientation
                        if let Ok((inv_score, inv_cigar)) = self.align_region(
                            seq1_idx, gap.query_start, gap.query_end,
                            seq2_idx, gap.target_start, gap.target_end,
                            true, scores, max_divergence
                        ) {
                            if inv_score >= 0 && inv_score < score / 2 { // Inversion alignment is good
                                if verbose {
                                    println!("    Inversion confirmed! Score: {}", inv_score);
                                }
                                
                                // Add the patched inversion alignment
                                alignments.push(Alignment {
                                    seq1_idx,
                                    seq2_idx,
                                    seq1_start: gap.query_start,
                                    seq1_end: gap.query_end,
                                    seq2_start: gap.target_start,
                                    seq2_end: gap.target_end,
                                    seq2_is_rc: true,
                                    cigar: inv_cigar,
                                    score: inv_score,
                                });
                            }
                        }
                    }
                }
                
                    // Add the main alignment (even if it has gaps)
                    alignments.push(alignment);
                }
            }
            Err(e) => {
                if verbose {
                    println!("  Forward-forward alignment failed: {}", e);
                }
            }
        }
        
        // For different sequences, also try full forward-reverse
        if seq1_idx != seq2_idx {
            if let Ok((score, cigar)) = self.align_region(
                seq1_idx, 0, self.sequences[seq1_idx].data.len(),
                seq2_idx, 0, self.sequences[seq2_idx].data.len(),
                true, scores, max_divergence
            ) {
                if score >= 0 {
                    alignments.push(Alignment {
                        seq1_idx,
                        seq2_idx,
                        seq1_start: 0,
                        seq1_end: self.sequences[seq1_idx].data.len(),
                        seq2_start: 0,
                        seq2_end: self.sequences[seq2_idx].data.len(),
                        seq2_is_rc: true,
                        cigar,
                        score,
                    });
                }
            }
        }
        
        Ok(alignments)
    }
    
    fn align_region(
        &self,
        seq1_idx: usize,
        seq1_start: usize,
        seq1_end: usize,
        seq2_idx: usize,
        seq2_start: usize,
        seq2_end: usize,
        seq2_is_rc: bool,
        scores: &AlignmentScores,
        max_divergence: Option<f64>,
    ) -> Result<(i32, String), Box<dyn std::error::Error>> {
        let seq1 = &self.sequences[seq1_idx].data[seq1_start..seq1_end];
        let seq2_region = &self.sequences[seq2_idx].data[seq2_start..seq2_end];
        
        let seq2 = if seq2_is_rc {
            reverse_complement(seq2_region)
        } else {
            seq2_region.to_vec()
        };
        
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
        
        let status = aligner.align(seq1, &seq2);
        
        match status {
            AlignmentStatus::Completed => {
                let score = aligner.score();
                
                if let Some(max_div) = max_divergence {
                    let max_allowed_score = scores.max_score_for_divergence(seq1.len().min(seq2.len()), max_div);
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
        alignment: &Alignment,
        min_match_length: usize,
        verbose: bool,
    ) {
        let seq1 = &self.sequences[alignment.seq1_idx];
        let seq2 = &self.sequences[alignment.seq2_idx];
        
        if verbose && (alignment.seq1_start > 0 || alignment.seq2_start > 0) {
            println!("  Processing patched alignment: {}[{}..{}] <-> {}[{}..{}] ({})",
                     seq1.id, alignment.seq1_start, alignment.seq1_end,
                     seq2.id, alignment.seq2_start, alignment.seq2_end,
                     if alignment.seq2_is_rc { "reverse" } else { "forward" });
        }
        
        // Parse CIGAR and process matches
        let mut operations = Vec::new();
        let mut count = 0;
        
        for ch in alignment.cigar.chars() {
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
        
        let mut pos1 = 0;
        let mut pos2 = 0;
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
                    if in_match_run && match_run_len >= min_match_length {
                        self.process_match_run(
                            alignment,
                            match_run_start1, match_run_start2, match_run_len,
                            verbose
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
        
        if in_match_run && match_run_len >= min_match_length {
            self.process_match_run(
                alignment,
                match_run_start1, match_run_start2, match_run_len,
                verbose
            );
        }
    }
    
    fn process_match_run(
        &mut self,
        alignment: &Alignment,
        start1: usize,
        start2: usize,
        length: usize,
        _verbose: bool,
    ) {
        let seq1 = &self.sequences[alignment.seq1_idx];
        let seq2 = &self.sequences[alignment.seq2_idx];
        
        for j in 0..length {
            let pos1_fwd = make_pos(seq1.offset + alignment.seq1_start + start1 + j, false);
            
            let pos2_oriented = if alignment.seq2_is_rc {
                let rev_offset = seq2.offset + (alignment.seq2_end - 1 - (start2 + j));
                make_pos(rev_offset, true)
            } else {
                make_pos(seq2.offset + alignment.seq2_start + start2 + j, false)
            };
            
            let idx1 = self.pos_to_uf_index(pos1_fwd);
            let idx2 = self.pos_to_uf_index(pos2_oriented);
            self.union_find.unite(idx1, idx2);
            
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
    
    fn build_bidirected_graph(&self, _verbose: bool) -> Result<BidirectedGraph, Box<dyn std::error::Error>> {
        let mut graph = BidirectedGraph::new();
        let mut union_to_node: HashMap<usize, usize> = HashMap::new();
        let mut next_node_id = 1;
        
        for seq in &self.sequences {
            let mut path = BiPath::new(seq.id.clone());
            
            for i in 0..seq.data.len() {
                let pos = make_pos(seq.offset + i, false);
                let uf_idx = self.pos_to_uf_index(pos);
                let union = self.union_find.find(uf_idx);
                
                let node_id = match union_to_node.get(&union) {
                    Some(&id) => id,
                    None => {
                        let id = next_node_id;
                        next_node_id += 1;
                        union_to_node.insert(union, id);
                        graph.add_node(id, vec![seq.data[i]]);
                        id
                    }
                };
                
                path.add_step(Handle::forward(node_id));
            }
            
            graph.paths.push(path);
        }
        
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

pub fn run_inversion_aware_seqrush(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    let sequences = read_fasta(&args.sequences)?;
    
    if args.verbose {
        println!("Loaded {} sequences", sequences.len());
        for seq in &sequences {
            println!("  {}: {} bp", seq.id, seq.data.len());
        }
    }
    
    let mut seqrush = InversionAwareSeqRush::new(sequences);
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