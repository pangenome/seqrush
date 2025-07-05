use clap::Parser;
use rayon::prelude::*;
use std::collections::{HashMap, BTreeMap};
use std::sync::atomic::{AtomicUsize, AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::fs::File;
use std::io::{BufWriter, Write, BufRead, BufReader};

/// Position with orientation encoding (like seqwish)
/// LSB = strand (0=forward, 1=reverse), remaining bits = offset
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct Position(u64);

impl Position {
    fn new(offset: usize, reverse: bool) -> Self {
        Position((offset as u64) << 1 | (reverse as u64))
    }
    
    fn offset(&self) -> usize {
        (self.0 >> 1) as usize
    }
    
    fn is_reverse(&self) -> bool {
        (self.0 & 1) != 0
    }
    
    fn increment(&mut self) {
        self.0 = ((self.offset() + 1) as u64) << 1 | (self.is_reverse() as u64);
    }
}

#[derive(Parser)]
#[command(name = "seqrush", version = "0.2.0")]
struct Args {
    #[arg(short, long)]
    sequences: String,
    
    #[arg(short, long, default_value = "output.gfa")]
    output: String,
    
    #[arg(short, long, default_value = "4")]
    threads: usize,
    
    #[arg(short = 'b', long, default_value = "1000000")]
    batch_size: usize,
    
    #[arg(short = 'k', long, default_value = "15")]
    min_match_length: usize,
    
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Clone)]
struct Sequence {
    id: String,
    data: Vec<u8>,
    offset: usize,
}

#[derive(Debug, Clone, Copy)]
struct Match {
    start: usize,    // Start in concatenated sequence
    end: usize,      // End in concatenated sequence  
    pos: Position,   // Position in other sequence
}

struct SeqRush {
    sequences: Vec<Sequence>,
    total_length: usize,
    matches: Vec<Match>,
    node_boundaries: Vec<bool>,  // Marks start of nodes
    seq_to_node: Vec<Option<usize>>,
    paths: Vec<Vec<usize>>,
}

impl SeqRush {
    fn new(sequences: Vec<Sequence>) -> Self {
        let total_length = sequences.iter().map(|s| s.data.len()).sum();
        Self {
            sequences,
            total_length,
            matches: Vec::new(),
            node_boundaries: vec![false; total_length + 1],
            seq_to_node: vec![None; total_length],
            paths: Vec::new(),
        }
    }
    
    fn build_graph(&mut self, args: &Args) {
        println!("Building graph with {} sequences (total length: {})", 
                 self.sequences.len(), self.total_length);
        
        // Phase 1: Collect matches (simplified - just all vs all alignments)
        self.collect_matches(args);
        
        // Phase 2: Build transitive closure in batches
        self.build_transitive_closure(args);
        
        // Phase 3: Extract paths
        self.extract_paths();
    }
    
    fn collect_matches(&mut self, args: &Args) {
        use lib_wfa2::affine_wavefront::{AffineWavefronts, MemoryMode};
        
        // Simple all-vs-all alignment for now
        for i in 0..self.sequences.len() {
            for j in i+1..self.sequences.len() {
                let seq1 = &self.sequences[i];
                let seq2 = &self.sequences[j];
                
                let mut wf = AffineWavefronts::with_penalties(
                    5, 8, 2,
                    MemoryMode::MemoryMed
                ).unwrap();
                
                let status = wf.align(&seq1.data, &seq2.data);
                if status.is_ok() {
                    // Extract matches from alignment
                    let cigar = wf.cigar_string();
                    self.extract_matches_from_cigar(
                        &cigar, 
                        seq1.offset, 
                        seq2.offset,
                        args.min_match_length
                    );
                }
            }
        }
        
        println!("Collected {} matches", self.matches.len());
    }
    
    fn extract_matches_from_cigar(&mut self, cigar: &str, offset1: usize, offset2: usize, min_len: usize) {
        let mut pos1 = 0;
        let mut pos2 = 0;
        let mut in_match = false;
        let mut match_start1 = 0;
        let mut match_start2 = 0;
        let mut match_len = 0;
        
        let mut count = 0;
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                count = count * 10 + (ch as usize - '0' as usize);
            } else {
                if count == 0 { count = 1; }
                
                match ch {
                    'M' | '=' => {
                        if !in_match {
                            in_match = true;
                            match_start1 = pos1;
                            match_start2 = pos2;
                            match_len = 0;
                        }
                        match_len += count;
                        pos1 += count;
                        pos2 += count;
                    }
                    _ => {
                        if in_match && match_len >= min_len {
                            // Add forward match
                            self.matches.push(Match {
                                start: offset1 + match_start1,
                                end: offset1 + match_start1 + match_len,
                                pos: Position::new(offset2 + match_start2, false),
                            });
                            // Add reverse match  
                            self.matches.push(Match {
                                start: offset2 + match_start2,
                                end: offset2 + match_start2 + match_len,
                                pos: Position::new(offset1 + match_start1, false),
                            });
                        }
                        in_match = false;
                        
                        match ch {
                            'X' => { pos1 += count; pos2 += count; }
                            'I' => { pos2 += count; }
                            'D' => { pos1 += count; }
                            _ => {}
                        }
                    }
                }
                count = 0;
            }
        }
        
        // Handle final match
        if in_match && match_len >= min_len {
            self.matches.push(Match {
                start: offset1 + match_start1,
                end: offset1 + match_start1 + match_len,
                pos: Position::new(offset2 + match_start2, false),
            });
            self.matches.push(Match {
                start: offset2 + match_start2,
                end: offset2 + match_start2 + match_len,
                pos: Position::new(offset1 + match_start1, false),
            });
        }
    }
    
    fn build_transitive_closure(&mut self, args: &Args) {
        use uf_rush::UFRush;
        
        let mut seen = vec![false; self.total_length];
        let mut seq_v_length = 0;
        
        // Process in batches
        let mut i = 0;
        while i < self.total_length {
            // Find next unseen position
            while i < self.total_length && seen[i] {
                i += 1;
            }
            if i >= self.total_length { break; }
            
            // Collect batch of unseen positions
            let batch_start = i;
            let mut batch_positions = Vec::new();
            let mut bases_in_batch = 0;
            
            while bases_in_batch < args.batch_size && i < self.total_length {
                if !seen[i] {
                    batch_positions.push(i);
                    bases_in_batch += 1;
                }
                i += 1;
            }
            
            if batch_positions.is_empty() { continue; }
            
            // Build union-find for this batch
            let mut position_map = HashMap::new();
            for (idx, &pos) in batch_positions.iter().enumerate() {
                position_map.insert(pos, idx);
            }
            
            let uf = Arc::new(UFRush::new(batch_positions.len()));
            
            // Process matches that overlap with batch
            for m in &self.matches {
                let mut united_any = false;
                for offset in 0..m.end.saturating_sub(m.start) {
                    let pos1 = m.start + offset;
                    let pos2 = m.pos.offset() + offset;
                    
                    if let (Some(&idx1), Some(&idx2)) = 
                        (position_map.get(&pos1), position_map.get(&pos2)) {
                        uf.unite(idx1, idx2);
                        united_any = true;
                    }
                }
            }
            
            // Extract components and create nodes
            let mut components: HashMap<usize, Vec<usize>> = HashMap::new();
            for (idx, &pos) in batch_positions.iter().enumerate() {
                let root = uf.find(idx);
                components.entry(root).or_insert_with(Vec::new).push(pos);
            }
            
            // Create nodes for each component
            for (_root, positions) in components {
                // Mark node boundary at first position
                if let Some(&first_pos) = positions.first() {
                    self.node_boundaries[first_pos] = true;
                }
                
                // Assign node ID to all positions
                for &pos in &positions {
                    self.seq_to_node[pos] = Some(seq_v_length);
                    seen[pos] = true;
                }
                
                seq_v_length += 1;
            }
        }
        
        println!("Created {} nodes", seq_v_length);
    }
    
    fn extract_paths(&mut self) {
        for seq in &self.sequences {
            let mut path = Vec::new();
            let mut last_node = None;
            
            for i in 0..seq.data.len() {
                let pos = seq.offset + i;
                if let Some(node_id) = self.seq_to_node[pos] {
                    // Only add if it's a node boundary or different from last
                    if self.node_boundaries[pos] || last_node != Some(node_id) {
                        path.push(node_id);
                        last_node = Some(node_id);
                    }
                }
            }
            
            self.paths.push(path);
        }
    }
    
    fn write_gfa(&self, output_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);
        
        writeln!(writer, "H\tVN:Z:1.0")?;
        
        // Collect unique nodes  
        let mut node_chars = BTreeMap::new();
        for i in 0..self.total_length {
            if self.node_boundaries[i] {
                if let Some(node_id) = self.seq_to_node[i] {
                    let seq_idx = self.sequences.iter()
                        .position(|s| i >= s.offset && i < s.offset + s.data.len())
                        .unwrap();
                    let char_idx = i - self.sequences[seq_idx].offset;
                    let ch = self.sequences[seq_idx].data[char_idx];
                    node_chars.entry(node_id).or_insert_with(Vec::new).push(ch);
                }
            }
        }
        
        // Write segments
        for (node_id, chars) in &node_chars {
            writeln!(writer, "S\t{}\t{}", 
                     node_id + 1,
                     String::from_utf8_lossy(&chars))?;
        }
        
        // Write paths
        for (seq_idx, path) in self.paths.iter().enumerate() {
            if !path.is_empty() {
                let path_str = path.iter()
                    .map(|&node_id| format!("{}+", node_id + 1))
                    .collect::<Vec<_>>()
                    .join(",");
                writeln!(writer, "P\t{}\t{}\t*", 
                         self.sequences[seq_idx].id, 
                         path_str)?;
            }
        }
        
        // Write edges
        let mut edges = std::collections::HashSet::new();
        for path in &self.paths {
            for window in path.windows(2) {
                if let [from, to] = window {
                    edges.insert((from + 1, to + 1));
                }
            }
        }
        
        for (from, to) in edges {
            writeln!(writer, "L\t{}\t+\t{}\t+\t0M", from, to)?;
        }
        
        Ok(())
    }
}

fn load_sequences(file_path: &str) -> Result<Vec<Sequence>, Box<dyn std::error::Error>> {
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();
    
    let sequences = load_sequences(&args.sequences)?;
    println!("Loaded {} sequences", sequences.len());
    
    let mut seqrush = SeqRush::new(sequences);
    seqrush.build_graph(&args);
    seqrush.write_gfa(&args.output)?;
    
    println!("Graph written to {}", args.output);
    Ok(())
}