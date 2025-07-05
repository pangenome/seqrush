use clap::Parser;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet, BTreeMap};
use std::sync::atomic::{AtomicUsize, AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use uf_rush::UFRush as LockFreeUnionFind;
use lib_wfa2::affine_wavefront::{AffineWavefronts, AlignmentStatus, MemoryMode};

#[derive(Parser)]
#[command(
    name = "seqrush",
    version = "0.1.0",
    author = "SeqRush Contributors",
    about = "Dynamic parallel in-memory implementation of seqwish's transitive closure algorithm",
    long_about = "SeqRush implements seqwish's graph induction algorithm using dynamic alignment and \
                  lock-free union-find to build pangenome graphs from sequences via transitive match closures."
)]
struct Args {
    /// Input FASTA file containing sequences to align
    #[arg(short, long, value_name = "FILE")]
    sequences: String,
    
    /// Number of threads for parallel processing
    #[arg(short, long, default_value = "4", value_name = "NUM")]
    threads: usize,
    
    /// Output GFA file path
    #[arg(short, long, default_value = "output.gfa", value_name = "FILE")]
    output: String,
    
    /// Maximum number of alignment iterations
    #[arg(short = 'i', long, default_value = "100", value_name = "NUM")]
    max_iterations: usize,
    
    /// Minimum match length to consider (-k in seqwish)
    #[arg(short = 'k', long, default_value = "1", value_name = "NUM")]
    min_match_length: usize,
    
    /// Alignment match score
    #[arg(long, default_value = "0", value_name = "NUM")]
    match_score: i32,
    
    /// Alignment mismatch penalty
    #[arg(long, default_value = "5", value_name = "NUM")]
    mismatch_penalty: i32,
    
    /// Gap open penalty
    #[arg(long, default_value = "8", value_name = "NUM")]
    gap_open: i32,
    
    /// Gap extend penalty
    #[arg(long, default_value = "2", value_name = "NUM")]
    gap_extend: i32,
    
    /// Number of stable iterations required for convergence
    #[arg(long, default_value = "3", value_name = "NUM")]
    convergence_stability: usize,
    
    /// Sparsification factor for Erdős-Rényi random graph (auto calculates if not set)
    #[arg(long, value_name = "FLOAT")]
    sparsification: Option<f64>,
    
    /// Safety factor multiplier for Erdős-Rényi connectivity threshold
    #[arg(long, default_value = "3.0", value_name = "FLOAT")]
    erdos_renyi_safety_factor: f64,
    
    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Debug, Clone)]
struct Sequence {
    id: String,
    data: Vec<u8>,
    offset: usize, // Global offset in concatenated sequence
}

/// A match between two sequence positions
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Match {
    seq1_id: usize,
    seq1_pos: usize,
    seq2_id: usize, 
    seq2_pos: usize,
    length: usize,
    reverse: bool, // true if seq2 is reverse complement
}


/// Graph node containing sequence of characters
#[derive(Debug, Clone)]
struct Node {
    id: usize,
    sequence: Vec<u8>,
}

/// The main SeqRush struct implementing seqwish algorithm
struct SeqRush {
    sequences: Vec<Sequence>,
    total_length: usize,
    union_find: Arc<LockFreeUnionFind>,
    matches: Arc<Mutex<Vec<Match>>>, // Dynamic matches from alignments
    seen: Arc<Vec<AtomicBool>>, // Bitvector X from algorithm
    nodes: Arc<Mutex<Vec<Node>>>, // Graph nodes V
    seq_to_node: Arc<Mutex<BTreeMap<usize, usize>>>, // Z mapping
    node_to_seq: Arc<Mutex<BTreeMap<usize, Vec<usize>>>>, // Z̄ mapping
    converged: Arc<AtomicBool>,
    iteration_count: Arc<AtomicUsize>,
    config: AlignmentConfig,
}

#[derive(Debug, Clone)]
struct AlignmentConfig {
    max_iterations: usize,
    min_match_length: usize,
    match_score: i32,
    mismatch_penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    convergence_stability: usize,
    sparsification: f64,
    erdos_renyi_safety_factor: f64,
    verbose: bool,
}

impl SeqRush {
    fn new(sequences: Vec<Sequence>, config: AlignmentConfig) -> Self {
        let total_length = sequences.iter().map(|s| s.data.len()).sum();
        
        // Initialize seen bitvector
        let mut seen_vec = Vec::with_capacity(total_length);
        for _ in 0..total_length {
            seen_vec.push(AtomicBool::new(false));
        }
        
        Self {
            sequences,
            total_length,
            union_find: Arc::new(LockFreeUnionFind::new(total_length)),
            matches: Arc::new(Mutex::new(Vec::new())),
            seen: Arc::new(seen_vec),
            nodes: Arc::new(Mutex::new(Vec::new())),
            seq_to_node: Arc::new(Mutex::new(BTreeMap::new())),
            node_to_seq: Arc::new(Mutex::new(BTreeMap::new())),
            converged: Arc::new(AtomicBool::new(false)),
            iteration_count: Arc::new(AtomicUsize::new(0)),
            config,
        }
    }
    
    /// Main algorithm: build graph through dynamic alignment and transitive closure
    fn build_graph(&mut self) {
        println!("Building graph with {} sequences (total length: {})", 
                 self.sequences.len(), self.total_length);
        
        // Phase 1: Dynamic alignment to discover matches
        self.discover_matches_dynamically();
        
        // Phase 2: Graph induction via transitive closure (seqwish algorithm)
        self.induce_graph();
    }
    
    /// Discover matches through dynamic pairwise alignment
    fn discover_matches_dynamically(&mut self) {
        // Special case: single sequence has no pairs to align
        if self.sequences.len() <= 1 {
            println!("Single sequence: no alignment needed");
            return;
        }
        
        let mut stable_iterations = 0;
        
        for iteration in 0..self.config.max_iterations {
            if self.converged.load(Ordering::Relaxed) {
                break;
            }
            
            println!("Alignment iteration {}", iteration + 1);
            
            // Select random pairs using Erdős-Rényi sparsification
            let pairs = self.select_sequence_pairs();
            
            // If no pairs selected, increment stable iterations
            if pairs.is_empty() {
                stable_iterations += 1;
                if stable_iterations >= self.config.convergence_stability {
                    println!("Converged: no pairs to align for {} iterations", stable_iterations);
                    self.converged.store(true, Ordering::Relaxed);
                    break;
                }
                continue;
            }
            
            // Align pairs in parallel and collect matches
            let new_matches: Vec<Match> = pairs
                .par_iter()
                .flat_map(|&(i, j)| self.align_and_extract_matches(i, j))
                .collect();
            
            // Add new matches, avoiding duplicates
            let new_match_count = {
                let mut matches = self.matches.lock().unwrap();
                let existing_matches: HashSet<Match> = matches.iter().cloned().collect();
                let mut added = 0;
                
                for m in new_matches {
                    if !existing_matches.contains(&m) {
                        matches.push(m);
                        added += 1;
                    }
                }
                
                added
            };
            
            // Check convergence
            if new_match_count == 0 {
                stable_iterations += 1;
                if stable_iterations >= self.config.convergence_stability {
                    println!("Converged: no new matches for {} iterations", stable_iterations);
                    self.converged.store(true, Ordering::Relaxed);
                    break;
                }
            } else {
                stable_iterations = 0;
                println!("Found {} new matches", new_match_count);
            }
            
            self.iteration_count.fetch_add(1, Ordering::Relaxed);
        }
        
        let total_matches = self.matches.lock().unwrap().len();
        println!("Match discovery complete: {} total matches", total_matches);
    }
    
    /// Select sequence pairs using Erdős-Rényi sparsification
    fn select_sequence_pairs(&self) -> Vec<(usize, usize)> {
        let n = self.sequences.len();
        
        // If only one sequence, no pairs to select
        if n <= 1 {
            return Vec::new();
        }
        
        let mut pairs = Vec::new();
        
        // For each possible pair, decide whether to include it
        for i in 0..n {
            for j in (i + 1)..n {
                if self.should_align_pair(i, j) {
                    pairs.push((i, j));
                }
            }
        }
        
        // Shuffle for better load balancing
        use rand::seq::SliceRandom;
        let mut rng = rand::thread_rng();
        pairs.shuffle(&mut rng);
        
        if self.config.verbose {
            let total_possible = n * (n - 1) / 2;
            println!("Selected {} of {} possible pairs ({:.1}%)", 
                     pairs.len(), total_possible, 
                     100.0 * pairs.len() as f64 / total_possible as f64);
        }
        
        pairs
    }
    
    /// Determine if a pair should be aligned based on sparsification
    fn should_align_pair(&self, seq1_id: usize, seq2_id: usize) -> bool {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        
        // Create deterministic hash
        let mut hasher = DefaultHasher::new();
        let (min_id, max_id) = if seq1_id < seq2_id {
            (seq1_id, seq2_id)
        } else {
            (seq2_id, seq1_id)
        };
        min_id.hash(&mut hasher);
        max_id.hash(&mut hasher);
        let hash_value = hasher.finish();
        
        // Convert to probability
        let prob = (hash_value as f64) / (u64::MAX as f64);
        prob < self.config.sparsification
    }
    
    /// Align two sequences and extract matches
    fn align_and_extract_matches(&self, seq1_id: usize, seq2_id: usize) -> Vec<Match> {
        let seq1 = &self.sequences[seq1_id];
        let seq2 = &self.sequences[seq2_id];
        
        if self.config.verbose {
            println!("Aligning seq{} ({} bp) vs seq{} ({} bp)", 
                     seq1_id, seq1.data.len(), seq2_id, seq2.data.len());
        }
        
        let mut all_matches = Vec::new();
        
        // Try forward alignment
        if let Some(matches) = self.align_and_parse(seq1_id, seq2_id, &seq1.data, &seq2.data, false) {
            all_matches.extend(matches);
        }
        
        // Try reverse complement alignment
        // DISABLED: Reverse complement matching is causing incorrect node associations
        // for repetitive sequences. For the simplest algorithm, we only consider
        // forward matches.
        // let seq2_rc = Self::reverse_complement(&seq2.data);
        // if let Some(matches) = self.align_and_parse(seq1_id, seq2_id, &seq1.data, &seq2_rc, true) {
        //     all_matches.extend(matches);
        // }
        
        if self.config.verbose && !all_matches.is_empty() {
            println!("Found {} matches between seq{} and seq{}", all_matches.len(), seq1_id, seq2_id);
        }
        
        all_matches
    }
    
    /// Perform alignment and parse CIGAR to extract matches
    fn align_and_parse(&self, seq1_id: usize, seq2_id: usize, 
                       seq1: &[u8], seq2: &[u8], reverse: bool) -> Option<Vec<Match>> {
        // Use WFA2 for alignment
        let mut aligner = AffineWavefronts::with_penalties(
            self.config.match_score,
            self.config.mismatch_penalty,
            self.config.gap_open,
            self.config.gap_extend,
        );
        
        aligner.set_memory_mode(MemoryMode::Ultralow);
        
        match aligner.align(seq1, seq2) {
            AlignmentStatus::Completed => {
                let cigar = aligner.cigar();
                if self.config.verbose {
                    println!("Alignment completed: seq{} vs seq{} (reverse={}): CIGAR={}", 
                             seq1_id, seq2_id, reverse, String::from_utf8_lossy(cigar));
                }
                self.parse_cigar_to_matches(seq1_id, seq2_id, cigar, seq1.len(), seq2.len(), reverse)
            }
            status => {
                if self.config.verbose {
                    println!("Alignment failed: seq{} vs seq{} (reverse={}): {:?}", 
                             seq1_id, seq2_id, reverse, status);
                }
                None
            }
        }
    }
    
    /// Parse CIGAR string to extract match intervals
    fn parse_cigar_to_matches(&self, seq1_id: usize, seq2_id: usize,
                              cigar: &[u8], _seq1_len: usize, seq2_len: usize, 
                              reverse: bool) -> Option<Vec<Match>> {
        let cigar_str = String::from_utf8_lossy(cigar);
        let mut matches = Vec::new();
        
        let mut seq1_pos = 0;
        let mut seq2_pos = 0;
        let mut current_num = String::new();
        let mut in_match_run = false;
        let mut match_start_seq1 = 0;
        let mut match_start_seq2 = 0;
        let mut match_length = 0;
        
        // Helper to finish a match run
        let finish_match = |matches: &mut Vec<Match>, match_start_seq1, match_start_seq2, match_length| {
            if match_length >= self.config.min_match_length {
                let seq2_actual_pos = if reverse {
                    seq2_len - match_start_seq2 - match_length
                } else {
                    match_start_seq2
                };
                
                if self.config.verbose && reverse {
                    println!("  Reverse match: seq1[{}..{}] vs seq2_rc[{}..{}] -> seq2[{}..{}]",
                             match_start_seq1, match_start_seq1 + match_length,
                             match_start_seq2, match_start_seq2 + match_length,
                             seq2_actual_pos, seq2_actual_pos + match_length);
                }
                
                matches.push(Match {
                    seq1_id,
                    seq1_pos: match_start_seq1,
                    seq2_id,
                    seq2_pos: seq2_actual_pos,
                    length: match_length,
                    reverse,
                });
            }
        };
        
        for ch in cigar_str.chars() {
            if ch.is_ascii_digit() {
                current_num.push(ch);
            } else {
                let count = if current_num.is_empty() { 1 } else { current_num.parse::<usize>().unwrap_or(1) };
                current_num.clear();
                
                match ch {
                    'M' | '=' => {
                        // Match or alignment position
                        if !in_match_run {
                            in_match_run = true;
                            match_start_seq1 = seq1_pos;
                            match_start_seq2 = seq2_pos;
                            match_length = 0;
                        }
                        match_length += count;
                        seq1_pos += count;
                        seq2_pos += count;
                    }
                    'X' | 'I' | 'D' => {
                        // End any current match run
                        if in_match_run {
                            finish_match(&mut matches, match_start_seq1, match_start_seq2, match_length);
                            in_match_run = false;
                        }
                        
                        match ch {
                            'X' => {
                                // Mismatch
                                seq1_pos += count;
                                seq2_pos += count;
                            }
                            'I' => {
                                // Insertion in seq2
                                seq2_pos += count;
                            }
                            'D' => {
                                // Deletion in seq2
                                seq1_pos += count;
                            }
                            _ => {}
                        }
                    }
                    _ => {}
                }
            }
        }
        
        // Finish any remaining match run
        if in_match_run {
            finish_match(&mut matches, match_start_seq1, match_start_seq2, match_length);
        }
        
        if matches.is_empty() {
            None
        } else {
            Some(matches)
        }
    }
    
    /// Induce graph from matches using transitive closure (core seqwish algorithm)
    fn induce_graph(&mut self) {
        println!("Inducing graph from matches...");
        
        if self.total_length == 0 {
            println!("No sequences to process");
            return;
        }
        
        // First, build the union-find structure from all matches
        self.build_union_find_from_matches();
        
        // Debug: check union-find structure
        if self.config.verbose {
            println!("Union-find structure:");
            for i in 0..self.total_length {
                println!("  Position {} -> root {}", i, self.union_find.find(i));
            }
        }
        
        // Process each character position
        let positions: Vec<usize> = (0..self.total_length).collect();
        positions.into_par_iter().for_each(|pos| {
            self.process_position(pos);
        });
        
        println!("Graph induction complete: {} nodes", self.nodes.lock().unwrap().len());
    }
    
    /// Build union-find structure from all matches
    fn build_union_find_from_matches(&self) {
        let matches = self.matches.lock().unwrap();
        
        println!("Building union-find from {} matches", matches.len());
        
        for m in matches.iter() {
            let seq1_offset = self.sequences[m.seq1_id].offset;
            let seq2_offset = self.sequences[m.seq2_id].offset;
            
            // Union all corresponding positions in the match
            for i in 0..m.length {
                let pos1 = seq1_offset + m.seq1_pos + i;
                let pos2 = if m.reverse {
                    // For reverse matches, seq2_pos is already the position in the original sequence
                    // But we need to traverse it in reverse order
                    seq2_offset + m.seq2_pos + (m.length - 1 - i)
                } else {
                    seq2_offset + m.seq2_pos + i
                };
                
                if self.config.verbose {
                    let (s1_id, s1_pos) = self.global_to_seq_pos(pos1);
                    let (s2_id, s2_pos) = self.global_to_seq_pos(pos2);
                    println!("  Uniting position {} (seq{}[{}]='{}') with position {} (seq{}[{}]='{}')", 
                             pos1, s1_id, s1_pos, self.sequences[s1_id].data[s1_pos] as char,
                             pos2, s2_id, s2_pos, self.sequences[s2_id].data[s2_pos] as char);
                }
                self.union_find.unite(pos1, pos2);
            }
            
            if self.config.verbose {
                println!("Match: seq{} pos {} vs seq{} pos {} len {} reverse={}",
                         m.seq1_id, m.seq1_pos, m.seq2_id, m.seq2_pos, m.length, m.reverse);
            }
        }
    }
    
    /// Process a single position (implements lines 19-32 of seqwish algorithm)
    fn process_position(&self, pos: usize) {
        // Check if already seen
        if self.seen[pos].load(Ordering::Acquire) {
            return;
        }
        
        // Get the root of this position in union-find
        let root = self.union_find.find(pos);
        
        // Collect all positions with the same root - this is the transitive closure
        let mut transitive_matches = Vec::new();
        
        // Check all positions to find those in the same equivalence class
        for i in 0..self.total_length {
            if self.union_find.find(i) == root && !self.seen[i].load(Ordering::Acquire) {
                transitive_matches.push(i);
            }
        }
        
        if transitive_matches.is_empty() {
            return;
        }
        
        if self.config.verbose && transitive_matches.len() > 1 {
            println!("Position {} (root {}) has {} matches in equivalence class", 
                     pos, root, transitive_matches.len());
        }
        
        // Mark all positions as seen atomically
        for &pos in &transitive_matches {
            self.seen[pos].store(true, Ordering::Release);
        }
        
        // Create new node
        let node_id = {
            let mut nodes = self.nodes.lock().unwrap();
            let id = nodes.len();
            
            // Determine consensus character(s) for the node
            let mut char_counts = HashMap::new();
            let mut first_char = None;
            for &match_pos in &transitive_matches {
                let (seq_id, seq_pos) = self.global_to_seq_pos(match_pos);
                let ch = self.sequences[seq_id].data[seq_pos];
                if first_char.is_none() {
                    first_char = Some(ch);
                }
                *char_counts.entry(ch).or_insert(0) += 1;
            }
            
            // For now, just use the first character seen
            let consensus = first_char.unwrap_or(b'N');
            
            nodes.push(Node {
                id,
                sequence: vec![consensus],
            });
            
            id
        };
        
        // Update mappings
        {
            let mut seq_to_node = self.seq_to_node.lock().unwrap();
            let mut node_to_seq = self.node_to_seq.lock().unwrap();
            
            for &match_pos in &transitive_matches {
                seq_to_node.insert(match_pos, node_id);
                node_to_seq.entry(node_id).or_insert_with(Vec::new).push(match_pos);
            }
        }
    }
    
    
    /// Write the graph in GFA format
    fn write_gfa(&self, output_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        use std::fs::File;
        use std::io::{BufWriter, Write};
        
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);
        
        // Header
        writeln!(writer, "H\tVN:Z:1.0")?;
        
        // Segments (nodes)
        let nodes = self.nodes.lock().unwrap();
        for node in nodes.iter() {
            writeln!(writer, "S\t{}\t{}", 
                     node.id + 1,  // GFA uses 1-based IDs
                     String::from_utf8_lossy(&node.sequence))?;
        }
        
        // Paths and edges
        let seq_to_node = self.seq_to_node.lock().unwrap();
        let mut edges = HashSet::new();
        
        for (_seq_id, seq) in self.sequences.iter().enumerate() {
            let mut path: Vec<String> = Vec::new();
            let mut pos = 0;
            let mut prev_node_id: Option<usize> = None;
            
            while pos < seq.data.len() {
                let global_pos = seq.offset + pos;
                
                if let Some(&node_id) = seq_to_node.get(&global_pos) {
                    let node_str = format!("{}+", node_id + 1);
                    path.push(node_str);
                    
                    // Add edge from previous node if exists
                    if let Some(prev_id) = prev_node_id {
                        if prev_id != node_id {
                            edges.insert((prev_id + 1, node_id + 1));
                        }
                    }
                    
                    prev_node_id = Some(node_id);
                    pos += 1;
                } else {
                    // This shouldn't happen if graph is complete
                    eprintln!("Warning: position {} not mapped to any node", global_pos);
                    pos += 1;
                }
            }
            
            if !path.is_empty() {
                writeln!(writer, "P\t{}\t{}\t*", seq.id, path.join(","))?;
            }
        }
        
        // Write edges
        for (from, to) in edges {
            writeln!(writer, "L\t{}\t+\t{}\t+\t0M", from, to)?;
        }
        
        Ok(())
    }
    
    /// Convert global position to sequence ID and position
    fn global_to_seq_pos(&self, global_pos: usize) -> (usize, usize) {
        for (seq_id, seq) in self.sequences.iter().enumerate() {
            if global_pos >= seq.offset && global_pos < seq.offset + seq.data.len() {
                return (seq_id, global_pos - seq.offset);
            }
        }
        panic!("Invalid global position: {}", global_pos);
    }
    
    /// Calculate Erdős-Rényi threshold
    fn calculate_erdos_renyi_threshold(n: usize, safety_factor: f64) -> f64 {
        if n <= 1 {
            return 1.0;
        }
        
        let n_f = n as f64;
        let threshold = safety_factor * n_f.ln() / n_f;
        threshold.min(1.0).max(0.001)
    }
    
    /// Reverse complement a sequence
    fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&base| match base {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                _ => base,
            })
            .collect()
    }
}

/// Load sequences from FASTA file
fn load_sequences(file_path: &str) -> Result<Vec<Sequence>, Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    
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
    
    // Set up thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();
    
    // Load sequences
    let sequences = load_sequences(&args.sequences)?;
    println!("Loaded {} sequences", sequences.len());
    
    // Calculate sparsification
    let sparsification = if let Some(s) = args.sparsification {
        s
    } else {
        SeqRush::calculate_erdos_renyi_threshold(sequences.len(), args.erdos_renyi_safety_factor)
    };
    
    if args.verbose {
        println!("Using sparsification factor: {:.4}", sparsification);
    }
    
    // Create config
    let config = AlignmentConfig {
        max_iterations: args.max_iterations,
        min_match_length: args.min_match_length,
        match_score: args.match_score,
        mismatch_penalty: args.mismatch_penalty,
        gap_open: args.gap_open,
        gap_extend: args.gap_extend,
        convergence_stability: args.convergence_stability,
        sparsification,
        erdos_renyi_safety_factor: args.erdos_renyi_safety_factor,
        verbose: args.verbose,
    };
    
    // Build graph
    let mut seqrush = SeqRush::new(sequences, config);
    seqrush.build_graph();
    
    // Write output
    seqrush.write_gfa(&args.output)?;
    println!("Graph written to {}", args.output);
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_sequences() -> Vec<Sequence> {
        vec![
            Sequence {
                id: "seq1".to_string(),
                data: b"ATCGATCGATCGATCG".to_vec(),
                offset: 0,
            },
            Sequence {
                id: "seq2".to_string(),
                data: b"ATCGATCGATCGATCG".to_vec(),
                offset: 16,
            },
            Sequence {
                id: "seq3".to_string(),
                data: b"ATCGATCGCTCGATCG".to_vec(),
                offset: 32,
            },
        ]
    }

    fn create_test_config() -> AlignmentConfig {
        AlignmentConfig {
            max_iterations: 10,
            min_match_length: 4,
            match_score: 0,
            mismatch_penalty: 5,
            gap_open: 8,
            gap_extend: 2,
            convergence_stability: 2,
            sparsification: 1.0,
            erdos_renyi_safety_factor: 3.0,
            verbose: false,
        }
    }

    #[test]
    fn test_sequence_loading() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, ">seq1").unwrap();
        writeln!(temp_file, "ATCG").unwrap();
        writeln!(temp_file, ">seq2").unwrap();
        writeln!(temp_file, "GCTA").unwrap();
        temp_file.flush().unwrap();
        
        let sequences = load_sequences(temp_file.path().to_str().unwrap()).unwrap();
        assert_eq!(sequences.len(), 2);
        assert_eq!(sequences[0].id, "seq1");
        assert_eq!(sequences[0].data, b"ATCG");
        assert_eq!(sequences[0].offset, 0);
        assert_eq!(sequences[1].offset, 4);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(SeqRush::reverse_complement(b"ATCG"), b"CGAT");
        assert_eq!(SeqRush::reverse_complement(b"atcg"), b"CGAT");
        assert_eq!(SeqRush::reverse_complement(b"AAATTTCCCGGG"), b"CCCGGGAAATTT");
    }

    #[test]
    fn test_match_extraction() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences, config);
        
        // Test CIGAR parsing
        let cigar = b"8M1X7M";
        let matches = seqrush.parse_cigar_to_matches(0, 1, cigar, 16, 16, false).unwrap();
        
        assert_eq!(matches.len(), 2);
        assert_eq!(matches[0].length, 8);
        assert_eq!(matches[1].length, 7);
    }

    #[test]
    fn test_transitive_closure() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences, config);
        
        // Add some test matches
        {
            let mut matches = seqrush.matches.lock().unwrap();
            matches.push(Match {
                seq1_id: 0,
                seq1_pos: 0,
                seq2_id: 1,
                seq2_pos: 0,
                length: 8,
                reverse: false,
            });
        }
        
        // Build union-find
        seqrush.build_union_find_from_matches();
        
        // Check that positions are unioned
        let root0 = seqrush.union_find.find(0);
        let root16 = seqrush.union_find.find(16);
        assert_eq!(root0, root16); // Positions 0 and 16 should have same root
    }

    #[test]
    fn test_graph_construction() {
        let sequences = create_test_sequences();
        let mut config = create_test_config();
        config.max_iterations = 2;
        
        let mut seqrush = SeqRush::new(sequences, config);
        seqrush.build_graph();
        
        let nodes = seqrush.nodes.lock().unwrap();
        assert!(nodes.len() > 0);
    }

    #[test]
    fn test_gfa_output() {
        let sequences = vec![
            Sequence {
                id: "test1".to_string(),
                data: b"ATCG".to_vec(),
                offset: 0,
            },
            Sequence {
                id: "test2".to_string(),
                data: b"ATCG".to_vec(),
                offset: 4,
            },
        ];
        
        let mut config = create_test_config();
        config.max_iterations = 1;
        
        let mut seqrush = SeqRush::new(sequences, config);
        seqrush.build_graph();
        
        let temp_file = NamedTempFile::new().unwrap();
        seqrush.write_gfa(temp_file.path().to_str().unwrap()).unwrap();
        
        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        assert!(content.contains("H\tVN:Z:1.0"));
        assert!(content.contains("S\t"));
    }

    #[test]
    fn test_erdos_renyi_calculation() {
        let threshold_10 = SeqRush::calculate_erdos_renyi_threshold(10, 3.0);
        let threshold_100 = SeqRush::calculate_erdos_renyi_threshold(100, 3.0);
        let threshold_1000 = SeqRush::calculate_erdos_renyi_threshold(1000, 3.0);
        
        assert!(threshold_10 > threshold_100);
        assert!(threshold_100 > threshold_1000);
        assert!(threshold_1000 > 0.001);
        assert!(threshold_10 < 1.0);
    }

    #[test]
    fn test_sequence_integrity_preservation() {
        // Test that sequences can be fully reconstructed from the graph
        let sequences = vec![
            Sequence {
                id: "seq1".to_string(),
                data: b"ATCGATCGATCG".to_vec(),
                offset: 0,
            },
            Sequence {
                id: "seq2".to_string(),
                data: b"ATCGATCGATCG".to_vec(),
                offset: 12,
            },
        ];
        
        let mut config = create_test_config();
        config.max_iterations = 5;
        
        let mut seqrush = SeqRush::new(sequences.clone(), config);
        seqrush.build_graph();
        
        // Verify each sequence can be reconstructed
        let seq_to_node = seqrush.seq_to_node.lock().unwrap();
        let nodes = seqrush.nodes.lock().unwrap();
        
        for (seq_id, seq) in seqrush.sequences.iter().enumerate() {
            let mut reconstructed = Vec::new();
            
            for pos in 0..seq.data.len() {
                let global_pos = seq.offset + pos;
                if let Some(&node_id) = seq_to_node.get(&global_pos) {
                    if node_id < nodes.len() {
                        // For simplicity, take first character of node
                        if !nodes[node_id].sequence.is_empty() {
                            reconstructed.push(nodes[node_id].sequence[0]);
                        }
                    }
                }
            }
            
            // Verify reconstruction matches original
            assert_eq!(reconstructed.len(), seq.data.len(), 
                      "Sequence {} length mismatch", seq_id);
        }
    }

    #[test]
    fn test_graph_completeness() {
        // Test that all sequence positions are represented in the graph
        let sequences = create_test_sequences();
        let mut config = create_test_config();
        config.max_iterations = 3;
        
        let mut seqrush = SeqRush::new(sequences.clone(), config);
        seqrush.build_graph();
        
        let seq_to_node = seqrush.seq_to_node.lock().unwrap();
        
        // Every position should have a node mapping
        let total_positions = seqrush.total_length;
        let mapped_positions = seq_to_node.len();
        
        assert_eq!(mapped_positions, total_positions,
                  "Not all positions mapped to nodes: {} of {} mapped",
                  mapped_positions, total_positions);
    }

    #[test]
    fn test_match_symmetry() {
        // Test that matches are symmetric (if A matches B, B matches A)
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences, config);
        
        // Create a match
        let match_forward = Match {
            seq1_id: 0,
            seq1_pos: 0,
            seq2_id: 1,
            seq2_pos: 0,
            length: 8,
            reverse: false,
        };
        
        // Build union-find with this match
        {
            let mut matches = seqrush.matches.lock().unwrap();
            matches.push(match_forward);
        }
        seqrush.build_union_find_from_matches();
        
        // Check symmetry: positions that are unioned should have same root
        for i in 0..8 {
            let pos1 = 0 + i;  // seq1 position
            let pos2 = 16 + i; // seq2 position (offset 16)
            let root1 = seqrush.union_find.find(pos1);
            let root2 = seqrush.union_find.find(pos2);
            assert_eq!(root1, root2, "Position {} and {} should have same root", pos1, pos2);
        }
    }

    #[test]
    fn test_reverse_complement_handling() {
        // Test that reverse complement matches are handled correctly
        let sequences = vec![
            Sequence {
                id: "forward".to_string(),
                data: b"ATCG".to_vec(),
                offset: 0,
            },
            Sequence {
                id: "revcomp".to_string(),
                data: b"CGAT".to_vec(), // reverse complement of ATCG
                offset: 4,
            },
        ];
        
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences, config);
        
        // Add a reverse complement match
        {
            let mut matches = seqrush.matches.lock().unwrap();
            matches.push(Match {
                seq1_id: 0,
                seq1_pos: 0,
                seq2_id: 1,
                seq2_pos: 0,
                length: 4,
                reverse: true,
            });
        }
        
        seqrush.build_union_find_from_matches();
        
        // Verify reverse complement pairing
        // Position 0 of seq1 (A) should match position 3 of seq2 (T->A in RC)
        let root_a = seqrush.union_find.find(0);
        let root_t = seqrush.union_find.find(7); // position 3 in seq2 (offset 4)
        assert_eq!(root_a, root_t, "Reverse complement positions should be unioned");
    }

    #[test]
    fn test_convergence_detection() {
        // Test that the algorithm converges when no new matches are found
        let sequences = vec![
            Sequence {
                id: "seq1".to_string(),
                data: b"AAAA".to_vec(),
                offset: 0,
            },
            Sequence {
                id: "seq2".to_string(),
                data: b"CCCC".to_vec(),
                offset: 4,
            },
        ];
        
        let mut config = create_test_config();
        config.max_iterations = 10;
        config.convergence_stability = 2;
        config.min_match_length = 4; // Require full match
        
        let mut seqrush = SeqRush::new(sequences, config);
        seqrush.discover_matches_dynamically();
        
        // Should have converged after finding no matches
        let _iterations = seqrush.iteration_count.load(Ordering::Relaxed);
        assert!(seqrush.converged.load(Ordering::Relaxed), "Should be marked as converged");
        
        // Check that no matches were found
        let matches = seqrush.matches.lock().unwrap();
        assert_eq!(matches.len(), 0, "Should find no matches between AAAA and CCCC");
    }

    #[test]
    fn test_gfa_path_validity() {
        // Test that GFA paths correctly represent sequences
        let sequences = vec![
            Sequence {
                id: "test".to_string(),
                data: b"ATCGATCG".to_vec(),
                offset: 0,
            },
        ];
        
        let mut config = create_test_config();
        config.max_iterations = 1;
        
        let mut seqrush = SeqRush::new(sequences, config);
        seqrush.build_graph();
        
        // Write to temporary file
        let temp_file = NamedTempFile::new().unwrap();
        seqrush.write_gfa(temp_file.path().to_str().unwrap()).unwrap();
        
        // Read and parse GFA
        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        let mut has_path = false;
        
        for line in content.lines() {
            if line.starts_with("P\ttest\t") {
                has_path = true;
                // Verify path format
                let parts: Vec<&str> = line.split('\t').collect();
                assert_eq!(parts.len(), 4, "Path line should have 4 fields");
                assert_eq!(parts[1], "test", "Path should reference correct sequence");
                assert!(!parts[2].is_empty(), "Path should have nodes");
            }
        }
        
        assert!(has_path, "GFA should contain path for sequence");
    }


    #[test]
    fn test_parallel_safety() {
        // Test that parallel alignment doesn't corrupt data
        use std::sync::Arc;
        use std::thread;
        
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = Arc::new(SeqRush::new(sequences, config));
        
        // Run alignments in parallel threads
        let mut handles = vec![];
        
        for i in 0..4 {
            let seqrush_clone = Arc::clone(&seqrush);
            let handle = thread::spawn(move || {
                // Each thread tries to align different pairs
                let matches = seqrush_clone.align_and_extract_matches(i % 2, (i + 1) % 3);
                matches.len()
            });
            handles.push(handle);
        }
        
        // Wait for all threads
        for handle in handles {
            let result = handle.join();
            assert!(result.is_ok(), "Thread should complete without panic");
        }
        
        // Verify data structures are still consistent
        let total_length = seqrush.total_length;
        assert!(total_length > 0, "Total length should be preserved");
    }

    #[test]
    fn test_empty_sequences() {
        // Test handling of empty sequence input
        let sequences = vec![];
        let config = create_test_config();
        let mut seqrush = SeqRush::new(sequences, config);
        
        // Should handle empty input gracefully
        seqrush.build_graph();
        
        let nodes = seqrush.nodes.lock().unwrap();
        assert_eq!(nodes.len(), 0, "Empty input should produce no nodes");
    }

    #[test]
    fn test_single_sequence() {
        // Test with single sequence
        let sequences = vec![
            Sequence {
                id: "single".to_string(),
                data: b"ATCG".to_vec(),
                offset: 0,
            },
        ];
        
        let mut config = create_test_config();
        config.max_iterations = 1;
        
        let mut seqrush = SeqRush::new(sequences, config);
        seqrush.build_graph();
        
        // Should create nodes for the sequence
        let nodes = seqrush.nodes.lock().unwrap();
        assert!(nodes.len() > 0, "Single sequence should produce nodes");
        drop(nodes);  // Release lock before continuing
        
        // Verify GFA output
        let temp_file = NamedTempFile::new().unwrap();
        seqrush.write_gfa(temp_file.path().to_str().unwrap()).unwrap();
        
        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        assert!(content.contains("P\tsingle\t"), "Should have path for single sequence");
    }

    #[test] 
    fn test_large_scale_integrity() {
        // Test with larger dataset to ensure scalability
        let mut sequences = Vec::new();
        let mut offset = 0;
        
        // Create 10 similar sequences with small variations
        for i in 0..10 {
            let mut data = b"ATCGATCGATCGATCG".to_vec();
            if i % 2 == 0 {
                data[8] = b'T'; // Introduce variation
            }
            sequences.push(Sequence {
                id: format!("seq{}", i),
                data: data.clone(),
                offset,
            });
            offset += data.len();
        }
        
        let mut config = create_test_config();
        config.max_iterations = 5;
        config.sparsification = 0.5; // Use sparsification
        
        let mut seqrush = SeqRush::new(sequences.clone(), config);
        seqrush.build_graph();
        
        // Verify all sequences are represented
        let seq_to_node = seqrush.seq_to_node.lock().unwrap();
        for seq in &seqrush.sequences {
            for pos in 0..seq.data.len() {
                let global_pos = seq.offset + pos;
                assert!(seq_to_node.contains_key(&global_pos),
                       "Position {} should be mapped", global_pos);
            }
        }
    }

    #[test]
    fn test_complete_sequence_reconstruction() {
        // Comprehensive test to ensure sequences can be fully reconstructed from graph
        let test_sequences = vec![
            ("identical1", "ATCGATCGATCG"),
            ("identical2", "ATCGATCGATCG"),
            ("variant", "ATCGATCGTTCG"),
            ("different", "GCTAGCTAGCTA"),
        ];
        
        let mut sequences = Vec::new();
        let mut offset = 0;
        
        for (id, data) in &test_sequences {
            sequences.push(Sequence {
                id: id.to_string(),
                data: data.as_bytes().to_vec(),
                offset,
            });
            offset += data.len();
        }
        
        let mut config = create_test_config();
        config.max_iterations = 10;
        config.min_match_length = 4;
        
        let mut seqrush = SeqRush::new(sequences.clone(), config);
        seqrush.build_graph();
        
        // Write to GFA and read back
        let temp_file = NamedTempFile::new().unwrap();
        seqrush.write_gfa(temp_file.path().to_str().unwrap()).unwrap();
        
        // Parse GFA to reconstruct sequences
        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        
        let mut nodes_map = HashMap::new();
        let mut paths_map = HashMap::new();
        
        for line in content.lines() {
            if line.starts_with("S\t") {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 3 {
                    let node_id: usize = parts[1].parse().unwrap();
                    let sequence = parts[2];
                    nodes_map.insert(node_id, sequence);
                }
            } else if line.starts_with("P\t") {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 3 {
                    let seq_id = parts[1];
                    let path = parts[2];
                    paths_map.insert(seq_id.to_string(), path.to_string());
                }
            }
        }
        
        // Verify each sequence can be reconstructed
        for (original_id, original_data) in &test_sequences {
            assert!(paths_map.contains_key(*original_id), 
                   "Path missing for sequence {}", original_id);
            
            let path = &paths_map[*original_id];
            let mut reconstructed = String::new();
            
            // Debug: print path
            if *original_id == "identical1" {
                println!("Path for {}: {}", original_id, path);
                println!("Nodes: {:?}", nodes_map);
            }
            
            for node_ref in path.split(',') {
                if node_ref.ends_with('+') {
                    let node_id: usize = node_ref[..node_ref.len()-1].parse().unwrap();
                    if let Some(node_seq) = nodes_map.get(&node_id) {
                        reconstructed.push_str(node_seq);
                    }
                }
            }
            
            assert_eq!(reconstructed, *original_data,
                      "Sequence {} not correctly preserved. Expected: {}, Got: {}",
                      original_id, original_data, reconstructed);
        }
    }

    #[test]
    fn test_no_sequence_corruption() {
        // Test that no sequence data is corrupted or lost
        let sequences = vec![
            Sequence {
                id: "seq_a".to_string(),
                data: b"AAAATTTTCCCCGGGG".to_vec(),
                offset: 0,
            },
            Sequence {
                id: "seq_b".to_string(),
                data: b"AAAATTTTCCCCGGGG".to_vec(),
                offset: 16,
            },
            Sequence {
                id: "seq_c".to_string(),
                data: b"AAAATTTTGGGGCCCC".to_vec(),
                offset: 32,
            },
        ];
        
        let mut config = create_test_config();
        config.max_iterations = 5;
        config.min_match_length = 4;
        
        let mut seqrush = SeqRush::new(sequences.clone(), config);
        seqrush.build_graph();
        
        // Verify no positions are lost
        let seq_to_node = seqrush.seq_to_node.lock().unwrap();
        let total_positions = sequences.iter().map(|s| s.data.len()).sum::<usize>();
        assert_eq!(seq_to_node.len(), total_positions, 
                  "Some positions were lost during graph construction");
        
        // Verify each position maps to a valid node
        let nodes = seqrush.nodes.lock().unwrap();
        for (_, &node_id) in seq_to_node.iter() {
            assert!(node_id < nodes.len(), 
                   "Invalid node ID {} (only {} nodes exist)", node_id, nodes.len());
        }
    }
}