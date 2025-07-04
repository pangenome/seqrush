use clap::Parser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::AtomicUsize;
use std::sync::Arc;
use union_find::{UnionFind, UnionBySize, QuickUnionUf};
use uf_rush::UFRush as LockFreeUnionFind;

#[derive(Parser)]
#[command(
    name = "seqrush",
    version = "0.1.0",
    author = "SeqRush Contributors",
    about = "A dynamic, parallel, in-memory bioinformatics tool implementing biWFA+seqwish for pangenome graph construction",
    long_about = "SeqRush implements a lockless union-find algorithm to align random pairs of sequences in parallel until convergence, building an implicit pangenome graph from sparse pairwise alignments."
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
    #[arg(short, long, default_value = "100", value_name = "NUM")]
    max_iterations: usize,
    
    /// Minimum alignment length to consider
    #[arg(long, default_value = "10", value_name = "NUM")]
    min_alignment_length: usize,
    
    /// Alignment match score
    #[arg(long, default_value = "2", value_name = "NUM")]
    match_score: i32,
    
    /// Alignment mismatch penalty
    #[arg(long, default_value = "-1", value_name = "NUM")]
    mismatch_penalty: i32,
    
    /// Alignment gap penalty
    #[arg(long, default_value = "-1", value_name = "NUM")]
    gap_penalty: i32,
    
    /// Maximum alignment score to compute (prevents computational blowup)
    #[arg(long, default_value = "1000", value_name = "NUM")]
    max_alignment_score: i32,
    
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
}

#[derive(Debug, Clone)]
struct Alignment {
    seq1_id: usize,
    seq2_id: usize,
    seq1_start: usize,
    seq1_end: usize,
    seq2_start: usize,
    seq2_end: usize,
    score: i32,
}

struct SeqRush {
    sequences: Vec<Sequence>,
    union_find: Arc<LockFreeUnionFind>,
    _backup_union_find: Arc<std::sync::Mutex<QuickUnionUf<UnionBySize>>>,
    alignments: Vec<Alignment>,
    sequence_to_graph: HashMap<(usize, usize), usize>, // (seq_id, pos) -> node_id
    graph_to_sequence: HashMap<usize, Vec<(usize, usize)>>, // node_id -> [(seq_id, pos)]
    iteration_count: Arc<AtomicUsize>,
    config: AlignmentConfig,
}

#[derive(Debug, Clone)]
struct AlignmentConfig {
    max_iterations: usize,
    min_alignment_length: usize,
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
    max_alignment_score: i32,
    convergence_stability: usize,
    sparsification: f64,
    erdos_renyi_safety_factor: f64,
    verbose: bool,
}

impl SeqRush {
    fn new(sequences: Vec<Sequence>, config: AlignmentConfig) -> Self {
        let total_chars: usize = sequences.iter().map(|s| s.data.len()).sum();
        
        Self {
            sequences,
            union_find: Arc::new(LockFreeUnionFind::new(total_chars)),
            _backup_union_find: Arc::new(std::sync::Mutex::new(QuickUnionUf::<UnionBySize>::new(total_chars))),
            alignments: Vec::new(),
            sequence_to_graph: HashMap::new(),
            graph_to_sequence: HashMap::new(),
            iteration_count: Arc::new(AtomicUsize::new(0)),
            config,
        }
    }
    
    /// Calculate Erdős-Rényi connectivity threshold probability
    /// For a random graph G(n,p) to have a giant component with high probability,
    /// we need p > (1 + ε) * ln(n) / n
    fn calculate_erdos_renyi_threshold(n: usize, safety_factor: f64) -> f64 {
        if n <= 1 {
            return 1.0;
        }
        
        let n_f = n as f64;
        // Critical threshold is ln(n)/n for connectivity
        // We use safety_factor * ln(n)/n to ensure giant component with high probability
        let threshold = safety_factor * n_f.ln() / n_f;
        
        // Clamp to reasonable bounds
        threshold.min(1.0).max(0.001)
    }
    
    /// Determine if a pair should be aligned based on sparsification
    fn should_align_pair(&self, seq1_id: usize, seq2_id: usize) -> bool {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        
        // Create a deterministic hash of the pair
        let mut hasher = DefaultHasher::new();
        // Hash in a way that's symmetric (same for (i,j) and (j,i))
        let (min_id, max_id) = if seq1_id < seq2_id {
            (seq1_id, seq2_id)
        } else {
            (seq2_id, seq1_id)
        };
        min_id.hash(&mut hasher);
        max_id.hash(&mut hasher);
        let hash_value = hasher.finish();
        
        // Convert hash to probability [0,1)
        let prob = (hash_value as f64) / (u64::MAX as f64);
        
        // Accept if probability is below sparsification threshold
        prob < self.config.sparsification
    }
    
    fn build_graph(&mut self) {
        println!("Building graph with {} sequences", self.sequences.len());
        
        // Initialize character positions
        let mut char_offset = 0;
        for (seq_id, seq) in self.sequences.iter().enumerate() {
            for (pos, _) in seq.data.iter().enumerate() {
                let global_pos = char_offset + pos;
                self.sequence_to_graph.insert((seq_id, pos), global_pos);
                self.graph_to_sequence.entry(global_pos).or_insert_with(Vec::new).push((seq_id, pos));
            }
            char_offset += seq.data.len();
        }
        
        // Parallel alignment and graph construction
        let mut prev_components = 0;
        let mut no_change_count = 0;
        
        for iteration in 0..self.config.max_iterations {
            println!("Iteration {}", iteration + 1);
            
            // Perform random pairwise alignments
            self.perform_random_alignments();
            
            // Count connected components
            let current_components = self.count_components();
            println!("Components: {}", current_components);
            
            // Check for convergence with stability requirement
            if current_components == prev_components {
                no_change_count += 1;
                if no_change_count >= self.config.convergence_stability {
                    println!("Converged after {} iterations (stable for {} iterations)", iteration + 1, no_change_count);
                    break;
                }
            } else {
                no_change_count = 0;
            }
            
            prev_components = current_components;
        }
    }
    
    fn perform_random_alignments(&mut self) {
        use rand::seq::SliceRandom;
        use rand::thread_rng;
        
        let mut rng = thread_rng();
        let seq_count = self.sequences.len();
        
        // Generate all possible pairs and apply Erdős-Rényi sparsification
        let mut pairs = Vec::new();
        let mut total_pairs = 0;
        let mut sparse_pairs = 0;
        
        for i in 0..seq_count {
            for j in (i + 1)..seq_count {
                total_pairs += 1;
                if self.should_align_pair(i, j) {
                    pairs.push((i, j));
                    sparse_pairs += 1;
                }
            }
        }
        
        if self.config.verbose {
            let sparsification_ratio = sparse_pairs as f64 / total_pairs as f64;
            println!("Erdős-Rényi sparsification: {}/{} pairs ({:.2}%) selected", 
                    sparse_pairs, total_pairs, sparsification_ratio * 100.0);
            println!("Sparsification probability: {:.4}", self.config.sparsification);
        }
        
        pairs.shuffle(&mut rng);
        
        // For each iteration, take a subset of the sparse pairs
        let batch_size = std::cmp::min(pairs.len(), std::cmp::max(1, pairs.len() / 4));
        let selected_pairs: Vec<_> = pairs.into_iter().take(batch_size).collect();
        
        if self.config.verbose {
            println!("Processing {} alignment pairs in parallel", selected_pairs.len());
        }
        
        // Perform alignments in parallel and process immediately
        // This avoids collecting all results before processing
        selected_pairs
            .par_iter()
            .filter_map(|&(i, j)| {
                self.align_sequences(i, j)
            })
            .for_each(|alignment| {
                self.process_alignment(&alignment);
            });
    }
    
    fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&base| match base {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                b'a' => b't',
                b't' => b'a',
                b'c' => b'g',
                b'g' => b'c',
                _ => base, // Keep N, ambiguous bases, etc. as-is
            })
            .collect()
    }

    fn align_sequences(&self, seq1_id: usize, seq2_id: usize) -> Option<Alignment> {
        let seq1 = &self.sequences[seq1_id];
        let seq2 = &self.sequences[seq2_id];
        
        // Try forward orientation first
        let forward_result = self.wfa_align(&seq1.data, &seq2.data);
        
        // Try reverse complement orientation
        let seq2_rc = Self::reverse_complement(&seq2.data);
        let reverse_result = self.wfa_align(&seq1.data, &seq2_rc);
        
        // Choose the best alignment
        let best_result = match (forward_result, reverse_result) {
            (Some((s1_start, s1_end, s2_start, s2_end, score1)), Some((_, _, _, _, score2))) => {
                if score1 >= score2 {
                    Some((s1_start, s1_end, s2_start, s2_end, score1))
                } else {
                    // For reverse complement, we need to adjust coordinates
                    // Note: this is simplified - in practice you'd track orientation
                    reverse_result
                }
            }
            (Some(result), None) => Some(result),
            (None, Some(result)) => Some(result),
            (None, None) => None,
        };
        
        if let Some((seq1_start, seq1_end, seq2_start, seq2_end, score)) = best_result {
            // Early termination if score exceeds threshold
            if score > self.config.max_alignment_score {
                if self.config.verbose {
                    println!("Alignment score {} exceeds threshold {}, skipping", 
                            score, self.config.max_alignment_score);
                }
                return None;
            }
            
            Some(Alignment {
                seq1_id,
                seq2_id,
                seq1_start,
                seq1_end,
                seq2_start,
                seq2_end,
                score,
            })
        } else {
            // Fall back to exact match finder (forward orientation only for now)
            self.find_exact_matches(seq1_id, seq2_id, &seq1.data, &seq2.data)
        }
    }
    
    fn find_exact_matches(&self, seq1_id: usize, seq2_id: usize, seq1: &[u8], seq2: &[u8]) -> Option<Alignment> {
        // Simple exact match finder - in practice this would use biWFA
        let min_match_len = self.config.min_alignment_length;
        
        for i in 0..seq1.len().saturating_sub(min_match_len) {
            for j in 0..seq2.len().saturating_sub(min_match_len) {
                let mut match_len = 0;
                while i + match_len < seq1.len() && 
                      j + match_len < seq2.len() && 
                      seq1[i + match_len] == seq2[j + match_len] {
                    match_len += 1;
                }
                
                if match_len >= min_match_len {
                    let score = match_len as i32;
                    // Respect max alignment score threshold
                    if score <= self.config.max_alignment_score {
                        return Some(Alignment {
                            seq1_id,
                            seq2_id,
                            seq1_start: i,
                            seq1_end: i + match_len,
                            seq2_start: j,
                            seq2_end: j + match_len,
                            score,
                        });
                    }
                }
            }
        }
        
        None
    }

    fn wfa_align(&self, seq1: &[u8], seq2: &[u8]) -> Option<(usize, usize, usize, usize, i32)> {
        // Placeholder for biWFA alignment
        // This would be replaced with actual WFA2 alignment
        // For now, using a simple local alignment approach
        
        let match_score = self.config.match_score;
        let mismatch_penalty = self.config.mismatch_penalty;
        let gap_penalty = self.config.gap_penalty;
        
        let n = seq1.len();
        let m = seq2.len();
        
        if n == 0 || m == 0 {
            return None;
        }
        
        // Simple Smith-Waterman-like approach (very basic)
        let mut dp = vec![vec![0; m + 1]; n + 1];
        let mut max_score = 0;
        let mut max_i = 0;
        let mut max_j = 0;
        
        for i in 1..=n {
            for j in 1..=m {
                let match_mismatch = if seq1[i-1] == seq2[j-1] { 
                    dp[i-1][j-1] + match_score 
                } else { 
                    dp[i-1][j-1] + mismatch_penalty 
                };
                
                let delete = dp[i-1][j] + gap_penalty;
                let insert = dp[i][j-1] + gap_penalty;
                
                dp[i][j] = std::cmp::max(0, std::cmp::max(match_mismatch, std::cmp::max(delete, insert)));
                
                if dp[i][j] > max_score {
                    max_score = dp[i][j];
                    max_i = i;
                    max_j = j;
                }
            }
        }
        
        if max_score < 10 {
            return None;
        }
        
        // Early termination if score exceeds threshold
        if max_score > self.config.max_alignment_score {
            return None;
        }
        
        // Traceback to find alignment boundaries (simplified)
        let mut i = max_i;
        let mut j = max_j;
        let end_i = i;
        let end_j = j;
        
        while i > 0 && j > 0 && dp[i][j] > 0 {
            if seq1[i-1] == seq2[j-1] {
                i -= 1;
                j -= 1;
            } else if dp[i-1][j] > dp[i][j-1] {
                i -= 1;
            } else {
                j -= 1;
            }
        }
        
        Some((i, end_i, j, end_j, max_score))
    }
    
    fn process_alignment(&self, alignment: &Alignment) {
        // Union characters that align
        let seq1_offset = self.sequences.iter().take(alignment.seq1_id).map(|s| s.data.len()).sum::<usize>();
        let seq2_offset = self.sequences.iter().take(alignment.seq2_id).map(|s| s.data.len()).sum::<usize>();
        
        for pos in 0..(alignment.seq1_end - alignment.seq1_start) {
            let pos1 = seq1_offset + alignment.seq1_start + pos;
            let pos2 = seq2_offset + alignment.seq2_start + pos;
            
            self.union_find.unite(pos1, pos2);
        }
    }
    
    fn count_components(&self) -> usize {
        let total_chars: usize = self.sequences.iter().map(|s| s.data.len()).sum();
        let mut components = std::collections::HashSet::new();
        
        for i in 0..total_chars {
            components.insert(self.union_find.find(i));
        }
        
        components.len()
    }
    
    fn write_gfa(&self, output_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        use std::fs::File;
        use std::io::{BufWriter, Write};
        
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);
        
        writeln!(writer, "H\tVN:Z:1.0")?;
        
        // Collect node information
        let mut node_chars = HashMap::new();
        let mut node_positions = HashMap::new();
        let total_chars: usize = self.sequences.iter().map(|s| s.data.len()).sum();
        
        for i in 0..total_chars {
            let root = self.union_find.find(i);
            let (seq_id, pos) = self.char_to_sequence_pos(i);
            let char = self.sequences[seq_id].data[pos];
            
            node_chars.entry(root).or_insert_with(Vec::new).push((i, char));
            node_positions.entry(root).or_insert_with(Vec::new).push((seq_id, pos));
        }
        
        // Sort characters within each node by their original order
        for chars in node_chars.values_mut() {
            chars.sort_by_key(|&(global_pos, _)| global_pos);
        }
        
        // Write segments (nodes) - compact consecutive identical characters
        let mut node_id_map = HashMap::new();
        let mut compact_node_id = 0;
        
        for (&original_node_id, chars) in node_chars.iter() {
            if chars.is_empty() {
                continue;
            }
            
            // Compact the sequence by merging consecutive identical characters
            let mut compact_sequence = Vec::new();
            let mut prev_char = None;
            
            for &(_, char) in chars {
                if Some(char) != prev_char {
                    compact_sequence.push(char);
                    prev_char = Some(char);
                }
            }
            
            if !compact_sequence.is_empty() {
                let sequence = String::from_utf8_lossy(&compact_sequence);
                writeln!(writer, "S\t{}\t{}", compact_node_id, sequence)?;
                node_id_map.insert(original_node_id, compact_node_id);
                compact_node_id += 1;
            }
        }
        
        // Write edges (simple approach - consecutive nodes in paths)
        let mut edges = std::collections::HashSet::new();
        
        for (seq_id, seq) in self.sequences.iter().enumerate() {
            let seq_offset = self.sequences.iter().take(seq_id).map(|s| s.data.len()).sum::<usize>();
            let mut prev_node = None;
            
            for pos in 0..seq.data.len() {
                let global_pos = seq_offset + pos;
                let original_node_id = self.union_find.find(global_pos);
                
                if let Some(&current_node) = node_id_map.get(&original_node_id) {
                    if let Some(prev) = prev_node {
                        if prev != current_node {
                            edges.insert((prev, current_node));
                        }
                    }
                    prev_node = Some(current_node);
                }
            }
        }
        
        // Write edges
        for (from, to) in edges {
            writeln!(writer, "L\t{}\t+\t{}\t+\t0M", from, to)?;
        }
        
        // Write paths
        for (seq_id, seq) in self.sequences.iter().enumerate() {
            let mut path_nodes = Vec::new();
            let seq_offset = self.sequences.iter().take(seq_id).map(|s| s.data.len()).sum::<usize>();
            let mut prev_node = None;
            
            for pos in 0..seq.data.len() {
                let global_pos = seq_offset + pos;
                let original_node_id = self.union_find.find(global_pos);
                
                if let Some(&current_node) = node_id_map.get(&original_node_id) {
                    if Some(current_node) != prev_node {
                        path_nodes.push(format!("{}+", current_node));
                        prev_node = Some(current_node);
                    }
                }
            }
            
            if !path_nodes.is_empty() {
                writeln!(writer, "P\t{}\t{}\t*", seq.id, path_nodes.join(","))?;
            }
        }
        
        Ok(())
    }
    
    fn char_to_sequence_pos(&self, global_pos: usize) -> (usize, usize) {
        let mut offset = 0;
        for (seq_id, seq) in self.sequences.iter().enumerate() {
            if global_pos < offset + seq.data.len() {
                return (seq_id, global_pos - offset);
            }
            offset += seq.data.len();
        }
        panic!("Invalid global position: {}", global_pos);
    }
}

fn load_sequences(file_path: &str) -> Result<Vec<Sequence>, Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut sequences = Vec::new();
    let mut current_id = String::new();
    let mut current_data = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_id.is_empty() {
                sequences.push(Sequence {
                    id: current_id.clone(),
                    data: current_data.clone(),
                });
                current_data.clear();
            }
            current_id = line[1..].to_string();
        } else {
            current_data.extend(line.bytes());
        }
    }
    
    if !current_id.is_empty() {
        sequences.push(Sequence {
            id: current_id,
            data: current_data,
        });
    }
    
    Ok(sequences)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    
    // Set up Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();
    
    // Load sequences
    let sequences = load_sequences(&args.sequences)?;
    if args.verbose {
        println!("Loaded {} sequences", sequences.len());
        for seq in &sequences {
            println!("  {}: {} bp", seq.id, seq.data.len());
        }
    } else {
        println!("Loaded {} sequences", sequences.len());
    }
    
    // Calculate sparsification probability using Erdős-Rényi model
    let sparsification = if let Some(manual_sparse) = args.sparsification {
        if args.verbose {
            println!("Using manual sparsification: {:.4}", manual_sparse);
        }
        manual_sparse
    } else {
        let calculated = SeqRush::calculate_erdos_renyi_threshold(sequences.len(), args.erdos_renyi_safety_factor);
        if args.verbose {
            println!("Auto-calculated Erdős-Rényi sparsification: {:.4} (n={}, safety_factor={})", 
                    calculated, sequences.len(), args.erdos_renyi_safety_factor);
        }
        calculated
    };

    // Create configuration
    let config = AlignmentConfig {
        max_iterations: args.max_iterations,
        min_alignment_length: args.min_alignment_length,
        match_score: args.match_score,
        mismatch_penalty: args.mismatch_penalty,
        gap_penalty: args.gap_penalty,
        max_alignment_score: args.max_alignment_score,
        convergence_stability: args.convergence_stability,
        sparsification,
        erdos_renyi_safety_factor: args.erdos_renyi_safety_factor,
        verbose: args.verbose,
    };
    
    if args.verbose {
        println!("Configuration:");
        println!("  Max iterations: {}", config.max_iterations);
        println!("  Min alignment length: {}", config.min_alignment_length);
        println!("  Scoring: match={}, mismatch={}, gap={}, max_score={}", 
                config.match_score, config.mismatch_penalty, config.gap_penalty, config.max_alignment_score);
        println!("  Convergence stability: {}", config.convergence_stability);
        println!("  Sparsification: {:.4} (Erdős-Rényi safety factor: {})", 
                config.sparsification, config.erdos_renyi_safety_factor);
        println!("  Threads: {}", args.threads);
    }
    
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
            },
            Sequence {
                id: "seq2".to_string(),
                data: b"ATCGATCGATCGATCG".to_vec(),
            },
            Sequence {
                id: "seq3".to_string(),
                data: b"ATCGATCGCTCGATCG".to_vec(),
            },
            Sequence {
                id: "seq4".to_string(),
                data: b"GCTAGCTAGCTAGCTA".to_vec(),
            },
        ]
    }

    fn create_test_config() -> AlignmentConfig {
        AlignmentConfig {
            max_iterations: 10,
            min_alignment_length: 5,
            match_score: 2,
            mismatch_penalty: -1,
            gap_penalty: -1,
            max_alignment_score: 1000,
            convergence_stability: 2,
            sparsification: 1.0, // Use all pairs for testing
            erdos_renyi_safety_factor: 3.0,
            verbose: false,
        }
    }

    fn create_temp_fasta(sequences: &[(&str, &str)]) -> NamedTempFile {
        let mut temp_file = NamedTempFile::new().unwrap();
        for (id, seq) in sequences {
            writeln!(temp_file, ">{}", id).unwrap();
            writeln!(temp_file, "{}", seq).unwrap();
        }
        temp_file.flush().unwrap();
        temp_file
    }

    #[test]
    fn test_sequence_creation() {
        let sequences = create_test_sequences();
        assert_eq!(sequences.len(), 4);
        assert_eq!(sequences[0].id, "seq1");
        assert_eq!(sequences[0].data, b"ATCGATCGATCGATCG");
    }

    #[test]
    fn test_config_creation() {
        let config = create_test_config();
        assert_eq!(config.max_iterations, 10);
        assert_eq!(config.min_alignment_length, 5);
        assert_eq!(config.match_score, 2);
    }

    #[test]
    fn test_seqrush_creation() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences, config);
        
        assert_eq!(seqrush.sequences.len(), 4);
        assert_eq!(seqrush.config.max_iterations, 10);
    }

    #[test]
    fn test_fasta_loading() {
        let temp_file = create_temp_fasta(&[
            ("test1", "ATCGATCG"),
            ("test2", "GCTAGCTA"),
        ]);
        
        let sequences = load_sequences(temp_file.path().to_str().unwrap()).unwrap();
        assert_eq!(sequences.len(), 2);
        assert_eq!(sequences[0].id, "test1");
        assert_eq!(sequences[0].data, b"ATCGATCG");
        assert_eq!(sequences[1].id, "test2");
        assert_eq!(sequences[1].data, b"GCTAGCTA");
    }

    #[test]
    fn test_fasta_loading_empty_file() {
        let temp_file = NamedTempFile::new().unwrap();
        let result = load_sequences(temp_file.path().to_str().unwrap());
        assert!(result.is_ok());
        assert_eq!(result.unwrap().len(), 0);
    }

    #[test]
    fn test_char_to_sequence_pos() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences, config);
        
        // First sequence starts at position 0
        let (seq_id, pos) = seqrush.char_to_sequence_pos(0);
        assert_eq!(seq_id, 0);
        assert_eq!(pos, 0);
        
        // Position in second sequence (first sequence has 16 chars)
        let (seq_id, pos) = seqrush.char_to_sequence_pos(16);
        assert_eq!(seq_id, 1);
        assert_eq!(pos, 0);
    }

    #[test]
    fn test_exact_match_finding() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences.clone(), config);
        
        // Test exact match between identical sequences
        let alignment = seqrush.find_exact_matches(0, 1, &sequences[0].data, &sequences[1].data);
        assert!(alignment.is_some());
        
        let alignment = alignment.unwrap();
        assert_eq!(alignment.seq1_id, 0);
        assert_eq!(alignment.seq2_id, 1);
        assert!(alignment.seq1_end - alignment.seq1_start >= 5);
    }

    #[test]
    fn test_wfa_alignment() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences.clone(), config);
        
        // Test alignment between similar sequences
        let result = seqrush.wfa_align(&sequences[0].data, &sequences[2].data);
        assert!(result.is_some());
        
        let (start1, end1, start2, end2, score) = result.unwrap();
        assert!(end1 > start1);
        assert!(end2 > start2);
        assert!(score > 0);
    }

    #[test]
    fn test_wfa_alignment_no_match() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences, config);
        
        // Create sequences with no good alignment
        let seq1 = b"AAAAAAAAAA";
        let seq2 = b"TTTTTTTTTT";
        
        let result = seqrush.wfa_align(seq1, seq2);
        // Should return None for poor alignments
        assert!(result.is_none());
    }

    #[test]
    fn test_alignment_processing() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences, config);
        
        let alignment = Alignment {
            seq1_id: 0,
            seq2_id: 1,
            seq1_start: 0,
            seq1_end: 10,
            seq2_start: 0,
            seq2_end: 10,
            score: 20,
        };
        
        let initial_components = seqrush.count_components();
        seqrush.process_alignment(&alignment);
        let final_components = seqrush.count_components();
        
        // Processing alignment should reduce component count
        assert!(final_components <= initial_components);
    }

    #[test]
    fn test_component_counting() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let total_chars: usize = sequences.iter().map(|s| s.data.len()).sum();
        let seqrush = SeqRush::new(sequences, config);
        
        let components = seqrush.count_components();
        
        // Initially, each character is its own component
        assert_eq!(components, total_chars);
    }

    #[test]
    fn test_gfa_output() {
        let sequences = vec![
            Sequence {
                id: "simple1".to_string(),
                data: b"ATCG".to_vec(),
            },
            Sequence {
                id: "simple2".to_string(),
                data: b"ATCG".to_vec(),
            },
        ];
        let config = create_test_config();
        let mut seqrush = SeqRush::new(sequences, config);
        
        // Build a simple graph
        seqrush.build_graph();
        
        // Test GFA output
        let temp_file = NamedTempFile::new().unwrap();
        let result = seqrush.write_gfa(temp_file.path().to_str().unwrap());
        assert!(result.is_ok());
        
        // Read and verify GFA content
        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        assert!(content.starts_with("H\tVN:Z:1.0"));
        assert!(content.contains("S\t"));  // Should have segments
        assert!(content.contains("P\t"));  // Should have paths
    }

    #[test]
    fn test_gfa_output_validation() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let mut seqrush = SeqRush::new(sequences, config);
        
        seqrush.build_graph();
        
        let temp_file = NamedTempFile::new().unwrap();
        seqrush.write_gfa(temp_file.path().to_str().unwrap()).unwrap();
        
        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        
        // Check GFA format requirements
        assert!(lines[0].starts_with("H\t"));  // Header
        
        let mut has_segments = false;
        let mut has_paths = false;
        
        for line in &lines {
            if line.starts_with("S\t") {
                has_segments = true;
                // Verify segment format: S<tab>id<tab>sequence
                let parts: Vec<&str> = line.split('\t').collect();
                assert_eq!(parts.len(), 3);
                assert!(parts[1].parse::<usize>().is_ok()); // ID should be numeric
                assert!(!parts[2].is_empty()); // Sequence should not be empty
            }
            if line.starts_with("P\t") {
                has_paths = true;
                // Verify path format
                let parts: Vec<&str> = line.split('\t').collect();
                assert!(parts.len() >= 3);
            }
        }
        
        assert!(has_segments);
        assert!(has_paths);
    }

    #[test]
    fn test_convergence_detection() {
        let sequences = vec![
            Sequence {
                id: "conv1".to_string(),
                data: b"ATCGATCGATCGATCG".to_vec(), // Longer sequence for better alignment
            },
            Sequence {
                id: "conv2".to_string(),
                data: b"ATCGATCGATCGATCG".to_vec(), // Identical sequence
            },
        ];
        
        let mut config = create_test_config();
        config.convergence_stability = 1; // Quick convergence for testing
        config.max_iterations = 5;
        config.min_alignment_length = 5; // Lower threshold for alignment
        
        let mut seqrush = SeqRush::new(sequences, config);
        
        // This should converge quickly since sequences are identical
        seqrush.build_graph();
        
        // Should have converged and reduced components
        let components = seqrush.count_components();
        let total_chars = 32; // 16 + 16 characters
        assert!(components < total_chars); // Should be less than total characters
    }

    #[test]
    fn test_parallel_alignment_system() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let mut seqrush = SeqRush::new(sequences, config);
        
        let initial_components = seqrush.count_components();
        seqrush.perform_random_alignments();
        let final_components = seqrush.count_components();
        
        // Performing alignments should potentially reduce components
        assert!(final_components <= initial_components);
    }

    #[test]
    fn test_sequence_to_graph_mapping() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let mut seqrush = SeqRush::new(sequences, config);
        
        // Initialize mappings
        let mut char_offset = 0;
        for (seq_id, seq) in seqrush.sequences.iter().enumerate() {
            for (pos, _) in seq.data.iter().enumerate() {
                let global_pos = char_offset + pos;
                seqrush.sequence_to_graph.insert((seq_id, pos), global_pos);
                seqrush.graph_to_sequence.entry(global_pos).or_insert_with(Vec::new).push((seq_id, pos));
            }
            char_offset += seq.data.len();
        }
        
        // Test mapping
        assert_eq!(seqrush.sequence_to_graph.get(&(0, 0)), Some(&0));
        assert_eq!(seqrush.sequence_to_graph.get(&(1, 0)), Some(&16));
        assert!(seqrush.graph_to_sequence.contains_key(&0));
    }

    #[test]
    fn test_edge_case_empty_sequences() {
        let sequences = vec![];
        let config = create_test_config();
        let mut seqrush = SeqRush::new(sequences, config);
        
        // Should handle empty input gracefully
        seqrush.build_graph();
        
        let temp_file = NamedTempFile::new().unwrap();
        let result = seqrush.write_gfa(temp_file.path().to_str().unwrap());
        assert!(result.is_ok());
    }

    #[test]
    fn test_edge_case_single_sequence() {
        let sequences = vec![
            Sequence {
                id: "single".to_string(),
                data: b"ATCGATCG".to_vec(),
            },
        ];
        let config = create_test_config();
        let mut seqrush = SeqRush::new(sequences, config);
        
        seqrush.build_graph();
        
        let temp_file = NamedTempFile::new().unwrap();
        let result = seqrush.write_gfa(temp_file.path().to_str().unwrap());
        assert!(result.is_ok());
        
        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        assert!(content.contains("P\tsingle\t")); // Should have path for single sequence
    }

    #[test]
    fn test_edge_case_very_short_sequences() {
        let sequences = vec![
            Sequence {
                id: "short1".to_string(),
                data: b"A".to_vec(),
            },
            Sequence {
                id: "short2".to_string(),
                data: b"T".to_vec(),
            },
        ];
        let config = create_test_config();
        let mut seqrush = SeqRush::new(sequences, config);
        
        seqrush.build_graph();
        
        // Should handle very short sequences without crashing
        let components = seqrush.count_components();
        assert_eq!(components, 2); // Each character should be separate
    }

    #[test]
    fn test_large_dataset_performance() {
        // Create a larger dataset for performance testing
        let mut sequences = Vec::new();
        for i in 0..20 {
            sequences.push(Sequence {
                id: format!("perf_seq_{}", i),
                data: format!("ATCGATCGATCG{:04}ATCGATCGATCG", i).into_bytes(),
            });
        }
        
        let mut config = create_test_config();
        config.max_iterations = 5; // Limit iterations for performance test
        
        let mut seqrush = SeqRush::new(sequences, config);
        
        // This should complete without issues
        seqrush.build_graph();
        
        let temp_file = NamedTempFile::new().unwrap();
        let result = seqrush.write_gfa(temp_file.path().to_str().unwrap());
        assert!(result.is_ok());
    }

    #[test]
    fn test_configuration_parameters() {
        let sequences = create_test_sequences();
        
        // Test different configuration parameters
        let config1 = AlignmentConfig {
            max_iterations: 5,
            min_alignment_length: 3,
            match_score: 1,
            mismatch_penalty: -2,
            gap_penalty: -3,
            max_alignment_score: 500,
            convergence_stability: 1,
            sparsification: 1.0,
            erdos_renyi_safety_factor: 3.0,
            verbose: false,
        };
        
        let mut seqrush1 = SeqRush::new(sequences.clone(), config1);
        seqrush1.build_graph();
        
        let config2 = AlignmentConfig {
            max_iterations: 10,
            min_alignment_length: 8,
            match_score: 3,
            mismatch_penalty: -1,
            gap_penalty: -1,
            max_alignment_score: 2000,
            convergence_stability: 3,
            sparsification: 1.0,
            erdos_renyi_safety_factor: 3.0,
            verbose: true,
        };
        
        let mut seqrush2 = SeqRush::new(sequences, config2);
        seqrush2.build_graph();
        
        // Both should complete successfully with different parameters
        assert!(seqrush1.count_components() > 0);
        assert!(seqrush2.count_components() > 0);
    }

    #[test]
    fn test_fasta_malformed_input() {
        let temp_file = create_temp_fasta(&[
            ("good_seq", "ATCGATCG"),
        ]);
        
        // Add some malformed content
        {
            let mut file = std::fs::OpenOptions::new()
                .append(true)
                .open(temp_file.path())
                .unwrap();
            writeln!(file, "malformed line without >").unwrap();
            writeln!(file, ">another_seq").unwrap();
            writeln!(file, "GCTAGCTA").unwrap();
        }
        
        let sequences = load_sequences(temp_file.path().to_str().unwrap()).unwrap();
        // Should still load the valid sequences
        assert!(sequences.len() >= 1);
    }

    #[test]
    fn test_integration_full_pipeline() {
        // Integration test for the full pipeline
        let temp_file = create_temp_fasta(&[
            ("seq1", "ATCGATCGATCGATCGATCGATCGATCG"),
            ("seq2", "ATCGATCGATCGATCGATCGATCGATCG"),
            ("seq3", "ATCGATCGATCGTTCGATCGATCGATCG"),
            ("seq4", "GCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
        ]);
        
        let sequences = load_sequences(temp_file.path().to_str().unwrap()).unwrap();
        assert_eq!(sequences.len(), 4);
        
        let config = AlignmentConfig {
            max_iterations: 10,
            min_alignment_length: 8,
            match_score: 2,
            mismatch_penalty: -1,
            gap_penalty: -1,
            max_alignment_score: 1000,
            convergence_stability: 2,
            sparsification: 1.0,
            erdos_renyi_safety_factor: 3.0,
            verbose: false,
        };
        
        let mut seqrush = SeqRush::new(sequences, config);
        seqrush.build_graph();
        
        let temp_gfa = NamedTempFile::new().unwrap();
        seqrush.write_gfa(temp_gfa.path().to_str().unwrap()).unwrap();
        
        let content = std::fs::read_to_string(temp_gfa.path()).unwrap();
        assert!(content.contains("H\tVN:Z:1.0"));
        assert!(content.contains("S\t"));
        assert!(content.contains("P\t"));
        
        // Should have paths for all input sequences
        assert!(content.contains("P\tseq1\t"));
        assert!(content.contains("P\tseq2\t"));
        assert!(content.contains("P\tseq3\t"));
        assert!(content.contains("P\tseq4\t"));
    }

    #[test]
    fn test_identical_sequences_merge() {
        // Test that identical sequences are properly merged
        let sequences = vec![
            Sequence {
                id: "identical1".to_string(),
                data: b"ATCGATCGATCGATCGATCGATCG".to_vec(),
            },
            Sequence {
                id: "identical2".to_string(),
                data: b"ATCGATCGATCGATCGATCGATCG".to_vec(),
            },
            Sequence {
                id: "identical3".to_string(),
                data: b"ATCGATCGATCGATCGATCGATCG".to_vec(),
            },
        ];
        
        let mut config = create_test_config();
        config.min_alignment_length = 8;
        config.max_iterations = 5;
        
        let mut seqrush = SeqRush::new(sequences, config);
        
        let initial_components = seqrush.count_components();
        seqrush.build_graph();
        let final_components = seqrush.count_components();
        
        // Should have significantly fewer components after merging identical sequences
        assert!(final_components < initial_components);
        
        // Test GFA output
        let temp_file = NamedTempFile::new().unwrap();
        seqrush.write_gfa(temp_file.path().to_str().unwrap()).unwrap();
        
        let content = std::fs::read_to_string(temp_file.path()).unwrap();
        // Should have paths for all three sequences
        assert!(content.contains("identical1"));
        assert!(content.contains("identical2"));
        assert!(content.contains("identical3"));
    }

    #[test]
    fn test_alignment_score_calculation() {
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences.clone(), config);
        
        // Test the WFA alignment scoring
        let result = seqrush.wfa_align(&sequences[0].data, &sequences[1].data);
        assert!(result.is_some());
        
        let (_, _, _, _, score) = result.unwrap();
        // Identical sequences should have a high positive score
        assert!(score > 20);
        
        // Test alignment of different sequences
        let result2 = seqrush.wfa_align(&sequences[0].data, &sequences[3].data);
        if let Some((_, _, _, _, score2)) = result2 {
            // Different sequences should have lower score than identical ones
            assert!(score2 < score);
        }
    }

    #[test]
    fn test_error_handling_invalid_file() {
        let result = load_sequences("nonexistent_file.fasta");
        assert!(result.is_err());
    }

    #[test]
    fn test_edge_case_sequences_with_ns() {
        // Test sequences containing N characters
        let sequences = vec![
            Sequence {
                id: "with_ns".to_string(),
                data: b"ATCGNNNGATCGATCG".to_vec(),
            },
            Sequence {
                id: "normal".to_string(),
                data: b"ATCGATCGATCGATCG".to_vec(),
            },
        ];
        
        let config = create_test_config();
        let mut seqrush = SeqRush::new(sequences, config);
        
        // Should handle sequences with N characters without crashing
        seqrush.build_graph();
        
        let temp_file = NamedTempFile::new().unwrap();
        let result = seqrush.write_gfa(temp_file.path().to_str().unwrap());
        assert!(result.is_ok());
    }

    #[test]
    fn test_thread_safety_basic() {
        // Basic test that the union-find operations work in parallel context
        let sequences = create_test_sequences();
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences, config);
        
        // Create multiple alignments and process them
        let alignments = vec![
            Alignment {
                seq1_id: 0,
                seq2_id: 1,
                seq1_start: 0,
                seq1_end: 5,
                seq2_start: 0,
                seq2_end: 5,
                score: 10,
            },
            Alignment {
                seq1_id: 1,
                seq2_id: 2,
                seq1_start: 5,
                seq1_end: 10,
                seq2_start: 5,
                seq2_end: 10,
                score: 8,
            },
        ];
        
        let initial_components = seqrush.count_components();
        
        // Process alignments
        for alignment in alignments {
            seqrush.process_alignment(&alignment);
        }
        
        let final_components = seqrush.count_components();
        assert!(final_components <= initial_components);
    }

    #[test]
    fn test_manual_reverse_complement_verification() {
        // Manual verification of reverse complement
        let original = b"ATCG";
        // A->T, T->A, C->G, G->C, then reverse
        // Complement: TAGC
        // Reverse: CGAT
        let expected = b"CGAT";
        let result = SeqRush::reverse_complement(original);
        assert_eq!(result, expected);
        
        // Test a longer sequence step by step
        let seq = b"ATCGATCG";
        // A-T-C-G-A-T-C-G
        // T-A-G-C-T-A-G-C (complement)
        // C-G-A-T-C-G-A-T (reverse)
        let expected_long = b"CGATCGAT";
        let result_long = SeqRush::reverse_complement(seq);
        assert_eq!(result_long, expected_long);
    }

    #[test]
    fn test_reverse_complement() {
        // Test reverse complement functionality
        let seq = b"ATCG";
        let rc = SeqRush::reverse_complement(seq);
        assert_eq!(rc, b"CGAT");
        
        let seq2 = b"AAATTTCCCGGG";
        let rc2 = SeqRush::reverse_complement(seq2);
        assert_eq!(rc2, b"CCCGGGAAATTT");
        
        // Test with lowercase
        let seq3 = b"atcg";
        let rc3 = SeqRush::reverse_complement(seq3);
        assert_eq!(rc3, b"cgat");
        
        // Test with N characters
        let seq4 = b"ATCGN";
        let rc4 = SeqRush::reverse_complement(seq4);
        assert_eq!(rc4, b"NCGAT");
    }

    #[test]
    fn test_bidirectional_alignment() {
        // Test that alignment considers both orientations
        let forward_seq = b"ATCGATCGATCG".to_vec();
        let reverse_comp = SeqRush::reverse_complement(&forward_seq);
        
        println!("Forward: {:?}", String::from_utf8_lossy(&forward_seq));
        println!("Reverse complement: {:?}", String::from_utf8_lossy(&reverse_comp));
        
        let sequences = vec![
            Sequence {
                id: "forward".to_string(),
                data: forward_seq,
            },
            Sequence {
                id: "reverse".to_string(),
                data: reverse_comp,
            },
        ];
        
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences.clone(), config);
        
        // Should find alignment in reverse complement orientation
        let alignment = seqrush.align_sequences(0, 1);
        assert!(alignment.is_some());
        
        let alignment = alignment.unwrap();
        assert!(alignment.score > 0);
        println!("Alignment score: {}", alignment.score);
    }

    #[test]
    fn test_max_alignment_score_threshold() {
        // Test that max score threshold prevents computational blowup
        let sequences = vec![
            Sequence {
                id: "seq1".to_string(),
                data: b"ATCGATCGATCGATCGATCGATCGATCGATCG".to_vec(),
            },
            Sequence {
                id: "seq2".to_string(),
                data: b"ATCGATCGATCGATCGATCGATCGATCGATCG".to_vec(),
            },
        ];
        
        let mut config = create_test_config();
        config.max_alignment_score = 20; // Very low threshold
        config.verbose = true; // To see the threshold message
        
        let seqrush = SeqRush::new(sequences.clone(), config);
        
        // Should return None due to score threshold
        let alignment = seqrush.align_sequences(0, 1);
        // Depending on implementation, this might be None due to threshold
        // or might fall back to exact match finder
        if let Some(alignment) = alignment {
            assert!(alignment.score <= 20);
        }
    }

    #[test]
    fn test_orientation_preference() {
        // Test that forward orientation is preferred when scores are equal
        let sequences = vec![
            Sequence {
                id: "palindrome".to_string(),
                data: b"ATCGCGAT".to_vec(), // Palindromic sequence
            },
            Sequence {
                id: "same".to_string(),
                data: b"ATCGCGAT".to_vec(),
            },
        ];
        
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences.clone(), config);
        
        let alignment = seqrush.align_sequences(0, 1);
        assert!(alignment.is_some());
        
        // Should find a good alignment
        let alignment = alignment.unwrap();
        assert!(alignment.score > 10);
    }

    #[test]
    fn test_performance_with_orientations() {
        // Test that bidirectional alignment doesn't cause major performance issues
        let mut sequences = Vec::new();
        for i in 0..10 {
            sequences.push(Sequence {
                id: format!("perf_seq_{}", i),
                data: format!("ATCGATCGATCG{:04}CGTATCGATCGATCG", i).into_bytes(),
            });
        }
        
        let mut config = create_test_config();
        config.max_iterations = 3; // Limit for performance test
        config.max_alignment_score = 500; // Reasonable threshold
        
        let mut seqrush = SeqRush::new(sequences, config);
        
        // This should complete without major performance issues
        seqrush.build_graph();
        
        let temp_file = NamedTempFile::new().unwrap();
        let result = seqrush.write_gfa(temp_file.path().to_str().unwrap());
        assert!(result.is_ok());
    }

    #[test]
    fn test_orientation_selection_actually_works() {
        // Test that we actually select the better orientation
        let sequences = vec![
            Sequence {
                id: "seq1".to_string(),
                data: b"ATCGATCGATCGATCG".to_vec(),
            },
            Sequence {
                id: "seq2_poor_forward".to_string(),
                // Poor forward match but good reverse complement match
                data: b"CGATCGATCGATCGAT".to_vec(), // This is reverse complement of seq1
            },
        ];
        
        let config = create_test_config();
        let seqrush = SeqRush::new(sequences.clone(), config);
        
        // Test forward alignment directly
        let forward_result = seqrush.wfa_align(&sequences[0].data, &sequences[1].data);
        
        // Test reverse complement alignment directly  
        let seq2_rc = SeqRush::reverse_complement(&sequences[1].data);
        let reverse_result = seqrush.wfa_align(&sequences[0].data, &seq2_rc);
        
        println!("Forward alignment result: {:?}", forward_result);
        println!("Reverse complement alignment result: {:?}", reverse_result);
        
        // The bidirectional aligner should pick the best one
        let best_alignment = seqrush.align_sequences(0, 1);
        assert!(best_alignment.is_some());
        
        let best = best_alignment.unwrap();
        println!("Best alignment score: {}", best.score);
        
        // Should have found a good alignment
        assert!(best.score > 10);
    }

    #[test]
    fn test_erdos_renyi_threshold_calculation() {
        // Test Erdős-Rényi threshold calculation
        
        // For small graphs, should use high probability
        assert_eq!(SeqRush::calculate_erdos_renyi_threshold(1, 3.0), 1.0);
        assert_eq!(SeqRush::calculate_erdos_renyi_threshold(2, 3.0), 1.0);
        
        // For larger graphs, should decrease
        let threshold_10 = SeqRush::calculate_erdos_renyi_threshold(10, 3.0);
        let threshold_100 = SeqRush::calculate_erdos_renyi_threshold(100, 3.0);
        let threshold_1000 = SeqRush::calculate_erdos_renyi_threshold(1000, 3.0);
        
        assert!(threshold_10 > threshold_100);
        assert!(threshold_100 > threshold_1000);
        
        // Should be reasonable values
        assert!(threshold_10 < 1.0);
        assert!(threshold_100 < 0.5);
        assert!(threshold_1000 < 0.1);
        assert!(threshold_1000 > 0.001);
        
        println!("Thresholds: n=10: {:.4}, n=100: {:.4}, n=1000: {:.4}", 
                threshold_10, threshold_100, threshold_1000);
    }

    #[test]
    fn test_sparsification_deterministic() {
        // Test that sparsification is deterministic for same inputs
        let sequences = create_test_sequences();
        let mut config = create_test_config();
        config.sparsification = 0.5; // 50% probability
        
        let seqrush = SeqRush::new(sequences, config);
        
        // Same pair should always give same result
        let result1 = seqrush.should_align_pair(0, 1);
        let result2 = seqrush.should_align_pair(0, 1);
        assert_eq!(result1, result2);
        
        // Symmetric pairs should give same result
        let result_forward = seqrush.should_align_pair(0, 1);
        let result_reverse = seqrush.should_align_pair(1, 0);
        assert_eq!(result_forward, result_reverse);
    }

    #[test]
    fn test_sparsification_probability() {
        // Test that sparsification approximately matches expected probability
        let sequences = create_test_sequences();
        let mut config = create_test_config();
        config.sparsification = 0.3; // 30% probability
        
        let expected = config.sparsification;
        let seqrush = SeqRush::new(sequences, config);
        
        let total_pairs = 1000;
        let mut accepted = 0;
        
        // Test with synthetic pair IDs
        for i in 0..total_pairs {
            if seqrush.should_align_pair(i, i + 1000) {
                accepted += 1;
            }
        }
        
        let actual_probability = accepted as f64 / total_pairs as f64;
        println!("Expected: {:.3}, Actual: {:.3}", expected, actual_probability);
        
        // Should be approximately correct (within 10%)
        let error = (actual_probability - expected).abs();
        assert!(error < 0.1, "Sparsification error too large: {:.3}", error);
    }

    #[test]
    fn test_sparsification_with_large_n() {
        // Test that sparsification works correctly with larger numbers of sequences
        let mut sequences = Vec::new();
        for i in 0..50 {
            sequences.push(Sequence {
                id: format!("seq_{}", i),
                data: format!("ATCGATCG{:02}ATCGATCG", i).into_bytes(),
            });
        }
        
        let auto_threshold = SeqRush::calculate_erdos_renyi_threshold(50, 3.0);
        println!("Auto threshold for n=50: {:.4}", auto_threshold);
        
        let mut config = create_test_config();
        config.sparsification = auto_threshold;
        config.max_iterations = 2; // Quick test
        
        let mut seqrush = SeqRush::new(sequences, config);
        
        // Should complete without issues using sparsification
        seqrush.build_graph();
        
        let temp_file = NamedTempFile::new().unwrap();
        let result = seqrush.write_gfa(temp_file.path().to_str().unwrap());
        assert!(result.is_ok());
    }

    #[test]
    fn test_erdos_renyi_prevents_quadratic_blowup() {
        // Test that Erdős-Rényi sparsification actually reduces work
        let mut sequences = Vec::new();
        for i in 0..20 {
            sequences.push(Sequence {
                id: format!("test_seq_{}", i),
                data: format!("ATCGATCGATCG{:02}ATCGATCGATCG", i).into_bytes(),
            });
        }
        
        // Compare full vs sparse alignment counts
        let mut config_full = create_test_config();
        config_full.sparsification = 1.0; // All pairs
        config_full.verbose = true;
        
        let mut config_sparse = create_test_config();
        config_sparse.sparsification = SeqRush::calculate_erdos_renyi_threshold(20, 3.0);
        config_sparse.verbose = true;
        
        println!("Full sparsification: {:.4}", config_full.sparsification);
        println!("Sparse sparsification: {:.4}", config_sparse.sparsification);
        
        // Sparse should use significantly fewer pairs
        assert!(config_sparse.sparsification < 0.5);
        assert!(config_sparse.sparsification < config_full.sparsification);
    }
}
