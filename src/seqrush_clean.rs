use crate::graph_ops::{Edge, Graph, Node};
use allwave::{AlignmentParams, AllPairIterator, SparsificationStrategy};
use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::Arc;
use std::sync::Mutex;
use uf_rush::UFRush;

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

    /// Output GFA file
    #[arg(short, long, default_value = "output.gfa")]
    pub output: String,

    /// Number of threads
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,

    /// Minimum match length
    #[arg(short = 'k', long, default_value_t = 0)]
    pub min_match_length: usize,

    /// Alignment scores (match,mismatch,gap_open,gap_extend[,gap2_open,gap2_extend])
    #[arg(short = 'S', long, default_value = "0,5,8,2,24,1")]
    pub scores: String,

    /// Orientation check alignment scores for fast edit distance (match,mismatch,gap_open,gap_extend)
    #[arg(long, default_value = "0,1,1,1")]
    pub orientation_scores: String,

    /// Maximum divergence threshold (0.0-1.0, e.g., 0.1 = 10% divergence)
    #[arg(short = 'd', long)]
    pub max_divergence: Option<f64>,

    /// Verbose output
    #[arg(short, long)]
    pub verbose: bool,

    /// Output alignments to PAF file
    #[arg(long = "output-alignments")]
    pub output_alignments: Option<String>,
}

#[derive(Clone, Debug)]
pub struct Sequence {
    pub id: String,
    pub data: Vec<u8>,
    pub offset: usize,
}

pub struct SeqRush {
    pub sequences: Vec<Sequence>,
    pub total_length: usize,
    pub union_find: Arc<UFRush>,
}

impl SeqRush {
    pub fn new(sequences: Vec<Sequence>, total_length: usize) -> Self {
        let union_find = Arc::new(UFRush::new(total_length));

        Self {
            sequences,
            total_length,
            union_find,
        }
    }

    pub fn build_graph(&mut self, args: &Args) {
        println!("Loaded {} sequences", self.sequences.len());
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
        // Parse alignment parameters
        let params = parse_alignment_params(&args.scores, args.max_divergence);
        let orientation_params = parse_alignment_params(&args.orientation_scores, None);

        // Convert sequences to allwave format
        let allwave_sequences: Vec<allwave::Sequence> = self
            .sequences
            .iter()
            .map(|s| allwave::Sequence {
                id: s.id.clone(),
                seq: s.data.clone(),
            })
            .collect();

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

        // Create aligner excluding self-alignments
        let aligner = AllPairIterator::with_options(
            &allwave_sequences,
            params,
            true, // exclude self
            SparsificationStrategy::None,
        )
        .with_orientation_params(orientation_params);

        let n = self.sequences.len();
        println!(
            "Total sequence pairs: {} (sparsification: None)",
            n * (n - 1)
        );

        // Process alignments sequentially for debugging
        let alignments: Vec<_> = aligner.collect();
        for alignment in alignments {
            // Write to PAF if requested
            if let Some(writer) = &paf_writer {
                let paf_record = allwave::alignment_to_paf(&alignment, &allwave_sequences);
                if let Ok(mut w) = writer.lock() {
                    let _ = writeln!(w, "{}", paf_record);
                }
            }

            // Process alignment for union-find
            let cigar = allwave::cigar_bytes_to_string(&alignment.cigar_bytes);
            self.process_alignment(
                &cigar,
                &self.sequences[alignment.query_idx],
                &self.sequences[alignment.target_idx],
                args.min_match_length,
                alignment.is_reverse,
            );
        }

        // Force all union-find operations to complete
        std::sync::atomic::fence(std::sync::atomic::Ordering::SeqCst);

        // Flush PAF writer
        if let Some(writer) = &paf_writer {
            if let Ok(mut w) = writer.lock() {
                let _ = w.flush();
            }
        }
    }

    fn process_alignment(
        &self,
        cigar: &str,
        seq1: &Sequence,
        seq2: &Sequence,
        min_match_len: usize,
        seq2_is_rc: bool,
    ) {
        let mut pos1 = 0;
        let mut pos2 = 0;
        let mut count = 0;

        // Track matches
        let mut in_match = false;
        let mut match_start1 = 0;
        let mut match_start2 = 0;
        let mut match_len = 0;

        // Get seq2 data (potentially reverse complemented)
        let seq2_data = if seq2_is_rc {
            reverse_complement(&seq2.data)
        } else {
            seq2.data.clone()
        };

        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                count = count * 10 + (ch as usize - '0' as usize);
            } else {
                if count == 0 {
                    count = 1;
                }

                match ch {
                    'M' | '=' => {
                        // Check each position in this M operation
                        for k in 0..count {
                            if pos1 + k < seq1.data.len() && pos2 + k < seq2_data.len() {
                                if seq1.data[pos1 + k] == seq2_data[pos2 + k] {
                                    // Bases match
                                    if !in_match {
                                        in_match = true;
                                        match_start1 = pos1 + k;
                                        match_start2 = pos2 + k;
                                        match_len = 1;
                                    } else {
                                        match_len += 1;
                                    }
                                } else {
                                    // Mismatch - process any accumulated matches
                                    if in_match && match_len >= min_match_len {
                                        self.unite_match(
                                            seq1,
                                            seq2,
                                            match_start1,
                                            match_start2,
                                            match_len,
                                            seq2_is_rc,
                                        );
                                    }
                                    in_match = false;
                                    match_len = 0;
                                }
                            }
                        }
                        pos1 += count;
                        pos2 += count;
                    }
                    _ => {
                        // Not a match - process any accumulated matches
                        if in_match && match_len >= min_match_len {
                            self.unite_match(
                                seq1,
                                seq2,
                                match_start1,
                                match_start2,
                                match_len,
                                seq2_is_rc,
                            );
                        }
                        in_match = false;
                        match_len = 0;

                        match ch {
                            'X' => {
                                pos1 += count;
                                pos2 += count;
                            }
                            'I' => {
                                pos1 += count;
                            }
                            'D' => {
                                pos2 += count;
                            }
                            _ => {}
                        }
                    }
                }
                count = 0;
            }
        }

        // Process final match if any
        if in_match && match_len >= min_match_len {
            self.unite_match(
                seq1,
                seq2,
                match_start1,
                match_start2,
                match_len,
                seq2_is_rc,
            );
        }
    }

    fn unite_match(
        &self,
        seq1: &Sequence,
        seq2: &Sequence,
        start1: usize,
        start2: usize,
        length: usize,
        seq2_is_rc: bool,
    ) {
        static DEBUG_COUNT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let count = DEBUG_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);

        if count < 5 {
            println!(
                "DEBUG unite_match {}: {} -> {} start1={} start2={} len={} rc={}",
                count, seq1.id, seq2.id, start1, start2, length, seq2_is_rc
            );
        }

        for i in 0..length {
            let global_pos1 = seq1.offset + start1 + i;

            let global_pos2 = if seq2_is_rc {
                // Map RC position back to forward strand coordinates
                let rc_pos = start2 + i;
                let forward_pos = seq2.data.len() - 1 - rc_pos;
                seq2.offset + forward_pos
            } else {
                seq2.offset + start2 + i
            };

            if count < 5 && i < 3 && seq2_is_rc {
                println!("    RC mapping: rc_pos={} forward_pos={} seq2.len={} seq2.offset={} -> global_pos2={}",
                         start2 + i, seq2.data.len() - 1 - (start2 + i), seq2.data.len(), seq2.offset, global_pos2);
            }

            // Following seqwish: skip self-mappings
            if global_pos1 != global_pos2 {
                self.union_find.unite(global_pos1, global_pos2);
                if count < 5 && i < 3 {
                    println!("  united pos {} with pos {}", global_pos1, global_pos2);
                }
            }
        }
    }

    fn write_gfa(&self, args: &Args) -> Result<(), Box<dyn std::error::Error>> {
        let output_path = &args.output;

        // Build graph from union-find
        let mut graph = Graph::new();
        let mut union_to_node: HashMap<usize, usize> = HashMap::new();
        let mut next_node_id = 1;

        // First, ensure all finds are fully compressed
        let mut all_positions = Vec::new();
        for seq in &self.sequences {
            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                all_positions.push((global_pos, seq.data[i]));
            }
        }

        // Do multiple rounds of find to ensure full path compression
        for _ in 0..3 {
            for &(pos, _) in &all_positions {
                self.union_find.find(pos);
            }
        }

        // Now collect positions by their union representative
        let mut union_positions: HashMap<usize, Vec<(usize, u8)>> = HashMap::new();
        for (pos, base) in all_positions {
            let union_rep = self.union_find.find(pos);
            union_positions
                .entry(union_rep)
                .or_default()
                .push((pos, base));
        }

        if args.verbose {
            println!("Found {} union components", union_positions.len());
            // Debug: check if positions are properly united
            let pos0_rep = self.union_find.find(0);
            let pos12_rep = self.union_find.find(12);
            let pos24_rep = self.union_find.find(24);
            println!(
                "DEBUG: pos 0 rep = {}, pos 12 rep = {}, pos 24 rep = {}",
                pos0_rep, pos12_rep, pos24_rep
            );
            println!("DEBUG: Are 0 and 12 in same set? {}", pos0_rep == pos12_rep);
            println!("DEBUG: Are 0 and 24 in same set? {}", pos0_rep == pos24_rep);
        }

        // Debug the union positions
        if args.verbose {
            println!("\nDEBUG union_positions:");
            for (rep, positions) in &union_positions {
                println!("  Rep {}: {} positions", rep, positions.len());
                if positions.len() > 1 {
                    println!(
                        "    Positions: {:?}",
                        positions.iter().map(|(p, _)| p).take(5).collect::<Vec<_>>()
                    );
                }
            }
        }

        // Create nodes for each union component
        for (union_rep, positions) in &union_positions {
            let node_id = next_node_id;
            next_node_id += 1;

            // Use the first base we see (they should all be the same)
            let base = positions[0].1;

            // Map the union representative to this node
            union_to_node.insert(*union_rep, node_id);

            let node = Node {
                id: node_id,
                sequence: vec![base],
                rank: node_id as f64,
            };
            graph.nodes.insert(node_id, node);
        }

        // Second pass: build paths
        for seq in &self.sequences {
            let mut path = Vec::new();

            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                let union_rep = self.union_find.find(global_pos);
                let node_id = union_to_node[&union_rep];
                path.push(node_id);
            }

            graph.paths.push((seq.id.clone(), path));
        }

        // Build edges from paths
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

        // Write GFA
        let file = File::create(output_path)?;
        let mut writer = BufWriter::new(file);

        // Header
        writeln!(writer, "H\tVN:Z:1.0")?;

        // Segments
        for (id, node) in &graph.nodes {
            writeln!(
                writer,
                "S\t{}\t{}",
                id,
                String::from_utf8_lossy(&node.sequence)
            )?;
        }

        // Links
        for edge in &graph.edges {
            writeln!(writer, "L\t{}\t+\t{}\t+\t0M", edge.from, edge.to)?;
        }

        // Paths
        for (name, path) in &graph.paths {
            let path_str = path
                .iter()
                .map(|&id| format!("{}+", id))
                .collect::<Vec<_>>()
                .join(",");
            writeln!(writer, "P\t{}\t{}\t*", name, path_str)?;
        }

        writer.flush()?;

        println!(
            "Graph written to {}: {} nodes, {} edges",
            output_path,
            graph.nodes.len(),
            graph.edges.len()
        );

        Ok(())
    }
}

fn parse_alignment_params(scores_str: &str, max_divergence: Option<f64>) -> AlignmentParams {
    let parts: Vec<i32> = scores_str
        .split(',')
        .filter_map(|s| s.parse().ok())
        .collect();

    AlignmentParams {
        match_score: parts.get(0).copied().unwrap_or(0),
        mismatch_penalty: parts.get(1).copied().unwrap_or(5),
        gap_open: parts.get(2).copied().unwrap_or(8),
        gap_extend: parts.get(3).copied().unwrap_or(2),
        gap2_open: parts.get(4).copied(),
        gap2_extend: parts.get(5).copied(),
        max_divergence,
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
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
    let sequences = load_sequences(&args.sequences)?;
    let total_length: usize = sequences.iter().map(|s| s.data.len()).sum();

    let mut seqrush = SeqRush::new(sequences, total_length);
    seqrush.build_graph(&args);

    Ok(())
}
