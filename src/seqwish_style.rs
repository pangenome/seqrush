use crate::bidirected_graph::{BiEdge, BiNode, BiPath, Handle};
use crate::bidirected_ops::BidirectedGraph;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Compute reverse complement of a base
fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => base,
    }
}

/// Union-Find data structure for transitive closure
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    fn unite(&mut self, x: usize, y: usize) {
        let root_x = self.find(x);
        let root_y = self.find(y);

        if root_x != root_y {
            if self.rank[root_x] < self.rank[root_y] {
                self.parent[root_x] = root_y;
            } else if self.rank[root_x] > self.rank[root_y] {
                self.parent[root_y] = root_x;
            } else {
                self.parent[root_y] = root_x;
                self.rank[root_x] += 1;
            }
        }
    }
}

/// Seqwish-style graph builder
pub struct SeqwishStyleBuilder {
    sequences: Vec<(String, Vec<u8>)>,
    seq_offsets: Vec<usize>,
    total_length: usize,
}

impl SeqwishStyleBuilder {
    pub fn new() -> Self {
        Self {
            sequences: Vec::new(),
            seq_offsets: Vec::new(),
            total_length: 0,
        }
    }

    pub fn add_sequence(&mut self, name: String, data: Vec<u8>) {
        self.seq_offsets.push(self.total_length);
        self.total_length += data.len();
        self.sequences.push((name, data));
    }

    /// Build graph from PAF alignments
    pub fn build_from_paf(
        self,
        paf_path: &str,
        min_match_len: usize,
        verbose: bool,
    ) -> Result<BidirectedGraph, Box<dyn std::error::Error>> {
        // Step 1: Build union-find from alignments
        let mut uf = UnionFind::new(self.total_length);

        // CRITICAL: Unite forward and reverse orientations of each position
        // This ensures we work with a single "strand" where alignments can happen in both orientations
        for i in 0..self.total_length {
            // In this simple union-find, we need to double the size to handle orientations
            // Position i forward = i * 2, Position i reverse = i * 2 + 1
            uf.unite(i * 2, i * 2 + 1);
        }

        let alignment_count =
            self.process_alignments_paf(paf_path, min_match_len, &mut uf, verbose)?;

        if verbose {
            eprintln!("[seqwish-style] Processed {} alignments", alignment_count);
        }

        // Step 2: Build graph from union-find components
        let graph = self.build_graph_from_union_find(&mut uf, verbose)?;

        Ok(graph)
    }

    /// Process PAF alignments and build union-find
    fn process_alignments_paf(
        &self,
        paf_path: &str,
        min_match_len: usize,
        uf: &mut UnionFind,
        _verbose: bool,
    ) -> Result<usize, Box<dyn std::error::Error>> {
        let file = File::open(paf_path)?;
        let reader = BufReader::new(file);
        let mut alignment_count = 0;
        let mut total_unites = 0;

        for line in reader.lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 12 {
                continue;
            }

            let query_name = fields[0];
            let query_len: usize = fields[1].parse()?;
            let query_start: usize = fields[2].parse()?;
            let query_end: usize = fields[3].parse()?;
            let query_strand = fields[4];
            let target_name = fields[5];
            let _target_len: usize = fields[6].parse()?;
            let target_start: usize = fields[7].parse()?;
            let target_end: usize = fields[8].parse()?;

            // DON'T skip self-alignments - they establish the backbone!
            // if query_name == target_name && query_start == target_start && query_end == target_end {
            //     continue;
            // }

            // Debug first few alignments
            if alignment_count < 3 && _verbose {
                eprintln!(
                    "DEBUG: Alignment {} -> {}: {}..{} to {}..{} (strand {})",
                    query_name,
                    target_name,
                    query_start,
                    query_end,
                    target_start,
                    target_end,
                    query_strand
                );

                // Also debug the offsets
                let query_offset = self
                    .sequences
                    .iter()
                    .position(|(name, _)| name == query_name)
                    .map(|i| self.seq_offsets[i])
                    .unwrap_or(0);
                let target_offset = self
                    .sequences
                    .iter()
                    .position(|(name, _)| name == target_name)
                    .map(|i| self.seq_offsets[i])
                    .unwrap_or(0);
                eprintln!(
                    "DEBUG: query_offset={}, target_offset={}",
                    query_offset, target_offset
                );
            }

            // Find CIGAR string
            let cigar = fields
                .iter()
                .find(|f| f.starts_with("cg:Z:"))
                .map(|f| &f[5..])
                .unwrap_or("");

            if !cigar.is_empty() {
                // Debug CIGAR
                if alignment_count < 3 && _verbose {
                    eprintln!("DEBUG: CIGAR: {}", cigar);
                }
                // Find sequence offsets - match by prefix since FASTA has full description
                let query_offset = self
                    .sequences
                    .iter()
                    .position(|(name, _)| name.starts_with(query_name) || name == query_name)
                    .map(|i| self.seq_offsets[i])
                    .unwrap_or(0);
                let target_offset = self
                    .sequences
                    .iter()
                    .position(|(name, _)| name.starts_with(target_name) || name == target_name)
                    .map(|i| self.seq_offsets[i])
                    .unwrap_or(0);

                // Process CIGAR alignment
                let unites = self.process_cigar_alignment(
                    cigar,
                    query_offset + query_start,
                    target_offset + target_start,
                    query_strand == "-",
                    query_len,
                    min_match_len,
                    uf,
                );
                total_unites += unites;
                alignment_count += 1;
            }
        }

        if _verbose {
            eprintln!("[seqwish-style] Total unite operations: {}", total_unites);
        }
        Ok(alignment_count)
    }

    /// Process CIGAR string to unite matching positions
    fn process_cigar_alignment(
        &self,
        cigar: &str,
        query_start: usize,
        target_start: usize,
        query_is_rev: bool,
        query_len: usize,
        min_match_len: usize,
        uf: &mut UnionFind,
    ) -> usize {
        let mut query_pos = 0;
        let mut target_pos = 0;
        let mut i = 0;
        let cigar_bytes = cigar.as_bytes();
        let mut unite_count = 0;

        while i < cigar_bytes.len() {
            // Parse count
            let mut count = 0;
            while i < cigar_bytes.len() && cigar_bytes[i].is_ascii_digit() {
                count = count * 10 + (cigar_bytes[i] - b'0') as usize;
                i += 1;
            }

            if i >= cigar_bytes.len() {
                break;
            }

            let op = cigar_bytes[i] as char;
            i += 1;

            // Debug: track operation types
            if unite_count == 0 && (op == '=' || op == 'M') {
                eprintln!(
                    "DEBUG: CIGAR op '{}' with count {}, min_match_len={}",
                    op, count, min_match_len
                );
            }

            match op {
                '=' => {
                    // Exact match - unite the positions
                    if count >= min_match_len {
                        if unite_count == 0 {
                            eprintln!(
                                "DEBUG: Inside '=' case, count={}, min_match_len={}",
                                count, min_match_len
                            );
                        }
                        for j in 0..count {
                            let q_pos = if query_is_rev {
                                query_start + query_len - query_pos - count + j
                            } else {
                                query_start + query_pos + j
                            };
                            let t_pos = target_start + target_pos + j;

                            // Debug first few unites
                            if unite_count < 5 {
                                eprintln!(
                                    "DEBUG: Unite #{}: q_pos={} t_pos={}",
                                    unite_count, q_pos, t_pos
                                );
                            }

                            // Unite query and target positions
                            uf.unite(q_pos, t_pos);
                            unite_count += 1;
                        }
                    }
                    query_pos += count;
                    target_pos += count;
                }
                'M' => {
                    // For 'M' operations, we should check if bases match
                    // For now, let's treat all M as matches
                    if count >= min_match_len {
                        for j in 0..count {
                            let q_pos = if query_is_rev {
                                query_start + query_len - query_pos - count + j
                            } else {
                                query_start + query_pos + j
                            };
                            let t_pos = target_start + target_pos + j;

                            // Check if bases match
                            let q_base = self.get_base_at_position(q_pos);
                            let t_base = self.get_base_at_position(t_pos);

                            if q_base == t_base || (query_is_rev && q_base == complement(t_base)) {
                                uf.unite(q_pos, t_pos);
                            }
                        }
                    }
                    query_pos += count;
                    target_pos += count;
                }
                'X' => {
                    // Mismatch - skip
                    query_pos += count;
                    target_pos += count;
                }
                'I' => {
                    // Insertion in query
                    query_pos += count;
                }
                'D' => {
                    // Deletion in query
                    target_pos += count;
                }
                _ => {
                    // Skip unknown operations
                }
            }
        }
        unite_count
    }

    /// Build graph from union-find components
    fn build_graph_from_union_find(
        &self,
        uf: &mut UnionFind,
        verbose: bool,
    ) -> Result<BidirectedGraph, Box<dyn std::error::Error>> {
        let mut graph = BidirectedGraph::new();

        // Step 1: Collect all components and their positions
        let mut components: HashMap<usize, Vec<usize>> = HashMap::new();
        for pos in 0..self.total_length {
            let root = uf.find(pos);
            components.entry(root).or_default().push(pos);
        }

        if verbose {
            eprintln!(
                "[seqwish-style] Found {} union-find components",
                components.len()
            );
            let single_pos_components = components.values().filter(|v| v.len() == 1).count();
            let multi_pos_components = components.values().filter(|v| v.len() > 1).count();
            eprintln!(
                "[seqwish-style]   {} single-position components",
                single_pos_components
            );
            eprintln!(
                "[seqwish-style]   {} multi-position components",
                multi_pos_components
            );

            // Debug: check total positions covered
            let total_positions: usize = components.values().map(|v| v.len()).sum();
            eprintln!(
                "[seqwish-style] Total positions in components: {} (out of {})",
                total_positions, self.total_length
            );
        }

        // Step 2: Sort components by their minimum position
        let mut sorted_components: Vec<(usize, Vec<usize>)> = components.into_iter().collect();
        sorted_components.sort_by_key(|(_, positions)| *positions.iter().min().unwrap());

        // Step 3: Create one node per component (single-base nodes for now)
        let mut position_to_node: HashMap<usize, Handle> = HashMap::new();
        let mut node_id = 1u64; // Start from 1

        for (root, positions) in &sorted_components {
            // Get the base for this component (they should all be the same)
            let base = self.get_base_at_position(positions[0]);

            // Create node
            let node = BiNode {
                id: node_id as usize,
                sequence: vec![base],
                rank: Some(0),
            };
            let id = node_id as usize;
            // Ensure the vec is large enough
            if id >= graph.nodes.len() {
                graph.nodes.resize(id + 1, None);
            }
            graph.nodes[id] = Some(node);

            let handle = Handle::new(node_id as usize, false);

            // Map all positions in this component to this node
            for &pos in positions {
                position_to_node.insert(pos, handle);
            }

            node_id += 1;
        }

        if verbose {
            eprintln!("[seqwish-style] Created {} nodes", node_id - 1);
        }

        // Step 4: Create paths
        for (seq_idx, (name, data)) in self.sequences.iter().enumerate() {
            let seq_start = self.seq_offsets[seq_idx];
            let mut path = BiPath {
                name: name.clone(),
                steps: Vec::new(),
            };

            for i in 0..data.len() {
                let abs_pos = seq_start + i;
                if let Some(&handle) = position_to_node.get(&abs_pos) {
                    // Only add if different from previous step (avoid consecutive duplicates)
                    if path.steps.is_empty() || path.steps.last() != Some(&handle) {
                        path.steps.push(handle);
                    }
                }
            }

            graph.paths.push(path);
        }

        // Step 5: Create edges based on path connectivity
        let mut edge_set = std::collections::HashSet::new();
        for path in &graph.paths {
            for window in path.steps.windows(2) {
                let from = window[0];
                let to = window[1];
                if from != to {
                    // Don't create self-loops
                    let edge = BiEdge { from, to };
                    edge_set.insert(edge);
                }
            }
        }

        graph.edges = edge_set.into_iter().collect();

        if verbose {
            eprintln!("[seqwish-style] Graph construction complete");
        }

        Ok(graph)
    }

    /// Get the base at a given absolute position
    fn get_base_at_position(&self, pos: usize) -> u8 {
        let mut offset = 0;
        for (_, seq) in &self.sequences {
            if pos < offset + seq.len() {
                return seq[pos - offset];
            }
            offset += seq.len();
        }
        b'N' // Should never reach here
    }
}
