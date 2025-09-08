use crate::bidirected_graph::{BiNode, BiPath, Handle};
use crate::bidirected_ops::BidirectedGraph;
use std::collections::{HashMap, HashSet};

/// Represents an alignment range between sequences
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AlignmentRange {
    pub seq1_start: usize,
    pub seq1_end: usize,
    pub seq2_start: usize,
    pub seq2_end: usize,
    pub seq2_is_rc: bool,
}

/// Builds a graph using alignment ranges (like seqwish)
pub struct RangeBasedGraphBuilder {
    /// All alignment ranges from PAF
    ranges: Vec<AlignmentRange>,
    /// Sequence data
    sequences: Vec<(String, Vec<u8>)>,
}

impl RangeBasedGraphBuilder {
    pub fn new() -> Self {
        Self {
            ranges: Vec::new(),
            sequences: Vec::new(),
        }
    }

    pub fn add_sequence(&mut self, name: String, data: Vec<u8>) {
        self.sequences.push((name, data));
    }

    pub fn add_alignment_range(&mut self, range: AlignmentRange) {
        self.ranges.push(range);
    }

    /// Build graph like seqwish does
    pub fn build_graph(
        &self,
        verbose: bool,
    ) -> Result<BidirectedGraph, Box<dyn std::error::Error>> {
        if verbose {
            eprintln!("Building graph from {} alignment ranges", self.ranges.len());
        }

        // Step 1: Create a "graph sequence" by laying out all sequences
        let mut graph_seq = Vec::new();
        let mut seq_offsets = Vec::new();
        let mut total_offset = 0;

        for (_, data) in &self.sequences {
            seq_offsets.push(total_offset);
            graph_seq.extend_from_slice(data);
            total_offset += data.len();
        }

        // Step 2: Build interval trees for alignments (like seqwish's node_iitree and path_iitree)
        // For now, we'll use a simpler approach with HashMaps
        let mut position_to_ranges: HashMap<usize, Vec<usize>> = HashMap::new();

        // Add self-alignments first (each sequence maps to itself)
        let mut all_ranges = self.ranges.clone();
        for (seq_idx, (_, data)) in self.sequences.iter().enumerate() {
            let offset = seq_offsets[seq_idx];
            // Add one self-alignment for the entire sequence
            all_ranges.push(AlignmentRange {
                seq1_start: offset,
                seq1_end: offset + data.len(),
                seq2_start: offset,
                seq2_end: offset + data.len(),
                seq2_is_rc: false,
            });
        }

        // Step 3: Mark node boundaries (like seqwish's compact_nodes)
        let mut node_boundaries = HashSet::new();
        node_boundaries.insert(0); // First position is always a boundary

        // For each alignment range, mark start and end as boundaries
        for (range_idx, range) in all_ranges.iter().enumerate() {
            // Mark boundaries in graph sequence
            node_boundaries.insert(range.seq1_start);
            node_boundaries.insert(range.seq1_end);

            // Track which positions are covered by this range
            for pos in range.seq1_start..range.seq1_end {
                position_to_ranges.entry(pos).or_default().push(range_idx);
            }
        }

        node_boundaries.insert(graph_seq.len()); // Last position is always a boundary

        // Convert boundaries to sorted vector
        let mut boundaries: Vec<usize> = node_boundaries.into_iter().collect();
        boundaries.sort();

        if verbose {
            eprintln!("Found {} node boundaries", boundaries.len());
        }

        // Step 4: Create nodes from boundary segments
        let mut graph = BidirectedGraph::new();
        let mut position_to_node: HashMap<usize, (usize, bool)> = HashMap::new();
        let mut node_id = 1;

        for i in 0..boundaries.len() - 1 {
            let start = boundaries[i];
            let end = boundaries[i + 1];

            if start >= end {
                continue;
            }

            // Extract sequence for this node
            let node_seq = graph_seq[start..end].to_vec();

            // Create the node
            let bi_node = BiNode {
                id: node_id,
                sequence: node_seq,
                rank: Some(node_id as u64),
            };
            graph.nodes.insert(node_id, bi_node);

            // Map positions to this node
            for pos in start..end {
                position_to_node.insert(pos, (node_id, false));
            }

            if verbose && node_id <= 5 {
                eprintln!(
                    "Created node {} spanning positions {}..{} with sequence '{}'",
                    node_id,
                    start,
                    end,
                    String::from_utf8_lossy(&graph.nodes[&node_id].sequence)
                );
            }

            node_id += 1;
        }

        // Step 5: Build paths through the graph
        for (seq_idx, (seq_name, seq_data)) in self.sequences.iter().enumerate() {
            let offset = seq_offsets[seq_idx];
            let mut path_handles = Vec::new();
            let mut last_node_id = None;

            for i in 0..seq_data.len() {
                let pos = offset + i;
                if let Some(&(node_id, is_rev)) = position_to_node.get(&pos) {
                    // Only add to path if it's a different node than the last one
                    if last_node_id != Some(node_id) {
                        path_handles.push(Handle::new(node_id, is_rev));
                        last_node_id = Some(node_id);
                    }
                }
            }

            let path = BiPath {
                name: seq_name.clone(),
                steps: path_handles,
            };
            graph.paths.push(path);
        }

        // Step 6: Add edges between consecutive nodes in paths
        let edges_to_add: Vec<(Handle, Handle)> = graph
            .paths
            .iter()
            .flat_map(|path| {
                if path.steps.len() > 1 {
                    (0..path.steps.len() - 1)
                        .map(|i| (path.steps[i], path.steps[i + 1]))
                        .collect::<Vec<_>>()
                } else {
                    Vec::new()
                }
            })
            .collect();

        for (from, to) in edges_to_add {
            graph.add_edge(from, to);
        }

        if verbose {
            eprintln!("Built graph with {} nodes", graph.nodes.len());
        }

        Ok(graph)
    }
}
