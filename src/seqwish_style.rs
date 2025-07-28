// Implement seqwish-style graph construction

use crate::pos::{Pos, make_pos};
use crate::bidirected_union_find::BidirectedUnionFind;
use std::collections::{HashMap, HashSet, BinaryHeap};
use std::cmp::Ordering;

#[derive(Clone, Eq, PartialEq)]
struct DSU {
    parent: Pos,
    size: usize,
}

impl Ord for DSU {
    fn cmp(&self, other: &Self) -> Ordering {
        // Order by size (descending) then by parent position
        other.size.cmp(&self.size)
            .then_with(|| self.parent.cmp(&other.parent))
    }
}

impl PartialOrd for DSU {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Build a graph sequence in the style of seqwish
pub fn build_graph_sequence(
    sequences: &[(String, Vec<u8>, usize)], // (id, data, offset)
    union_find: &BidirectedUnionFind,
) -> (Vec<u8>, Vec<Pos>) {
    // Step 1: Collect all disjoint sets
    let mut dsets: HashMap<Pos, HashSet<Pos>> = HashMap::new();
    let total_positions: usize = sequences.iter().map(|(_, data, _)| data.len()).sum();
    
    // For each position, find its representative and group
    for (_, seq_data, offset) in sequences {
        for i in 0..seq_data.len() {
            let pos = make_pos(offset + i, false);
            let rep = union_find.find(pos);
            dsets.entry(rep).or_insert_with(HashSet::new).insert(pos);
        }
    }
    
    // Step 2: Sort disjoint sets by size (largest first) for better layout
    let mut sorted_dsets: BinaryHeap<DSU> = BinaryHeap::new();
    for (parent, positions) in dsets {
        sorted_dsets.push(DSU {
            parent,
            size: positions.len(),
        });
    }
    
    // Step 3: Build the graph sequence
    let mut graph_sequence = Vec::new();
    let mut position_to_graph_offset = Vec::new();
    position_to_graph_offset.resize(total_positions * 2, 0); // *2 for both orientations
    
    while let Some(dsu) = sorted_dsets.pop() {
        let positions = &dsets[&dsu.parent];
        
        // Get the consensus base for this component
        let mut base_counts: HashMap<u8, usize> = HashMap::new();
        for &pos in positions {
            let seq_idx = find_sequence_for_position(pos, sequences);
            if let Some((_, data, offset)) = seq_idx {
                let local_pos = (pos >> 1) - offset;
                if local_pos < data.len() {
                    *base_counts.entry(data[local_pos]).or_insert(0) += 1;
                }
            }
        }
        
        // Use the most common base
        if let Some((&base, _)) = base_counts.iter().max_by_key(|(_, &count)| count) {
            let graph_offset = graph_sequence.len();
            graph_sequence.push(base);
            
            // Record where each position maps in the graph sequence
            for &pos in positions {
                position_to_graph_offset[pos as usize] = graph_offset;
            }
        }
    }
    
    (graph_sequence, position_to_graph_offset)
}

fn find_sequence_for_position(pos: Pos, sequences: &[(String, Vec<u8>, usize)]) -> Option<&(String, Vec<u8>, usize)> {
    let offset = pos >> 1; // Remove orientation bit
    sequences.iter().find(|(_, data, seq_offset)| {
        offset >= *seq_offset && offset < seq_offset + data.len()
    })
}