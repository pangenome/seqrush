use handlegraph::handle::{Handle, NodeId, Edge};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::mutablehandlegraph::MutableHandleGraph;
use handlegraph::pathhandlegraph::embedded_paths::{GraphPaths, MutableGraphPaths};
use handlegraph::hashgraph::HashGraph;

use crate::seqrush::SeqRush;
use crate::pos::{make_pos, offset, is_rev};
use std::collections::HashMap;

/// Build a handlegraph from union-find results
impl SeqRush {
    pub fn build_handlegraph(&self, verbose: bool) -> Result<HashGraph, Box<dyn std::error::Error>> {
        let mut graph = HashGraph::new();
        
        // Track which union representatives we've seen and their node IDs
        let mut union_to_node: HashMap<u64, NodeId> = HashMap::new();
        let mut next_node_id = 1;
        
        // First pass: create all nodes
        if verbose {
            eprintln!("Creating nodes from union-find components...");
        }
        
        // Collect unique union representatives
        let mut unique_unions = std::collections::HashSet::new();
        for seq in &self.sequences {
            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                let pos_fwd = make_pos(global_pos, false);
                let union_fwd = self.union_find.find(pos_fwd);
                unique_unions.insert(union_fwd);
            }
        }
        
        if verbose {
            eprintln!("Found {} unique union components", unique_unions.len());
        }
        
        // Create a node for each union component
        for &union_rep in &unique_unions {
            // Find the base for this union
            let base = if let Some(source_seq) = self.find_sequence_for_position_hg(union_rep) {
                let local_offset = offset(union_rep) - source_seq.offset;
                source_seq.data[local_offset]
            } else {
                b'N' // Fallback
            };
            
            let node_id = NodeId::from(next_node_id);
            next_node_id += 1;
            
            graph.create_handle(&[base], node_id);
            union_to_node.insert(union_rep, node_id);
        }
        
        if verbose {
            eprintln!("Created {} nodes", graph.node_count());
        }
        
        // Second pass: build paths
        if verbose {
            eprintln!("Building paths through the graph...");
        }
        
        for (seq_idx, seq) in self.sequences.iter().enumerate() {
            let path_name = seq.id.as_bytes().to_vec();
            let path = graph.create_path(&path_name, false).unwrap();
            
            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                
                // For each position, find which union it belongs to
                let pos_fwd = make_pos(global_pos, false);
                let pos_rev = make_pos(global_pos, true);
                
                let union_fwd = self.union_find.find(pos_fwd);
                let union_rev = self.union_find.find(pos_rev);
                
                // Determine which orientation to use
                let (union_rep, is_reverse) = if union_to_node.contains_key(&union_fwd) {
                    (union_fwd, false)
                } else if union_to_node.contains_key(&union_rev) {
                    (union_rev, true)
                } else {
                    // Check if either is connected to an existing union
                    let mut found = None;
                    for &existing_union in union_to_node.keys() {
                        if self.union_find.same(union_fwd, existing_union) {
                            found = Some((existing_union, false));
                            break;
                        }
                        if self.union_find.same(union_rev, existing_union) {
                            found = Some((existing_union, true));
                            break;
                        }
                    }
                    found.unwrap_or((union_fwd, false))
                };
                
                // Get the node ID
                let node_id = union_to_node[&union_rep];
                let handle = Handle::new(node_id, is_reverse);
                
                // Append to path
                graph.append_step(&path, handle);
            }
        }
        
        if verbose {
            eprintln!("Created {} paths", graph.path_count());
        }
        
        // Third pass: create edges from path traversals
        if verbose {
            eprintln!("Creating edges from path traversals...");
        }
        
        let mut edges_added = 0;
        for path_id in graph.paths_iter() {
            let mut prev_handle: Option<Handle> = None;
            
            for handle in graph.steps_iter(&path_id) {
                if let Some(prev) = prev_handle {
                    // Create edge from prev to current
                    graph.create_edge(&Edge(prev, handle));
                    edges_added += 1;
                }
                prev_handle = Some(handle);
            }
        }
        
        if verbose {
            eprintln!("Added {} edges", edges_added);
            eprintln!("Final graph: {} nodes, {} edges, {} paths",
                     graph.node_count(), graph.edge_count(), graph.path_count());
        }
        
        Ok(graph)
    }
    
    /// Find which sequence contains a given position (handlegraph version)
    fn find_sequence_for_position_hg(&self, pos: u64) -> Option<&crate::seqrush::Sequence> {
        let offset_val = offset(pos);
        self.sequences.iter().find(|seq| {
            offset_val >= seq.offset && offset_val < seq.offset + seq.data.len()
        })
    }
}