use handlegraph::handle::{Handle, NodeId, Edge};
use handlegraph::mutablehandlegraph::{MutableHandleGraph, AdditiveHandleGraph};
use handlegraph::pathhandlegraph::embedded_paths::MutableGraphPaths;
use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::HandleGraph;

use crate::seqrush::SeqRush;
use crate::pos::{make_pos, offset};
use std::collections::HashMap;

impl SeqRush {
    /// Build a handlegraph from union-find results
    pub fn build_handlegraph_v2(&self, verbose: bool) -> Result<HashGraph, Box<dyn std::error::Error>> {
        let mut graph = HashGraph::new();
        
        // Track which union representatives we've seen and their node IDs
        let mut union_to_node: HashMap<u64, NodeId> = HashMap::new();
        let mut next_node_id = 1u64;
        
        if verbose {
            eprintln!("Building handlegraph from union-find results...");
        }
        
        // First, create all nodes from union components
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
        
        // Create nodes
        for &union_rep in &unique_unions {
            let base = self.get_base_for_union_v2(union_rep);
            let node_id = NodeId::from(next_node_id);
            next_node_id += 1;
            
            graph.create_handle(&[base], node_id);
            union_to_node.insert(union_rep, node_id);
        }
        
        if verbose {
            eprintln!("Created {} nodes", graph.node_count());
        }
        
        // Build paths and collect edges
        let mut edges_to_add = Vec::new();
        
        for seq in &self.sequences {
            let path_name = seq.id.as_bytes();
            let path_id = graph.create_path(path_name, false)?;
            
            let mut prev_handle: Option<Handle> = None;
            
            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                
                // Find which union this position belongs to
                let pos_fwd = make_pos(global_pos, false);
                let pos_rev = make_pos(global_pos, true);
                
                let union_fwd = self.union_find.find(pos_fwd);
                let union_rev = self.union_find.find(pos_rev);
                
                // Determine orientation
                let (union_rep, is_reverse) = if union_to_node.contains_key(&union_fwd) {
                    (union_fwd, false)
                } else if union_to_node.contains_key(&union_rev) {
                    (union_rev, true)
                } else {
                    // Find which union this is connected to
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
                
                let node_id = union_to_node[&union_rep];
                let handle = Handle::new(node_id, is_reverse);
                
                // Add to path
                graph.path_append_step(&path_id, handle);
                
                // Record edge
                if let Some(prev) = prev_handle {
                    edges_to_add.push(Edge(prev, handle));
                }
                prev_handle = Some(handle);
            }
        }
        
        // Add edges (deduplicating)
        let mut seen_edges = std::collections::HashSet::new();
        for edge in edges_to_add {
            let key = (edge.0.as_integer(), edge.1.as_integer());
            if seen_edges.insert(key) {
                graph.create_edge(&edge);
            }
        }
        
        if verbose {
            eprintln!("Final graph: {} nodes, {} edges, {} paths",
                     graph.node_count(), graph.edge_count(), graph.path_count());
        }
        
        Ok(graph)
    }
    
    /// Get the base for a union representative
    fn get_base_for_union_v2(&self, union_rep: u64) -> u8 {
        let offset_val = offset(union_rep);
        
        // Find which sequence contains this position
        for seq in &self.sequences {
            if offset_val >= seq.offset && offset_val < seq.offset + seq.data.len() {
                let local_offset = offset_val - seq.offset;
                return seq.data[local_offset];
            }
        }
        
        b'N' // Fallback
    }
}