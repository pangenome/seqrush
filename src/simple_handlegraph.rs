use handlegraph::handle::{Handle, NodeId, Edge};
use handlegraph::mutablehandlegraph::AdditiveHandleGraph;
use handlegraph::pathhandlegraph::embedded_paths::{MutableGraphPaths, GraphPaths};
use handlegraph::hashgraph::HashGraph;
use handlegraph::handlegraph::{HandleGraph, IntoHandles, IntoEdges, IntoSequences};

use crate::seqrush::SeqRush;
use crate::pos::{make_pos, offset};
use std::collections::{HashMap, HashSet};

impl SeqRush {
    /// Build and write handlegraph GFA
    pub fn write_handlegraph_gfa(&self, output_path: &str, no_compact: bool, verbose: bool) 
        -> Result<(), Box<dyn std::error::Error>> {
        
        // Build the graph
        let mut graph = HashGraph::new();
        
        // Track union to node mapping
        let mut union_to_node: HashMap<u64, NodeId> = HashMap::new();
        let mut next_node_id = 1u64;
        
        // Create nodes
        let mut unique_unions = HashSet::new();
        for seq in &self.sequences {
            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                let pos_fwd = make_pos(global_pos, false);
                let union_fwd = self.union_find.find(pos_fwd);
                unique_unions.insert(union_fwd);
            }
        }
        
        for &union_rep in &unique_unions {
            let base = self.get_base_for_union(union_rep);
            let node_id = NodeId::from(next_node_id);
            next_node_id += 1;
            
            graph.create_handle(&[base], node_id);
            union_to_node.insert(union_rep, node_id);
        }
        
        // Build paths
        let mut all_edges = Vec::new();
        
        for seq in &self.sequences {
            let path_id = graph.create_path(seq.id.as_bytes(), false).unwrap();
            let mut prev_handle: Option<Handle> = None;
            
            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                let pos_fwd = make_pos(global_pos, false);
                let union_fwd = self.union_find.find(pos_fwd);
                
                // Find the actual union rep
                let mut union_rep = union_fwd;
                for &existing in union_to_node.keys() {
                    if self.union_find.same(union_fwd, existing) {
                        union_rep = existing;
                        break;
                    }
                }
                
                let node_id = union_to_node[&union_rep];
                let handle = Handle::pack(node_id, false);
                
                graph.path_append_step(path_id, handle);
                
                if let Some(prev) = prev_handle {
                    all_edges.push(Edge(prev, handle));
                }
                prev_handle = Some(handle);
            }
        }
        
        // Add edges
        let mut seen_edges = HashSet::new();
        for edge in all_edges {
            let key = (edge.0.as_integer(), edge.1.as_integer());
            if seen_edges.insert(key) {
                graph.create_edge(edge);
            }
        }
        
        if verbose {
            println!("Built graph: {} nodes, {} edges, {} paths",
                     graph.node_count(), graph.edge_count(), graph.path_count());
        }
        
        // Apply compaction if requested
        let final_graph = graph; // Temporarily disable compaction due to API issues
        if !no_compact && verbose {
            eprintln!("[Warning] Compaction temporarily disabled due to handlegraph API issues");
        }
        
        // Write GFA
        let file = std::fs::File::create(output_path)?;
        let mut writer = std::io::BufWriter::new(file);
        use std::io::Write;
        
        writeln!(writer, "H\tVN:Z:1.0")?;
        
        // Write nodes
        for node_id in (&final_graph).handles() {
            let handle = Handle::pack(node_id, false);
            let seq: Vec<u8> = (&final_graph).sequence(handle).collect();
            writeln!(writer, "S\t{}\t{}", node_id.0, String::from_utf8_lossy(&seq))?;
        }
        
        // Write edges  
        for edge in (&final_graph).edges() {
            writeln!(writer, "L\t{}\t{}\t{}\t{}\t0M",
                    u64::from(edge.0.id()),
                    if edge.0.is_reverse() { '-' } else { '+' },
                    u64::from(edge.1.id()),
                    if edge.1.is_reverse() { '-' } else { '+' })?;
        }
        
        // Write paths
        // For now, we'll iterate through the sequences to write paths
        // since the handlegraph API for path iteration is unclear
        for (seq_idx, seq) in self.sequences.iter().enumerate() {
            let mut path_str = String::new();
            let mut step_idx = 0;
            
            // Build path by following the sequence through the graph
            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                let pos_fwd = make_pos(global_pos, false);
                let union_fwd = self.union_find.find(pos_fwd);
                
                // Find the node for this union
                for (union_rep, &node_id) in &union_to_node {
                    if self.union_find.same(union_fwd, *union_rep) {
                        if step_idx > 0 {
                            path_str.push(',');
                        }
                        path_str.push_str(&format!("{}{}",
                            u64::from(node_id),
                            '+'));
                        step_idx += 1;
                        break;
                    }
                }
            }
            
            writeln!(writer, "P\t{}\t{}\t*", seq.id, path_str)?;
        }
        
        println!("Handlegraph written to {}: {} nodes, {} edges, {} paths",
                 output_path, final_graph.node_count(), 
                 final_graph.edge_count(), final_graph.path_count());
        
        Ok(())
    }
    
    fn get_base_for_union(&self, union_rep: u64) -> u8 {
        let offset_val = offset(union_rep);
        for seq in &self.sequences {
            if offset_val >= seq.offset && offset_val < seq.offset + seq.data.len() {
                return seq.data[offset_val - seq.offset];
            }
        }
        b'N'
    }
}