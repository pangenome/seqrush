use crate::bidirected_graph::Handle;
use crate::bidirected_ops::BidirectedGraph;
use crate::pos::{Pos, make_pos, offset};
use crate::graph_ops::{Graph, Node, Edge};
use crate::seqrush::{SeqRush, Sequence};
use std::collections::HashMap;

impl SeqRush {
    /// Build a bidirected graph from the union-find results
    pub fn build_bidirected_graph(&self, verbose: bool) -> Result<BidirectedGraph, Box<dyn std::error::Error>> {
        let mut graph = BidirectedGraph::new();
        
        // Track which union representatives we've seen and their node IDs
        let mut union_to_node: HashMap<Pos, usize> = HashMap::new();
        let mut next_node_id = 1;
        
        // Build paths and discover nodes
        for seq in &self.sequences {
            let mut path_handles = Vec::new();
            
            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                
                // Since forward and reverse orientations are already united,
                // we just need to find the union representative for the forward orientation
                let pos_fwd = make_pos(global_pos, false);
                let union_rep = self.union_find.find(pos_fwd);
                
                // The path orientation is always forward since we're following the sequence
                let _path_orientation = false;
                
                if verbose && i < 5 {
                    eprintln!("  [BIDIRECTED] {} pos {} - union: {}", 
                        seq.id, i, union_rep);
                }
                
                // Use the union representative as the canonical form
                let canonical_union = union_rep;
                
                // Get or create node ID
                let node_id = match union_to_node.get(&canonical_union) {
                    Some(&id) => id,
                    None => {
                        let id = next_node_id;
                        next_node_id += 1;
                        union_to_node.insert(canonical_union, id);
                        
                        // Create node with the base from the canonical position
                        let base = if let Some(source_seq) = self.find_sequence_for_position(canonical_union) {
                            let local_offset = offset(canonical_union) - source_seq.offset;
                            source_seq.data[local_offset]
                        } else {
                            seq.data[i] // Fallback to current sequence's base
                        };
                        
                        graph.add_node(id, vec![base]);
                        id
                    }
                };
                
                // Add handle with forward orientation to path
                let handle = Handle::new(node_id, false);
                path_handles.push(handle);
            }
            
            // Build path with handles
            graph.build_path(seq.id.clone(), 
                path_handles.iter().map(|h| (h.node_id(), h.is_reverse())).collect());
        }
        
        // Build edges from paths
        let edges_to_add: Vec<(Handle, Handle)> = graph.paths.iter()
            .flat_map(|path| {
                (0..path.steps.len().saturating_sub(1))
                    .map(move |i| (path.steps[i], path.steps[i + 1]))
            })
            .collect();
            
        for (from, to) in edges_to_add {
            graph.add_edge(from, to);
        }
        
        if verbose {
            println!("Built bidirected graph: {} nodes, {} edges, {} paths", 
                     graph.nodes.len(), graph.edges.len(), graph.paths.len());
            
            // Debug: count unique union representatives
            let mut unique_unions = std::collections::HashSet::new();
            for seq in &self.sequences {
                for i in 0..seq.data.len() {
                    let global_pos = seq.offset + i;
                    let pos_fwd = make_pos(global_pos, false);
                    let union_fwd = self.union_find.find(pos_fwd);
                    unique_unions.insert(union_fwd);
                }
            }
            println!("Debug: {} unique union representatives from {} total positions", 
                     unique_unions.len(), self.sequences.iter().map(|s| s.data.len()).sum::<usize>());
        }
        
        Ok(graph)
    }
    
    /// Convert bidirected graph back to simple graph (temporary compatibility)
    pub fn bidirected_to_simple_graph(&self, bi_graph: BidirectedGraph) -> Graph {
        let mut graph = Graph::new();
        
        // Convert nodes
        for (id, bi_node) in bi_graph.nodes {
            graph.nodes.insert(id, Node {
                id,
                sequence: bi_node.sequence,
                rank: bi_node.rank.unwrap_or(id as u64) as f64,
            });
        }
        
        // Convert edges (lose orientation information)
        for edge in bi_graph.edges {
            graph.edges.insert(Edge {
                from: edge.from.node_id(),
                to: edge.to.node_id(),
            });
        }
        
        // Convert paths (lose orientation information)  
        for path in bi_graph.paths {
            let simple_path: Vec<usize> = path.steps.iter()
                .map(|handle| handle.node_id())
                .collect();
            graph.paths.push((path.name, simple_path));
        }
        
        graph
    }
    
    /// Find which sequence contains a given position
    fn find_sequence_for_position(&self, pos: Pos) -> Option<&Sequence> {
        let offset_val = offset(pos);
        self.sequences.iter().find(|seq| {
            offset_val >= seq.offset && offset_val < seq.offset + seq.data.len()
        })
    }
}