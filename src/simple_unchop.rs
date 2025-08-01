use handlegraph::handle::{Handle, NodeId, Edge};
use handlegraph::handlegraph::{HandleGraph, IntoHandles, IntoEdges, IntoSequences};
use handlegraph::mutablehandlegraph::AdditiveHandleGraph;
use handlegraph::pathhandlegraph::embedded_paths::{MutableGraphPaths, GraphPathNames, IntoPathIds, GraphPathsRef};
use handlegraph::pathhandlegraph::path::PathSteps;
use handlegraph::hashgraph::HashGraph;
use std::collections::{HashMap, HashSet};

/// Simple unchop that finds linear chains and merges them
pub fn simple_unchop(graph: &HashGraph, verbose: bool) -> Result<HashGraph, String> {
    if verbose {
        eprintln!("[unchop] Starting with {} nodes", graph.node_count());
    }
    
    // Find simple chains: nodes with single in/out edges
    let mut chains: Vec<Vec<NodeId>> = Vec::new();
    let mut in_chain = HashSet::new();
    
    // Get node degrees
    let mut in_degree: HashMap<NodeId, usize> = HashMap::new();
    let mut out_degree: HashMap<NodeId, usize> = HashMap::new();
    let mut predecessors: HashMap<NodeId, NodeId> = HashMap::new();
    let mut successors: HashMap<NodeId, NodeId> = HashMap::new();
    
    // Count edges and track neighbors
    for edge in graph.edges() {
        let from_id = edge.0.id();
        let to_id = edge.1.id();
        
        // Skip self-loops and reverse edges for now
        if from_id == to_id || edge.0.is_reverse() || edge.1.is_reverse() {
            continue;
        }
        
        *out_degree.entry(from_id).or_insert(0) += 1;
        *in_degree.entry(to_id).or_insert(0) += 1;
        
        // Only track if single predecessor/successor
        if out_degree[&from_id] == 1 {
            successors.insert(from_id, to_id);
        } else {
            successors.remove(&from_id);
        }
        
        if in_degree[&to_id] == 1 {
            predecessors.insert(to_id, from_id);
        } else {
            predecessors.remove(&to_id);
        }
    }
    
    // Find chains by starting from nodes with in_degree=0 or out_degree=0
    for handle in graph.handles() {
        let node_id = handle.id();
        if in_chain.contains(&node_id) {
            continue;
        }
        
        let in_deg = in_degree.get(&node_id).copied().unwrap_or(0);
        let out_deg = out_degree.get(&node_id).copied().unwrap_or(0);
        
        // Start a chain from nodes that are at the beginning of a potential chain
        if (in_deg == 0 || in_deg == 1) && out_deg == 1 {
            let mut chain = vec![node_id];
            in_chain.insert(node_id);
            
            // Extend the chain forward
            let mut current = node_id;
            while let Some(&next) = successors.get(&current) {
                if in_chain.contains(&next) {
                    break;
                }
                
                let next_in = in_degree.get(&next).copied().unwrap_or(0);
                let next_out = out_degree.get(&next).copied().unwrap_or(0);
                
                if next_in == 1 && next_out <= 1 {
                    chain.push(next);
                    in_chain.insert(next);
                    current = next;
                } else {
                    break;
                }
            }
            
            if chain.len() > 1 {
                chains.push(chain);
            }
        }
    }
    
    if verbose {
        eprintln!("[unchop] Found {} chains to merge", chains.len());
        let total_reduction: usize = chains.iter().map(|c| c.len() - 1).sum();
        eprintln!("[unchop] Will reduce node count by {}", total_reduction);
    }
    
    if chains.is_empty() {
        // Nothing to merge, rebuild graph as-is
        let mut new_graph = HashGraph::new();
        
        // Copy nodes
        for handle in graph.handles() {
            let node_id = handle.id();
            let seq: Vec<u8> = graph.sequence(handle).collect();
            new_graph.create_handle(&seq, node_id);
        }
        
        // Copy edges
        for edge in graph.edges() {
            new_graph.create_edge(edge);
        }
        
        // Copy paths
        for path_id in graph.path_ids() {
            let path_name: Vec<u8> = graph.get_path_name(path_id).unwrap().collect();
            let new_path_id = new_graph.create_path(&path_name, false)
                .ok_or_else(|| "Failed to create path".to_string())?;
            if let Some(path_ref) = graph.get_path_ref(path_id) {
                for step in path_ref.steps() {
                    new_graph.path_append_step(new_path_id, step.1);
                }
            }
        }
        
        return Ok(new_graph);
    }
    
    // Create mapping from old nodes to new nodes
    let mut node_mapping: HashMap<NodeId, NodeId> = HashMap::new();
    let mut next_new_id = 1u64;
    
    // Assign new IDs to merged chains
    for chain in &chains {
        let new_id = NodeId::from(next_new_id);
        next_new_id += 1;
        
        for &old_id in chain {
            node_mapping.insert(old_id, new_id);
        }
    }
    
    // Assign new IDs to unmerged nodes
    for handle in graph.handles() {
        let node_id = handle.id();
        if !node_mapping.contains_key(&node_id) {
            let new_id = NodeId::from(next_new_id);
            next_new_id += 1;
            node_mapping.insert(node_id, new_id);
        }
    }
    
    // Build new graph
    let mut new_graph = HashGraph::new();
    let mut created_nodes = HashSet::new();
    
    // Create merged nodes
    for chain in &chains {
        let new_id = node_mapping[&chain[0]];
        if created_nodes.insert(new_id) {
            // Concatenate sequences
            let mut seq = Vec::new();
            for &node_id in chain {
                let handle = Handle::pack(node_id, false);
                seq.extend(graph.sequence(handle));
            }
            new_graph.create_handle(&seq, new_id);
        }
    }
    
    // Create unmerged nodes
    for handle in graph.handles() {
        let node_id = handle.id();
        let new_id = node_mapping[&node_id];
        if created_nodes.insert(new_id) {
            let handle = Handle::pack(node_id, false);
            let seq: Vec<u8> = graph.sequence(handle).collect();
            new_graph.create_handle(&seq, new_id);
        }
    }
    
    // Create edges
    let mut created_edges = HashSet::new();
    for edge in graph.edges() {
        let old_from = edge.0.id();
        let old_to = edge.1.id();
        let new_from = node_mapping[&old_from];
        let new_to = node_mapping[&old_to];
        
        // Skip internal edges within merged chains
        if new_from == new_to && !edge.0.is_reverse() && !edge.1.is_reverse() {
            continue;
        }
        
        let new_edge = Edge(
            Handle::pack(new_from, edge.0.is_reverse()),
            Handle::pack(new_to, edge.1.is_reverse())
        );
        
        let key = (new_edge.0.as_integer(), new_edge.1.as_integer());
        if created_edges.insert(key) {
            new_graph.create_edge(new_edge);
        }
    }
    
    // Recreate paths
    for path_id in graph.path_ids() {
        let path_name: Vec<u8> = graph.get_path_name(path_id).unwrap().collect();
        let new_path_id = new_graph.create_path(&path_name, false)
            .ok_or_else(|| "Failed to create path".to_string())?;
        
        let mut last_added: Option<(NodeId, bool)> = None;
        
        if let Some(path_ref) = graph.get_path_ref(path_id) {
            for step in path_ref.steps() {
                let handle = step.1;  // Step is a tuple (StepIx, Handle)
                let old_id = handle.id();
                let new_id = node_mapping[&old_id];
                let is_rev = handle.is_reverse();
                
                // Skip consecutive steps that map to same node
                if let Some((last_id, last_rev)) = last_added {
                    if last_id == new_id && last_rev == is_rev {
                        continue;
                    }
                }
                
                new_graph.path_append_step(new_path_id, Handle::pack(new_id, is_rev));
                last_added = Some((new_id, is_rev));
            }
        }
    }
    
    if verbose {
        eprintln!("[unchop] Created new graph with {} nodes", new_graph.node_count());
    }
    
    Ok(new_graph)
}