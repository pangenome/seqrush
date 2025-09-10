/// ODGI topological sort adapted for HashGraph
use handlegraph::handle::{Edge, Handle, NodeId};
use handlegraph::handlegraph::*;
use handlegraph::hashgraph::HashGraph;
use handlegraph::mutablehandlegraph::*;
use handlegraph::pathhandlegraph::*;
use std::collections::{HashMap, HashSet};

/// Find all nodes with no edges on their left sides (heads)
pub fn head_nodes(graph: &HashGraph) -> Vec<Handle> {
    let mut heads = Vec::new();
    
    for handle in graph.handles() {
        if !handle.is_reverse() {  // Only check forward orientation
            let mut has_left_edges = false;
            
            // Check all edges to see if any come to this handle's left side
            // Left side of forward handle = edges coming TO the forward orientation
            for edge in graph.edges() {
                if edge.1 == handle {  // edge.1 is the target
                    has_left_edges = true;
                    break;
                }
            }
            
            if !has_left_edges {
                heads.push(handle);
            }
        }
    }
    
    // Sort for deterministic behavior
    heads.sort_by_key(|h| h.id().0);
    heads
}

/// Find all nodes with no edges on their right sides (tails)
pub fn tail_nodes(graph: &HashGraph) -> Vec<Handle> {
    let mut tails = Vec::new();
    
    for handle in graph.handles() {
        if !handle.is_reverse() {  // Only check forward orientation
            let mut has_right_edges = false;
            
            // Check all edges to see if any go from this handle's right side
            // Right side of forward handle = edges going FROM the forward orientation
            for edge in graph.edges() {
                if edge.0 == handle {  // edge.0 is the source
                    has_right_edges = true;
                    break;
                }
            }
            
            if !has_right_edges {
                tails.push(handle);
            }
        }
    }
    
    // Sort for deterministic behavior
    tails.sort_by_key(|h| h.id().0);
    tails
}

/// ODGI's topological_order algorithm adapted for HashGraph
pub fn odgi_topological_order(
    graph: &HashGraph,
    use_heads: bool,
    use_tails: bool,
    verbose: bool,
) -> Vec<Handle> {
    let mut sorted = Vec::new();
    
    if graph.node_count() == 0 {
        return sorted;
    }
    
    // Track which node IDs have been added to the result
    let mut emitted_nodes = HashSet::new();
    
    // S - set of oriented handles ready to be processed
    let mut s = HashSet::new();
    
    // Unvisited - track which handles haven't been processed yet
    let mut unvisited = HashMap::new();
    for handle in graph.handles() {
        // Add both orientations
        unvisited.insert(handle, true);
        unvisited.insert(handle.flip(), true);
    }
    
    // Seeds for breaking cycles - maps handle to priority
    let mut seeds: HashMap<Handle, i32> = HashMap::new();
    
    // Track masked (logically removed) edges
    let mut masked_edges = HashSet::new();
    
    // Initialize with heads or tails if requested
    if use_heads {
        for handle in head_nodes(graph) {
            s.insert(handle);
            unvisited.remove(&handle);
            unvisited.remove(&handle.flip());
        }
        if verbose {
            eprintln!("[odgi_topo] Starting with {} head nodes", s.len());
        }
    } else if use_tails {
        for handle in tail_nodes(graph) {
            s.insert(handle);
            unvisited.remove(&handle);
            unvisited.remove(&handle.flip());
        }
        if verbose {
            eprintln!("[odgi_topo] Starting with {} tail nodes", s.len());
        }
    }
    
    // Collect all edges once for efficiency
    let all_edges: Vec<Edge> = graph.edges().collect();
    
    // Main loop - continue until all nodes are visited
    while !unvisited.is_empty() || !s.is_empty() {
        
        // If S is empty, need to pick a seed to break into a cycle
        if s.is_empty() {
            // First try previously identified seeds
            let mut found_seed = false;
            
            // Sort seeds for deterministic selection
            let mut seed_handles: Vec<_> = seeds.keys().cloned().collect();
            seed_handles.sort_by_key(|h| (h.id().0, h.is_reverse()));
            
            for handle in seed_handles {
                if unvisited.get(&handle).copied().unwrap_or(false) {
                    s.insert(handle);
                    unvisited.remove(&handle);
                    unvisited.remove(&handle.flip());
                    found_seed = true;
                    if verbose {
                        eprintln!("[odgi_topo] Using seed: node {} orient {}", 
                                 handle.id().0, 
                                 if handle.is_reverse() { "rev" } else { "fwd" });
                    }
                    break;
                }
            }
            
            // If no seeds available, pick arbitrary unvisited handle
            if !found_seed {
                // Get minimum unvisited handle for deterministic behavior
                let mut handles: Vec<_> = unvisited.keys()
                    .filter(|h| !h.is_reverse()) // Prefer forward orientation
                    .cloned()
                    .collect();
                handles.sort_by_key(|h| h.id().0);
                
                if let Some(&handle) = handles.first() {
                    s.insert(handle);
                    unvisited.remove(&handle);
                    unvisited.remove(&handle.flip());
                    if verbose {
                        eprintln!("[odgi_topo] Using arbitrary: node {} orient fwd", 
                                 handle.id().0);
                    }
                }
            }
        }
        
        // Process handles in S
        while !s.is_empty() {
            // Get minimum handle for deterministic behavior
            let handle = *s.iter().min_by_key(|h| (h.id().0, h.is_reverse())).unwrap();
            s.remove(&handle);
            
            // Only emit each node once (in forward orientation)
            if !handle.is_reverse() && emitted_nodes.insert(handle.id()) {
                sorted.push(handle);
            }
            
            // Look at edges coming into the left side of this handle (going backward)
            // These are edges we're "consuming" by placing this handle
            for edge in &all_edges {
                let edge_key = (edge.0.as_integer(), edge.1.as_integer());
                if edge.1 == handle && !masked_edges.contains(&edge_key) {
                    masked_edges.insert(edge_key);
                }
            }
            
            // Look at edges going out from the right side of this handle (going forward)
            // These lead to potential next handles to process
            for edge in &all_edges {
                let edge_key = (edge.0.as_integer(), edge.1.as_integer());
                if edge.0 == handle && !masked_edges.contains(&edge_key) {
                    masked_edges.insert(edge_key);
                    
                    let next_handle = edge.1;
                    
                    // Only process if not yet visited
                    if unvisited.get(&next_handle).copied().unwrap_or(false) {
                        // Check if next_handle has any other unmasked incoming edges
                        let mut has_unmasked_incoming = false;
                        
                        for other_edge in &all_edges {
                            let other_key = (other_edge.0.as_integer(), other_edge.1.as_integer());
                            if other_edge.1 == next_handle && 
                               !masked_edges.contains(&other_key) {
                                has_unmasked_incoming = true;
                                break;
                            }
                        }
                        
                        if !has_unmasked_incoming {
                            // No more incoming edges, ready to process
                            s.insert(next_handle);
                            unvisited.remove(&next_handle);
                            unvisited.remove(&next_handle.flip());
                        } else {
                            // Still has dependencies, mark as potential seed for cycle breaking
                            seeds.insert(next_handle, 1);
                        }
                    }
                }
            }
        }
    }
    
    if verbose {
        eprintln!("[odgi_topo] Topological sort complete: {} nodes", sorted.len());
    }
    
    sorted
}

/// Apply ODGI topological sort to a HashGraph and create a new sorted graph
pub fn apply_odgi_sort(
    graph: &HashGraph,
    verbose: bool,
) -> Result<HashGraph, Box<dyn std::error::Error>> {
    if verbose {
        eprintln!("[odgi_sort] Starting ODGI topological sort...");
    }
    
    let ordering = odgi_topological_order(graph, true, false, verbose);
    
    if ordering.is_empty() {
        // Return empty graph if no nodes
        let empty = HashGraph::new();
        return Ok(empty);
    }
    
    // Create old to new ID mapping (1-based numbering)
    let mut old_to_new = HashMap::new();
    for (new_idx, handle) in ordering.iter().enumerate() {
        old_to_new.insert(handle.id(), NodeId::from((new_idx + 1) as u64));
    }
    
    if verbose {
        eprintln!("[odgi_sort] Renumbering {} nodes", old_to_new.len());
    }
    
    // Build new graph with renumbered nodes
    let mut new_graph = HashGraph::new();
    
    // Create nodes with new IDs
    for (idx, &handle) in ordering.iter().enumerate() {
        let new_id = NodeId::from((idx + 1) as u64);
        let seq: Vec<u8> = graph.sequence(handle).collect();
        new_graph.create_handle(&seq, new_id);
    }
    
    // Add edges with remapped IDs
    let mut seen_edges = HashSet::new();
    for edge in graph.edges() {
        let old_from_id = edge.0.id();
        let old_to_id = edge.1.id();
        
        // Both nodes must exist in the mapping
        if let (Some(&new_from_id), Some(&new_to_id)) = 
            (old_to_new.get(&old_from_id), old_to_new.get(&old_to_id)) {
            
            let new_edge = Edge(
                Handle::pack(new_from_id, edge.0.is_reverse()),
                Handle::pack(new_to_id, edge.1.is_reverse())
            );
            
            // Avoid duplicate edges
            let edge_key = (new_edge.0.as_integer(), new_edge.1.as_integer());
            if seen_edges.insert(edge_key) {
                new_graph.create_edge(new_edge);
            }
        }
    }
    
    // Copy paths with remapped node IDs
    for path_id in graph.path_ids() {
        let path_name: Vec<u8> = graph.get_path_name(path_id).unwrap().collect();
        let new_path_id = new_graph.create_path(&path_name, false).unwrap();
        
        if let Some(path_ref) = graph.get_path_ref(path_id) {
            for step in path_ref.steps() {
                let handle = step.1;
                let old_id = handle.id();
                
                if let Some(&new_id) = old_to_new.get(&old_id) {
                    let new_handle = Handle::pack(new_id, handle.is_reverse());
                    new_graph.path_append_step(new_path_id, new_handle);
                }
            }
        }
    }
    
    if verbose {
        eprintln!(
            "[odgi_sort] Sort complete: {} nodes, {} edges, {} paths",
            new_graph.node_count(),
            new_graph.edge_count(),
            new_graph.path_count()
        );
    }
    
    Ok(new_graph)
}