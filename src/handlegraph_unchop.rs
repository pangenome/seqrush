use handlegraph::handle::{Handle, NodeId, Direction, Edge};
use handlegraph::handlegraph::{HandleGraph, IntoHandles, IntoNeighbors};
use handlegraph::mutablehandlegraph::{MutableHandleGraph, AdditiveHandleGraph};
use handlegraph::pathhandlegraph::embedded_paths::{GraphPaths, GraphPathNames, MutableGraphPaths};
use handlegraph::pathhandlegraph::GraphPaths as PathTrait;
use handlegraph::hashgraph::HashGraph;
use std::collections::{HashMap, HashSet, VecDeque};

/// Find which nodes should be merged together
pub fn find_merge_groups(graph: &HashGraph) -> Vec<Vec<NodeId>> {
    let mut groups = Vec::new();
    let mut visited = HashSet::new();
    
    // For each node, try to extend a group in both directions
    for node_id in graph.handles() {
        if visited.contains(&node_id) {
            continue;
        }
        
        // Start a new group
        let mut group = VecDeque::new();
        group.push_back(node_id);
        visited.insert(node_id);
        
        // Extend forward
        let mut current = node_id;
        loop {
            let handle = Handle::pack(current, false);
            let next_handles: Vec<_> = graph.neighbors(handle, Direction::Right).collect();
            
            if next_handles.len() == 1 {
                let next = next_handles[0];
                if !next.is_reverse() {
                    let next_id = next.id();
                    
                    // Check if next only has one predecessor
                    let prev_handles: Vec<_> = graph.neighbors(Handle::pack(next_id, false), Direction::Left).collect();
                    if prev_handles.len() == 1 && prev_handles[0] == handle {
                        // Check paths are compatible
                        if paths_compatible(graph, handle, next) {
                            group.push_back(next_id);
                            visited.insert(next_id);
                            current = next_id;
                            continue;
                        }
                    }
                }
            }
            break;
        }
        
        // Extend backward
        current = node_id;
        loop {
            let handle = Handle::pack(current, false);
            let prev_handles: Vec<_> = graph.neighbors(handle, Direction::Left).collect();
            
            if prev_handles.len() == 1 {
                let prev = prev_handles[0];
                if !prev.is_reverse() {
                    let prev_id = prev.id();
                    
                    // Check if prev only has one successor
                    let next_handles: Vec<_> = graph.neighbors(Handle::pack(prev_id, false), Direction::Right).collect();
                    if next_handles.len() == 1 && next_handles[0] == handle {
                        // Check paths are compatible
                        if paths_compatible(graph, prev, handle) {
                            group.push_front(prev_id);
                            visited.insert(prev_id);
                            current = prev_id;
                            continue;
                        }
                    }
                }
            }
            break;
        }
        
        // Only keep groups with more than one node
        if group.len() > 1 {
            groups.push(group.into_iter().collect());
        }
    }
    
    groups
}

/// Check if all paths that go from one handle to another are compatible
fn paths_compatible(graph: &HashGraph, from: Handle, to: Handle) -> bool {
    // Get all paths on the from handle
    let from_paths: HashSet<_> = graph.steps_on_handle(from)
        .map(|(path_id, _)| path_id)
        .collect();
    
    // Get all paths on the to handle  
    let to_paths: HashSet<_> = graph.steps_on_handle(to)
        .map(|(path_id, _)| path_id)
        .collect();
    
    // They must have the same paths
    if from_paths != to_paths {
        return false;
    }
    
    // For each path, check that from is followed by to
    for path_id in from_paths {
        let path_handles: Vec<_> = graph.path(&path_id).collect();
        
        // Find occurrences of from and to
        let mut from_positions = Vec::new();
        let mut to_positions = Vec::new();
        
        for (i, &h) in path_handles.iter().enumerate() {
            if h == from {
                from_positions.push(i);
            }
            if h == to {
                to_positions.push(i);
            }
        }
        
        // Check that each from is followed by to
        for &from_pos in &from_positions {
            if from_pos + 1 >= path_handles.len() || path_handles[from_pos + 1] != to {
                return false;
            }
        }
    }
    
    true
}

/// Create a new compacted graph by merging the specified groups
pub fn create_compacted_graph(
    graph: &HashGraph,
    merge_groups: &[Vec<NodeId>]
) -> Result<HashGraph, String> {
    let mut new_graph = HashGraph::new();
    
    // Create mapping from old node IDs to new node IDs
    let mut node_mapping: HashMap<NodeId, NodeId> = HashMap::new();
    let mut next_node_id = 1u64;
    
    // First, handle merged groups
    for group in merge_groups {
        // Concatenate sequences
        let mut merged_seq = Vec::new();
        for &node_id in group {
            let handle = Handle::pack(node_id, false);
            merged_seq.extend_from_slice(graph.sequence(handle));
        }
        
        // Create new node
        let new_node_id = NodeId::from(next_node_id);
        next_node_id += 1;
        new_graph.create_handle(&merged_seq, new_node_id);
        
        // Map all nodes in group to the new node
        for &node_id in group {
            node_mapping.insert(node_id, new_node_id);
        }
    }
    
    // Then, handle unmerged nodes
    for node_id in graph.handles() {
        if !node_mapping.contains_key(&node_id) {
            let handle = Handle::pack(node_id, false);
            let seq = graph.sequence(handle);
            
            let new_node_id = NodeId::from(next_node_id);
            next_node_id += 1;
            new_graph.create_handle(seq, new_node_id);
            
            node_mapping.insert(node_id, new_node_id);
        }
    }
    
    // Create edges in the new graph
    let mut edges_added = HashSet::new();
    
    for edge in graph.edges() {
        let from = edge.0;
        let to = edge.1;
        
        // Map to new node IDs
        let new_from_id = node_mapping[&from.id()];
        let new_to_id = node_mapping[&to.id()];
        
        // Skip if this edge is internal to a merged group
        if new_from_id == new_to_id && !from.is_reverse() && !to.is_reverse() {
            continue;
        }
        
        let new_from = Handle::pack(new_from_id, from.is_reverse());
        let new_to = Handle::pack(new_to_id, to.is_reverse());
        let new_edge = Edge(new_from, new_to);
        
        // Avoid duplicate edges
        let edge_key = (new_from.as_integer(), new_to.as_integer());
        if !edges_added.contains(&edge_key) {
            new_graph.create_edge(&new_edge);
            edges_added.insert(edge_key);
        }
    }
    
    // Create paths in the new graph
    for path_id in graph.path_ids() {
        let path_name = graph.path_handle_to_name(&path_id);
        let new_path_id = new_graph.create_path(&path_name, false)?;
        
        let mut last_added: Option<Handle> = None;
        
        for handle in graph.path(&path_id) {
            let new_node_id = node_mapping[&handle.id()];
            let new_handle = Handle::pack(new_node_id, handle.is_reverse());
            
            // Skip consecutive identical handles (from merging)
            if last_added != Some(new_handle) {
                new_graph.path_append_step(&new_path_id, new_handle);
                last_added = Some(new_handle);
            }
        }
    }
    
    Ok(new_graph)
}

/// Unchop the graph by finding and merging simple chains
pub fn unchop(graph: &HashGraph, verbose: bool) -> Result<HashGraph, String> {
    if verbose {
        eprintln!("[unchop] Starting with {} nodes", graph.node_count());
    }
    
    // Find groups of nodes to merge
    let merge_groups = find_merge_groups(graph);
    
    if verbose {
        eprintln!("[unchop] Found {} groups to merge", merge_groups.len());
        for group in &merge_groups {
            eprintln!("  Group of {} nodes: {:?}", group.len(), group);
        }
    }
    
    // Create new compacted graph
    let new_graph = create_compacted_graph(graph, &merge_groups)?;
    
    if verbose {
        eprintln!("[unchop] Created new graph with {} nodes", new_graph.node_count());
    }
    
    Ok(new_graph)
}