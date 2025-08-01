use handlegraph::handle::{Handle, NodeId, Direction, Edge};
use handlegraph::handlegraph::{HandleGraph, IntoHandles, IntoNeighbors};
use handlegraph::mutablehandlegraph::{MutableHandleGraph, AdditiveHandleGraph};
use handlegraph::pathhandlegraph::embedded_paths::{GraphPaths, GraphPathNames, MutableGraphPaths};
use handlegraph::hashgraph::HashGraph;
use std::collections::{HashMap, HashSet, VecDeque};

/// Represents a group of nodes that should be merged into one
#[derive(Debug, Clone)]
struct MergeGroup {
    nodes: Vec<NodeId>,        // Nodes in order
    new_sequence: Vec<u8>,     // Concatenated sequence
}

/// Find chains of nodes that can be merged
pub fn find_merge_groups(graph: &HashGraph) -> Vec<MergeGroup> {
    let mut groups = Vec::new();
    let mut visited = HashSet::new();
    
    // For each node, try to extend a chain in both directions
    for node_id in graph.handles() {
        if visited.contains(&node_id) {
            continue;
        }
        
        // Start a new chain
        let mut chain = VecDeque::new();
        chain.push_back(node_id);
        visited.insert(node_id);
        
        // Try to extend forward
        let mut current = node_id;
        loop {
            let handle = Handle::pack(current, false);
            let next_handles: Vec<_> = (&graph).neighbors(handle, Direction::Right).collect();
            
            if next_handles.len() == 1 {
                let next = next_handles[0];
                if !next.is_reverse() {
                    let next_id = next.id();
                    
                    // Check if next only has one predecessor
                    let prev_handles: Vec<_> = (&graph).neighbors(Handle::pack(next_id, false), Direction::Left).collect();
                    if prev_handles.len() == 1 && prev_handles[0].id() == current {
                        // Check if all paths that visit current also visit next in sequence
                        if can_merge_nodes(graph, current, next_id) {
                            chain.push_back(next_id);
                            visited.insert(next_id);
                            current = next_id;
                            continue;
                        }
                    }
                }
            }
            break;
        }
        
        // Try to extend backward
        current = node_id;
        loop {
            let handle = Handle::pack(current, false);
            let prev_handles: Vec<_> = (&graph).neighbors(handle, Direction::Left).collect();
            
            if prev_handles.len() == 1 {
                let prev = prev_handles[0];
                if !prev.is_reverse() {
                    let prev_id = prev.id();
                    
                    // Check if prev only has one successor
                    let next_handles: Vec<_> = (&graph).neighbors(Handle::pack(prev_id, false), Direction::Right).collect();
                    if next_handles.len() == 1 && next_handles[0].id() == current {
                        // Check if all paths that visit prev also visit current in sequence
                        if can_merge_nodes(graph, prev_id, current) {
                            chain.push_front(prev_id);
                            visited.insert(prev_id);
                            current = prev_id;
                            continue;
                        }
                    }
                }
            }
            break;
        }
        
        // If we found a chain of more than one node, create a merge group
        if chain.len() > 1 {
            // Build the concatenated sequence
            let mut new_sequence = Vec::new();
            for &node_id in &chain {
                let handle = Handle::pack(node_id, false);
                new_sequence.extend((&graph).sequence(handle));
            }
            
            groups.push(MergeGroup {
                nodes: chain.into_iter().collect(),
                new_sequence,
            });
        }
    }
    
    groups
}

/// Check if two adjacent nodes can be merged
fn can_merge_nodes(graph: &HashGraph, from_id: NodeId, to_id: NodeId) -> bool {
    // Get all paths on both nodes
    let from_handle = Handle::pack(from_id, false);
    let to_handle = Handle::pack(to_id, false);
    
    let from_paths: HashSet<_> = graph.occurrences_iter(from_handle)
        .map(|(path_id, _)| path_id)
        .collect();
    
    let to_paths: HashSet<_> = graph.occurrences_iter(to_handle)
        .map(|(path_id, _)| path_id)
        .collect();
    
    // They must have the same paths
    if from_paths != to_paths {
        return false;
    }
    
    // For each path, check that from is always followed by to
    for path_id in from_paths {
        let path_handles: Vec<_> = graph.get_path(&path_id).unwrap().collect();
        
        // Find all occurrences of from_handle and to_handle
        for i in 0..path_handles.len() {
            if path_handles[i].id() == from_id {
                // Check if the next handle is to_handle with matching orientation
                if i + 1 >= path_handles.len() || path_handles[i + 1].id() != to_id {
                    return false;
                }
                // Also check orientation consistency
                if path_handles[i].is_reverse() != path_handles[i + 1].is_reverse() {
                    return false;
                }
            }
        }
    }
    
    true
}

/// Create a projection mapping from old node IDs to new node IDs
fn create_node_projection(graph: &HashGraph, merge_groups: &[MergeGroup]) -> HashMap<NodeId, NodeId> {
    let mut projection = HashMap::new();
    let mut next_new_id = 1u64;
    
    // First, assign new IDs to merged groups
    for group in merge_groups {
        let new_id = NodeId::from(next_new_id);
        next_new_id += 1;
        
        // Map all nodes in the group to the new ID
        for &old_id in &group.nodes {
            projection.insert(old_id, new_id);
        }
    }
    
    // Then, assign new IDs to unmerged nodes
    for node_id in graph.handles() {
        if !projection.contains_key(&node_id) {
            let new_id = NodeId::from(next_new_id);
            next_new_id += 1;
            projection.insert(node_id, new_id);
        }
    }
    
    projection
}

/// Build a new compacted graph using the projection
pub fn build_compacted_graph(
    graph: &HashGraph,
    merge_groups: &[MergeGroup],
    projection: &HashMap<NodeId, NodeId>,
    verbose: bool
) -> Result<HashGraph, String> {
    let mut new_graph = HashGraph::new();
    
    // Create nodes in the new graph
    let mut created_nodes = HashSet::new();
    
    // First, create merged nodes
    for group in merge_groups {
        let new_id = projection[&group.nodes[0]];
        if created_nodes.insert(new_id) {
            new_graph.create_handle(&group.new_sequence, new_id);
        }
    }
    
    // Then, create unmerged nodes
    for node_id in graph.handles() {
        let new_id = projection[&node_id];
        if created_nodes.insert(new_id) {
            let handle = Handle::pack(node_id, false);
            let seq: Vec<u8> = graph.sequence(handle).collect();
            new_graph.create_handle(&seq, new_id);
        }
    }
    
    if verbose {
        eprintln!("[unchop] Created {} nodes in new graph", created_nodes.len());
    }
    
    // Create edges in the new graph
    let mut created_edges = HashSet::new();
    
    for edge in graph.edges() {
        let from = edge.0;
        let to = edge.1;
        
        // Map to new node IDs
        let new_from_id = projection[&from.id()];
        let new_to_id = projection[&to.id()];
        
        // Skip self-loops within merged groups
        if new_from_id == new_to_id && !from.is_reverse() && !to.is_reverse() {
            // This edge is internal to a merged group
            continue;
        }
        
        // Create new handles with same orientation
        let new_from = Handle::pack(new_from_id, from.is_reverse());
        let new_to = Handle::pack(new_to_id, to.is_reverse());
        let new_edge = Edge(new_from, new_to);
        
        // Avoid duplicate edges
        let edge_key = (new_from.as_integer(), new_to.as_integer());
        if created_edges.insert(edge_key) {
            new_graph.create_edge(new_edge);
        }
    }
    
    if verbose {
        eprintln!("[unchop] Created {} edges in new graph", created_edges.len());
    }
    
    // Rebuild paths using the projection
    for path_id in graph.path_ids() {
        let path_name = graph.path_handle_to_name(&path_id);
        let new_path_id = new_graph.create_path(&path_name, false)?;
        
        let mut last_added: Option<(NodeId, bool)> = None;
        
        for handle in graph.get_path(&path_id).unwrap() {
            let old_id = handle.id();
            let new_id = projection[&old_id];
            let is_reverse = handle.is_reverse();
            
            // Skip consecutive steps that map to the same new node
            if let Some((last_id, last_rev)) = last_added {
                if last_id == new_id && last_rev == is_reverse {
                    continue;
                }
            }
            
            let new_handle = Handle::pack(new_id, is_reverse);
            new_graph.path_append_step(&new_path_id, new_handle);
            last_added = Some((new_id, is_reverse));
        }
    }
    
    if verbose {
        eprintln!("[unchop] Rebuilt {} paths in new graph", new_graph.path_count());
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
        let total_nodes_to_remove: usize = merge_groups.iter()
            .map(|g| g.nodes.len() - 1)  // Each group reduces node count by len-1
            .sum();
        eprintln!("[unchop] Will remove {} nodes", total_nodes_to_remove);
    }
    
    if merge_groups.is_empty() {
        // Nothing to merge, return a copy of the original graph
        if verbose {
            eprintln!("[unchop] No nodes to merge, returning original graph");
        }
        return Ok(graph.clone());
    }
    
    // Create the projection mapping
    let projection = create_node_projection(graph, &merge_groups);
    
    // Build the new compacted graph
    let new_graph = build_compacted_graph(graph, &merge_groups, &projection, verbose)?;
    
    if verbose {
        eprintln!("[unchop] Created new graph with {} nodes", new_graph.node_count());
    }
    
    Ok(new_graph)
}