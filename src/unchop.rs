use handlegraph::handle::{Handle, NodeId, Direction, Edge};
use handlegraph::handlegraph::HandleGraph;
use handlegraph::mutablehandlegraph::MutableHandleGraph;
use handlegraph::pathhandlegraph::embedded_paths::{GraphPaths, MutableGraphPaths, IntoNodeOccurrences};
use handlegraph::hashgraph::HashGraph;
use std::collections::{HashMap, HashSet};

/// Check if two nodes are perfect path neighbors
/// (all paths that visit the left node continue directly to the right node)
fn nodes_are_perfect_path_neighbors(graph: &HashGraph, left: Handle, right: Handle) -> bool {
    let mut expected_next = 0;
    let mut ok = true;
    
    // Check all path steps on the left handle
    for (path_id, step_index) in graph.steps_iter_on_handle(left) {
        // Get the actual handle at this step (might be reverse of what we're checking)
        let steps: Vec<_> = graph.steps_iter(&path_id).collect();
        if step_index >= steps.len() {
            continue;
        }
        
        let step_handle = steps[step_index];
        let step_is_reverse = step_handle != left;
        
        // Check if there's a next step
        let next_idx = if step_is_reverse {
            if step_index == 0 { 
                ok = false;
                break;
            }
            step_index - 1
        } else {
            step_index + 1
        };
        
        if next_idx >= steps.len() {
            ok = false;
            break;
        }
        
        let next_handle = steps[next_idx];
        let expected_handle = if step_is_reverse {
            right.flip()
        } else {
            right
        };
        
        if next_handle != expected_handle {
            ok = false;
            break;
        }
        
        expected_next += 1;
    }
    
    if !ok {
        return false;
    }
    
    // Count steps on right handle
    let observed_next = graph.steps_on_handle(right).count();
    
    observed_next == expected_next
}

/// Find simple components that can be merged
fn find_simple_components(graph: &HashGraph, min_size: usize) -> Vec<Vec<Handle>> {
    let mut dset = union_find::UnionFind::new(graph.node_count() * 2);
    let mut handle_to_idx: HashMap<Handle, usize> = HashMap::new();
    let mut idx_to_handle: Vec<Handle> = Vec::new();
    
    // Map handles to indices
    let mut idx = 0;
    for node_id in graph.handles_iter() {
        for &is_reverse in &[false, true] {
            let handle = Handle::new(node_id, is_reverse);
            handle_to_idx.insert(handle, idx);
            idx_to_handle.push(handle);
            idx += 1;
        }
    }
    
    // Find mergeable neighbors
    for node_id in graph.handles_iter() {
        for &is_reverse in &[false, true] {
            let handle = Handle::new(node_id, is_reverse);
            
            // Check forward direction
            let forward_neighbors: Vec<_> = graph.neighbors(handle, Direction::Right).collect();
            if forward_neighbors.len() == 1 {
                let next = forward_neighbors[0];
                if handle.id() != next.id() {
                    let backward_neighbors: Vec<_> = graph.neighbors(next, Direction::Left).collect();
                    if backward_neighbors.len() == 1 && nodes_are_perfect_path_neighbors(graph, handle, next) {
                        let from_idx = handle_to_idx[&handle];
                        let to_idx = handle_to_idx[&next];
                        dset.unite(from_idx, to_idx);
                    }
                }
            }
        }
    }
    
    // Collect components
    let mut components: HashMap<usize, Vec<Handle>> = HashMap::new();
    for node_id in graph.handles_iter() {
        let handle = Handle::new(node_id, false);
        let idx = handle_to_idx[&handle];
        let root = dset.find(idx);
        components.entry(root).or_default().push(handle);
    }
    
    // Filter and order components
    let mut result = Vec::new();
    for (_, mut handles) in components {
        if handles.len() >= min_size {
            // Order the handles in the component
            let ordered = order_component(graph, &handles);
            if !ordered.is_empty() {
                result.push(ordered);
            }
        }
    }
    
    result
}

/// Order handles in a component from start to end
fn order_component(graph: &HashGraph, handles: &[Handle]) -> Vec<Handle> {
    if handles.is_empty() {
        return Vec::new();
    }
    
    let handle_set: HashSet<_> = handles.iter().map(|h| h.id()).collect();
    
    // Find a start handle (one with no predecessor in the component)
    let mut start = None;
    for &handle in handles {
        let mut has_pred_in_comp = false;
        for pred in graph.neighbors(handle, Direction::Left) {
            if handle_set.contains(&pred.id()) {
                has_pred_in_comp = true;
                break;
            }
        }
        if !has_pred_in_comp {
            start = Some(handle);
            break;
        }
    }
    
    let mut current = match start {
        Some(h) => h,
        None => return Vec::new(), // Circular component, skip
    };
    
    let mut ordered = vec![current];
    let mut visited = HashSet::new();
    visited.insert(current.id());
    
    // Follow the chain
    loop {
        let mut found_next = None;
        for next in graph.neighbors(current, Direction::Right) {
            if handle_set.contains(&next.id()) && !visited.contains(&next.id()) {
                found_next = Some(next);
                break;
            }
        }
        
        match found_next {
            Some(next) => {
                ordered.push(next);
                visited.insert(next.id());
                current = next;
            }
            None => break,
        }
    }
    
    ordered
}

/// Concatenate nodes into a single new node
fn concat_nodes(graph: &mut HashGraph, nodes: &[Handle]) -> Handle {
    // Build the concatenated sequence
    let mut seq = Vec::new();
    for &handle in nodes {
        let node_seq = graph.sequence(handle);
        seq.extend_from_slice(&node_seq);
    }
    
    // Create new node
    let new_id = {
        let max_id = graph.handles_iter().map(|h| h.as_integer()).max().unwrap_or(0);
        NodeId::from((max_id + 1) as u64)
    };
    let new_handle = graph.create_handle(&seq, new_id);
    
    // Get edges from first and last nodes
    let mut left_neighbors = Vec::new();
    for neighbor in graph.neighbors(nodes[0], Direction::Left) {
        if neighbor == nodes.last().copied().unwrap() {
            left_neighbors.push(new_handle);
        } else if neighbor == nodes[0].flip() {
            left_neighbors.push(new_handle.flip());
        } else {
            left_neighbors.push(neighbor);
        }
    }
    
    let mut right_neighbors = Vec::new();
    for neighbor in graph.neighbors(nodes.last().copied().unwrap(), Direction::Right) {
        if neighbor == nodes[0] {
            // Already handled from other direction
        } else if neighbor == nodes.last().copied().unwrap().flip() {
            right_neighbors.push(new_handle.flip());
        } else {
            right_neighbors.push(neighbor);
        }
    }
    
    // Create edges
    for &neighbor in &left_neighbors {
        graph.create_edge(&Edge(neighbor, new_handle));
    }
    for &neighbor in &right_neighbors {
        graph.create_edge(&Edge(new_handle, neighbor));
    }
    
    // Update paths
    let paths_to_update: Vec<_> = graph.steps_iter_on_handle(nodes[0])
        .map(|(p, _)| p)
        .collect();
    
    for path_id in paths_to_update {
        // Find the range of steps to rewrite
        let steps: Vec<_> = graph.steps_iter(&path_id).collect();
        let mut start_idx = None;
        let mut end_idx = None;
        let mut is_reverse = false;
        
        for (i, &step) in steps.iter().enumerate() {
            if step.id() == nodes[0].id() {
                start_idx = Some(i);
                is_reverse = step != nodes[0];
                
                // Find end
                let target = if is_reverse { 
                    nodes.last().copied().unwrap().flip() 
                } else { 
                    nodes.last().copied().unwrap() 
                };
                
                for j in i..steps.len() {
                    if steps[j] == target {
                        end_idx = Some(j);
                        break;
                    }
                }
                break;
            }
        }
        
        if let (Some(start), Some(end)) = (start_idx, end_idx) {
            // Rebuild the path
            let mut new_steps = Vec::new();
            
            // Before the merged section
            new_steps.extend_from_slice(&steps[..start]);
            
            // The merged node
            new_steps.push(if is_reverse { new_handle.flip() } else { new_handle });
            
            // After the merged section
            if end + 1 < steps.len() {
                new_steps.extend_from_slice(&steps[end + 1..]);
            }
            
            // Clear and rebuild the path
            graph.path_clear(&path_id);
            for &handle in &new_steps {
                graph.append_step(&path_id, handle);
            }
        }
    }
    
    // Delete old nodes
    for &handle in nodes {
        graph.destroy_handle(handle);
    }
    
    new_handle
}

/// Unchop the graph by merging simple paths
pub fn unchop(graph: &mut HashGraph, verbose: bool) -> Result<usize, Box<dyn std::error::Error>> {
    let initial_nodes = graph.node_count();
    
    if verbose {
        eprintln!("[unchop] Starting with {} nodes", initial_nodes);
    }
    
    // Find components to merge
    let components = find_simple_components(graph, 2);
    
    if verbose {
        eprintln!("[unchop] Found {} components to merge", components.len());
    }
    
    let mut nodes_removed = 0;
    for component in components {
        if component.len() >= 2 {
            concat_nodes(graph, &component);
            nodes_removed += component.len() - 1;
        }
    }
    
    let final_nodes = graph.node_count();
    
    if verbose {
        eprintln!("[unchop] Merged {} nodes, final count: {}", nodes_removed, final_nodes);
    }
    
    Ok(nodes_removed)
}