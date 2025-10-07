use handlegraph::handle::{Edge, Handle, NodeId};
use handlegraph::handlegraph::{HandleGraph, IntoEdges, IntoHandles, IntoSequences};
use handlegraph::hashgraph::HashGraph;
use handlegraph::mutablehandlegraph::AdditiveHandleGraph;
use handlegraph::pathhandlegraph::embedded_paths::{
    GraphPathNames, GraphPathsRef, IntoPathIds, MutableGraphPaths,
};
use handlegraph::pathhandlegraph::path::PathSteps;
use std::collections::{HashMap, HashSet};

/// Bidirected unchop that properly handles all orientations
pub fn bidirected_unchop(graph: &HashGraph, verbose: bool) -> Result<HashGraph, String> {
    if verbose {
        eprintln!(
            "[unchop] Starting bidirected unchop with {} nodes",
            graph.node_count()
        );
    }

    // A simple chain is a sequence of nodes where each internal node has exactly one
    // predecessor and one successor, considering all orientations

    // First, compute degree information for all handles (both orientations)
    let mut out_degree: HashMap<Handle, usize> = HashMap::new();
    let mut in_degree: HashMap<Handle, usize> = HashMap::new();
    let mut successors: HashMap<Handle, Handle> = HashMap::new();
    let mut predecessors: HashMap<Handle, Handle> = HashMap::new();

    // Count edges
    for edge in graph.edges() {
        let Edge(from, to) = edge;

        // Outgoing edge from 'from'
        *out_degree.entry(from).or_insert(0) += 1;
        // Incoming edge to 'to'
        *in_degree.entry(to).or_insert(0) += 1;

        // Track single successors/predecessors
        if out_degree[&from] == 1 {
            successors.insert(from, to);
        } else {
            successors.remove(&from);
        }

        if in_degree[&to] == 1 {
            predecessors.insert(to, from);
        } else {
            predecessors.remove(&to);
        }
    }

    if verbose {
        let total_handles = out_degree.len() + in_degree.len();
        let single_out: usize = out_degree.values().filter(|&&d| d == 1).count();
        let single_in: usize = in_degree.values().filter(|&&d| d == 1).count();
        eprintln!(
            "[unchop] Handles with single successor: {}/{}",
            single_out, total_handles
        );
        eprintln!(
            "[unchop] Handles with single predecessor: {}/{}",
            single_in, total_handles
        );
    }

    // Find simple chains
    let mut chains: Vec<Vec<Handle>> = Vec::new();
    let mut used_handles = HashSet::new();

    // Start chains from handles that are at the beginning (in_degree 0 or not simple)
    for handle in graph.handles() {
        // Check both orientations
        for &start_handle in &[handle, handle.flip()] {
            if used_handles.contains(&start_handle) {
                continue;
            }

            let in_deg = in_degree.get(&start_handle).copied().unwrap_or(0);
            let out_deg = out_degree.get(&start_handle).copied().unwrap_or(0);

            // Start a chain if this looks like a chain beginning
            if out_deg == 1 && (in_deg == 0 || in_deg > 1) {
                let mut chain = vec![start_handle];
                used_handles.insert(start_handle);

                // Extend the chain
                let mut current = start_handle;
                while let Some(&next) = successors.get(&current) {
                    if used_handles.contains(&next) {
                        break;
                    }

                    let next_in = in_degree.get(&next).copied().unwrap_or(0);
                    let next_out = out_degree.get(&next).copied().unwrap_or(0);

                    // Continue chain if next node is simple (1 in, 1 out) or terminal (1 in, 0 out)
                    if next_in == 1 && next_out <= 1 {
                        // Check path consistency before adding
                        if !check_path_consistency(graph, current, next) {
                            break;
                        }

                        chain.push(next);
                        used_handles.insert(next);
                        current = next;
                    } else {
                        break;
                    }
                }

                if chain.len() > 1 {
                    chains.push(chain);
                } else {
                    // Return unused handle
                    used_handles.remove(&start_handle);
                }
            }
        }
    }

    // Also look for cycles (all nodes have degree 1-1)
    for handle in graph.handles() {
        for &start_handle in &[handle, handle.flip()] {
            if used_handles.contains(&start_handle) {
                continue;
            }

            let in_deg = in_degree.get(&start_handle).copied().unwrap_or(0);
            let out_deg = out_degree.get(&start_handle).copied().unwrap_or(0);

            if in_deg == 1 && out_deg == 1 {
                let mut chain = vec![start_handle];
                used_handles.insert(start_handle);

                let mut current = start_handle;
                loop {
                    let Some(&next) = successors.get(&current) else {
                        break;
                    };

                    if next == start_handle {
                        // Found a cycle
                        break;
                    }

                    if used_handles.contains(&next) {
                        break;
                    }

                    let next_in = in_degree.get(&next).copied().unwrap_or(0);
                    let next_out = out_degree.get(&next).copied().unwrap_or(0);

                    if next_in == 1 && next_out == 1 {
                        if !check_path_consistency(graph, current, next) {
                            break;
                        }

                        chain.push(next);
                        used_handles.insert(next);
                        current = next;
                    } else {
                        break;
                    }
                }

                if chain.len() > 1 {
                    chains.push(chain);
                } else {
                    // Return unused handle
                    used_handles.remove(&start_handle);
                }
            }
        }
    }

    if verbose {
        eprintln!("[unchop] Found {} chains to merge", chains.len());
        let total_handles_in_chains: usize = chains.iter().map(|c| c.len()).sum();
        eprintln!(
            "[unchop] Total handles in chains: {}",
            total_handles_in_chains
        );

        // Show chain length distribution
        let mut chain_lengths: HashMap<usize, usize> = HashMap::new();
        for chain in &chains {
            *chain_lengths.entry(chain.len()).or_insert(0) += 1;
        }
        let mut sorted_lengths: Vec<_> = chain_lengths.into_iter().collect();
        sorted_lengths.sort_by_key(|&(len, _)| len);

        eprintln!("[unchop] Chain length distribution:");
        for (len, count) in sorted_lengths.iter().take(10) {
            eprintln!("  Length {}: {} chains", len, count);
        }
    }

    // Create node mapping - each chain becomes one node
    let mut node_mapping: HashMap<NodeId, NodeId> = HashMap::new();
    let mut next_new_id = 1u64;

    // Map chains to new nodes
    for chain in &chains {
        let new_id = NodeId::from(next_new_id);
        next_new_id += 1;

        for &handle in chain {
            node_mapping.insert(handle.id(), new_id);
        }
    }

    // Map unmerged nodes
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

    // Create merged nodes with concatenated sequences
    for chain in &chains {
        let new_id = node_mapping[&chain[0].id()];
        if created_nodes.insert(new_id) {
            // Concatenate sequences following the chain
            let mut seq = Vec::new();
            for &handle in chain {
                let node_handle = Handle::pack(handle.id(), false);
                if handle.is_reverse() {
                    // Need reverse complement
                    let fwd_seq: Vec<u8> = graph.sequence(node_handle).collect();
                    seq.extend(crate::utils::reverse_complement(&fwd_seq));
                } else {
                    seq.extend(graph.sequence(node_handle));
                }
            }
            new_graph.create_handle(&seq, new_id);
        }
    }

    // Create unmerged nodes
    for handle in graph.handles() {
        let node_id = handle.id();
        let new_id = node_mapping[&node_id];
        if created_nodes.insert(new_id) {
            let seq: Vec<u8> = graph.sequence(handle).collect();
            new_graph.create_handle(&seq, new_id);
        }
    }

    // Recreate edges, handling orientation changes due to merging
    let mut created_edges = HashSet::new();
    for edge in graph.edges() {
        let Edge(from_handle, to_handle) = edge;

        // Skip edges that are internal to merged chains
        let from_new = node_mapping[&from_handle.id()];
        let to_new = node_mapping[&to_handle.id()];

        if from_new == to_new {
            // This edge is now internal to a merged node
            // Check if it's a self-loop that should be preserved
            let chain = chains
                .iter()
                .find(|c| c.iter().any(|h| h.id() == from_handle.id()));
            if let Some(chain) = chain {
                // Only keep if it connects the ends of the chain
                let first = chain.first().unwrap();
                let last = chain.last().unwrap();
                if !(from_handle == *last && to_handle == *first) {
                    continue;
                }
            }
        }

        // Create the edge with proper orientations
        let new_edge = Edge(
            Handle::pack(from_new, from_handle.is_reverse()),
            Handle::pack(to_new, to_handle.is_reverse()),
        );

        let key = (new_edge.0.as_integer(), new_edge.1.as_integer());
        if created_edges.insert(key) {
            new_graph.create_edge(new_edge);
        }
    }

    // Recreate paths
    for path_id in graph.path_ids() {
        let path_name: Vec<u8> = graph.get_path_name(path_id).unwrap().collect();
        let new_path_id = new_graph
            .create_path(&path_name, false)
            .ok_or_else(|| "Failed to create path".to_string())?;

        if let Some(path_ref) = graph.get_path_ref(path_id) {
            let mut last_added: Option<(NodeId, bool)> = None;

            for step in path_ref.steps() {
                let handle = step.1;
                let new_id = node_mapping[&handle.id()];
                let is_rev = handle.is_reverse();

                // Skip duplicates
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
        eprintln!(
            "[unchop] Created new graph with {} nodes",
            new_graph.node_count()
        );
    }

    Ok(new_graph)
}

/// Check if two adjacent handles have consistent paths
fn check_path_consistency(graph: &HashGraph, from: Handle, to: Handle) -> bool {
    // For each path that visits 'from', check that it continues to 'to'
    for path_id in graph.path_ids() {
        if let Some(path_ref) = graph.get_path_ref(path_id) {
            let steps: Vec<_> = path_ref.steps().collect();

            for (i, step) in steps.iter().enumerate() {
                if step.1 == from {
                    // Found 'from' in path - check if next is 'to'
                    if i + 1 >= steps.len() || steps[i + 1].1 != to {
                        return false;
                    }
                }
            }
        }
    }

    true
}
