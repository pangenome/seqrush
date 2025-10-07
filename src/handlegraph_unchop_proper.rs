use handlegraph::handle::{Edge, Handle, NodeId};
use handlegraph::handlegraph::{HandleGraph, IntoEdges, IntoHandles, IntoSequences};
use handlegraph::hashgraph::HashGraph;
use handlegraph::mutablehandlegraph::AdditiveHandleGraph;
use handlegraph::pathhandlegraph::embedded_paths::{
    GraphPathNames, GraphPathsRef, IntoPathIds, MutableGraphPaths,
};
use handlegraph::pathhandlegraph::path::PathSteps;
use std::collections::{HashMap, HashSet};

/// Proper unchop that handles bidirected graphs correctly and checks path consistency
pub fn proper_unchop(graph: &HashGraph, verbose: bool) -> Result<HashGraph, String> {
    if verbose {
        eprintln!(
            "[unchop] Starting proper unchop with {} nodes",
            graph.node_count()
        );
    }

    // Find nodes that can potentially be merged
    // A node can be merged with its successor if:
    // 1. It has exactly one successor (considering both orientations)
    // 2. That successor has exactly one predecessor (considering both orientations)
    // 3. All paths that visit one also visit the other consecutively

    let mut mergeable_pairs: Vec<(Handle, Handle)> = Vec::new();
    let mut nodes_in_pairs = HashSet::new();

    // For each node, check if it can be merged with a successor
    for handle in graph.handles() {
        // Check both orientations of this node
        for &node_handle in &[handle, handle.flip()] {
            // Get all edges from this handle
            let mut successors: Vec<Handle> = Vec::new();
            for edge in graph.edges() {
                if edge.0 == node_handle {
                    successors.push(edge.1);
                }
            }

            // Can only merge if exactly one successor
            if successors.len() != 1 {
                continue;
            }

            let successor = successors[0];

            // Check if successor has exactly one predecessor
            let mut predecessor_count = 0;
            for edge in graph.edges() {
                if edge.1 == successor {
                    predecessor_count += 1;
                }
            }

            if predecessor_count != 1 {
                continue;
            }

            // Now check path consistency - all paths visiting node_handle must continue to successor
            let mut path_consistent = true;

            // Check all paths for consistency
            for path_id in graph.path_ids() {
                if let Some(path_ref) = graph.get_path_ref(path_id) {
                    let steps: Vec<_> = path_ref.steps().collect();

                    // Find occurrences of node_handle in the path
                    for (i, step) in steps.iter().enumerate() {
                        let handle = step.1; // Step is a tuple (StepIx, Handle)

                        if handle == node_handle {
                            // Check if this is the last step
                            if i + 1 >= steps.len() {
                                // Path ends here, can't merge
                                path_consistent = false;
                                break;
                            }

                            // Check if next step is the successor
                            let next_handle = steps[i + 1].1;
                            if next_handle != successor {
                                // Path goes somewhere else
                                path_consistent = false;
                                break;
                            }
                        }

                        if handle == successor {
                            // Check if this is the first step
                            if i == 0 {
                                // Path starts here, can't merge
                                path_consistent = false;
                                break;
                            }

                            // Check if previous step is node_handle
                            let prev_handle = steps[i - 1].1;
                            if prev_handle != node_handle {
                                // Path comes from somewhere else
                                path_consistent = false;
                                break;
                            }
                        }
                    }

                    if !path_consistent {
                        break;
                    }
                }
            }

            if path_consistent
                && !nodes_in_pairs.contains(&node_handle.id())
                && !nodes_in_pairs.contains(&successor.id())
            {
                // Found a valid merge pair!
                mergeable_pairs.push((node_handle, successor));
                nodes_in_pairs.insert(node_handle.id());
                nodes_in_pairs.insert(successor.id());

                if verbose && mergeable_pairs.len() < 10 {
                    eprintln!(
                        "[unchop] Can merge {} -> {}",
                        format!(
                            "{:?}{}",
                            node_handle.id(),
                            if node_handle.is_reverse() { "-" } else { "+" }
                        ),
                        format!(
                            "{:?}{}",
                            successor.id(),
                            if successor.is_reverse() { "-" } else { "+" }
                        )
                    );
                }

                // Don't check other orientations if we found a merge
                break;
            }
        }
    }

    if verbose {
        eprintln!(
            "[unchop] Found {} node pairs to merge",
            mergeable_pairs.len()
        );

        // Analyze why some nodes couldn't be merged
        let mut single_successor_nodes = 0;
        let mut single_predecessor_nodes = 0;
        let mut path_inconsistent_nodes = 0;

        for handle in graph.handles() {
            let mut out_count = 0;
            let mut in_count = 0;

            for edge in graph.edges() {
                if edge.0.id() == handle.id() && !edge.0.is_reverse() && !edge.1.is_reverse() {
                    out_count += 1;
                }
                if edge.1.id() == handle.id() && !edge.0.is_reverse() && !edge.1.is_reverse() {
                    in_count += 1;
                }
            }

            if out_count == 1 {
                single_successor_nodes += 1;
            }
            if in_count == 1 {
                single_predecessor_nodes += 1;
            }
        }

        eprintln!(
            "[unchop] Nodes with single successor: {}",
            single_successor_nodes
        );
        eprintln!(
            "[unchop] Nodes with single predecessor: {}",
            single_predecessor_nodes
        );
    }

    // Build chains from pairs
    let mut chains: Vec<Vec<Handle>> = Vec::new();
    let mut used_nodes = HashSet::new();

    for &(from, to) in &mergeable_pairs {
        if used_nodes.contains(&from.id()) {
            continue;
        }

        // Start a new chain
        let mut chain = vec![from, to];
        used_nodes.insert(from.id());
        used_nodes.insert(to.id());

        // Extend forward
        let mut current = to;
        loop {
            let mut found_next = false;
            for &(f, t) in &mergeable_pairs {
                if f == current && !used_nodes.contains(&t.id()) {
                    chain.push(t);
                    used_nodes.insert(t.id());
                    current = t;
                    found_next = true;
                    break;
                }
            }
            if !found_next {
                break;
            }
        }

        chains.push(chain);
    }

    if verbose {
        eprintln!("[unchop] Built {} chains from pairs", chains.len());
        let total_nodes_in_chains: usize = chains.iter().map(|c| c.len()).sum();
        eprintln!("[unchop] Total nodes in chains: {}", total_nodes_in_chains);
    }

    // Create node mapping
    let mut node_mapping: HashMap<NodeId, NodeId> = HashMap::new();
    let mut next_new_id = 1u64;

    // Assign new IDs to chains
    for chain in &chains {
        let new_id = NodeId::from(next_new_id);
        next_new_id += 1;

        for &handle in chain {
            node_mapping.insert(handle.id(), new_id);
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
        let new_id = node_mapping[&chain[0].id()];
        if created_nodes.insert(new_id) {
            // Concatenate sequences in the correct order
            let mut seq = Vec::new();
            for &handle in chain {
                let node_handle = Handle::pack(handle.id(), false);
                if handle.is_reverse() {
                    // If the handle is reverse, we need to reverse complement the sequence
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

    // Create edges - need to handle orientation correctly
    let mut created_edges = HashSet::new();
    for edge in graph.edges() {
        let Edge(from_handle, to_handle) = edge;

        let new_from_id = node_mapping[&from_handle.id()];
        let new_to_id = node_mapping[&to_handle.id()];

        // Skip edges within merged chains
        if new_from_id == new_to_id && !from_handle.is_reverse() && !to_handle.is_reverse() {
            continue;
        }

        let new_edge = Edge(
            Handle::pack(new_from_id, from_handle.is_reverse()),
            Handle::pack(new_to_id, to_handle.is_reverse()),
        );

        let key = (new_edge.0.as_integer(), new_edge.1.as_integer());
        if created_edges.insert(key) {
            new_graph.create_edge(new_edge);
        }
    }

    // Recreate paths with proper deduplication
    for path_id in graph.path_ids() {
        let path_name: Vec<u8> = graph.get_path_name(path_id).unwrap().collect();
        let new_path_id = new_graph
            .create_path(&path_name, false)
            .ok_or_else(|| "Failed to create path".to_string())?;

        let mut last_added: Option<(NodeId, bool)> = None;

        if let Some(path_ref) = graph.get_path_ref(path_id) {
            for step in path_ref.steps() {
                let handle = step.1;
                let old_id = handle.id();
                let new_id = node_mapping[&old_id];
                let is_rev = handle.is_reverse();

                // Skip consecutive steps that map to same node with same orientation
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
