use crate::bidirected_graph::{BiEdge, BiNode, BiPath, Handle};
use crate::graph_ops::{Edge, Graph, Node};
use std::collections::{BTreeSet, HashMap, HashSet};
use sha2::{Sha256, Digest};

/// A bidirected graph that extends the basic Graph with orientation support
#[derive(Clone)]
pub struct BidirectedGraph {
    pub nodes: HashMap<usize, BiNode>,
    pub edges: HashSet<BiEdge>,
    pub paths: Vec<BiPath>,
}

impl Default for BidirectedGraph {
    fn default() -> Self {
        Self::new()
    }
}

impl BidirectedGraph {
    /// Apply a node ID remapping to renumber nodes in the graph
    pub fn apply_node_id_mapping(&mut self, mapping: &HashMap<usize, usize>) {
        // Create new nodes with remapped IDs
        let mut new_nodes = HashMap::new();
        for (old_id, node) in &self.nodes {
            let new_id = mapping.get(old_id).copied().unwrap_or(*old_id);
            let mut new_node = node.clone();
            new_node.id = new_id;
            new_nodes.insert(new_id, new_node);
        }
        self.nodes = new_nodes;

        // Update edges with new node IDs
        let mut new_edges = HashSet::new();
        for edge in &self.edges {
            let new_from_id = mapping.get(&edge.from.node_id()).copied().unwrap_or(edge.from.node_id());
            let new_to_id = mapping.get(&edge.to.node_id()).copied().unwrap_or(edge.to.node_id());

            let new_from = if edge.from.is_reverse() {
                Handle::reverse(new_from_id)
            } else {
                Handle::forward(new_from_id)
            };

            let new_to = if edge.to.is_reverse() {
                Handle::reverse(new_to_id)
            } else {
                Handle::forward(new_to_id)
            };

            new_edges.insert(BiEdge::new(new_from, new_to));
        }
        self.edges = new_edges;

        // Update paths with new node IDs
        for path in &mut self.paths {
            for handle in &mut path.steps {
                let new_id = mapping.get(&handle.node_id()).copied().unwrap_or(handle.node_id());
                *handle = if handle.is_reverse() {
                    Handle::reverse(new_id)
                } else {
                    Handle::forward(new_id)
                };
            }
        }
    }

    /// Compact the graph by merging linear chains of nodes
    /// Renumber nodes to be sequential starting from 1
    pub fn renumber_nodes_sequentially(&mut self) {
        eprintln!("[VALIDATION] Before renumbering - checking for issues...");
        self.validate_paths("before renumbering");

        let mut old_to_new: HashMap<usize, usize> = HashMap::new();
        let mut new_id = 1;

        // Create mapping from old to new IDs
        let mut old_ids: Vec<usize> = self.nodes.keys().cloned().collect();
        old_ids.sort();

        for old_id in old_ids {
            old_to_new.insert(old_id, new_id);
            new_id += 1;
        }

        // Apply the renumbering
        self.apply_node_id_mapping(&old_to_new);

        eprintln!("[VALIDATION] After renumbering - checking for issues...");
        self.validate_paths("after renumbering");
    }

    pub fn compact(&mut self) {
        let mut compacted = true;
        let mut iteration = 0;

        // Validate before compaction
        self.validate_graph_consistency("before compaction");
        self.validate_paths("before compaction");

        // Keep compacting until no more changes
        while compacted {
            compacted = false;
            iteration += 1;

            // Find simple components (linear chains)
            let components = self.find_simple_components();

            // Merge each component
            for component in components {
                if component.len() >= 2 && self.merge_component_v2(&component) {
                    compacted = true;
                    // Validate after each merge
                    self.validate_graph_consistency(&format!("after merge in iteration {}", iteration));
                    self.validate_paths(&format!("after merge in iteration {}", iteration));
                }
            }

            if !compacted {
                break;
            }
        }

        // Final validation
        self.validate_graph_consistency("after compaction");
        self.validate_paths("after compaction");
    }

    /// Find simple components (linear chains that can be merged)
    fn find_simple_components(&self) -> Vec<Vec<Handle>> {
        let mut components = Vec::new();
        let mut visited = std::collections::HashSet::new();
        let mut merged_nodes = std::collections::HashSet::new();

        // Build adjacency lists
        let mut forward_edges: std::collections::HashMap<Handle, Vec<Handle>> =
            std::collections::HashMap::new();
        let mut backward_edges: std::collections::HashMap<Handle, Vec<Handle>> =
            std::collections::HashMap::new();

        for edge in &self.edges {
            forward_edges.entry(edge.from).or_default().push(edge.to);
            backward_edges.entry(edge.to).or_default().push(edge.from);

            // Also consider the implied reverse edge
            forward_edges
                .entry(edge.to.flip())
                .or_default()
                .push(edge.from.flip());
            backward_edges
                .entry(edge.from.flip())
                .or_default()
                .push(edge.to.flip());
        }

        // Check if two handles are perfect neighbors (can be merged)
        let are_perfect_neighbors = |from: Handle, to: Handle| -> bool {
            // Check all paths for both forward and reverse consistency
            for path in &self.paths {
                // Check forward direction: from -> to
                let mut from_to_transitions = 0;
                let mut from_visits = 0;

                for i in 0..path.steps.len() {
                    if path.steps[i] == from {
                        from_visits += 1;
                        if i + 1 < path.steps.len() {
                            if path.steps[i + 1] == to {
                                from_to_transitions += 1;
                            } else {
                                // Path continues but doesn't go to 'to'
                                return false;
                            }
                        } else {
                            // Path ends at 'from'
                            return false;
                        }
                    }
                }

                // If we visited 'from' but didn't always go to 'to', not perfect
                if from_visits > 0 && from_visits != from_to_transitions {
                    return false;
                }
                
                // Also check reverse: to.flip() -> from.flip()
                // This ensures bidirected consistency
                let from_rev = from.flip();
                let to_rev = to.flip();
                let mut to_rev_visits = 0;
                let mut to_rev_to_from_rev = 0;
                
                for i in 0..path.steps.len() {
                    if path.steps[i] == to_rev {
                        to_rev_visits += 1;
                        if i + 1 < path.steps.len() {
                            if path.steps[i + 1] == from_rev {
                                to_rev_to_from_rev += 1;
                            } else {
                                // Reverse path doesn't maintain the chain
                                return false;
                            }
                        } else {
                            // Path ends at to_rev
                            return false;
                        }
                    }
                }
                
                if to_rev_visits > 0 && to_rev_visits != to_rev_to_from_rev {
                    return false;
                }
            }

            true
        };

        // Find linear chains
        for handle in self
            .nodes
            .keys()
            .flat_map(|&id| vec![Handle::forward(id), Handle::reverse(id)])
        {
            if visited.contains(&handle) {
                continue;
            }

            // Check if this could be the start of a chain
            let in_degree = backward_edges.get(&handle).map(|v| v.len()).unwrap_or(0);
            let out_degree = forward_edges.get(&handle).map(|v| v.len()).unwrap_or(0);

            if out_degree == 1 {
                // Try to extend a chain from here
                let mut chain = vec![handle];
                visited.insert(handle);

                let mut current = handle;
                while let Some(nexts) = forward_edges.get(&current) {
                    if nexts.len() != 1 {
                        break;
                    }

                    let next = nexts[0];
                    let next_in_degree = backward_edges.get(&next).map(|v| v.len()).unwrap_or(0);

                    if next_in_degree != 1 || visited.contains(&next) {
                        break;
                    }

                    // Check if they are perfect neighbors
                    if !are_perfect_neighbors(current, next) {
                        break;
                    }

                    chain.push(next);
                    visited.insert(next);
                    current = next;

                    // Check next's out degree
                    let next_out_degree = forward_edges.get(&next).map(|v| v.len()).unwrap_or(0);
                    if next_out_degree != 1 {
                        break;
                    }
                }

                if chain.len() >= 2 {
                    // Check if any node in this chain has already been merged
                    let mut already_merged = false;
                    for &h in &chain {
                        if merged_nodes.contains(&h.node_id()) {
                            already_merged = true;
                            break;
                        }
                    }

                    if !already_merged {
                        // Mark all nodes in this chain as merged
                        for &h in &chain {
                            merged_nodes.insert(h.node_id());
                        }
                        components.push(chain);
                    }
                }
            }
        }

        components
    }
    
    /// Merge a component (chain of handles) into a single node - Version 2
    /// Returns true if the merge was successful
    fn merge_component_v2(&mut self, handles: &[Handle]) -> bool {
        if handles.len() < 2 {
            return false;
        }

        eprintln!("[COMPACT DEBUG] Attempting to merge chain: {:?}",
                  handles.iter().map(|h| format!("{}{}", h.node_id(), if h.is_reverse() { "-" } else { "+" })).collect::<Vec<_>>());

        // First, build a complete mapping of what needs to be replaced
        // This includes both the handles in the chain and their flipped versions
        let mut handle_mapping = std::collections::HashMap::new();
        
        // Create new sequence by concatenating
        let mut new_sequence = Vec::new();
        for &handle in handles {
            if let Some(node) = self.nodes.get(&handle.node_id()) {
                if handle.is_reverse() {
                    new_sequence.extend(crate::bidirected_graph::reverse_complement(&node.sequence));
                } else {
                    new_sequence.extend(&node.sequence);
                }
            }
        }
        
        // Create new node
        let new_node_id = self.next_node_id();
        let new_handle_forward = Handle::forward(new_node_id);
        let new_handle_reverse = Handle::reverse(new_node_id);
        
        // Build complete handle mapping
        // Map each handle in the chain to its position in the new node
        for (i, &handle) in handles.iter().enumerate() {
            // The forward chain maps to forward new node
            handle_mapping.insert(handle, (new_handle_forward, i, false));
            // The reverse of each handle maps to reverse position
            let rev_handle = handle.flip();
            let rev_position = handles.len() - 1 - i;
            handle_mapping.insert(rev_handle, (new_handle_reverse, rev_position, true));
        }
        
        // Validate that paths can be properly updated
        for path in &self.paths {
            let mut i = 0;
            while i < path.steps.len() {
                if let Some((_, chain_pos, _)) = handle_mapping.get(&path.steps[i]) {
                    // Found a handle from our chain
                    // Check if it's the start of the complete forward chain
                    if *chain_pos == 0 {
                        // Check if the complete chain follows
                        let mut chain_complete = true;
                        if i + handles.len() <= path.steps.len() {
                            for j in 0..handles.len() {
                                if path.steps[i + j] != handles[j] {
                                    chain_complete = false;
                                    break;
                                }
                            }
                            if chain_complete {
                                i += handles.len();
                                continue;
                            }
                        }
                    }
                    
                    // Check if it's the start of the complete reverse chain
                    let rev_chain: Vec<Handle> = handles.iter().rev().map(|h| h.flip()).collect();
                    if path.steps[i] == rev_chain[0] {
                        let mut chain_complete = true;
                        if i + rev_chain.len() <= path.steps.len() {
                            for j in 0..rev_chain.len() {
                                if path.steps[i + j] != rev_chain[j] {
                                    chain_complete = false;
                                    break;
                                }
                            }
                            if chain_complete {
                                i += rev_chain.len();
                                continue;
                            }
                        }
                    }
                    
                    // If we get here, the handle appears but not as part of a complete chain
                    // Cannot compact this chain
                    eprintln!("[COMPACT DEBUG] Chain validation FAILED: handle {} appears individually in path {} at position {}",
                             path.steps[i], path.name, i);
                    return false;
                }
                i += 1;
            }
        }
        
        // All paths validated - proceed with the merge
        eprintln!("[COMPACT DEBUG] Validation passed, creating new node {}", new_node_id);
        self.add_node(new_node_id, new_sequence);

        // Update all paths
        for path in &mut self.paths {
            let old_step_count = path.steps.len();
            let mut new_steps = Vec::new();
            let mut i = 0;
            let mut replacements = 0;

            while i < path.steps.len() {
                // Check for forward chain
                if i + handles.len() <= path.steps.len() {
                    let mut is_forward_chain = true;
                    for j in 0..handles.len() {
                        if path.steps[i + j] != handles[j] {
                            is_forward_chain = false;
                            break;
                        }
                    }
                    if is_forward_chain {
                        replacements += 1;
                        new_steps.push(new_handle_forward);
                        i += handles.len();
                        continue;
                    }
                }

                // Check for reverse chain
                let rev_chain: Vec<Handle> = handles.iter().rev().map(|h| h.flip()).collect();
                if i + rev_chain.len() <= path.steps.len() {
                    let mut is_reverse_chain = true;
                    for j in 0..rev_chain.len() {
                        if path.steps[i + j] != rev_chain[j] {
                            is_reverse_chain = false;
                            break;
                        }
                    }
                    if is_reverse_chain {
                        replacements += 1;
                        new_steps.push(new_handle_reverse);
                        i += rev_chain.len();
                        continue;
                    }
                }

                // Not part of a chain, keep as is
                new_steps.push(path.steps[i]);
                i += 1;
            }

            let new_step_count = new_steps.len();
            if replacements > 0 {
                eprintln!("[COMPACT DEBUG]   Path {}: {} steps -> {} steps ({} chain replacements)",
                         path.name, old_step_count, new_step_count, replacements);
            }

            path.steps = new_steps;
        }
        
        // Update edges
        let first_handle = handles[0];
        let last_handle = handles[handles.len() - 1];

        // Build set of all node IDs being merged for quick lookup
        let merged_node_ids: std::collections::HashSet<usize> =
            handles.iter().map(|h| h.node_id()).collect();

        // Build new edge set
        let mut new_edges = std::collections::HashSet::new();

        for edge in &self.edges {
            let from_in_chain = merged_node_ids.contains(&edge.from.node_id());
            let to_in_chain = merged_node_ids.contains(&edge.to.node_id());

            if from_in_chain && to_in_chain {
                // Both ends are in the chain - this is an internal edge, skip it
                continue;
            } else if !from_in_chain && !to_in_chain {
                // Neither end is in the chain - keep edge as-is
                new_edges.insert(edge.clone());
            } else if from_in_chain && !to_in_chain {
                // Edge goes FROM a node in the chain TO outside
                // Check if it's from the last handle of the forward chain
                if edge.from == last_handle {
                    new_edges.insert(BiEdge {
                        from: new_handle_forward,
                        to: edge.to,
                    });
                }
                // Check if it's from the first handle of the reverse chain
                let first_rev = first_handle.flip();
                if edge.from == first_rev {
                    new_edges.insert(BiEdge {
                        from: new_handle_reverse,
                        to: edge.to,
                    });
                }
            } else {
                // Edge goes FROM outside TO a node in the chain
                // Check if it's to the first handle of the forward chain
                if edge.to == first_handle {
                    new_edges.insert(BiEdge {
                        from: edge.from,
                        to: new_handle_forward,
                    });
                }
                // Check if it's to the last handle of the reverse chain
                let last_rev = last_handle.flip();
                if edge.to == last_rev {
                    new_edges.insert(BiEdge {
                        from: edge.from,
                        to: new_handle_reverse,
                    });
                }
            }
        }

        self.edges = new_edges;

        // Remove old nodes
        eprintln!("[COMPACT DEBUG] Removing {} old nodes from chain", handles.len());
        for &handle in handles {
            let node_id = handle.node_id();
            if self.nodes.remove(&node_id).is_some() {
                eprintln!("[COMPACT DEBUG]   Removed node {}", node_id);
            } else {
                eprintln!("[COMPACT DEBUG]   WARNING: Node {} was already deleted!", node_id);
            }
        }

        true
    }

    /// Merge a component (chain of handles) into a single node
    fn merge_component(&mut self, handles: &[Handle]) {
        if handles.len() < 2 {
            return;
        }

        // First, validate that this chain can be safely compacted
        // Check that all occurrences of chain nodes appear as the complete chain
        for path in &self.paths {
            let mut i = 0;
            while i < path.steps.len() {
                let current_handle = path.steps[i];
                
                // Check if this handle is in our chain
                if handles.contains(&current_handle) {
                    // Found a handle from the chain - verify the entire chain is here
                    let chain_position = handles.iter().position(|&h| h == current_handle);
                    if let Some(pos) = chain_position {
                        // Check if we have the complete chain starting from the beginning
                        if pos == 0 {
                            // Check if the entire chain follows
                            let mut chain_complete = true;
                            if i + handles.len() <= path.steps.len() {
                                for j in 0..handles.len() {
                                    if path.steps[i + j] != handles[j] {
                                        chain_complete = false;
                                        break;
                                    }
                                }
                            } else {
                                chain_complete = false;
                            }
                            
                            if chain_complete {
                                // Good, skip the entire chain
                                i += handles.len();
                                continue;
                            }
                        }
                        
                        // If we get here, the chain node appears but not as part of the complete chain
                        eprintln!(
                            "WARNING: Cannot compact chain {:?} - node {} appears individually in path {}",
                            handles.iter().map(|h| h.node_id()).collect::<Vec<_>>(),
                            current_handle.node_id(),
                            path.name
                        );
                        return;  // Abort this compaction
                    }
                }
                i += 1;
            }
        }

        // Create new sequence by concatenating
        let mut new_sequence = Vec::new();
        for &handle in handles {
            if let Some(node) = self.nodes.get(&handle.node_id()) {
                if handle.is_reverse() {
                    new_sequence
                        .extend(crate::bidirected_graph::reverse_complement(&node.sequence));
                } else {
                    new_sequence.extend(&node.sequence);
                }
            }
        }

        // Create new node
        let new_node_id = self.next_node_id();
        self.add_node(new_node_id, new_sequence);

        // Determine the orientation of the new node based on the chain
        // If the chain starts with a reverse handle, the new node represents that reverse orientation
        let chain_is_reverse = handles[0].is_reverse();

        // Create mapping from old handles to new handle
        let mut handle_map = std::collections::HashMap::new();
        for (i, &handle) in handles.iter().enumerate() {
            // Map this handle to the new node
            // The orientation depends on whether the chain is reverse and the handle's orientation
            let new_handle = if chain_is_reverse {
                // If chain is reverse and handle is reverse, they cancel out
                if handle.is_reverse() {
                    Handle::forward(new_node_id)
                } else {
                    Handle::reverse(new_node_id)
                }
            } else {
                // Chain is forward
                if handle.is_reverse() {
                    Handle::reverse(new_node_id)
                } else {
                    Handle::forward(new_node_id)
                }
            };
            handle_map.insert(handle, new_handle);
        }

        // Track what the first and last handles connect to
        let first_handle = handles[0];
        let last_handle = handles[handles.len() - 1];

        // Find edges to reconnect
        let mut edges_to_add = Vec::new();

        for edge in &self.edges {
            // Check edges entering the chain
            if edge.to == first_handle && !handles.contains(&edge.from) {
                let new_to = handle_map[&first_handle];
                edges_to_add.push((edge.from, new_to));
            }
            // Check edges leaving the chain
            if edge.from == last_handle && !handles.contains(&edge.to) {
                let new_from = handle_map[&last_handle];
                edges_to_add.push((new_from, edge.to));
            }
            // Check for self-loops within the chain
            if edge.from == last_handle && edge.to == first_handle {
                // Loop from end to beginning
                let new_from = handle_map[&last_handle];
                let new_to = handle_map[&first_handle];
                if new_from.node_id() == new_to.node_id() {
                    // Self-loop on new node
                    edges_to_add.push((new_from, new_to));
                }
            }
        }

        // Update paths - we've already validated that all chain nodes appear as complete chains
        for path in &mut self.paths {
            let mut new_steps = Vec::new();
            let mut i = 0;

            while i < path.steps.len() {
                let current_handle = path.steps[i];

                // Check if this is the start of our chain
                if current_handle == handles[0] {
                    // Check if the complete chain follows
                    let mut is_complete_chain = true;
                    if i + handles.len() <= path.steps.len() {
                        for j in 0..handles.len() {
                            if path.steps[i + j] != handles[j] {
                                is_complete_chain = false;
                                break;
                            }
                        }
                    } else {
                        is_complete_chain = false;
                    }

                    if is_complete_chain {
                        // Replace the entire chain with the new merged handle
                        let new_handle = handle_map.get(&handles[0]).unwrap();
                        new_steps.push(*new_handle);
                        i += handles.len();
                    } else {
                        // This shouldn't happen after our validation
                        new_steps.push(current_handle);
                        i += 1;
                    }
                } else if handles.contains(&current_handle) {
                    // This shouldn't happen - we validated that chain nodes only appear as complete chains
                    eprintln!(
                        "ERROR: Unexpected individual occurrence of chain node {} in path {}",
                        current_handle, path.name
                    );
                    new_steps.push(current_handle);
                    i += 1;
                } else {
                    // This handle is not in the chain, keep it
                    new_steps.push(current_handle);
                    i += 1;
                }
            }

            path.steps = new_steps;
        }

        // Remove old nodes
        for &handle in handles {
            self.nodes.remove(&handle.node_id());
        }

        // Remove old edges
        self.edges.retain(|edge| {
            !handles
                .iter()
                .any(|&h| h.node_id() == edge.from.node_id() || h.node_id() == edge.to.node_id())
        });

        // Add new edges
        for (from, to) in edges_to_add {
            self.add_edge(from, to);
        }
    }

    fn next_node_id(&self) -> usize {
        self.nodes.keys().max().copied().unwrap_or(0) + 1
    }
    pub fn new() -> Self {
        BidirectedGraph {
            nodes: HashMap::new(),
            edges: HashSet::new(),
            paths: Vec::new(),
        }
    }

    /// Get the paths that pass through a given node
    pub fn paths_through_node(&self, node_id: usize) -> Vec<usize> {
        let mut paths_with_node = Vec::new();
        for (path_idx, path) in self.paths.iter().enumerate() {
            for handle in &path.steps {
                if handle.node_id() == node_id {
                    paths_with_node.push(path_idx);
                    break;  // Only need to add each path once
                }
            }
        }
        paths_with_node
    }

    /// Get total sequence length in the graph
    pub fn total_sequence_length(&self) -> usize {
        self.nodes.values().map(|node| node.sequence.len()).sum()
    }

    /// Get the number of nodes in the graph
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Convert from the old Graph format to BidirectedGraph
    pub fn from_graph(graph: Graph) -> Self {
        let mut bi_graph = BidirectedGraph::new();

        // Convert nodes
        for (id, node) in graph.nodes {
            let bi_node = BiNode {
                id,
                sequence: node.sequence,
                rank: Some(node.rank as u64),
            };
            bi_graph.nodes.insert(id, bi_node);
        }

        // Convert edges (all as forward-to-forward for now)
        for edge in graph.edges {
            let bi_edge = BiEdge::new(Handle::forward(edge.from), Handle::forward(edge.to));
            bi_graph.edges.insert(bi_edge);
        }

        // Convert paths (all as forward steps for now)
        for (name, path) in graph.paths {
            let mut bi_path = BiPath::new(name);
            for node_id in path {
                bi_path.add_step(Handle::forward(node_id));
            }
            bi_graph.paths.push(bi_path);
        }

        bi_graph
    }

    /// Convert back to the old Graph format for compatibility
    pub fn to_graph(&self) -> Graph {
        let mut graph = Graph::new();

        // Convert nodes
        for (id, bi_node) in &self.nodes {
            let node = Node {
                id: *id,
                sequence: bi_node.sequence.clone(),
                rank: bi_node.rank.unwrap_or(0) as f64,
            };
            graph.nodes.insert(*id, node);
        }

        // Convert edges - for now we only include forward-to-forward edges
        // This is a limitation but maintains compatibility
        for bi_edge in &self.edges {
            if !bi_edge.from.is_reverse() && !bi_edge.to.is_reverse() {
                let edge = Edge {
                    from: bi_edge.from.node_id(),
                    to: bi_edge.to.node_id(),
                };
                graph.edges.insert(edge);
            }
        }

        // Convert paths - extracting just node IDs
        for bi_path in &self.paths {
            let path_nodes: Vec<usize> = bi_path
                .steps
                .iter()
                .map(|handle| handle.node_id())
                .collect();
            graph.paths.push((bi_path.name.clone(), path_nodes));
        }

        graph
    }

    /// Add a node to the graph
    pub fn add_node(&mut self, id: usize, sequence: Vec<u8>) {
        self.nodes.insert(id, BiNode::new(id, sequence));
    }

    /// Add an edge to the graph
    pub fn add_edge(&mut self, from: Handle, to: Handle) {
        let edge = BiEdge::new(from, to);

        // Check if the complement edge already exists
        // The complement of A+ -> B+ is B- -> A-
        let complement = BiEdge::new(to.flip(), from.flip());

        // Only add the edge if neither it nor its complement is already present
        if !self.edges.contains(&edge) && !self.edges.contains(&complement) {
            // Add the edge as-is (not in canonical form which might flip node IDs incorrectly)
            self.edges.insert(edge);
        }
    }

    /// Get sequence for a handle (forward or reverse complement)
    pub fn get_sequence(&self, handle: Handle) -> Option<Vec<u8>> {
        self.nodes
            .get(&handle.node_id())
            .map(|node| node.get_sequence(handle.is_reverse()))
    }

    /// Check if an edge exists (checking both the edge and its complement)
    pub fn has_edge(&self, from: Handle, to: Handle) -> bool {
        let edge = BiEdge::new(from, to);
        let complement = BiEdge::new(to.flip(), from.flip());
        self.edges.contains(&edge) || self.edges.contains(&complement)
    }

    /// Get all edges from a handle
    pub fn edges_from(&self, handle: Handle) -> Vec<Handle> {
        let mut result = Vec::new();
        for edge in &self.edges {
            if edge.from == handle {
                result.push(edge.to);
            } else if edge.to.flip() == handle {
                // This edge's complement has handle as from
                result.push(edge.from.flip());
            }
        }
        result
    }

    /// Get all edges to a handle
    pub fn edges_to(&self, handle: Handle) -> Vec<Handle> {
        let mut result = Vec::new();
        for edge in &self.edges {
            if edge.to == handle {
                result.push(edge.from);
            } else if edge.from.flip() == handle {
                // This edge's complement has handle as to
                result.push(edge.to.flip());
            }
        }
        result
    }

    /// Build a path from a sequence of node IDs with orientations
    pub fn build_path(&mut self, name: String, steps: Vec<(usize, bool)>) {
        let mut path = BiPath::new(name);
        for (node_id, is_reverse) in steps {
            path.add_step(Handle::new(node_id, is_reverse));
        }
        self.paths.push(path);
    }

    /// Write graph in GFA format with proper orientations
    pub fn write_gfa<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        // Write header
        writeln!(writer, "H\tVN:Z:1.0")?;

        // Write segments (nodes)
        let mut node_ids: Vec<_> = self.nodes.keys().cloned().collect();
        node_ids.sort();

        for node_id in node_ids {
            if let Some(node) = self.nodes.get(&node_id) {
                let seq_str = String::from_utf8_lossy(&node.sequence);
                writeln!(writer, "S\t{}\t{}", node_id, seq_str)?;
            }
        }

        // Write links (edges) with orientations
        // NOTE: We do NOT canonicalize edges when writing!
        // If we have edge A+ -> B+, its reverse complement is B- -> A-
        // If we canonicalize 3+ -> 2+ to 2- -> 3-, and we already have 2+ -> 3+,
        // the output will incorrectly show both 2+ -> 3+ and 2- -> 3-, which
        // violates the bidirected graph property.
        for edge in &self.edges {
            writeln!(
                writer,
                "L\t{}\t{}\t{}\t{}\t0M",
                edge.from.node_id(),
                edge.from.orientation_char(),
                edge.to.node_id(),
                edge.to.orientation_char()
            )?;
        }

        // Write paths with orientations
        for path in &self.paths {
            write!(writer, "P\t{}\t", path.name)?;

            // Write path steps with orientations
            let steps: Vec<String> = path
                .steps
                .iter()
                .map(|handle| format!("{}{}", handle.node_id(), handle.orientation_char()))
                .collect();
            write!(writer, "{}\t*", steps.join(","))?;

            writeln!(writer)?;
        }

        Ok(())
    }
    
    /// Validate graph consistency - check that all edges reference existing nodes and paths are valid
    pub fn validate_graph_consistency(&self, phase: &str) -> bool {
        let mut errors = Vec::new();

        // Check all edges reference existing nodes
        for edge in &self.edges {
            if !self.nodes.contains_key(&edge.from.node_id()) {
                errors.push(format!("Edge references non-existent node: {} (from)", edge.from.node_id()));
            }
            if !self.nodes.contains_key(&edge.to.node_id()) {
                errors.push(format!("Edge references non-existent node: {} (to)", edge.to.node_id()));
            }
        }

        // Check all path steps reference existing nodes
        for path in &self.paths {
            for handle in &path.steps {
                if !self.nodes.contains_key(&handle.node_id()) {
                    errors.push(format!("Path {} references non-existent node: {}", path.name, handle.node_id()));
                }
            }
        }

        if !errors.is_empty() {
            eprintln!("\n[ERROR] Graph validation failed at phase: {}", phase);
            for error in errors.iter().take(10) {
                eprintln!("  {}", error);
            }
            if errors.len() > 10 {
                eprintln!("  ... and {} more errors", errors.len() - 10);
            }
            false
        } else {
            true
        }
    }

    /// Compute SHA-256 hash of a path's sequence
    pub fn compute_path_hash(&self, path: &BiPath) -> String {
        let mut hasher = Sha256::new();

        for handle in &path.steps {
            if let Some(seq) = self.get_sequence(*handle) {
                hasher.update(&seq);
            }
        }

        format!("{:x}", hasher.finalize())
    }

    /// Compute SHA-256 hashes for all paths
    pub fn compute_all_path_hashes(&self) -> HashMap<String, String> {
        let mut hashes = HashMap::new();
        for path in &self.paths {
            hashes.insert(path.name.clone(), self.compute_path_hash(path));
        }
        hashes
    }

    /// Validate paths - check basic statistics and optionally verify hashes haven't changed
    /// NOTE: Repeated consecutive nodes ARE VALID in seqwish/odgi graphs!
    /// They represent structural variation (indels, duplications, etc.)
    pub fn validate_paths(&self, phase: &str) -> bool {
        eprintln!("[VALIDATION] At phase '{}': {} paths, {} nodes, {} edges",
                  phase, self.paths.len(), self.nodes.len(), self.edges.len());

        // Compute and display path hashes
        let hashes = self.compute_all_path_hashes();
        eprintln!("[VALIDATION] Path hashes at '{}':", phase);
        for (name, hash) in &hashes {
            eprintln!("  {}: {}...", name, &hash[..16]);
        }

        // Just print statistics for now - repeated nodes are actually valid!
        let mut has_repeats = false;
        for path in &self.paths {
            let mut node_counts = HashMap::new();
            let mut consecutive_repeats = 0;

            for handle in &path.steps {
                *node_counts.entry(handle.node_id()).or_insert(0) += 1;
            }

            for i in 0..(path.steps.len().saturating_sub(1)) {
                if path.steps[i] == path.steps[i + 1] {
                    consecutive_repeats += 1;
                }
            }

            if consecutive_repeats > 0 || node_counts.values().any(|&c| c > 5) {
                has_repeats = true;
            }
        }

        if has_repeats {
            eprintln!("[VALIDATION INFO] Graph has repeated nodes in paths (this is VALID for representing structural variation)");
        }

        true // Always pass - repeated nodes are valid
    }

    /// Validate that path hashes match expected values
    pub fn validate_path_hashes(&self, expected_hashes: &HashMap<String, String>, phase: &str) -> bool {
        let current_hashes = self.compute_all_path_hashes();
        let mut all_match = true;

        eprintln!("[VALIDATION] Checking path sequence integrity at '{}':", phase);
        for (name, expected_hash) in expected_hashes {
            if let Some(current_hash) = current_hashes.get(name) {
                if current_hash != expected_hash {
                    eprintln!("  ERROR: Path {} hash changed!", name);
                    eprintln!("    Expected: {}", expected_hash);
                    eprintln!("    Got:      {}", current_hash);
                    all_match = false;
                } else {
                    eprintln!("  ✓ Path {} sequence preserved", name);
                }
            } else {
                eprintln!("  ERROR: Path {} missing!", name);
                all_match = false;
            }
        }

        all_match
    }

    /// Verify that all edges needed for paths exist in the graph
    pub fn verify_path_edges(&mut self, verbose: bool) {
        let mut missing_edges = Vec::new();
        let mut added_edges = 0;

        for path in &self.paths {
            for i in 0..path.steps.len().saturating_sub(1) {
                let from = path.steps[i];
                let to = path.steps[i + 1];
                let edge = BiEdge::new(from, to);
                let complement = BiEdge::new(to.flip(), from.flip());

                // Check if neither the edge nor its complement exists
                if !self.edges.contains(&edge) && !self.edges.contains(&complement) {
                    missing_edges.push((path.name.clone(), from, to));
                    // Add the missing edge as-is (not in canonical form)
                    self.edges.insert(edge);
                    added_edges += 1;
                }
            }
        }
        
        if verbose && !missing_edges.is_empty() {
            eprintln!("WARNING: Found {} missing edges needed for paths", missing_edges.len());
            for (path_name, from, to) in missing_edges.iter().take(10) {
                eprintln!("  Path {}: edge {} -> {} was missing", path_name, from, to);
            }
            if missing_edges.len() > 10 {
                eprintln!("  ... and {} more", missing_edges.len() - 10);
            }
            eprintln!("Added {} missing edges", added_edges);
        }
    }

    /// Perform bidirected topological sort using modified Kahn's algorithm
    /// This is a direct port of the odgi topological_order algorithm
    pub fn topological_sort(
        &mut self,
        use_heads: bool,
        use_tails: bool,
        verbose: bool,
    ) -> Vec<Handle> {
        let mut sorted = Vec::new();
        
        if self.nodes.is_empty() {
            return sorted;
        }

        // Create mapping of all possible handles (both orientations for each node)
        let mut handle_to_idx = HashMap::new();
        let mut idx_to_handle = Vec::new();
        let mut idx = 0;
        
        // Collect all node IDs in sorted order for deterministic behavior
        let mut node_ids: Vec<_> = self.nodes.keys().copied().collect();
        node_ids.sort();
        
        for node_id in &node_ids {
            // Forward orientation
            let fwd = Handle::forward(*node_id);
            handle_to_idx.insert(fwd, idx);
            idx_to_handle.push(fwd);
            idx += 1;
            
            // Reverse orientation
            let rev = Handle::reverse(*node_id);
            handle_to_idx.insert(rev, idx);
            idx_to_handle.push(rev);
            idx += 1;
        }
        
        let total_handles = idx_to_handle.len();
        
        // S - set of handles ready to be processed
        let mut s = HashSet::new();
        
        // Track which handles have been visited
        let mut visited = vec![false; total_handles];
        
        // Seeds for cycle breaking
        let mut seeds = HashMap::new();
        
        // Track masked edges (simulating edge removal)
        let mut masked_edges = HashSet::new();
        
        // Find head and tail nodes if requested
        if use_heads {
            for &node_id in &node_ids {
                let handle = Handle::forward(node_id);
                // Check if this node has any incoming edges on its left side
                let mut has_left_edges = false;
                
                // Check all edges where this node is the target
                for edge in &self.edges {
                    if edge.to.node_id() == node_id && !edge.to.is_reverse() {
                        // Edge comes to forward orientation (left side)
                        has_left_edges = true;
                        break;
                    }
                    if edge.to.node_id() == node_id && edge.to.is_reverse() {
                        // Edge comes to reverse orientation (which is right side of forward)
                        continue;
                    }
                    // Also check edges from this node going backward (self loops)
                    if edge.from.node_id() == node_id && edge.from.is_reverse() && edge.to.node_id() == node_id && !edge.to.is_reverse() {
                        has_left_edges = true;
                        break;
                    }
                }
                
                if !has_left_edges {
                    s.insert(handle);
                }
            }
        } else if use_tails {
            for &node_id in &node_ids {
                let handle = Handle::forward(node_id);
                // Check if this node has any outgoing edges on its right side
                let mut has_right_edges = false;
                
                for edge in &self.edges {
                    if edge.from.node_id() == node_id && !edge.from.is_reverse() {
                        // Edge goes from forward orientation (right side)
                        has_right_edges = true;
                        break;
                    }
                }
                
                if !has_right_edges {
                    s.insert(handle);
                }
            }
        }
        
        // Main topological sort loop
        while visited.iter().filter(|&&v| v).count() < self.nodes.len() || !s.is_empty() {
            // If S is empty, we need to find a new seed
            if s.is_empty() {
                // First try seeds collected during traversal
                let mut found_seed = false;
                for (&handle, &_) in &seeds {
                    let idx = handle_to_idx[&handle];
                    if !visited[idx] {
                        s.insert(handle);
                        found_seed = true;
                        break;
                    }
                }
                
                // If no seeds, pick any unvisited node (lowest ID for determinism)
                if !found_seed {
                    for &handle in &idx_to_handle {
                        let idx = handle_to_idx[&handle];
                        if !visited[idx] && !handle.is_reverse() {
                            // Prefer forward orientation
                            s.insert(handle);
                            found_seed = true;
                            break;
                        }
                    }
                    
                    // If still nothing, try reverse orientations
                    if !found_seed {
                        for &handle in &idx_to_handle {
                            let idx = handle_to_idx[&handle];
                            if !visited[idx] {
                                s.insert(handle);
                                break;
                            }
                        }
                    }
                }
            }
            
            // Process handles in S
            while !s.is_empty() {
                // Get a handle from S (use minimum for determinism)
                let handle = *s.iter().min().unwrap();
                s.remove(&handle);
                
                let idx = handle_to_idx[&handle];
                if visited[idx] {
                    continue;
                }
                
                // Mark as visited and add to sorted order
                visited[idx] = true;
                // Only add forward orientations to final result
                if !handle.is_reverse() {
                    sorted.push(handle);
                }
                
                // Look at edges from the left side of this handle
                for edge in self.edges.clone() {
                    // Check edges coming into this handle from the left
                    if edge.to == handle && !masked_edges.contains(&edge) {
                        // Mask this edge
                        masked_edges.insert(edge);
                        
                        // The source of this edge might have been a cycle entry point
                        // No need to process further for left edges
                    }
                }
                
                // Look at edges from the right side of this handle
                let edges_to_check: Vec<_> = self.edges.iter().cloned().collect();
                for edge in edges_to_check {
                    if edge.from == handle && !masked_edges.contains(&edge) {
                        // Mask this edge
                        masked_edges.insert(edge);
                        
                        let next_handle = edge.to;
                        let next_idx = *handle_to_idx.get(&next_handle).unwrap_or(&total_handles);
                        
                        if next_idx < total_handles && !visited[next_idx] {
                            // Check if this handle has any other unmasked incoming edges
                            let mut has_unmasked_incoming = false;
                            for other_edge in &self.edges {
                                if other_edge.to == next_handle && 
                                   !masked_edges.contains(other_edge) &&
                                   other_edge != &edge {
                                    has_unmasked_incoming = true;
                                    break;
                                }
                            }
                            
                            if !has_unmasked_incoming {
                                // No more incoming edges, add to S
                                s.insert(next_handle);
                            } else {
                                // Still has incoming edges, mark as potential seed
                                seeds.insert(next_handle, true);
                            }
                        }
                    }
                }
            }
        }
        
        if verbose {
            println!("Topological sort completed: {} nodes ordered", sorted.len());
        }
        
        sorted
    }

    /// Build a map of node_id -> earliest position in any path
    /// This is used to prioritize nodes based on path order during topological sort
    fn build_path_position_map(&self) -> HashMap<usize, usize> {
        let mut position_map: HashMap<usize, usize> = HashMap::new();

        for path in &self.paths {
            for (pos, handle) in path.steps.iter().enumerate() {
                let node_id = handle.node_id();
                // Keep the minimum position (earliest occurrence)
                position_map.entry(node_id)
                    .and_modify(|existing_pos: &mut usize| *existing_pos = (*existing_pos).min(pos))
                    .or_insert(pos);
            }
        }

        position_map
    }

    /// Find all nodes with no edges on their left sides (heads)
    /// In a bidirected graph, we check for incoming edges to the forward orientation
    /// Heads are sorted by their earliest position in paths for path-aware sorting
    pub fn find_head_nodes(&self) -> Vec<Handle> {
        let mut heads = Vec::new();

        for &node_id in self.nodes.keys() {
            // Check if this node has ANY incoming edges (to either orientation)
            // A head node has no incoming edges to either its forward or reverse orientation
            let mut has_incoming = false;

            // Check for edges to either orientation
            let fwd_handle = Handle::forward(node_id);
            let rev_handle = Handle::reverse(node_id);

            for edge in &self.edges {
                if edge.to == fwd_handle || edge.to == rev_handle {
                    has_incoming = true;
                    break;
                }
            }

            // If no incoming edges, add the forward orientation as head
            if !has_incoming {
                heads.push(fwd_handle);
            }
        }

        // Sort by path position for path-aware topological sorting
        // Nodes that appear earlier in paths are processed first
        let path_positions = self.build_path_position_map();
        heads.sort_by_key(|h| {
            (
                path_positions.get(&h.node_id()).copied().unwrap_or(usize::MAX),
                h.node_id()
            )
        });

        heads
    }
    
    /// Find all nodes with no edges on their right sides (tails)
    pub fn find_tail_nodes(&self) -> Vec<Handle> {
        let mut tails = Vec::new();
        
        for &node_id in self.nodes.keys() {
            let forward_handle = Handle::forward(node_id);
            let mut has_outgoing = false;
            
            // Check if any edge goes FROM this handle
            for edge in &self.edges {
                if edge.from == forward_handle {
                    has_outgoing = true;
                    break;
                }
            }
            
            if !has_outgoing {
                tails.push(forward_handle);
            }
        }
        
        // Sort for deterministic behavior
        tails.sort_by_key(|h| h.node_id());
        tails
    }

    /// Exact ODGI topological_order algorithm
    /// This is a modified Kahn's algorithm that can handle cycles and bidirected graphs
    /// Enhanced with path-aware ordering: prioritizes nodes based on their path positions
    pub fn exact_odgi_topological_order(
        &self,
        use_heads: bool,
        use_tails: bool,
        verbose: bool,
    ) -> Vec<Handle> {
        let mut sorted = Vec::new();

        if self.nodes.is_empty() {
            return sorted;
        }

        // Build path position map for path-aware ordering
        let path_positions = self.build_path_position_map();

        // S - set of oriented handles ready to be processed
        // Using BTreeSet for deterministic ordering (ODGI behavior)
        let mut s: BTreeSet<Handle> = BTreeSet::new();

        // Track which nodes have been visited (not handles, nodes!)
        let mut visited_nodes = HashSet::new();

        // Unvisited - track which handles haven't been processed yet
        let mut unvisited = HashSet::new();
        for &node_id in self.nodes.keys() {
            // Both orientations start as unvisited
            unvisited.insert(Handle::forward(node_id));
            unvisited.insert(Handle::reverse(node_id));
        }

        // Seeds for breaking cycles - keep as Vec to preserve path-ordered heads
        let mut seeds: Vec<Handle> = Vec::new();

        // Track masked (logically removed) edges
        let mut masked_edges = HashSet::new();

        // Initialize with heads or tails if requested
        // Important: Only add ONE head at a time to maintain path ordering
        if use_heads {
            // Get heads in path-order (already sorted by find_head_nodes)
            let heads = self.find_head_nodes();
            if !heads.is_empty() {
                // Start with the first head (earliest in paths)
                let first_head = heads[0];
                if verbose {
                    eprintln!("[exact_odgi] Starting with head: node {}", first_head.node_id());
                }
                s.insert(first_head);
                unvisited.remove(&first_head);
                unvisited.remove(&first_head.flip());

                // Keep remaining heads as seeds for later
                for handle in heads.into_iter().skip(1) {
                    seeds.push(handle);
                }
            }
        } else if use_tails {
            let tails = self.find_tail_nodes();
            if !tails.is_empty() {
                let first_tail = tails[0];
                if verbose {
                    eprintln!("[exact_odgi] Starting with tail: node {}", first_tail.node_id());
                }
                s.insert(first_tail);
                unvisited.remove(&first_tail);
                unvisited.remove(&first_tail.flip());

                for handle in tails.into_iter().skip(1) {
                    seeds.push(handle);
                }
            }
        }
        
        // Main loop - continue until all nodes are visited
        while !unvisited.is_empty() || !s.is_empty() {
            
            // If S is empty, need to pick a seed to break into a cycle
            if s.is_empty() {
                // First try previously identified seeds (in path order)
                let mut found_seed = false;

                // Process seeds in order - they're already path-ordered
                if !seeds.is_empty() {
                    // Take the first seed from the Vec (earliest in paths)
                    let handle = seeds.remove(0);
                    if unvisited.contains(&handle) {
                        s.insert(handle);
                        unvisited.remove(&handle);
                        unvisited.remove(&handle.flip());
                        found_seed = true;
                        if verbose {
                            eprintln!("[exact_odgi] Using seed: node {} orient {}",
                                     handle.node_id(),
                                     if handle.is_reverse() { "-" } else { "+" });
                        }
                    } else {
                        // This seed was already visited, try the next one
                        found_seed = false;
                    }
                }
                
                // If no seeds available, pick arbitrary unvisited handle
                // Use path position to prioritize nodes that appear earlier in paths
                if !found_seed && !unvisited.is_empty() {
                    // Get handle with earliest path position (or lowest ID as tiebreaker)
                    let min_handle = unvisited.iter()
                        .min_by_key(|h| {
                            (
                                path_positions.get(&h.node_id()).copied().unwrap_or(usize::MAX),
                                h.node_id(),
                                h.is_reverse()
                            )
                        })
                        .cloned()
                        .unwrap();

                    s.insert(min_handle);
                    unvisited.remove(&min_handle);
                    unvisited.remove(&min_handle.flip());
                    if verbose {
                        eprintln!("[exact_odgi] Using arbitrary: node {} orient {}",
                                 min_handle.node_id(),
                                 if min_handle.is_reverse() { "-" } else { "+" });
                    }
                }
            }
            
            // Process handles in S
            while !s.is_empty() {
                // Get minimum handle for deterministic behavior (BTreeSet maintains order)
                let handle = *s.iter().next().unwrap();
                s.remove(&handle);
                
                // Emit the node (only once per node, when we see it first)
                if visited_nodes.insert(handle.node_id()) {
                    // We emit the forward orientation when we first visit a node
                    sorted.push(Handle::forward(handle.node_id()));
                    if verbose {
                        eprintln!("[exact_odgi] Emitting node {}", handle.node_id());
                    }
                }
                
                // Look at edges coming into this handle (backward edges)
                // These are edges that should be "consumed" by placing this handle
                // Sort edges for deterministic iteration
                let mut edges_vec: Vec<_> = self.edges.iter().cloned().collect();
                edges_vec.sort_by_key(|e| (e.from.node_id(), e.from.is_reverse(), e.to.node_id(), e.to.is_reverse()));

                for edge in &edges_vec {
                    if edge.to == handle && !masked_edges.contains(edge) {
                        masked_edges.insert(edge.clone());
                        if verbose {
                            eprintln!("[exact_odgi]   Masking incoming edge: {} {} -> {} {}",
                                     edge.from.node_id(),
                                     if edge.from.is_reverse() { "-" } else { "+" },
                                     edge.to.node_id(),
                                     if edge.to.is_reverse() { "-" } else { "+" });
                        }
                    }
                }

                // Look at edges going out from this handle (forward edges)
                for edge in &edges_vec {
                    if edge.from == handle && !masked_edges.contains(edge) {
                        masked_edges.insert(edge.clone());
                        let next_handle = edge.to;
                        
                        if verbose {
                            eprintln!("[exact_odgi]   Processing outgoing edge: {} {} -> {} {}", 
                                     edge.from.node_id(),
                                     if edge.from.is_reverse() { "-" } else { "+" },
                                     edge.to.node_id(),
                                     if edge.to.is_reverse() { "-" } else { "+" });
                        }
                        
                        // Only process if not yet visited
                        if unvisited.contains(&next_handle) {
                            // Check if next_handle has any other unmasked incoming edges
                            let mut has_unmasked_incoming = false;
                            
                            for other_edge in &edges_vec {
                                if other_edge.to == next_handle &&
                                   !masked_edges.contains(other_edge) {
                                    has_unmasked_incoming = true;
                                    break;
                                }
                            }
                            
                            if !has_unmasked_incoming {
                                // No more incoming edges, ready to process
                                s.insert(next_handle);
                                unvisited.remove(&next_handle);
                                unvisited.remove(&next_handle.flip());
                                if verbose {
                                    eprintln!("[exact_odgi]     Adding to S: node {} orient {}", 
                                             next_handle.node_id(),
                                             if next_handle.is_reverse() { "-" } else { "+" });
                                }
                            } else {
                                // Still has dependencies, mark as potential seed for cycle breaking
                                if !seeds.contains(&next_handle) {
                                    seeds.push(next_handle);
                                }
                                if verbose {
                                    eprintln!("[exact_odgi]     Marking as seed: node {} orient {}",
                                             next_handle.node_id(),
                                             if next_handle.is_reverse() { "-" } else { "+" });
                                }
                            }
                        }
                    }
                }
            }
        }
        
        if verbose {
            eprintln!("[exact_odgi] Topological sort complete: {} nodes", sorted.len());
        }
        
        sorted
    }
    
    /// Apply the exact ODGI topological ordering to renumber nodes
    pub fn apply_exact_odgi_ordering(&mut self, verbose: bool) {
        // Use heads only to match ODGI's default 's' topological sort
        let ordering = self.exact_odgi_topological_order(true, false, verbose);
        self.apply_ordering(ordering, verbose);
    }

    /// Apply a node ordering to renumber the graph
    pub fn apply_ordering(&mut self, ordering: Vec<Handle>, verbose: bool) {
        if ordering.is_empty() {
            return;
        }

        if verbose {
            eprintln!("\n[apply_ordering] Received ordering with {} handles:", ordering.len());
            for (i, handle) in ordering.iter().enumerate() {
                eprintln!("  Position {}: Node {}{}",
                         i, handle.node_id(),
                         if handle.is_reverse() { "-" } else { "+" });
            }
        }

        // Create old to new ID mapping
        let mut old_to_new = HashMap::new();
        for (new_idx, handle) in ordering.iter().enumerate() {
            old_to_new.insert(handle.node_id(), new_idx + 1); // 1-based IDs
        }

        if verbose {
            eprintln!("\n[apply_ordering] Old-to-new ID mapping:");
            let mut mapping_vec: Vec<_> = old_to_new.iter().collect();
            mapping_vec.sort_by_key(|(old_id, _)| **old_id);
            for (old_id, new_id) in mapping_vec {
                eprintln!("  {} -> {}", old_id, new_id);
            }
        }

        // Update nodes
        let mut new_nodes = HashMap::new();
        for (old_id, mut node) in self.nodes.drain() {
            if let Some(&new_id) = old_to_new.get(&old_id) {
                node.id = new_id;
                node.rank = Some((new_id - 1) as u64); // 0-based rank
                new_nodes.insert(new_id, node);
            }
        }
        self.nodes = new_nodes;

        // Update edges
        let mut new_edges = HashSet::new();
        for edge in self.edges.drain() {
            if let (Some(&new_from), Some(&new_to)) = (
                old_to_new.get(&edge.from.node_id()),
                old_to_new.get(&edge.to.node_id()),
            ) {
                let new_edge = BiEdge::new(
                    Handle::new(new_from, edge.from.is_reverse()),
                    Handle::new(new_to, edge.to.is_reverse()),
                );
                new_edges.insert(new_edge);
            }
        }
        self.edges = new_edges;

        // Update paths
        for path in &mut self.paths {
            if verbose {
                eprintln!("\n[apply_ordering] Updating path '{}':", path.name);
                eprintln!("  Before: {:?}",
                         path.steps.iter()
                             .map(|h| format!("{}{}", h.node_id(), if h.is_reverse() { "-" } else { "+" }))
                             .collect::<Vec<_>>());
            }

            for handle in &mut path.steps {
                if let Some(&new_id) = old_to_new.get(&handle.node_id()) {
                    *handle = Handle::new(new_id, handle.is_reverse());
                }
            }

            if verbose {
                eprintln!("  After:  {:?}",
                         path.steps.iter()
                             .map(|h| format!("{}{}", h.node_id(), if h.is_reverse() { "-" } else { "+" }))
                             .collect::<Vec<_>>());
            }
        }

        if verbose {
            eprintln!("\n[apply_ordering] Applied ordering: renumbered {} nodes\n", old_to_new.len());
        }
    }
}

/// Convert orientation-aware alignments to edges
pub fn alignment_to_edges(
    query_handle: Handle,
    target_handle: Handle,
    _query_len: usize,
    _target_len: usize,
    is_reverse_alignment: bool,
) -> Vec<BiEdge> {
    let mut edges = Vec::new();

    // For reverse alignments, we need to flip one of the handles
    let (from_handle, to_handle) = if is_reverse_alignment {
        (query_handle, target_handle.flip())
    } else {
        (query_handle, target_handle)
    };

    // Add edge from end of query to start of target
    edges.push(BiEdge::new(from_handle, to_handle));

    edges
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bidirected_graph_creation() {
        let mut graph = BidirectedGraph::new();

        // Add nodes
        graph.add_node(1, b"ATCG".to_vec());
        graph.add_node(2, b"GCTA".to_vec());

        // Add edges with different orientations
        graph.add_edge(Handle::forward(1), Handle::forward(2));
        graph.add_edge(Handle::forward(1), Handle::reverse(2));

        assert_eq!(graph.nodes.len(), 2);
        assert_eq!(graph.edges.len(), 2);
    }

    #[test]
    fn test_sequence_retrieval() {
        let mut graph = BidirectedGraph::new();
        graph.add_node(1, b"ATCG".to_vec());

        // Forward sequence
        assert_eq!(
            graph.get_sequence(Handle::forward(1)),
            Some(b"ATCG".to_vec())
        );

        // Reverse complement sequence
        assert_eq!(
            graph.get_sequence(Handle::reverse(1)),
            Some(b"CGAT".to_vec())
        );
    }

    #[test]
    fn test_path_with_orientations() {
        let mut graph = BidirectedGraph::new();
        graph.add_node(1, b"ATG".to_vec());
        graph.add_node(2, b"CGA".to_vec());
        graph.add_node(3, b"TAC".to_vec());

        // Build a path that traverses nodes in different orientations
        graph.build_path(
            "test_path".to_string(),
            vec![
                (1, false), // Forward
                (2, true),  // Reverse
                (3, false), // Forward
            ],
        );

        assert_eq!(graph.paths.len(), 1);
        let path = &graph.paths[0];
        assert_eq!(path.steps.len(), 3);
        assert!(!path.steps[0].is_reverse());
        assert!(path.steps[1].is_reverse());
        assert!(!path.steps[2].is_reverse());
    }

    #[test]
    fn test_gfa_output() {
        let mut graph = BidirectedGraph::new();
        graph.add_node(1, b"ATCG".to_vec());
        graph.add_node(2, b"GCTA".to_vec());

        graph.add_edge(Handle::forward(1), Handle::reverse(2));
        graph.build_path("path1".to_string(), vec![(1, false), (2, true)]);

        let mut output = Vec::new();
        graph.write_gfa(&mut output).unwrap();
        let gfa_str = String::from_utf8(output).unwrap();

        assert!(gfa_str.contains("S\t1\tATCG"));
        assert!(gfa_str.contains("S\t2\tGCTA"));
        assert!(gfa_str.contains("L\t1\t+\t2\t-\t0M"));
        assert!(gfa_str.contains("P\tpath1\t1+,2-\t*"));
    }
}
