use crate::bidirected_graph::{BiEdge, BiNode, BiPath, Handle};
use crate::graph_ops::{Edge, Graph, Node};
use bitvec::prelude::*;
use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};

/// A bidirected graph that extends the basic Graph with orientation support
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
    /// Compact the graph by merging linear chains of nodes
    pub fn compact(&mut self) {
        // Find simple components (linear chains)
        let components = self.find_simple_components();

        // Merge each component
        for component in components {
            if component.len() >= 2 {
                self.merge_component(&component);
            }
        }
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
            // Check all paths
            for path in &self.paths {
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

    /// Merge a component (chain of handles) into a single node
    fn merge_component(&mut self, handles: &[Handle]) {
        if handles.len() < 2 {
            return;
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

        // Update paths
        for path in &mut self.paths {
            let mut new_steps = Vec::new();
            let mut i = 0;

            while i < path.steps.len() {
                let current_handle = path.steps[i];

                // Check if this handle is in our chain
                if let Some(&new_handle) = handle_map.get(&current_handle) {
                    // This handle is being replaced
                    // Check if we're at the start of the complete chain
                    let mut is_chain_start = true;
                    if i + handles.len() <= path.steps.len() {
                        for j in 0..handles.len() {
                            if path.steps[i + j] != handles[j] {
                                is_chain_start = false;
                                break;
                            }
                        }
                    } else {
                        is_chain_start = false;
                    }

                    if is_chain_start {
                        // We found the complete chain, replace it with the new handle
                        new_steps.push(new_handle);
                        i += handles.len();
                    } else {
                        // This handle is part of the chain but we're not at the chain start
                        // This shouldn't happen if paths are coherent
                        eprintln!(
                            "WARNING: Found handle {} in chain but not at chain start in path {}",
                            current_handle, path.name
                        );
                        new_steps.push(current_handle);
                        i += 1;
                    }
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
        self.edges.insert(BiEdge::new(from, to));
    }

    /// Get sequence for a handle (forward or reverse complement)
    pub fn get_sequence(&self, handle: Handle) -> Option<Vec<u8>> {
        self.nodes
            .get(&handle.node_id())
            .map(|node| node.get_sequence(handle.is_reverse()))
    }

    /// Check if an edge exists
    pub fn has_edge(&self, from: Handle, to: Handle) -> bool {
        self.edges.contains(&BiEdge::new(from, to))
    }

    /// Get all edges from a handle
    pub fn edges_from(&self, handle: Handle) -> Vec<Handle> {
        self.edges
            .iter()
            .filter(|edge| edge.from == handle)
            .map(|edge| edge.to)
            .collect()
    }

    /// Get all edges to a handle
    pub fn edges_to(&self, handle: Handle) -> Vec<Handle> {
        self.edges
            .iter()
            .filter(|edge| edge.to == handle)
            .map(|edge| edge.from)
            .collect()
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

    /// Perform bidirected topological sort using modified Kahn's algorithm
    /// This handles cycles and bidirected edges properly
    pub fn topological_sort(
        &mut self,
        use_heads: bool,
        use_tails: bool,
        verbose: bool,
    ) -> Vec<Handle> {
        let mut ordering = Vec::new();
        let num_nodes = self.nodes.len();

        if num_nodes == 0 {
            return ordering;
        }

        // Create handle set for all possible orientations
        let mut all_handles = BTreeSet::new();
        for &node_id in self.nodes.keys() {
            all_handles.insert(Handle::forward(node_id));
            all_handles.insert(Handle::reverse(node_id));
        }

        // Track visited handles and potential seeds for cycle breaking
        let mut visited: BitVec<u64, Lsb0> = BitVec::repeat(false, all_handles.len());
        let mut seeds: BitVec<u64, Lsb0> = BitVec::repeat(false, all_handles.len());

        // Create handle to index mapping
        let mut handle_to_idx = HashMap::new();
        for (idx, &handle) in all_handles.iter().enumerate() {
            handle_to_idx.insert(handle, idx);
        }

        // Build edge index for efficient queries
        let mut edges_by_from: HashMap<Handle, Vec<(Handle, usize)>> = HashMap::new();
        let mut edges_by_to: HashMap<Handle, Vec<(Handle, usize)>> = HashMap::new();
        let mut edge_list = Vec::new();

        for (edge_idx, edge) in self.edges.iter().enumerate() {
            edge_list.push(*edge);

            // Direct edge: from → to
            edges_by_from
                .entry(edge.from)
                .or_default()
                .push((edge.to, edge_idx));
            edges_by_to
                .entry(edge.to)
                .or_default()
                .push((edge.from, edge_idx));

            // Implied reverse edge: to.flip() → from.flip()
            edges_by_from
                .entry(edge.to.flip())
                .or_default()
                .push((edge.from.flip(), edge_idx));
            edges_by_to
                .entry(edge.from.flip())
                .or_default()
                .push((edge.to.flip(), edge_idx));
        }

        // Track masked edges (simulating edge removal without modifying graph)
        let mut masked_edges: BitVec<u64, Lsb0> = BitVec::repeat(false, edge_list.len());

        // Count unmasked in-edges for each handle
        let mut in_degree: HashMap<Handle, usize> = HashMap::new();
        for &handle in &all_handles {
            in_degree.insert(handle, 0);
        }

        for edge in &self.edges {
            // Direct edge: from → to
            *in_degree.get_mut(&edge.to).unwrap() += 1;
            // Implied reverse edge: to.flip() → from.flip()
            *in_degree.get_mut(&edge.from.flip()).unwrap() += 1;
        }

        // Initialize queue with appropriate seeds
        let mut queue = VecDeque::new();

        // Function to add seeds based on criteria
        let add_seeds = |queue: &mut VecDeque<Handle>,
                         in_degree: &HashMap<Handle, usize>,
                         visited: &BitVec<u64, Lsb0>,
                         handle_to_idx: &HashMap<Handle, usize>,
                         use_heads: bool,
                         use_tails: bool| {
            // Add head nodes (no incoming edges)
            if use_heads {
                for &handle in &all_handles {
                    let idx = handle_to_idx[&handle];
                    if !visited[idx] && in_degree[&handle] == 0 {
                        queue.push_back(handle);
                    }
                }
            }

            // Add tail nodes (no outgoing edges)
            if use_tails && queue.is_empty() {
                for &handle in &all_handles {
                    let idx = handle_to_idx[&handle];
                    if !visited[idx]
                        && edges_by_from
                            .get(&handle)
                            .map(|v| v.is_empty())
                            .unwrap_or(true)
                    {
                        queue.push_back(handle);
                    }
                }
            }
        };

        // Start with initial seeds - but prioritize forward orientations
        add_seeds(
            &mut queue,
            &in_degree,
            &visited,
            &handle_to_idx,
            use_heads,
            use_tails,
        );

        // Sort queue to prioritize forward orientations
        let mut queue_vec: Vec<_> = queue.drain(..).collect();
        queue_vec.sort_by_key(|h| (h.is_reverse(), h.node_id()));
        queue.extend(queue_vec);

        // Main topological sort loop
        while ordering.len() < all_handles.len() {
            // If queue is empty, we need to break a cycle
            if queue.is_empty() {
                // Try seeds collected during traversal (good cycle entry points)
                let mut found_seed = false;
                for &handle in &all_handles {
                    let idx = handle_to_idx[&handle];
                    if seeds[idx] && !visited[idx] {
                        queue.push_back(handle);
                        seeds.set(idx, false);
                        found_seed = true;
                        break;
                    }
                }

                // If no seeds, pick any unvisited node
                if !found_seed {
                    for &handle in &all_handles {
                        let idx = handle_to_idx[&handle];
                        if !visited[idx] {
                            queue.push_back(handle);
                            if verbose {
                                println!("Breaking cycle at handle {}", handle);
                            }
                            break;
                        }
                    }
                }

                // If still empty, we're done
                if queue.is_empty() {
                    break;
                }
            }

            // Process next handle
            let handle = queue.pop_front().unwrap();
            let idx = handle_to_idx[&handle];

            // Skip if already visited
            if visited[idx] {
                continue;
            }

            // Mark as visited and add to ordering
            visited.set(idx, true);
            ordering.push(handle);

            // Process outgoing edges
            if let Some(out_edges) = edges_by_from.get(&handle) {
                for &(to_handle, edge_idx) in out_edges {
                    // Skip if edge is masked
                    if masked_edges[edge_idx] {
                        continue;
                    }

                    // Mask this edge
                    masked_edges.set(edge_idx, true);

                    // Update in-degree of target
                    if let Some(degree) = in_degree.get_mut(&to_handle) {
                        if *degree > 0 {
                            *degree -= 1;

                            // If in-degree becomes 0, add to queue
                            if *degree == 0 {
                                let to_idx = handle_to_idx[&to_handle];
                                if !visited[to_idx] {
                                    queue.push_back(to_handle);
                                }
                            }
                        }
                    }

                    // If target still has incoming edges but we've seen it,
                    // mark it as a potential seed for cycle breaking
                    let to_idx = handle_to_idx[&to_handle];
                    if !visited[to_idx] && in_degree[&to_handle] > 0 {
                        seeds.set(to_idx, true);
                    }
                }
            }
        }

        // Filter ordering to only include forward orientations of nodes that exist
        let mut node_ordering = Vec::new();
        let mut seen_nodes = HashSet::new();

        for handle in ordering {
            let node_id = handle.node_id();
            if self.nodes.contains_key(&node_id) && seen_nodes.insert(node_id) {
                node_ordering.push(Handle::forward(node_id));
            }
        }

        if verbose {
            println!(
                "Topological sort completed: {} nodes ordered",
                node_ordering.len()
            );
        }

        node_ordering
    }

    /// Apply a node ordering to renumber the graph
    pub fn apply_ordering(&mut self, ordering: Vec<Handle>, verbose: bool) {
        if ordering.is_empty() {
            return;
        }

        // Create old to new ID mapping
        let mut old_to_new = HashMap::new();
        for (new_idx, handle) in ordering.iter().enumerate() {
            old_to_new.insert(handle.node_id(), new_idx + 1); // 1-based IDs
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
            for handle in &mut path.steps {
                if let Some(&new_id) = old_to_new.get(&handle.node_id()) {
                    *handle = Handle::new(new_id, handle.is_reverse());
                }
            }
        }

        if verbose {
            println!("Applied ordering: renumbered {} nodes", old_to_new.len());
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
