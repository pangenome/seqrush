use crate::bidirected_graph::{BiEdge, Handle};
use crate::bidirected_ops::BidirectedGraph;
use std::collections::{HashMap, HashSet};

impl BidirectedGraph {
    /// Perform compaction using the perfect neighbor algorithm
    /// Returns the number of nodes that were compacted
    pub fn compact_perfect_neighbors(&mut self) -> Result<usize, String> {
        // Step 1: Validate paths before compaction
        let original_path_sequences = self.validate_and_save_paths()?;

        // Track unique original node IDs that are compacted away
        let original_nodes: HashSet<usize> = self.nodes.iter()
            .enumerate()
            .filter_map(|(id, n)| if n.is_some() { Some(id) } else { None })
            .collect();
        let mut compacted_original_nodes: HashSet<usize> = HashSet::new();
        let mut iteration = 0;

        // Iteratively compact perfect neighbor pairs
        loop {
            iteration += 1;
            eprintln!("Compaction iteration {}", iteration);

            // Find perfect neighbor pairs
            let perfect_pairs = self.find_perfect_neighbor_pairs();

            if perfect_pairs.is_empty() {
                eprintln!("No more perfect pairs found, stopping compaction");
                break;
            }

            eprintln!("Found {} perfect neighbor pairs", perfect_pairs.len());

            // Apply compaction for these pairs only (no transitive grouping)
            let nodes_compacted = self.apply_direct_compaction(&perfect_pairs)?;
            // Track which original nodes disappeared after this iteration
            for (from, to) in &perfect_pairs {
                let f = from.node_id();
                let t = to.node_id();
                if original_nodes.contains(&f) {
                    compacted_original_nodes.insert(f);
                }
                if original_nodes.contains(&t) {
                    compacted_original_nodes.insert(t);
                }
            }
            let _ = nodes_compacted; // progress indicator still used for logging

            eprintln!("Compacted {} nodes in this iteration", nodes_compacted);

            // If we found pairs but didn't compact any, we're stuck
            if !perfect_pairs.is_empty() && nodes_compacted == 0 {
                eprintln!(
                    "ERROR: Found pairs but made no progress, stopping to prevent infinite loop"
                );
                break;
            }
        }

        // Validate paths after compaction
        self.validate_paths_match(&original_path_sequences)?;

        Ok(compacted_original_nodes.len())
    }

    /// Save current path sequences for validation
    fn validate_and_save_paths(&self) -> Result<HashMap<String, Vec<u8>>, String> {
        let mut path_sequences = HashMap::new();

        for path in &self.paths {
            let mut sequence = Vec::new();

            for &handle in &path.steps {
                let node_id = handle.node_id();
                let node = self.nodes.get(node_id)
                    .and_then(|n| n.as_ref())
                    .ok_or_else(|| {
                        format!(
                            "Path {} references non-existent node {}",
                            path.name, node_id
                        )
                    })?;

                if handle.is_reverse() {
                    // Reverse complement for reverse orientation
                    let rc = crate::bidirected_graph::reverse_complement(&node.sequence);
                    sequence.extend_from_slice(&rc);
                } else {
                    sequence.extend_from_slice(&node.sequence);
                }
            }

            path_sequences.insert(path.name.clone(), sequence);
        }

        Ok(path_sequences)
    }

    /// Find all perfect neighbor pairs in the graph
    fn find_perfect_neighbor_pairs(&self) -> Vec<(Handle, Handle)> {
        let mut perfect_pairs = Vec::new();

        // Build connections from path traversals - the paths ARE the links!
        let (inbound, outbound) = self.build_connections_from_paths();

        // Check each handle for perfect neighbors
        for (node_id, node_opt) in self.nodes.iter().enumerate() {
            if node_opt.is_none() {
                continue;
            }
            for is_reverse in [false, true] {
                let handle = Handle::new(node_id, is_reverse);

                // Get outbound connections from this handle
                let outgoing = outbound
                    .get(&handle)
                    .map(|set| set.iter().cloned().collect::<Vec<_>>())
                    .unwrap_or_default();

                // A has only one outbound connection
                if outgoing.len() == 1 {
                    let next_handle = outgoing[0];

                    // Get inbound connections to the next handle
                    let next_incoming = inbound
                        .get(&next_handle)
                        .map(|set| set.iter().cloned().collect::<Vec<_>>())
                        .unwrap_or_default();

                    // B has only one inbound connection, and it's from A
                    if next_incoming.len() == 1 && next_incoming[0] == handle {
                        // That's it! A→B is a perfect pair
                        // If A only goes to B and B only comes from A,
                        // they MUST be adjacent in all paths
                        perfect_pairs.push((handle, next_handle));
                    }
                }
            }
        }

        perfect_pairs
    }

    /// Build connection maps from path traversals
    fn build_connections_from_paths(
        &self,
    ) -> (
        HashMap<Handle, HashSet<Handle>>,
        HashMap<Handle, HashSet<Handle>>,
    ) {
        let mut inbound: HashMap<Handle, HashSet<Handle>> = HashMap::new();
        let mut outbound: HashMap<Handle, HashSet<Handle>> = HashMap::new();

        // The paths define the links!
        for path in &self.paths {
            for i in 0..path.steps.len() - 1 {
                let from = path.steps[i];
                let to = path.steps[i + 1];

                outbound.entry(from).or_default().insert(to);
                inbound.entry(to).or_default().insert(from);
            }
        }

        (inbound, outbound)
    }

    /// Use union-find to group nodes that should be compacted together
    fn build_compaction_groups(&self, perfect_pairs: &[(Handle, Handle)]) -> Vec<Vec<Handle>> {
        // Create a mapping from handles to indices
        let mut handle_to_idx = HashMap::new();
        let mut idx_to_handle = Vec::new();
        let mut idx = 0;

        for &(from, to) in perfect_pairs {
            if !handle_to_idx.contains_key(&from) {
                handle_to_idx.insert(from, idx);
                idx_to_handle.push(from);
                idx += 1;
            }
            if !handle_to_idx.contains_key(&to) {
                handle_to_idx.insert(to, idx);
                idx_to_handle.push(to);
                idx += 1;
            }
        }

        // Build union-find manually
        let mut parent: Vec<usize> = (0..idx).collect();
        let mut rank: Vec<usize> = vec![0; idx];

        fn find(parent: &mut Vec<usize>, x: usize) -> usize {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }

        fn unite(parent: &mut Vec<usize>, rank: &mut Vec<usize>, x: usize, y: usize) {
            let root_x = find(parent, x);
            let root_y = find(parent, y);

            if root_x != root_y {
                if rank[root_x] < rank[root_y] {
                    parent[root_x] = root_y;
                } else if rank[root_x] > rank[root_y] {
                    parent[root_y] = root_x;
                } else {
                    parent[root_y] = root_x;
                    rank[root_x] += 1;
                }
            }
        }

        for &(from, to) in perfect_pairs {
            let from_idx = handle_to_idx[&from];
            let to_idx = handle_to_idx[&to];
            unite(&mut parent, &mut rank, from_idx, to_idx);
        }

        // Extract groups
        let mut groups: HashMap<usize, Vec<Handle>> = HashMap::new();
        for i in 0..idx {
            let root = find(&mut parent, i);
            groups.entry(root).or_default().push(idx_to_handle[i]);
        }

        groups.into_values().collect()
    }

    /// Order nodes within each compaction group based on path traversal
    fn order_compaction_groups(&self, groups: &[Vec<Handle>]) -> Result<Vec<Vec<Handle>>, String> {
        let mut ordered_groups = Vec::new();

        for group in groups {
            if group.len() < 2 {
                continue; // Skip single-node groups
            }

            // Find the ordering by looking at any path that traverses this group
            let ordered = self.order_single_group(group)?;
            ordered_groups.push(ordered);
        }

        Ok(ordered_groups)
    }

    /// Order a single compaction group using path information
    fn order_single_group(&self, group: &[Handle]) -> Result<Vec<Handle>, String> {
        let group_set: HashSet<_> = group.iter().cloned().collect();

        // Find a path that contains nodes from this group
        for path in &self.paths {
            let mut group_nodes_in_path = Vec::new();

            for &handle in &path.steps {
                if group_set.contains(&handle) {
                    group_nodes_in_path.push(handle);
                }
            }

            if group_nodes_in_path.len() == group.len() {
                // This path contains all nodes in the group, use its ordering
                return Ok(group_nodes_in_path);
            }
        }

        // If no path contains all nodes in order, try to find any valid ordering
        // This can happen in complex graphs with cycles
        // For now, just use the order from the group
        // In a more sophisticated implementation, we'd analyze the graph structure
        Ok(group.to_vec())
    }

    /// Apply direct compaction for perfect neighbor pairs without transitive grouping
    fn apply_direct_compaction(
        &mut self,
        perfect_pairs: &[(Handle, Handle)],
    ) -> Result<usize, String> {
        let mut nodes_compacted = 0;
        let mut handle_remapping: HashMap<Handle, Handle> = HashMap::new();
        let mut new_node_id = self.nodes.len();
        let mut actually_merged_pairs: Vec<(Handle, Handle)> = Vec::new();

        // First, let's see what pairs we're trying to compact
        if perfect_pairs.len() < 10 {
            eprintln!("Perfect pairs to compact:");
            for &(from, to) in perfect_pairs {
                eprintln!("  {:?} -> {:?}", from, to);
            }
        }

        // Process each pair individually
        for &(from, to) in perfect_pairs {
            // Skip if either handle was already remapped
            if handle_remapping.contains_key(&from) || handle_remapping.contains_key(&to) {
                if perfect_pairs.len() < 10 {
                    eprintln!("  Skipping {:?} -> {:?} (already remapped)", from, to);
                }
                continue;
            }

            // Get the nodes
            let from_node = match self.nodes.get(from.node_id()).and_then(|n| n.as_ref()) {
                Some(node) => node,
                None => {
                    if perfect_pairs.len() < 10 {
                        eprintln!(
                            "  Skipping {:?} -> {:?} (from node {} doesn't exist)",
                            from,
                            to,
                            from.node_id()
                        );
                    }
                    continue; // Node was already removed
                }
            };
            let to_node = match self.nodes.get(to.node_id()).and_then(|n| n.as_ref()) {
                Some(node) => node,
                None => {
                    if perfect_pairs.len() < 10 {
                        eprintln!(
                            "  Skipping {:?} -> {:?} (to node {} doesn't exist)",
                            from,
                            to,
                            to.node_id()
                        );
                    }
                    continue; // Node was already removed
                }
            };

            // Create merged sequence
            let mut merged_sequence = Vec::new();

            // Add sequence from 'from' handle
            if from.is_reverse() {
                merged_sequence.extend_from_slice(&crate::bidirected_graph::reverse_complement(
                    &from_node.sequence,
                ));
            } else {
                merged_sequence.extend_from_slice(&from_node.sequence);
            }

            // Add sequence from 'to' handle
            if to.is_reverse() {
                merged_sequence.extend_from_slice(&crate::bidirected_graph::reverse_complement(
                    &to_node.sequence,
                ));
            } else {
                merged_sequence.extend_from_slice(&to_node.sequence);
            }

            // Create new node
            let new_handle = Handle::forward(new_node_id);
            self.add_node(new_node_id, merged_sequence);

            // Map the handles used in this merge
            handle_remapping.insert(from, new_handle);
            handle_remapping.insert(to, new_handle);

            // Also map their complements appropriately
            // The complement of 'from' maps to the complement of the new handle
            handle_remapping.insert(from.flip(), new_handle.flip());
            // The complement of 'to' maps to the complement of the new handle
            handle_remapping.insert(to.flip(), new_handle.flip());

            new_node_id += 1;
            nodes_compacted += 2; // We merged 2 nodes
        }

        // Store which pairs were actually merged
        let merged_pairs: HashSet<(Handle, Handle)> = perfect_pairs
            .iter()
            .filter(|(from, to)| {
                handle_remapping.contains_key(from) && handle_remapping.contains_key(to)
            })
            .cloned()
            .collect();

        // Update paths - only replace exact A→B sequences with C
        for path in &mut self.paths {
            let mut new_steps = Vec::new();
            let mut i = 0;

            while i < path.steps.len() {
                if i + 1 < path.steps.len() {
                    let current = path.steps[i];
                    let next = path.steps[i + 1];

                    // Check if this is a merged pair
                    if merged_pairs.contains(&(current, next)) {
                        // This is a perfect pair that was merged
                        if let Some(&new_handle) = handle_remapping.get(&current) {
                            new_steps.push(new_handle);
                            i += 2; // Skip both handles
                            continue;
                        }
                    }
                }

                // Not part of a merged pair, keep the original handle
                new_steps.push(path.steps[i]);
                i += 1;
            }

            path.steps = new_steps;
        }

        // Recompute which nodes are still used by any path
        let mut used_node_ids: HashSet<usize> = HashSet::new();
        for path in &self.paths {
            for &h in &path.steps {
                used_node_ids.insert(h.node_id());
            }
        }
        // Remove only nodes that are not used anymore
        let nodes_to_remove: HashSet<_> = handle_remapping
            .keys()
            .map(|h| h.node_id())
            .filter(|nid| !used_node_ids.contains(nid))
            .collect();
        for node_id in nodes_to_remove {
            if node_id < self.nodes.len() {
                self.nodes[node_id] = None;
            }
        }

        // Rebuild edges from updated paths to reflect current structure
        let mut rebuilt_edges: HashSet<BiEdge> = HashSet::new();
        for path in &self.paths {
            for w in path.steps.windows(2) {
                if let [from, to] = w {
                    if from != to {
                        rebuilt_edges.insert(BiEdge {
                            from: *from,
                            to: *to,
                        });
                    }
                }
            }
        }
        self.edges = rebuilt_edges;

        Ok(nodes_compacted)
    }

    /// Apply compaction by creating new merged nodes and updating the graph
    fn apply_compaction(&mut self, ordered_groups: &[Vec<Handle>]) -> Result<usize, String> {
        let mut nodes_compacted = 0;
        let mut handle_remapping: HashMap<Handle, Handle> = HashMap::new();
        let mut new_node_id = self.nodes.len();

        for group in ordered_groups {
            if group.len() < 2 {
                continue;
            }

            // Create merged sequence
            let mut merged_sequence = Vec::new();
            for &handle in group {
                let node = self.nodes.get(handle.node_id())
                    .and_then(|n| n.as_ref())
                    .ok_or_else(|| format!("Node {} not found", handle.node_id()))?;
                if handle.is_reverse() {
                    merged_sequence.extend_from_slice(
                        &crate::bidirected_graph::reverse_complement(&node.sequence),
                    );
                } else {
                    merged_sequence.extend_from_slice(&node.sequence);
                }
            }

            // Debug: log compaction groups that might cause issues
            if new_node_id > 4100 && new_node_id < 4200 {
                eprintln!(
                    "Creating merged node {} ({} bp) from {} nodes:",
                    new_node_id,
                    merged_sequence.len(),
                    group.len()
                );
                for &handle in group {
                    if let Some(node) = self.nodes.get(handle.node_id()).and_then(|n| n.as_ref()) {
                        eprintln!(
                            "  - node {} ({} bp, {})",
                            handle.node_id(),
                            node.sequence.len(),
                            if handle.is_reverse() { "rev" } else { "fwd" }
                        );
                    }
                }
            }

            // Create new node
            let new_handle = Handle::forward(new_node_id);
            self.add_node(new_node_id, merged_sequence);

            // Map all handles in the group to the new handle
            // Also map their reverse complements
            for &handle in group {
                // Map the forward orientation
                handle_remapping.insert(handle, new_handle);
                // Map the reverse orientation of the same node
                let reverse_handle = Handle::new(handle.node_id(), !handle.is_reverse());
                handle_remapping.insert(
                    reverse_handle,
                    Handle::new(new_handle.node_id(), !new_handle.is_reverse()),
                );
            }

            new_node_id += 1;
            nodes_compacted += group.len();
        }

        // Remove old nodes
        let nodes_to_remove: HashSet<_> = handle_remapping.keys().map(|h| h.node_id()).collect();

        for node_id in nodes_to_remove {
            if node_id < self.nodes.len() {
                self.nodes[node_id] = None;
            }
        }

        // Update edges
        let mut new_edges = HashSet::new();
        for edge in &self.edges {
            let from = handle_remapping
                .get(&edge.from)
                .copied()
                .unwrap_or(edge.from);
            let to = handle_remapping.get(&edge.to).copied().unwrap_or(edge.to);

            // Don't add self-loops from compaction
            if from != to {
                new_edges.insert(BiEdge { from, to });
            }
        }
        self.edges = new_edges;

        // Update paths
        for path in &mut self.paths {
            let mut new_steps = Vec::new();
            let mut i = 0;

            while i < path.steps.len() {
                let handle = path.steps[i];

                if let Some(&new_handle) = handle_remapping.get(&handle) {
                    // This handle was compacted, add the new handle
                    new_steps.push(new_handle);

                    // Skip all subsequent handles that were compacted into the same node
                    i += 1;
                    while i < path.steps.len()
                        && handle_remapping.get(&path.steps[i]) == Some(&new_handle)
                    {
                        i += 1;
                    }
                } else {
                    // This handle wasn't compacted
                    new_steps.push(handle);
                    i += 1;
                }
            }

            path.steps = new_steps;
        }

        Ok(nodes_compacted)
    }

    /// Validate that paths still produce the same sequences after compaction
    fn validate_paths_match(
        &self,
        original_sequences: &HashMap<String, Vec<u8>>,
    ) -> Result<(), String> {
        for path in &self.paths {
            let mut current_sequence = Vec::new();

            for &handle in &path.steps {
                let node = self.nodes.get(handle.node_id())
                    .and_then(|n| n.as_ref())
                    .ok_or_else(|| {
                        format!(
                            "Path {} references non-existent node {}",
                            path.name,
                            handle.node_id()
                        )
                    })?;

                if handle.is_reverse() {
                    current_sequence.extend_from_slice(
                        &crate::bidirected_graph::reverse_complement(&node.sequence),
                    );
                } else {
                    current_sequence.extend_from_slice(&node.sequence);
                }
            }

            let original = original_sequences
                .get(&path.name)
                .ok_or_else(|| format!("Path {} not found in original sequences", path.name))?;

            if current_sequence != *original {
                // Debug: show what changed
                eprintln!("Path {} validation failed:", path.name);
                eprintln!("  Original length: {}", original.len());
                eprintln!("  Current length:  {}", current_sequence.len());
                eprintln!("  Path has {} steps", path.steps.len());

                // Show first few steps and last few steps
                eprintln!("  First 10 steps:");
                for (i, &handle) in path.steps.iter().enumerate().take(10) {
                    if let Some(node) = self.nodes.get(handle.node_id()).and_then(|n| n.as_ref()) {
                        eprintln!(
                            "    Step {}: node {} ({} bp, {})",
                            i,
                            handle.node_id(),
                            node.sequence.len(),
                            if handle.is_reverse() { "rev" } else { "fwd" }
                        );
                    }
                }

                eprintln!("  Last 10 steps:");
                let start = path.steps.len().saturating_sub(10);
                for (i, &handle) in path.steps.iter().enumerate().skip(start) {
                    if let Some(node) = self.nodes.get(handle.node_id()).and_then(|n| n.as_ref()) {
                        eprintln!(
                            "    Step {}: node {} ({} bp, {})",
                            i,
                            handle.node_id(),
                            node.sequence.len(),
                            if handle.is_reverse() { "rev" } else { "fwd" }
                        );
                    }
                }

                // Check if sequences differ at start or end
                let min_len = original.len().min(current_sequence.len());
                let mut first_diff = None;
                for i in 0..min_len {
                    if original[i] != current_sequence[i] {
                        first_diff = Some(i);
                        break;
                    }
                }

                if let Some(pos) = first_diff {
                    eprintln!("  First difference at position {}", pos);
                    eprintln!(
                        "    Original: {}",
                        String::from_utf8_lossy(
                            &original
                                [pos.saturating_sub(5)..pos.saturating_add(5).min(original.len())]
                        )
                    );
                    eprintln!(
                        "    Current:  {}",
                        String::from_utf8_lossy(
                            &current_sequence[pos.saturating_sub(5)
                                ..pos.saturating_add(5).min(current_sequence.len())]
                        )
                    );
                } else if original.len() != current_sequence.len() {
                    eprintln!(
                        "  Sequences match up to position {} but have different lengths",
                        min_len
                    );
                    if current_sequence.len() > original.len() {
                        eprintln!(
                            "  Extra bases: {}",
                            String::from_utf8_lossy(&current_sequence[original.len()..])
                        );
                    }
                }

                return Err(format!(
                    "Path {} sequence changed after compaction! Original length: {}, Current length: {}", 
                    path.name, original.len(), current_sequence.len()
                ));
            }
        }

        Ok(())
    }
}
