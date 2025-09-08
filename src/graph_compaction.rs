use crate::graph_ops::{Edge, Graph};
use std::collections::{HashMap, HashSet};

impl Graph {
    /// Find nodes that can be merged (simple components)
    /// Following ODGI's simple_components approach
    pub fn find_simple_components(&self) -> Vec<Vec<usize>> {
        let debug = std::env::var("SEQRUSH_DEBUG_COMPACT").is_ok();
        let mut components = Vec::new();

        // Build adjacency lists
        let mut forward_edges: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut backward_edges: HashMap<usize, Vec<usize>> = HashMap::new();

        for edge in &self.edges {
            forward_edges.entry(edge.from).or_default().push(edge.to);
            backward_edges.entry(edge.to).or_default().push(edge.from);
        }

        // Check if two nodes are perfect path neighbors (following ODGI's logic)
        let are_perfect_neighbors = |left: usize, right: usize| -> bool {
            // Count visits to left node and check if they all continue to right
            let mut left_visits = 0;
            let mut left_to_right_visits = 0;

            // Check all paths that visit left node
            for (_path_name, path) in &self.paths {
                // Find all occurrences of left node in path
                for (i, &node) in path.iter().enumerate() {
                    if node == left {
                        left_visits += 1;
                        // Check if next node in path is right
                        if i + 1 < path.len() {
                            if path[i + 1] == right {
                                left_to_right_visits += 1;
                            } else {
                                // Path continues but doesn't go to right - not perfect neighbors
                                return false;
                            }
                        } else {
                            // Path ends at left node - not perfect neighbors
                            return false;
                        }
                    }
                }
            }

            // Count actual visits to right node
            let mut right_visits = 0;
            for (_, path) in &self.paths {
                for &node in path {
                    if node == right {
                        right_visits += 1;
                    }
                }
            }

            // They are perfect neighbors if:
            // 1. All paths visiting left continue to right
            // 2. The number of visits to right equals the number of transitions from left
            left_visits > 0
                && left_to_right_visits == left_visits
                && left_to_right_visits == right_visits
        };

        // Find simple chains using union-find approach
        let mut parent: HashMap<usize, usize> = HashMap::new();
        for &node_id in self.nodes.keys() {
            parent.insert(node_id, node_id);
        }

        let find = |parent: &mut HashMap<usize, usize>, x: usize| -> usize {
            let mut root = x;
            while parent[&root] != root {
                root = parent[&root];
            }
            // Path compression
            let mut curr = x;
            while parent[&curr] != root {
                let next = parent[&curr];
                parent.insert(curr, root);
                curr = next;
            }
            root
        };

        let unite = |parent: &mut HashMap<usize, usize>, x: usize, y: usize| {
            let root_x = find(parent, x);
            let root_y = find(parent, y);
            if root_x != root_y {
                parent.insert(root_x, root_y);
            }
        };

        // Count degree statistics
        let mut degree_stats = HashMap::new();
        let mut eligible_pairs = 0;
        let mut perfect_neighbor_pairs = 0;

        // Unite nodes that form simple chains
        for &node_id in self.nodes.keys() {
            let out_degree = forward_edges.get(&node_id).map(|v| v.len()).unwrap_or(0);
            let in_degree = backward_edges.get(&node_id).map(|v| v.len()).unwrap_or(0);

            *degree_stats.entry((in_degree, out_degree)).or_insert(0) += 1;

            // Check forward direction
            if out_degree == 1 {
                if let Some(nexts) = forward_edges.get(&node_id) {
                    let next = nexts[0];
                    if node_id != next {
                        let next_in_degree =
                            backward_edges.get(&next).map(|v| v.len()).unwrap_or(0);
                        if next_in_degree == 1 {
                            eligible_pairs += 1;
                            if are_perfect_neighbors(node_id, next) {
                                perfect_neighbor_pairs += 1;
                                unite(&mut parent, node_id, next);
                            }
                        }
                    }
                }
            }
        }

        if debug {
            println!("DEBUG: Node degree statistics:");
            let mut sorted_degrees: Vec<_> = degree_stats.iter().collect();
            sorted_degrees.sort_by_key(|&(_, count)| -{ *count });
            for ((in_deg, out_deg), count) in sorted_degrees.iter().take(10) {
                println!("  (in={}, out={}): {} nodes", in_deg, out_deg, count);
            }
            println!("  Eligible pairs (degree 1->1): {}", eligible_pairs);
            println!("  Perfect neighbor pairs: {}", perfect_neighbor_pairs);
        }

        // Collect components
        let mut component_map: HashMap<usize, Vec<usize>> = HashMap::new();
        for &node_id in self.nodes.keys() {
            let root = find(&mut parent, node_id);
            component_map.entry(root).or_default().push(node_id);
        }

        if debug {
            println!("DEBUG: Total components found: {}", component_map.len());
            let mut component_sizes: HashMap<usize, usize> = HashMap::new();
            for comp in component_map.values() {
                *component_sizes.entry(comp.len()).or_insert(0) += 1;
            }
            println!("  Component size distribution:");
            let mut sorted_sizes: Vec<_> = component_sizes.iter().collect();
            sorted_sizes.sort_by_key(|&(size, _)| -(*size as i32));
            for (size, count) in sorted_sizes.iter().take(10) {
                println!("    Size {}: {} components", size, count);
            }
        }

        // For simple chains, we need proper ordering
        // But for complex components, we should still compact them
        let mut ordered_components = 0;
        let mut unordered_components = 0;

        for (_, mut comp) in component_map {
            if comp.len() >= 2 {
                // Try to find a simple linear ordering
                let mut start = None;
                for &node in &comp {
                    let in_degree = backward_edges.get(&node).map(|v| v.len()).unwrap_or(0);
                    // Check if this could be the start of a linear chain
                    if in_degree == 0
                        || (in_degree == 1
                            && backward_edges
                                .get(&node)
                                .and_then(|v| v.first())
                                .map(|&prev| !comp.contains(&prev))
                                .unwrap_or(false))
                    {
                        start = Some(node);
                        break;
                    }
                }

                let mut use_ordering = false;
                if let Some(start_node) = start {
                    // Try to order the component by following edges
                    let mut ordered = vec![start_node];
                    let mut current = start_node;
                    let mut visited = HashSet::new();
                    visited.insert(start_node);

                    while ordered.len() < comp.len() {
                        let mut found = false;
                        if let Some(nexts) = forward_edges.get(&current) {
                            // Look for a next node that's in the component and not visited
                            for &next in nexts {
                                if comp.contains(&next) && !visited.contains(&next) {
                                    // Check if this next node has only one incoming edge from within component
                                    let next_in_from_comp = backward_edges
                                        .get(&next)
                                        .map(|prevs| {
                                            prevs.iter().filter(|&&p| comp.contains(&p)).count()
                                        })
                                        .unwrap_or(0);

                                    if next_in_from_comp == 1 {
                                        ordered.push(next);
                                        visited.insert(next);
                                        current = next;
                                        found = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if !found {
                            break;
                        }
                    }

                    if ordered.len() == comp.len() {
                        components.push(ordered);
                        use_ordering = true;
                        ordered_components += 1;
                    }
                }

                // If we couldn't order it as a simple chain, still compact it
                // but use the original order (which might not be path-coherent)
                if !use_ordering {
                    // For complex components, we need a different approach
                    // For now, just sort by node ID to have a consistent order
                    comp.sort_unstable();
                    let comp_len = comp.len();
                    components.push(comp);
                    unordered_components += 1;

                    if debug && comp_len > 10 {
                        println!(
                            "DEBUG: Using unordered compaction for component of size {}",
                            comp_len
                        );
                    }
                }
            }
        }

        if debug {
            println!(
                "DEBUG: {} ordered components, {} unordered components",
                ordered_components, unordered_components
            );
        }

        components
    }

    /// Compact nodes using ODGI-style algorithm
    pub fn compact_nodes_odgi(&mut self) -> usize {
        // Enable debug for troubleshooting
        let debug = std::env::var("SEQRUSH_DEBUG_COMPACT").is_ok();
        self.compact_nodes_odgi_internal(debug)
    }

    pub fn compact_nodes_odgi_internal(&mut self, debug: bool) -> usize {
        // Save original path sequences for validation
        let original_path_sequences: HashMap<String, Vec<u8>> = self
            .paths
            .iter()
            .map(|(name, path)| {
                let mut seq = Vec::new();
                for &node_id in path {
                    if let Some(node) = self.nodes.get(&node_id) {
                        seq.extend_from_slice(&node.sequence);
                    }
                }
                (name.clone(), seq)
            })
            .collect();

        let components = self.find_simple_components();
        let mut compacted_count = 0;

        if debug {
            println!("DEBUG: Found {} components to compact", components.len());
            for (i, comp) in components.iter().enumerate() {
                if i < 5 {
                    // Show first 5 components
                    println!("  Component {}: {} nodes: {:?}", i, comp.len(), comp);
                }
            }
            if components.len() > 5 {
                println!("  ... and {} more components", components.len() - 5);
            }
        }

        for component in components {
            if component.len() < 2 {
                continue;
            }

            // Merge this component
            let new_id = component[0]; // Use first node's ID
            let mut new_sequence = Vec::new();

            // Concatenate sequences
            for &node_id in &component {
                if let Some(node) = self.nodes.get(&node_id) {
                    new_sequence.extend_from_slice(&node.sequence);
                }
            }

            // Update the first node
            if let Some(node) = self.nodes.get_mut(&new_id) {
                node.sequence = new_sequence;
            }

            // Remove other nodes
            for &node_id in &component[1..] {
                self.nodes.remove(&node_id);
            }

            // Update paths - this is the critical part
            // According to ODGI, components should only be merged if paths traverse them
            // in the exact same order
            for (_path_name, path) in &mut self.paths {
                let mut new_path = Vec::new();
                let mut i = 0;

                while i < path.len() {
                    // Check if we're at the start of this component in the path
                    if path[i] == component[0] && i + component.len() <= path.len() {
                        // Check if the entire component appears consecutively
                        let mut matches = true;
                        for j in 0..component.len() {
                            if path[i + j] != component[j] {
                                matches = false;
                                break;
                            }
                        }

                        if matches {
                            // Replace the entire component with the merged node
                            new_path.push(new_id);
                            i += component.len();
                        } else {
                            // Component doesn't appear in full - this path doesn't traverse
                            // the component in the expected order, so we shouldn't have
                            // included these nodes in the component. Keep original.
                            new_path.push(path[i]);
                            i += 1;
                        }
                    } else {
                        // Not at component start - keep original
                        new_path.push(path[i]);
                        i += 1;
                    }
                }

                *path = new_path;
            }

            // Update edges
            let component_set: HashSet<usize> = component.iter().cloned().collect();
            let mut new_edges = HashSet::new();

            for edge in &self.edges {
                let mut from = edge.from;
                let mut to = edge.to;

                // Map nodes in component to the merged node
                if component_set.contains(&from) {
                    from = new_id;
                }
                if component_set.contains(&to) {
                    to = new_id;
                }

                // Skip internal edges within the component
                if from == new_id && to == new_id {
                    // Check if this was an internal edge
                    let mut is_internal = false;
                    for i in 0..component.len() - 1 {
                        if edge.from == component[i] && edge.to == component[i + 1] {
                            is_internal = true;
                            break;
                        }
                    }
                    if is_internal {
                        continue;
                    }
                }

                new_edges.insert(Edge { from, to });
            }

            self.edges = new_edges;
            compacted_count += component.len() - 1;
        }

        // Validate that path sequences are preserved
        for (path_name, path) in &self.paths {
            let mut reconstructed = Vec::new();
            for &node_id in path {
                if let Some(node) = self.nodes.get(&node_id) {
                    reconstructed.extend_from_slice(&node.sequence);
                }
            }

            if let Some(original) = original_path_sequences.get(path_name) {
                if reconstructed != *original {
                    eprintln!(
                        "WARNING: Path {} sequence changed during compaction!",
                        path_name
                    );
                    eprintln!(
                        "  Original length: {}, Reconstructed length: {}",
                        original.len(),
                        reconstructed.len()
                    );
                }
            }
        }

        compacted_count
    }
}
