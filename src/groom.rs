/// Grooming algorithm to remove spurious inverting links from the graph
/// by exploring from the orientation supported by the most paths.
///
/// Based on ODGI's grooming algorithm
use crate::bidirected_graph::Handle;
use crate::bidirected_ops::BidirectedGraph;
use std::collections::{HashSet, VecDeque};

impl BidirectedGraph {
    /// Analyze orientation preferences from paths
    /// Returns (forward_count, reverse_count) for each node
    fn analyze_orientation_preferences(&self) -> std::collections::HashMap<usize, (usize, usize)> {
        let mut preferences = std::collections::HashMap::new();

        for path in &self.paths {
            for &handle in &path.steps {
                let entry = preferences.entry(handle.node_id()).or_insert((0, 0));
                if handle.is_reverse() {
                    entry.1 += 1; // reverse count
                } else {
                    entry.0 += 1; // forward count
                }
            }
        }

        preferences
    }

    /// Count how many paths use each edge (coverage)
    /// Returns map of (from_handle, to_handle) -> path_count
    fn count_edge_coverage(&self) -> std::collections::HashMap<(Handle, Handle), usize> {
        let mut coverage = std::collections::HashMap::new();

        for path in &self.paths {
            for i in 0..path.steps.len().saturating_sub(1) {
                let from = path.steps[i];
                let to = path.steps[i + 1];
                *coverage.entry((from, to)).or_insert(0) += 1;
            }
        }

        coverage
    }

    /// Apply grooming to orient the graph consistently
    /// This explores the graph via BFS/DFS from seed nodes and flips nodes
    /// as needed to maintain consistent orientation along paths
    /// Returns handles in their CURRENT order (like ODGI), just with flipped orientations
    pub fn groom(&mut self, use_bfs: bool, verbose: bool) -> Vec<Handle> {
        self.groom_with_mode(use_bfs, false, verbose)
    }

    /// Apply grooming with explicit mode selection
    /// use_bfs: if true, use BFS; if false and use_coverage_dfs is false, use regular DFS
    /// use_coverage_dfs: if true, use coverage-weighted DFS (overrides use_bfs)
    pub fn groom_with_mode(&mut self, use_bfs: bool, use_coverage_dfs: bool, verbose: bool) -> Vec<Handle> {
        if verbose {
            eprintln!("[groom] Starting grooming with {} nodes", self.nodes.len());
            if use_coverage_dfs {
                eprintln!("[groom] Using coverage-weighted DFS mode");
            } else if use_bfs {
                eprintln!("[groom] Using BFS mode");
            } else {
                eprintln!("[groom] Using regular DFS mode");
            }
        }

        // Analyze orientation preferences from paths (majority-based approach)
        let preferences = self.analyze_orientation_preferences();

        // Count edge coverage if using coverage-weighted DFS
        let edge_coverage = if use_coverage_dfs {
            self.count_edge_coverage()
        } else {
            std::collections::HashMap::new()
        };

        if verbose {
            let total_nodes = preferences.len();
            let unanimous = preferences.values()
                .filter(|(fwd, rev)| *fwd == 0 || *rev == 0)
                .count();
            let conflicts = total_nodes - unanimous;
            eprintln!("[groom] Orientation analysis: {} unanimous, {} conflicts",
                     unanimous, conflicts);

            // Show details of conflicts
            if conflicts > 0 && conflicts <= 10 {
                eprintln!("[groom] Conflict nodes:");
                let mut conflict_nodes: Vec<_> = preferences.iter()
                    .filter(|(_, (fwd, rev))| *fwd > 0 && *rev > 0)
                    .collect();
                conflict_nodes.sort_by_key(|(node_id, _)| **node_id);
                for (node_id, (fwd, rev)) in conflict_nodes {
                    eprintln!("[groom]   Node {}: {} forward, {} reverse (majority: {})",
                             node_id, fwd, rev,
                             if fwd > rev { "forward" } else { "reverse" });
                }
            }

            if use_coverage_dfs {
                eprintln!("[groom] Counted coverage for {} edges", edge_coverage.len());
            }
        }

        // Find seed nodes - use graph structure as-is (like ODGI)
        let seeds = self.find_head_nodes();

        if verbose {
            eprintln!("[groom] Found {} head nodes as seeds", seeds.len());
        }

        // Track visited and flipped status
        let mut visited = HashSet::new();
        let mut flipped = HashSet::new();

        // If we have no heads, pick first unvisited node (like ODGI)
        let mut current_seeds = if seeds.is_empty() {
            // Pick first node in forward orientation
            if let Some(node_id) = self.nodes.iter().enumerate()
                .filter_map(|(id, n)| if n.is_some() { Some(id) } else { None })
                .min() {
                if verbose {
                    eprintln!("[groom] No heads found, using first node {}", node_id);
                }
                vec![Handle::forward(node_id)]
            } else {
                vec![]
            }
        } else {
            seeds
        };

        // Process all components
        let mut component_count = 0;
        while visited.len() < self.nodes.len() {
            if current_seeds.is_empty() {
                // Find an unvisited node as new seed (forward orientation, like ODGI)
                for (node_id, node_opt) in self.nodes.iter().enumerate() {
                    if node_opt.is_none() {
                        continue;
                    }
                    if !visited.contains(&node_id) {
                        let seed = Handle::forward(node_id);  // Always use forward (like ODGI)

                        current_seeds.push(seed);
                        if verbose {
                            eprintln!("[groom] Starting new component with node {}+", node_id);
                        }
                        component_count += 1;
                        break;
                    }
                }
                if current_seeds.is_empty() {
                    break; // All nodes visited
                }
            }

            let visited_before = visited.len();
            if use_coverage_dfs {
                self.groom_coverage_weighted_dfs(&current_seeds, &mut visited, &mut flipped, &edge_coverage, verbose);
            } else if use_bfs {
                self.groom_bfs_majority(&current_seeds, &mut visited, &mut flipped, &preferences, verbose);
            } else {
                self.groom_dfs(&current_seeds, &mut visited, &mut flipped);
            }

            if verbose && visited.len() > visited_before {
                eprintln!("[groom] Component {} visited {} nodes",
                         component_count, visited.len() - visited_before);
            }

            current_seeds.clear();
        }

        // Build the final handle order in CURRENT NODE ORDER (like ODGI)
        // Get nodes in sorted order
        let mut node_ids: Vec<_> = self.nodes.iter().enumerate()
            .filter_map(|(id, n)| if n.is_some() { Some(id) } else { None })
            .collect();
        node_ids.sort();

        let mut final_order = Vec::new();
        for node_id in node_ids {
            let handle = if flipped.contains(&node_id) {
                Handle::reverse(node_id)
            } else {
                Handle::forward(node_id)
            };
            final_order.push(handle);
        }

        if verbose {
            eprintln!("[groom] Flipped {} nodes", flipped.len());
        }

        final_order
    }

    /// Find head nodes and orient them according to majority preference
    fn find_head_nodes_with_majority(&self, preferences: &std::collections::HashMap<usize, (usize, usize)>) -> Vec<Handle> {
        let heads = self.find_head_nodes();

        // Re-orient based on majority preference
        heads.into_iter().map(|h| {
            let node_id = h.node_id();
            let majority_forward = preferences.get(&node_id)
                .map(|(fwd, rev)| fwd >= rev)
                .unwrap_or(true);

            if majority_forward {
                Handle::forward(node_id)
            } else {
                Handle::reverse(node_id)
            }
        }).collect()
    }

    /// Pick seeds by strongest orientation preference
    fn pick_seeds_by_preference(&self, preferences: &std::collections::HashMap<usize, (usize, usize)>,
                                verbose: bool) -> Vec<Handle> {
        // Find node with strongest preference (largest difference between fwd and rev)
        let best = preferences.iter()
            .max_by_key(|(_, (fwd, rev))| {
                if fwd > rev { fwd - rev } else { rev - fwd }
            });

        if let Some((node_id, (fwd, rev))) = best {
            let majority_forward = fwd >= rev;
            if verbose {
                eprintln!("[groom] No heads found, using node {} with strongest preference ({} fwd, {} rev)",
                         node_id, fwd, rev);
            }
            vec![if majority_forward {
                Handle::forward(*node_id)
            } else {
                Handle::reverse(*node_id)
            }]
        } else if let Some(node_id) = self.nodes.iter().enumerate()
            .filter_map(|(id, n)| if n.is_some() { Some(id) } else { None })
            .min() {
            if verbose {
                eprintln!("[groom] No preferences found, using first node {}", node_id);
            }
            vec![Handle::forward(node_id)]
        } else {
            vec![]
        }
    }

    /// BFS-based grooming - ODGI greedy algorithm
    /// Just uses edge orientation on first visit, no fancy majority voting
    fn groom_bfs_majority(&self, seeds: &[Handle], visited: &mut HashSet<usize>,
                         flipped: &mut HashSet<usize>,
                         preferences: &std::collections::HashMap<usize, (usize, usize)>,
                         verbose: bool) {
        let mut queue = VecDeque::new();
        let verbose_groom = std::env::var("SEQRUSH_GROOM_DEBUG").is_ok();

        // Initialize queue with seeds
        for &seed in seeds {
            if !visited.contains(&seed.node_id()) {
                queue.push_back(seed);
                visited.insert(seed.node_id());

                if verbose_groom {
                    eprintln!("[groom_bfs] Seed node {} orientation {}",
                             seed.node_id(), if seed.is_reverse() { "REVERSE" } else { "FORWARD" });
                }

                // ODGI logic: flip if we visit via reverse orientation
                if seed.is_reverse() {
                    flipped.insert(seed.node_id());
                    if verbose_groom {
                        eprintln!("[groom_bfs]   -> Marking {} as FLIPPED (seed was reverse)", seed.node_id());
                    }
                }
            }
        }

        // BFS traversal - pure ODGI greedy algorithm
        while let Some(current) = queue.pop_front() {
            // Get edges from this handle and sort them for deterministic iteration (like ODGI)
            let mut matching_edges: Vec<_> = self.edges.iter()
                .filter(|e| e.from == current)
                .collect();
            matching_edges.sort();  // Use BiEdge's Ord trait for deterministic order

            for edge in matching_edges {
                let next = edge.to;
                if !visited.contains(&next.node_id()) {
                    visited.insert(next.node_id());

                    if verbose_groom {
                        eprintln!("[groom_bfs] From {}{}  ->  {}{}",
                                 current.node_id(), if current.is_reverse() { "-" } else { "+" },
                                 next.node_id(), if next.is_reverse() { "-" } else { "+" });
                    }

                    // ODGI greedy logic: mark as flipped if we reach it via reverse orientation
                    if next.is_reverse() {
                        flipped.insert(next.node_id());
                        if verbose_groom {
                            eprintln!("[groom_bfs]   -> Marking {} as FLIPPED (reached via reverse)", next.node_id());
                        }
                    }

                    // Continue BFS from the handle we arrived at
                    queue.push_back(next);
                }
            }
        }
    }

    /// BFS-based grooming with traversal order tracking
    fn groom_bfs_with_order(&self, seeds: &[Handle], visited: &mut HashSet<usize>,
                            flipped: &mut HashSet<usize>, order: &mut Vec<usize>) {
        let mut queue = VecDeque::new();
        let verbose_groom = std::env::var("SEQRUSH_GROOM_DEBUG").is_ok();

        // Initialize queue with seeds
        for &seed in seeds {
            if !visited.contains(&seed.node_id()) {
                queue.push_back(seed);
                visited.insert(seed.node_id());
                order.push(seed.node_id());  // Track traversal order

                if verbose_groom {
                    eprintln!("[groom_bfs] Seed node {} orientation {}",
                             seed.node_id(), if seed.is_reverse() { "REVERSE" } else { "FORWARD" });
                }

                // ODGI logic: flip if we visit via reverse orientation
                if seed.is_reverse() {
                    flipped.insert(seed.node_id());
                    if verbose_groom {
                        eprintln!("[groom_bfs]   -> Marking {} as FLIPPED (seed was reverse)", seed.node_id());
                    }
                }
            }
        }

        // BFS traversal
        while let Some(current) = queue.pop_front() {
            // Get edges from this handle
            for edge in &self.edges {
                if edge.from == current {
                    let next = edge.to;
                    if !visited.contains(&next.node_id()) {
                        visited.insert(next.node_id());
                        order.push(next.node_id());  // Track traversal order

                        if verbose_groom {
                            eprintln!("[groom_bfs] From {}{}  ->  {}{}",
                                     current.node_id(), if current.is_reverse() { "-" } else { "+" },
                                     next.node_id(), if next.is_reverse() { "-" } else { "+" });
                        }

                        // ODGI logic: mark as flipped if we reach it via reverse orientation
                        if next.is_reverse() {
                            flipped.insert(next.node_id());
                            if verbose_groom {
                                eprintln!("[groom_bfs]   -> Marking {} as FLIPPED (reached via reverse)", next.node_id());
                            }
                        }

                        // Continue BFS from the handle we arrived at
                        queue.push_back(next);
                    }
                }
            }
        }
    }

    /// DFS-based grooming with traversal order tracking
    fn groom_dfs_with_order(&self, seeds: &[Handle], visited: &mut HashSet<usize>,
                            flipped: &mut HashSet<usize>, order: &mut Vec<usize>) {
        let mut stack = Vec::new();

        // Initialize stack with seeds
        for &seed in seeds {
            stack.push(seed);
        }

        // DFS traversal
        while let Some(current) = stack.pop() {
            if visited.contains(&current.node_id()) {
                continue;
            }

            visited.insert(current.node_id());
            order.push(current.node_id());  // Track traversal order

            if current.is_reverse() {
                flipped.insert(current.node_id());
            }

            // Get edges from this handle
            for edge in &self.edges {
                if edge.from == current {
                    let next = edge.to;
                    if !visited.contains(&next.node_id()) {
                        stack.push(next);
                    }
                }
            }
        }
    }

    /// Keep old BFS for compatibility
    fn groom_bfs(&self, seeds: &[Handle], visited: &mut HashSet<usize>, flipped: &mut HashSet<usize>) {
        let mut _order = Vec::new();
        self.groom_bfs_with_order(seeds, visited, flipped, &mut _order);
    }

    /// Keep old DFS for compatibility
    fn groom_dfs(&self, seeds: &[Handle], visited: &mut HashSet<usize>, flipped: &mut HashSet<usize>) {
        let mut _order = Vec::new();
        self.groom_dfs_with_order(seeds, visited, flipped, &mut _order);
    }

    /// Coverage-weighted DFS grooming
    /// Prioritizes edges with higher coverage (used by more paths)
    /// This naturally follows the "main bundle" and handles hairpins better
    fn groom_coverage_weighted_dfs(&self, seeds: &[Handle], visited: &mut HashSet<usize>,
                                   flipped: &mut HashSet<usize>,
                                   edge_coverage: &std::collections::HashMap<(Handle, Handle), usize>,
                                   verbose: bool) {
        let mut stack = Vec::new();
        let verbose_groom = std::env::var("SEQRUSH_GROOM_DEBUG").is_ok();

        // Initialize stack with seeds
        for &seed in seeds {
            stack.push(seed);
        }

        // DFS traversal
        while let Some(current) = stack.pop() {
            if visited.contains(&current.node_id()) {
                continue;
            }

            visited.insert(current.node_id());

            if verbose_groom {
                eprintln!("[groom_cov_dfs] Visiting node {}{}",
                         current.node_id(),
                         if current.is_reverse() { "-" } else { "+" });
            }

            // ODGI logic: flip if we visit via reverse orientation
            if current.is_reverse() {
                flipped.insert(current.node_id());
                if verbose_groom {
                    eprintln!("[groom_cov_dfs]   -> Marking {} as FLIPPED", current.node_id());
                }
            }

            // Get outgoing edges and sort by coverage (highest first)
            let mut outgoing: Vec<_> = self.edges.iter()
                .filter(|e| e.from == current)
                .map(|e| {
                    let coverage = edge_coverage.get(&(e.from, e.to)).copied().unwrap_or(0);
                    (e.to, coverage)
                })
                .collect();

            // Sort by coverage descending, then by node ID for determinism
            outgoing.sort_by_key(|(handle, cov)| (std::cmp::Reverse(*cov), handle.node_id(), handle.is_reverse()));

            if verbose_groom && !outgoing.is_empty() {
                eprintln!("[groom_cov_dfs]   Outgoing edges (sorted by coverage):");
                for (next, cov) in &outgoing {
                    eprintln!("[groom_cov_dfs]     -> {}{} (coverage: {})",
                             next.node_id(),
                             if next.is_reverse() { "-" } else { "+" },
                             cov);
                }
            }

            // Push in reverse order so highest-coverage is processed first (DFS stack)
            for (next_handle, _) in outgoing.into_iter().rev() {
                if !visited.contains(&next_handle.node_id()) {
                    stack.push(next_handle);
                }
            }
        }
    }

    /// Apply grooming then topological sort
    pub fn groom_and_sort(&mut self, verbose: bool) {
        if verbose {
            eprintln!("[groom_and_sort] Starting groom and sort process");
        }

        // Step 1: Apply grooming
        let groomed_order = self.groom(true, verbose); // Use BFS

        // Apply the grooming orientation changes
        self.apply_grooming(groomed_order.clone(), verbose);

        // Step 2: Apply topological sort
        if verbose {
            eprintln!("[groom_and_sort] Applying topological sort after grooming");
        }
        self.apply_exact_odgi_ordering(verbose);

        if verbose {
            eprintln!("[groom_and_sort] Groom and sort complete");
        }
    }

    /// Alternative: Sort first, then groom, then sort again
    pub fn sort_groom_sort(&mut self, verbose: bool) {
        if verbose {
            eprintln!("[sort_groom_sort] Starting sort-groom-sort process");
        }

        // Step 1: Initial topological sort
        if verbose {
            eprintln!("[sort_groom_sort] Initial topological sort");
        }
        self.apply_exact_odgi_ordering(false);

        // Step 2: Apply grooming with coverage-weighted DFS
        if verbose {
            eprintln!("[sort_groom_sort] Grooming after first sort (using coverage-weighted DFS)");
        }
        let groomed_order = self.groom_with_mode(false, true, false); // Use coverage-weighted DFS
        self.apply_grooming(groomed_order, false);

        // Step 3: Final topological sort
        if verbose {
            eprintln!("[sort_groom_sort] Final topological sort");
        }
        self.apply_exact_odgi_ordering(false);

        if verbose {
            eprintln!("[sort_groom_sort] Sort-groom-sort complete");
        }
    }

    /// Iterative sort-groom-sort until stabilization or max iterations
    pub fn iterative_groom(&mut self, max_iterations: usize, verbose: bool) -> usize {
        if verbose {
            eprintln!("[iterative_groom] Starting iterative grooming (max {} iterations)", max_iterations);
        }

        let mut iteration = 0;
        let mut prev_flipped_count = usize::MAX;

        while iteration < max_iterations {
            iteration += 1;

            if verbose {
                eprintln!("[iterative_groom] === Iteration {} ===", iteration);
            }

            // Step 1: Sort
            if verbose {
                eprintln!("[iterative_groom] Sorting...");
            }
            self.apply_exact_odgi_ordering(false);

            // Step 2: Groom and count how many nodes were flipped
            if verbose {
                eprintln!("[iterative_groom] Grooming...");
            }
            let groomed_order = self.groom(true, false);

            // Count how many nodes need flipping
            let flipped_count = groomed_order.iter()
                .filter(|h| h.is_reverse())
                .count();

            if verbose {
                eprintln!("[iterative_groom] Iteration {} flipped {} nodes", iteration, flipped_count);
            }

            // Apply the grooming
            self.apply_grooming(groomed_order, false);

            // Step 3: Final sort
            if verbose {
                eprintln!("[iterative_groom] Final sort...");
            }
            self.apply_exact_odgi_ordering(false);

            // Check for stabilization
            if flipped_count == prev_flipped_count || flipped_count == 0 {
                if verbose {
                    eprintln!("[iterative_groom] Stabilized after {} iterations", iteration);
                }
                break;
            }

            prev_flipped_count = flipped_count;
        }

        if verbose && iteration == max_iterations {
            eprintln!("[iterative_groom] Reached maximum iterations ({})", max_iterations);
        }

        iteration
    }

    /// Apply grooming changes (flip nodes as needed and optionally reorder)
    fn apply_grooming(&mut self, groomed_handles: Vec<Handle>, verbose: bool) {
        self.apply_grooming_with_reorder(groomed_handles, false, verbose);
    }

    /// Apply grooming changes with optional node reordering (like ODGI does)
    pub fn apply_grooming_with_reorder(&mut self, groomed_handles: Vec<Handle>, reorder: bool, verbose: bool) {
        // Build a map of which nodes need to be flipped
        let mut node_flips = HashSet::new();

        for handle in &groomed_handles {
            if handle.is_reverse() {
                node_flips.insert(handle.node_id());
            }
        }

        if verbose && !node_flips.is_empty() {
            eprintln!("[apply_grooming] Flipping {} nodes", node_flips.len());
        }

        // Flip the sequences of nodes that need flipping
        for &node_id in &node_flips {
            if let Some(Some(node)) = self.nodes.get_mut(node_id) {
                // Reverse complement the sequence
                node.sequence = crate::bidirected_graph::reverse_complement(&node.sequence);
            }
        }

        // Update edges: when a node is flipped, we XOR the orientation of handles to/from it
        let mut new_edges = HashSet::new();
        for edge in self.edges.drain() {
            // XOR orientations if the node is flipped
            let new_from = if node_flips.contains(&edge.from.node_id()) {
                edge.from.flip()
            } else {
                edge.from
            };

            let new_to = if node_flips.contains(&edge.to.node_id()) {
                edge.to.flip()
            } else {
                edge.to
            };

            new_edges.insert(crate::bidirected_graph::BiEdge::new(new_from, new_to));
        }
        self.edges = new_edges;

        // Update paths: XOR the orientation when a node is flipped
        for path in &mut self.paths {
            for handle in &mut path.steps {
                if node_flips.contains(&handle.node_id()) {
                    *handle = handle.flip();
                }
            }
        }

        // Apply node reordering if requested (like ODGI does)
        if reorder {
            if verbose {
                eprintln!("[apply_grooming] Applying node reordering based on traversal order");
            }

            // Create mapping from old node IDs to new ones based on traversal order
            let mut id_mapping = std::collections::HashMap::new();
            for (new_id, handle) in groomed_handles.iter().enumerate() {
                let old_id = handle.node_id();
                let new_node_id = new_id + 1; // Node IDs start from 1
                id_mapping.insert(old_id, new_node_id);
            }

            // Apply the reordering
            self.apply_node_id_mapping(&id_mapping);
        }

        if verbose {
            eprintln!("[apply_grooming] Grooming complete");
        }
    }
}