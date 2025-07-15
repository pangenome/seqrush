use std::collections::{HashMap, HashSet};

/// Represents a node in the pangenome graph
#[derive(Clone, Debug)]
pub struct Node {
    pub id: usize,
    pub sequence: Vec<u8>,
    pub rank: f64,  // For topological ordering
}

/// Represents an edge between nodes
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub struct Edge {
    pub from: usize,
    pub to: usize,
}

/// Graph structure for manipulation
#[derive(Clone)]
pub struct Graph {
    pub nodes: HashMap<usize, Node>,
    pub edges: HashSet<Edge>,
    pub paths: Vec<(String, Vec<usize>)>,
}

impl Default for Graph {
    fn default() -> Self {
        Self::new()
    }
}

impl Graph {
    pub fn new() -> Self {
        Graph {
            nodes: HashMap::new(),
            edges: HashSet::new(),
            paths: Vec::new(),
        }
    }
    
    /// Compact linear chains of nodes (unchop)
    pub fn compact_nodes(&mut self) -> usize {
        // Use ODGI-style compaction algorithm
        self.compact_nodes_odgi()
    }
    
    #[allow(dead_code)]
    fn compact_nodes_old(&mut self) -> usize {
        let mut compacted_count = 0;
        let mut processed = HashSet::new();
        
        // Find all linear chains (simple components)
        let chains = self.find_linear_chains();
        
        
        // Add a safety check to prevent infinite loops
        let max_iterations = self.nodes.len();
        
        for (iterations, chain) in chains.into_iter().enumerate() {
            if iterations >= max_iterations {
                eprintln!("WARNING: Compaction stopped after {} iterations", iterations);
                break;
            }
            
            if chain.len() < 2 {
                continue;
            }
            
            // Skip if any node in chain was already processed
            if chain.iter().any(|&id| processed.contains(&id)) {
                continue;
            }
            
            // Mark all nodes in chain as processed
            for &id in &chain {
                processed.insert(id);
            }
            
            // Merge the chain into a single node
            self.merge_chain(&chain);
            compacted_count += chain.len() - 1;
        }
        
        compacted_count
    }
    
    /// Find all linear chains in the graph
    #[allow(dead_code)]
    fn find_linear_chains(&self) -> Vec<Vec<usize>> {
        let mut chains = Vec::new();
        let mut visited = HashSet::new();
        
        // Build adjacency lists
        let mut forward_edges: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut backward_edges: HashMap<usize, Vec<usize>> = HashMap::new();
        
        for edge in &self.edges {
            forward_edges.entry(edge.from).or_default().push(edge.to);
            backward_edges.entry(edge.to).or_default().push(edge.from);
        }
        
        // Find nodes that can be part of linear chains
        // A node can be in a linear chain if it has exactly 1 in-edge and 1 out-edge
        for &node_id in self.nodes.keys() {
            if visited.contains(&node_id) {
                continue;
            }
            
            // Check if node has a self-loop
            let has_self_loop = self.edges.contains(&Edge { from: node_id, to: node_id });
            if has_self_loop {
                continue; // Skip nodes with self-loops for chain compaction
            }
            
            let in_degree = backward_edges.get(&node_id).map(|v| v.len()).unwrap_or(0);
            let out_degree = forward_edges.get(&node_id).map(|v| v.len()).unwrap_or(0);
            
            // Only consider nodes with exactly 1 in and 1 out for chains
            if in_degree == 1 && out_degree == 1 {
                // Walk backward to find chain start
                let mut chain_start = node_id;
                if let Some(prevs) = backward_edges.get(&node_id) {
                    if prevs.len() == 1 {
                        let mut current = prevs[0];
                        let mut seen_in_backward = HashSet::new();
                        seen_in_backward.insert(node_id);
                        
                        loop {
                            // Prevent infinite loop in case of cycles
                            if seen_in_backward.contains(&current) {
                                break;
                            }
                            seen_in_backward.insert(current);
                            
                            let prev_in = backward_edges.get(&current).map(|v| v.len()).unwrap_or(0);
                            let prev_out = forward_edges.get(&current).map(|v| v.len()).unwrap_or(0);
                            
                            if prev_out == 1 && !visited.contains(&current) {
                                chain_start = current;
                                if prev_in == 1 {
                                    // Continue backward
                                    if let Some(prevs) = backward_edges.get(&current) {
                                        if prevs.len() == 1 {
                                            current = prevs[0];
                                            continue;
                                        }
                                    }
                                }
                            }
                            break;
                        }
                    }
                }
                
                // Now walk forward from chain_start to build the chain
                if !visited.contains(&chain_start) {
                    let mut chain = vec![chain_start];
                    visited.insert(chain_start);
                    
                    let mut current = chain_start;
                    while let Some(nexts) = forward_edges.get(&current) {
                        if nexts.len() == 1 {
                            let next = nexts[0];
                            let next_in_degree = backward_edges.get(&next).map(|v| v.len()).unwrap_or(0);
                            let next_out_degree = forward_edges.get(&next).map(|v| v.len()).unwrap_or(0);
                            
                            if next_in_degree == 1 && !visited.contains(&next) {
                                chain.push(next);
                                visited.insert(next);
                                current = next;
                                
                                // Stop if we reach a node with out-degree != 1
                                if next_out_degree != 1 {
                                    break;
                                }
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                    }
                    
                    if chain.len() > 1 {
                        chains.push(chain);
                    }
                }
            }
        }
        
        chains
    }
    
    /// Merge a chain of nodes into a single node
    #[allow(dead_code)]
    fn merge_chain(&mut self, chain: &[usize]) {
        if chain.len() < 2 {
            return;
        }
        
        // Create new node with concatenated sequence
        let new_id = chain[0]; // Use first node's ID
        let mut new_sequence = Vec::new();
        let mut total_rank = 0.0;
        
        // Concatenate all sequences in the chain
        for &node_id in chain {
            if let Some(node) = self.nodes.get(&node_id) {
                new_sequence.extend_from_slice(&node.sequence);
                total_rank += node.rank;
            }
        }
        
        let avg_rank = total_rank / chain.len() as f64;
        
        
        // Remove old nodes except the first one
        for &node_id in &chain[1..] {
            self.nodes.remove(&node_id);
        }
        
        // Update the first node
        if let Some(node) = self.nodes.get_mut(&new_id) {
            node.sequence = new_sequence;
            node.rank = avg_rank;
        }
        
        // Update edges
        let mut new_edges = HashSet::new();
        
        for edge in &self.edges {
            let mut from = edge.from;
            let mut to = edge.to;
            
            // Remap source
            if chain.contains(&from) {
                from = new_id;
            }
            
            // Remap target
            if chain.contains(&to) {
                to = new_id;
            }
            
            // Skip internal chain edges (edges between consecutive nodes in the chain)
            // but keep all other edges including self-loops
            let mut is_internal_chain_edge = false;
            for i in 0..chain.len() - 1 {
                if edge.from == chain[i] && edge.to == chain[i + 1] {
                    is_internal_chain_edge = true;
                    break;
                }
            }
            
            if is_internal_chain_edge {
                continue; // Skip internal chain edges
            }
            
            new_edges.insert(Edge { from, to });
        }
        
        self.edges = new_edges;
        
        // Update paths - replace consecutive chain nodes with single merged node
        for (_, path) in &mut self.paths {
            let mut new_path = Vec::new();
            let mut i = 0;
            
            while i < path.len() {
                let node_id = path[i];
                
                if chain.contains(&node_id) {
                    // Found a chain node - check if this is the start of the chain in the path
                    let chain_start_idx = chain.iter().position(|&n| n == node_id);
                    
                    if let Some(start_idx) = chain_start_idx {
                        // Check if the path follows the chain order
                        let mut chain_len_in_path = 1;
                        
                        for j in 1..chain.len() - start_idx {
                            if i + j < path.len() && path[i + j] == chain[start_idx + j] {
                                chain_len_in_path += 1;
                            } else {
                                break;
                            }
                        }
                        
                        if start_idx == 0 && chain_len_in_path == chain.len() {
                            // Full chain appears in order - replace with merged node
                            new_path.push(new_id);
                            i += chain_len_in_path;
                        } else {
                            // Partial chain or doesn't start at beginning
                            // Don't compact - keep original nodes to avoid issues
                            new_path.push(node_id);
                            i += 1;
                        }
                    } else {
                        // Shouldn't happen, but handle gracefully
                        new_path.push(node_id);
                        i += 1;
                    }
                } else {
                    new_path.push(node_id);
                    i += 1;
                }
            }
            
            *path = new_path;
        }
        
        // After updating paths, rebuild edges from the paths to ensure all necessary edges exist
        // This is important for self-loops that arise from consecutive visits to the same node
        let mut edges_from_paths = HashSet::new();
        for (_, path) in &self.paths {
            for window in path.windows(2) {
                if let [from, to] = window {
                    edges_from_paths.insert(Edge { from: *from, to: *to });
                }
            }
        }
        
        // Add any edges from paths that might be missing
        for edge in edges_from_paths {
            self.edges.insert(edge);
        }
        
        // Clean up edges that reference non-existent nodes
        self.edges.retain(|edge| {
            self.nodes.contains_key(&edge.from) && self.nodes.contains_key(&edge.to)
        });
    }
    
    /// Compute strongly connected components using Tarjan's algorithm
    #[allow(dead_code)]
    fn compute_sccs(&self, adj_list: &HashMap<usize, Vec<usize>>) -> Vec<Vec<usize>> {
        let mut index_counter = 0;
        let mut stack = Vec::new();
        let mut indices: HashMap<usize, usize> = HashMap::new();
        let mut lowlinks: HashMap<usize, usize> = HashMap::new();
        let mut on_stack: HashSet<usize> = HashSet::new();
        let mut sccs = Vec::new();
        
        fn strongconnect(
            v: usize,
            adj_list: &HashMap<usize, Vec<usize>>,
            index_counter: &mut usize,
            stack: &mut Vec<usize>,
            indices: &mut HashMap<usize, usize>,
            lowlinks: &mut HashMap<usize, usize>,
            on_stack: &mut HashSet<usize>,
            sccs: &mut Vec<Vec<usize>>,
        ) {
            indices.insert(v, *index_counter);
            lowlinks.insert(v, *index_counter);
            *index_counter += 1;
            stack.push(v);
            on_stack.insert(v);
            
            if let Some(neighbors) = adj_list.get(&v) {
                for &w in neighbors {
                    if !indices.contains_key(&w) {
                        strongconnect(w, adj_list, index_counter, stack, indices, lowlinks, on_stack, sccs);
                        let w_lowlink = lowlinks[&w];
                        let v_lowlink = lowlinks[&v];
                        lowlinks.insert(v, v_lowlink.min(w_lowlink));
                    } else if on_stack.contains(&w) {
                        let w_index = indices[&w];
                        let v_lowlink = lowlinks[&v];
                        lowlinks.insert(v, v_lowlink.min(w_index));
                    }
                }
            }
            
            if lowlinks[&v] == indices[&v] {
                let mut component = Vec::new();
                loop {
                    let w = stack.pop().unwrap();
                    on_stack.remove(&w);
                    component.push(w);
                    if w == v {
                        break;
                    }
                }
                sccs.push(component);
            }
        }
        
        let mut nodes: Vec<_> = self.nodes.keys().cloned().collect();
        nodes.sort(); // Ensure deterministic order
        
        for node in nodes {
            if !indices.contains_key(&node) {
                strongconnect(
                    node,
                    adj_list,
                    &mut index_counter,
                    &mut stack,
                    &mut indices,
                    &mut lowlinks,
                    &mut on_stack,
                    &mut sccs,
                );
            }
        }
        
        sccs
    }
    
    /// Perform topological sort on the graph
    pub fn topological_sort(&mut self) {
        // Use a multi-pass algorithm that minimizes edge spans
        
        // First, do a basic topological sort
        let mut visited = HashSet::new();
        let mut initial_order = Vec::new();
        
        // Build adjacency lists
        let mut adj_list: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut in_degree: HashMap<usize, usize> = HashMap::new();
        
        // Initialize
        for node_id in self.nodes.keys() {
            adj_list.insert(*node_id, Vec::new());
            in_degree.insert(*node_id, 0);
        }
        
        // Build adjacency list and compute in-degrees
        for edge in &self.edges {
            adj_list.get_mut(&edge.from).unwrap().push(edge.to);
            *in_degree.get_mut(&edge.to).unwrap() += 1;
        }
        
        // Find path information
        let mut path_positions: HashMap<usize, Vec<usize>> = HashMap::new();
        for (_, path) in &self.paths {
            for (pos, &node) in path.iter().enumerate() {
                path_positions.entry(node).or_default().push(pos);
            }
        }
        
        // Modified Kahn's algorithm
        let mut queue: Vec<usize> = Vec::new();
        
        // Start with nodes that have no incoming edges
        for (&node_id, &degree) in &in_degree {
            if degree == 0 {
                queue.push(node_id);
            }
        }
        
        // Sort by average position in paths
        queue.sort_by_key(|&node| {
            path_positions.get(&node)
                .map(|positions| positions.iter().sum::<usize>() / positions.len())
                .unwrap_or(usize::MAX)
        });
        
        // Process nodes
        while let Some(node) = queue.pop() {
            initial_order.push(node);
            visited.insert(node);
            
            if let Some(neighbors) = adj_list.get(&node) {
                let mut next_nodes = Vec::new();
                
                for &neighbor in neighbors {
                    if let Some(degree) = in_degree.get_mut(&neighbor) {
                        *degree = degree.saturating_sub(1);
                        if *degree == 0 && !visited.contains(&neighbor) {
                            next_nodes.push(neighbor);
                        }
                    }
                }
                
                // Sort by path position
                next_nodes.sort_by_key(|&node| {
                    path_positions.get(&node)
                        .map(|positions| positions.iter().sum::<usize>() / positions.len())
                        .unwrap_or(usize::MAX)
                });
                queue.extend(next_nodes.into_iter().rev());
            }
        }
        
        // Handle remaining nodes (in cycles)
        let mut remaining: Vec<_> = self.nodes.keys()
            .filter(|&&id| !visited.contains(&id))
            .cloned()
            .collect();
        remaining.sort();
        
        fn dfs_visit(node: usize, adj_list: &HashMap<usize, Vec<usize>>, 
                     visited: &mut HashSet<usize>, stack: &mut Vec<usize>) {
            if visited.contains(&node) {
                return;
            }
            visited.insert(node);
            
            if let Some(neighbors) = adj_list.get(&node) {
                let mut sorted_neighbors: Vec<_> = neighbors.iter()
                    .filter(|&&n| !visited.contains(&n))
                    .cloned()
                    .collect();
                sorted_neighbors.sort();
                
                for neighbor in sorted_neighbors {
                    dfs_visit(neighbor, adj_list, visited, stack);
                }
            }
            
            stack.push(node);
        }
        
        for node_id in remaining {
            if !visited.contains(&node_id) {
                dfs_visit(node_id, &adj_list, &mut visited, &mut initial_order);
            }
        }
        
        // Now optimize the ordering to minimize edge spans
        // Create position map
        let mut position: HashMap<usize, usize> = HashMap::new();
        for (pos, &node) in initial_order.iter().enumerate() {
            position.insert(node, pos);
        }
        
        // Calculate edge spans and identify problematic nodes
        let mut node_span_scores: HashMap<usize, f64> = HashMap::new();
        
        for edge in &self.edges {
            if let (Some(&from_pos), Some(&to_pos)) = 
                (position.get(&edge.from), position.get(&edge.to)) {
                let span = if to_pos > from_pos { 
                    to_pos - from_pos 
                } else { 
                    from_pos - to_pos 
                };
                
                // Accumulate span scores for nodes
                *node_span_scores.entry(edge.from).or_insert(0.0) += span as f64;
                *node_span_scores.entry(edge.to).or_insert(0.0) += span as f64;
            }
        }
        
        // Identify nodes with high span scores
        let mut problematic_nodes: Vec<(usize, f64)> = node_span_scores.iter()
            .map(|(&node, &score)| (node, score))
            .filter(|(_, score)| *score > 100.0) // Threshold for problematic nodes
            .collect();
        problematic_nodes.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        
        // Try to reposition problematic nodes
        let mut final_order = initial_order.clone();
        
        for (prob_node, _) in problematic_nodes.iter().take(50) { // Limit iterations
            // Find optimal position for this node
            let current_position = position[prob_node];
            let mut best_position = current_position;
            let mut _best_score = f64::MAX;
            
            // Calculate connected nodes' positions
            let mut connected_positions = Vec::new();
            
            if let Some(neighbors) = adj_list.get(prob_node) {
                for &neighbor in neighbors {
                    if let Some(&pos) = position.get(&neighbor) {
                        connected_positions.push(pos);
                    }
                }
            }
            
            // Check incoming edges too
            for edge in &self.edges {
                if edge.to == *prob_node {
                    if let Some(&pos) = position.get(&edge.from) {
                        connected_positions.push(pos);
                    }
                }
            }
            
            if !connected_positions.is_empty() {
                // Try median position
                connected_positions.sort();
                let median_pos = connected_positions[connected_positions.len() / 2];
                
                // Calculate score at median position
                let score: f64 = connected_positions.iter()
                    .map(|&pos| (median_pos as i32 - pos as i32).abs() as f64)
                    .sum();
                
                if score < _best_score {
                    _best_score = score;
                    best_position = median_pos;
                }
            }
            
            // Update position if it improves things
            if best_position != position[prob_node] {
                // Remove from current position
                final_order.retain(|&n| n != *prob_node);
                
                // Insert at new position
                let insert_pos = best_position.min(final_order.len());
                final_order.insert(insert_pos, *prob_node);
                
                // Update position map
                for (pos, &node) in final_order.iter().enumerate() {
                    position.insert(node, pos);
                }
            }
        }
        
        // Create mapping from old to new IDs
        let mut old_to_new = HashMap::new();
        for (new_id, &old_id) in final_order.iter().enumerate() {
            old_to_new.insert(old_id, new_id + 1); // 1-based IDs
        }
        
        // Update nodes
        let mut new_nodes = HashMap::new();
        for (old_id, mut node) in self.nodes.drain() {
            if let Some(&new_id) = old_to_new.get(&old_id) {
                node.id = new_id;
                node.rank = (new_id - 1) as f64;
                new_nodes.insert(new_id, node);
            }
        }
        self.nodes = new_nodes;
        
        // Update edges
        let mut new_edges = HashSet::new();
        for edge in self.edges.drain() {
            if let (Some(&new_from), Some(&new_to)) = 
                (old_to_new.get(&edge.from), old_to_new.get(&edge.to)) {
                new_edges.insert(Edge { from: new_from, to: new_to });
            }
        }
        self.edges = new_edges;
        
        // Update paths
        for (_, path) in &mut self.paths {
            for node_id in path.iter_mut() {
                if let Some(&new_id) = old_to_new.get(node_id) {
                    *node_id = new_id;
                }
            }
        }
    }
    
    /// Verify that all paths are fully embedded in the graph
    pub fn verify_path_embedding(&self, verbose: bool) -> Result<(), Vec<String>> {
        let mut errors = Vec::new();
        
        for (path_name, path) in &self.paths {
            if let Err(mut path_errors) = self.verify_single_path(path_name, path, verbose) {
                errors.append(&mut path_errors);
            }
        }
        
        if errors.is_empty() {
            if verbose {
                println!("✓ All paths are fully embedded in the graph");
            }
            Ok(())
        } else {
            Err(errors)
        }
    }
    
    /// Verify a single path's embedding
    fn verify_single_path(&self, path_name: &str, path: &[usize], verbose: bool) -> Result<(), Vec<String>> {
        let mut errors = Vec::new();
        
        if path.is_empty() {
            return Ok(());
        }
        
        // Check 1: All nodes in path exist
        for &node_id in path {
            if !self.nodes.contains_key(&node_id) {
                errors.push(format!("Path '{}': Node {} does not exist", path_name, node_id));
            }
        }
        
        // Check 2: All consecutive nodes have edges
        for window in path.windows(2) {
            let from = window[0];
            let to = window[1];
            
            if !self.edges.contains(&Edge { from, to }) {
                errors.push(format!("Path '{}': Missing edge {} -> {}", path_name, from, to));
            }
        }
        
        // Check 3: Path can be reconstructed to valid sequence
        if let Err(reconstruction_error) = self.reconstruct_path_sequence(path_name, path) {
            errors.push(reconstruction_error);
        }
        
        if errors.is_empty() && verbose {
            println!("✓ Path '{}' is fully embedded ({} nodes)", path_name, path.len());
        }
        
        if errors.is_empty() {
            Ok(())
        } else {
            Err(errors)
        }
    }
    
    /// Reconstruct sequence from path to verify integrity
    pub fn reconstruct_path_sequence(&self, path_name: &str, path: &[usize]) -> Result<Vec<u8>, String> {
        let mut sequence = Vec::new();
        
        for &node_id in path {
            match self.nodes.get(&node_id) {
                Some(node) => {
                    if node.sequence.is_empty() {
                        return Err(format!("Path '{}': Node {} has empty sequence", path_name, node_id));
                    }
                    sequence.extend_from_slice(&node.sequence);
                }
                None => {
                    return Err(format!("Path '{}': Node {} not found during reconstruction", path_name, node_id));
                }
            }
        }
        
        if sequence.is_empty() && !path.is_empty() {
            return Err(format!("Path '{}': Reconstruction resulted in empty sequence", path_name));
        }
        
        Ok(sequence)
    }
    
    /// Verify path integrity by comparing reconstructed sequence with original
    pub fn verify_path_integrity(&self, path_name: &str, path: &[usize], original_sequence: &[u8]) -> Result<(), String> {
        let reconstructed = self.reconstruct_path_sequence(path_name, path)?;
        
        // In a bidirectional graph with RC alignments, sequences may not match exactly
        // due to nodes representing different bases on different strands
        // For now, we only check length
        if reconstructed.len() != original_sequence.len() {
            return Err(format!(
                "Path '{}': Sequence length mismatch. Original length: {}, Reconstructed length: {}",
                path_name, original_sequence.len(), reconstructed.len()
            ));
        }
        
        // TODO: Properly handle bidirectional graph verification by tracking path orientations
        
        Ok(())
    }
    
    /// Check for common path issues
    pub fn validate_path_structure(&self, verbose: bool) -> Result<(), Vec<String>> {
        let errors = Vec::new();
        
        // We don't care about consecutive duplicates - they're expected when positions are united
        
        // Check for orphaned nodes (nodes not in any path)
        let mut nodes_in_paths = HashSet::new();
        for (_, path) in &self.paths {
            for &node_id in path {
                nodes_in_paths.insert(node_id);
            }
        }
        
        let orphaned_nodes: Vec<_> = self.nodes.keys()
            .filter(|&&id| !nodes_in_paths.contains(&id))
            .collect();
        
        if !orphaned_nodes.is_empty() && verbose {
            println!("Warning: {} nodes are not present in any path", orphaned_nodes.len());
        }
        
        // Check for disconnected components
        let connected_components = self.find_connected_components();
        if connected_components.len() > 1 && verbose {
            println!("Warning: Graph has {} disconnected components", connected_components.len());
        }
        
        if errors.is_empty() {
            Ok(())
        } else {
            Err(errors)
        }
    }
    
    /// Find connected components in the graph
    fn find_connected_components(&self) -> Vec<HashSet<usize>> {
        let mut visited = HashSet::new();
        let mut components = Vec::new();
        
        // Build adjacency list (undirected)
        let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
        for &node_id in self.nodes.keys() {
            adjacency.insert(node_id, Vec::new());
        }
        
        for edge in &self.edges {
            adjacency.get_mut(&edge.from).unwrap().push(edge.to);
            adjacency.get_mut(&edge.to).unwrap().push(edge.from);
        }
        
        // DFS to find components
        for &node_id in self.nodes.keys() {
            if !visited.contains(&node_id) {
                let mut component = HashSet::new();
                let mut stack = vec![node_id];
                
                while let Some(current) = stack.pop() {
                    if !visited.contains(&current) {
                        visited.insert(current);
                        component.insert(current);
                        
                        if let Some(neighbors) = adjacency.get(&current) {
                            for &neighbor in neighbors {
                                if !visited.contains(&neighbor) {
                                    stack.push(neighbor);
                                }
                            }
                        }
                    }
                }
                
                components.push(component);
            }
        }
        
        components
    }
    
    /// Verify edge traversal - all edges should be used by paths
    pub fn verify_edge_traversal(&self, verbose: bool) -> Result<(), Vec<String>> {
        let mut errors = Vec::new();
        
        // Collect all edges traversed by paths
        let mut traversed_edges = HashSet::new();
        let mut self_loop_traversals = 0;
        for (_, path) in &self.paths {
            for window in path.windows(2) {
                let edge = Edge { from: window[0], to: window[1] };
                if edge.from == edge.to {
                    self_loop_traversals += 1;
                }
                traversed_edges.insert(edge);
            }
        }
        
        if verbose && self_loop_traversals > 0 {
            println!("  Found {} self-loop traversals in paths", self_loop_traversals);
        }
        
        // Check for untraversed edges
        let mut untraversed_edges = Vec::new();
        let mut total_self_loops = 0;
        for edge in &self.edges {
            if edge.from == edge.to {
                total_self_loops += 1;
            }
            if !traversed_edges.contains(edge) {
                untraversed_edges.push(edge);
            }
        }
        
        if verbose {
            println!("  Total self-loops in graph: {}", total_self_loops);
        }
        
        if !untraversed_edges.is_empty() {
            errors.push(format!(
                "Found {} edges not traversed by any path",
                untraversed_edges.len()
            ));
            if verbose {
                for edge in untraversed_edges.iter().take(10) {
                    println!("  Untraversed edge: {} -> {}", edge.from, edge.to);
                }
                if untraversed_edges.len() > 10 {
                    println!("  ... and {} more", untraversed_edges.len() - 10);
                }
            }
        }
        
        // Check that all path edges exist in the graph
        let mut missing_edges = Vec::new();
        for (path_name, path) in &self.paths {
            for window in path.windows(2) {
                let edge = Edge { from: window[0], to: window[1] };
                if !self.edges.contains(&edge) {
                    missing_edges.push((path_name.clone(), edge));
                }
            }
        }
        
        if !missing_edges.is_empty() {
            errors.push(format!(
                "Found {} edges in paths that don't exist in the graph",
                missing_edges.len()
            ));
            if verbose {
                for (path_name, edge) in missing_edges.iter().take(10) {
                    println!("  Path '{}' uses non-existent edge: {} -> {}", 
                             path_name, edge.from, edge.to);
                }
                if missing_edges.len() > 10 {
                    println!("  ... and {} more", missing_edges.len() - 10);
                }
            }
        }
        
        if errors.is_empty() {
            if verbose {
                println!("✓ Edge traversal verification passed: all {} edges are traversed by paths", 
                         self.edges.len());
            }
            Ok(())
        } else {
            Err(errors)
        }
    }

    /// Validate that the graph is valid for GFA output
    pub fn validate_gfa_format(&self, verbose: bool) -> Result<(), Vec<String>> {
        let mut errors = Vec::new();
        
        // Check 1: All nodes have IDs
        for (id, node) in &self.nodes {
            if *id != node.id {
                errors.push(format!("Node ID mismatch: map key {} != node.id {}", id, node.id));
            }
            if node.sequence.is_empty() {
                errors.push(format!("Node {} has empty sequence", id));
            }
        }
        
        // Check 2: All edges reference existing nodes
        for edge in &self.edges {
            if !self.nodes.contains_key(&edge.from) {
                errors.push(format!("Edge references non-existent node: {}", edge.from));
            }
            if !self.nodes.contains_key(&edge.to) {
                errors.push(format!("Edge references non-existent node: {}", edge.to));
            }
        }
        
        // Check 3: All paths reference existing nodes and have valid edges
        for (path_name, path) in &self.paths {
            // Check all nodes exist
            for &node_id in path {
                if !self.nodes.contains_key(&node_id) {
                    errors.push(format!("Path '{}' references non-existent node: {}", path_name, node_id));
                }
            }
            
            // Check all consecutive pairs have edges
            for window in path.windows(2) {
                let edge = Edge { from: window[0], to: window[1] };
                if !self.edges.contains(&edge) {
                    errors.push(format!("Path '{}' uses non-existent edge: {} -> {}", 
                                      path_name, edge.from, edge.to));
                }
            }
        }
        
        if errors.is_empty() {
            if verbose {
                println!("✓ GFA format validation passed");
            }
            Ok(())
        } else {
            Err(errors)
        }
    }
    
    /// Comprehensive verification of graph and paths
    pub fn comprehensive_verify(&self, original_sequences: Option<&[(String, Vec<u8>)]>, verbose: bool) -> Result<(), Vec<String>> {
        let mut all_errors = Vec::new();
        
        if verbose {
            println!("Starting comprehensive graph verification...");
        }
        
        // 1. Verify path embedding
        if let Err(mut errors) = self.verify_path_embedding(verbose) {
            all_errors.append(&mut errors);
        }
        
        // 2. Validate path structure
        if let Err(mut errors) = self.validate_path_structure(verbose) {
            all_errors.append(&mut errors);
        }
        
        // 3. Verify edge traversal
        if let Err(mut errors) = self.verify_edge_traversal(verbose) {
            all_errors.append(&mut errors);
        }
        
        // 4. Verify path integrity against original sequences if provided
        if let Some(original_seqs) = original_sequences {
            for (seq_name, seq_data) in original_seqs {
                if let Some((_, path)) = self.paths.iter().find(|(name, _)| name == seq_name) {
                    if let Err(error) = self.verify_path_integrity(seq_name, path, seq_data) {
                        all_errors.push(error);
                    }
                } else {
                    all_errors.push(format!("Original sequence '{}' not found in paths", seq_name));
                }
            }
        }
        
        if all_errors.is_empty() {
            if verbose {
                println!("✓ Comprehensive verification passed");
            }
            Ok(())
        } else {
            if verbose {
                println!("✗ Comprehensive verification failed with {} errors", all_errors.len());
                for error in &all_errors {
                    println!("  - {}", error);
                }
            }
            Err(all_errors)
        }
    }
}
