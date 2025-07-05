use std::collections::{HashMap, HashSet, VecDeque};

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
pub struct Graph {
    pub nodes: HashMap<usize, Node>,
    pub edges: HashSet<Edge>,
    pub paths: Vec<(String, Vec<usize>)>,
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
        let mut compacted_count = 0;
        let mut processed = HashSet::new();
        
        // Find all linear chains (simple components)
        let chains = self.find_linear_chains();
        
        // Add a safety check to prevent infinite loops
        let max_iterations = self.nodes.len();
        let mut iterations = 0;
        
        for chain in chains {
            if iterations >= max_iterations {
                eprintln!("WARNING: Compaction stopped after {} iterations", iterations);
                break;
            }
            iterations += 1;
            
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
    fn find_linear_chains(&self) -> Vec<Vec<usize>> {
        let mut chains = Vec::new();
        let mut visited = HashSet::new();
        
        // Build adjacency lists
        let mut forward_edges: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut backward_edges: HashMap<usize, Vec<usize>> = HashMap::new();
        
        for edge in &self.edges {
            forward_edges.entry(edge.from).or_insert_with(Vec::new).push(edge.to);
            backward_edges.entry(edge.to).or_insert_with(Vec::new).push(edge.from);
        }
        
        // Find nodes that can be part of linear chains
        // A node can be in a linear chain if it has exactly 1 in-edge and 1 out-edge
        for &node_id in self.nodes.keys() {
            if visited.contains(&node_id) {
                continue;
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
    fn merge_chain(&mut self, chain: &[usize]) {
        if chain.len() < 2 {
            return;
        }
        
        // Create new node with concatenated sequence
        let new_id = chain[0]; // Use first node's ID
        let mut new_sequence = Vec::new();
        let mut total_rank = 0.0;
        
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
            
            // Skip self-loops unless they were already present
            if from != to || (edge.from == edge.to && chain.contains(&edge.from)) {
                new_edges.insert(Edge { from, to });
            }
        }
        
        // Remove internal edges of the chain
        for i in 0..chain.len() - 1 {
            new_edges.remove(&Edge { from: chain[i], to: chain[i + 1] });
        }
        
        self.edges = new_edges;
        
        // Update paths
        for (_, path) in &mut self.paths {
            let mut new_path = Vec::new();
            let mut i = 0;
            
            while i < path.len() {
                let node_id = path[i];
                
                if node_id == chain[0] {
                    // Check if this starts the chain in the path
                    let mut matches_chain = true;
                    for (j, &chain_node) in chain.iter().enumerate() {
                        if i + j >= path.len() || path[i + j] != chain_node {
                            matches_chain = false;
                            break;
                        }
                    }
                    
                    if matches_chain {
                        // Replace the chain with the merged node
                        new_path.push(new_id);
                        i += chain.len();
                    } else {
                        new_path.push(new_id);
                        i += 1;
                    }
                } else if chain.contains(&node_id) {
                    // This node was part of a chain but not at the expected position
                    // This shouldn't happen in a well-formed path
                    new_path.push(new_id);
                    i += 1;
                } else {
                    new_path.push(node_id);
                    i += 1;
                }
            }
            
            *path = new_path;
        }
    }
    
    /// Perform topological sort on the graph
    pub fn topological_sort(&mut self) {
        // Build adjacency lists
        let mut in_degree: HashMap<usize, usize> = HashMap::new();
        let mut out_edges: HashMap<usize, Vec<usize>> = HashMap::new();
        
        // Initialize all nodes with 0 in-degree
        for &node_id in self.nodes.keys() {
            in_degree.insert(node_id, 0);
            out_edges.insert(node_id, Vec::new());
        }
        
        // Count in-degrees and build adjacency list
        for edge in &self.edges {
            *in_degree.get_mut(&edge.to).unwrap() += 1;
            out_edges.get_mut(&edge.from).unwrap().push(edge.to);
        }
        
        // Find all nodes with no incoming edges
        let mut queue: VecDeque<usize> = VecDeque::new();
        for (&node_id, &degree) in &in_degree {
            if degree == 0 {
                queue.push_back(node_id);
            }
        }
        
        // Process nodes in topological order
        let mut rank = 0.0;
        let mut sorted_ids = Vec::new();
        
        while let Some(node_id) = queue.pop_front() {
            sorted_ids.push(node_id);
            
            // Assign rank
            if let Some(node) = self.nodes.get_mut(&node_id) {
                node.rank = rank;
                rank += 1.0;
            }
            
            // Process all neighbors
            if let Some(neighbors) = out_edges.get(&node_id) {
                for &neighbor in neighbors {
                    if let Some(degree) = in_degree.get_mut(&neighbor) {
                        *degree -= 1;
                        if *degree == 0 {
                            queue.push_back(neighbor);
                        }
                    }
                }
            }
        }
        
        // Handle cycles by assigning remaining nodes arbitrary ranks
        let mut remaining: Vec<usize> = self.nodes.keys()
            .filter(|&&id| !sorted_ids.contains(&id))
            .cloned()
            .collect();
        
        // Sort remaining by current rank to maintain some order
        remaining.sort_by(|&a, &b| {
            let rank_a = self.nodes.get(&a).map(|n| n.rank).unwrap_or(0.0);
            let rank_b = self.nodes.get(&b).map(|n| n.rank).unwrap_or(0.0);
            rank_a.partial_cmp(&rank_b).unwrap()
        });
        
        for node_id in remaining {
            if let Some(node) = self.nodes.get_mut(&node_id) {
                node.rank = rank;
                rank += 1.0;
            }
            sorted_ids.push(node_id);
        }
        
        // Renumber nodes based on topological order
        let mut old_to_new: HashMap<usize, usize> = HashMap::new();
        for (new_id, &old_id) in sorted_ids.iter().enumerate() {
            old_to_new.insert(old_id, new_id + 1); // 1-based IDs
        }
        
        // Create new nodes map with updated IDs
        let mut new_nodes = HashMap::new();
        for (old_id, mut node) in self.nodes.drain() {
            if let Some(&new_id) = old_to_new.get(&old_id) {
                node.id = new_id;
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
}