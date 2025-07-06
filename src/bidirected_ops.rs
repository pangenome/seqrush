use crate::bidirected_graph::{Handle, BiNode, BiPath, BiEdge};
use crate::graph_ops::{Node, Edge, Graph};
use std::collections::{HashMap, HashSet, VecDeque, BTreeSet};
use bitvec::prelude::*;

/// A bidirected graph that extends the basic Graph with orientation support
pub struct BidirectedGraph {
    pub nodes: HashMap<usize, BiNode>,
    pub edges: HashSet<BiEdge>,
    pub paths: Vec<BiPath>,
}

impl BidirectedGraph {
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
            let bi_edge = BiEdge::new(
                Handle::forward(edge.from),
                Handle::forward(edge.to),
            );
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
            let path_nodes: Vec<usize> = bi_path.steps.iter()
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
        self.nodes.get(&handle.node_id())
            .map(|node| node.get_sequence(handle.is_reverse()))
    }
    
    /// Check if an edge exists
    pub fn has_edge(&self, from: Handle, to: Handle) -> bool {
        self.edges.contains(&BiEdge::new(from, to))
    }
    
    /// Get all edges from a handle
    pub fn edges_from(&self, handle: Handle) -> Vec<Handle> {
        self.edges.iter()
            .filter(|edge| edge.from == handle)
            .map(|edge| edge.to)
            .collect()
    }
    
    /// Get all edges to a handle
    pub fn edges_to(&self, handle: Handle) -> Vec<Handle> {
        self.edges.iter()
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
        let mut written_edges = HashSet::new();
        for edge in &self.edges {
            // Get canonical form to avoid duplicates
            let canonical = edge.canonical();
            if written_edges.insert(canonical) {
                writeln!(writer, "L\t{}\t{}\t{}\t{}\t0M",
                    canonical.from.node_id(),
                    canonical.from.orientation_char(),
                    canonical.to.node_id(),
                    canonical.to.orientation_char()
                )?;
            }
        }
        
        // Write paths with orientations
        for path in &self.paths {
            write!(writer, "P\t{}\t", path.name)?;
            
            // Write path steps with orientations
            let steps: Vec<String> = path.steps.iter()
                .map(|handle| format!("{}{}", handle.node_id(), handle.orientation_char()))
                .collect();
            write!(writer, "{}\t*", steps.join(","))?;
            
            writeln!(writer)?;
        }
        
        Ok(())
    }
    
    /// Perform bidirected topological sort using modified Kahn's algorithm
    /// This handles cycles and bidirected edges properly
    pub fn topological_sort(&mut self, use_heads: bool, use_tails: bool, verbose: bool) -> Vec<Handle> {
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
            edges_by_from.entry(edge.from).or_insert_with(Vec::new).push((edge.to, edge_idx));
            edges_by_to.entry(edge.to).or_insert_with(Vec::new).push((edge.from, edge_idx));
            
            // Implied reverse edge: to.flip() → from.flip()
            edges_by_from.entry(edge.to.flip()).or_insert_with(Vec::new).push((edge.from.flip(), edge_idx));
            edges_by_to.entry(edge.from.flip()).or_insert_with(Vec::new).push((edge.to.flip(), edge_idx));
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
                    if !visited[idx] && edges_by_from.get(&handle).map(|v| v.is_empty()).unwrap_or(true) {
                        queue.push_back(handle);
                    }
                }
            }
        };
        
        // Start with initial seeds - but prioritize forward orientations
        add_seeds(&mut queue, &in_degree, &visited, &handle_to_idx, use_heads, use_tails);
        
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
            println!("Topological sort completed: {} nodes ordered", node_ordering.len());
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
            if let (Some(&new_from), Some(&new_to)) = 
                (old_to_new.get(&edge.from.node_id()), old_to_new.get(&edge.to.node_id())) {
                let new_edge = BiEdge::new(
                    Handle::new(new_from, edge.from.is_reverse()),
                    Handle::new(new_to, edge.to.is_reverse())
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
        assert_eq!(graph.get_sequence(Handle::forward(1)), Some(b"ATCG".to_vec()));
        
        // Reverse complement sequence
        assert_eq!(graph.get_sequence(Handle::reverse(1)), Some(b"CGAT".to_vec()));
    }
    
    #[test]
    fn test_path_with_orientations() {
        let mut graph = BidirectedGraph::new();
        graph.add_node(1, b"ATG".to_vec());
        graph.add_node(2, b"CGA".to_vec());
        graph.add_node(3, b"TAC".to_vec());
        
        // Build a path that traverses nodes in different orientations
        graph.build_path("test_path".to_string(), vec![
            (1, false), // Forward
            (2, true),  // Reverse
            (3, false), // Forward
        ]);
        
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