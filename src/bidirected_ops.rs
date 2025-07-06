use crate::bidirected_graph::{Handle, BiNode, BiPath, BiEdge};
use crate::graph_ops::{Node, Edge, Graph};
use std::collections::{HashMap, HashSet};

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