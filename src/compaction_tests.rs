#[cfg(test)]
mod compaction_tests {
    use crate::graph_ops::{Graph, Node, Edge};
    use std::collections::{HashMap, HashSet};
    
    fn create_test_graph() -> Graph {
        let mut graph = Graph::new();
        
        // Create a simple linear chain: 1 -> 2 -> 3 -> 4
        for i in 1..=4 {
            graph.nodes.insert(i, Node {
                id: i,
                sequence: vec![b'A' + (i as u8 - 1)],
                rank: i as f64,
            });
        }
        
        // Add edges
        graph.edges.insert(Edge { from: 1, to: 2 });
        graph.edges.insert(Edge { from: 2, to: 3 });
        graph.edges.insert(Edge { from: 3, to: 4 });
        
        // Add paths
        graph.paths.push(("path1".to_string(), vec![1, 2, 3, 4]));
        
        graph
    }
    
    #[test]
    fn test_basic_compaction() {
        let mut graph = create_test_graph();
        
        // Before compaction
        assert_eq!(graph.nodes.len(), 4);
        assert_eq!(graph.edges.len(), 3);
        
        let compacted = graph.compact_nodes();
        
        // After compaction - all 4 nodes should merge into 1
        assert_eq!(compacted, 3); // 3 nodes were merged away
        assert_eq!(graph.nodes.len(), 1);
        assert_eq!(graph.edges.len(), 0);
        
        // Check the merged node
        let node = graph.nodes.values().next().unwrap();
        assert_eq!(node.sequence, vec![b'A', b'B', b'C', b'D']);
        
        // Check path is updated correctly
        assert_eq!(graph.paths[0].1, vec![1]); // All nodes merged into node 1
    }
    
    #[test]
    fn test_branching_graph() {
        let mut graph = Graph::new();
        
        // Create a branching graph:
        //     2 -> 4
        //   /      \
        // 1          6
        //   \      /
        //     3 -> 5
        
        for i in 1..=6 {
            graph.nodes.insert(i, Node {
                id: i,
                sequence: vec![b'A' + (i as u8 - 1)],
                rank: i as f64,
            });
        }
        
        graph.edges.insert(Edge { from: 1, to: 2 });
        graph.edges.insert(Edge { from: 1, to: 3 });
        graph.edges.insert(Edge { from: 2, to: 4 });
        graph.edges.insert(Edge { from: 3, to: 5 });
        graph.edges.insert(Edge { from: 4, to: 6 });
        graph.edges.insert(Edge { from: 5, to: 6 });
        
        graph.paths.push(("path1".to_string(), vec![1, 2, 4, 6]));
        graph.paths.push(("path2".to_string(), vec![1, 3, 5, 6]));
        
        let initial_nodes = graph.nodes.len();
        let compacted = graph.compact_nodes();
        
        // Two chains can be compacted: 2->4 and 3->5
        assert_eq!(compacted, 2); // 2 nodes merged away
        assert_eq!(graph.nodes.len(), 4); // 6 - 2 = 4 nodes remain
    }
    
    #[test]
    fn test_partial_chain() {
        let mut graph = Graph::new();
        
        // Create a graph with a partial chain:
        // 1 -> 2 -> 3 -> 4 -> 5
        //      ^              |
        //      |              v
        //      7 -----------> 6
        
        for i in 1..=7 {
            graph.nodes.insert(i, Node {
                id: i,
                sequence: vec![b'A' + (i as u8 - 1)],
                rank: i as f64,
            });
        }
        
        graph.edges.insert(Edge { from: 1, to: 2 });
        graph.edges.insert(Edge { from: 2, to: 3 });
        graph.edges.insert(Edge { from: 3, to: 4 });
        graph.edges.insert(Edge { from: 4, to: 5 });
        graph.edges.insert(Edge { from: 5, to: 6 });
        graph.edges.insert(Edge { from: 7, to: 2 });
        graph.edges.insert(Edge { from: 7, to: 6 });
        
        graph.paths.push(("path1".to_string(), vec![1, 2, 3, 4, 5, 6]));
        graph.paths.push(("path2".to_string(), vec![7, 2, 3, 4, 5, 6]));
        graph.paths.push(("path3".to_string(), vec![7, 6]));
        
        let compacted = graph.compact_nodes();
        
        // 2->3->4->5 can be compacted (linear chain)
        assert_eq!(compacted, 3); // 3 nodes merged away
        assert_eq!(graph.nodes.len(), 4); // 7 - 3 = 4 nodes remain
        
        // Check that the compacted node has the right sequence (BCDE from nodes 2,3,4,5)
        let has_merged_sequence = graph.nodes.values()
            .any(|n| n.sequence == vec![b'B', b'C', b'D', b'E']);
        assert!(has_merged_sequence);
    }
    
    #[test]
    fn test_self_loop() {
        let mut graph = Graph::new();
        
        // Create a graph with a self-loop
        graph.nodes.insert(1, Node {
            id: 1,
            sequence: vec![b'A'],
            rank: 1.0,
        });
        graph.nodes.insert(2, Node {
            id: 2,
            sequence: vec![b'B'],
            rank: 2.0,
        });
        
        graph.edges.insert(Edge { from: 1, to: 2 });
        graph.edges.insert(Edge { from: 2, to: 2 }); // Self-loop
        
        graph.paths.push(("path1".to_string(), vec![1, 2, 2, 2]));
        
        let compacted = graph.compact_nodes();
        
        // No compaction should occur due to self-loop
        assert_eq!(compacted, 0);
        assert_eq!(graph.nodes.len(), 2);
    }
    
    #[test]
    fn test_complex_paths() {
        let mut graph = Graph::new();
        
        // Create a linear chain but with complex path usage
        for i in 1..=5 {
            graph.nodes.insert(i, Node {
                id: i,
                sequence: vec![b'A' + (i as u8 - 1)],
                rank: i as f64,
            });
        }
        
        graph.edges.insert(Edge { from: 1, to: 2 });
        graph.edges.insert(Edge { from: 2, to: 3 });
        graph.edges.insert(Edge { from: 3, to: 4 });
        graph.edges.insert(Edge { from: 4, to: 5 });
        
        // Paths that use the chain in different ways
        graph.paths.push(("full".to_string(), vec![1, 2, 3, 4, 5]));
        graph.paths.push(("skip".to_string(), vec![1, 2, 4, 5])); // Invalid path that skips 3
        
        let compacted = graph.compact_nodes();
        
        // The skip path prevents full compaction
        // Only nodes that are consistently used can be merged
        assert!(graph.nodes.len() < 5); // Some compaction should occur
    }
    
    #[test]
    fn test_empty_graph() {
        let mut graph = Graph::new();
        let compacted = graph.compact_nodes();
        
        assert_eq!(compacted, 0);
        assert_eq!(graph.nodes.len(), 0);
        assert_eq!(graph.edges.len(), 0);
    }
    
    #[test]
    fn test_single_node() {
        let mut graph = Graph::new();
        
        graph.nodes.insert(1, Node {
            id: 1,
            sequence: vec![b'A'],
            rank: 1.0,
        });
        
        graph.paths.push(("path1".to_string(), vec![1]));
        
        let compacted = graph.compact_nodes();
        
        assert_eq!(compacted, 0);
        assert_eq!(graph.nodes.len(), 1);
    }
    
    #[test]
    fn test_preserve_sequences() {
        let mut graph = Graph::new();
        
        // Create chain with specific sequences
        graph.nodes.insert(1, Node {
            id: 1,
            sequence: vec![b'A', b'T', b'C'],
            rank: 1.0,
        });
        graph.nodes.insert(2, Node {
            id: 2,
            sequence: vec![b'G', b'A', b'T'],
            rank: 2.0,
        });
        graph.nodes.insert(3, Node {
            id: 3,
            sequence: vec![b'C', b'G'],
            rank: 3.0,
        });
        
        graph.edges.insert(Edge { from: 1, to: 2 });
        graph.edges.insert(Edge { from: 2, to: 3 });
        
        graph.paths.push(("path1".to_string(), vec![1, 2, 3]));
        
        // Store original sequence
        let original_seq: Vec<u8> = vec![1, 2, 3].iter()
            .flat_map(|&id| graph.nodes[&id].sequence.clone())
            .collect();
        
        let compacted = graph.compact_nodes();
        
        assert_eq!(compacted, 2); // nodes 2 and 3 merged into node 1
        assert_eq!(graph.nodes.len(), 1);
        
        // Check merged sequence
        let merged_node = graph.nodes.values().next().unwrap();
        assert_eq!(merged_node.sequence, original_seq);
        assert_eq!(merged_node.sequence, vec![b'A', b'T', b'C', b'G', b'A', b'T', b'C', b'G']);
    }
}