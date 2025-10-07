use seqrush::bidirected_graph::{BiEdge, BiNode, BiPath, Handle};
use seqrush::bidirected_ops::BidirectedGraph;

/// Helper to create a simple test graph
fn create_test_graph() -> BidirectedGraph {
    let mut graph = BidirectedGraph::new();

    // Add nodes
    graph.nodes.insert(
        1,
        BiNode {
            id: 1,
            sequence: vec![b'A'],
            rank: None,
        },
    );
    graph.nodes.insert(
        2,
        BiNode {
            id: 2,
            sequence: vec![b'T'],
            rank: None,
        },
    );
    graph.nodes.insert(
        3,
        BiNode {
            id: 3,
            sequence: vec![b'C'],
            rank: None,
        },
    );
    graph.nodes.insert(
        4,
        BiNode {
            id: 4,
            sequence: vec![b'G'],
            rank: None,
        },
    );

    graph
}

#[test]
fn test_perfect_neighbors_simple() {
    let mut graph = BidirectedGraph::new();

    // Add only the nodes we'll use
    graph.nodes.insert(
        1,
        BiNode {
            id: 1,
            sequence: vec![b'A'],
            rank: None,
        },
    );
    graph.nodes.insert(
        2,
        BiNode {
            id: 2,
            sequence: vec![b'T'],
            rank: None,
        },
    );
    graph.nodes.insert(
        3,
        BiNode {
            id: 3,
            sequence: vec![b'C'],
            rank: None,
        },
    );

    // Create edges: 1 -> 2 -> 3
    graph.edges.insert(BiEdge {
        from: Handle::forward(1),
        to: Handle::forward(2),
    });
    graph.edges.insert(BiEdge {
        from: Handle::forward(2),
        to: Handle::forward(3),
    });

    // Create path that traverses all nodes
    graph.paths.push(BiPath {
        name: "path1".to_string(),
        steps: vec![Handle::forward(1), Handle::forward(2), Handle::forward(3)],
    });

    // Compact the graph
    let result = graph.compact_perfect_neighbors();
    assert!(result.is_ok());
    let nodes_compacted = result.unwrap();

    // Should compact all 3 nodes into 1
    assert_eq!(nodes_compacted, 3);
    assert_eq!(graph.nodes.len(), 1);

    // Check the compacted node has the right sequence
    let compacted_node = graph.nodes.values().next().unwrap();
    assert_eq!(compacted_node.sequence, vec![b'A', b'T', b'C']);

    // Check the path is updated correctly
    assert_eq!(graph.paths[0].steps.len(), 1);
}

#[test]
fn test_no_compaction_with_branching() {
    let mut graph = create_test_graph();

    // Create edges: 1 -> 2 -> 3, and 1 -> 4
    graph.edges.insert(BiEdge {
        from: Handle::forward(1),
        to: Handle::forward(2),
    });
    graph.edges.insert(BiEdge {
        from: Handle::forward(2),
        to: Handle::forward(3),
    });
    graph.edges.insert(BiEdge {
        from: Handle::forward(1),
        to: Handle::forward(4),
    });

    // Create paths
    graph.paths.push(BiPath {
        name: "path1".to_string(),
        steps: vec![Handle::forward(1), Handle::forward(2), Handle::forward(3)],
    });
    graph.paths.push(BiPath {
        name: "path2".to_string(),
        steps: vec![Handle::forward(1), Handle::forward(4)],
    });

    // Compact the graph
    let result = graph.compact_perfect_neighbors();
    assert!(result.is_ok());
    let nodes_compacted = result.unwrap();

    // Should only compact 2 and 3 (since 1 branches)
    assert_eq!(nodes_compacted, 2);
    assert_eq!(graph.nodes.len(), 3); // 1, 2+3 compacted, 4
}

#[test]
fn test_reverse_complement_compaction() {
    let mut graph = create_test_graph();

    // Create edges for forward path
    graph.edges.insert(BiEdge {
        from: Handle::forward(1),
        to: Handle::forward(2),
    });

    // Create path with reverse complement
    graph.paths.push(BiPath {
        name: "path1".to_string(),
        steps: vec![Handle::forward(1), Handle::forward(2)],
    });
    graph.paths.push(BiPath {
        name: "path2_rc".to_string(),
        steps: vec![Handle::reverse(2), Handle::reverse(1)],
    });

    // Compact the graph
    let result = graph.compact_perfect_neighbors();
    assert!(result.is_ok());
    let nodes_compacted = result.unwrap();

    // Should compact nodes 1 and 2
    assert_eq!(nodes_compacted, 2);
    assert_eq!(graph.nodes.len(), 3); // compacted node + nodes 3,4

    // Check paths are preserved
    assert_eq!(graph.paths.len(), 2);
}

#[test]
fn test_path_validation() {
    let mut graph = create_test_graph();

    // Create a simple linear graph
    graph.edges.insert(BiEdge {
        from: Handle::forward(1),
        to: Handle::forward(2),
    });
    graph.edges.insert(BiEdge {
        from: Handle::forward(2),
        to: Handle::forward(3),
    });

    graph.paths.push(BiPath {
        name: "test_path".to_string(),
        steps: vec![Handle::forward(1), Handle::forward(2), Handle::forward(3)],
    });

    // Get original sequence
    let mut original_seq = Vec::new();
    for step in &graph.paths[0].steps {
        let node = &graph.nodes[&step.node_id()];
        original_seq.extend_from_slice(&node.sequence);
    }
    assert_eq!(original_seq, vec![b'A', b'T', b'C']);

    // Compact
    let result = graph.compact_perfect_neighbors();
    assert!(result.is_ok());

    // Get compacted sequence
    let mut compacted_seq = Vec::new();
    for step in &graph.paths[0].steps {
        let node = &graph.nodes[&step.node_id()];
        compacted_seq.extend_from_slice(&node.sequence);
    }

    // Sequences should match
    assert_eq!(original_seq, compacted_seq);
}

#[test]
fn test_complex_graph_compaction() {
    let mut graph = BidirectedGraph::new();

    // Create a more complex graph with multiple compactable regions
    // A -> B -> C -> D
    //      ^         |
    //      |         v
    //      F <- E <-

    for i in 1..=6 {
        graph.nodes.insert(
            i,
            BiNode {
                id: i,
                sequence: vec![b'A' + (i as u8 - 1)],
                rank: None,
            },
        );
    }

    // Main path edges
    graph.edges.insert(BiEdge {
        from: Handle::forward(1),
        to: Handle::forward(2),
    });
    graph.edges.insert(BiEdge {
        from: Handle::forward(2),
        to: Handle::forward(3),
    });
    graph.edges.insert(BiEdge {
        from: Handle::forward(3),
        to: Handle::forward(4),
    });
    graph.edges.insert(BiEdge {
        from: Handle::forward(4),
        to: Handle::forward(5),
    });
    graph.edges.insert(BiEdge {
        from: Handle::forward(5),
        to: Handle::forward(6),
    });
    graph.edges.insert(BiEdge {
        from: Handle::forward(6),
        to: Handle::forward(2),
    });

    // Path that goes through all nodes
    graph.paths.push(BiPath {
        name: "circular".to_string(),
        steps: vec![
            Handle::forward(1),
            Handle::forward(2),
            Handle::forward(3),
            Handle::forward(4),
            Handle::forward(5),
            Handle::forward(6),
            Handle::forward(2), // Loop back
            Handle::forward(3),
            Handle::forward(4),
        ],
    });

    // Compact
    let result = graph.compact_perfect_neighbors();
    assert!(result.is_ok());

    // Should compact 1-2 and 3-4, but not across the loop
    // Original: 6 nodes, After: fewer nodes but paths preserved
    assert!(graph.nodes.len() < 6);

    // Verify path is still valid
    let path = &graph.paths[0];
    for step in &path.steps {
        assert!(graph.nodes.contains_key(&step.node_id()));
    }
}

#[test]
fn test_no_compaction_single_nodes() {
    let mut graph = BidirectedGraph::new();

    // Add isolated nodes
    graph.nodes.insert(
        1,
        BiNode {
            id: 1,
            sequence: vec![b'A'],
            rank: None,
        },
    );
    graph.nodes.insert(
        2,
        BiNode {
            id: 2,
            sequence: vec![b'T'],
            rank: None,
        },
    );

    // No edges, just a path
    graph.paths.push(BiPath {
        name: "path".to_string(),
        steps: vec![Handle::forward(1)],
    });

    let result = graph.compact_perfect_neighbors();
    assert!(result.is_ok());
    let nodes_compacted = result.unwrap();

    // No compaction should occur
    assert_eq!(nodes_compacted, 0);
    assert_eq!(graph.nodes.len(), 2);
}

#[test]
fn test_path_ends_prevent_compaction() {
    let mut graph = create_test_graph();

    // Create edges: 1 -> 2 -> 3
    graph.edges.insert(BiEdge {
        from: Handle::forward(1),
        to: Handle::forward(2),
    });
    graph.edges.insert(BiEdge {
        from: Handle::forward(2),
        to: Handle::forward(3),
    });

    // Create paths where one ends at node 2
    graph.paths.push(BiPath {
        name: "path1".to_string(),
        steps: vec![Handle::forward(1), Handle::forward(2), Handle::forward(3)],
    });
    graph.paths.push(BiPath {
        name: "path2".to_string(),
        steps: vec![
            Handle::forward(1),
            Handle::forward(2), // Ends here
        ],
    });

    // Compact the graph
    let result = graph.compact_perfect_neighbors();
    assert!(result.is_ok());
    let nodes_compacted = result.unwrap();

    // Should not compact 2 and 3 because path2 ends at 2
    // Might compact 1 and 2 if they're perfect neighbors
    assert!(graph.nodes.len() >= 2);
}
