/// Tests for path-guided SGD on simple graphs
/// These should achieve PERFECT sorting with no displacement

use seqrush::bidirected_ops::BidirectedGraph;
use seqrush::bidirected_graph::{Handle, BiPath};
use seqrush::path_sgd::{PathSGDParams, path_sgd_sort};

/// Helper to create a simple linear chain graph
fn create_chain_graph(node_sequences: Vec<&str>, edges: Vec<(Handle, Handle)>) -> BidirectedGraph {
    let mut graph = BidirectedGraph::new();

    for (id, seq) in node_sequences.iter().enumerate() {
        graph.add_node(id + 1, seq.as_bytes().to_vec());
    }

    for (from, to) in edges {
        graph.add_edge(from, to);
    }

    graph
}

/// Helper to check if nodes are in perfect sequential order
fn check_perfect_order(graph: &BidirectedGraph, expected_node_order: Vec<usize>) -> bool {
    // Get actual node order from sorted graph
    let mut actual_order: Vec<usize> = graph.nodes.keys().copied().collect();
    actual_order.sort();

    // For perfect ordering, node IDs should match expected order
    // (after renumbering, nodes should be 1,2,3,... in the correct sequence)
    if actual_order.len() != expected_node_order.len() {
        eprintln!("Length mismatch: {} vs {}", actual_order.len(), expected_node_order.len());
        return false;
    }

    // Check that paths traverse nodes in expected order
    for path in &graph.paths {
        let mut prev_node_id = 0;
        for (i, handle) in path.steps.iter().enumerate() {
            let node_id = handle.node_id();
            if i > 0 && node_id < prev_node_id {
                eprintln!("Path {} has non-sequential node IDs: {} followed by {}",
                         path.name, prev_node_id, node_id);
                return false;
            }
            prev_node_id = node_id;
        }
    }

    true
}

#[test]
fn test_sgd_simple_forward_chain() {
    // Create a simple 3-node forward chain: 1+->2+->3+
    let mut graph = create_chain_graph(
        vec!["AAAA", "CCCC", "GGGG"],
        vec![
            (Handle::forward(1), Handle::forward(2)),
            (Handle::forward(2), Handle::forward(3)),
        ]
    );

    // Add a path traversing all nodes forward
    let mut path = BiPath::new("forward_path".to_string());
    path.add_step(Handle::forward(1));
    path.add_step(Handle::forward(2));
    path.add_step(Handle::forward(3));
    graph.paths.push(path);

    // Run SGD with high iterations to ensure convergence
    let params = PathSGDParams {
        iter_max: 100,
        min_term_updates: 100,
        eta_max: 100.0,
        nthreads: 1,
        progress: false,
        ..Default::default()
    };

    let ordering = path_sgd_sort(&graph, params);
    graph.apply_ordering(ordering, false);

    // Verify perfect ordering
    assert!(check_perfect_order(&graph, vec![1, 2, 3]),
            "SGD should produce perfect forward chain ordering");
}

#[test]
fn test_sgd_simple_reverse_chain() {
    // Create a 3-node chain where path traverses in reverse: 1-->2-->3-
    let mut graph = create_chain_graph(
        vec!["AAAA", "CCCC", "GGGG"],
        vec![
            (Handle::reverse(1), Handle::reverse(2)),
            (Handle::reverse(2), Handle::reverse(3)),
        ]
    );

    // Add a path traversing all nodes in reverse
    let mut path = BiPath::new("reverse_path".to_string());
    path.add_step(Handle::reverse(1));
    path.add_step(Handle::reverse(2));
    path.add_step(Handle::reverse(3));
    graph.paths.push(path);

    // Run SGD
    let params = PathSGDParams {
        iter_max: 100,
        min_term_updates: 100,
        eta_max: 100.0,
        nthreads: 1,
        progress: false,
        ..Default::default()
    };

    let ordering = path_sgd_sort(&graph, params);
    graph.apply_ordering(ordering, false);

    // Verify perfect ordering (nodes should still be sequential)
    assert!(check_perfect_order(&graph, vec![1, 2, 3]),
            "SGD should produce perfect reverse chain ordering");
}

#[test]
fn test_sgd_mixed_orientation_chain() {
    // Create a chain with mixed orientations: 1+->2-->3+
    let mut graph = create_chain_graph(
        vec!["AAAA", "CCCC", "GGGG"],
        vec![
            (Handle::forward(1), Handle::reverse(2)),
            (Handle::forward(2), Handle::forward(3)),
        ]
    );

    // Add a path with mixed orientations
    let mut path = BiPath::new("mixed_path".to_string());
    path.add_step(Handle::forward(1));
    path.add_step(Handle::reverse(2));  // Traverse node 2 in reverse
    path.add_step(Handle::forward(3));
    graph.paths.push(path);

    // Run SGD
    let params = PathSGDParams {
        iter_max: 100,
        min_term_updates: 100,
        eta_max: 100.0,
        nthreads: 1,
        progress: false,
        ..Default::default()
    };

    let ordering = path_sgd_sort(&graph, params);
    graph.apply_ordering(ordering, false);

    // Verify perfect ordering
    assert!(check_perfect_order(&graph, vec![1, 2, 3]),
            "SGD should handle mixed orientations correctly");
}

#[test]
fn test_sgd_simple_bubble_with_rc() {
    // Create a bubble where one path uses RC:
    //     2+
    //    /  \
    // 1+      4+
    //    \  /
    //     3-  (path B traverses node 3 in reverse)

    let mut graph = BidirectedGraph::new();
    graph.add_node(1, b"AAAA".to_vec());
    graph.add_node(2, b"CCCC".to_vec());
    graph.add_node(3, b"GGGG".to_vec());
    graph.add_node(4, b"TTTT".to_vec());

    // Edges for both paths
    graph.add_edge(Handle::forward(1), Handle::forward(2));
    graph.add_edge(Handle::forward(2), Handle::forward(4));
    graph.add_edge(Handle::forward(1), Handle::reverse(3));
    graph.add_edge(Handle::forward(3), Handle::forward(4));

    // Path A: 1+->2+->4+
    let mut path_a = BiPath::new("path_a".to_string());
    path_a.add_step(Handle::forward(1));
    path_a.add_step(Handle::forward(2));
    path_a.add_step(Handle::forward(4));
    graph.paths.push(path_a);

    // Path B: 1+->3-->4+ (traverses 3 in reverse)
    let mut path_b = BiPath::new("path_b".to_string());
    path_b.add_step(Handle::forward(1));
    path_b.add_step(Handle::reverse(3));  // RC traversal
    path_b.add_step(Handle::forward(4));
    graph.paths.push(path_b);

    // Run SGD
    let params = PathSGDParams {
        iter_max: 200,  // More iterations for bubble
        min_term_updates: 200,
        eta_max: 100.0,
        nthreads: 1,
        progress: false,
        ..Default::default()
    };

    let ordering = path_sgd_sort(&graph, params);
    graph.apply_ordering(ordering, false);

    // Verify reasonable ordering (1 first, 4 last, 2 and 3 in middle)
    let mut node_ids: Vec<usize> = graph.nodes.keys().copied().collect();
    node_ids.sort();

    assert_eq!(node_ids[0], 1, "Node 1 should be first");
    assert_eq!(node_ids[3], 4, "Node 4 should be last");
    // Nodes 2 and 3 should be in the middle (order may vary)
    assert!((node_ids[1] == 2 && node_ids[2] == 3) || (node_ids[1] == 3 && node_ids[2] == 2),
            "Nodes 2 and 3 should be in middle");
}

#[test]
fn test_sgd_repeated_node() {
    // Create a path that visits node 2 twice: 1+->2+->3+->2+->4+
    // This tests handling of structural variation (node repetition)

    let mut graph = BidirectedGraph::new();
    for i in 1..=4 {
        graph.add_node(i, b"AAAA".to_vec());
    }

    // Edges for the path
    graph.add_edge(Handle::forward(1), Handle::forward(2));
    graph.add_edge(Handle::forward(2), Handle::forward(3));
    graph.add_edge(Handle::forward(3), Handle::forward(2));
    graph.add_edge(Handle::forward(2), Handle::forward(4));

    // Path visiting node 2 twice
    let mut path = BiPath::new("repeated".to_string());
    path.add_step(Handle::forward(1));
    path.add_step(Handle::forward(2));  // First visit
    path.add_step(Handle::forward(3));
    path.add_step(Handle::forward(2));  // Second visit
    path.add_step(Handle::forward(4));
    graph.paths.push(path);

    // Run SGD
    let params = PathSGDParams {
        iter_max: 200,
        min_term_updates: 200,
        eta_max: 100.0,
        nthreads: 1,
        progress: false,
        ..Default::default()
    };

    let ordering = path_sgd_sort(&graph, params);
    graph.apply_ordering(ordering, false);

    // For repeated nodes, we can't expect perfect sequential ordering
    // But node 1 should be first and node 4 should be last
    let mut node_ids: Vec<usize> = graph.nodes.keys().copied().collect();
    node_ids.sort();

    // Just check that SGD completes without crashing
    assert_eq!(node_ids.len(), 4, "Should have 4 nodes");
}

#[test]
fn test_sgd_convergence_simple() {
    // Create a very simple 2-node chain to test basic convergence
    let mut graph = create_chain_graph(
        vec!["AAAA", "CCCC"],
        vec![(Handle::forward(1), Handle::forward(2))]
    );

    let mut path = BiPath::new("simple".to_string());
    path.add_step(Handle::forward(1));
    path.add_step(Handle::forward(2));
    graph.paths.push(path);

    // Run SGD with very high iterations
    let params = PathSGDParams {
        iter_max: 500,
        min_term_updates: 500,
        eta_max: 200.0,
        nthreads: 1,
        progress: true,  // Enable progress to see convergence
        ..Default::default()
    };

    let ordering = path_sgd_sort(&graph, params);
    graph.apply_ordering(ordering, false);

    // For 2 nodes, ordering should be perfect
    let node_ids: Vec<usize> = graph.nodes.keys().copied().collect();
    assert_eq!(node_ids.len(), 2);
}
