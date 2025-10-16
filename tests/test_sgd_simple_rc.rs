/// Test SGD with simple RC traversal where node IDs match path order
use seqrush::bidirected_ops::BidirectedGraph;
use seqrush::bidirected_graph::{Handle, BiPath};
use seqrush::path_sgd::{PathSGDParams, path_sgd_sort};

#[test]
fn test_sgd_with_nodes_already_in_path_order() {
    // Simpler test: nodes in path order but with one RC traversal
    let mut graph = BidirectedGraph::new();

    // Use IDs that are ALREADY in path order: 1 -> 2 -> 3
    graph.add_node(1, b"AAAAAAAAAA".to_vec());       // 10bp
    graph.add_node(2, b"CCCCCCCCCCCCCCCC".to_vec()); // 16bp
    graph.add_node(3, b"GGGGGGGGGGGGGGGGGGGGGGGG".to_vec()); // 24bp

    // Edges
    graph.add_edge(Handle::forward(1), Handle::reverse(2));
    graph.add_edge(Handle::forward(2), Handle::forward(3));

    // Path: 1+ -> 2- -> 3+
    let mut path = BiPath::new("test_path".to_string());
    path.add_step(Handle::forward(1));
    path.add_step(Handle::reverse(2));  // RC traversal
    path.add_step(Handle::forward(3));
    graph.paths.push(path);

    eprintln!("\n=== Testing SGD with nodes already in path order ===");
    eprintln!("Node IDs: 1, 2, 3 (sorted)");
    eprintln!("Path: 1+->2-->3+");
    eprintln!("Node IDs already match path order, so SGD should keep them in order\n");

    let params = PathSGDParams {
        iter_max: 100,
        min_term_updates: 100,
        eta_max: 100.0,
        nthreads: 1,
        progress: true,
        ..Default::default()
    };

    let ordering = path_sgd_sort(&graph, params);

    eprintln!("\nSGD ordering:");
    for (i, handle) in ordering.iter().enumerate() {
        eprintln!("  Position {}: Node {}{}",
                 i,
                 handle.node_id(),
                 if handle.is_reverse() { "-" } else { "+" });
    }

    let node_ids: Vec<usize> = ordering.iter().map(|h| h.node_id()).collect();
    eprintln!("\nFinal node order: {:?}", node_ids);

    assert_eq!(node_ids, vec![1, 2, 3],
               "Nodes should remain in path order [1, 2, 3] since IDs already match path");

    eprintln!("✓ PASS: Nodes correctly remain in order despite RC traversal");
}
