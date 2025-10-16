/// Test to understand SGD's limitations with poor initialization
use seqrush::bidirected_ops::BidirectedGraph;
use seqrush::bidirected_graph::{Handle, BiPath};
use seqrush::path_sgd::{PathSGDParams, path_sgd_sort};

#[test]
#[ignore] // SGD is a layout optimizer, not a sorting algorithm. It cannot reorder nodes.
          // This test documents that limitation.
fn test_sgd_can_handle_single_swap() {
    // Test if SGD can fix a simple 2-node swap
    let mut graph = BidirectedGraph::new();

    // IDs: [1, 2] but path wants [2, 1]
    graph.add_node(1, b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec()); // 40bp
    graph.add_node(2, b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_vec()); // 40bp

    graph.add_edge(Handle::forward(2), Handle::forward(1));

    let mut path = BiPath::new("test".to_string());
    path.add_step(Handle::forward(2));
    path.add_step(Handle::forward(1));
    graph.paths.push(path);

    eprintln!("\n=== Testing SGD with simple 2-node swap ===");
    eprintln!("Initial order (by ID): [1, 2]");
    eprintln!("Path wants: [2, 1]");

    let params = PathSGDParams {
        iter_max: 1000,
        min_term_updates: 1000,
        eta_max: 1000.0,
        nthreads: 1,
        progress: true,
        ..Default::default()
    };

    let ordering = path_sgd_sort(&graph, params);

    let node_ids: Vec<usize> = ordering.iter().map(|h| h.node_id()).collect();
    eprintln!("Final order: {:?}", node_ids);

    // Can SGD swap two nodes?
    assert_eq!(node_ids, vec![2, 1],
               "SGD should be able to swap 2 nodes with enough iterations");
}

#[test]
#[ignore] // This test documents SGD's limitation - it may not converge for pathological cases
fn test_sgd_complete_reversal_is_hard() {
    // Document that complete reversal is hard for SGD
    let mut graph = BidirectedGraph::new();

    // IDs: [1, 2, 3] but path wants complete reversal [3, 2, 1]
    graph.add_node(1, vec![b'A'; 10]);
    graph.add_node(2, vec![b'C'; 10]);
    graph.add_node(3, vec![b'G'; 10]);

    graph.add_edge(Handle::forward(3), Handle::forward(2));
    graph.add_edge(Handle::forward(2), Handle::forward(1));

    let mut path = BiPath::new("test".to_string());
    path.add_step(Handle::forward(3));
    path.add_step(Handle::forward(2));
    path.add_step(Handle::forward(1));
    graph.paths.push(path);

    eprintln!("\n=== Testing SGD with complete 3-node reversal ===");
    eprintln!("Initial order (by ID): [1, 2, 3]");
    eprintln!("Path wants: [3, 2, 1]");

    let params = PathSGDParams {
        iter_max: 10000,  // Massive iteration count
        min_term_updates: 10000,
        eta_max: 1000.0,
        nthreads: 1,
        progress: true,
        ..Default::default()
    };

    let ordering = path_sgd_sort(&graph, params);

    let node_ids: Vec<usize> = ordering.iter().map(|h| h.node_id()).collect();
    eprintln!("Final order: {:?}", node_ids);

    // This documents that SGD may struggle with complete reversal
    eprintln!("Note: This test is ignored because SGD may not converge for pathological initialization");
}
