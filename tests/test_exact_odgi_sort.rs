use seqrush::bidirected_graph::Handle;
use seqrush::bidirected_ops::BidirectedGraph;

#[test]
fn test_exact_odgi_heads_and_tails() {
    let mut graph = BidirectedGraph::new();
    
    // Create linear graph: 1 -> 2 -> 3
    graph.add_node(1, b"A".to_vec());
    graph.add_node(2, b"C".to_vec());
    graph.add_node(3, b"G".to_vec());
    
    graph.add_edge(Handle::forward(1), Handle::forward(2));
    graph.add_edge(Handle::forward(2), Handle::forward(3));
    
    let heads = graph.find_head_nodes();
    assert_eq!(heads.len(), 1);
    assert_eq!(heads[0].node_id(), 1);
    
    let tails = graph.find_tail_nodes();
    assert_eq!(tails.len(), 1);
    assert_eq!(tails[0].node_id(), 3);
}

#[test]
fn test_exact_odgi_topological_sort_dag() {
    let mut dag = BidirectedGraph::new();
    
    // Create DAG:
    //   1 -> 2
    //   1 -> 3
    //   2 -> 4
    //   3 -> 4
    dag.add_node(1, b"A".to_vec());
    dag.add_node(2, b"C".to_vec());
    dag.add_node(3, b"G".to_vec());
    dag.add_node(4, b"T".to_vec());
    
    dag.add_edge(Handle::forward(1), Handle::forward(2));
    dag.add_edge(Handle::forward(1), Handle::forward(3));
    dag.add_edge(Handle::forward(2), Handle::forward(4));
    dag.add_edge(Handle::forward(3), Handle::forward(4));
    
    let ordering = dag.exact_odgi_topological_order(true, false, false);
    
    // Should have all 4 nodes
    assert_eq!(ordering.len(), 4);
    
    // Node 1 should come first (it's the only head)
    assert_eq!(ordering[0].node_id(), 1);
    
    // Node 4 should come last (it's the only tail)
    assert_eq!(ordering[3].node_id(), 4);
    
    // Nodes 2 and 3 should be in between
    let middle_nodes: Vec<usize> = ordering[1..3].iter().map(|h| h.node_id()).collect();
    assert!(middle_nodes.contains(&2));
    assert!(middle_nodes.contains(&3));
}

#[test]
fn test_exact_odgi_topological_sort_with_cycle() {
    let mut graph = BidirectedGraph::new();
    
    // Create graph with cycle:
    // 3 -> 1 -> 2 -> 4
    //      ^---------/
    graph.add_node(3, b"A".to_vec());
    graph.add_node(1, b"C".to_vec());
    graph.add_node(2, b"G".to_vec());
    graph.add_node(4, b"T".to_vec());
    
    graph.add_edge(Handle::forward(3), Handle::forward(1));
    graph.add_edge(Handle::forward(1), Handle::forward(2));
    graph.add_edge(Handle::forward(2), Handle::forward(4));
    graph.add_edge(Handle::forward(4), Handle::forward(1)); // Creates cycle
    
    let ordering = graph.exact_odgi_topological_order(true, false, false);
    
    // Should still order all nodes
    assert_eq!(ordering.len(), 4);
    
    // Node 3 is the only head, should come first
    assert_eq!(ordering[0].node_id(), 3);
}

#[test]
fn test_exact_odgi_with_bidirected_edges() {
    let mut graph = BidirectedGraph::new();
    
    // Create graph with bidirected edges
    graph.add_node(1, b"ATG".to_vec());
    graph.add_node(2, b"CGA".to_vec());
    graph.add_node(3, b"TAC".to_vec());
    
    // Mix of forward and reverse edges
    graph.add_edge(Handle::forward(1), Handle::reverse(2));
    graph.add_edge(Handle::forward(2), Handle::forward(3));
    
    let ordering = graph.exact_odgi_topological_order(true, false, false);
    
    // Should order all nodes
    assert_eq!(ordering.len(), 3);
    
    // All handles in ordering should be forward orientation
    for handle in &ordering {
        assert!(!handle.is_reverse());
    }
}

#[test]
fn test_exact_odgi_apply_ordering() {
    let mut graph = BidirectedGraph::new();
    
    // Create graph with non-sequential IDs
    graph.add_node(10, b"A".to_vec());
    graph.add_node(5, b"C".to_vec());
    graph.add_node(15, b"G".to_vec());
    
    graph.add_edge(Handle::forward(10), Handle::forward(5));
    graph.add_edge(Handle::forward(5), Handle::forward(15));
    
    // Apply topological ordering
    graph.apply_exact_odgi_ordering(false);
    
    // Nodes should be renumbered 1, 2, 3
    let node_ids: Vec<usize> = graph.nodes.keys().cloned().collect();
    assert_eq!(node_ids.len(), 3);
    assert!(node_ids.contains(&1));
    assert!(node_ids.contains(&2));
    assert!(node_ids.contains(&3));
}