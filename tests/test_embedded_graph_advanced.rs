use seqrush::embedded_graph::*;
use seqrush::bidirected_graph::Handle;

#[test]
fn test_empty_graph() {
    let graph = EmbeddedGraph::new();
    assert_eq!(graph.nodes.len(), 0);
    assert_eq!(graph.paths.len(), 0);
    
    // Finding perfect pairs in empty graph should return empty
    let pairs = graph.find_perfect_pairs();
    assert_eq!(pairs.len(), 0);
}

#[test]
fn test_single_node_graph() {
    let mut graph = EmbeddedGraph::new();
    let n1 = graph.add_node(b"ACGT".to_vec());
    
    // Single node has no connections
    let next = graph.get_next_steps(Handle::forward(n1));
    let prev = graph.get_prev_steps(Handle::forward(n1));
    assert_eq!(next.len(), 0);
    assert_eq!(prev.len(), 0);
    
    // No perfect pairs
    let pairs = graph.find_perfect_pairs();
    assert_eq!(pairs.len(), 0);
}

#[test]
fn test_path_without_nodes() {
    let mut graph = EmbeddedGraph::new();
    let path = graph.add_path("empty_path".to_string());
    
    // Empty path should have empty sequence
    let seq = graph.get_path_sequence(path).unwrap();
    assert_eq!(seq.len(), 0);
}

#[test]
fn test_nonexistent_node_error() {
    let mut graph = EmbeddedGraph::new();
    let path = graph.add_path("test".to_string());
    
    // Try to extend path with non-existent node
    let result = graph.extend_path(path, Handle::forward(999));
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("Node 999 not found"));
}

#[test]
fn test_nonexistent_path_error() {
    let mut graph = EmbeddedGraph::new();
    graph.add_node(b"ACGT".to_vec());
    
    // Try to get sequence of non-existent path
    let result = graph.get_path_sequence(PathId(999));
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("Path PathId(999) not found"));
}

#[test]
fn test_merge_non_neighbors_error() {
    let mut graph = EmbeddedGraph::new();
    
    // Create non-adjacent nodes
    let a = graph.add_node(b"A".to_vec());
    let b = graph.add_node(b"B".to_vec());
    let c = graph.add_node(b"C".to_vec());
    
    let path = graph.add_path("test".to_string());
    graph.extend_path(path, Handle::forward(a)).unwrap();
    graph.extend_path(path, Handle::forward(c)).unwrap(); // Skip B
    
    // Try to merge non-neighbors
    let result = graph.merge_perfect_neighbors(Handle::forward(a), Handle::forward(b));
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("are not perfect neighbors"));
}

#[test]
fn test_self_loop_handling() {
    let mut graph = EmbeddedGraph::new();
    
    // Create a self-loop: A -> A -> B
    let a = graph.add_node(b"A".to_vec());
    let b = graph.add_node(b"B".to_vec());
    
    let path = graph.add_path("test".to_string());
    graph.extend_path(path, Handle::forward(a)).unwrap();
    graph.extend_path(path, Handle::forward(a)).unwrap(); // Self-loop
    graph.extend_path(path, Handle::forward(b)).unwrap();
    
    // A should have itself as both next and prev
    let a_next = graph.get_next_steps(Handle::forward(a));
    let a_prev = graph.get_prev_steps(Handle::forward(a));
    
    // A has two next steps: itself and B
    assert_eq!(a_next.len(), 2);
    assert!(a_next.contains(&Handle::forward(a)));
    assert!(a_next.contains(&Handle::forward(b)));
    
    // A has itself as prev
    assert!(a_prev.contains(&Handle::forward(a)));
    
    // A->A is not a perfect pair (A has multiple next steps)
    assert!(!graph.are_perfect_neighbors(Handle::forward(a), Handle::forward(a)));
}

#[test]
fn test_multiple_paths_same_nodes() {
    let mut graph = EmbeddedGraph::new();
    
    // Create nodes
    let a = graph.add_node(b"A".to_vec());
    let b = graph.add_node(b"B".to_vec());
    let c = graph.add_node(b"C".to_vec());
    
    // Path 1: A -> B -> C
    let path1 = graph.add_path("path1".to_string());
    graph.extend_path(path1, Handle::forward(a)).unwrap();
    graph.extend_path(path1, Handle::forward(b)).unwrap();
    graph.extend_path(path1, Handle::forward(c)).unwrap();
    
    // Path 2: A -> B -> C (same nodes)
    let path2 = graph.add_path("path2".to_string());
    graph.extend_path(path2, Handle::forward(a)).unwrap();
    graph.extend_path(path2, Handle::forward(b)).unwrap();
    graph.extend_path(path2, Handle::forward(c)).unwrap();
    
    // Nodes should track both paths
    let node_a = &graph.nodes[&a];
    assert_eq!(node_a.path_steps.len(), 2);
    assert!(node_a.path_steps.contains_key(&(path1, false)));
    assert!(node_a.path_steps.contains_key(&(path2, false)));
    
    // Still perfect neighbors (all paths agree)
    assert!(graph.are_perfect_neighbors(Handle::forward(a), Handle::forward(b)));
    assert!(graph.are_perfect_neighbors(Handle::forward(b), Handle::forward(c)));
}

#[test]
fn test_mixed_orientations() {
    let mut graph = EmbeddedGraph::new();
    
    // Create nodes
    let a = graph.add_node(b"ACGT".to_vec());
    let b = graph.add_node(b"TGCA".to_vec());
    
    // Path 1: A+ -> B-
    let path1 = graph.add_path("mixed1".to_string());
    graph.extend_path(path1, Handle::forward(a)).unwrap();
    graph.extend_path(path1, Handle::reverse(b)).unwrap();
    
    // Check connections
    let a_fwd_next = graph.get_next_steps(Handle::forward(a));
    assert_eq!(a_fwd_next.len(), 1);
    assert!(a_fwd_next.contains(&Handle::reverse(b)));
    
    let b_rev_prev = graph.get_prev_steps(Handle::reverse(b));
    assert_eq!(b_rev_prev.len(), 1);
    assert!(b_rev_prev.contains(&Handle::forward(a)));
    
    // They should be perfect neighbors
    assert!(graph.are_perfect_neighbors(Handle::forward(a), Handle::reverse(b)));
    
    // Verify sequence
    let seq = graph.get_path_sequence(path1).unwrap();
    assert_eq!(seq, b"ACGTTGCA"); // ACGT + reverse_complement(TGCA)
}

#[test]
fn test_complex_bubble_structure() {
    let mut graph = EmbeddedGraph::new();
    
    // Create a more complex bubble:
    // A -> B -> C -> E
    //  \        /
    //   -> D ---
    
    let a = graph.add_node(b"START".to_vec());
    let b = graph.add_node(b"TOP1".to_vec());
    let c = graph.add_node(b"TOP2".to_vec());
    let d = graph.add_node(b"BOTTOM".to_vec());
    let e = graph.add_node(b"END".to_vec());
    
    // Path 1: A -> B -> C -> E
    let path1 = graph.add_path("top_path".to_string());
    graph.extend_path(path1, Handle::forward(a)).unwrap();
    graph.extend_path(path1, Handle::forward(b)).unwrap();
    graph.extend_path(path1, Handle::forward(c)).unwrap();
    graph.extend_path(path1, Handle::forward(e)).unwrap();
    
    // Path 2: A -> D -> E
    let path2 = graph.add_path("bottom_path".to_string());
    graph.extend_path(path2, Handle::forward(a)).unwrap();
    graph.extend_path(path2, Handle::forward(d)).unwrap();
    graph.extend_path(path2, Handle::forward(e)).unwrap();
    
    // Check perfect neighbors
    assert!(graph.are_perfect_neighbors(Handle::forward(b), Handle::forward(c)));
    assert!(!graph.are_perfect_neighbors(Handle::forward(a), Handle::forward(b)));
    assert!(!graph.are_perfect_neighbors(Handle::forward(a), Handle::forward(d)));
    assert!(!graph.are_perfect_neighbors(Handle::forward(c), Handle::forward(e)));
    assert!(!graph.are_perfect_neighbors(Handle::forward(d), Handle::forward(e)));
    
    // Compact and verify
    let merged = graph.compact().unwrap();
    assert_eq!(merged, 1); // Only B->C should merge
    
    // Check final structure
    assert_eq!(graph.nodes.len(), 4); // 5 - 1 = 4
    
    // Verify paths still work
    let seq1 = graph.get_path_sequence(path1).unwrap();
    let seq2 = graph.get_path_sequence(path2).unwrap();
    assert_eq!(seq1, b"STARTTOP1TOP2END");
    assert_eq!(seq2, b"STARTBOTTOMEND");
}

#[test]
fn test_iterative_compaction() {
    let mut graph = EmbeddedGraph::new();
    
    // Create a long chain: A -> B -> C -> D -> E -> F
    let nodes: Vec<_> = (0..6)
        .map(|i| graph.add_node(vec![b'A' + i as u8]))
        .collect();
    
    let path = graph.add_path("chain".to_string());
    for &node in &nodes {
        graph.extend_path(path, Handle::forward(node)).unwrap();
    }
    
    // All should be perfect neighbors initially
    for i in 0..5 {
        assert!(graph.are_perfect_neighbors(
            Handle::forward(nodes[i]), 
            Handle::forward(nodes[i + 1])
        ));
    }
    
    // Compact should merge everything into one node
    let merged = graph.compact().unwrap();
    assert_eq!(merged, 5); // 5 merges for 6 nodes
    assert_eq!(graph.nodes.len(), 1);
    
    // Sequence should be preserved
    let seq = graph.get_path_sequence(path).unwrap();
    assert_eq!(seq, b"ABCDEF");
}

#[test]
fn test_circular_path() {
    let mut graph = EmbeddedGraph::new();
    
    // Create a circular structure: A -> B -> C -> A
    let a = graph.add_node(b"A".to_vec());
    let b = graph.add_node(b"B".to_vec());
    let c = graph.add_node(b"C".to_vec());
    
    let path = graph.add_path("circular".to_string());
    graph.extend_path(path, Handle::forward(a)).unwrap();
    graph.extend_path(path, Handle::forward(b)).unwrap();
    graph.extend_path(path, Handle::forward(c)).unwrap();
    graph.extend_path(path, Handle::forward(a)).unwrap(); // Back to A
    
    // No perfect neighbors in a cycle
    assert!(!graph.are_perfect_neighbors(Handle::forward(a), Handle::forward(b)));
    assert!(!graph.are_perfect_neighbors(Handle::forward(b), Handle::forward(c)));
    assert!(!graph.are_perfect_neighbors(Handle::forward(c), Handle::forward(a)));
    
    // No compaction possible
    let merged = graph.compact().unwrap();
    assert_eq!(merged, 0);
    assert_eq!(graph.nodes.len(), 3);
}

#[test]
fn test_handle_flip_consistency() {
    let h1 = Handle::forward(42);
    let h2 = Handle::reverse(42);
    
    assert_eq!(h1.node_id(), 42);
    assert_eq!(h2.node_id(), 42);
    assert!(!h1.is_reverse());
    assert!(h2.is_reverse());
    assert_eq!(h1.flip(), h2);
    assert_eq!(h2.flip(), h1);
    assert_eq!(h1.flip().flip(), h1);
}

#[test]
fn test_large_graph_performance() {
    let mut graph = EmbeddedGraph::new();
    
    // Create a graph with 1000 nodes in 10 paths
    let nodes_per_path = 100;
    let num_paths = 10;
    
    // Create shared start and end nodes
    let start = graph.add_node(b"START".to_vec());
    let end = graph.add_node(b"END".to_vec());
    
    for i in 0..num_paths {
        let path = graph.add_path(format!("path{}", i));
        
        // Start node
        graph.extend_path(path, Handle::forward(start)).unwrap();
        
        // Middle nodes unique to this path
        for j in 0..nodes_per_path {
            let node = graph.add_node(vec![b'A' + (i as u8 % 26), b'0' + (j as u8 % 10)]);
            graph.extend_path(path, Handle::forward(node)).unwrap();
        }
        
        // End node
        graph.extend_path(path, Handle::forward(end)).unwrap();
    }
    
    // Should have 2 shared + (100 * 10) unique = 1002 nodes
    assert_eq!(graph.nodes.len(), 1002);
    
    // Find perfect pairs
    let pairs = graph.find_perfect_pairs();
    
    // Each path has 99 perfect pairs in its middle section
    assert_eq!(pairs.len(), 99 * num_paths);
    
    // Compact
    let merged = graph.compact().unwrap();
    assert_eq!(merged, 99 * num_paths);
    
    // Should have 2 shared + 10 compacted middle nodes = 12 nodes
    assert_eq!(graph.nodes.len(), 12);
}