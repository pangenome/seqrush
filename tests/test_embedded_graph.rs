use seqrush::bidirected_graph::Handle;
use seqrush::embedded_graph::*;

#[test]
fn test_embedded_graph_creation() {
    let mut graph = EmbeddedGraph::new();

    // Add nodes
    let n1 = graph.add_node(b"ACGT".to_vec());
    let n2 = graph.add_node(b"TGCA".to_vec());
    let n3 = graph.add_node(b"GGCC".to_vec());

    // Create a path
    let path1 = graph.add_path("path1".to_string());

    // Build path: n1 -> n2 -> n3
    graph.extend_path(path1, Handle::forward(n1)).unwrap();
    graph.extend_path(path1, Handle::forward(n2)).unwrap();
    graph.extend_path(path1, Handle::forward(n3)).unwrap();

    // Check connections
    let n1_next = graph.get_next_steps(Handle::forward(n1));
    assert_eq!(n1_next.len(), 1);
    assert!(n1_next.contains(&Handle::forward(n2)));

    let n2_prev = graph.get_prev_steps(Handle::forward(n2));
    assert_eq!(n2_prev.len(), 1);
    assert!(n2_prev.contains(&Handle::forward(n1)));

    // Verify path sequence
    let seq = graph.get_path_sequence(path1).unwrap();
    assert_eq!(seq, b"ACGTTGCAGGCC");
}

#[test]
fn test_perfect_neighbors() {
    let mut graph = EmbeddedGraph::new();

    // Create a simple linear graph: A -> B -> C
    let a = graph.add_node(b"AAA".to_vec());
    let b = graph.add_node(b"BBB".to_vec());
    let c = graph.add_node(b"CCC".to_vec());

    let path = graph.add_path("test".to_string());
    graph.extend_path(path, Handle::forward(a)).unwrap();
    graph.extend_path(path, Handle::forward(b)).unwrap();
    graph.extend_path(path, Handle::forward(c)).unwrap();

    // A->B and B->C should be perfect neighbors
    assert!(graph.are_perfect_neighbors(Handle::forward(a), Handle::forward(b)));
    assert!(graph.are_perfect_neighbors(Handle::forward(b), Handle::forward(c)));

    // But A->C should not be
    assert!(!graph.are_perfect_neighbors(Handle::forward(a), Handle::forward(c)));
}

#[test]
fn test_branching_graph() {
    let mut graph = EmbeddedGraph::new();

    // Create a branching graph:
    //     B
    //   /   \
    // A       D
    //   \   /
    //     C

    let a = graph.add_node(b"A".to_vec());
    let b = graph.add_node(b"B".to_vec());
    let c = graph.add_node(b"C".to_vec());
    let d = graph.add_node(b"D".to_vec());

    // Path 1: A -> B -> D
    let path1 = graph.add_path("path1".to_string());
    graph.extend_path(path1, Handle::forward(a)).unwrap();
    graph.extend_path(path1, Handle::forward(b)).unwrap();
    graph.extend_path(path1, Handle::forward(d)).unwrap();

    // Path 2: A -> C -> D
    let path2 = graph.add_path("path2".to_string());
    graph.extend_path(path2, Handle::forward(a)).unwrap();
    graph.extend_path(path2, Handle::forward(c)).unwrap();
    graph.extend_path(path2, Handle::forward(d)).unwrap();

    // A has two outgoing edges, so no perfect neighbors from A
    assert!(!graph.are_perfect_neighbors(Handle::forward(a), Handle::forward(b)));
    assert!(!graph.are_perfect_neighbors(Handle::forward(a), Handle::forward(c)));

    // Similarly, D has two incoming edges
    assert!(!graph.are_perfect_neighbors(Handle::forward(b), Handle::forward(d)));
    assert!(!graph.are_perfect_neighbors(Handle::forward(c), Handle::forward(d)));
}

#[test]
fn test_compaction() {
    let mut graph = EmbeddedGraph::new();

    // Create a simple graph that can be compacted: A -> B -> C -> D
    let a = graph.add_node(b"AA".to_vec());
    let b = graph.add_node(b"BB".to_vec());
    let c = graph.add_node(b"CC".to_vec());
    let d = graph.add_node(b"DD".to_vec());

    let path = graph.add_path("test".to_string());
    graph.extend_path(path, Handle::forward(a)).unwrap();
    graph.extend_path(path, Handle::forward(b)).unwrap();
    graph.extend_path(path, Handle::forward(c)).unwrap();
    graph.extend_path(path, Handle::forward(d)).unwrap();

    // Before compaction
    assert_eq!(graph.nodes.len(), 4);
    let seq_before = graph.get_path_sequence(path).unwrap();
    assert_eq!(seq_before, b"AABBCCDD");

    // Compact
    let merged_count = graph.compact().unwrap();
    assert_eq!(merged_count, 3); // Should merge 3 pairs

    // After compaction
    assert_eq!(graph.nodes.len(), 1); // Should have one big node
    let seq_after = graph.get_path_sequence(path).unwrap();
    assert_eq!(seq_after, b"AABBCCDD"); // Sequence should be preserved
}

#[test]
fn test_reverse_complement_path() {
    let mut graph = EmbeddedGraph::new();

    // Create nodes
    let n1 = graph.add_node(b"ACGT".to_vec());
    let n2 = graph.add_node(b"TGCA".to_vec());

    // Forward path: n1+ -> n2+
    let fwd_path = graph.add_path("forward".to_string());
    graph.extend_path(fwd_path, Handle::forward(n1)).unwrap();
    graph.extend_path(fwd_path, Handle::forward(n2)).unwrap();

    // Reverse path: n2- -> n1-
    let rev_path = graph.add_path("reverse".to_string());
    graph.extend_path(rev_path, Handle::reverse(n2)).unwrap();
    graph.extend_path(rev_path, Handle::reverse(n1)).unwrap();

    // Check sequences
    let fwd_seq = graph.get_path_sequence(fwd_path).unwrap();
    let rev_seq = graph.get_path_sequence(rev_path).unwrap();

    assert_eq!(fwd_seq, b"ACGTTGCA");
    assert_eq!(rev_seq, b"TGCAACGT"); // RC of forward path
}

#[test]
fn test_compaction_with_branches() {
    let mut graph = EmbeddedGraph::new();

    // Create a graph with a bubble:
    // A -> B -> D -> E
    //  \-> C ->/

    let a = graph.add_node(b"A".to_vec());
    let b = graph.add_node(b"B".to_vec());
    let c = graph.add_node(b"C".to_vec());
    let d = graph.add_node(b"D".to_vec());
    let e = graph.add_node(b"E".to_vec());

    // Path 1: A -> B -> D -> E
    let path1 = graph.add_path("path1".to_string());
    graph.extend_path(path1, Handle::forward(a)).unwrap();
    graph.extend_path(path1, Handle::forward(b)).unwrap();
    graph.extend_path(path1, Handle::forward(d)).unwrap();
    graph.extend_path(path1, Handle::forward(e)).unwrap();

    // Path 2: A -> C -> D -> E
    let path2 = graph.add_path("path2".to_string());
    graph.extend_path(path2, Handle::forward(a)).unwrap();
    graph.extend_path(path2, Handle::forward(c)).unwrap();
    graph.extend_path(path2, Handle::forward(d)).unwrap();
    graph.extend_path(path2, Handle::forward(e)).unwrap();

    // Before compaction
    assert_eq!(graph.nodes.len(), 5);

    // Compact - only D->E should be compactable
    let merged_count = graph.compact().unwrap();
    assert_eq!(merged_count, 1); // Only D->E can be merged

    // After compaction
    assert_eq!(graph.nodes.len(), 4); // 5 - 1 = 4

    // Verify paths still work
    let seq1 = graph.get_path_sequence(path1).unwrap();
    let seq2 = graph.get_path_sequence(path2).unwrap();
    assert_eq!(seq1, b"ABDE");
    assert_eq!(seq2, b"ACDE");
}
