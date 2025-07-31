use seqrush::embedded_graph::EmbeddedGraph;
use seqrush::bidirected_graph::Handle;

fn main() {
    println!("Testing embedded graph with circular paths...\n");
    
    // Test 1: Simple circular path
    println!("Test 1: Simple circular path (A->B->C->A)");
    let mut graph = EmbeddedGraph::new();
    
    // Add nodes
    let node_a = graph.add_node(b"A".to_vec());
    let node_b = graph.add_node(b"B".to_vec());
    let node_c = graph.add_node(b"C".to_vec());
    
    // Create circular path
    let path_id = graph.add_path("circular".to_string());
    graph.extend_path(path_id, Handle::forward(node_a)).unwrap();
    graph.extend_path(path_id, Handle::forward(node_b)).unwrap();
    graph.extend_path(path_id, Handle::forward(node_c)).unwrap();
    graph.extend_path(path_id, Handle::forward(node_a)).unwrap(); // Back to A
    
    // Verify connections
    let next_from_a = graph.get_next_steps(Handle::forward(node_a));
    let next_from_b = graph.get_next_steps(Handle::forward(node_b));
    let next_from_c = graph.get_next_steps(Handle::forward(node_c));
    
    println!("  Next from A: {:?}", next_from_a);
    println!("  Next from B: {:?}", next_from_b);
    println!("  Next from C: {:?}", next_from_c);
    
    // Test 2: Path visiting same node multiple times
    println!("\nTest 2: Path visiting same node multiple times (A->B->A->C->B->C)");
    let mut graph2 = EmbeddedGraph::new();
    
    let node_a2 = graph2.add_node(b"A".to_vec());
    let node_b2 = graph2.add_node(b"B".to_vec());
    let node_c2 = graph2.add_node(b"C".to_vec());
    
    let path_id2 = graph2.add_path("multi_visit".to_string());
    graph2.extend_path(path_id2, Handle::forward(node_a2)).unwrap();
    graph2.extend_path(path_id2, Handle::forward(node_b2)).unwrap();
    graph2.extend_path(path_id2, Handle::forward(node_a2)).unwrap(); // Visit A again
    graph2.extend_path(path_id2, Handle::forward(node_c2)).unwrap();
    graph2.extend_path(path_id2, Handle::forward(node_b2)).unwrap(); // Visit B again
    graph2.extend_path(path_id2, Handle::forward(node_c2)).unwrap(); // Visit C again
    
    // Check node A has two occurrences
    if let Some(node_a_obj) = graph2.nodes.get(&node_a2) {
        if let Some(steps) = node_a_obj.path_steps.get(&path_id2) {
            println!("  Node A has {} occurrences in path", steps.len());
            for (i, step) in steps.iter().enumerate() {
                println!("    Occurrence {}: prev={:?}, next={:?}", 
                         i, step.prev, step.next);
            }
        }
    }
    
    // Check node B has two occurrences
    if let Some(node_b_obj) = graph2.nodes.get(&node_b2) {
        if let Some(steps) = node_b_obj.path_steps.get(&path_id2) {
            println!("  Node B has {} occurrences in path", steps.len());
            for (i, step) in steps.iter().enumerate() {
                println!("    Occurrence {}: prev={:?}, next={:?}", 
                         i, step.prev, step.next);
            }
        }
    }
    
    // Test 3: Figure-8 path
    println!("\nTest 3: Figure-8 path (A->B->C->D->B->E->F->D->G)");
    let mut graph3 = EmbeddedGraph::new();
    
    let nodes: Vec<_> = (0..7).map(|i| {
        graph3.add_node(vec![b'A' + i as u8])
    }).collect();
    
    let path_id3 = graph3.add_path("figure8".to_string());
    let path_nodes = vec![0, 1, 2, 3, 1, 4, 5, 3, 6]; // A->B->C->D->B->E->F->D->G
    
    for &node_idx in &path_nodes {
        graph3.extend_path(path_id3, Handle::forward(nodes[node_idx])).unwrap();
    }
    
    // Check nodes B(1) and D(3) which are visited twice
    for &node_idx in &[1, 3] {
        if let Some(node) = graph3.nodes.get(&nodes[node_idx]) {
            if let Some(steps) = node.path_steps.get(&path_id3) {
                let node_name = (b'A' + node_idx as u8) as char;
                println!("  Node {} has {} occurrences", node_name, steps.len());
            }
        }
    }
    
    // Test 4: Self-loop
    println!("\nTest 4: Self-loop (A->A->B)");
    let mut graph4 = EmbeddedGraph::new();
    
    let node_a4 = graph4.add_node(b"A".to_vec());
    let node_b4 = graph4.add_node(b"B".to_vec());
    
    let path_id4 = graph4.add_path("self_loop".to_string());
    graph4.extend_path(path_id4, Handle::forward(node_a4)).unwrap();
    graph4.extend_path(path_id4, Handle::forward(node_a4)).unwrap(); // Self-loop
    graph4.extend_path(path_id4, Handle::forward(node_b4)).unwrap();
    
    if let Some(node_a_obj) = graph4.nodes.get(&node_a4) {
        if let Some(steps) = node_a_obj.path_steps.get(&path_id4) {
            println!("  Node A has {} occurrences", steps.len());
            for (i, step) in steps.iter().enumerate() {
                println!("    Occurrence {}: prev={:?}, next={:?}", 
                         i, step.prev, step.next);
            }
        }
    }
    
    // Test path reconstruction
    println!("\nTest 5: Path reconstruction");
    match graph2.get_path_sequence(path_id2) {
        Ok(seq) => {
            println!("  Reconstructed path: {}", String::from_utf8_lossy(&seq));
            println!("  Expected: ABACBC");
            if seq == b"ABACBC" {
                println!("  ✓ Path reconstruction correct!");
            } else {
                println!("  ✗ Path reconstruction incorrect!");
            }
        }
        Err(e) => println!("  Error: {}", e),
    }
    
    println!("\nAll tests completed!");
}