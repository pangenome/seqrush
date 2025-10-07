use seqrush::bidirected_graph::Handle;
use seqrush::bidirected_ops::BidirectedGraph;
use seqrush::embedded_builder::convert_to_embedded;

fn main() {
    println!("Testing bidirected to embedded graph conversion...\n");

    // Test 1: Simple linear path
    println!("Test 1: Simple linear path");
    let mut bi_graph = BidirectedGraph::new();
    bi_graph.add_node(1, b"A".to_vec());
    bi_graph.add_node(2, b"B".to_vec());
    bi_graph.add_node(3, b"C".to_vec());
    bi_graph.build_path("seq1".to_string(), vec![(1, false), (2, false), (3, false)]);

    match convert_to_embedded(&bi_graph) {
        Ok(embedded) => {
            println!("  ✓ Conversion successful");
            println!(
                "  Nodes: {}, Paths: {}",
                embedded.nodes.len(),
                embedded.paths.len()
            );

            // Check path reconstruction
            let path_id = seqrush::embedded_graph::PathId(0);
            match embedded.get_path_sequence(path_id) {
                Ok(seq) => println!("  Path sequence: {}", String::from_utf8_lossy(&seq)),
                Err(e) => println!("  Error getting path: {}", e),
            }
        }
        Err(e) => println!("  ✗ Conversion failed: {}", e),
    }

    // Test 2: Path with reverse complement
    println!("\nTest 2: Path with reverse complement");
    let mut bi_graph2 = BidirectedGraph::new();
    bi_graph2.add_node(1, b"ACGT".to_vec());
    bi_graph2.add_node(2, b"TTTT".to_vec());
    bi_graph2.build_path("seq_rc".to_string(), vec![(1, false), (2, true)]);

    match convert_to_embedded(&bi_graph2) {
        Ok(embedded) => {
            println!("  ✓ Conversion successful");
            let path_id = seqrush::embedded_graph::PathId(0);
            match embedded.get_path_sequence(path_id) {
                Ok(seq) => {
                    println!("  Path sequence: {}", String::from_utf8_lossy(&seq));
                    println!("  Expected: ACGTAAAA");
                }
                Err(e) => println!("  Error getting path: {}", e),
            }
        }
        Err(e) => println!("  ✗ Conversion failed: {}", e),
    }

    // Test 3: Multiple paths with shared nodes
    println!("\nTest 3: Multiple paths with shared nodes");
    let mut bi_graph3 = BidirectedGraph::new();
    bi_graph3.add_node(1, b"A".to_vec());
    bi_graph3.add_node(2, b"B".to_vec());
    bi_graph3.add_node(3, b"C".to_vec());
    bi_graph3.add_node(4, b"D".to_vec());

    // Path 1: A->B->D
    bi_graph3.build_path(
        "path1".to_string(),
        vec![(1, false), (2, false), (4, false)],
    );
    // Path 2: A->C->D
    bi_graph3.build_path(
        "path2".to_string(),
        vec![(1, false), (3, false), (4, false)],
    );

    match convert_to_embedded(&bi_graph3) {
        Ok(embedded) => {
            println!("  ✓ Conversion successful");
            println!(
                "  Nodes: {}, Paths: {}",
                embedded.nodes.len(),
                embedded.paths.len()
            );

            for i in 0..2 {
                let path_id = seqrush::embedded_graph::PathId(i);
                match embedded.get_path_sequence(path_id) {
                    Ok(seq) => println!("  Path {}: {}", i, String::from_utf8_lossy(&seq)),
                    Err(e) => println!("  Error getting path {}: {}", i, e),
                }
            }
        }
        Err(e) => println!("  ✗ Conversion failed: {}", e),
    }

    // Test 4: Circular path with multiple visits
    println!("\nTest 4: Circular path with multiple visits");
    let mut bi_graph4 = BidirectedGraph::new();
    bi_graph4.add_node(1, b"A".to_vec());
    bi_graph4.add_node(2, b"B".to_vec());
    bi_graph4.add_node(3, b"C".to_vec());

    // Path: A->B->C->B->A (visits B twice, circular)
    bi_graph4.build_path(
        "circular".to_string(),
        vec![(1, false), (2, false), (3, false), (2, false), (1, false)],
    );

    match convert_to_embedded(&bi_graph4) {
        Ok(embedded) => {
            println!("  ✓ Conversion successful");

            // Check multiple visits
            for (node_id, node) in &embedded.nodes {
                if node.sequence == b"B" {
                    let path_id = seqrush::embedded_graph::PathId(0);
                    if let Some(steps) = node.path_steps.get(&path_id) {
                        println!("  Node B (id={}) has {} visits", node_id, steps.len());
                    }
                }
            }

            let path_id = seqrush::embedded_graph::PathId(0);
            match embedded.get_path_sequence(path_id) {
                Ok(seq) => {
                    println!("  Path sequence: {}", String::from_utf8_lossy(&seq));
                    println!("  Expected: ABCBA");
                }
                Err(e) => println!("  Error getting path: {}", e),
            }
        }
        Err(e) => println!("  ✗ Conversion failed: {}", e),
    }

    println!("\nAll conversion tests completed!");
}
