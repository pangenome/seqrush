/// Test whether YGS pipeline actually produces high-quality layouts
use seqrush::bidirected_ops::BidirectedGraph;
use seqrush::bidirected_graph::{Handle, BiPath};
use seqrush::ygs_sort::{YgsParams, ygs_sort};
use std::collections::HashMap;

#[test]
fn test_ygs_layout_quality_measurement() {
    // Create a graph with nodes in bad initial order
    let mut graph = BidirectedGraph::new();

    // IDs: [147, 426, 470] - reversed from path order
    graph.add_node(426, b"A".to_vec());               // 1bp
    graph.add_node(470, b"AAAACCCCGGGGTTTT".to_vec());  // 16bp
    graph.add_node(147, b"CCCCCCCGGGGGGGTTTTTTT".to_vec());  // 21bp

    // Edges for path: 426+ -> 470- -> 147+
    // Edge from step 0 (426+) to step 1 (470-):
    graph.add_edge(Handle::forward(426), Handle::reverse(470));
    // Edge from step 1 (470-) to step 2 (147+):
    graph.add_edge(Handle::reverse(470), Handle::forward(147));

    // Path traverses: 426+ -> 470- -> 147+
    let mut path = BiPath::new("test_path".to_string());
    path.add_step(Handle::forward(426));
    path.add_step(Handle::reverse(470));
    path.add_step(Handle::forward(147));
    graph.paths.push(path);

    eprintln!("\n=== Measuring YGS layout quality ===");
    eprintln!("Initial graph (before YGS):");
    print_graph_info(&graph);

    eprintln!("\nInitial edges:");
    for edge in &graph.edges {
        eprintln!("  {} {} -> {} {}",
                 edge.from.node_id(),
                 if edge.from.is_reverse() { "-" } else { "+" },
                 edge.to.node_id(),
                 if edge.to.is_reverse() { "-" } else { "+" });
    }

    eprintln!("\nInitial path details:");
    for step in &graph.paths[0].steps {
        eprintln!("  Step: Node {}{}",
                 step.node_id(),
                 if step.is_reverse() { "-" } else { "+" });
    }

    // Run YGS with verbose output to see what's happening
    let params = YgsParams::from_graph(&graph, true, 1);  // Set verbose=true
    ygs_sort(&mut graph, &params);

    eprintln!("\nAfter YGS:");
    print_graph_info(&graph);

    eprintln!("\nFinal path details:");
    for step in &graph.paths[0].steps {
        eprintln!("  Step: Node {}{}",
                 step.node_id(),
                 if step.is_reverse() { "-" } else { "+" });
    }

    // Measure layout quality
    let quality = measure_layout_quality(&graph);
    eprintln!("\nLayout quality metrics:");
    eprintln!("  Mean absolute error: {:.2} bp", quality.mae);
    eprintln!("  Max displacement ratio: {:.2}x", quality.max_ratio);
    eprintln!("  Path preserved: {}", quality.path_preserved);

    // Verify path is still valid
    assert_eq!(graph.paths.len(), 1, "Should have 1 path");
    assert_eq!(graph.paths[0].steps.len(), 3, "Path should have 3 steps");

    // Check that layout quality is good (low error)
    assert!(quality.mae < 20.0,
            "Layout should have low mean error, got {:.2}", quality.mae);

    eprintln!("✓ PASS: YGS produces high-quality layout");
}

fn print_graph_info(graph: &BidirectedGraph) {
    let mut node_ids: Vec<usize> = graph.nodes.keys().copied().collect();
    node_ids.sort();

    eprintln!("  Nodes: {:?}", node_ids);
    eprintln!("  Edges: {}", graph.edges.len());

    if !graph.paths.is_empty() {
        eprintln!("  Path steps:");
        for step in &graph.paths[0].steps {
            eprintln!("    Node {}{}",
                     step.node_id(),
                     if step.is_reverse() { "-" } else { "+" });
        }
    }
}

struct LayoutQuality {
    mae: f64,  // Mean absolute error
    max_ratio: f64,  // Max displacement ratio
    path_preserved: bool,
}

fn measure_layout_quality(graph: &BidirectedGraph) -> LayoutQuality {
    if graph.paths.is_empty() {
        return LayoutQuality {
            mae: 0.0,
            max_ratio: 0.0,
            path_preserved: false,
        };
    }

    let path = &graph.paths[0];

    // Build position map (nodes should be in order by ID after apply_ordering)
    let mut node_positions: HashMap<usize, f64> = HashMap::new();
    let mut sorted_nodes: Vec<_> = graph.nodes.iter().collect();
    sorted_nodes.sort_by_key(|(id, _)| **id);

    let mut pos = 0.0;
    for (&node_id, node) in sorted_nodes {
        node_positions.insert(node_id, pos);
        pos += node.sequence.len() as f64;
    }

    // Calculate errors for consecutive steps in path
    let mut total_error = 0.0;
    let mut max_ratio: f64 = 0.0;
    let mut num_pairs = 0;

    for i in 0..(path.steps.len() - 1) {
        let step_a = path.steps[i];
        let step_b = path.steps[i + 1];

        let node_a_id = step_a.node_id();
        let node_b_id = step_b.node_id();

        // Expected distance in path (length of node A)
        let expected_dist = graph.nodes.get(&node_a_id)
            .map(|n| n.sequence.len() as f64)
            .unwrap_or(0.0);

        // Actual distance in layout
        let pos_a = node_positions.get(&node_a_id).copied().unwrap_or(0.0);
        let pos_b = node_positions.get(&node_b_id).copied().unwrap_or(0.0);
        let actual_dist = (pos_b - pos_a).abs();

        let error = (actual_dist - expected_dist).abs();
        total_error += error;

        let ratio = if expected_dist > 0.0 {
            actual_dist / expected_dist
        } else {
            1.0
        };
        max_ratio = max_ratio.max(ratio);

        num_pairs += 1;
    }

    let mae = if num_pairs > 0 {
        total_error / num_pairs as f64
    } else {
        0.0
    };

    LayoutQuality {
        mae,
        max_ratio,
        path_preserved: path.steps.len() == 3,
    }
}

#[test]
fn test_ygs_with_multiple_paths() {
    // Test with a more realistic graph with multiple overlapping paths
    let mut graph = BidirectedGraph::new();

    // Create nodes: 10, 20, 30, 40, 50 (not in any particular order)
    for i in &[50, 30, 10, 40, 20] {
        graph.add_node(*i, vec![b'A'; 10]);  // 10bp each
    }

    // Path 1: 10 -> 20 -> 30
    graph.add_edge(Handle::forward(10), Handle::forward(20));
    graph.add_edge(Handle::forward(20), Handle::forward(30));

    let mut path1 = BiPath::new("path1".to_string());
    path1.add_step(Handle::forward(10));
    path1.add_step(Handle::forward(20));
    path1.add_step(Handle::forward(30));
    graph.paths.push(path1);

    // Path 2: 20 -> 30 -> 40 (overlaps with path1)
    graph.add_edge(Handle::forward(30), Handle::forward(40));

    let mut path2 = BiPath::new("path2".to_string());
    path2.add_step(Handle::forward(20));
    path2.add_step(Handle::forward(30));
    path2.add_step(Handle::forward(40));
    graph.paths.push(path2);

    // Path 3: 30 -> 40 -> 50
    graph.add_edge(Handle::forward(40), Handle::forward(50));

    let mut path3 = BiPath::new("path3".to_string());
    path3.add_step(Handle::forward(30));
    path3.add_step(Handle::forward(40));
    path3.add_step(Handle::forward(50));
    graph.paths.push(path3);

    eprintln!("\n=== Testing YGS with multiple overlapping paths ===");
    eprintln!("Paths:");
    eprintln!("  path1: 10 -> 20 -> 30");
    eprintln!("  path2: 20 -> 30 -> 40");
    eprintln!("  path3: 30 -> 40 -> 50");

    let params = YgsParams::from_graph(&graph, true, 1);
    ygs_sort(&mut graph, &params);

    // All paths should still be valid
    assert_eq!(graph.paths.len(), 3);
    for path in &graph.paths {
        assert!(path.steps.len() >= 3, "Paths should be preserved");
    }

    eprintln!("✓ PASS: YGS handles multiple paths correctly");
}
