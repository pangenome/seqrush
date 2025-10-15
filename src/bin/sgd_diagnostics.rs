/// Diagnostic tool to analyze SGD behavior and identify problematic node pairs
/// Particularly focusing on reverse complement paths
use seqrush::bidirected_ops::BidirectedGraph;
use seqrush::bidirected_graph::Handle;
use std::collections::HashMap;

#[derive(Debug)]
struct NodePairDiagnostic {
    node_a: usize,
    node_b: usize,
    handle_a: Handle,
    handle_b: Handle,
    _path_name: String,
    pos_in_path_a: usize,
    pos_in_path_b: usize,
    path_distance: f64,
    sgd_distance: f64,
    ratio: f64,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <input.gfa>", args[0]);
        eprintln!("Analyzes SGD behavior to find problematic node pairs");
        std::process::exit(1);
    }

    let input_path = &args[1];

    eprintln!("[sgd_diagnostics] Loading graph from {}", input_path);
    let gfa_content = std::fs::read_to_string(input_path)
        .expect("Failed to read GFA file");

    let mut graph = BidirectedGraph::new();

    // Parse S lines (segments/nodes)
    for line in gfa_content.lines() {
        if line.starts_with('S') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let node_id: usize = parts[1].parse().expect("Invalid node ID");
                let sequence = parts[2].as_bytes().to_vec();
                graph.add_node(node_id, sequence);
            }
        }
    }

    // Parse L lines (links/edges)
    for line in gfa_content.lines() {
        if line.starts_with('L') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 5 {
                let from_id: usize = parts[1].parse().expect("Invalid from node ID");
                let from_orient = parts[2];
                let to_id: usize = parts[3].parse().expect("Invalid to node ID");
                let to_orient = parts[4];

                let from = if from_orient == "+" {
                    Handle::forward(from_id)
                } else {
                    Handle::reverse(from_id)
                };

                let to = if to_orient == "+" {
                    Handle::forward(to_id)
                } else {
                    Handle::reverse(to_id)
                };

                graph.add_edge(from, to);
            }
        }
    }

    // Parse P lines (paths)
    use seqrush::bidirected_graph::BiPath;
    for line in gfa_content.lines() {
        if line.starts_with('P') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let path_name = parts[1].to_string();
                let mut path = BiPath::new(path_name);

                let steps_str = parts[2];
                for step_str in steps_str.split(',') {
                    let node_str = step_str.trim_end_matches('+').trim_end_matches('-');
                    let node_id: usize = node_str.parse().expect("Invalid node ID in path");
                    let is_reverse = step_str.ends_with('-');

                    let handle = if is_reverse {
                        Handle::reverse(node_id)
                    } else {
                        Handle::forward(node_id)
                    };

                    path.add_step(handle);
                }

                graph.paths.push(path);
            }
        }
    }

    eprintln!("[sgd_diagnostics] Graph loaded:");
    eprintln!("  Nodes: {}", graph.nodes.len());
    eprintln!("  Paths: {}", graph.paths.len());
    eprintln!("  Edges: {}", graph.edges.len());

    // Analyze path orientations
    eprintln!("\n[sgd_diagnostics] Path orientation analysis:");
    for path in &graph.paths {
        let total_steps = path.steps.len();
        let reverse_steps = path.steps.iter()
            .filter(|h| h.is_reverse())
            .count();
        let forward_steps = total_steps - reverse_steps;
        let pct_reverse = (reverse_steps as f64 / total_steps as f64) * 100.0;

        eprintln!("  {}: {} steps, {} forward, {} reverse ({:.1}% reverse)",
                 path.name, total_steps, forward_steps, reverse_steps, pct_reverse);
    }

    // Build sorted node list for SGD positioning
    let mut sorted_nodes: Vec<_> = graph.nodes.iter().collect();
    sorted_nodes.sort_by_key(|(node_id, _)| *node_id);

    let mut node_id_to_sgd_pos: HashMap<usize, f64> = HashMap::new();
    let mut pos = 0.0;
    for (&node_id, node) in &sorted_nodes {
        node_id_to_sgd_pos.insert(node_id, pos);
        pos += node.sequence.len() as f64;
    }

    // Analyze each path for problematic pairs
    eprintln!("\n[sgd_diagnostics] Analyzing path step pairs:");

    for path in &graph.paths {
        eprintln!("\n  Path: {}", path.name);

        // Calculate positions for each step
        let mut step_positions: Vec<usize> = Vec::new();
        let mut cumulative_pos = 0;
        for handle in &path.steps {
            step_positions.push(cumulative_pos);
            if let Some(node) = graph.nodes.get(&handle.node_id()) {
                cumulative_pos += node.sequence.len();
            }
        }

        // Look at adjacent step pairs
        let mut problematic_pairs = Vec::new();
        for i in 0..path.steps.len().saturating_sub(1) {
            let handle_a = path.steps[i];
            let handle_b = path.steps[i + 1];

            let pos_in_path_a = step_positions[i];
            let pos_in_path_b = step_positions[i + 1];
            let path_distance = (pos_in_path_b as f64 - pos_in_path_a as f64).abs();

            let sgd_pos_a = node_id_to_sgd_pos.get(&handle_a.node_id()).copied().unwrap_or(0.0);
            let sgd_pos_b = node_id_to_sgd_pos.get(&handle_b.node_id()).copied().unwrap_or(0.0);
            let sgd_distance = (sgd_pos_b - sgd_pos_a).abs();

            // Flag if SGD distance is much larger than path distance
            let ratio = if path_distance > 0.0 {
                sgd_distance / path_distance
            } else {
                0.0
            };

            if ratio > 10.0 {  // SGD distance is 10x larger than expected
                problematic_pairs.push(NodePairDiagnostic {
                    node_a: handle_a.node_id(),
                    node_b: handle_b.node_id(),
                    handle_a,
                    handle_b,
                    _path_name: path.name.clone(),
                    pos_in_path_a,
                    pos_in_path_b,
                    path_distance,
                    sgd_distance,
                    ratio,
                });
            }
        }

        // Report problematic pairs
        if problematic_pairs.is_empty() {
            eprintln!("    No problematic adjacent pairs found");
        } else {
            eprintln!("    Found {} problematic adjacent pairs:", problematic_pairs.len());
            for diag in &problematic_pairs {
                eprintln!("      Node {}{}->{}{}:",
                         diag.node_a,
                         if diag.handle_a.is_reverse() { "-" } else { "+" },
                         diag.node_b,
                         if diag.handle_b.is_reverse() { "-" } else { "+" });
                eprintln!("        Path positions: {} -> {} (dist={:.0}bp)",
                         diag.pos_in_path_a, diag.pos_in_path_b, diag.path_distance);
                eprintln!("        SGD positions: {:.0} -> {:.0} (dist={:.0})",
                         node_id_to_sgd_pos.get(&diag.node_a).unwrap_or(&0.0),
                         node_id_to_sgd_pos.get(&diag.node_b).unwrap_or(&0.0),
                         diag.sgd_distance);
                eprintln!("        Ratio: {:.1}x (SGD dist / path dist)", diag.ratio);
            }
        }
    }

    eprintln!("\n[sgd_diagnostics] Analysis complete");
}
