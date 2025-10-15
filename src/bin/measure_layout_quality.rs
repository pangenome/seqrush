/// Measure the quality of a 1D graph layout by calculating path entropy/loss
///
/// This measures how well the 1D node positions match the actual genomic distances
/// along paths. Lower values indicate better layout quality.
use seqrush::bidirected_ops::BidirectedGraph;
use seqrush::bidirected_graph::Handle;
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: {} <graph.gfa>", args[0]);
        eprintln!();
        eprintln!("Calculates layout quality metrics:");
        eprintln!("  - Mean squared error between 1D distances and path distances");
        eprintln!("  - Mean absolute error");
        eprintln!("  - Normalized metrics per base pair");
        std::process::exit(1);
    }

    let gfa_path = &args[1];

    // Load the graph
    eprintln!("Loading graph from {}...", gfa_path);
    let gfa_content = std::fs::read_to_string(gfa_path)?;

    // Parse GFA into BidirectedGraph
    let mut graph = BidirectedGraph::new();

    // Parse S lines (segments/nodes)
    for line in gfa_content.lines() {
        if line.starts_with('S') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let node_id: usize = parts[1].parse()?;
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
                let from_id: usize = parts[1].parse()?;
                let from_orient = parts[2];
                let to_id: usize = parts[3].parse()?;
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
                    let node_id: usize = node_str.parse()?;
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

    eprintln!("Graph has {} nodes and {} paths", graph.nodes.len(), graph.paths.len());

    // Get node positions from node IDs (assuming sequential ordering represents layout)
    let mut node_positions: HashMap<usize, f64> = HashMap::new();
    let mut pos = 0.0;
    for (&node_id, node) in graph.nodes.iter() {
        node_positions.insert(node_id, pos);
        pos += node.sequence.len() as f64;
    }

    let total_length = pos;

    // Calculate metrics for each path
    let mut total_squared_error = 0.0;
    let mut total_absolute_error = 0.0;
    let mut total_steps = 0usize;
    let mut total_path_length = 0.0;

    println!("\nPer-path metrics:");
    println!("{:<40} {:>12} {:>12} {:>12}", "Path", "MSE", "MAE", "Length");
    println!("{}", "-".repeat(80));

    for path in &graph.paths {
        if path.steps.len() < 2 {
            continue;
        }

        let mut path_squared_error = 0.0;
        let mut path_absolute_error = 0.0;
        let mut path_length = 0.0;
        let mut path_steps = 0usize;

        // For each consecutive pair of nodes in the path
        for i in 0..path.steps.len() - 1 {
            let handle_a = path.steps[i];
            let handle_b = path.steps[i + 1];

            let node_a_id = handle_a.node_id();
            let node_b_id = handle_b.node_id();

            // Get genomic distance (sum of node lengths between them in the path)
            // For consecutive nodes, this is just the length of node A
            let genomic_distance = if let Some(node_a) = graph.nodes.get(&node_a_id) {
                node_a.sequence.len() as f64
            } else {
                continue;
            };

            path_length += genomic_distance;

            // Get 1D layout distance
            let pos_a = node_positions.get(&node_a_id).copied().unwrap_or(0.0);
            let pos_b = node_positions.get(&node_b_id).copied().unwrap_or(0.0);
            let layout_distance = (pos_b - pos_a).abs();

            // Calculate error
            let error = layout_distance - genomic_distance;
            path_squared_error += error * error;
            path_absolute_error += error.abs();
            path_steps += 1;
        }

        total_squared_error += path_squared_error;
        total_absolute_error += path_absolute_error;
        total_steps += path_steps;
        total_path_length += path_length;

        if path_steps > 0 {
            let mse = path_squared_error / path_steps as f64;
            let mae = path_absolute_error / path_steps as f64;

            let display_name = if path.name.len() > 40 {
                format!("{}...", &path.name[..37])
            } else {
                path.name.clone()
            };

            println!("{:<40} {:>12.2} {:>12.2} {:>12.0}",
                     display_name, mse, mae, path_length);
        }
    }

    println!("{}", "-".repeat(80));

    // Calculate overall metrics
    if total_steps > 0 {
        let overall_mse = total_squared_error / total_steps as f64;
        let overall_mae = total_absolute_error / total_steps as f64;
        let overall_rmse = overall_mse.sqrt();

        // Normalized metrics (per base pair of total length)
        let normalized_mse = total_squared_error / total_length;
        let normalized_mae = total_absolute_error / total_length;

        println!("\nOverall metrics:");
        println!("  Total steps analyzed: {}", total_steps);
        println!("  Total path length: {:.0} bp", total_path_length);
        println!("  Total graph length: {:.0} bp", total_length);
        println!();
        println!("  Mean Squared Error (MSE): {:.2}", overall_mse);
        println!("  Root Mean Squared Error (RMSE): {:.2}", overall_rmse);
        println!("  Mean Absolute Error (MAE): {:.2}", overall_mae);
        println!();
        println!("  Normalized MSE (per bp): {:.6}", normalized_mse);
        println!("  Normalized MAE (per bp): {:.6}", normalized_mae);
        println!();
        println!("Layout quality score (lower is better):");
        println!("  RMSE: {:.2} bp displacement per step", overall_rmse);
        println!("  Relative error: {:.2}%", (overall_mae / (total_path_length / total_steps as f64)) * 100.0);
    }

    Ok(())
}
