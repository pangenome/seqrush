/// Test program for Ygs sorting
/// Usage: cargo run --bin test_ygs_sort -- input.gfa output.gfa

use seqrush::bidirected_ops::BidirectedGraph;
use seqrush::ygs_sort::{YgsParams, ygs_sort};
use std::env;
use std::fs;
use std::io::{Read, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: {} <input.gfa> <output.gfa> [--verbose] [--threads N]", args[0]);
        eprintln!("\nApplies the exact 'odgi sort -p Ygs' algorithm:");
        eprintln!("  Y = Path-guided SGD");
        eprintln!("  g = Grooming");
        eprintln!("  s = Topological sort (heads only)");
        std::process::exit(1);
    }

    let input_file = &args[1];
    let output_file = &args[2];
    let verbose = args.contains(&"--verbose".to_string()) || args.contains(&"-v".to_string());

    // Parse threads argument
    let mut nthreads = 1;
    for i in 0..args.len() {
        if args[i] == "--threads" || args[i] == "-t" {
            if i + 1 < args.len() {
                nthreads = args[i + 1].parse().unwrap_or(1);
            }
        }
    }

    // Read input GFA
    if verbose {
        eprintln!("[test_ygs_sort] Reading input GFA from {}", input_file);
    }

    let gfa_content = fs::read_to_string(input_file)?;

    // Parse GFA into BidirectedGraph
    // For now, we'll use a simple parser
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
    use seqrush::bidirected_graph::Handle;
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

    if verbose {
        eprintln!("[test_ygs_sort] Loaded graph: {} nodes, {} edges, {} paths",
                 graph.nodes.len(), graph.edges.len(), graph.paths.len());
    }

    // Calculate parameters and apply Ygs sort
    let params = YgsParams::from_graph(&graph, verbose, nthreads);

    ygs_sort(&mut graph, &params);

    // Write output GFA
    if verbose {
        eprintln!("[test_ygs_sort] Writing sorted GFA to {}", output_file);
    }

    let mut output = fs::File::create(output_file)?;
    graph.write_gfa(&mut output)?;

    if verbose {
        eprintln!("[test_ygs_sort] Complete! Output graph: {} nodes, {} edges",
                 graph.nodes.len(), graph.edges.len());
    }

    Ok(())
}
