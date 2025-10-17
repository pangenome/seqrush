/// Standalone GFA sorting tool using Ygs pipeline
use seqrush::bidirected_ops::BidirectedGraph;
use seqrush::bidirected_graph::{Handle, BiNode, BiPath};
use seqrush::ygs_sort::{YgsParams, ygs_sort};
use clap::Parser;
use std::process;

#[derive(Parser)]
#[command(name = "sort_gfa")]
#[command(about = "Sort a GFA file using the Ygs pipeline (SGD + grooming + topological sort)")]
struct Args {
    /// Input GFA file
    #[arg(short = 'i', long)]
    input: String,

    /// Output GFA file
    #[arg(short = 'o', long)]
    output: String,

    /// Number of SGD iterations
    #[arg(long, default_value = "100")]
    iter_max: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value = "1")]
    threads: usize,

    /// Verbose output
    #[arg(short = 'v', long)]
    verbose: bool,
}

fn parse_gfa(content: &str) -> Result<BidirectedGraph, String> {
    let mut graph = BidirectedGraph::new();

    // Parse S lines (segments/nodes)
    for line in content.lines() {
        if line.starts_with('S') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let node_id: usize = parts[1].parse()
                    .map_err(|e| format!("Failed to parse node ID: {}", e))?;
                let sequence = parts[2].as_bytes().to_vec();
                let node = BiNode {
                    id: node_id,
                    sequence,
                    rank: None,
                };
                // Ensure the vec is large enough
                if node_id >= graph.nodes.len() {
                    graph.nodes.resize(node_id + 1, None);
                }
                graph.nodes[node_id] = Some(node);
            }
        }
    }

    // Parse L lines (links/edges)
    for line in content.lines() {
        if line.starts_with('L') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 5 {
                let from_id: usize = parts[1].parse()
                    .map_err(|e| format!("Failed to parse from ID: {}", e))?;
                let from_orient = parts[2];
                let to_id: usize = parts[3].parse()
                    .map_err(|e| format!("Failed to parse to ID: {}", e))?;
                let to_orient = parts[4];

                let from_handle = if from_orient == "+" {
                    Handle::forward(from_id)
                } else {
                    Handle::reverse(from_id)
                };

                let to_handle = if to_orient == "+" {
                    Handle::forward(to_id)
                } else {
                    Handle::reverse(to_id)
                };

                graph.add_edge(from_handle, to_handle);
            }
        }
    }

    // Parse P lines (paths)
    for line in content.lines() {
        if line.starts_with('P') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let path_name = parts[1].to_string();
                let mut path = BiPath::new(path_name);

                // Parse path steps (e.g., "1+,2-,3+")
                for step_str in parts[2].split(',') {
                    let step = step_str.trim();
                    if step.is_empty() {
                        continue;
                    }

                    let orient = step.chars().last().unwrap();
                    let node_id: usize = step[..step.len()-1].parse()
                        .map_err(|e| format!("Failed to parse path node ID: {}", e))?;

                    let handle = if orient == '+' {
                        Handle::forward(node_id)
                    } else {
                        Handle::reverse(node_id)
                    };

                    path.steps.push(handle);
                }

                graph.paths.push(path);
            }
        }
    }

    Ok(graph)
}

fn main() {
    let args = Args::parse();

    if args.verbose {
        eprintln!("Reading GFA from: {}", args.input);
    }

    // Read and parse the GFA file
    let content = match std::fs::read_to_string(&args.input) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error reading file: {}", e);
            process::exit(1);
        }
    };

    let mut graph = match parse_gfa(&content) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("Error parsing GFA: {}", e);
            process::exit(1);
        }
    };

    if args.verbose {
        eprintln!("Loaded graph: {} nodes, {} edges, {} paths",
                 graph.nodes.len(), graph.edges.len(), graph.paths.len());
    }

    // Calculate parameters based on graph structure
    let mut params = YgsParams::from_graph(&graph, args.verbose, args.threads);
    params.path_sgd.iter_max = args.iter_max as u64;

    if args.verbose {
        eprintln!("Starting Ygs sorting...");
    }

    // Apply Ygs sorting
    ygs_sort(&mut graph, &params);

    if args.verbose {
        eprintln!("Writing sorted graph to: {}", args.output);
    }

    // Write the sorted graph
    let mut output_buffer = Vec::new();
    if let Err(e) = graph.write_gfa(&mut output_buffer) {
        eprintln!("Error writing GFA: {}", e);
        process::exit(1);
    }

    if let Err(e) = std::fs::write(&args.output, output_buffer) {
        eprintln!("Error writing output file: {}", e);
        process::exit(1);
    }

    if args.verbose {
        eprintln!("Done!");
    }
}
