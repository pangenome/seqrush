/// Test program to output graph at each Ygs stage
use seqrush::bidirected_ops::BidirectedGraph;
use seqrush::path_sgd_exact::{path_sgd_sort, PathSGDParams, XPIndex};
use std::env;
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <input.gfa>", args[0]);
        std::process::exit(1);
    }

    let input_file = &args[1];
    let verbose = true;

    eprintln!("[test_ygs_stages] Reading input GFA from {}", input_file);
    let gfa_content = fs::read_to_string(input_file)?;

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

    eprintln!("[test_ygs_stages] Loaded graph: {} nodes, {} edges, {} paths",
             graph.nodes.len(), graph.edges.len(), graph.paths.len());

    // Calculate parameters
    let path_index = XPIndex::from_graph(&graph);
    let mut sum_path_step_count = 0u64;
    let mut max_path_step_count = 0usize;
    let mut max_path_length = 0usize;

    for path_id in 0..graph.paths.len() {
        let step_count = path_index.get_path_step_count(path_id);
        sum_path_step_count += step_count as u64;
        max_path_step_count = max_path_step_count.max(step_count);
        max_path_length = max_path_length.max(path_index.get_path_length(path_id));
    }

    let mut params = PathSGDParams::default();
    params.nthreads = 1;
    params.progress = verbose;
    params.min_term_updates = sum_path_step_count;
    params.eta_max = (max_path_step_count * max_path_step_count) as f64;
    params.space = max_path_length as u64;

    eprintln!("[test_ygs_stages] SGD parameters:");
    eprintln!("  eta_max: {}", params.eta_max);
    eprintln!("  min_term_updates: {}", params.min_term_updates);
    eprintln!("  space: {}", params.space);

    // Stage 1: Y (SGD only)
    eprintln!("\n[test_ygs_stages] === Stage 1: Y (SGD) ===");
    let sgd_ordering = path_sgd_sort(&graph, params);
    graph.apply_ordering(sgd_ordering, verbose);

    let mut output = fs::File::create("stage1_Y.gfa")?;
    graph.write_gfa(&mut output)?;
    eprintln!("[test_ygs_stages] Written stage1_Y.gfa");

    // Stage 2: Yg (SGD + groom)
    eprintln!("\n[test_ygs_stages] === Stage 2: g (groom) ===");
    let groomed_order = graph.groom(true, verbose);
    graph.apply_grooming_with_reorder(groomed_order, false, verbose);

    let mut output = fs::File::create("stage2_Yg.gfa")?;
    graph.write_gfa(&mut output)?;
    eprintln!("[test_ygs_stages] Written stage2_Yg.gfa");

    // Stage 3: Ygs (SGD + groom + topo)
    eprintln!("\n[test_ygs_stages] === Stage 3: s (topological sort) ===");
    let topo_order = graph.exact_odgi_topological_order(true, false, verbose);
    graph.apply_ordering(topo_order, verbose);

    let mut output = fs::File::create("stage3_Ygs.gfa")?;
    graph.write_gfa(&mut output)?;
    eprintln!("[test_ygs_stages] Written stage3_Ygs.gfa");

    eprintln!("\n[test_ygs_stages] Complete!");

    Ok(())
}
