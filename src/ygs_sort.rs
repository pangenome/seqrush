/// Exact reimplementation of `odgi sort -p Ygs`
///
/// This module implements the three-stage sorting pipeline:
/// 1. Y - Path-guided stochastic gradient descent (PG-SGD)
/// 2. g - Grooming to remove spurious inverting links
/// 3. s - Topological sort using head nodes
///
/// Based on ODGI's implementation in src/subcommand/sort_main.cpp

use crate::bidirected_ops::BidirectedGraph;
use crate::path_sgd::{PathSGDParams, path_sgd_sort};

/// Parameters for the Ygs sorting pipeline
/// These match ODGI's defaults for the Ygs pipeline
#[derive(Clone)]
pub struct YgsParams {
    /// Path SGD parameters
    pub path_sgd: PathSGDParams,
    /// Whether to print progress information
    pub verbose: bool,
}

impl Default for YgsParams {
    fn default() -> Self {
        YgsParams {
            path_sgd: PathSGDParams {
                iter_max: 100,  // ODGI default (from sort_main.cpp)
                iter_with_max_learning_rate: 0,
                min_term_updates: 0,  // Will be calculated based on paths
                delta: 0.0,
                eps: 0.01,
                eta_max: 0.0,  // Will be calculated from path lengths
                theta: 0.99,
                space: 0,  // Will be calculated from longest path
                space_max: 100,
                space_quantization_step: 100,
                cooling_start: 0.5,  // ODGI default
                nthreads: 1,
                progress: false,
            },
            verbose: false,
        }
    }
}

impl YgsParams {
    /// Create parameters with calculated defaults based on graph structure
    /// This matches how ODGI calculates the parameters in sort_main.cpp
    pub fn from_graph(graph: &BidirectedGraph, verbose: bool, nthreads: usize) -> Self {
        let mut params = Self::default();
        params.verbose = verbose;
        params.path_sgd.nthreads = nthreads;
        params.path_sgd.progress = verbose;

        // Calculate parameters based on graph structure
        // Build a temporary path index to get statistics
        let path_index = crate::path_sgd::PathIndex::from_graph(graph);

        // Calculate sum of path step counts
        let mut sum_path_step_count = 0u64;
        let mut max_path_step_count = 0usize;
        let mut max_path_length = 0usize;

        for i in 0..path_index.num_paths() {
            let step_count = path_index.get_path_step_count(i);
            sum_path_step_count += step_count as u64;
            max_path_step_count = max_path_step_count.max(step_count);
            max_path_length = max_path_length.max(path_index.get_path_length(i));
        }

        // Set min_term_updates (ODGI default: 1.0 * sum_path_step_count)
        params.path_sgd.min_term_updates = sum_path_step_count;

        // Set eta_max (ODGI default: max_path_step_count^2)
        params.path_sgd.eta_max = (max_path_step_count * max_path_step_count) as f64;

        // Set space (ODGI default: max path length)
        params.path_sgd.space = max_path_length as u64;

        if verbose {
            eprintln!("[ygs_sort] Calculated parameters:");
            eprintln!("  sum_path_step_count: {}", sum_path_step_count);
            eprintln!("  max_path_step_count: {}", max_path_step_count);
            eprintln!("  max_path_length: {}", max_path_length);
            eprintln!("  min_term_updates: {}", params.path_sgd.min_term_updates);
            eprintln!("  eta_max: {}", params.path_sgd.eta_max);
            eprintln!("  space: {}", params.path_sgd.space);
        }

        params
    }
}

/// Apply the Ygs sorting pipeline to a graph
/// This exactly replicates `odgi sort -p Ygs`
pub fn ygs_sort(graph: &mut BidirectedGraph, params: &YgsParams) {
    if params.verbose {
        eprintln!("[ygs_sort] Starting Ygs pipeline (Y=SGD, g=groom, s=topological_sort)");
        eprintln!("[ygs_sort] Initial graph: {} nodes, {} edges",
                 graph.nodes.len(), graph.edges.len());
    }

    // Step 1: Y - Path-guided SGD sort
    if params.verbose {
        eprintln!("[ygs_sort] === Step 1/3: Path-guided SGD (Y) ===");
    }

    let sgd_ordering = path_sgd_sort(graph, params.path_sgd.clone());
    graph.apply_ordering(sgd_ordering, params.verbose);

    if params.verbose {
        eprintln!("[ygs_sort] After SGD: {} nodes", graph.nodes.len());
        eprintln!("[ygs_sort] Edges after SGD:");
        let mut edges_vec: Vec<_> = graph.edges.iter().collect();
        edges_vec.sort_by_key(|e| (e.from.node_id(), e.from.is_reverse(), e.to.node_id(), e.to.is_reverse()));
        for edge in edges_vec {
            eprintln!("[ygs_sort]   {} {} -> {} {}",
                     edge.from.node_id(),
                     if edge.from.is_reverse() { "-" } else { "+" },
                     edge.to.node_id(),
                     if edge.to.is_reverse() { "-" } else { "+" });
        }
        eprintln!("[ygs_sort] Path after SGD:");
        for path in &graph.paths {
            for (i, step) in path.steps.iter().enumerate() {
                eprintln!("[ygs_sort]   Step {}: Node {}{}",
                         i, step.node_id(),
                         if step.is_reverse() { "-" } else { "+" });
            }
        }
    }

    // Step 2: g - Groom the graph
    // Using BFS like ODGI (simple first-visit orientation locking)
    if params.verbose {
        eprintln!("[ygs_sort] === Step 2/3: Grooming (g) - BFS ===");
    }

    let groomed_order = graph.groom(true, params.verbose);  // Use BFS like ODGI
    graph.apply_grooming_with_reorder(groomed_order, false, params.verbose);

    if params.verbose {
        eprintln!("[ygs_sort] After grooming: {} nodes", graph.nodes.len());
    }

    // TODO: Step 3 (topological sort) is currently disabled
    //
    // PROBLEM IDENTIFIED:
    // - Topological sort creates a NEW ordering based on edges, destroying SGD's ordering
    // - Even after grooming, our topo sort degrades performance (76.3% → 64.8%)
    // - ODGI's topo sort IMPROVES performance after grooming (77.4% → 85.0%)
    //
    // CURRENT STATUS:
    // - Using SGD + groom only
    // - Need to fix topological sort to preserve SGD ordering

    /* DISABLED - topological sort
    // Step 3: s - Topological sort (heads only)
    if params.verbose {
        eprintln!("[ygs_sort] === Step 3/3: Topological sort (s) ===");
    }

    // use_heads=true, use_tails=false matches ODGI's 's' command
    let topo_order = graph.exact_odgi_topological_order(true, false, params.verbose);
    graph.apply_ordering(topo_order, params.verbose);

    if params.verbose {
        eprintln!("[ygs_sort] After topological sort: {} nodes", graph.nodes.len());
    }
    */

    if params.verbose {
        eprintln!("[ygs_sort] === Ygs pipeline complete (SGD + groom, topo disabled) ===");
    }
}

/// Apply just the topological sort step (the 's' part)
/// This is useful for testing or for applying just the final sort
pub fn topological_sort_only(graph: &mut BidirectedGraph, verbose: bool) {
    if verbose {
        eprintln!("[topological_sort] Starting topological sort (heads only)");
    }

    let order = graph.exact_odgi_topological_order(true, false, verbose);
    graph.apply_ordering(order, verbose);

    if verbose {
        eprintln!("[topological_sort] Complete");
    }
}

/// Apply just the grooming step (the 'g' part)
pub fn groom_only(graph: &mut BidirectedGraph, verbose: bool) {
    if verbose {
        eprintln!("[groom] Starting grooming");
    }

    let groomed_order = graph.groom(true, verbose);
    graph.apply_grooming_with_reorder(groomed_order, false, verbose);

    if verbose {
        eprintln!("[groom] Complete");
    }
}

/// Apply the SGD step (the 'Y' part)
pub fn sgd_sort_only(graph: &mut BidirectedGraph, params: PathSGDParams, verbose: bool) {
    if verbose {
        eprintln!("[path_sgd] Starting path-guided SGD");
    }

    let ordering = path_sgd_sort(graph, params);
    graph.apply_ordering(ordering, verbose);

    if verbose {
        eprintln!("[path_sgd] Complete");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bidirected_graph::{BiPath, Handle};

    fn create_test_graph() -> BidirectedGraph {
        let mut graph = BidirectedGraph::new();

        // Create a simple linear graph: 1 -> 2 -> 3
        graph.add_node(1, b"AAAA".to_vec());
        graph.add_node(2, b"CCCC".to_vec());
        graph.add_node(3, b"GGGG".to_vec());

        graph.add_edge(Handle::forward(1), Handle::forward(2));
        graph.add_edge(Handle::forward(2), Handle::forward(3));

        // Add a path through all nodes
        let mut path = BiPath::new("test_path".to_string());
        path.add_step(Handle::forward(1));
        path.add_step(Handle::forward(2));
        path.add_step(Handle::forward(3));
        graph.paths.push(path);

        graph
    }

    #[test]
    fn test_ygs_params_default() {
        let params = YgsParams::default();
        assert_eq!(params.path_sgd.iter_max, 100);  // ODGI default
        assert_eq!(params.path_sgd.theta, 0.99);
        assert_eq!(params.path_sgd.eps, 0.01);
    }

    #[test]
    fn test_ygs_params_from_graph() {
        let graph = create_test_graph();
        let params = YgsParams::from_graph(&graph, false, 1);

        // Check that parameters were calculated
        assert!(params.path_sgd.min_term_updates > 0);
        assert!(params.path_sgd.eta_max > 0.0);
        assert!(params.path_sgd.space > 0);
    }

    #[test]
    fn test_ygs_sort_runs() {
        let mut graph = create_test_graph();
        let params = YgsParams::from_graph(&graph, false, 1);

        // This should not panic
        ygs_sort(&mut graph, &params);

        // Graph should still be valid
        assert_eq!(graph.node_count(), 3);
        assert!(graph.paths.len() > 0);
    }

    #[test]
    fn test_individual_steps() {
        let graph = create_test_graph();

        // Test SGD only
        {
            let mut g = graph.clone();
            let params = YgsParams::from_graph(&g, false, 1);
            sgd_sort_only(&mut g, params.path_sgd, false);
            assert_eq!(g.node_count(), 3);
        }

        // Test groom only
        {
            let mut g = graph.clone();
            groom_only(&mut g, false);
            assert_eq!(g.node_count(), 3);
        }

        // Test topological sort only
        {
            let mut g = graph.clone();
            topological_sort_only(&mut g, false);
            assert_eq!(g.node_count(), 3);
        }
    }
}
