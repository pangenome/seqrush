/// Write BidirectedGraph directly to GFA without HashGraph conversion
use crate::ygs_sort::{YgsParams, ygs_sort, sgd_sort_only, groom_only, topological_sort_only};
use crate::seqrush::SeqRush;

impl SeqRush {
    /// Write BidirectedGraph directly to GFA
    /// By default, uses the full Ygs pipeline (Y=SGD, g=groom, s=topological sort)
    pub fn write_bidirected_gfa(
        &self,
        output_path: &str,
        no_compact: bool,
        no_sort: bool,
        skip_sgd: bool,
        skip_groom: bool,
        skip_topo: bool,
        sgd_iter_max: u64,
        sgd_eta_max: Option<f64>,
        sgd_theta: f64,
        sgd_eps: f64,
        sgd_cooling_start: f64,
        threads: usize,
        verbose: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Build the BidirectedGraph without initial sorting (Ygs will sort it)
        if verbose {
            if no_sort {
                eprintln!("[bidirected_gfa] Building BidirectedGraph without sorting...");
            } else {
                eprintln!("[bidirected_gfa] Building BidirectedGraph (Ygs pipeline will sort)...");
            }
        }
        let mut bi_graph = self.build_bidirected_graph_with_options(verbose, false)?;

        if verbose {
            eprintln!("[bidirected_gfa] BidirectedGraph has {} nodes before sorting", bi_graph.nodes.len());
        }

        // Compact BEFORE sorting (to match ODGI/seqwish behavior)
        if !no_compact {
            if verbose {
                eprintln!("[bidirected_gfa] Compacting graph before sorting...");
                let nodes_before = bi_graph.nodes.len();
                bi_graph.compact();
                let nodes_after = bi_graph.nodes.len();
                eprintln!("[bidirected_gfa] Compacted from {} to {} nodes", nodes_before, nodes_after);
                bi_graph.renumber_nodes_sequentially();
            } else {
                bi_graph.compact();
                bi_graph.renumber_nodes_sequentially();
            }
        }

        // Apply sorting if requested
        if !no_sort && !bi_graph.nodes.is_empty() {
            // Determine which phases to run
            let use_sgd = !skip_sgd;
            let use_groom = !skip_groom;
            let use_topo = !skip_topo;

            if use_sgd || use_groom || use_topo {
                if verbose {
                    eprintln!("[bidirected_gfa] Applying Ygs pipeline:");
                    eprintln!("  Y (path-guided SGD): {}", if use_sgd { "ENABLED" } else { "SKIP" });
                    eprintln!("  g (grooming): {}", if use_groom { "ENABLED" } else { "SKIP" });
                    eprintln!("  s (topological sort): {}", if use_topo { "ENABLED" } else { "SKIP" });
                }

                // Create YGS parameters
                let mut params = YgsParams::from_graph(&bi_graph, verbose, threads);
                params.path_sgd.iter_max = sgd_iter_max;
                params.path_sgd.theta = sgd_theta;
                params.path_sgd.eps = sgd_eps;
                params.path_sgd.cooling_start = sgd_cooling_start;

                // Override eta_max if specified
                if let Some(eta_max) = sgd_eta_max {
                    params.path_sgd.eta_max = eta_max;
                }

                if verbose {
                    eprintln!("[bidirected_gfa] SGD parameters:");
                    eprintln!("  iter_max: {}", params.path_sgd.iter_max);
                    eprintln!("  eta_max: {}", params.path_sgd.eta_max);
                    eprintln!("  theta: {}", params.path_sgd.theta);
                    eprintln!("  eps: {}", params.path_sgd.eps);
                    eprintln!("  cooling_start: {}", params.path_sgd.cooling_start);
                    eprintln!("  nthreads: {}", params.path_sgd.nthreads);
                }

                // Apply selected phases
                if use_sgd && use_groom && use_topo {
                    // Use full Ygs pipeline
                    ygs_sort(&mut bi_graph, &params);
                } else {
                    // Apply individual phases
                    if use_sgd {
                        if verbose {
                            eprintln!("[bidirected_gfa] Applying path-guided SGD (Y)...");
                        }
                        sgd_sort_only(&mut bi_graph, params.path_sgd.clone(), verbose);
                    }

                    if use_groom {
                        if verbose {
                            eprintln!("[bidirected_gfa] Applying grooming (g)...");
                        }
                        groom_only(&mut bi_graph, verbose);
                    }

                    if use_topo {
                        if verbose {
                            eprintln!("[bidirected_gfa] Applying topological sort (s)...");
                        }
                        topological_sort_only(&mut bi_graph, verbose);
                    }
                }
            }
        } else if verbose {
            if bi_graph.nodes.is_empty() {
                eprintln!("[bidirected_gfa] Skipping sorting - no nodes in graph");
            } else {
                eprintln!("[bidirected_gfa] Skipping sorting (--no-sort specified)");
            }
        }

        // Write GFA directly from BidirectedGraph
        if verbose {
            eprintln!("[bidirected_gfa] Writing GFA to {}", output_path);
            let mut sorted_ids: Vec<usize> = bi_graph.nodes.keys().cloned().collect();
            sorted_ids.sort();
            eprintln!("[bidirected_gfa] Final node count: {}", bi_graph.nodes.len());
            if !sorted_ids.is_empty() {
                eprintln!("[bidirected_gfa] Node ID range: {} to {}", sorted_ids[0], sorted_ids[sorted_ids.len() - 1]);
            }
        }

        let file = std::fs::File::create(output_path)?;
        let mut writer = std::io::BufWriter::new(file);

        // Validate that all paths correctly represent their sequences before writing
        if let Err(e) = self.validate_paths_match_sequences(&bi_graph) {
            eprintln!("ERROR: Path validation failed!");
            eprintln!("{}", e);
            return Err(e);
        }

        if verbose {
            println!("✓ All paths correctly represent their input sequences");
        }

        bi_graph.write_gfa(&mut writer)?;

        if verbose {
            println!(
                "BidirectedGraph written to {}: {} nodes, {} edges, {} paths",
                output_path,
                bi_graph.nodes.len(),
                bi_graph.edges.len(),
                bi_graph.paths.len()
            );
        }

        Ok(())
    }
}
