/// Write BidirectedGraph directly to GFA without HashGraph conversion
use crate::bidirected_graph::Handle;
use crate::path_sgd_exact::{path_sgd_sort, PathSGDParams};
use crate::seqrush::SeqRush;

impl SeqRush {
    /// Write BidirectedGraph directly to GFA
    pub fn write_bidirected_gfa(
        &self,
        output_path: &str,
        no_compact: bool,
        no_sort: bool,
        groom: bool,
        sort_groom_sort: bool,
        iterative_groom: Option<usize>,
        odgi_style_groom: bool,
        sgd_sort: bool,
        sgd_iter_max: u64,
        threads: usize,
        verbose: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Build the BidirectedGraph
        // When using SGD, skip initial topological sort since SGD will provide the ordering
        let should_sort_initially = !no_sort && !sgd_sort;
        if verbose {
            if sgd_sort {
                eprintln!("[bidirected_gfa] Building BidirectedGraph without initial sorting (SGD will sort)...");
            } else if no_sort {
                eprintln!("[bidirected_gfa] Building BidirectedGraph without sorting...");
            } else {
                eprintln!("[bidirected_gfa] Building BidirectedGraph with exact ODGI sort...");
            }
        }
        let mut bi_graph = self.build_bidirected_graph_with_options(verbose, should_sort_initially)?;
        
        if verbose {
            eprintln!("[bidirected_gfa] BidirectedGraph has {} nodes (sorted)", bi_graph.nodes.len());
            eprintln!("[bidirected_gfa] After build - has node 1? {}, has node 2? {}",
                bi_graph.nodes.contains_key(&1), bi_graph.nodes.contains_key(&2));
        }

        // Compact BEFORE sorting if using SGD (to match ODGI/seqwish behavior)
        if sgd_sort && !no_compact {
            if verbose {
                eprintln!("[bidirected_gfa] Compacting graph BEFORE SGD (like seqwish does)...");
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

        // Apply grooming and/or sorting strategies
        if sgd_sort {
            // Path-guided SGD sorting (like odgi sort -p Ygs)
            if verbose {
                eprintln!("[bidirected_gfa] Using path-guided SGD sorting (like odgi sort -p Ygs)");
                eprintln!("[bidirected_gfa] Graph has {} nodes before SGD", bi_graph.nodes.len());
            }

            // Only apply SGD if we have nodes
            if !bi_graph.nodes.is_empty() {
                // Calculate parameters from graph structure (like ODGI does)
                use crate::path_sgd_exact::XPIndex;
                let path_index = XPIndex::from_graph(&bi_graph);

                let mut sum_path_step_count = 0u64;
                let mut max_path_step_count = 0usize;
                let mut max_path_length = 0usize;

                for path_id in 0..bi_graph.paths.len() {
                    let step_count = path_index.get_path_step_count(path_id);
                    sum_path_step_count += step_count as u64;
                    max_path_step_count = max_path_step_count.max(step_count);
                    max_path_length = max_path_length.max(path_index.get_path_length(path_id));
                }

                // Run exact ODGI path-guided SGD with calculated parameters
                let mut params = PathSGDParams::default();
                params.iter_max = sgd_iter_max;
                params.nthreads = threads;
                params.progress = verbose;
                params.min_term_updates = sum_path_step_count;
                params.eta_max = (max_path_step_count * max_path_step_count) as f64;
                params.space = max_path_length as u64;

                // Calculate space_quantization_step dynamically like ODGI does
                // This ensures we have approximately 100-102 zipf distributions instead of 488
                const MAX_NUMBER_OF_ZIPF_DISTRIBUTIONS: u64 = 100;
                let space_max = params.space_max;
                // ODGI uses max(space_max + 1, MAX_NUMBER_OF_ZIPF_DISTRIBUTIONS)
                let max_num_distributions = (space_max + 1).max(MAX_NUMBER_OF_ZIPF_DISTRIBUTIONS);

                if params.space > space_max && max_num_distributions > space_max {
                    // Dynamic calculation to achieve approximately MAX_NUMBER_OF_ZIPF_DISTRIBUTIONS
                    params.space_quantization_step = 2u64.max(
                        ((params.space - space_max) as f64 / (max_num_distributions - space_max) as f64).ceil() as u64
                    );
                } else {
                    // Fallback to ODGI's default
                    params.space_quantization_step = 100;
                }

                if verbose {
                    eprintln!("[bidirected_gfa] Calculated SGD parameters:");
                    eprintln!("  sum_path_step_count: {}", sum_path_step_count);
                    eprintln!("  max_path_step_count: {}", max_path_step_count);
                    eprintln!("  max_path_length: {}", max_path_length);
                    eprintln!("  min_term_updates: {}", params.min_term_updates);
                    eprintln!("  eta_max: {}", params.eta_max);
                    eprintln!("  space: {}", params.space);
                    eprintln!("  space_max: {}", params.space_max);
                    eprintln!("  space_quantization_step: {} (dynamically calculated)", params.space_quantization_step);
                    eprintln!("[bidirected_gfa] Running exact ODGI path_linear_sgd with {} iterations", params.iter_max);
                }

                // Apply full Ygs pipeline: Y (SGD) + g (groom) + s (topological sort)

                // Step 1: Y - Path-guided SGD
                let sorted_handles = path_sgd_sort(&bi_graph, params);

                if verbose {
                    eprintln!("[bidirected_gfa] Step 1/3: Path SGD produced ordering of {} nodes", sorted_handles.len());
                }

                // Convert to forward handles only for node ordering
                let forward_order: Vec<Handle> = sorted_handles.iter()
                    .map(|h| Handle::new(h.node_id(), false))
                    .collect();

                // Apply the SGD ordering
                bi_graph.apply_ordering(forward_order, verbose);

                // Validate after SGD ordering
                eprintln!("[VALIDATION] After SGD ordering:");
                bi_graph.validate_paths("after SGD ordering");

                // DEBUG: Write stage 1 output
                if verbose {
                    let stage1_file = std::fs::File::create("stage1_seqrush_Y.gfa").unwrap();
                    let mut stage1_writer = std::io::BufWriter::new(stage1_file);
                    bi_graph.write_gfa(&mut stage1_writer).unwrap();
                    eprintln!("[bidirected_gfa] DEBUG: Written stage1_seqrush_Y.gfa");
                }

                // Step 2: g - Grooming
                if verbose {
                    eprintln!("[bidirected_gfa] Step 2/3: Applying grooming after SGD...");
                }
                let groomed_order = bi_graph.groom(true, verbose);
                bi_graph.apply_grooming_with_reorder(groomed_order, false, verbose);

                // Validate after grooming
                eprintln!("[VALIDATION] After grooming:");
                bi_graph.validate_paths("after grooming");

                // DEBUG: Write stage 2 output
                if verbose {
                    let stage2_file = std::fs::File::create("stage2_seqrush_Yg.gfa").unwrap();
                    let mut stage2_writer = std::io::BufWriter::new(stage2_file);
                    bi_graph.write_gfa(&mut stage2_writer).unwrap();
                    eprintln!("[bidirected_gfa] DEBUG: Written stage2_seqrush_Yg.gfa");
                }

                // CRITICAL: Renumber nodes sequentially after grooming to preserve Y+g ordering
                if verbose {
                    eprintln!("[bidirected_gfa] Renumbering nodes sequentially to preserve Y+g order...");
                }
                bi_graph.renumber_nodes_sequentially();

                // NOTE: Skipping topological sort reordering because it would destroy the Y+g layout
                // The topological sort was causing massive jumps in node connectivity because it
                // traverses the graph topology (following edges) and then renumbers nodes based on
                // that traversal, which doesn't respect the linear Y+g layout.
                // TODO: Investigate if we need a different kind of topological adjustment that
                // only fixes orientations or local ordering without global reordering
            } else if verbose {
                eprintln!("[bidirected_gfa] Skipping SGD - no nodes in graph yet");
            }

        } else if let Some(max_iterations) = iterative_groom {
            // Iterative grooming
            if verbose {
                eprintln!("[bidirected_gfa] Using iterative grooming (max {} iterations)", max_iterations);
            }
            let iterations = bi_graph.iterative_groom(max_iterations, verbose);
            if verbose {
                eprintln!("[bidirected_gfa] Iterative grooming completed in {} iterations", iterations);
            }
        } else if odgi_style_groom {
            // ODGI-style grooming (orientation flips only, no reordering)
            if verbose {
                eprintln!("[bidirected_gfa] Using ODGI-style grooming");
            }
            let groomed_order = bi_graph.groom(true, verbose);  // Use BFS like ODGI
            bi_graph.apply_grooming_with_reorder(groomed_order, false, verbose);  // No reordering, just flips
        } else if sort_groom_sort {
            // Strategy: sort -> groom -> sort
            if verbose {
                eprintln!("[bidirected_gfa] Using sort-groom-sort strategy");
            }
            bi_graph.sort_groom_sort(verbose);
        } else if groom {
            // Strategy: groom -> sort
            if verbose {
                eprintln!("[bidirected_gfa] Using groom-then-sort strategy");
            }
            bi_graph.groom_and_sort(verbose);
        }
        // Note: Regular sorting is already done in build_bidirected_graph

        // Apply compaction if requested (but skip if we already compacted before SGD)
        let already_compacted = sgd_sort && !no_compact;
        if !no_compact && !already_compacted {
            if verbose {
                eprintln!("[bidirected_gfa] Applying compaction...");
                let nodes_before = bi_graph.nodes.len();
                bi_graph.compact();
                let nodes_after = bi_graph.nodes.len();
                eprintln!("[bidirected_gfa] Compacted from {} to {} nodes", nodes_before, nodes_after);

                // Renumber nodes to be sequential for ODGI compatibility
                bi_graph.renumber_nodes_sequentially();

                // Verify edges after compaction
                bi_graph.verify_path_edges(verbose);

                // Re-apply sorting/grooming after compaction to fix node numbering
                if sgd_sort {
                    // For SGD, only re-apply grooming, not sorting (preserve SGD order)
                    eprintln!("[bidirected_gfa] Re-applying grooming after compaction (preserving SGD order)...");
                    let groomed_order = bi_graph.groom(true, false);
                    bi_graph.apply_grooming_with_reorder(groomed_order, false, false);
                } else if sort_groom_sort {
                    eprintln!("[bidirected_gfa] Re-applying sort-groom-sort after compaction...");
                    bi_graph.sort_groom_sort(false);
                } else if odgi_style_groom {
                    eprintln!("[bidirected_gfa] Re-applying ODGI-style grooming after compaction...");
                    let groomed_order = bi_graph.groom(true, false);
                    bi_graph.apply_grooming_with_reorder(groomed_order, false, false);
                } else if groom {
                    eprintln!("[bidirected_gfa] Re-applying groom-and-sort after compaction...");
                    bi_graph.groom_and_sort(false);
                } else {
                    eprintln!("[bidirected_gfa] Re-applying topological sort after compaction...");
                    bi_graph.apply_exact_odgi_ordering(false);
                }

                // Verify edges again after sorting
                bi_graph.verify_path_edges(verbose);
            } else {
                bi_graph.compact();
                bi_graph.renumber_nodes_sequentially();
                bi_graph.verify_path_edges(false);

                // Re-apply sorting/grooming after compaction
                if sgd_sort {
                    // For SGD, only re-apply grooming, not sorting (preserve SGD order)
                    let groomed_order = bi_graph.groom(true, false);
                    bi_graph.apply_grooming_with_reorder(groomed_order, false, false);
                } else if sort_groom_sort {
                    bi_graph.sort_groom_sort(false);
                } else if odgi_style_groom {
                    let groomed_order = bi_graph.groom(true, false);
                    bi_graph.apply_grooming_with_reorder(groomed_order, false, false);
                } else if groom {
                    bi_graph.groom_and_sort(false);
                } else {
                    bi_graph.apply_exact_odgi_ordering(false);
                }
                bi_graph.verify_path_edges(false);
            }
        }
        
        // Write GFA directly from BidirectedGraph
        if verbose {
            eprintln!("[bidirected_gfa] Writing GFA to {}", output_path);
            eprintln!("[bidirected_gfa] Before write - has node 1? {}, has node 2? {}", 
                bi_graph.nodes.contains_key(&1), bi_graph.nodes.contains_key(&2));
            let mut sorted_ids: Vec<usize> = bi_graph.nodes.keys().cloned().collect();
            sorted_ids.sort();
            eprintln!("[bidirected_gfa] First 10 sorted node IDs: {:?}", &sorted_ids[..10.min(sorted_ids.len())]);
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