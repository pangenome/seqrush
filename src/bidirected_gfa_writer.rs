/// Write BidirectedGraph directly to GFA without HashGraph conversion

use crate::bidirected_ops::BidirectedGraph;
use crate::bidirected_graph::Handle;
use crate::path_sgd_exact::{path_sgd_sort, PathSGDParams};
use crate::seqrush::SeqRush;
use std::io::Write;

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
        
        // Apply grooming and/or sorting strategies
        if sgd_sort {
            // Path-guided SGD sorting (like odgi sort -p Ygs)
            if verbose {
                eprintln!("[bidirected_gfa] Using path-guided SGD sorting (like odgi sort -p Ygs)");
                eprintln!("[bidirected_gfa] Graph has {} nodes before SGD", bi_graph.nodes.len());
            }

            // Only apply SGD if we have nodes
            if !bi_graph.nodes.is_empty() {
                // Run exact ODGI path-guided SGD with same parameters as odgi sort -p Ygs
                let mut params = PathSGDParams::default();
                params.nthreads = threads;
                params.progress = verbose;

                if verbose {
                    eprintln!("[bidirected_gfa] Running exact ODGI path_linear_sgd with {} iterations", params.iter_max);
                }

                let sorted_handles = path_sgd_sort(&bi_graph, params);

                if verbose {
                    eprintln!("[bidirected_gfa] Path SGD produced ordering of {} nodes", sorted_handles.len());
                }

                // Convert to forward handles only for node ordering
                let forward_order: Vec<Handle> = sorted_handles.iter()
                    .map(|h| Handle::new(h.node_id(), false))
                    .collect();

                // Apply the SGD ordering
                bi_graph.apply_ordering(forward_order, verbose);

                // Then apply grooming (degroom) - only flip orientations, don't reorder!
                if verbose {
                    eprintln!("[bidirected_gfa] Applying grooming (orientation flips only) after SGD...");
                }
                let groomed_order = bi_graph.groom(true, false);
                bi_graph.apply_grooming_with_reorder(groomed_order, false, false);

                // NO topological sort after SGD! That would destroy the SGD ordering
            } else {
                if verbose {
                    eprintln!("[bidirected_gfa] Skipping SGD - no nodes in graph yet");
                }
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

        // Apply compaction if requested
        if !no_compact {
            if verbose {
                eprintln!("[bidirected_gfa] Applying compaction...");
                let nodes_before = bi_graph.nodes.len();
                bi_graph.compact();
                let nodes_after = bi_graph.nodes.len();
                eprintln!("[bidirected_gfa] Compacted from {} to {} nodes", nodes_before, nodes_after);

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