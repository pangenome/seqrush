/// Write BidirectedGraph directly to GFA without HashGraph conversion

use crate::bidirected_ops::BidirectedGraph;
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
        verbose: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Build the BidirectedGraph
        if verbose {
            if no_sort {
                eprintln!("[bidirected_gfa] Building BidirectedGraph without sorting...");
            } else {
                eprintln!("[bidirected_gfa] Building BidirectedGraph with exact ODGI sort...");
            }
        }
        let mut bi_graph = self.build_bidirected_graph_with_options(verbose, !no_sort)?;
        
        if verbose {
            eprintln!("[bidirected_gfa] BidirectedGraph has {} nodes (sorted)", bi_graph.nodes.len());
            eprintln!("[bidirected_gfa] After build - has node 1? {}, has node 2? {}", 
                bi_graph.nodes.contains_key(&1), bi_graph.nodes.contains_key(&2));
        }
        
        // Apply grooming and/or sorting strategies
        if let Some(max_iterations) = iterative_groom {
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
                if sort_groom_sort {
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
                if sort_groom_sort {
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