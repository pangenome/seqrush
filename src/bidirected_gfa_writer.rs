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
        _no_sort: bool,  // Sorting already done in build_bidirected_graph
        verbose: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Build the BidirectedGraph (which includes exact ODGI topological sort)
        if verbose {
            eprintln!("[bidirected_gfa] Building BidirectedGraph with exact ODGI sort...");
        }
        let mut bi_graph = self.build_bidirected_graph(verbose)?;
        
        if verbose {
            eprintln!("[bidirected_gfa] BidirectedGraph has {} nodes (sorted)", bi_graph.nodes.len());
            eprintln!("[bidirected_gfa] After build - has node 1? {}, has node 2? {}", 
                bi_graph.nodes.contains_key(&1), bi_graph.nodes.contains_key(&2));
        }
        
        // Apply compaction if requested
        if !no_compact {
            if verbose {
                eprintln!("[bidirected_gfa] Applying compaction...");
                let nodes_before = bi_graph.nodes.len();
                bi_graph.compact();
                let nodes_after = bi_graph.nodes.len();
                eprintln!("[bidirected_gfa] Compacted from {} to {} nodes", nodes_before, nodes_after);
                
                // Re-apply topological sort after compaction to fix node numbering
                eprintln!("[bidirected_gfa] Re-applying topological sort after compaction...");
                bi_graph.apply_exact_odgi_ordering(false);
            } else {
                bi_graph.compact();
                // Re-apply topological sort after compaction to fix node numbering
                bi_graph.apply_exact_odgi_ordering(false);
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