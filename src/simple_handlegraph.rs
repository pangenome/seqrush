use handlegraph::handle::{Edge, Handle, NodeId};
use handlegraph::handlegraph::{HandleGraph, IntoEdges, IntoHandles, IntoSequences};
use handlegraph::hashgraph::HashGraph;
use handlegraph::mutablehandlegraph::AdditiveHandleGraph;
use handlegraph::pathhandlegraph::embedded_paths::{
    GraphPathNames, GraphPaths, GraphPathsRef, IntoPathIds, MutableGraphPaths,
};
use handlegraph::pathhandlegraph::path::PathSteps;

use crate::pos::{make_pos, offset};
use crate::seqrush::SeqRush;
use std::collections::{HashMap, HashSet, VecDeque};

impl SeqRush {
    /// Build and write handlegraph GFA
    pub fn write_handlegraph_gfa(
        &self,
        output_path: &str,
        no_compact: bool,
        no_sort: bool,
        verbose: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Build the graph
        let mut graph = HashGraph::new();

        // Track union to node mapping
        let mut union_to_node: HashMap<u64, NodeId> = HashMap::new();
        let mut next_node_id = 1u64;

        // Create nodes
        let mut unique_unions = HashSet::new();
        for seq in &self.sequences {
            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                let pos_fwd = make_pos(global_pos, false);
                let union_fwd = self.union_find.find(pos_fwd);
                unique_unions.insert(union_fwd);
            }
        }

        let mut debug_count = 0;
        for &union_rep in &unique_unions {
            let base = self.get_base_for_union(union_rep);
            let node_id = NodeId::from(next_node_id);

            if verbose && debug_count < 5 {
                println!(
                    "  Creating node {} for union_rep {}",
                    next_node_id, union_rep
                );
                debug_count += 1;
            }

            next_node_id += 1;

            graph.create_handle(&[base], node_id);
            union_to_node.insert(union_rep, node_id);
        }

        if verbose {
            println!(
                "Created {} nodes with IDs 1..{}",
                unique_unions.len(),
                next_node_id
            );
        }

        // Build paths with correct orientation per step (seqwish-style)
        let mut all_edges = Vec::new();

        for seq in &self.sequences {
            let path_id = graph.create_path(seq.id.as_bytes(), false).unwrap();
            let mut prev_handle: Option<Handle> = None;

            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;
                let pos_fwd = make_pos(global_pos, false);
                let pos_rev = make_pos(global_pos, true);
                let union_fwd = self.union_find.find(pos_fwd);
                let union_rev = self.union_find.find(pos_rev);

                // Determine component representative and orientation for this step
                let mut union_rep = union_fwd;
                let mut is_reverse = false;
                // Prefer direct key hits first
                if !union_to_node.contains_key(&union_rep) {
                    if union_to_node.contains_key(&union_rev) {
                        union_rep = union_rev;
                        is_reverse = true;
                    } else {
                        // Fall back: scan existing reps in same component
                        let mut found = None;
                        for &existing in union_to_node.keys() {
                            if self.union_find.same(union_fwd, existing) {
                                found = Some((existing, false));
                                break;
                            }
                            if self.union_find.same(union_rev, existing) {
                                found = Some((existing, true));
                                break;
                            }
                        }
                        if let Some((existing, rev)) = found {
                            union_rep = existing;
                            is_reverse = rev;
                        } else {
                            panic!(
                                "No node found for union component of pos {} (fwd={}, rev={})",
                                global_pos, union_fwd, union_rev
                            );
                        }
                    }
                }

                let node_id = union_to_node[&union_rep];
                let handle = Handle::pack(node_id, is_reverse);

                if verbose && i < 5 {
                    println!(
                        "  Seq {} pos {}: union_rep={}, node_id={}, orient={}",
                        seq.id,
                        i,
                        union_rep,
                        node_id.0,
                        if is_reverse { '-' } else { '+' }
                    );
                }

                graph.path_append_step(path_id, handle);

                if let Some(prev) = prev_handle {
                    all_edges.push(Edge(prev, handle));
                }
                prev_handle = Some(handle);
            }
        }

        // Add edges
        let mut seen_edges = HashSet::new();
        let mut edge_count = 0;
        for edge in &all_edges {
            let key = (edge.0.as_integer(), edge.1.as_integer());
            if seen_edges.insert(key) {
                graph.create_edge(*edge);
                edge_count += 1;

                if verbose && edge_count <= 5 {
                    println!(
                        "  Edge {}: {} -> {} (handles: {:?} -> {:?})",
                        edge_count,
                        edge.0.id().0,
                        edge.1.id().0,
                        edge.0,
                        edge.1
                    );
                }
            }
        }

        if verbose {
            println!(
                "Built graph: {} nodes, {} edges, {} paths",
                graph.node_count(),
                graph.edge_count(),
                graph.path_count()
            );
        }

        // Apply compaction if requested
        let compacted_graph = if !no_compact {
            if verbose {
                eprintln!("[unchop] Applying graph compaction...");
            }
            match crate::bidirected_unchop::bidirected_unchop(&graph, verbose) {
                Ok(compacted) => compacted,
                Err(e) => {
                    eprintln!("Warning: Compaction failed: {}", e);
                    graph
                }
            }
        } else {
            graph
        };

        // Apply topological sorting if requested (default is to sort)
        let final_graph = if !no_sort {
            if verbose {
                eprintln!("[sort] Applying topological sort to renumber nodes 1..N...");
            }
            match self.topological_sort_handlegraph(&compacted_graph, verbose) {
                Ok(sorted) => sorted,
                Err(e) => {
                    eprintln!("Warning: Topological sort failed: {}", e);
                    compacted_graph
                }
            }
        } else {
            compacted_graph
        };

        // Write GFA
        let file = std::fs::File::create(output_path)?;
        let mut writer = std::io::BufWriter::new(file);
        use std::io::Write;

        writeln!(writer, "H\tVN:Z:1.0")?;

        // Write nodes in sorted order (by node ID)
        let mut node_handles: Vec<_> = (&final_graph).handles().collect();
        node_handles.sort_by_key(|h| h.id().0);
        
        let mut node_count = 0;
        for handle in node_handles {
            let node_id = handle.id();
            let seq: Vec<u8> = (&final_graph).sequence(handle).collect();
            writeln!(
                writer,
                "S\t{}\t{}",
                node_id.0,
                String::from_utf8_lossy(&seq)
            )?;
            node_count += 1;

            if verbose && node_count <= 5 {
                println!("  Writing node {}: handle={:?}", node_id.0, handle);
            }
        }

        // Write edges
        for edge in (&final_graph).edges() {
            writeln!(
                writer,
                "L\t{}\t{}\t{}\t{}\t0M",
                u64::from(edge.0.id()),
                if edge.0.is_reverse() { '-' } else { '+' },
                u64::from(edge.1.id()),
                if edge.1.is_reverse() { '-' } else { '+' }
            )?;
        }

        // Write paths from the final graph in deterministic (name) order
        let mut path_ids: Vec<_> = (&final_graph).path_ids().collect();
        path_ids.sort_by_key(|pid| {
            let name: Vec<u8> = (&final_graph).get_path_name(*pid).unwrap().collect();
            String::from_utf8_lossy(&name).to_string()
        });
        for path_id in path_ids {
            let path_name: Vec<u8> = (&final_graph).get_path_name(path_id).unwrap().collect();
            let mut path_str = String::new();
            let mut step_idx = 0;

            if let Some(path_ref) = (&final_graph).get_path_ref(path_id) {
                for step in path_ref.steps() {
                    let handle = step.1; // Step is a tuple (StepIx, Handle)
                    if step_idx > 0 {
                        path_str.push(',');
                    }
                    path_str.push_str(&format!(
                        "{}{}",
                        u64::from(handle.id()),
                        if handle.is_reverse() { '-' } else { '+' }
                    ));
                    step_idx += 1;
                }
            }

            writeln!(
                writer,
                "P\t{}\t{}\t*",
                String::from_utf8_lossy(&path_name),
                path_str
            )?;
        }

        println!(
            "Handlegraph written to {}: {} nodes, {} edges, {} paths",
            output_path,
            final_graph.node_count(),
            final_graph.edge_count(),
            final_graph.path_count()
        );

        Ok(())
    }

    fn get_base_for_union(&self, union_rep: u64) -> u8 {
        let offset_val = offset(union_rep);
        for seq in &self.sequences {
            if offset_val >= seq.offset && offset_val < seq.offset + seq.data.len() {
                return seq.data[offset_val - seq.offset];
            }
        }
        b'N'
    }

    /// Apply topological sort to a handlegraph, renumbering nodes from 1 to N
    fn topological_sort_handlegraph(
        &self,
        graph: &HashGraph,
        verbose: bool,
    ) -> Result<HashGraph, Box<dyn std::error::Error>> {
        if verbose {
            eprintln!("[sort] Starting optimized topological sort...");
        }
        
        // Collect all nodes
        let mut all_nodes: Vec<NodeId> = graph.handles().map(|h| h.id()).collect();
        
        // Build adjacency information for both directions
        let mut forward_edges: HashMap<NodeId, Vec<NodeId>> = HashMap::new();
        let mut backward_edges: HashMap<NodeId, Vec<NodeId>> = HashMap::new();
        
        for handle in graph.handles() {
            let node_id = handle.id();
            forward_edges.insert(node_id, Vec::new());
            backward_edges.insert(node_id, Vec::new());
        }
        
        // Build edge lists considering all edge types
        for edge in graph.edges() {
            let from_id = edge.0.id();
            let to_id = edge.1.id();
            
            // Track forward edges (ignoring orientation for now)
            forward_edges.entry(from_id).or_default().push(to_id);
            backward_edges.entry(to_id).or_default().push(from_id);
        }
        
        // Find nodes that appear in paths and use them as anchors
        let mut path_nodes = HashSet::new();
        let mut path_order = Vec::new();
        
        for path_id in graph.path_ids() {
            if let Some(path_ref) = graph.get_path_ref(path_id) {
                for step in path_ref.steps() {
                    let node_id = step.1.id();
                    if path_nodes.insert(node_id) {
                        path_order.push(node_id);
                    }
                }
            }
        }
        
        // Start with path nodes in the order they appear
        let mut sorted_order = Vec::new();
        let mut visited = HashSet::new();
        
        // Add path nodes first
        for node_id in path_order {
            if visited.insert(node_id) {
                sorted_order.push(node_id);
            }
        }
        
        // Add remaining nodes based on connectivity
        let mut remaining: Vec<_> = all_nodes
            .into_iter()
            .filter(|id| !visited.contains(id))
            .collect();
        
        // Sort remaining nodes by their connectivity to already sorted nodes
        remaining.sort_by_key(|&node_id| {
            let mut min_distance = usize::MAX;
            
            // Find minimum distance to any sorted node
            for &sorted_id in &sorted_order {
                // Check if directly connected
                if forward_edges.get(&sorted_id).map_or(false, |v| v.contains(&node_id)) ||
                   backward_edges.get(&sorted_id).map_or(false, |v| v.contains(&node_id)) {
                    min_distance = min_distance.min(1);
                }
            }
            
            min_distance
        });
        
        sorted_order.extend(remaining);
        
        if verbose {
            eprintln!("[sort] Sorted {} nodes", sorted_order.len());
            eprintln!("[sort] Original graph had {} nodes", graph.node_count());
            if sorted_order.len() != graph.node_count() {
                eprintln!("[sort] WARNING: Node count mismatch!");
            }
        }
        
        // Create mapping from old IDs to new IDs (1 to N)
        let mut old_to_new: HashMap<NodeId, NodeId> = HashMap::new();
        for (new_idx, &old_id) in sorted_order.iter().enumerate() {
            old_to_new.insert(old_id, NodeId::from((new_idx + 1) as u64));
        }
        
        // Build new graph with renumbered nodes
        let mut new_graph = HashGraph::new();
        
        // Create nodes with new IDs
        for (idx, &old_id) in sorted_order.iter().enumerate() {
            let new_id = NodeId::from((idx + 1) as u64);
            let old_handle = Handle::pack(old_id, false);
            let seq: Vec<u8> = graph.sequence(old_handle).collect();
            let created_handle = new_graph.create_handle(&seq, new_id);
            if verbose && idx < 5 {
                eprintln!("[sort] Created node with intended ID {} -> actual ID {}", 
                         new_id.0, created_handle.id().0);
            }
        }
        
        // Add edges with remapped IDs
        let mut seen_edges = HashSet::new();
        for edge in graph.edges() {
            let old_from_id = edge.0.id();
            let old_to_id = edge.1.id();
            let new_from_id = old_to_new[&old_from_id];
            let new_to_id = old_to_new[&old_to_id];
            
            let new_edge = Edge(
                Handle::pack(new_from_id, edge.0.is_reverse()),
                Handle::pack(new_to_id, edge.1.is_reverse())
            );
            
            // Avoid duplicate edges
            let edge_key = (new_edge.0.as_integer(), new_edge.1.as_integer());
            if seen_edges.insert(edge_key) {
                new_graph.create_edge(new_edge);
            }
        }
        
        // Copy paths with remapped node IDs
        for path_id in graph.path_ids() {
            let path_name: Vec<u8> = graph.get_path_name(path_id).unwrap().collect();
            let new_path_id = new_graph.create_path(&path_name, false).unwrap();
            
            if let Some(path_ref) = graph.get_path_ref(path_id) {
                for step in path_ref.steps() {
                    let handle = step.1; // Step is a tuple (StepIx, Handle)
                    let old_id = handle.id();
                    let new_id = if let Some(&mapped_id) = old_to_new.get(&old_id) {
                        mapped_id
                    } else {
                        if verbose {
                            eprintln!("[sort] WARNING: Node {} not found in mapping!", old_id.0);
                        }
                        old_id // Keep original if not found
                    };
                    let new_handle = Handle::pack(new_id, handle.is_reverse());
                    new_graph.path_append_step(new_path_id, new_handle);
                }
            }
        }
        
        if verbose {
            eprintln!(
                "[sort] Created sorted graph: {} nodes, {} edges, {} paths",
                new_graph.node_count(),
                new_graph.edge_count(),
                new_graph.path_count()
            );
            eprintln!("[sort] Nodes renumbered from 1 to {}", sorted_order.len());
            
            // Debug: check a few node IDs to verify they're sequential
            let mut all_ids: Vec<u64> = Vec::new();
            for handle in new_graph.handles() {
                all_ids.push(handle.id().0);
            }
            all_ids.sort();
            let first_10: Vec<u64> = all_ids.iter().take(10).copied().collect();
            let last_10: Vec<u64> = all_ids.iter().rev().take(10).rev().copied().collect();
            eprintln!("[sort] First 10 node IDs (sorted): {:?}", first_10);
            eprintln!("[sort] Last 10 node IDs (sorted): {:?}", last_10);
            eprintln!("[sort] Min ID: {}, Max ID: {}", all_ids.first().unwrap_or(&0), all_ids.last().unwrap_or(&0));
        }
        
        Ok(new_graph)
    }
}
