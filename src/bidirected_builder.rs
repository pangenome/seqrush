use crate::bidirected_graph::Handle;
use crate::bidirected_ops::BidirectedGraph;
use crate::graph_ops::{Edge, Graph, Node};
use crate::pos::{make_pos, offset, Pos};
use crate::seqrush::{SeqRush, Sequence};
use std::collections::HashMap;

impl SeqRush {
    /// Build a bidirected graph from the union-find results
    pub fn build_bidirected_graph(
        &self,
        verbose: bool,
    ) -> Result<BidirectedGraph, Box<dyn std::error::Error>> {
        self.build_bidirected_graph_with_options(verbose, true)
    }

    pub fn build_bidirected_graph_with_options(
        &self,
        verbose: bool,
        apply_sort: bool,
    ) -> Result<BidirectedGraph, Box<dyn std::error::Error>> {
        let mut graph = BidirectedGraph::new();

        // Track which union representatives we've seen and their node IDs
        let mut union_to_node: HashMap<Pos, usize> = HashMap::new();
        let mut next_node_id = 1;

        // Build paths and discover nodes
        for (seq_idx, seq) in self.sequences.iter().enumerate() {
            if verbose {
                eprintln!(
                    "Building path for sequence {} of {}: {}",
                    seq_idx + 1,
                    self.sequences.len(),
                    seq.id
                );
            }
            let mut path_handles = Vec::new();

            for i in 0..seq.data.len() {
                let global_pos = seq.offset + i;

                // For each position in the sequence, we need to determine which orientation to use
                let pos_fwd = make_pos(global_pos, false);
                let pos_rev = make_pos(global_pos, true);

                // Find which orientation belongs to which union
                let union_fwd = self.union_find.find(pos_fwd);
                let union_rev = self.union_find.find(pos_rev);

                // For this position, check which orientation (if any) has already been
                // mapped to a node.
                let union_rep = if union_to_node.contains_key(&union_fwd) {
                    // Forward orientation already has a node
                    union_fwd
                } else if union_to_node.contains_key(&union_rev) {
                    // Reverse orientation already has a node
                    union_rev
                } else {
                    // Neither exact union is in the map. But they might be united with
                    // other unions that ARE in the map. We need to find the canonical
                    // representative of this union component.

                    // Get all representatives in the union-find
                    let fwd_root = self.union_find.find(pos_fwd);
                    let rev_root = self.union_find.find(pos_rev);

                    // Check if these roots are already mapped
                    if union_to_node.contains_key(&fwd_root) {
                        fwd_root
                    } else if union_to_node.contains_key(&rev_root) {
                        rev_root
                    } else {
                        // Still not found. One more check: maybe this position's union
                        // was united with something else that's already in the map.
                        // This happens with RC alignments where pos_fwd might be united
                        // with some other_pos_rev.

                        let mut found = None;

                        if verbose && seq.id == "seq2" && i == 0 {
                            eprintln!(
                                "  Searching for existing unions connected to fwd={} or rev={}",
                                union_fwd, union_rev
                            );
                            eprintln!(
                                "  Testing: is union 8 same as union 6? {}",
                                self.union_find.same(make_pos(4, false), make_pos(3, false))
                            );
                            eprintln!(
                                "  Testing: is union 8 same as existing union 6? {}",
                                self.union_find.same(union_fwd, make_pos(3, false))
                            );
                        }

                        for &existing_union in union_to_node.keys() {
                            // Check if current position's forward union is in same component
                            if self.union_find.same(union_fwd, existing_union) {
                                if verbose && seq.id == "seq2" && i == 0 {
                                    eprintln!(
                                        "    Found! union {} is same as existing union {}",
                                        union_fwd, existing_union
                                    );
                                }
                                found = Some(existing_union);
                                break;
                            }
                            // Check if current position's reverse union is in same component
                            if self.union_find.same(union_rev, existing_union) {
                                if verbose && seq.id == "seq2" && i == 0 {
                                    eprintln!(
                                        "    Found! union {} is same as existing union {}",
                                        union_rev, existing_union
                                    );
                                }
                                found = Some(existing_union);
                                break;
                            }
                        }

                        if let Some(existing) = found {
                            existing
                        } else {
                            // Truly new union - use forward by default
                            union_fwd
                        }
                    }
                };

                if verbose && (i < 5 || seq.id == "seq2") {
                    eprintln!("  [BIDIRECTED] {} pos {} - fwd_union: {}, rev_union: {}, chosen: {}",
                        seq.id, i, union_fwd, union_rev, union_rep);
                    // Debug: show what's already in union_to_node
                    if seq.id == "seq2" && i == 0 {
                        eprintln!("  Current union_to_node mappings:");
                        let mut mappings: Vec<_> = union_to_node.iter().collect();
                        mappings.sort_by_key(|(k, _)| *k);
                        for (union, node) in mappings {
                            eprintln!("    union {} -> node {}", union, node);
                        }
                    }
                }

                // Get or create node ID and determine its stored base
                let (node_id, node_base) = match union_to_node.get(&union_rep) {
                    Some(&id) => {
                        // Node already exists, get its base
                        let node_base = graph.nodes.get(id)
                            .and_then(|n| n.as_ref())
                            .map(|n| n.sequence[0])
                            .unwrap_or(seq.data[i]);
                        (id, node_base)
                    },
                    None => {
                        let id = next_node_id;
                        next_node_id += 1;
                        union_to_node.insert(union_rep, id);

                        // CRITICAL: Also map all other unions in this component to the same node!
                        // This ensures that when we encounter the same component through a different
                        // orientation or position, we'll find the existing node.
                        // Check both orientations of the current position
                        if self.union_find.same(pos_fwd, union_rep)
                            && !union_to_node.contains_key(&union_fwd)
                        {
                            union_to_node.insert(union_fwd, id);
                        }
                        if self.union_find.same(pos_rev, union_rep)
                            && !union_to_node.contains_key(&union_rev)
                        {
                            union_to_node.insert(union_rev, id);
                        }

                        // Create node with the base from the union representative position
                        // Store the forward base - paths will handle orientation
                        let base =
                            if let Some(source_seq) = self.find_sequence_for_position(union_rep) {
                                let local_offset = offset(union_rep) - source_seq.offset;
                                source_seq.data[local_offset]
                            } else {
                                seq.data[i] // Fallback to current sequence's base
                            };

                        graph.add_node(id, vec![base]);
                        (id, base)
                    }
                };

                // Determine path orientation based on whether we need to flip to get the correct base
                let expected_base = seq.data[i].to_ascii_uppercase();
                let node_base_upper = node_base.to_ascii_uppercase();

                let need_reverse = match (node_base_upper, expected_base) {
                    (b'A', b'T') => true,
                    (b'T', b'A') => true,
                    (b'C', b'G') => true,
                    (b'G', b'C') => true,
                    _ => false, // Same base or ambiguous - use forward
                };

                // Add handle with correct orientation to path
                let handle = Handle::new(node_id, need_reverse);
                path_handles.push(handle);
            }

            // Build path with handles
            graph.build_path(
                seq.id.clone(),
                path_handles
                    .iter()
                    .map(|h| (h.node_id(), h.is_reverse()))
                    .collect(),
            );
        }

        // Build edges from paths
        let edges_to_add: Vec<(Handle, Handle)> = graph
            .paths
            .iter()
            .flat_map(|path| {
                (0..path.steps.len().saturating_sub(1))
                    .map(move |i| (path.steps[i], path.steps[i + 1]))
            })
            .collect();

        for (from, to) in edges_to_add {
            graph.add_edge(from, to);
        }

        if verbose {
            println!(
                "Built bidirected graph: {} nodes, {} edges, {} paths",
                graph.nodes.len(),
                graph.edges.len(),
                graph.paths.len()
            );

            // Debug: count unique union representatives
            let mut unique_unions = std::collections::HashSet::new();
            for seq in &self.sequences {
                for i in 0..seq.data.len() {
                    let global_pos = seq.offset + i;
                    let pos_fwd = make_pos(global_pos, false);
                    let union_fwd = self.union_find.find(pos_fwd);
                    unique_unions.insert(union_fwd);
                }
            }
            println!(
                "Debug: {} unique union representatives from {} total positions",
                unique_unions.len(),
                self.sequences.iter().map(|s| s.data.len()).sum::<usize>()
            );
        }
        
        if verbose {
            eprintln!("[bidirected_builder] Graph built with {} nodes", graph.node_count());
            let first_10: Vec<usize> = graph.nodes.iter()
                .enumerate()
                .filter_map(|(id, n)| if n.is_some() { Some(id) } else { None })
                .take(10)
                .collect();
            eprintln!("[bidirected_builder] First 10 node IDs: {:?}", first_10);
            eprintln!("[bidirected_builder] Does graph have node 1? {}",
                1 < graph.nodes.len() && graph.nodes[1].is_some());
            eprintln!("[bidirected_builder] Does graph have node 2? {}",
                2 < graph.nodes.len() && graph.nodes[2].is_some());
        }

        // Verify all path edges exist before sorting
        if verbose {
            eprintln!("[bidirected_builder] Verifying path edges...");
        }
        graph.verify_path_edges(verbose);
        
        // Apply exact ODGI topological sort if requested
        if apply_sort {
            if verbose {
                eprintln!("[bidirected_builder] Applying exact ODGI topological sort...");
            }
            graph.apply_exact_odgi_ordering(verbose);
        } else if verbose {
            eprintln!("[bidirected_builder] Skipping topological sort (--no-sort)");
        }
        
        // Verify edges again after sorting
        graph.verify_path_edges(verbose);
        
        Ok(graph)
    }

    /// Convert bidirected graph back to simple graph (temporary compatibility)
    pub fn bidirected_to_simple_graph(&self, bi_graph: BidirectedGraph) -> Graph {
        let mut graph = Graph::new();

        // Convert nodes
        for (id, node_opt) in bi_graph.nodes.iter().enumerate() {
            if let Some(bi_node) = node_opt {
                graph.nodes.insert(
                    id,
                    Node {
                        id,
                        sequence: bi_node.sequence.clone(),
                        rank: bi_node.rank.unwrap_or(id as u64) as f64,
                    },
                );
            }
        }

        // Convert edges (lose orientation information)
        for edge in bi_graph.edges {
            graph.edges.insert(Edge {
                from: edge.from.node_id(),
                to: edge.to.node_id(),
            });
        }

        // Convert paths (lose orientation information)
        for path in bi_graph.paths {
            let simple_path: Vec<usize> =
                path.steps.iter().map(|handle| handle.node_id()).collect();
            graph.paths.push((path.name, simple_path));
        }

        graph
    }

    /// Find which sequence contains a given position
    fn find_sequence_for_position(&self, pos: Pos) -> Option<&Sequence> {
        let offset_val = offset(pos);
        self.sequences
            .iter()
            .find(|seq| offset_val >= seq.offset && offset_val < seq.offset + seq.data.len())
    }
}
