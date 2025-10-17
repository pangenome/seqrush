use crate::bidirected_graph::Handle;
use crate::bidirected_ops::BidirectedGraph;
use crate::embedded_graph::EmbeddedGraph;
use std::collections::HashMap;

/// Convert a BidirectedGraph to an EmbeddedGraph
pub fn convert_to_embedded(bi_graph: &BidirectedGraph) -> Result<EmbeddedGraph, String> {
    let mut embedded = EmbeddedGraph::new();

    // First, add all nodes
    let mut node_mapping = HashMap::new();
    for (old_id, node_opt) in bi_graph.nodes.iter().enumerate() {
        if let Some(node) = node_opt {
            let new_id = embedded.add_node(node.sequence.clone());
            node_mapping.insert(old_id, new_id);
        }
    }

    // Then, add all paths
    let mut path_mapping = HashMap::new();
    for (i, path) in bi_graph.paths.iter().enumerate() {
        let path_id = embedded.add_path(path.name.clone());
        path_mapping.insert(i, path_id);

        // Add each step in the path
        for &handle in &path.steps {
            let old_node_id = handle.node_id();
            let new_node_id = node_mapping
                .get(&old_node_id)
                .ok_or_else(|| format!("Node {} not found in mapping", old_node_id))?;

            let new_handle = Handle::new(*new_node_id, handle.is_reverse());
            embedded.extend_path(path_id, new_handle)?;
        }
    }

    Ok(embedded)
}

/// Build an EmbeddedGraph directly from sequences and union-find results
pub struct EmbeddedGraphBuilder {
    sequences: Vec<Vec<u8>>,
    sequence_names: Vec<String>,
}

impl EmbeddedGraphBuilder {
    pub fn new(sequences: Vec<Vec<u8>>, sequence_names: Vec<String>) -> Self {
        EmbeddedGraphBuilder {
            sequences,
            sequence_names,
        }
    }

    /// Build embedded graph from union-find components
    pub fn build_from_components(
        &self,
        components: &HashMap<usize, Vec<(usize, usize, bool)>>,
    ) -> Result<EmbeddedGraph, String> {
        let mut graph = EmbeddedGraph::new();

        // Create nodes from union components
        let mut component_to_node = HashMap::new();
        for (component_id, positions) in components {
            // Get the base from the first position
            let (seq_idx, pos, _) = positions[0];
            let base = self.sequences[seq_idx][pos];

            // Verify all positions have the same base
            for &(seq_idx, pos, is_reverse) in positions {
                let actual_base = if is_reverse {
                    match self.sequences[seq_idx][pos] {
                        b'A' => b'T',
                        b'T' => b'A',
                        b'C' => b'G',
                        b'G' => b'C',
                        b'N' => b'N',
                        _ => self.sequences[seq_idx][pos],
                    }
                } else {
                    self.sequences[seq_idx][pos]
                };

                if actual_base != base {
                    return Err(format!(
                        "Component {} has inconsistent bases: {} vs {}",
                        component_id, base as char, actual_base as char
                    ));
                }
            }

            let node_id = graph.add_node(vec![base]);
            component_to_node.insert(*component_id, node_id);
        }

        // Build paths through the graph
        for (seq_idx, sequence) in self.sequences.iter().enumerate() {
            let path_id = graph.add_path(self.sequence_names[seq_idx].clone());

            // For each position in the sequence, find its component and add to path
            for pos in 0..sequence.len() {
                // Find which component this position belongs to
                let mut found_component = None;
                for (component_id, positions) in components {
                    if positions.iter().any(|&(s, p, _)| s == seq_idx && p == pos) {
                        found_component = Some(*component_id);
                        break;
                    }
                }

                if let Some(component_id) = found_component {
                    let node_id = component_to_node[&component_id];

                    // Determine orientation based on whether this position was united
                    // with reverse orientation
                    let is_reverse = components[&component_id]
                        .iter()
                        .find(|&&(s, p, _)| s == seq_idx && p == pos)
                        .map(|&(_, _, is_rev)| is_rev)
                        .unwrap_or(false);

                    let handle = Handle::new(node_id, is_reverse);
                    graph.extend_path(path_id, handle)?;
                } else {
                    return Err(format!(
                        "Position {} in sequence {} not found in any component",
                        pos, seq_idx
                    ));
                }
            }
        }

        Ok(graph)
    }
}
