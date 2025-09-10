/// Direct port of ODGI's topological sort algorithm
/// Based on odgi/src/algorithms/topological_sort.cpp

use crate::bidirected_graph::Handle;
use crate::bidirected_ops::BidirectedGraph;
use std::collections::{HashMap, HashSet};

impl BidirectedGraph {
    /// Find all nodes with no edges on their left sides (heads)
    pub fn head_nodes(&self) -> Vec<Handle> {
        let mut heads = Vec::new();
        
        for &node_id in self.nodes.keys() {
            let handle = Handle::forward(node_id);
            let mut has_left_edges = false;
            
            // Check if this handle has any incoming edges on its left side
            // In a bidirected graph, left side of forward handle receives edges to the forward orientation
            for edge in &self.edges {
                // Check if edge comes TO this handle's left side
                if edge.to == handle {
                    has_left_edges = true;
                    break;
                }
            }
            
            if !has_left_edges {
                heads.push(handle);
            }
        }
        
        // Sort for deterministic behavior
        heads.sort_by_key(|h| h.node_id());
        heads
    }
    
    /// Find all nodes with no edges on their right sides (tails)
    pub fn tail_nodes(&self) -> Vec<Handle> {
        let mut tails = Vec::new();
        
        for &node_id in self.nodes.keys() {
            let handle = Handle::forward(node_id);
            let mut has_right_edges = false;
            
            // Check if this handle has any outgoing edges on its right side
            for edge in &self.edges {
                // Check if edge goes FROM this handle's right side
                if edge.from == handle {
                    has_right_edges = true;
                    break;
                }
            }
            
            if !has_right_edges {
                tails.push(handle);
            }
        }
        
        // Sort for deterministic behavior
        tails.sort_by_key(|h| h.node_id());
        tails
    }

    /// ODGI's topological_order algorithm
    /// This is a bidirected adaptation of Kahn's topological sort that can handle cycles
    pub fn odgi_topological_order(
        &self,
        use_heads: bool,
        use_tails: bool,
        verbose: bool,
    ) -> Vec<Handle> {
        let mut sorted = Vec::new();
        
        if self.nodes.is_empty() {
            return sorted;
        }
        
        // Track which node IDs have been added to the result
        let mut emitted_nodes = HashSet::new();
        
        // S - set of oriented handles ready to be processed  
        let mut s = HashSet::new();
        
        // Unvisited - track which handles haven't been processed yet
        let mut unvisited = HashMap::new();
        for &node_id in self.nodes.keys() {
            unvisited.insert(Handle::forward(node_id), true);
            unvisited.insert(Handle::reverse(node_id), true);
        }
        
        // Seeds for breaking cycles - maps handle to priority
        let mut seeds: HashMap<Handle, i32> = HashMap::new();
        
        // Track masked (logically removed) edges
        let mut masked_edges = HashSet::new();
        
        // Initialize with heads or tails if requested
        if use_heads {
            for handle in self.head_nodes() {
                s.insert(handle);
                unvisited.remove(&handle);
                unvisited.remove(&handle.flip());
            }
            if verbose {
                eprintln!("[odgi_topo] Starting with {} head nodes", s.len());
            }
        } else if use_tails {
            for handle in self.tail_nodes() {
                s.insert(handle);
                unvisited.remove(&handle);
                unvisited.remove(&handle.flip());
            }
            if verbose {
                eprintln!("[odgi_topo] Starting with {} tail nodes", s.len());
            }
        }
        
        // Main loop - continue until all nodes are visited
        while !unvisited.is_empty() || !s.is_empty() {
            
            // If S is empty, need to pick a seed to break into a cycle
            if s.is_empty() {
                // First try previously identified seeds
                let mut found_seed = false;
                
                // Sort seeds for deterministic selection
                let mut seed_handles: Vec<_> = seeds.keys().cloned().collect();
                seed_handles.sort_by_key(|h| (h.node_id(), h.is_reverse()));
                
                for handle in seed_handles {
                    if unvisited.get(&handle).copied().unwrap_or(false) {
                        s.insert(handle);
                        unvisited.remove(&handle);
                        unvisited.remove(&handle.flip());
                        found_seed = true;
                        if verbose {
                            eprintln!("[odgi_topo] Using seed: node {} orient {}", 
                                     handle.node_id(), 
                                     if handle.is_reverse() { "rev" } else { "fwd" });
                        }
                        break;
                    }
                }
                
                // If no seeds available, pick arbitrary unvisited handle
                if !found_seed {
                    // Get minimum unvisited handle for deterministic behavior
                    let mut handles: Vec<_> = unvisited.keys()
                        .filter(|h| !h.is_reverse()) // Prefer forward orientation
                        .cloned()
                        .collect();
                    handles.sort_by_key(|h| h.node_id());
                    
                    if let Some(&handle) = handles.first() {
                        s.insert(handle);
                        unvisited.remove(&handle);
                        unvisited.remove(&handle.flip());
                        if verbose {
                            eprintln!("[odgi_topo] Using arbitrary: node {} orient fwd", 
                                     handle.node_id());
                        }
                    }
                }
            }
            
            // Process handles in S
            while !s.is_empty() {
                // Get minimum handle for deterministic behavior
                let handle = *s.iter().min().unwrap();
                s.remove(&handle);
                
                // Only emit each node once (in forward orientation)
                if !handle.is_reverse() && emitted_nodes.insert(handle.node_id()) {
                    sorted.push(handle);
                }
                
                // Look at edges coming into the left side of this handle (going backward)
                // These are edges we're "consuming" by placing this handle
                for edge in self.edges.clone() {
                    if edge.to == handle && !masked_edges.contains(&edge) {
                        masked_edges.insert(edge.clone());
                        // The source might be a cycle entry point, but we don't need to track it further
                    }
                }
                
                // Look at edges going out from the right side of this handle (going forward)
                // These lead to potential next handles to process
                for edge in self.edges.clone() {
                    if edge.from == handle && !masked_edges.contains(&edge) {
                        masked_edges.insert(edge.clone());
                        
                        let next_handle = edge.to;
                        
                        // Only process if not yet visited
                        if unvisited.get(&next_handle).copied().unwrap_or(false) {
                            // Check if next_handle has any other unmasked incoming edges
                            let mut has_unmasked_incoming = false;
                            
                            for other_edge in &self.edges {
                                if other_edge.to == next_handle && 
                                   !masked_edges.contains(other_edge) {
                                    has_unmasked_incoming = true;
                                    break;
                                }
                            }
                            
                            if !has_unmasked_incoming {
                                // No more incoming edges, ready to process
                                s.insert(next_handle);
                                unvisited.remove(&next_handle);
                                unvisited.remove(&next_handle.flip());
                            } else {
                                // Still has dependencies, mark as potential seed for cycle breaking
                                seeds.insert(next_handle, 1);
                            }
                        }
                    }
                }
            }
        }
        
        if verbose {
            eprintln!("[odgi_topo] Topological sort complete: {} nodes", sorted.len());
        }
        
        sorted
    }
    
    /// Apply the ODGI topological ordering to renumber nodes
    pub fn apply_odgi_ordering(&mut self, verbose: bool) {
        let ordering = self.odgi_topological_order(true, false, verbose);
        
        if ordering.is_empty() {
            return;
        }
        
        // Create old to new ID mapping (1-based numbering)
        let mut old_to_new = HashMap::new();
        for (new_idx, handle) in ordering.iter().enumerate() {
            old_to_new.insert(handle.node_id(), new_idx + 1);
        }
        
        if verbose {
            eprintln!("[odgi_topo] Renumbering {} nodes", old_to_new.len());
        }
        
        // Apply the renumbering
        self.apply_ordering(ordering, verbose);
    }
}