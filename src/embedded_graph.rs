use std::collections::{HashMap, HashSet};
use crate::bidirected_graph::{Handle, reverse_complement};

/// A path identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PathId(pub usize);

/// A step in a path through the graph
#[derive(Debug, Clone)]
pub struct PathStep {
    pub path_id: PathId,
    pub handle: Handle,  // The handle used to visit this node
    pub prev: Option<Handle>,
    pub next: Option<Handle>,
}

/// A node in the embedded graph where paths are stored as linked lists
#[derive(Debug, Clone)]
pub struct EmbeddedNode {
    pub id: usize,
    pub sequence: Vec<u8>,
    /// For each (path, orientation) that visits this node, store the path step
    /// Key is (PathId, is_reverse)
    pub path_steps: HashMap<(PathId, bool), PathStep>,
}

impl EmbeddedNode {
    pub fn new(id: usize, sequence: Vec<u8>) -> Self {
        EmbeddedNode {
            id,
            sequence,
            path_steps: HashMap::new(),
        }
    }
    
    /// Get the sequence in the specified orientation
    pub fn get_sequence(&self, is_reverse: bool) -> Vec<u8> {
        if is_reverse {
            reverse_complement(&self.sequence)
        } else {
            self.sequence.clone()
        }
    }
}

/// Path metadata (name and terminal nodes)
#[derive(Debug, Clone)]
pub struct PathMeta {
    pub name: String,
    pub start: Option<Handle>,
    pub end: Option<Handle>,
}

/// An embedded bidirected graph where paths are stored as linked lists through nodes
#[derive(Debug)]
pub struct EmbeddedGraph {
    /// The nodes in the graph
    pub nodes: HashMap<usize, EmbeddedNode>,
    /// Path metadata
    pub paths: HashMap<PathId, PathMeta>,
    /// Next available node ID
    next_node_id: usize,
    /// Next available path ID
    next_path_id: usize,
}

impl EmbeddedGraph {
    pub fn new() -> Self {
        EmbeddedGraph {
            nodes: HashMap::new(),
            paths: HashMap::new(),
            next_node_id: 1,
            next_path_id: 0,
        }
    }
    
    /// Add a node to the graph
    pub fn add_node(&mut self, sequence: Vec<u8>) -> usize {
        let id = self.next_node_id;
        self.next_node_id += 1;
        self.nodes.insert(id, EmbeddedNode::new(id, sequence));
        id
    }
    
    /// Create a new path
    pub fn add_path(&mut self, name: String) -> PathId {
        let id = PathId(self.next_path_id);
        self.next_path_id += 1;
        self.paths.insert(id, PathMeta {
            name,
            start: None,
            end: None,
        });
        id
    }
    
    /// Add a step to a path
    pub fn extend_path(&mut self, path_id: PathId, handle: Handle) -> Result<(), String> {
        let node_id = handle.node_id();
        
        // Get the path metadata
        let path_meta = self.paths.get_mut(&path_id)
            .ok_or_else(|| format!("Path {:?} not found", path_id))?;
        
        // Get the node
        let node = self.nodes.get_mut(&node_id)
            .ok_or_else(|| format!("Node {} not found", node_id))?;
        
        // Determine the previous handle
        let prev_handle = path_meta.end;
        
        // Create the path step
        let step = PathStep {
            path_id,
            handle,
            prev: prev_handle,
            next: None,
        };
        
        // Insert the step into the node, keyed by (path_id, orientation)
        node.path_steps.insert((path_id, handle.is_reverse()), step);
        
        // Update the previous node's next pointer if it exists
        if let Some(prev) = prev_handle {
            let prev_node = self.nodes.get_mut(&prev.node_id())
                .ok_or_else(|| format!("Previous node {} not found", prev.node_id()))?;
            
            // The previous step was stored with its orientation
            if let Some(prev_step) = prev_node.path_steps.get_mut(&(path_id, prev.is_reverse())) {
                prev_step.next = Some(handle);
            }
        }
        
        // Update path metadata
        if path_meta.start.is_none() {
            path_meta.start = Some(handle);
        }
        path_meta.end = Some(handle);
        
        Ok(())
    }
    
    /// Get all handles that come after this handle in any path
    pub fn get_next_steps(&self, handle: Handle) -> HashSet<Handle> {
        let mut next_handles = HashSet::new();
        
        if let Some(node) = self.nodes.get(&handle.node_id()) {
            // Look for steps that use this specific handle
            for ((_, is_reverse), step) in &node.path_steps {
                if *is_reverse == handle.is_reverse() {
                    // This step uses the same orientation we're querying
                    if let Some(next) = step.next {
                        next_handles.insert(next);
                    }
                }
            }
        }
        
        next_handles
    }
    
    /// Get all handles that come before this handle in any path
    pub fn get_prev_steps(&self, handle: Handle) -> HashSet<Handle> {
        let mut prev_handles = HashSet::new();
        
        if let Some(node) = self.nodes.get(&handle.node_id()) {
            // Look for steps that use this specific handle
            for ((_, is_reverse), step) in &node.path_steps {
                if *is_reverse == handle.is_reverse() {
                    // This step uses the same orientation we're querying
                    if let Some(prev) = step.prev {
                        prev_handles.insert(prev);
                    }
                }
            }
        }
        
        prev_handles
    }
    
    /// Check if two handles are perfect neighbors
    pub fn are_perfect_neighbors(&self, from: Handle, to: Handle) -> bool {
        let next_steps = self.get_next_steps(from);
        let prev_steps = self.get_prev_steps(to);
        
        // Perfect neighbors if:
        // 1. from only goes to `to`
        // 2. to only comes from `from`
        next_steps.len() == 1 && next_steps.contains(&to) &&
        prev_steps.len() == 1 && prev_steps.contains(&from)
    }
    
    /// Find all perfect neighbor pairs in the graph
    pub fn find_perfect_pairs(&self) -> Vec<(Handle, Handle)> {
        let mut pairs = Vec::new();
        
        for &node_id in self.nodes.keys() {
            for is_reverse in [false, true] {
                let handle = Handle::new(node_id, is_reverse);
                let next_steps = self.get_next_steps(handle);
                
                if next_steps.len() == 1 {
                    let next_handle = *next_steps.iter().next().unwrap();
                    if self.are_perfect_neighbors(handle, next_handle) {
                        pairs.push((handle, next_handle));
                    }
                }
            }
        }
        
        pairs
    }
    
    /// Merge two nodes that are perfect neighbors
    pub fn merge_perfect_neighbors(&mut self, from: Handle, to: Handle) -> Result<usize, String> {
        if !self.are_perfect_neighbors(from, to) {
            return Err(format!("{} and {} are not perfect neighbors", from, to));
        }
        
        // Get both nodes
        let from_node = self.nodes.get(&from.node_id())
            .ok_or_else(|| format!("Node {} not found", from.node_id()))?
            .clone();
        let to_node = self.nodes.get(&to.node_id())
            .ok_or_else(|| format!("Node {} not found", to.node_id()))?
            .clone();
        
        // Create merged sequence
        let mut merged_seq = Vec::new();
        merged_seq.extend_from_slice(&from_node.get_sequence(from.is_reverse()));
        merged_seq.extend_from_slice(&to_node.get_sequence(to.is_reverse()));
        
        // Create new node
        let new_id = self.add_node(merged_seq);
        let new_handle = Handle::forward(new_id);
        
        // For each path that goes through from->to, update it to go through the new node
        for ((path_id, from_is_rev), from_step) in &from_node.path_steps {
            if from_step.next == Some(to) && *from_is_rev == from.is_reverse() {
                // This path goes from->to, update it
                let new_node = self.nodes.get_mut(&new_id).unwrap();
                
                // Get the step from 'to' node
                let to_key = (*path_id, to.is_reverse());
                let to_step = to_node.path_steps.get(&to_key)
                    .ok_or_else(|| format!("Expected to find path {:?} in to_node", path_id))?;
                
                // Create new step in the merged node
                let new_step = PathStep {
                    path_id: *path_id,
                    handle: new_handle,
                    prev: from_step.prev,
                    next: to_step.next,
                };
                
                new_node.path_steps.insert((*path_id, false), new_step);
                
                // Update previous node's next pointer
                if let Some(prev_handle) = from_step.prev {
                    if let Some(prev_node) = self.nodes.get_mut(&prev_handle.node_id()) {
                        let prev_key = (*path_id, prev_handle.is_reverse());
                        if let Some(prev_step) = prev_node.path_steps.get_mut(&prev_key) {
                            prev_step.next = Some(new_handle);
                        }
                    }
                }
                
                // Update next node's prev pointer
                if let Some(next_handle) = to_step.next {
                    if let Some(next_node) = self.nodes.get_mut(&next_handle.node_id()) {
                        let next_key = (*path_id, next_handle.is_reverse());
                        if let Some(next_step) = next_node.path_steps.get_mut(&next_key) {
                            next_step.prev = Some(new_handle);
                        }
                    }
                }
                
                // Update path metadata if needed
                if let Some(path_meta) = self.paths.get_mut(path_id) {
                    if path_meta.start == Some(from) {
                        path_meta.start = Some(new_handle);
                    }
                    if path_meta.end == Some(to) {
                        path_meta.end = Some(new_handle);
                    }
                }
            }
        }
        
        // Remove the old nodes
        self.nodes.remove(&from.node_id());
        self.nodes.remove(&to.node_id());
        
        Ok(new_id)
    }
    
    /// Perform iterative compaction
    pub fn compact(&mut self) -> Result<usize, String> {
        let mut total_merged = 0;
        let mut iteration = 0;
        
        loop {
            iteration += 1;
            eprintln!("Compaction iteration {}", iteration);
            
            let pairs = self.find_perfect_pairs();
            if pairs.is_empty() {
                eprintln!("No more perfect pairs found");
                break;
            }
            
            eprintln!("Found {} perfect pairs", pairs.len());
            
            let mut merged_in_iteration = 0;
            let mut processed = HashSet::new();
            
            for (from, to) in pairs {
                // Skip if we already processed either node
                if processed.contains(&from.node_id()) || processed.contains(&to.node_id()) {
                    continue;
                }
                
                // Try to merge
                match self.merge_perfect_neighbors(from, to) {
                    Ok(_new_id) => {
                        processed.insert(from.node_id());
                        processed.insert(to.node_id());
                        merged_in_iteration += 1;
                    }
                    Err(e) => {
                        eprintln!("Failed to merge {} -> {}: {}", from, to, e);
                    }
                }
            }
            
            eprintln!("Merged {} pairs in iteration {}", merged_in_iteration, iteration);
            total_merged += merged_in_iteration;
            
            if merged_in_iteration == 0 {
                break;
            }
        }
        
        Ok(total_merged)
    }
    
    /// Get the sequence for a path
    pub fn get_path_sequence(&self, path_id: PathId) -> Result<Vec<u8>, String> {
        let path_meta = self.paths.get(&path_id)
            .ok_or_else(|| format!("Path {:?} not found", path_id))?;
        
        let mut sequence = Vec::new();
        let mut current = path_meta.start;
        
        while let Some(handle) = current {
            let node = self.nodes.get(&handle.node_id())
                .ok_or_else(|| format!("Node {} not found", handle.node_id()))?;
            
            sequence.extend_from_slice(&node.get_sequence(handle.is_reverse()));
            
            // Get next handle - need to find the step for this path
            let mut found_next = None;
            for ((pid, is_rev), step) in &node.path_steps {
                if *pid == path_id && *is_rev == handle.is_reverse() {
                    found_next = step.next;
                    break;
                }
            }
            current = found_next;
        }
        
        Ok(sequence)
    }
}