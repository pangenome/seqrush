use crate::bidirected_graph::{reverse_complement, Handle};
use std::collections::{HashMap, HashSet};
use std::io::Write;

/// A path identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PathId(pub usize);

/// A unique step identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StepId {
    pub path_id: PathId,
    pub node_id: usize,
    pub occurrence: usize, // Which occurrence of this path on this node (0-based)
}

/// A step in a path through the graph
#[derive(Debug, Clone)]
pub struct PathStep {
    pub path_id: PathId,
    pub handle: Handle,                // The handle used to visit this node
    pub prev: Option<(Handle, usize)>, // (prev_handle, occurrence_on_that_node)
    pub next: Option<(Handle, usize)>, // (next_handle, occurrence_on_that_node)
}

/// A node in the embedded graph where paths are stored as linked lists
#[derive(Debug, Clone)]
pub struct EmbeddedNode {
    pub id: usize,
    pub sequence: Vec<u8>,
    /// For each path that visits this node, store ordered list of steps
    /// A path can visit a node multiple times
    pub path_steps: HashMap<(PathId, bool), Vec<PathStep>>,
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
        self.paths.insert(
            id,
            PathMeta {
                name,
                start: None,
                end: None,
            },
        );
        id
    }

    /// Add a step to a path
    pub fn extend_path(&mut self, path_id: PathId, handle: Handle) -> Result<(), String> {
        let node_id = handle.node_id();

        // Get the path metadata
        let path_meta = self
            .paths
            .get_mut(&path_id)
            .ok_or_else(|| format!("Path {:?} not found", path_id))?;

        // Determine the previous handle and its occurrence
        let prev_info = if let Some(prev_handle) = path_meta.end {
            // Find the last occurrence of this path on the previous node
            let prev_node = self
                .nodes
                .get(&prev_handle.node_id())
                .ok_or_else(|| format!("Previous node {} not found", prev_handle.node_id()))?;

            if let Some(prev_steps) = prev_node
                .path_steps
                .get(&(path_id, prev_handle.is_reverse()))
            {
                // Find the step with matching handle
                let matching_steps: Vec<_> = prev_steps
                    .iter()
                    .enumerate()
                    .filter(|(_, s)| s.handle == prev_handle)
                    .collect();

                if let Some((idx, _)) = matching_steps.last() {
                    Some((prev_handle, *idx))
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        };

        // Create the path step
        let step = PathStep {
            path_id,
            handle,
            prev: prev_info,
            next: None,
        };

        // Get the node and determine current occurrence
        let node = self
            .nodes
            .get_mut(&node_id)
            .ok_or_else(|| format!("Node {} not found", node_id))?;

        let steps = node
            .path_steps
            .entry((path_id, handle.is_reverse()))
            .or_default();
        let current_occurrence = steps.iter().filter(|s| s.handle == handle).count();

        steps.push(step);

        // Update the previous node's next pointer if it exists
        if let Some((prev_handle, prev_occurrence)) = prev_info {
            let prev_node = self
                .nodes
                .get_mut(&prev_handle.node_id())
                .ok_or_else(|| format!("Previous node {} not found", prev_handle.node_id()))?;

            if let Some(prev_steps) = prev_node
                .path_steps
                .get_mut(&(path_id, prev_handle.is_reverse()))
            {
                if let Some(prev_step) = prev_steps
                    .iter_mut()
                    .filter(|s| s.handle == prev_handle)
                    .nth(prev_occurrence)
                {
                    prev_step.next = Some((handle, current_occurrence));
                }
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
            // Look at all steps for all paths on this node
            for (_path_id, steps) in &node.path_steps {
                for step in steps {
                    // Only consider steps that use the same orientation
                    if step.handle == handle {
                        if let Some((next_handle, _)) = step.next {
                            next_handles.insert(next_handle);
                        }
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
            // Look at all steps for all paths on this node
            for (_path_id, steps) in &node.path_steps {
                for step in steps {
                    // Only consider steps that use the same orientation
                    if step.handle == handle {
                        if let Some((prev_handle, _)) = step.prev {
                            prev_handles.insert(prev_handle);
                        }
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
        next_steps.len() == 1
            && next_steps.contains(&to)
            && prev_steps.len() == 1
            && prev_steps.contains(&from)
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
        let from_node = self
            .nodes
            .get(&from.node_id())
            .ok_or_else(|| format!("Node {} not found", from.node_id()))?
            .clone();
        let to_node = self
            .nodes
            .get(&to.node_id())
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
        for ((path_id, _path_is_rev), from_steps) in &from_node.path_steps {
            // Find steps that use the correct orientation and go to 'to'
            for (from_idx, from_step) in from_steps.iter().enumerate() {
                if from_step.handle == from {
                    if let Some((next_handle, next_occurrence)) = from_step.next {
                        if next_handle == to {
                            // This path goes from->to, update it
                            let new_node = self.nodes.get_mut(&new_id).unwrap();

                            // Find the matching step in 'to' node
                            if let Some(to_steps) =
                                to_node.path_steps.get(&(*path_id, to.is_reverse()))
                            {
                                if let Some(to_step) = to_steps
                                    .iter()
                                    .filter(|s| s.handle == to)
                                    .nth(next_occurrence)
                                {
                                    // Create new step in the merged node
                                    let new_step = PathStep {
                                        path_id: *path_id,
                                        handle: new_handle,
                                        prev: from_step.prev,
                                        next: to_step.next,
                                    };

                                    let new_steps = new_node
                                        .path_steps
                                        .entry((*path_id, new_handle.is_reverse()))
                                        .or_default();
                                    let new_occurrence = new_steps.len();
                                    new_steps.push(new_step);

                                    // Update previous node's next pointer
                                    if let Some((prev_handle, prev_occurrence)) = from_step.prev {
                                        if let Some(prev_node) =
                                            self.nodes.get_mut(&prev_handle.node_id())
                                        {
                                            if let Some(prev_steps) = prev_node
                                                .path_steps
                                                .get_mut(&(*path_id, prev_handle.is_reverse()))
                                            {
                                                if let Some(prev_step) = prev_steps
                                                    .iter_mut()
                                                    .filter(|s| s.handle == prev_handle)
                                                    .nth(prev_occurrence)
                                                {
                                                    prev_step.next =
                                                        Some((new_handle, new_occurrence));
                                                }
                                            }
                                        }
                                    }

                                    // Update next node's prev pointer
                                    if let Some((next_handle, next_occurrence_on_next)) =
                                        to_step.next
                                    {
                                        if let Some(next_node) =
                                            self.nodes.get_mut(&next_handle.node_id())
                                        {
                                            if let Some(next_steps) = next_node
                                                .path_steps
                                                .get_mut(&(*path_id, next_handle.is_reverse()))
                                            {
                                                if let Some(next_step) = next_steps
                                                    .iter_mut()
                                                    .filter(|s| s.handle == next_handle)
                                                    .nth(next_occurrence_on_next)
                                                {
                                                    next_step.prev =
                                                        Some((new_handle, new_occurrence));
                                                }
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
                        }
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

            eprintln!(
                "Merged {} pairs in iteration {}",
                merged_in_iteration, iteration
            );
            total_merged += merged_in_iteration;

            if merged_in_iteration == 0 {
                break;
            }
        }

        Ok(total_merged)
    }

    /// Write the graph to GFA format
    pub fn write_gfa(&self, writer: &mut impl Write) -> Result<(), std::io::Error> {
        // Write header
        writeln!(writer, "H\tVN:Z:1.0")?;

        // Write nodes (segments)
        for (node_id, node) in &self.nodes {
            writeln!(
                writer,
                "S\t{}\t{}",
                node_id,
                String::from_utf8_lossy(&node.sequence)
            )?;
        }

        // Write edges (links) by traversing paths
        let mut written_edges = HashSet::new();

        for (_path_key, steps_by_path) in self
            .nodes
            .values()
            .flat_map(|n| &n.path_steps)
            .flat_map(|(pid, steps)| steps.iter().map(move |s| (pid, s)))
        {
            if let Some((next_handle, _)) = steps_by_path.next {
                let from = steps_by_path.handle;
                let to = next_handle;

                // Create edge key to avoid duplicates
                let edge_key = (
                    from.node_id(),
                    from.is_reverse(),
                    to.node_id(),
                    to.is_reverse(),
                );

                if !written_edges.contains(&edge_key) {
                    written_edges.insert(edge_key);

                    let from_orient = if from.is_reverse() { '-' } else { '+' };
                    let to_orient = if to.is_reverse() { '-' } else { '+' };

                    writeln!(
                        writer,
                        "L\t{}\t{}\t{}\t{}\t0M",
                        from.node_id(),
                        from_orient,
                        to.node_id(),
                        to_orient
                    )?;
                }
            }
        }

        // Write paths
        for (path_id, path_meta) in &self.paths {
            if let Some(start_handle) = path_meta.start {
                let mut path_string = String::new();
                let mut current = Some((start_handle, 0usize));
                let mut visited = HashSet::new();

                while let Some((handle, occurrence)) = current {
                    // Prevent infinite loops
                    let visit_key = (handle, occurrence);
                    if visited.contains(&visit_key) {
                        break; // Circular path
                    }
                    visited.insert(visit_key);

                    // Add to path string
                    if !path_string.is_empty() {
                        path_string.push(',');
                    }
                    path_string.push_str(&format!(
                        "{}{}",
                        handle.node_id(),
                        if handle.is_reverse() { '-' } else { '+' }
                    ));

                    // Get next handle
                    let mut found_next = None;
                    if let Some(node) = self.nodes.get(&handle.node_id()) {
                        if let Some(steps) = node.path_steps.get(&(*path_id, handle.is_reverse())) {
                            if let Some(step) =
                                steps.iter().filter(|s| s.handle == handle).nth(occurrence)
                            {
                                found_next = step.next;
                            }
                        }
                    }
                    current = found_next;
                }

                writeln!(writer, "P\t{}\t{}\t*", path_meta.name, path_string)?;
            }
        }

        Ok(())
    }

    /// Get the sequence for a path
    pub fn get_path_sequence(&self, path_id: PathId) -> Result<Vec<u8>, String> {
        let path_meta = self
            .paths
            .get(&path_id)
            .ok_or_else(|| format!("Path {:?} not found", path_id))?;

        let mut sequence = Vec::new();

        if let Some(start_handle) = path_meta.start {
            let mut current = Some((start_handle, 0usize));
            let mut visited = HashSet::new();

            while let Some((handle, occurrence)) = current {
                // Prevent infinite loops
                let visit_key = (handle, occurrence);
                if visited.contains(&visit_key) {
                    return Err(format!(
                        "Cycle detected in path {:?} at {:?}",
                        path_id, visit_key
                    ));
                }
                visited.insert(visit_key);

                let node = self
                    .nodes
                    .get(&handle.node_id())
                    .ok_or_else(|| format!("Node {} not found", handle.node_id()))?;

                sequence.extend_from_slice(&node.get_sequence(handle.is_reverse()));

                // Get next handle - find the specific occurrence
                let mut found_next = None;
                if let Some(steps) = node.path_steps.get(&(path_id, handle.is_reverse())) {
                    if let Some(step) = steps.iter().filter(|s| s.handle == handle).nth(occurrence)
                    {
                        found_next = step.next;
                    }
                }
                current = found_next;
            }
        }

        Ok(sequence)
    }
}
