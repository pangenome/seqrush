use std::fmt;

/// A handle represents an oriented reference to a node in the graph.
/// The least significant bit (LSB) indicates orientation:
/// - 0 = forward strand
/// - 1 = reverse strand
/// The remaining bits store the node ID.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Handle(u64);

impl Handle {
    /// Create a new handle with the given node ID and orientation
    pub fn new(node_id: usize, is_reverse: bool) -> Self {
        let mut value = (node_id as u64) << 1;
        if is_reverse {
            value |= 1;
        }
        Handle(value)
    }

    /// Create a forward handle for the given node ID
    pub fn forward(node_id: usize) -> Self {
        Self::new(node_id, false)
    }

    /// Create a reverse handle for the given node ID
    pub fn reverse(node_id: usize) -> Self {
        Self::new(node_id, true)
    }

    /// Get the node ID from this handle
    pub fn node_id(&self) -> usize {
        (self.0 >> 1) as usize
    }

    /// Check if this handle is in reverse orientation
    pub fn is_reverse(&self) -> bool {
        (self.0 & 1) == 1
    }

    /// Get the orientation sign as a char ('+' or '-')
    pub fn orientation_char(&self) -> char {
        if self.is_reverse() {
            '-'
        } else {
            '+'
        }
    }

    /// Flip the orientation of this handle
    pub fn flip(&self) -> Self {
        Handle(self.0 ^ 1)
    }

    /// Get the raw u64 value
    pub fn as_u64(&self) -> u64 {
        self.0
    }

    /// Create from raw u64 value
    pub fn from_u64(value: u64) -> Self {
        Handle(value)
    }
}

impl fmt::Display for Handle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.node_id(), self.orientation_char())
    }
}

/// Compute the reverse complement of a DNA sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&base| match base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            _ => base, // Keep any other characters unchanged
        })
        .collect()
}

/// A bidirected graph node containing a DNA sequence
#[derive(Debug, Clone)]
pub struct BiNode {
    pub id: usize,
    pub sequence: Vec<u8>,
    pub rank: Option<u64>,
}

impl BiNode {
    /// Create a new node
    pub fn new(id: usize, sequence: Vec<u8>) -> Self {
        BiNode {
            id,
            sequence,
            rank: None,
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

    /// Get sequence as string in the specified orientation
    pub fn get_sequence_string(&self, is_reverse: bool) -> String {
        String::from_utf8_lossy(&self.get_sequence(is_reverse)).to_string()
    }
}

/// A path through the bidirected graph
#[derive(Debug, Clone)]
pub struct BiPath {
    pub name: String,
    pub steps: Vec<Handle>,
}

impl BiPath {
    /// Create a new path
    pub fn new(name: String) -> Self {
        BiPath {
            name,
            steps: Vec::new(),
        }
    }

    /// Add a step to the path
    pub fn add_step(&mut self, handle: Handle) {
        self.steps.push(handle);
    }

    /// Get the full sequence of this path given a node lookup
    pub fn get_sequence<'a, F>(&self, get_node: F) -> Vec<u8>
    where
        F: Fn(usize) -> Option<&'a BiNode>,
    {
        let mut result = Vec::new();
        for &handle in &self.steps {
            if let Some(node) = get_node(handle.node_id()) {
                result.extend(node.get_sequence(handle.is_reverse()));
            }
        }
        result
    }
}

/// An edge in the bidirected graph
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct BiEdge {
    pub from: Handle,
    pub to: Handle,
}

impl BiEdge {
    /// Create a new edge
    pub fn new(from: Handle, to: Handle) -> Self {
        BiEdge { from, to }
    }

    /// Get the canonical form of this edge (smaller handle first)
    pub fn canonical(&self) -> Self {
        if self.from <= self.to {
            *self
        } else {
            BiEdge {
                from: self.to.flip(),
                to: self.from.flip(),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_handle_creation() {
        let h1 = Handle::forward(42);
        assert_eq!(h1.node_id(), 42);
        assert!(!h1.is_reverse());
        assert_eq!(h1.orientation_char(), '+');

        let h2 = Handle::reverse(42);
        assert_eq!(h2.node_id(), 42);
        assert!(h2.is_reverse());
        assert_eq!(h2.orientation_char(), '-');
    }

    #[test]
    fn test_handle_flip() {
        let h1 = Handle::forward(10);
        let h2 = h1.flip();
        assert_eq!(h2.node_id(), 10);
        assert!(h2.is_reverse());

        let h3 = h2.flip();
        assert_eq!(h3, h1);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GCTA"), b"TAGC");
        assert_eq!(reverse_complement(b"N"), b"N");
    }

    #[test]
    fn test_binode_sequences() {
        let node = BiNode::new(1, b"ATCG".to_vec());
        assert_eq!(node.get_sequence(false), b"ATCG");
        assert_eq!(node.get_sequence(true), b"CGAT");
    }

    #[test]
    fn test_path_sequence() {
        let node1 = BiNode::new(1, b"ATG".to_vec());
        let node2 = BiNode::new(2, b"CGA".to_vec());

        let mut path = BiPath::new("test".to_string());
        path.add_step(Handle::forward(1));
        path.add_step(Handle::reverse(2));

        let nodes = [None, Some(&node1), Some(&node2)];
        let get_node = |id: usize| nodes.get(id).and_then(|n| *n);

        let seq = path.get_sequence(get_node);
        assert_eq!(seq, b"ATGTCG"); // ATG + reverse_complement(CGA)
    }
}
