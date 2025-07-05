# DESIGN.md - SeqRush Implementation Design

## Overview

SeqRush will be a high-performance Rust implementation of seqwish's transitive closure algorithm. This design document outlines the architecture and implementation strategy.

## Core Components

### 1. Lock-Free Union-Find

**Requirements:**
- Support for up to 2^63 positions (with orientation bit)
- Lock-free concurrent operations
- Path compression and union by rank
- Performance comparable to seqwish's implementation

**Implementation Options:**
1. **Use uf_rush** (current) - Verify it supports lock-free operations
2. **Custom implementation** - Port seqwish's dset64-gccAtomic.hpp to Rust
3. **Alternative crate** - Evaluate other union-find crates

**Design Decision:** 
- First, benchmark uf_rush against a simple test case
- If performance is inadequate, implement custom lock-free union-find using Rust's atomic types

### 2. Position Representation

```rust
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
struct Position {
    value: u64,  // LSB = strand, remaining 63 bits = offset
}

impl Position {
    fn new(offset: usize, reverse: bool) -> Self {
        Position {
            value: (offset as u64) << 1 | (reverse as u64)
        }
    }
    
    fn offset(&self) -> usize {
        (self.value >> 1) as usize
    }
    
    fn is_reverse(&self) -> bool {
        (self.value & 1) != 0
    }
}
```

### 3. Batch Processing Architecture

```rust
struct BatchProcessor {
    batch_size: usize,  // Default 1MB worth of sequence
    sequences: Vec<Sequence>,
    union_find: Arc<UnionFind>,
    processed: Arc<Vec<AtomicBool>>,
}

impl BatchProcessor {
    fn process_batch(&self, start: usize, end: usize) {
        // Process a chunk of positions
        // Apply transitive closure within batch
        // Update global union-find structure
    }
}
```

### 4. Match Processing Pipeline

```rust
struct MatchProcessor {
    min_match_length: usize,
    sparsification_factor: Option<f64>,
}

impl MatchProcessor {
    fn process_alignment(&self, alignment: &Alignment) -> Vec<Match> {
        // Extract matches from CIGAR
        // Filter by minimum length
        // Apply sparsification if enabled
        // Return validated matches
    }
}
```

### 5. Transitive Closure Algorithm

```rust
impl SeqRush {
    fn build_transitive_closure(&mut self) {
        // Phase 1: Process all matches
        for match in &self.matches {
            self.process_match(match);
        }
        
        // Phase 2: Process positions in batches
        let batch_size = self.calculate_batch_size();
        let batches: Vec<_> = (0..self.total_length)
            .step_by(batch_size)
            .collect();
            
        batches.par_iter().for_each(|&start| {
            let end = (start + batch_size).min(self.total_length);
            self.process_batch(start, end);
        });
    }
    
    fn process_batch(&self, start: usize, end: usize) {
        // Create work queue for this batch
        let queue = Arc::new(Mutex::new(VecDeque::new()));
        
        // Initialize with unprocessed positions
        for pos in start..end {
            if !self.processed[pos].load(Ordering::Acquire) {
                queue.lock().unwrap().push_back(pos);
            }
        }
        
        // Process queue in parallel
        (0..num_threads).into_par_iter().for_each(|_| {
            loop {
                let pos = {
                    let mut q = queue.lock().unwrap();
                    q.pop_front()
                };
                
                match pos {
                    Some(p) => self.process_position(p, &queue),
                    None => break,
                }
            }
        });
    }
}
```

### 6. Node Creation and Deduplication

```rust
struct NodeManager {
    nodes: Arc<Mutex<Vec<Node>>>,
    position_to_node: Arc<DashMap<usize, usize>>,
    node_to_positions: Arc<DashMap<usize, Vec<usize>>>,
}

impl NodeManager {
    fn create_node_for_component(&self, positions: Vec<usize>) -> usize {
        // Check if any position already has a node
        for &pos in &positions {
            if let Some(node_id) = self.position_to_node.get(&pos) {
                // Extend existing node with new positions
                self.extend_node(*node_id, positions);
                return *node_id;
            }
        }
        
        // Create new node
        let node_id = self.create_new_node(positions);
        node_id
    }
}
```

### 7. Path Construction

```rust
impl PathBuilder {
    fn build_path(&self, sequence: &Sequence) -> Vec<usize> {
        let mut path = Vec::new();
        let mut last_node = None;
        
        for pos in 0..sequence.len() {
            let global_pos = sequence.offset + pos;
            let node_id = self.position_to_node[&global_pos];
            
            // Avoid duplicate consecutive nodes
            if last_node != Some(node_id) {
                path.push(node_id);
                last_node = Some(node_id);
            }
        }
        
        path
    }
}
```

## Performance Optimizations

### 1. Memory Layout
- Use cache-aligned data structures
- Group frequently accessed data together
- Minimize pointer chasing

### 2. Parallelization Strategy
- Batch processing with work stealing
- Thread-local caches for hot data
- Atomic operations only where necessary

### 3. I/O Optimization
- Buffered writing for GFA output
- Memory-mapped files for large inputs
- Streaming processing where possible

## Error Handling

### 1. Duplicate Nodes
- Track which positions have been processed
- Ensure each position maps to exactly one node
- Validate paths during construction

### 2. Memory Management
- Use Arena allocators for temporary data
- Clear intermediate structures after each batch
- Monitor memory usage and adjust batch size

## Testing Strategy

### 1. Unit Tests
- Union-find operations
- Position encoding/decoding
- Match extraction from CIGAR

### 2. Integration Tests
- Small test cases with known outputs
- Comparison with seqwish results
- Performance benchmarks

### 3. Validation
- Use odgi validate on output
- Check for duplicate nodes in paths
- Verify sequence reconstruction

## Key Differences from Current Implementation

### 1. Algorithm Flow
**Current SeqRush:**
- Collects all matches upfront
- Processes all positions in one pass
- Creates nodes on-the-fly

**Seqwish:**
- Processes in batches to manage memory
- Uses atomic queues for work distribution
- Sorts and compresses disjoint sets before node creation

### 2. Data Structures
**Current SeqRush:**
- Simple Match struct
- Direct position-to-node mapping
- No interval trees

**Seqwish:**
- Interval trees for efficient range queries
- Bit vectors for node boundaries
- Separate graph sequence from input sequences

### 3. Performance Issues
- No batch processing (processes entire sequence at once)
- Sequential match processing
- Inefficient transitive closure computation
- No memory management for large inputs

## Recommended Implementation

### Phase 1: Core Infrastructure
- [ ] Implement Position type with orientation bit
- [ ] Add interval tree data structure (or use existing crate)
- [ ] Implement atomic queue for work distribution
- [ ] Create proper match storage with bidirectional indexing

### Phase 2: Batch Processing Algorithm
- [ ] Implement chunk-based processing (default 1MB)
- [ ] Add parallel overlap exploration
- [ ] Implement proper union-find with compression
- [ ] Add node boundary tracking

### Phase 3: Graph Construction
- [ ] Build graph sequence incrementally
- [ ] Track node boundaries with bit vector
- [ ] Implement path extraction with deduplication
- [ ] Add proper edge derivation

### Phase 4: Optimizations
- [ ] Memory-mapped files for large inputs
- [ ] Parallel sorting with IPS⁴o
- [ ] Lock-free data structures throughout
- [ ] Cache-friendly memory layout

## Critical Fixes Needed

1. **Duplicate Nodes**: Track node boundaries properly
2. **Performance**: Implement batch processing
3. **Memory**: Don't load everything into memory at once
4. **Correctness**: Follow exact transitive closure algorithm