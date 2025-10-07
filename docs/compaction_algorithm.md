# SeqRush Compaction Algorithm

## Overview

The compaction algorithm merges adjacent nodes in a bidirected graph that can be safely combined without changing the sequences represented by any path through the graph.

## Algorithm

### 1. Build Connection Maps from Path Traversals

The path traversals ARE the links\! For each path in the graph:
- Traverse the path step by step
- Each consecutive pair of nodes (A, B) in a path represents a link/edge
- Record these connections:
  - A has an outbound connection to B
  - B has an inbound connection from A

### 2. Perfect Neighbor Detection

For each pair of adjacent nodes (A, B) found in the paths:

1. **Check connections**: For nodes A and B to be perfect neighbors:
   - A's outbound connections must equal {B} (A only connects to B in all paths)
   - B's inbound connections must equal {A} (B only receives from A in all paths)
   - If we find ANY difference in these lists, we bail out - they cannot be merged

2. **Unite compatible nodes**: If A and B pass the above checks, they can be safely merged

### 3. Key Properties

- The path traversals define the graph connectivity
- The algorithm operates from one strand perspective of the graph
- It's mathematically sound: if two nodes have identical connection patterns (except for their connection to each other), they can be merged
- The algorithm is conservative: any difference in connection patterns prevents merging

### 4. Implementation Details

```rust
// Build connections from path traversals
let mut inbound: HashMap<Handle, HashSet<Handle>> = HashMap::new();
let mut outbound: HashMap<Handle, HashSet<Handle>> = HashMap::new();

// The paths define the links\!
for path in &graph.paths {
    for i in 0..path.steps.len() - 1 {
        let from = path.steps[i];
        let to = path.steps[i + 1];
        
        outbound.entry(from).or_default().insert(to);
        inbound.entry(to).or_default().insert(from);
    }
}

// Two nodes are perfect neighbors if:
// 1. A only points to B (outbound[A] == {B})
// 2. B only receives from A (inbound[B] == {A})
```

### 5. Correctness

This algorithm ensures:
- No path sequences are changed
- No connectivity is lost  
- The graph remains valid
- All biological sequences are preserved

## References

This algorithm is based on the mathematical properties of directed graphs and the specific requirements of pangenome graphs where sequence preservation is critical. The key insight is that in variation graphs, the paths through the graph define the connectivity - there are no edges independent of paths.
EOF < /dev/null