# ODGI Unchop Algorithm Documentation

## Overview
The unchop algorithm merges linear chains of nodes in a bidirected graph while preserving all path information. It's designed to reduce graph complexity without changing the sequences represented by the paths.

## Key Components

### 1. Simple Components Discovery (`simple_components.cpp`)

The algorithm finds "simple components" - chains of nodes that can be safely merged:

1. **Handle Creation**: For each node, create both forward and reverse handles
2. **Perfect Hash Function**: Build a MPHF (Minimal Perfect Hash Function) for O(1) handle lookups
3. **Disjoint Set Union**: Initialize a union-find structure for grouping handles
4. **Chain Detection**: For each handle:
   - If in-degree = 1: Check backward neighbor
     - Neighbor must have out-degree = 1
     - Must be perfect path neighbors
     - Different node IDs (no self-loops)
   - If out-degree = 1: Check forward neighbor
     - Same constraints apply
   - Unite handles in the disjoint set if conditions met

### 2. Perfect Path Neighbors (`perfect_neighbors.cpp`)

Two handles are "perfect path neighbors" if:
1. **All paths through left handle continue to right handle**
   - No path ends at left handle
   - No path from left goes to a different handle
2. **Path step counts match**
   - Number of steps leaving left = number entering right
   - This ensures no path "appears" or "disappears"

Algorithm:
```
for each path step on left_handle:
    if path ends here: return false
    next_step = follow path forward/backward
    if next_step != right_handle: return false
    expected_next++

count steps on right_handle
return (observed_next == expected_next)
```

### 3. Component Ordering

After finding components via union-find:
1. **Group by Component**: Collect all handles with same union-find root
2. **Find Start**: For each component ≥ min_size:
   - Find handle with in-degree = 0 OR
   - Handle with in-degree = 1 but predecessor outside component
   - If circular, use first handle
3. **Walk Component**: Starting from identified handle:
   - Follow edges forward within component
   - Build ordered list of handles
   - Fail if can't reach all handles (indicates complex topology)

### 4. Node Concatenation (`concat_nodes` in `unchop.cpp`)

For each ordered component:
1. **Create New Node**:
   - Concatenate sequences from all handles in order
   - Create node with combined sequence
2. **Reconnect Edges**:
   - Find all edges entering first handle (from outside component)
   - Find all edges leaving last handle (to outside component)
   - Handle self-loops specially:
     - Loop from last back to first → self-loop on new node
     - Loop from first back to last → reverse self-loop
3. **Update Paths**:
   - For each path step through component handles:
     - Rewrite to single step on new node
     - Preserve path continuity
4. **Cleanup**: Delete original nodes

### 5. Validation

After unchop:
1. **Path Sequence Preservation**: 
   - Extract sequence from each path before/after
   - Verify sequences are identical
2. **Node Ordering**: 
   - Maintain approximate topological order
   - Use average rank of merged nodes

## Key Invariants

1. **Path Preservation**: Every path's sequence must remain unchanged
2. **No Information Loss**: Can reconstruct original graph from paths
3. **Degree Constraints**: Only merge nodes forming simple chains (degree 1-1)
4. **Path Coherence**: All paths must traverse merged nodes identically

## Differences from Simple Merging

1. **Bidirected**: Handles both orientations of each node
2. **Path-Aware**: Only merges if ALL paths agree on traversal
3. **Atomic**: Either merges entire component or none
4. **Self-Loop Handling**: Special cases for loops within components

## Common Pitfalls

1. **Orientation Confusion**: Must track handle orientations carefully
2. **Path Step Counting**: Missing reverse complement paths
3. **Edge Reconnection**: Forgetting to handle self-loops
4. **Component Ordering**: Circular components need special handling