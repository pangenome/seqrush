# Algorithm Comparison: Our Implementation vs ODGI

## Key Differences Found

### 1. **Handle Representation & Packing**
**ODGI:**
- Uses `number_bool_packing::pack/unpack` to encode node ID + orientation in a single integer
- The packed integer encodes both the node and its orientation efficiently
- Uses `unpack_number()` to get the rank, `unpack_bit()` for orientation

**Ours:**
- Uses `Handle` struct directly with separate `id()` and `is_reverse()` methods
- Less efficient packing/unpacking

### 2. **Edge Masking System**
**ODGI:**
- Uses succinct bitvectors for efficient edge masking
- Delta encoding for edge storage (encodes edge as offset from source)
- `is_masked_edge()` checks if edge was already processed

**Ours:**
- Simple HashSet for processed handles
- No edge masking - we just track processed handles
- **CRITICAL MISS:** We don't mask edges, so we might process them multiple times

### 3. **The Core Loop Structure**
**ODGI:**
```cpp
while(unvisited.rank1(unvisited.size())!=0 || s.rank1(s.size())!=0) {
    // First: Try to pull from seeds if s is empty
    while(s.rank1(s.size())==0 && seeds.rank1(seeds.size())!=0) {
        // Get seed and add to s
    }
    
    // Second: If still empty, grab arbitrary unvisited node
    if(s.rank1(s.size())==0) {
        // Take first unvisited node
    }
    
    // Third: Process all nodes in s
    while (s.rank1(s.size())!=0) {
        // Process node
        // Follow edges backward to mask cycle edges
        // Follow edges forward to find next nodes
    }
}
```

**Ours:**
```rust
while !queue.is_empty() || !seeds.is_empty() {
    // Only processes seeds when queue is empty
    // Doesn't handle the two-phase edge following
}
```

### 4. **Critical Algorithm Difference: Two-Phase Edge Following**
**ODGI does TWO passes:**
1. `follow_edges(n, true, ...)` - Looks BACKWARD to find and mask edges from cycle entry points
2. `follow_edges(n, false, ...)` - Looks FORWARD to process outgoing edges

**We only look forward!** This is probably why we miss proper cycle handling.

### 5. **Node Output**
**ODGI:**
- Outputs handles (node + orientation) to the sorted vector
- Only outputs forward orientation of each node once

**Ours:**
- We track nodes, not handles
- We mark both orientations as processed when we see one

### 6. **Seed Management**
**ODGI:**
- Seeds are handles with specific orientations
- Only adds seeds if that node hasn't been seeded before
- Uses the exact orientation that was encountered

**Ours:**
- Seeds are handles but we don't track orientation properly
- We process both orientations together

### 7. **The "Last Incoming Edge" Check**
**ODGI's critical logic:**
```cpp
// After masking current edge, check if next_node has any other unmasked incoming edges
bool unmasked_incoming_edge = false;
g->follow_edges(next_node, true, [&](const handle_t& prev_node) {
    auto prev_edge = g->edge_handle(prev_node, next_node);
    if (!is_masked_edge(prev_edge)) {
        unmasked_incoming_edge = true;
        return false;
    }
    return true;
});

if(!unmasked_incoming_edge) {
    // This was the last incoming edge - node is ready to process
    s.set(next_node_rank, 1);
    unvisited.set(next_node_rank, 0);
} else {
    // Still has incoming edges - make it a seed for later
    seeds[next_node_rank] = 1;
    seeds_rev[next_node_rank] = orientation;
}
```

**Our version doesn't do this!** We use simple in-degree counting, but don't track masked edges properly.

## The Core Problem

Our implementation treats the graph as simple directed when it's actually bidirected. The key issues:

1. **We don't mask edges** - ODGI masks each edge as it's traversed to ensure it's only processed once
2. **We don't look backward** - ODGI checks backward edges to handle cycle entry points
3. **We process nodes, not handles** - ODGI processes specific orientations
4. **We don't check for "last incoming edge"** - ODGI only adds a node when ALL its incoming edges are masked

## What Needs to Change

To achieve parity with ODGI:

1. Implement proper edge masking (track which edges have been processed)
2. Add backward edge checking phase
3. Process handles (node+orientation) not just nodes
4. Implement the "last incoming edge" check using masked edges, not simple in-degree
5. Use proper orientation tracking for seeds
6. Output handles to maintain orientation information

The fundamental issue is that we're doing a simplified Kahn's algorithm when ODGI is doing a sophisticated bidirected topological sort with virtual edge removal.