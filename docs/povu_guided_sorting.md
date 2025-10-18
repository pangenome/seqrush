# Povu-Guided Topological Sorting

## Problem Statement

Path-Guided SGD (the "Y" in Ygs sorting) has a catastrophic failure mode on certain graph structures where nodes that are adjacent in path traversals get placed far apart in the final ordering.

### Example: A-3105 Catastrophic Edges

**Graph:** HLA-zoo A-3105 (1,549 nodes)

**Catastrophic Edge:** Node 1367 → Node 687
- **Jump distance:** -680 positions (43.9% of graph!)
- **Path reality:** These nodes are **adjacent** in path `gi|528476637` (steps 575→576)
- **Graph ordering:** Node 1367 at position 1366 (88%), Node 687 at position 686 (44%)

**Root Cause:**
- SGD initializes nodes by ID order (not path order)
- Complex bubble structure creates competing constraints from multiple paths
- SGD gets trapped in local minimum even with 100 iterations
- **ODGI has the exact same problem** (A-3105: max backward jump -774)

### Why This Matters

- **Path coherence destroyed:** 33-86% of path steps go backward in graph
- **RMSE degraded:** 21,959 bp vs ~750 bp for well-ordered graphs (29x worse!)
- **Visualization broken:** Long backward edges make graph layouts unreadable
- **This is an inherent limitation of pure SGD**, not a bug

## Solution: Povu-Guided Hybrid Sorting

Use **povu** (https://github.com/pangenome/povu) to decompose the graph into bubble chains, then combine:

1. **Topological sorting within bubbles** - Perfect for local structure (deterministic, no catastrophic edges)
2. **Path-Guided SGD between bubble chains** - Good for global layout (handles complex path relationships)

### The Povu Decomposition

Povu identifies "flubbles" (bubbles) and organizes them into linear chains:

```bash
~/povu/bin/povu decompose -i graph.gfa -o regions/
```

**Output:** `.pvst` file describing bubble chains

**Example from A-3105:**
- **513 bubbles** identified in linear chains
- **Bubble 171:** Has 12 alternative paths (complex structural variant)
  - Contains nodes 687, 689 (start) and 1367-1548 (end)
  - These are the exact nodes causing catastrophic -680 position jumps!

### Algorithm Design

#### Current Ygs Pipeline:
```
Y (SGD) → g (groom) → s (topological sort heads only)
```

#### Proposed Povu-Guided Pipeline:
```
1. Decompose: Run povu to identify bubble chains
2. Local topo sort: Sort nodes within each bubble topologically
3. Global SGD: Use Path-Guided SGD to order bubble chains
4. Flatten: Convert bubble chain order to final node order
5. Groom: Apply grooming (local adjustments)
6. Topo heads: Final topological sort of heads (DAG cleanup)
```

### Detailed Algorithm

```rust
// 1. DECOMPOSE: Parse povu output to get bubble structure
struct Bubble {
    id: usize,
    nodes: Vec<usize>,           // Node IDs in this bubble
    alternatives: Vec<Vec<usize>>, // Alternative paths through bubble
}

struct BubbleChain {
    bubbles: Vec<Bubble>,
    start_node: usize,
    end_node: usize,
}

fn parse_povu_decomposition(pvst_file: &str) -> Vec<BubbleChain> {
    // Parse .pvst format:
    // H	0.0.3	.	.	.
    // D	0	.	[node list]	.
    // F	bubble_id	>from>to	.	L  (or alternatives)
    // ...
}

// 2. LOCAL TOPO SORT: Order nodes within each bubble
fn topologically_sort_bubble(bubble: &Bubble, graph: &BidirectedGraph) -> Vec<usize> {
    // Within bubble, nodes form a DAG (by definition of bubble)
    // Use standard topological sort (Kahn's algorithm or DFS)
    // This ensures no backward edges within bubble structure

    let mut in_degree = HashMap::new();
    let mut adj_list = HashMap::new();

    // Build subgraph for this bubble
    for &node in &bubble.nodes {
        for edge in graph.edges_from(node) {
            if bubble.nodes.contains(&edge.to) {
                *in_degree.entry(edge.to).or_insert(0) += 1;
                adj_list.entry(edge.from).or_default().push(edge.to);
            }
        }
    }

    // Kahn's algorithm
    let mut queue: VecDeque<_> = bubble.nodes.iter()
        .filter(|&&n| in_degree.get(&n).unwrap_or(&0) == &0)
        .copied()
        .collect();

    let mut sorted = Vec::new();
    while let Some(node) = queue.pop_front() {
        sorted.push(node);
        for &neighbor in adj_list.get(&node).unwrap_or(&vec![]) {
            let deg = in_degree.get_mut(&neighbor).unwrap();
            *deg -= 1;
            if *deg == 0 {
                queue.push_back(neighbor);
            }
        }
    }

    sorted
}

// 3. GLOBAL SGD: Order bubble chains using Path-Guided SGD
fn order_bubble_chains(chains: &[BubbleChain], graph: &BidirectedGraph) -> Vec<usize> {
    // Create "super-graph" where each bubble is a super-node
    // Run Path-Guided SGD on super-nodes
    // This handles global path relationships without bubble-internal conflicts

    let super_graph = create_bubble_chain_graph(chains, graph);

    // Use existing path_linear_sgd but on bubble chains
    let mut sgd_params = PathSGDParams::default();
    sgd_params.nthreads = 4;

    let bubble_order = path_linear_sgd(Arc::new(super_graph), sgd_params);

    bubble_order
}

// 4. FLATTEN: Convert bubble chain order to final node order
fn flatten_bubble_order(
    bubble_order: Vec<usize>,
    bubbles: &[Bubble],
    bubble_node_orders: &HashMap<usize, Vec<usize>>
) -> Vec<usize> {
    let mut final_order = Vec::new();

    for bubble_id in bubble_order {
        let bubble = &bubbles[bubble_id];
        let node_order = &bubble_node_orders[&bubble_id];
        final_order.extend(node_order);
    }

    final_order
}

// Main entry point
pub fn povu_guided_ygs_sort(
    graph: &mut BidirectedGraph,
    povu_decomposition: &str,
    params: &YgsParams
) {
    // Step 1: Parse povu decomposition
    let chains = parse_povu_decomposition(povu_decomposition);

    // Step 2: Topologically sort within each bubble
    let mut bubble_node_orders = HashMap::new();
    for chain in &chains {
        for bubble in &chain.bubbles {
            let sorted_nodes = topologically_sort_bubble(bubble, graph);
            bubble_node_orders.insert(bubble.id, sorted_nodes);
        }
    }

    // Step 3: Use SGD to order bubble chains globally
    let bubble_chain_order = order_bubble_chains(&chains, graph);

    // Step 4: Flatten to final node order
    let final_order = flatten_bubble_order(bubble_chain_order, &all_bubbles, &bubble_node_orders);

    // Step 5: Apply final ordering to graph
    apply_ordering(graph, &final_order);

    // Step 6: Groom (existing)
    groom(graph, params);

    // Step 7: Topological sort heads (existing)
    topological_sort_heads(graph);
}
```

## Expected Impact

### A-3105 Fix:

**Before (SGD only):**
- Max backward jump: -680 positions (43.9% of graph)
- 5 catastrophic edges
- 33-86% backward path steps
- RMSE: 21,959 bp

**After (Povu-guided):**
- Bubble 171 nodes kept together (687, 689, 1367-1548)
- No inter-bubble catastrophic edges (topologically sorted within)
- Expected max jump: <20 positions
- Expected RMSE: ~750 bp (same as working graphs)

### Benefits:

1. **Eliminates catastrophic failures** - Local topo sort prevents bubble-internal backward edges
2. **Preserves SGD strengths** - Global path coherence maintained via bubble-chain SGD
3. **Deterministic within bubbles** - No random initialization issues for local structure
4. **Better than ODGI** - ODGI doesn't use bubble decomposition, has same failures

### Potential Drawbacks:

1. **Dependency on povu** - Need to integrate povu library or call binary
2. **Performance overhead** - Additional decomposition step (should be fast)
3. **Complexity** - More moving parts than pure SGD

## Implementation Plan

### Phase 1: Prototype (Python script)
- Parse povu `.pvst` format
- Manually order A-3105 using bubble structure
- Measure if catastrophic edges eliminated
- **Validation metric:** Max edge jump should drop from 680 to <20

### Phase 2: Rust Integration
- Add povu as dependency (C++ library, need FFI bindings)
  - Alternative: Call povu binary as subprocess
- Implement `.pvst` parser in Rust
- Implement bubble-aware topological sort
- Integrate with existing Ygs pipeline

### Phase 3: Testing
- Run on full HLA-zoo suite
- Compare RMSE vs pure Ygs
- Ensure working graphs (B-3106, C-3107) not degraded
- Benchmark performance overhead

### Phase 4: Optimization
- Cache povu decomposition (only decompose once per graph)
- Parallelize bubble-internal topo sorts
- Tune SGD parameters for bubble-chain ordering

## File Locations

- **Povu repository:** `~/povu/`
- **Povu binary:** `~/povu/bin/povu`
- **Test GFA:** `/tmp/A-3105_before_odgi.gfa`
- **Povu output:** `/tmp/povu_test/1.pvst`
- **SeqRush Ygs impl:** `/home/erik/seqrush/src/ygs_sort.rs`
- **Path-SGD impl:** `/home/erik/seqrush/src/path_sgd.rs`

## References

- **Povu paper:** (TBD - check povu README for citations)
- **ODGI Path-SGD:** `~/odgi/src/algorithms/path_sgd.cpp`
- **Catastrophic edge analysis:** `/tmp/long_edge_fix_summary.md`
- **Detailed investigation:** `/tmp/analyze_a3105_detailed.py`

## Alternative Approaches Considered

### 1. Path-based initialization
**Idea:** Initialize SGD from average path positions instead of node IDs
**Rejected:** User correctly noted "isn't the Pathguided SGD like learning the embedding?" - SGD should learn correct layout regardless of initialization. The problem is local minima, not initialization.

### 2. Increase SGD iterations
**Tried:** 100 iterations with full convergence
**Result:** Still produces catastrophic edges - trapped in local minimum

### 3. Better topological sort
**Tried:** Sort all nodes (not just heads)
**Result:** Doesn't fix bubble-internal conflicts from SGD misplacement

### 4. Multi-temperature SGD
**Idea:** Simulated annealing approach
**Not tried:** Would add complexity; povu approach more principled

## Conclusion

The povu-guided hybrid approach combines the best of both worlds:
- **Topological sorting** ensures correctness within local structure (bubbles)
- **Path-Guided SGD** handles global layout and path coherence

This addresses a fundamental limitation of ODGI's pure SGD approach and could become a significant improvement for pangenome graph sorting.
