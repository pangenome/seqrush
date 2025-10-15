# SGD Path Position Calculation Bug - Reverse Handles

## Problem

Path-guided SGD produces suboptimal layouts with "bubbles" - nodes that should be adjacent in paths are placed far apart. Topological sort would eliminate these, indicating SGD is not calculating correct positions.

## Root Cause

**Location**: `src/path_sgd_exact.rs` lines 528-574

The bug: When calculating `pos_in_path_b`, the code searches for `handle_b` by comparing full handles:

```rust
let pos_in_path_b = {
    let mut pos = 0;
    for (pid, h) in &path_index_steps {
        if *pid == path_id {
            if *h == handle_b {  // ← Compares node ID + orientation
                break;
            }
            if let Some(node) = graph_clone.nodes.get(&h.node_id()) {
                pos += node.sequence.len();
            }
        }
    }
    pos
};
```

**The Problem**: This finds the FIRST occurrence of `handle_b` in the path. If a node appears multiple times (valid for structural variation), or if we're looking for the wrong rank, we calculate the wrong position.

**However**, my attempted fix (using `b_rank` instead) made things WORSE, suggesting the current implementation is closer to correct.

## Alternative Hypothesis

The real issue may be that **path positions are calculated correctly**, but the SGD is not converging because:

1. Initial positions (seeded by node ID order) create massive separations
2. Not enough iterations to overcome poor initialization
3. Conflicting constraints from paths with mixed orientations

## Diagnostic Evidence

Created `sgd_diagnostics` tool showing:
- Adjacent nodes in paths (1-40bp apart) are placed 100x-3976x apart in final layout
- Example: Node 1 and 205 are 1bp apart in path, but 845 positions apart in layout
- ALL paths show this problem (not just RC paths)

## Next Steps

1. Create focused tests for simple cases (chains with RC edges)
2. Verify SGD can achieve perfect layout on trivial graphs
3. If tests fail, the bug is confirmed in position calculation or SGD convergence
4. Fix and verify

## Test Cases Needed

1. **Simple chain forward**: `1+->2+->3+` should place nodes sequentially
2. **Simple chain reverse**: `1-->2-->3-` should place nodes sequentially
3. **Mixed orientation chain**: `1+->2-->3+` should place nodes sequentially
4. **Simple bubble with RC**:
   ```
   Path A: 1+->2+->3+
   Path B: 1+->2-->3+  (traverses node 2 in reverse)
   ```
5. **Repeated node**: `1+->2+->1+` should handle revisiting correctly

Each test should verify SGD produces perfect ordering (no displacement).
