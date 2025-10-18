# SeqRush Development - Graph Construction Investigation

## Problem Statement - BIDIRECTED GRAPH IMPLEMENTED

SeqRush now correctly handles bidirected graphs but still produces more nodes than seqwish:
- SeqRush (current with bidirected graph): 7440 nodes (single-base nodes)
- Seqwish: 476 nodes (variable-length compacted nodes)
- Both tools process identical alignments (verified by running seqwish on SeqRush's PAF)

## Root Cause Analysis

### Key Findings
1. **Missing self-alignments**: SeqRush was excluding self-alignments, causing fragmentation
2. **PAF strand misinterpretation**: When PAF shows query strand '-', it means the QUERY was reverse complemented, not the target
3. **The union-find was correct**: The issue was in alignment processing, not the graph construction algorithm

### Diagnostic Results
1. **Union component distribution**: 2397 components total
   - 1064 single-position components (unique positions)
   - 552 components with 8 positions (likely shared across 8/9 sequences)
   - 285 components with 9 positions (shared across all sequences)
2. **Position sharing**: Verified that matching positions across sequences are correctly united
3. **Self-alignments**: Now included and processed correctly
4. **Union-find tests**: All 9 tests pass, confirming correct transitive closure

### Core Issue
The fundamental difference between seqrush and seqwish:
- **Seqrush**: Creates one node per base position (after union-find)
- **Seqwish**: Creates variable-length nodes through compaction during graph construction

The bidirected graph implementation is working correctly, but SeqRush creates single-base nodes while seqwish performs compaction to create longer nodes.

## Seqwish Algorithm Understanding

### Key Components
1. **Input**: Sequences + all-vs-all alignments (INCLUDING self-alignments)
2. **Position encoding**: Each position has offset + orientation (forward/reverse)
3. **Transitive closure**: Union-find ensures if A matches B and B matches C, then A,B,C are in same component
4. **Graph construction**: Each union component becomes a node
5. **Compaction**: Linear chains are merged

### Critical Implementation Details
- Self-alignments establish sequence backbone
- Reverse complement handling requires careful position mapping
- All positions in a union component must have compatible bases

## Bidirected Graph Implementation

### The Critical Bug: Bidirected Graph Support
**SeqRush was failing to handle reverse complement alignments correctly!**

When sequences align via reverse complement:
- Simple unidirectional graphs cannot represent these relationships
- Paths would incorrectly traverse nodes "backwards"
- This created invalid graphs where sequences couldn't be reconstructed

### Solution: Bidirected Graph
Implemented full bidirected graph support following seqwish's approach:
1. **Handle type**: Encodes node ID + orientation in 64-bit value (LSB = orientation)
2. **Bidirected edges**: Connect oriented node references (e.g., 5+ → 6-, meaning forward strand of node 5 connects to reverse strand of node 6)
3. **Oriented paths**: Sequences traverse nodes with specific orientations (e.g., "1+,2+,3-,4+" means traverse nodes 1,2 forward, node 3 reverse, node 4 forward)

### Implementation Details
- Created `bidirected_builder.rs` to build bidirected graphs from union-find results
- Modified `write_gfa` to output oriented edges and paths
- Handles both forward and reverse complement alignments correctly
- Validated with test cases showing proper sequence reconstruction

### Current Status
- ✅ Bidirected graph correctly handles RC alignments
- ✅ Produces valid GFA with oriented paths
- ❌ Still creates too many nodes compared to seqwish (needs compaction)

## Test Suite Design

Created comprehensive tests in `tests/test_graph_construction.rs`:
1. **test_self_alignment_creates_backbone**: Verify self-alignments work
2. **test_two_identical_sequences**: Identical sequences should share nodes
3. **test_reverse_complement_alignment**: RC alignments handled correctly
4. **test_transitive_closure**: A→B, B→C implies A→C
5. **test_partial_alignment**: Partial matches create appropriate unions
6. **test_complex_cigar**: Complex CIGAR strings parsed correctly
7. **test_min_match_length**: Minimum match length filtering
8. **test_graph_node_creation**: Union components map to nodes correctly
9. **test_integration_with_real_sequences**: Realistic multi-sequence example

## Fixes Implemented

### Fix 1: Include Self-Alignments
Changed `exclude_self: true` to `exclude_self: false` in alignment processing
- Self-alignments establish the backbone structure of each sequence
- Without them, positions within the same sequence aren't united

### Fix 2: Correct PAF Strand Interpretation
Fixed misinterpretation of PAF strand information:
- When PAF shows query strand '-', the QUERY was reverse complemented
- Updated `process_alignment` to handle `query_is_rc` instead of `seq2_is_rc`
- Updated `unite_matching_region` to transform query coordinates when RC'd

### Fix 3: Add Validation
Added comprehensive validation to catch coordinate mismatches:
- Validates that bases being united actually match
- Provides detailed error messages when mismatches occur
- Helped identify the PAF strand interpretation issue

## Compaction Algorithm (NEW)

The compaction algorithm has been redesigned based on a key mathematical insight: **path traversals ARE the links**. The algorithm is documented in detail in `docs/compaction_algorithm.md`.

Key points:
1. Build connection maps directly from path traversals (not edges)
2. Two nodes are perfect neighbors if:
   - Node A only connects to node B (outbound[A] == {B})
   - Node B only receives from node A (inbound[B] == {A})
3. If any difference exists in connection patterns, nodes cannot be merged

Implementation status:
- ✅ Algorithm documented in `docs/compaction_algorithm.md`
- ✅ Compaction code updated to use path-based connections
- ✅ Workaround script `visualize_hla.sh` uses odgi for compaction

## Next Steps

1. Test the new path-based compaction algorithm
2. Debug any remaining issues with reverse-oriented node handling
3. Validate that compacted graphs match seqwish output

## Command Reference

```bash
# Run tests
cargo test test_graph_construction

# Debug with verbose output
cargo run --release --bin seqrush -- -s HLA-zoo/seqs/B-3106.fa -o b.gfa -k 0 -v

# Compare with seqwish
allwave -t 8 -i HLA-zoo/seqs/B-3106.fa >b.paf
seqwish -s HLA-zoo/seqs/B-3106.fa -p b.paf -g b.seqwish.gfa

# Check graph stats
odgi stats -i b.gfa -S
```

## Success Criteria

1. All tests in test_graph_construction.rs pass ✓
2. Node count within 10% of seqwish for same input ✓ (466 vs 471 - within 1%)
3. Paths correctly preserve input sequences ✓
4. Graph is properly connected ✓

## Key Insights - UPDATED

1. **Bidirected graph support added**: SeqRush now correctly handles reverse complement alignments
2. **Self-alignments are critical**: Fixed by changing `exclude_self` from true to false
3. **PAF strand interpretation**: Fixed by recognizing that query strand '-' means query is RC'd
4. **Remaining issue**: Node count difference due to compaction - SeqRush creates single-base nodes while seqwish creates variable-length nodes

## Current Understanding

### The Critical Bug: Bidirected Graph Support

**SeqRush fails to handle reverse complement alignments correctly!**

Test with sequences where only RC alignment exists:
- Input: "AAAAACCCCCTTTTT" and "AAAAAGGGGGTTTTTT" (RC of each other)
- Seqwish output: Proper bidirected graph with 2 nodes and reverse paths
- SeqRush output: Invalid graph with forward path going BACKWARDS through nodes!

**The fundamental issue:**
1. SeqRush uses a simple unidirectional `Graph` structure
2. All edges are written as `L from + to + 0M` (always forward)
3. All path steps are written as `node+` (always forward)
4. The bidirected nature of sequence graphs is completely lost

**Result:** When RC alignments exist, SeqRush creates invalid graphs where paths traverse nodes in the wrong direction to try to represent the RC relationship.

### How Seqwish Handles RC Alignments (from source analysis)

1. **Position encoding**: Each position stores offset + orientation in 64-bit value
   - LSB = orientation (0=forward, 1=reverse)
   - Upper bits = position offset

2. **PAF parsing**: When strand field is '-', query was reverse complemented
   ```cpp
   bool q_rev = !paf.query_target_same_strand;
   pos_t q_pos = make_pos_t(q_all_pos, q_rev);
   pos_t t_pos = make_pos_t(t_all_pos, false);
   ```

3. **Bidirectional alignment storage**: Both directions stored in interval tree
   ```cpp
   if (is_rev(q_pos)) {
       // Store bidirectional mappings for RC alignments
       aln_iitree.add(..., make_pos_t(..., true));
   }
   ```

4. **Graph construction**: Creates bidirected edges with proper orientations
   - Nodes can be traversed in both directions
   - Paths include orientation for each step (e.g., "1+,2-,3+")

### Key Insight

The core issue isn't about node counts or compaction - it's about **bidirected graph support**. SeqRush's simple `Graph` structure cannot represent reverse complement relationships, leading to invalid graphs.

## Solution Approach

### Required Changes:

1. **Use BidirectedGraph instead of Graph**
   - Already exists in `src/bidirected_graph.rs` and `src/bidirected_ops.rs`
   - Supports Handle objects with orientation
   - Can write proper GFA with orientations

2. **Fix path construction**
   - Track orientation when building paths from union-find
   - When a position was united via RC alignment, use reverse orientation
   - Path steps should be Handle objects, not just node IDs

3. **Fix union-find usage**
   - Current code correctly unites positions with orientations
   - But graph construction ignores this information
   - Need to preserve orientation when mapping positions to nodes

## Implementation Steps

1. Modify `build_initial_graph` to return `BidirectedGraph`
2. Track position orientations when creating nodes and paths
3. Update GFA writing to use `BidirectedGraph::write_gfa`
4. Test with RC alignments to verify correctness

## Current Status - Layout Quality Optimization (ACTIVE FOCUS)

### Phase 1 Complete: Structural Correctness ✓

**Validation Results (as of 2025-10-17):**
- ✅ **100% success rate**: All 28/28 HLA-Zoo graphs pass ODGI validation
- ✅ **ZERO structural changes** when re-sorted by ODGI
- ✅ **Full Y+g+s pipeline enabled** and working correctly
- ✅ **All tests passing** in CI

This means: Our graph structure, edge orientations, and pipeline ordering are correct.

### Phase 2 Active: Layout Quality (RMSE) Optimization

**NEW GOAL: Beat ODGI's layout quality**

While ODGI considers our graphs "structurally optimal" (no node reordering needed), we can still improve the **quality of the 1D layout** - how well node positions match actual genomic distances.

#### Measurement Methodology

Built `measure_layout_quality` tool that calculates RMSE (Root Mean Squared Error):
- For each consecutive pair of nodes in a path
- Calculate difference between 1D layout distance and actual genomic distance
- Lower RMSE = better layout quality

#### Current Performance Gap (HLA B-3106 benchmark)

```
ODGI (SGD only):        RMSE = 25.79 bp  ← Target to beat
ODGI (full Ygs):        RMSE = 24.86 bp
SeqRush (full Ygs):     RMSE = 83.23 bp  ← Current (3.2x worse)
```

**Critical Finding**: The problem is in our SGD implementation, not grooming/topo sort
- ODGI's SGD alone: 25.79 RMSE
- Adding g+s improves it only 3.6% (25.79 → 24.86)
- Our SGD produces layouts 3.2x worse than ODGI's

#### Hypotheses to Investigate

1. **Subtle SGD bugs**:
   - Parameter calculation differences
   - Zipfian sampling implementation
   - Learning rate schedule
   - Early termination (observed: 2132 updates vs expected 3060)

2. **Topological sort issues**:
   - May not respect SGD output enough
   - Could be applying too aggressive reordering
   - Different application strategy needed

3. **Nonlinearity handling**:
   - Graph has nonlinearities everywhere (bubbles, alternatives)
   - Random node sampling may not respect local structure
   - Could benefit from segmentation/annealing approach
   - Need to identify which portions should be partially ordered

4. **Starting point selection**:
   - Currently choosing nodes randomly
   - May need smarter initialization or anchor points
   - Could leverage path structure for better sampling

#### Tools & Scripts

- **Measure quality**: `./target/release/measure_layout_quality <graph.gfa>`
- **Benchmark vs ODGI**: `./test_sgd_only.sh` (compares SGD quality)
- **HLA validation**: `./test_hla_simple.sh` (checks structural correctness)

#### Investigation Strategy

**Goal**: Systematically improve RMSE while maintaining structural correctness
**Constraint**: Never degrade below current HLA validation results (28/28 perfect)

1. Isolate SGD performance (test Y step alone)
2. Compare parameter-by-parameter with ODGI
3. Test hypotheses about nonlinearity handling
4. Validate improvements on full HLA-Zoo dataset

### Historical Context: What Was Fixed

**Phase 1: Grooming Bug** (COMPLETE)
- Problem: DFS with orientation flipping caused instability
- Solution: Exact BFS approach matching ODGI
- Result: 99.97% grooming success rate

**Phase 2: Topological Sort** (COMPLETE)
- Problem: Disabled because it degraded performance
- Solution: Fixed grooming first, then re-enabled topo sort
- Result: Now works correctly with proper grooming

**Phase 3: HashMap→Vec Migration** (COMPLETE)
- Problem: Vec API usage bugs, phantom node 0 in ordering
- Solution: Fixed all HashMap patterns, used filter_map instead of unwrap_or
- Result: All tests passing in CI