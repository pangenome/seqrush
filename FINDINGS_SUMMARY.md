# Seqwish vs Seqrush Parity Investigation - Findings Summary

## Problem
- **Target**: seqrush should produce 453 nodes from x.clean.paf (matching seqwish)
- **Current**: seqrush produces 2715 nodes from the same input
- **Gap**: 6x more nodes than seqwish

## Investigation Results

### Phase 1: Algorithm Analysis

#### Hypothesis 1: Self-Mapping Guard Missing ❌ NOT THE ROOT CAUSE
**Finding**: Seqwish has a guard `offset(q_pos) != offset(t_pos)` that seqrush was missing.

**Fix Implemented**: Added the guard to seqrush (src/seqrush.rs lines 777-810):
```rust
// Calculate global offsets for self-mapping guard
let global_offset1 = if query_is_rc {
    let rc_local_pos = pos1 + k;
    let forward_local_pos = seq1.data.len() - 1 - rc_local_pos;
    seq1.offset + forward_local_pos
} else {
    seq1.offset + pos1 + k
};
let global_offset2 = seq2.offset + pos2 + k;

// CRITICAL SEQWISH PARITY FIX:
// Only unite if:
// 1. Bases match
// 2. Base is not 'N' (seqwish filters N bases)
// 3. Global offsets are different (guard against self-mapping)
if base1 == base2
    && base1 != b'N'
    && global_offset1 != global_offset2 {
```

**Result**: Still produces 2715 nodes after the fix.

**Conclusion**: While this fix is correct and matches seqwish's logic, it doesn't explain the 6x difference in node count.

#### Hypothesis 2: N-Base Filtering ❌ NOT THE ISSUE
**Finding**: Seqwish filters N bases (`query_base != 'N'`)

**Test**: Checked input sequences - they contain ZERO N bases.

**Conclusion**: Not applicable to this dataset.

#### Hypothesis 3: Different Graph Compaction ✅ **LIKELY ROOT CAUSE**
**Finding**: Fundamental difference in node structure:

**Seqwish output** (x.seqwish.gfa):
```
S	1	AT
S	2	TCTGGAA
S	3	G
S	4	GTTCTCAGGTCTTTATTTGCTCT
```
Nodes are **multi-base sequences** (variable length).

**Seqrush output** (test_fixed.gfa):
```
S	1	C
S	2	G
S	3	G
S	4	G
```
Nodes are **single bases** (always length 1).

**Node Count Math**:
- Total sequence length: 30,751 bp
- Seqwish nodes: 453 (average ~68 bp per node)
- Seqrush nodes: 2,715 (much smaller nodes)

### Phase 2: Architectural Difference

**Seqwish uses a TWO-PHASE approach**:
1. **Phase 1** (`alignments.cpp`): Store matches in interval tree
2. **Phase 2** (`transclosure.cpp`): Apply transitive closure with union-find, then compact

**Seqrush uses a ONE-PHASE approach**:
1. **Single phase**: Unite positions directly during alignment processing
2. Compaction happens after graph construction (but may not be working correctly)

### Phase 3: Root Cause Analysis

The 2715 nodes in seqrush suggests that:
1. **Union-find is creating too many separate components** - positions that should be united into the same node are remaining separate
2. **Compaction is not running or not working** - even if union-find created the right components, compaction should merge linear chains

**Evidence**:
- Seqrush with `--no-compact` flag: 2715 nodes ❌
- Seqrush WITHOUT flag (default=compact): 2715 nodes ❌
- Seqwish (always compacts): 453 nodes ✅

This indicates **compaction is not running at all** in seqrush, OR the union-find is creating far too many separate components.

## Key Code Locations

### Seqwish
- PAF parsing: `seqwish-local/src/alignments.cpp` lines 19-116
- Transitive closure: `seqwish-local/src/transclosure.cpp` lines 490-540
- Union-find usage: `transclosure.cpp` line 504: `disjoint_sets.unite(q_curr_rank(j), q_curr_rank(offset(p)))`

### Seqrush
- PAF parsing: `src/seqrush.rs` lines 641-998 (`process_alignment`)
- Union-find: `src/bidirected_union_find.rs` lines 60-98 (`unite_matching_region`)
- Graph building: `src/seqrush_bidirected.rs` or similar

## Next Steps

1. **Verify union-find is working**:
   - Add debug output to count union-find components
   - Compare with seqwish's transitive closure output

2. **Check compaction**:
   - Verify bidirected compaction is actually running
   - Check if compaction algorithm has bugs

3. **Compare intermediate results**:
   - Dump union-find state from both tools
   - Compare which positions are united

4. **Simplify test case**:
   - Create minimal example (2-3 sequences) where difference is clear
   - Debug step-by-step through both tools

## Status

- ✅ Self-mapping guard implemented (correct but not sufficient)
- ✅ N-base filtering implemented (correct but not applicable)
- ❌ Node count still 2715 instead of 453
- 🔍 Need to investigate union-find component count and compaction

## Files Modified

- `/home/erik/seqrush/src/seqrush.rs` - Added self-mapping guard and N-filtering
- `/home/erik/seqrush/ANALYSIS_SEQWISH_PARITY.md` - Initial analysis
- `/home/erik/seqrush/tests/seqwish_parity_tests.rs` - Test cases (not all compilable yet)
- `/home/erik/seqrush/FINDINGS_SUMMARY.md` - This document
