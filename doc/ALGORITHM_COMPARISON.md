# SeqRush vs Seqwish Algorithm Comparison

## Key Finding

The core issue was that SeqRush was **excluding self-alignments**, which are critical for establishing the backbone structure of each sequence. Without self-alignments, positions within the same sequence don't get united, leading to excessive fragmentation.

## Algorithm Comparison

### Both SeqRush and Seqwish:

1. **Process alignments** to build a union-find data structure
2. **Create one initial position/node per union component** 
3. **Compact linear chains** into single nodes

### The Critical Difference: Self-Alignments

**SeqRush (before fix):**
- Excluded self-alignments (`exclude_self = true`)
- Result: 2397 union components → 1785 nodes after compaction

**SeqRush (after fix):**
- Includes self-alignments (`exclude_self = false`)
- Result: 5208 union components → 1108 nodes after compaction

**Seqwish:**
- Always includes self-alignments
- Result: ~5000 components → 471 nodes after compaction

## Why SeqRush Still Has More Nodes (1108 vs 471)

Even with self-alignments included, SeqRush produces more nodes than seqwish. This suggests differences in:

1. **Compaction algorithm effectiveness**
2. **Alignment parameters or filtering**
3. **CIGAR processing details**
4. **Handling of complex regions**

## The Graph Sequence Misconception

Initially, I thought seqwish used a fundamentally different "graph sequence" approach. However, examining the code reveals:

1. Seqwish DOES create one position in the graph sequence per union component
2. The "graph sequence" is just the linearized set of union components
3. Both tools use the same basic approach: union components → nodes → compaction

## Recommendations for Further Investigation

1. **Compare compaction algorithms** - SeqRush's compaction may be less aggressive
2. **Check alignment parameters** - Different scoring might produce different alignments
3. **Examine complex regions** - The HLA data has repetitive regions that might be handled differently
4. **Profile union component sizes** - Understanding the distribution might reveal insights

## Test Results Summary

| Test Case | SeqRush Initial | SeqRush Compacted | Seqwish |
|-----------|----------------|-------------------|---------|
| 3 identical 12bp sequences | 12 nodes | 1 node | 1 node |
| 9 HLA sequences (~3.3kb each) | 5208 nodes | 1108 nodes | 471 nodes |

The simple test case shows SeqRush's algorithm is fundamentally correct. The difference in the HLA case suggests the issue is in the details of alignment processing or compaction effectiveness rather than a fundamental algorithmic difference.