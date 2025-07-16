# SeqRush Development - Graph Construction Investigation

## Problem Statement - RESOLVED

SeqRush was producing graphs with significantly more nodes than seqwish when processing the same input:
- SeqRush (original): 1785 nodes (after compaction)
- SeqRush (after fixes): 466 nodes ✓
- Seqwish: 471 nodes
- Both tools use identical PAF alignments

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
- **Seqrush**: Creates one node per union component (2397 nodes)
- **Seqwish**: Builds a "graph sequence" through transitive closure, then identifies nodes based on graph topology (471 nodes)

The union-find is working correctly, but we need a different approach to node creation that considers the graph structure, not just union components.

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

## Next Steps

1. Run the test suite to identify specific failures
2. Fix union-find implementation based on test results
3. Add more detailed logging to trace union operations
4. Compare intermediate results with seqwish at each stage
5. Iterate until all tests pass and node counts match seqwish

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

## Key Insights - RESOLVED

1. **Union-find was working correctly**: The algorithm was sound, issue was in alignment processing
2. **Self-alignments are critical**: Fixed by changing `exclude_self` from true to false
3. **PAF strand interpretation**: Fixed by recognizing that query strand '-' means query is RC'd
4. **Final result**: 466 nodes (SeqRush) vs 471 nodes (seqwish) - within 1% ✓

## Conclusion

The original hypothesis about needing a "graph sequence" approach was incorrect. SeqRush's algorithm of creating one node per union component is fundamentally sound and matches seqwish's approach. The issues were in the alignment processing:

1. **Missing self-alignments** caused excessive fragmentation
2. **Incorrect PAF strand interpretation** created spurious matches and a tangled graph

With these fixes, SeqRush now produces graphs nearly identical to seqwish (466 vs 471 nodes), demonstrating that the core algorithm was correct all along.