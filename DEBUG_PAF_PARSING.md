# PAF Parsing Bug Analysis

## Critical Bug Found: Incorrect RC Alignment Position Calculation

### The Problem

SeqRush is misinterpreting reverse complement alignments by reverse complementing the sequence data instead of just tracking orientation and adjusting position traversal.

### How Seqwish Handles RC Alignments (Correct)

From `seqwish-local/src/alignments.cpp` lines 45-50:

```cpp
bool q_rev = !paf.query_target_same_strand;
size_t q_all_pos = (q_rev ? seqidx.pos_in_all_seqs(query_idx, paf.query_end, false) - 1
                    : seqidx.pos_in_all_seqs(query_idx, paf.query_start, false));
size_t t_all_pos = seqidx.pos_in_all_seqs(target_idx, paf.target_start, false);
pos_t q_pos = make_pos_t(q_all_pos, q_rev);
pos_t t_pos = make_pos_t(t_all_pos, false);
```

Key points:
1. **Starting position**: When `q_rev` is true, start at `query_end - 1`, NOT `query_start`
2. **Position encoding**: Positions include orientation: `make_pos_t(offset, is_reverse)`
3. **Traversal direction**: `incr_pos()` moves BACKWARDS when position is reverse-oriented

From `seqwish-local/src/pos.cpp`:
```cpp
void incr_pos(pos_t& pos) {
    if (is_rev(pos)) {
        pos -= 2;  // Move BACKWARDS through sequence
    } else {
        pos += 2;  // Move forward normally
    }
}
```

### How SeqRush Currently Handles RC Alignments (WRONG)

From `src/seqrush.rs` lines 714-718:

```rust
// If query was reverse complemented for alignment, we need to compare against the RC of seq1
let seq1_data = if query_is_rc {
    seq1.reverse_complement()  // WRONG!
} else {
    seq1.data.clone()
};
```

Then at lines 778-785:
```rust
let global_offset1 = if query_is_rc {
    // For RC alignments, need to map RC coordinates back to forward
    let rc_local_pos = pos1 + k;
    let forward_local_pos = seq1.data.len() - 1 - rc_local_pos;
    seq1.offset + forward_local_pos
} else {
    seq1.offset + pos1 + k
};
```

**The problem**:
- We're RC'ing the entire sequence and then processing forward
- This creates DIFFERENT offsets than seqwish
- The position mapping is fundamentally wrong

### Correct Approach (from seqwish source analysis)

Seqwish uses oriented positions throughout:

1. **Position encoding**: `pos_t = (offset << 1) | orientation_bit`
2. **Base access**: When accessing a reverse-oriented position, automatically RC the base:
   ```cpp
   char c = at(offset(pos));
   if (is_rev(pos)) {
       c = dna_reverse_complement(c);
   }
   ```
3. **Position increment**: When position is reverse-oriented, decrement offset:
   ```cpp
   void incr_pos(pos_t& pos) {
       if (is_rev(pos)) {
           pos -= 2;  // Move backwards!
       } else {
           pos += 2;
       }
   }
   ```

When PAF shows query strand is '-':

1. **Don't RC the sequence data** - keep original sequence
2. **Start from query_end - 1** - the end of the aligned region
3. **Track orientation separately** - query positions are reverse-oriented
4. **RC bases on access** - when comparing, RC the query base if query_is_rc
5. **Move backwards** - when query_is_rc, decrement position instead of increment

### Example

Given:
- Query sequence: "ATCGATCG" (positions 0-7)
- PAF alignment: query_start=2, query_end=6, strand='-'
- CIGAR: 4M

Seqwish processing:
1. Start at position 5 (query_end - 1)
2. Process CIGAR positions: 5, 4, 3, 2 (backwards!)
3. Compare RC'd bases at these positions: RC(seq[5]), RC(seq[4]), RC(seq[3]), RC(seq[2])
4. Unite with target positions: 0, 1, 2, 3

SeqRush current (WRONG):
1. RC entire sequence: "CGATCGAT"
2. Start at position 0 in RC'd sequence
3. Process positions 0, 1, 2, 3 in RC'd sequence
4. These map to DIFFERENT original offsets!

### Impact

This bug causes incorrect position unification:
- Positions that should be united are NOT united
- Creates extra graph components (2715 nodes instead of 453)
- Graph structure doesn't match seqwish

### Fix Required

Modify `process_alignment` in `src/seqrush.rs`:
1. Remove the `seq1.reverse_complement()` call
2. When `query_is_rc` is true:
   - Access query sequence in REVERSE order
   - RC each base on-the-fly for comparison
   - Calculate correct original offsets
3. Ensure position offsets match seqwish exactly

## Test Case

Create minimal test with two sequences where one is RC of the other:
- Forward: "ATCGATCG"
- Reverse: "CGATCGAT"

Expected: Both sequences should share the same nodes (with different orientations)
Current: Creates separate nodes due to position offset mismatch
