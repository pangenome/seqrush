# Seqwish vs Seqrush Discrepancy Analysis

## Problem Statement
- Seqwish produces 453 nodes from x.clean.paf (72 alignments)
- Seqrush produces 2715 nodes from the same PAF
- Both tools should produce identical graphs from identical inputs

## Key Findings

### 1. Position Encoding (IDENTICAL)
Both seqwish and seqrush use the same position encoding scheme:
- Position = (offset << 1) | orientation_bit
- LSB=0 means forward strand (+)
- LSB=1 means reverse strand (-)

This is confirmed by comparing:
- seqwish: `seqwish-local/src/pos.cpp` lines 13-20
- seqrush: `src/pos.rs` lines 10-12

### 2. CIGAR Parsing (IDENTICAL)
Both tools parse CIGAR strings identically:
- Match operations: M, =, X
- Insertions: I (advance query only)
- Deletions: D (advance target only)

This is confirmed by comparing:
- seqwish: `seqwish-local/src/alignments.cpp` lines 52-114
- seqrush: `src/seqrush.rs` lines 760-950

### 3. CRITICAL DIFFERENCE: Self-Mapping Guard

**Seqwish has a critical guard that seqrush is missing:**

```cpp
// seqwish-local/src/alignments.cpp line 84
if (query_base == target_base
    && query_base != 'N'
    && offset(q_pos) != offset(t_pos)) {  // <-- CRITICAL GUARD
```

**Seqrush does NOT have this guard:**

```rust
// src/seqrush.rs line 788
if base1 == base2 {
    // No offset check - unites even when offsets are equal!
```

### 4. What This Guard Does

The guard `offset(q_pos) != offset(t_pos)` prevents uniting a position with itself in the concatenated sequence space.

**Why this matters:**
1. All sequences are concatenated into a single sequence space
2. Each position has a global offset in this concatenated space
3. During alignment processing, the same global offset could theoretically appear as both query and target
4. Without the guard, we'd unite a position with itself (which is a no-op but might affect the union-find structure)

**However**, after examining x.clean.paf:
- There are NO self-alignments (where query_name == target_name)
- All 72 alignments are between different sequences
- So this shouldn't be the issue... unless...

### 5. Hypothesis: N-base Filtering

**Seqwish also filters N bases:**

```cpp
// seqwish-local/src/alignments.cpp line 83
if (query_base == target_base
    && query_base != 'N'  // <-- Filters N bases
    && offset(q_pos) != offset(t_pos))
```

**Seqrush does NOT filter N bases - it only checks if bases match:**

```rust
// src/seqrush.rs line 788
if base1 == base2 {
    // No N filtering!
```

This could cause differences if sequences contain N bases that seqwish ignores but seqrush unites.

### 6. Testing the Hypothesis

To test if N-filtering is the issue:
1. Count N bases in the input sequences
2. Check if any alignments contain N bases in their match regions
3. Compare node counts with/without N filtering

## Recommended Fixes

### Fix 1: Add Self-Mapping Guard (CRITICAL)

Add the offset check to seqrush's alignment processing:

```rust
// In src/seqrush.rs around line 788
if base1 == base2 && base1 != b'N' {
    // Calculate global offsets
    let global_offset1 = seq1.offset + (if query_is_rc {
        seq1.data.len() - 1 - (pos1 + k)
    } else {
        pos1 + k
    });
    let global_offset2 = seq2.offset + pos2 + k;

    // Only unite if different global positions
    if global_offset1 != global_offset2 {
        // Continue with match processing...
    }
}
```

### Fix 2: Filter N Bases

Add N-base filtering:

```rust
if base1 == base2 && base1 != b'N' {
    // Process match...
}
```

## Next Steps

1. ✅ Document discrepancies (this file)
2. Create targeted test cases
3. Implement fixes
4. Verify node count matches seqwish (should be 453)
