# Seqwish Algorithm Analysis

## Overview

Seqwish implements a variation graph inducer that converts pairwise alignments between sequences into a variation graph. The algorithm is designed to handle very large inputs using disk-backed sorts and memory-efficient data structures.

## Key Components

### 1. Input Processing

- **Sequences (Q)**: Input sequences are concatenated into a single sequence Q
- **Sequence Index**: Built using a compressed suffix array (CSA) that maps sequence names to offsets in Q
- **Alignments (A)**: PAF format alignments are parsed and stored as pairs of bidirectional positions

### 2. PAF Alignment Processing (`alignments.cpp`)

The PAF parsing extracts exact matches from CIGAR strings:

```cpp
// Key features:
- Processes PAF alignments in parallel using multiple threads
- Extracts exact matches from CIGAR operations (M, =, X)
- Filters matches below minimum length threshold
- Supports match sparsification using hash-based sampling
- CRITICAL: Excludes self-alignments (offset(q_pos) != offset(t_pos))
- Stores bidirectional alignments in interval tree (aln_iitree)
```

Important observations:
- When query is on reverse strand (!paf.query_target_same_strand), positions are adjusted
- Each match is stored bidirectionally (query→target and target→query)
- Self-alignments are explicitly filtered out in the match processing

### 3. Transitive Closure (`transclosure.cpp`)

This is the core algorithm that builds the graph sequence:

#### Algorithm Steps:

1. **Batch Processing**: Processes input in chunks (default 1MB) to manage memory
2. **Union-Find**: Uses lock-free parallel union-find to compute transitive closures
3. **Graph Sequence Construction**: Builds output sequence S where each position represents a graph node

#### Key Data Structures:

- `q_seen_bv`: Tracks which positions in Q have been processed
- `q_curr_bv`: Current batch's positions being processed
- `DisjointSets`: Lock-free union-find implementation
- `node_iitree`: Maps graph sequence positions to input sequence positions
- `path_iitree`: Maps input sequence positions to graph sequence positions

#### Detailed Process:

1. **Chunk Selection**: Find next unprocessed chunk of bases
2. **Overlap Collection**: For each position in chunk, find all overlapping alignments
3. **Union-Find**: Unite all positions that align together
4. **Sort by Position**: Sort disjoint sets by minimum position
5. **Graph Emission**: Write bases to graph sequence S in sorted order

### 4. Union-Find Implementation (`dset64-gccAtomic.hpp`)

- Lock-free implementation using 128-bit atomic operations
- Supports concurrent find(), unite(), and same() operations
- Uses path compression and union by rank
- Stores parent (64 bits) and rank (64 bits) in single 128-bit atomic

### 5. Node Compaction (`compact.cpp`)

After building the graph sequence, nodes are compacted:

```cpp
// For each sequence:
- Find where it maps in the graph sequence
- Mark node boundaries at:
  - Start of each mapped region
  - End of each mapped region
  - Positions with multiple incoming/outgoing edges
```

This creates a bit vector `seq_id_bv` that marks node boundaries, allowing compression of unitigs into single nodes.

## Key Differences from SeqRush

### 1. Self-Alignment Handling
- **Seqwish**: Explicitly excludes self-alignments in alignment processing
- **SeqRush**: Was including self-alignments (fixed in CLAUDE.md)

### 2. Algorithm Approach
- **Seqwish**: Builds graph sequence S through transitive closure, then identifies nodes
- **SeqRush**: Creates nodes directly from union components

### 3. PAF Strand Interpretation
- **Seqwish**: Correctly handles query reverse complement when parsing PAF
- **SeqRush**: Had incorrect interpretation (fixed in CLAUDE.md)

### 4. Processing Model
- **Seqwish**: Batch-based processing with configurable batch size
- **SeqRush**: Processes all alignments at once

### 5. Memory Management
- **Seqwish**: Uses disk-backed data structures for scalability
- **SeqRush**: In-memory processing

## Critical Implementation Details

### 1. Position Encoding
- Positions encoded with offset + orientation (forward/reverse)
- `make_pos_t(offset, is_rev)` creates position
- `is_rev()`, `offset()` extract components

### 2. Bidirectional Storage
- Every alignment stored in both directions
- Ensures transitive closure captures all relationships

### 3. Batch-Based Union-Find
- Processes in chunks to limit memory usage
- Parallel overlap collection and union operations
- Careful synchronization using atomic operations

### 4. Graph Sequence Generation
- Positions sorted by disjoint set ID
- Within each set, sorted by minimum position
- Ensures deterministic output

## Important Discovery: Self-Alignment Exclusion

After analyzing the source code, I found that **seqwish explicitly excludes self-alignments** in its alignment processing:

```cpp
// In alignments.cpp, line 84:
if (query_base == target_base
    && query_base != 'N'
    && offset(q_pos) != offset(t_pos)) { // guard against self mappings
```

This contradicts the initial understanding that self-alignments were required. The algorithm establishes sequence structure through the batch-based transitive closure process, not through self-alignments.

## Conclusion

Seqwish's algorithm is fundamentally sound and matches SeqRush's corrected approach. The key insights are:

1. **One node per union component** is the correct approach
2. **Self-alignments are excluded** in alignment processing (contrary to initial belief)
3. **PAF strand interpretation** must account for query reverse complement
4. **Transitive closure** ensures all related positions are united
5. **Batch processing** allows handling of large datasets

The main difference is in implementation strategy: Seqwish uses a batch-based, disk-backed approach for scalability, while SeqRush uses in-memory processing. Both produce nearly identical results when implemented correctly.