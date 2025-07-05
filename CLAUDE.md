# CLAUDE.md - SeqRush Implementation Notes

## Project Overview

SeqRush is a Rust reimplementation of seqwish's transitive closure algorithm for building pangenome graphs. The goal is to create a dynamic, parallel, in-memory implementation that can efficiently build and write subgraphs for any partition.

## Key Concepts from Seqwish

### 1. Transitive Closure Algorithm

The core idea: if position A matches position B, and position B matches position C, then A, B, and C should all be represented by the same node in the graph. This transitive relationship is computed using a union-find (disjoint set) data structure.

### 2. Lock-Free Union-Find

Seqwish uses a sophisticated lock-free union-find implementation:
- **dset64-gccAtomic.hpp**: Primary implementation using GCC's `__sync_bool_compare_and_swap`
- Stores parent and rank in a single 128-bit atomic value
- Implements path compression and union by rank
- Wait-free parallel operations

### 3. Batch Processing Algorithm

The transitive closure works in chunks:

1. **Chunk Selection**: Start at first unseen position, collect `transclose_batch_size` unseen bases
2. **Overlap Collection**: 
   - Use atomic queues for work distribution
   - Multiple threads explore overlaps in parallel
   - Each overlap is added to the work queue
3. **Union-Find Phase**:
   - Create a DisjointSets structure for the chunk
   - Unite positions that overlap according to matches
4. **Compression and Sorting**:
   - Find representative for each position
   - Sort by disjoint set ID
   - Renumber sets by minimum position
5. **Graph Emission**:
   - Process disjoint sets in order
   - Create one node per equivalence class
   - Handle repeats and minimum distance constraints

### 4. Match Processing

Matches are extracted from alignments and stored in interval trees:
- Bidirectional indexing (query→target and target→query)
- Support for both forward and reverse strand alignments
- Minimum match length filtering
- Optional sparsification using hash-based sampling

### 5. Graph Construction Pipeline

1. **Parse alignments** → Extract matches
2. **Build union-find** → Group equivalent positions
3. **Process positions** → Create nodes for equivalence classes
4. **Compact nodes** → Merge linear paths (optional)
5. **Derive links** → Connect adjacent nodes
6. **Output GFA** → Write graph in standard format

## Critical Implementation Details

### Position Encoding
- 64-bit positions with orientation bit (LSB)
- 0 = forward strand, 1 = reverse strand
- Remaining 63 bits store the offset

### Duplicate Node Prevention
The error "the link 1775+,1775+ is missing" indicates duplicate nodes in paths. Seqwish handles this with:
- `repeat_max`: Limits copies of input bases in graph
- `min_repeat_dist`: Prevents closure of nearby repeats
- Careful tracking of which positions have been processed

### Performance Considerations

1. **Memory Access Patterns**
   - Sequential access when possible
   - Batch processing for cache efficiency
   - Memory-mapped files for large datasets

2. **Parallelization Strategy**
   - Thread-local work queues
   - Atomic operations for synchronization
   - Lock-free data structures throughout

3. **Union-Find Efficiency**
   - Path compression during find operations
   - Union by rank for balanced trees
   - Atomic 128-bit operations for thread safety

## Current Issues in SeqRush

### 1. Performance Problems
- Using uf_rush crate which may not be as optimized as seqwish's custom implementation
- Sequential processing in some critical paths
- Not utilizing batch processing effectively

### 2. Algorithm Differences
- Not implementing the exact transitive closure algorithm from seqwish
- Missing proper handling of repeat sequences
- Incorrect handling of position mappings

### 3. Missing Features
- No support for minimum repeat distance
- No batch processing for memory efficiency
- Limited parallelization compared to seqwish

## Required Improvements

### 1. Union-Find Implementation
- Need to verify uf_rush performance characteristics
- Consider implementing custom lock-free union-find if needed
- Ensure atomic operations are as efficient as seqwish

### 2. Batch Processing
- Implement chunked processing of sequences
- Use parallel iterators more effectively
- Optimize memory access patterns

### 3. Duplicate Node Handling
- Implement proper repeat handling logic
- Track processed positions more carefully
- Ensure each position creates at most one node

### 4. Algorithm Correctness
- Follow seqwish's exact algorithm for transitive closure
- Implement proper position mapping
- Handle edge cases like self-alignments

## Next Steps

1. **Analyze uf_rush performance** - Benchmark against seqwish's union-find
2. **Implement batch processing** - Process sequences in chunks like seqwish
3. **Fix duplicate node issue** - Ensure proper path construction
4. **Optimize parallel processing** - Use rayon more effectively
5. **Add missing features** - Repeat handling, sparsification, etc.

## References

- Seqwish paper: https://doi.org/10.1093/bioinformatics/btac453
- Seqwish code: https://github.com/ekg/seqwish
- Lock-free union-find: Based on "Wait-free Parallel Algorithms for the Union-Find Problem" by Anderson and Woll