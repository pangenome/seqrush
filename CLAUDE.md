# SeqRush Development Guide

## Current Focus: Iterative Alignment Feature

### Overview
Implementing iterative graph construction with stabilization detection to reduce computational work by processing alignment pairs in priority order and stopping early when the graph converges.

### Status (2025-10-18)
- ✅ `--iterative` CLI flag added to Args struct (src/seqrush.rs:149)
- ✅ Flag recognized in `--help` output
- ✅ Build integration added (src/seqrush.rs:346-353)
- ❌ **NOT YET FUNCTIONAL** - shows error message directing to implementation
- ❌ Implementation needs to be ported from seqrush_clean.rs

### Implementation Plan

**Core Algorithm** (from seqrush_clean.rs:218-359):
1. Use AllWave TreeSampling to get prioritized pair list (k-nearest + k-farthest + random)
2. Process pairs iteratively in priority order
3. After every 10 pairs, count union-find components
4. Stop early if component count doesn't change for N consecutive checks (N=10)

**Technical Challenge:**
- seqrush_clean.rs uses simple `UFRush` (lock-free union-find over usize positions)
- seqrush.rs uses `BidirectedUnionFind` (handles oriented positions with Pos type)
- Different APIs, coordinate systems, and RC handling

**Required Work:**
1. Port `count_components()` to work with BidirectedUnionFind's Pos type
2. Adapt `align_and_unite_iterative()` to use seqrush.rs's alignment infrastructure
3. Test on simple sequences (/tmp/test_iter.fa)
4. Benchmark: does it reduce work while producing equivalent graphs?

### Testing Command
```bash
# Currently shows error message:
./target/release/seqrush -s /tmp/test_iter.fa -o /tmp/test.gfa -k 0 --iterative

# Expected behavior after implementation:
# - Align pairs in priority order
# - Print progress showing component count evolution
# - Stop early when graph stabilizes
# - Produce same graph as non-iterative mode
```

---

## Production Status: Graph Construction

### Current Quality (2025-10-17)
✅ **Structural Validation**: 100% success (28/28 HLA-Zoo graphs pass ODGI validation with ZERO changes)
- Bidirected graph implementation is correct
- Full Ygs pipeline (SGD + groom + topo sort) working
- Edge orientations and path handling verified

⚠️ **Layout Quality Gap**: RMSE 3.2x worse than ODGI
```
ODGI (full Ygs):        RMSE = 24.86 bp  ← Target
SeqRush (full Ygs):     RMSE = 83.23 bp  ← Current (3.2x worse)
```

**Root cause**: SGD implementation differences (parameter calculation, sampling strategy, or early termination)

### Command Reference
```bash
# Build graph with full Ygs pipeline (default)
./target/release/seqrush -s input.fa -o output.gfa -k 0

# Skip specific pipeline stages
./target/release/seqrush -s input.fa -o output.gfa -k 0 --skip-sgd
./target/release/seqrush -s input.fa -o output.gfa -k 0 --skip-groom
./target/release/seqrush -s input.fa -o output.gfa -k 0 --skip-topo

# Validate against ODGI
./test_hla_simple.sh  # Tests all HLA-Zoo graphs
```

---

## Architecture Notes

### Key Modules
- **src/seqrush.rs** - Production implementation (uses BidirectedUnionFind, full pipeline)
- **src/seqrush_clean.rs** - Simple implementation (uses UFRush, reference for iterative mode)
- **src/bidirected_union_find.rs** - Union-find for oriented positions (Pos type)
- **src/bidirected_builder.rs** - Builds bidirected graphs from union-find results
- **src/aligner/allwave_impl.rs** - AllWave aligner wrapper with prioritized pairs

### Graph Construction Pipeline
```
Load FASTA
  ↓
Align all-vs-all (AllWave/SweepGA)
  ↓
Unite matching positions (BidirectedUnionFind)
  ↓
Build bidirected graph from union components
  ↓
Compact linear chains (if --no-compact not set)
  ↓
Apply Ygs sorting:
  - Y: Path-guided SGD layout
  - g: Grooming (flip orientations for consistency)
  - s: Topological sort
  ↓
Write GFA
```

### Important Implementation Details

**Bidirected Graph (CRITICAL)**:
- Handles forward AND reverse complement alignments
- Each node has two orientations (forward/reverse)
- Edges connect oriented nodes: `5+ → 6-` means "forward strand of node 5 connects to reverse strand of node 6"
- Paths include orientations: `1+,2+,3-,4+` means "traverse nodes 1,2 forward, 3 reverse, 4 forward"
- Without bidirected support, RC alignments create invalid graphs

**Union-Find**:
- Positions encoded as `Pos` type (offset + orientation bit)
- LSB encodes orientation: even=forward, odd=reverse
- `unite_matching_region()` handles RC coordinate transformation
- BidirectedUnionFind ensures forward/reverse of same position are united

**Alignment Processing**:
- CIGAR parsing: only unite exact matches ('M' ops checked base-by-base)
- Minimum match length filtering (default k=0)
- Self-alignments MUST be included (establish sequence backbone)
- When PAF strand='-', the QUERY was reverse complemented (not target)

---

## Historical Context (Collapsed Issues - RESOLVED)

<details>
<summary>Click to expand resolved issues</summary>

### Issue 1: Missing Self-Alignments (FIXED)
**Problem**: Excluding self-alignments caused graph fragmentation
**Solution**: Changed `exclude_self: true` → `exclude_self: false`
**Result**: Positions within sequences now properly united

### Issue 2: PAF Strand Misinterpretation (FIXED)
**Problem**: Incorrectly assumed target was RC'd when strand='-'
**Solution**: Recognize that strand='-' means QUERY was RC'd
**Result**: RC alignments now processed correctly

### Issue 3: Grooming Instability (FIXED)
**Problem**: DFS with orientation flipping caused non-deterministic results
**Solution**: Exact BFS implementation matching ODGI
**Result**: 99.97% grooming success rate

### Issue 4: Vec Migration Bugs (FIXED)
**Problem**: HashMap→Vec migration introduced phantom node 0
**Solution**: Fixed all HashMap patterns, used filter_map correctly
**Result**: All tests passing in CI

</details>

---

## Quick Start for Development

```bash
# Test with simple sequences
cat > /tmp/test.fa << 'EOF'
>seq1
ACGTACGTACGT
>seq2
ACGTACGTACGT
>seq3
ACGTACGTACGT
EOF

./target/release/seqrush -s /tmp/test.fa -o /tmp/test.gfa -k 0 -v

# Validate output
odgi stats -i /tmp/test.gfa -S

# Run full HLA-Zoo validation
./test_hla_simple.sh
```
