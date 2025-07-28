# Reverse Complement Alignment Issue in SeqRush

## Problem Description

When aligning sequences that have reverse complement relationships (common in real genomic data), SeqRush correctly detects and aligns them in the proper orientation but then forces all sequences into forward orientation in the output graph. This creates a messy, incorrect graph structure.

## Example

Using HLA sequences from DRB1-3123.fa:
```bash
cargo run --release --bin seqrush -- -s HLA-zoo/seqs/DRB1-3123.fa -o d.gfa -k 0
odgi viz -z -i d.gfa -o d.png
```

Running wfmash shows several sequences align better in reverse complement:
```bash
wfmash HLA-zoo/seqs/DRB1-3123.fa -t 4 -k 19 > alignments.paf
grep -E "\s-\s" alignments.paf  # Shows reverse strand alignments
```

## Root Cause

1. **Alignment Phase**: SeqRush correctly:
   - Tests both forward-forward and forward-reverse orientations
   - Chooses the better scoring orientation
   - Maps reverse complement positions back to forward strand coordinates

2. **Graph Construction Phase**: The problem occurs here:
   - All orientation information is lost
   - The union-find structure only tracks position equivalences, not orientations
   - Graph edges are all written as forward-to-forward (`L\t{}\t+\t{}\t+\t0M`)
   - Path nodes are all written with forward orientation (`{}+`)

## Code Locations

- Orientation detection: `src/seqrush.rs:275-296`
- Position mapping for RC: `src/seqrush.rs:412-417`
- Graph output (forced forward): `src/seqrush.rs:685,692-693`

## Impact

Sequences that should be reverse complemented are forced into forward orientation, creating:
- Incorrect graph topology
- Many unnecessary nodes and edges
- Poor visualization with many crossing paths
- Biologically incorrect representation

## Potential Solutions

### Short-term Workaround
- Use `--test-mode` flag to disable reverse complement detection (not recommended for real data)

### Proper Fix
The codebase already includes bidirected graph support in:
- `src/bidirected_graph.rs`
- `src/bidirected_ops.rs`
- `src/seqrush_bidirected.rs`

The fix would require:
1. Using bidirected graph representation throughout
2. Tracking orientation in union-find operations
3. Outputting proper orientations in GFA format
4. Possibly switching to a different main implementation that supports bidirected graphs

## Comparison with Other Tools

This is noted as a "long standing bug in odgi's toposort, and vg's too", suggesting it's a common issue in pangenome graph construction when forcing unidirectional representation of inherently bidirectional data.