# SeqRush

![CI](https://github.com/KristopherKubicki/seqrush/actions/workflows/ci.yml/badge.svg)
[![codecov](https://codecov.io/gh/KristopherKubicki/seqrush/branch/main/graph/badge.svg)](https://codecov.io/gh/KristopherKubicki/seqrush)

A high-performance, parallel pangenome graph construction tool that implements a simplified seqwish-like algorithm using lock-free union-find data structures.

## Overview

SeqRush builds pangenome graphs by:
1. Performing all-vs-all pairwise alignments using WFA2 (Wavefront Alignment)
2. Using a lock-free union-find data structure to merge matching positions
3. Constructing a graph where sequences are embedded as paths

The key innovation is the use of UFRush (lock-free union-find) to enable true parallel graph construction without synchronization overhead.

## Quick Start

```bash
# Create a test FASTA file
cat > test.fasta << EOF
>seq1
ATCGATCGATCGATCG
>seq2
ATCGATGGATCGATCG
>seq3
ATCGATCGATCGATGG
EOF

# Build the pangenome graph
seqrush -s test.fasta -o test.gfa

# View the output
cat test.gfa
```

## Features

- **Lock-free Parallel Processing**: Uses UFRush for wait-free union-find operations
- **Memory Efficient**: O(s) memory alignment using WFA2's UltraLow mode
- **All-vs-All Alignment**: Comprehensive pairwise alignment of all input sequences
- **Configurable Parameters**: Alignment scoring, minimum match length
- **Standard GFA Output**: Compatible with tools like `odgi` and `vg`
- **Path Integrity**: Sequences are perfectly reconstructible from the graph

## Installation

### Prerequisites

- Rust 1.70 or later
- Git

### Build from Source

```bash
git clone https://github.com/KristopherKubicki/seqrush.git
cd seqrush
cargo build --release
```

The binary will be available at `target/release/seqrush`.

## Usage

### Basic Usage

```bash
seqrush -s sequences.fasta -o graph.gfa
```

### Advanced Options

```bash
seqrush \
  -s sequences.fasta \        # Input FASTA file
  -o graph.gfa \              # Output GFA file
  -t 8 \                      # Number of threads (default: 4)
  -k 15 \                     # Minimum match length (default: 15)
  --match-score 0 \           # WFA match score (default: 0)
  --mismatch-penalty 5 \      # WFA mismatch penalty (default: 5)
  --gap-open 8 \              # WFA gap open penalty (default: 8)
  --gap-extend 2 \            # WFA gap extend penalty (default: 2)
  -v                          # Verbose output
```

### Example Workflow

```bash
# Build a pangenome graph
seqrush -s genomes.fasta -o pangenome.gfa

# Visualize with odgi
odgi build -g pangenome.gfa -o pangenome.og
odgi viz -i pangenome.og -o pangenome.png

# Check graph statistics
odgi stats -i pangenome.og -S
```

## Algorithm Details

### Core Algorithm

1. **Load Sequences**: Read FASTA file and assign global positions to each base
2. **Initialize Union-Find**: Create a UFRush instance with one element per base
3. **Pairwise Alignment**: Align all sequence pairs using WFA2
4. **Process Matches**: For each alignment, unite positions with matches ≥ min_match_length
5. **Build Graph**: Walk sequences through union-find to identify nodes and edges

### Key Implementation Details

- **CIGAR Processing**: Handles WFA2's fine-grained CIGAR output (individual M operations)
- **Match Accumulation**: Tracks consecutive matches across multiple CIGAR operations
- **Base Verification**: Confirms actual base matches since 'M' can represent match or mismatch
- **Path Construction**: Each sequence becomes a path through deduplicated nodes

### Performance Characteristics

- Time Complexity: O(n²×L) for n sequences of length L (pairwise alignment)
- Space Complexity: O(N) where N is total sequence length
- Parallel Scaling: Near-linear with thread count for alignment phase

## Output Format

SeqRush generates GFA 1.0 format with:
- **H**: Header with version
- **S**: Segments (nodes) with single-character sequences
- **P**: Paths representing input sequences
- **L**: Links between adjacent nodes in paths

Example output:
```
H	VN:Z:1.0
S	1	A
S	2	C
S	3	G
P	seq1	1+,2+,3+	*
L	1	+	2	+	0M
L	2	+	3	+	0M
```

## Development

### Running Tests

```bash
# Run all tests
cargo test

# Run with verbose output
cargo test -- --nocapture

# Run specific test
cargo test test_performance_scaling
```

### Building Documentation

```bash
cargo doc --open
```

### Debugging Tools

SeqRush includes several debugging utilities:

```bash
# Debug union-find operations
cargo run --bin debug_union_find

# Debug alignment processing
cargo run --bin debug_alignment

# Debug CIGAR string parsing
cargo run --bin debug_cigar

# Test match run tracking
cargo run --bin test_match_runs
```

### Project Structure

```
seqrush/
├── src/
│   ├── lib.rs          # Library interface
│   ├── main.rs         # CLI binary
│   └── seqrush.rs      # Core implementation
├── tests/
│   └── integration_tests.rs
└── Cargo.toml
```

## Limitations

- Input sequences must fit in memory
- Currently builds the entire graph at once (no streaming)
- Single-character nodes (no compaction)
- Limited to DNA sequences (ACGT alphabet)

## Known Issues

- Path integrity verification may fail for some sequences with complex indel patterns. The graph structure is correct, but path reconstruction needs improvement for certain edge cases.

## Citation

If you use SeqRush in your research, please cite:

```
SeqRush: Lock-free parallel pangenome graph construction
[Your Name], 2024
https://github.com/KristopherKubicki/seqrush
```

## License

[License information here]

## Acknowledgments

SeqRush is inspired by:
- [seqwish](https://github.com/ekg/seqwish) by Erik Garrison
- [WFA2-lib](https://github.com/smarco/WFA2-lib) by Santiago Marco-Sola
- [UFRush](https://crates.io/crates/uf_rush) lock-free union-find implementation