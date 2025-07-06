# SeqRush

![CI](https://github.com/KristopherKubicki/seqrush/actions/workflows/ci.yml/badge.svg)
[![codecov](https://codecov.io/gh/KristopherKubicki/seqrush/branch/main/graph/badge.svg)](https://codecov.io/gh/KristopherKubicki/seqrush)

SeqRush is a prototype pangenome graph construction tool inspired by seqwish. Planned features include lock-free union-find structures and WFA2-based alignments.

## Overview

SeqRush builds pangenome graphs by:
1. *Planned*: perform all-vs-all pairwise alignments using WFA2 (Wavefront Alignment)
2. *Planned*: use a lock-free union-find data structure to merge matching positions
3. Construct a graph where sequences are embedded as paths

The design aims to leverage UFRush (lock-free union-find) for true parallel graph construction once implemented.

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

# Build the project with CLI support
cargo build --release --features cli
# Build the pangenome graph using flags
seqrush -s test.fasta -o test.gfa

# View the output
cat test.gfa
```

If you build without `--features cli`, use positional arguments instead:

```bash
seqrush test.fasta test.gfa
```

## Features

- *Planned*: **Lock-free Parallel Processing** via UFRush
- *Planned*: **Memory Efficient** alignment using WFA2's UltraLow mode
- *Planned*: **All-vs-All Alignment** of all input sequences
- **Configurable Parameters**: Alignment scoring, minimum match length
- **Standard GFA Output**: Compatible with tools like `odgi` and `vg`
- **Path Integrity**: Sequences are perfectly reconstructible from the graph

## Installation

### Prerequisites

- Rust via `rustup` (install the version pinned in `rust-toolchain.toml` or
  provide it via an offline setup)
- Git

### Build from Source

```bash
git clone https://github.com/KristopherKubicki/seqrush.git
cd seqrush
# rustup will automatically install the toolchain defined in `rust-toolchain.toml`
cargo build --release
```

To use the CLI flags, build the binary with the `cli` feature enabled:

```bash
cargo build --release --features cli
```

The binary will be available at `target/release/seqrush`.

## Usage

### Basic Usage

```bash
seqrush -s sequences.fasta -o graph.gfa
```

Enable the optional `cli` feature to use command-line flags.

### Advanced Options

```bash
seqrush \
  -s sequences.fasta \        # Input FASTA file
  -o graph.gfa \              # Output GFA file
  -t 8 \                      # Number of threads (default: 1)
  -k 15 \                     # Minimum match length (default: 15)
  -S "0,5,8,2,24,1" \         # Alignment scores: match,mismatch,gap1_open,gap1_extend,gap2_open,gap2_extend
  -v                          # Verbose output
```

#### Alignment Scoring

The `-S/--scores` parameter accepts comma-separated values:

```bash
# Two-piece affine gap model (default)
-S "0,5,8,2,24,1"  # match=0, mismatch=5, gap1_open=8, gap1_extend=2, gap2_open=24, gap2_extend=1

# Single affine gap model
-S "0,5,8,2"  # match=0, mismatch=5, gap_open=8, gap_extend=2

# Custom scoring for high similarity sequences
-S "0,4,6,1"  # More permissive scoring
```

Note: Two-piece affine gap support requires compatible WFA2 library version.

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
2. *Planned*: initialize a UFRush instance with one element per base
3. *Planned*: align all sequence pairs using WFA2
4. *Planned*: process matches ≥ `min_match_length` and unite positions
5. **Build Graph**: Walk sequences to identify nodes and edges

### Key Implementation Details

- *Planned*: **CIGAR Processing** for WFA2's fine-grained output
- *Planned*: **Match Accumulation** across CIGAR operations
- *Planned*: **Base Verification** since 'M' may represent mismatch
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
cargo test --features cli

# Run with verbose output
cargo test -- --nocapture

# Run a specific test
cargo test run_seqrush_writes_output
```

### Building Documentation

```bash
cargo doc --open
```

### Project Structure

![Project Structure](docs/project_structure.svg)

```
seqrush/
├── src/
│   ├── lib.rs          # Library interface
│   └── main.rs         # CLI binary
├── tests/
│   └── integration_tests.rs
├── Cargo.toml
└── README.md
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
Kristopher Kubicki, 2025
https://github.com/KristopherKubicki/seqrush
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgments

SeqRush is inspired by:
- [seqwish](https://github.com/ekg/seqwish) by Erik Garrison
- [WFA2-lib](https://github.com/smarco/WFA2-lib) by Santiago Marco-Sola
- [UFRush](https://crates.io/crates/uf_rush) lock-free union-find implementation

## Future Work

- Integrate WFA2 for efficient pairwise alignments
- Implement lock-free union-find (UFRush) for parallel graph construction
- Support streaming graph output for large datasets
- Extend sequence alphabet beyond DNA
