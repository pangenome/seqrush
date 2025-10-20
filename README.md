# SeqRush

![CI](https://github.com/pangenome/seqrush/actions/workflows/ci.yml/badge.svg)

SeqRush is a fast pangenome graph construction tool that builds bidirected sequence graphs with integrated layout optimization.

## Overview

SeqRush constructs pangenome variation graphs by:
1. Performing all-vs-all pairwise sequence alignments using AllWave or SweepGA aligners
2. Using a lock-free bidirected union-find structure to merge matching positions
3. Building a bidirected graph with proper handling of forward and reverse complement alignments
4. Applying the Ygs pipeline (SGD layout + grooming + topological sort) for optimized graph layout
5. Outputting GFA format compatible with tools like ODGI and VG

## Features

- **Bidirected Graph Support**: Properly handles forward and reverse complement alignments
- **Multiple Alignment Backends**: AllWave (wavefront-based) and SweepGA (sweep-line genetic algorithm)
- **Sparsification Strategies**: TreeSampling for intelligent pair selection (k-nearest + k-farthest + random)
- **Iterative Alignment**: Optional two-phase alignment with early stopping for large datasets
- **Graph Layout Optimization**: Integrated Ygs pipeline (path-guided SGD + grooming + topological sort)
- **Node Compaction**: Linear chain compaction to reduce graph complexity
- **Standard GFA Output**: Compatible with downstream pangenome tools

## Quick Start

```bash
# Install Rust if needed
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Clone and build
git clone --recursive https://github.com/pangenome/seqrush.git
cd seqrush
cargo build --release --features use-allwave

# Build a graph from sequences
./target/release/seqrush -s input.fa -o output.gfa -k 0

# View statistics
odgi stats -i output.gfa -S
```

## Installation

### Prerequisites

- Rust 1.70+ (automatically managed via `rust-toolchain.toml`)
- Git with submodule support

### Build from Source

```bash
git clone --recursive https://github.com/pangenome/seqrush.git
cd seqrush
cargo build --release --features use-allwave
```

The binary will be available at `target/release/seqrush`.

## Usage

### Basic Usage

```bash
# Build graph with default AllWave aligner
./target/release/seqrush -s sequences.fa -o graph.gfa -k 0

# Use SweepGA aligner instead
./target/release/seqrush -s sequences.fa -o graph.gfa -k 0 --aligner sweepga

# Enable iterative alignment for large datasets
./target/release/seqrush -s sequences.fa -o graph.gfa -k 0 --iterative
```

### Common Options

```bash
seqrush \
  -s sequences.fa \          # Input FASTA file
  -o graph.gfa \             # Output GFA file
  -k 0 \                     # Minimum match length (0 = unite all exact matches)
  -t 4 \                     # Number of threads
  -v \                       # Verbose output
  --iterative \              # Use iterative alignment with early stopping
  --no-compact \             # Skip linear chain compaction
  --skip-sgd \               # Skip SGD layout step
  --skip-groom \             # Skip grooming step
  --skip-topo                # Skip topological sort step
```

### Alignment Scoring

```bash
# Default two-piece affine gap scoring
-S "0,5,8,2,24,1"  # match, mismatch, gap1_open, gap1_extend, gap2_open, gap2_extend

# Orientation refinement scoring
--orientation-scores "0,1,1,1"  # match, mismatch, gap_open, gap_extend
```

### Sparsification

Control which sequence pairs to align:

```bash
# Align all pairs (default)
--sparsification "1.0"

# TreeSampling: k-nearest + k-farthest + random fraction
# Format: "tree:k_nearest,k_farthest,random_fraction"
--sparsification "tree:3,3,0.1"
```

### Iterative Alignment

For large datasets, use iterative mode to reduce computational work:

```bash
./target/release/seqrush -s input.fa -o output.gfa -k 0 --iterative
```

Iterative mode processes alignments in two phases:
1. **Phase 1**: All tree pairs (k-nearest + k-farthest) - guarantees connectivity
2. **Phase 2**: Random pairs with early stopping when graph stabilizes

## Algorithm Details

### Graph Construction Pipeline

```
Load FASTA sequences
  ↓
Align sequence pairs (AllWave/SweepGA)
  ↓
Unite matching positions (BidirectedUnionFind)
  ↓
Build bidirected graph from union components
  ↓
Compact linear chains (optional)
  ↓
Apply Ygs sorting pipeline:
  - Y: Path-guided SGD layout
  - g: Grooming (orientation consistency)
  - s: Topological sort
  ↓
Write GFA
```

### Key Implementation Details

**Bidirected Graphs**: Sequences can align in forward or reverse complement orientation. Each node has two orientations (forward `+` and reverse `-`), and edges connect oriented nodes: `5+ → 6-` means "forward strand of node 5 connects to reverse strand of node 6".

**Union-Find**: Lock-free bidirected union-find efficiently merges matching positions while tracking orientation. Positions are encoded as `Pos` type with orientation bit.

**Alignment Processing**: Only exact matches from CIGAR strings are united. Self-alignments are included to establish the sequence backbone.

**Graph Layout**: The Ygs pipeline optimizes node ordering:
- **Y (SGD)**: Path-guided stochastic gradient descent for 1D layout
- **g (grooming)**: BFS-based orientation flipping for consistency
- **s (topological sort)**: Final ordering based on path structure

### Performance Characteristics

- Time Complexity: O(n²×L×A) for n sequences of length L with alignment cost A
- Space Complexity: O(N) where N is total sequence length
- Parallel Scaling: Near-linear for alignment phase

## Output Format

SeqRush generates GFA 1.0 format:

```
H	VN:Z:1.0
S	1	ACGT         # Segment (node) with sequence
S	2	GC
L	1	+	2	+	0M   # Link between oriented nodes
P	seq1	1+,2+	*    # Path representing input sequence
```

## Development

### Running Tests

```bash
# Run all tests
cargo test

# Run specific test suite
cargo test bidirected

# Run with output
cargo test -- --nocapture
```

### Project Structure

```
seqrush/
├── src/
│   ├── lib.rs                      # Library interface
│   ├── main.rs                     # CLI entry point
│   ├── seqrush.rs                  # Main graph construction logic
│   ├── bidirected_union_find.rs   # Lock-free oriented union-find
│   ├── bidirected_builder.rs      # Graph builder from union-find
│   ├── path_sgd.rs                 # Path-guided SGD layout
│   ├── grooming.rs                 # Orientation consistency
│   ├── topological_sort.rs         # Path-based topological ordering
│   └── aligner/                    # Alignment backends
│       ├── allwave_impl.rs
│       └── sweepga_impl.rs
├── tests/                          # Integration tests
├── allwave/                        # AllWave aligner (submodule)
├── WFA2-lib/                       # WFA2 library (submodule)
└── CLAUDE.md                       # Development guide
```

## Validation

SeqRush has been validated against ODGI on the HLA-Zoo dataset (28 HLA graphs):
- **Structural correctness**: 100% (all graphs pass ODGI validation with zero changes)
- **Layout quality**: RMSE currently 3.2x worse than ODGI (under investigation)

## Limitations

- Input sequences must fit in memory
- SGD layout quality needs improvement (see CLAUDE.md for details)
- Limited to DNA sequences (ACGT alphabet)

## Citation

If you use SeqRush in your research, please cite:

```
SeqRush: Fast bidirected pangenome graph construction
Erik Garrison, Kristopher Kubicki, 2025
https://github.com/pangenome/seqrush
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

SeqRush is inspired by and builds upon:
- [seqwish](https://github.com/ekg/seqwish) - Original pangenome graph construction algorithm
- [ODGI](https://github.com/pangenome/odgi) - Optimized dynamic genome/graph implementation
- [AllWave](https://github.com/urbanslug/allwave) - Wavefront-based sequence aligner
- [WFA2-lib](https://github.com/smarco/WFA2-lib) - Wavefront alignment algorithm library

Special thanks to Erik Garrison for the seqwish algorithm and graph construction concepts.
