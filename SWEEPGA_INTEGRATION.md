# SweepGA Integration for SeqRush

## Overview

SeqRush now supports **sweepga** as an alternative alignment backend to allwave. This integration uses sweepga's FastGA wrapper for fast, efficient genome alignment.

## Architecture

### Feature Flags

The integration uses Cargo features to select the aligner at compile time:

```toml
[features]
default = ["use-allwave"]       # Default: use allwave
use-allwave = ["allwave", "lib_wfa2"]  # Allwave with WFA2
use-sweepga = ["sweepga"]       # SweepGA with FastGA
```

### Aligner Trait

All aligners implement the `Aligner` trait:

```rust
pub trait Aligner {
    fn align_sequences(
        &self,
        sequences: &[AlignmentSequence],
    ) -> Result<Vec<AlignmentRecord>, Box<dyn Error>>;
}
```

This provides a unified interface regardless of which backend is used.

### Implementations

1. **AllwaveAligner** (`src/aligner/allwave_impl.rs`)
   - Uses allwave's `AllPairIterator` for in-memory alignment
   - Fast for smaller datasets
   - Uses WFA2 algorithm

2. **SweepgaAligner** (`src/aligner/sweepga_impl.rs`)
   - Uses sweepga's `FastGAIntegration` for file-based alignment
   - Better for larger genomes
   - Uses FastGA with plane-sweep filtering
   - Supports k-mer frequency tuning

## Building

### With Allwave (Default)

```bash
cargo build --release
```

### With SweepGA

```bash
cargo build --release --no-default-features --features use-sweepga
```

## Usage

The aligner is transparent to the user - SeqRush automatically uses whichever backend was compiled in.

### Example with SweepGA

```bash
# Build with sweepga
cargo build --release --no-default-features --features use-sweepga

# Run normally
./target/release/seqrush -s input.fa -o output.gfa
```

### Optional: K-mer Frequency Tuning

When using sweepga, you can tune the k-mer frequency threshold for better performance on multi-genome datasets. This will require adding a CLI parameter:

```bash
./target/release/seqrush -s input.fa -o output.gfa --frequency 10
```

## Implementation Details

### SweepGA Workflow

1. **Write FASTA**: Sequences are written to a temporary FASTA file
2. **Create Index**: GDB and GIX index files are created via FastGA
3. **Align**: FastGA performs all-vs-all alignment
4. **Parse PAF**: Output is parsed into `AlignmentRecord` structures
5. **Graph Construction**: Records are fed into SeqRush's union-find algorithm

### File Management

- Temporary files are automatically cleaned up via Rust's `NamedTempFile`
- GDB/GIX index files are created alongside the temp FASTA
- PAF output is also in a temp file

### Performance Characteristics

**Allwave:**
- Pros: Fast in-memory alignment, no file I/O overhead
- Cons: Higher memory usage for large datasets
- Best for: Small to medium genomes (< 1 Gbp)

**SweepGA:**
- Pros: Efficient for large genomes, good filtering, k-mer frequency control
- Cons: File I/O overhead, index creation time
- Best for: Large genomes (> 1 Gbp), many genomes

## Testing

### Unit Tests

Both implementations share the same test suite to ensure consistency:

```bash
cargo test --no-default-features --features use-sweepga
```

### Integration Test

Compare outputs from both backends:

```bash
# With allwave
cargo build --release
./target/release/seqrush -s test.fa -o test_allwave.gfa

# With sweepga
cargo build --release --no-default-features --features use-sweepga
./target/release/seqrush -s test.fa -o test_sweepga.gfa

# Compare (should be similar, may have minor differences due to alignment heuristics)
diff test_allwave.gfa test_sweepga.gfa
```

## Future Enhancements

1. **Runtime Selection**: Allow choosing aligner at runtime via CLI flag
2. **Hybrid Mode**: Use allwave for small datasets, sweepga for large ones
3. **Parameter Tuning**: Expose more FastGA parameters (min identity, chain parameters, etc.)
4. **Parallel Alignment**: Leverage sweepga's all-pairs mode for multi-genome datasets
5. **Streaming Interface**: Process alignments as they're generated rather than batch processing

## Dependencies

### SweepGA Path

The integration expects sweepga to be at `../sweepga` relative to seqrush. Adjust the path in `Cargo.toml` if needed:

```toml
sweepga = { path = "/path/to/sweepga", optional = true }
```

### Build Requirements

- Rust 1.70+
- C compiler (for FastGA)
- Make
- Standard build tools

## Troubleshooting

### Build Times

FastGA includes C/C++ code that takes time to compile. First builds may take 5-10 minutes. Use `cargo build -j 1` to reduce memory usage during compilation.

### Missing Symbols

If you see linker errors about missing FastGA symbols, ensure sweepga builds correctly:

```bash
cd ../sweepga
cargo build --release
```

### Performance Issues

If sweepga is slower than expected, try adjusting the frequency parameter (lower = faster but less sensitive):

```bash
# Fast but may miss some alignments
--frequency 5

# Slower but more sensitive
--frequency 20
```

## References

- [SweepGA Repository](https://github.com/pangenome/sweepga)
- [FastGA Paper](https://doi.org/10.1101/2024.01.18.576254)
- [Allwave Repository](https://github.com/ekg/allwave)
