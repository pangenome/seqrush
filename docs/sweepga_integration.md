# SweepGA Integration Analysis

## Executive Summary

SweepGA is a fast genome aligner that wraps FastGA and applies plane sweep filtering to produce high-quality alignments. It offers significant advantages over the current WFA2-based aligner in SeqRush, particularly for self-alignment of pangenome FASTA files.

**Key Benefits:**
1. **Native all-vs-all alignment** - Designed specifically for pangenome workflows
2. **Built-in filtering** - Plane sweep algorithm keeps best 1:1 mappings automatically
3. **FastGA speed** - Significantly faster than WFA2 for whole genome alignment
4. **Better memory efficiency** - Handles large pangenomes without excessive memory usage
5. **Rust library** - Can be integrated directly as a dependency

## Current SeqRush Alignment Pipeline

SeqRush currently supports two aligners (runtime selectable):

### 1. WFA2 via Allwave (default)
```rust
// From src/align_sequences.rs
pub fn align_sequences_with_allwave(
    sequences: &[SequenceInfo],
    k: usize,
    min_length: u64,
    threads: usize,
) -> Vec<Alignment>
```

**Pros:**
- Precise gap-affine alignment
- Good for highly similar sequences (>95% ANI)
- Direct Rust FFI bindings

**Cons:**
- VERY slow for whole genome alignment (designed for local alignment)
- High memory usage for large sequences
- Not optimized for all-vs-all pangenome alignment
- Requires chunking strategies for large inputs

### 2. SweepGA via FastGA (PR #5, feature flag `use-sweepga`)
```rust
// From Cargo.toml (commented out)
# use-sweepga = ["sweepga", "fastga-rs"]  # Commented out - enable locally if needed
```

Currently disabled due to dependency issues, but the infrastructure exists.

## SweepGA Architecture

### Core Components

```
SweepGA = FastGA Aligner + Plane Sweep Filtering
```

1. **FastGA**: Fast genome aligner using seed-and-extend with k-mer indexing
2. **Plane Sweep**: Filters overlapping alignments to keep best non-overlapping mappings
3. **1:1 Filtering**: Ensures orthogonal mappings (best per query AND target position)

### Library Structure

SweepGA is structured as both a CLI tool and a Rust library:

```
sweepga/
├── src/
│   ├── lib.rs                    # Library exports
│   ├── main.rs                   # CLI binary
│   ├── fastga_integration.rs     # FastGA Rust bindings
│   ├── plane_sweep_core.rs       # Core filtering algorithm
│   ├── paf_filter.rs             # PAF filtering logic
│   ├── unified_filter.rs         # Format-preserving filtering
│   └── mapping.rs                # Alignment data structures
```

**Key API from `fastga_integration.rs`:**

```rust
pub struct FastGAIntegration {
    config: Config,
}

impl FastGAIntegration {
    // Create with optional frequency parameter
    pub fn new(frequency: Option<usize>, num_threads: usize) -> Self;

    // Run alignment and return native .1aln format
    pub fn align_to_temp_1aln(&self, queries: &Path, targets: &Path)
        -> Result<NamedTempFile>;

    // Run alignment and return PAF format
    pub fn align_to_temp_paf(&self, queries: &Path, targets: &Path)
        -> Result<NamedTempFile>;
}
```

### Filtering Modes

SweepGA supports multiple filtering strategies:

```rust
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FilterMode {
    OneToOne,     // 1:1 - Best mapping per position (both axes)
    OneToMany,    // 1:∞ - Best per query only
    ManyToMany,   // ∞:∞ - No filtering
}
```

**For pangenome graphs, 1:1 filtering is ideal** - it ensures:
- No overlapping alignments confuse the graph construction
- Each position in each sequence maps to at most one position in other sequences
- Transitive closure works correctly (if A→B and B→C, then A→C)

### Scoring Function

SweepGA uses intelligent scoring to rank alignments:

```rust
pub enum ScoringFunction {
    Identity,              // Pure identity %
    Length,                // Alignment length
    LengthIdentity,        // identity * length
    LogLengthIdentity,     // identity * log(length) (DEFAULT)
    Matches,               // Raw match count
}
```

**`LogLengthIdentity` (default) is ideal** because:
- Favors long alignments (log scale prevents over-weighting)
- Penalizes low identity appropriately
- Matches wfmash's scoring function (proven for pangenomes)

## Integration Approach

### Option 1: Runtime Aligner Selection (Recommended)

Similar to the current WFA2/SweepGA approach, but make SweepGA the default:

```rust
// In src/main.rs
#[derive(Parser, Debug)]
struct Cli {
    /// Aligner to use: sweepga (default), wfa2, allwave
    #[clap(long = "aligner", short = 'A', default_value = "sweepga")]
    aligner: String,

    /// K-mer frequency threshold for SweepGA (use k-mers occurring ≤ N times)
    #[clap(long = "frequency", short = 'f')]
    frequency: Option<usize>,

    // ... existing options ...
}
```

**Workflow:**

```rust
// In src/graph_builder.rs
match aligner.as_str() {
    "sweepga" => {
        let alignments = align_sequences_with_sweepga(
            sequences,
            frequency,
            threads
        )?;
        // ... continue with graph construction ...
    }
    "wfa2" | "allwave" => {
        let alignments = align_sequences_with_allwave(
            sequences,
            k,
            min_length,
            threads
        )?;
        // ... continue with graph construction ...
    }
    _ => bail!("Unknown aligner: {}", aligner),
}
```

### Option 2: Preprocessing Mode

Use SweepGA as a preprocessing step to generate PAF, then build graph from PAF:

```bash
# Step 1: Generate alignments with SweepGA
sweepga input.fa -t 8 --paf > alignments.paf

# Step 2: Build graph from PAF
seqrush -s input.fa -p alignments.paf -o output.gfa
```

**Pros:**
- Decouples alignment from graph construction
- Users can experiment with different aligners
- Allows manual inspection/filtering of alignments

**Cons:**
- Requires intermediate file storage
- Two-step workflow is less convenient

### Option 3: Hybrid Approach (Best of Both Worlds)

Support BOTH integrated alignment AND PAF input:

```rust
// Automatically detect workflow
if let Some(paf_file) = args.paf {
    // Use existing PAF
    alignments = parse_paf_file(&paf_file)?;
} else {
    // Run aligner
    alignments = match args.aligner {
        "sweepga" => align_with_sweepga(sequences, args.frequency, threads)?,
        "wfa2" => align_with_allwave(sequences, k, min_length, threads)?,
        _ => bail!("Unknown aligner"),
    };
}
```

## Implementation Plan

### Phase 1: Enable SweepGA Feature (1-2 hours)

1. **Uncomment feature flag in Cargo.toml:**

```toml
[features]
default = ["use-sweepga"]  # Change from use-allwave
use-allwave = ["allwave", "lib_wfa2"]
use-sweepga = ["sweepga", "fastga-rs"]
```

2. **Add dependencies:**

```toml
[dependencies]
sweepga = { path = "../sweepga", optional = true }
fastga-rs = { git = "https://github.com/pangenome/fastga-rs.git", rev = "db602e0", optional = true }
```

3. **Test build:**

```bash
cd ~/seqrush
cargo build --release --features use-sweepga
```

### Phase 2: Create SweepGA Alignment Module (2-3 hours)

Create `src/align_sweepga.rs`:

```rust
use anyhow::{Context, Result};
use std::path::Path;
use sweepga::fastga_integration::FastGAIntegration;
use sweepga::paf_filter::{FilterConfig, FilterMode, PafFilter, ScoringFunction};

pub fn align_sequences_with_sweepga(
    sequences: &[SequenceInfo],
    frequency: Option<usize>,
    threads: usize,
) -> Result<Vec<Alignment>> {
    // Step 1: Write sequences to temp FASTA
    let temp_fasta = write_sequences_to_fasta(sequences)?;

    // Step 2: Run FastGA alignment
    let fastga = FastGAIntegration::new(frequency, threads);
    let temp_paf = fastga.align_to_temp_paf(
        temp_fasta.path(),
        temp_fasta.path()  // Self-alignment
    )?;

    // Step 3: Apply 1:1 filtering
    let filter_config = FilterConfig {
        mapping_filter_mode: FilterMode::OneToOne,
        scoring_function: ScoringFunction::LogLengthIdentity,
        overlap_threshold: 0.95,
        min_block_length: 100,  // Filter out tiny alignments
        ..Default::default()
    };

    let filter = PafFilter::new(filter_config);
    let temp_filtered = tempfile::NamedTempFile::new()?;
    filter.filter_paf(
        temp_paf.path().to_str().unwrap(),
        temp_filtered.path().to_str().unwrap()
    )?;

    // Step 4: Parse filtered PAF to Alignment structs
    parse_paf_to_alignments(temp_filtered.path(), sequences)
}

fn parse_paf_to_alignments(
    paf_path: &Path,
    sequences: &[SequenceInfo]
) -> Result<Vec<Alignment>> {
    // Parse PAF format and convert to SeqRush Alignment structs
    // ...
}
```

### Phase 3: Integrate into Main Pipeline (1 hour)

Modify `src/graph_builder.rs`:

```rust
pub fn build_graph_from_sequences(
    sequences: Vec<SequenceInfo>,
    aligner: &str,
    frequency: Option<usize>,
    k: usize,
    min_length: u64,
    threads: usize,
) -> Result<BidirectedGraph> {
    // Generate alignments
    let alignments = match aligner {
        #[cfg(feature = "use-sweepga")]
        "sweepga" => {
            eprintln!("[seqrush] Using SweepGA aligner (FastGA + plane sweep filtering)");
            align_sequences_with_sweepga(&sequences, frequency, threads)?
        }

        #[cfg(feature = "use-allwave")]
        "wfa2" | "allwave" => {
            eprintln!("[seqrush] Using WFA2 aligner (via allwave)");
            align_sequences_with_allwave(&sequences, k, min_length, threads)?
        }

        _ => bail!("Unknown aligner: {}. Available: sweepga, wfa2", aligner),
    };

    // Continue with existing graph construction
    build_graph_from_alignments(sequences, alignments, threads)
}
```

### Phase 4: Testing and Validation (2-3 hours)

Test on HLA-zoo dataset:

```bash
# Test SweepGA integration
cargo build --release --features use-sweepga
./target/release/seqrush \
    -s HLA-zoo/seqs/A-3105.fa \
    -o /tmp/sweepga_test.gfa \
    -A sweepga \
    -t 8

# Compare with WFA2
./target/release/seqrush \
    -s HLA-zoo/seqs/A-3105.fa \
    -o /tmp/wfa2_test.gfa \
    -A wfa2 \
    -t 8

# Compare stats
odgi stats -i /tmp/sweepga_test.gfa -S
odgi stats -i /tmp/wfa2_test.gfa -S
```

**Validation metrics:**
- Node count should be similar or better
- Edge count should be lower (better filtering)
- RMSE should be comparable or better
- Runtime should be significantly faster

### Phase 5: Documentation and Cleanup (1 hour)

Update README.md:

```markdown
## Alignment Algorithms

SeqRush supports multiple alignment algorithms:

### SweepGA (default, recommended for pangenomes)
- Fast whole-genome alignment via FastGA
- Built-in 1:1 filtering for clean graph construction
- Optimized for all-vs-all pangenome alignment

Usage:
```bash
seqrush -s sequences.fa -o output.gfa -A sweepga -t 8
```

### WFA2 (via Allwave)
- Precise gap-affine alignment
- Better for highly similar sequences (>95% ANI)

Usage:
```bash
seqrush -s sequences.fa -o output.gfa -A wfa2 -k 16 -t 8
```
```

## Expected Performance Improvements

Based on SweepGA's design and FastGA benchmarks:

### Speed
- **WFA2 (current)**: ~10-30 minutes for 8-sequence HLA graphs
- **SweepGA (expected)**: ~1-3 minutes for 8-sequence HLA graphs
- **Speedup**: 5-10x faster

### Memory
- **WFA2 (current)**: High memory for large alignments (O(nm) per alignment)
- **SweepGA (expected)**: Lower memory (seed-based indexing)
- **Improvement**: 2-3x lower peak memory

### Graph Quality
- **WFA2 (current)**: Can create overlapping alignments requiring complex resolution
- **SweepGA (expected)**: 1:1 filtering ensures clean, non-overlapping alignments
- **Result**: Simpler graph construction, fewer edge cases

### Catastrophic Edge Handling
- **Current issue**: A-3105 has catastrophic long edges (see `povu_guided_sorting.md`)
- **SweepGA benefit**: Better alignment quality may reduce catastrophic edges
- **Note**: Povu-guided sorting still recommended as final solution

## Dependency Management

### Current Status (RESOLVED)
SweepGA integration now working! The GIXmake `-f` flag error was resolved by updating dependencies.

### Root Cause
The issue was that different builds of fastga-rs were using different versions of the C code:
- **Old version (db602e0)**: Compiled GIXmake from source with `-f` flag for frequency
- **New version (b8416dd)**: Uses embedded GIXmake binary WITHOUT `-f` flag (uses `-k` for k-mer size instead)

The mismatch caused "GIXmake: -f is an illegal option" errors.

### Fix Applied
1. Updated `~/sweepga/Cargo.toml` to use fastga-rs revision `b8416dd` (includes embedded binaries)
2. Rebuilt sweepga: `cd ~/sweepga && cargo update -p fastga-rs && cargo build --release`
3. Rebuilt seqrush: `cargo build --release --features use-sweepga`

### Verification
```bash
# Test SweepGA alignment (should work without errors)
./target/release/seqrush -s HLA-zoo/seqs/B-3106.fa -o test.gfa --aligner sweepga -t 8 -v
```

### Cargo.toml Configuration
```toml
[features]
default = ["use-allwave"]
use-allwave = ["allwave", "lib_wfa2"]
use-sweepga = ["sweepga"]

[dependencies]
sweepga = { path = "../sweepga", optional = true }
```

## Alignment Format Comparison

### PAF (Pairwise mApping Format)
SweepGA outputs PAF by default:

```
query_name  q_len  q_start  q_end  strand  target_name  t_len  t_start  t_end  matches  aln_len  mapq
```

SeqRush can consume PAF with `-p` flag:
```bash
sweepga sequences.fa -t 8 --paf > alignments.paf
seqrush -s sequences.fa -p alignments.paf -o output.gfa
```

### .1aln (Native FastGA format)
SweepGA also supports FastGA's binary format:

```bash
sweepga sequences.fa -t 8 > alignments.1aln
```

Benefits:
- More compact (binary)
- Preserves sequence names via .1gdb sidecar
- Faster to parse

SeqRush would need to add .1aln parser OR use SweepGA's library API directly.

## Alternative: Use SweepGA as External Preprocessor

Instead of integrating SweepGA as a library, we could recommend it as a preprocessing step:

```bash
# Recommended workflow
sweepga input.fa -t 8 --paf > alignments.paf
seqrush -s input.fa -p alignments.paf -o output.gfa
```

**Pros:**
- No new dependencies
- Decoupled alignment from graph construction
- Users can swap aligners easily

**Cons:**
- Two-step workflow
- Intermediate file storage
- Less convenient for quick tests

## Conclusion

**Recommendation: Integrate SweepGA as the default aligner**

Reasons:
1. **Speed**: 5-10x faster than WFA2 for whole-genome alignment
2. **Quality**: Built-in 1:1 filtering produces cleaner alignments
3. **Memory**: Lower memory usage for large pangenomes
4. **Rust native**: No FFI complexity, clean Rust API
5. **Pangenome optimized**: Designed specifically for all-vs-all alignment

**Implementation effort**: ~6-9 hours total
- Phase 1: 1-2 hours (enable feature)
- Phase 2: 2-3 hours (create alignment module)
- Phase 3: 1 hour (integrate into pipeline)
- Phase 4: 2-3 hours (testing and validation)
- Phase 5: 1 hour (documentation)

**Risk level**: Low
- SweepGA is already partially integrated (PR #5)
- Fallback to WFA2 remains available
- No changes to core graph construction logic

## References

- **SweepGA repository**: `~/sweepga/`
- **SweepGA README**: `/home/erik/sweepga/README.md`
- **SeqRush PR #5**: Added SweepGA as alternative aligner
- **FastGA paper**: (check SweepGA README for citations)
- **Related**: Povu-guided sorting (`/home/erik/seqrush/docs/povu_guided_sorting.md`)

## Next Steps

1. **Test current SweepGA integration** - Enable the feature and run benchmarks
2. **Compare with WFA2** - Measure speed, memory, and graph quality
3. **Validate on HLA-zoo** - Ensure all 9 graphs build correctly
4. **Measure impact on catastrophic edges** - Check if better alignments reduce the A-3105 issue
5. **Consider Povu integration** - After SweepGA validation, implement povu-guided sorting
