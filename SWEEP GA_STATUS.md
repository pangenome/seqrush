# SweepGA Integration Status

## ✅ Completed

1. **Architecture**: Created trait-based aligner abstraction
   - `src/aligner.rs` - Core trait and types
   - `src/aligner/allwave_impl.rs` - Allwave implementation
   - `src/aligner/sweepga_impl.rs` - SweepGA implementation

2. **Feature Flags**: Both aligners can be built together
   - Default features: `["use-allwave", "use-sweepga"]`
   - Runtime selection via `--aligner` CLI flag

3. **CLI Integration**: Added `--aligner` parameter
   - `--aligner allwave` (default)
   - `--aligner sweepga`

4. **Allwave Testing**: ✅ WORKS PERFECTLY
   ```bash
   ./target/release/seqrush -s test_small.fa -o test_allwave.gfa --aligner allwave -v
   ```
   - Successfully builds graph
   - Produces valid GFA output
   - 3 sequences → 4 compacted nodes

5. **Build System**: ✅ COMPILES SUCCESSFULLY
   - FastGA C++ code compiles (~16 seconds)
   - All dependencies resolve correctly
   - No warnings or errors

## ⚠️ Known Issues

### SweepGA FastGA Integration Issue

**Problem**: FastGA's `FAtoGDB` wrapper in sweepga fails when called from seqrush

**Symptoms**:
```
[FastGA] prepare_gdb: Converting "/tmp/.tmpRENTKJ" to GDB
[FastGA] Calling FAtoGDB: /home/erik/seqrush/target/release/build/fastga-rs-a03906b4c279d5f4/out/FAtoGDB /tmp/.tmpRENTKJ
Alignment failed: Failed to prepare GDB index
```

**Root Cause**: The fastga-rs Rust wrapper for FAtoGDB appears to have issues with:
- Temporary file paths
- Error propagation from C binaries
- Binary discovery/pathing

**Evidence**:
- `FAtoGDB` binary exists and is executable
- Running `FAtoGDB test_small.fa` directly **WORKS** and creates `.1gdb` and `.bps` files
- The issue is in the Rust wrapper layer, not the C binary itself

### Workarounds

#### Option 1: Use allwave (Recommended for now)
```bash
./target/release/seqrush -s input.fa -o output.gfa --aligner allwave
```

#### Option 2: Use external PAF with seqrush
If you want to use sweepga's filtering, generate PAF externally:
```bash
# Generate alignments with sweepga CLI
sweepga input.fa -o alignments.paf

# Build graph with seqrush
./target/release/seqrush -s input.fa -p alignments.paf -o output.gfa
```

#### Option 3: Fix fastga-rs (Future work)
The issue is in `sweepga/src/fastga_integration.rs` around lines 162-186 where `FastGAIntegration::prepare_gdb()` calls the orchestrator. Potential fixes:
1. Better error handling/logging from C binary
2. Fix temp file path handling
3. Ensure binaries are in PATH or use absolute paths
4. Check stdout/stderr from FAtoGDB for actual error messages

## Files Modified

### In seqrush:
- `Cargo.toml` - Added features and dependencies
- `src/lib.rs` - Added aligner module
- `src/aligner.rs` - Created trait abstraction (NEW)
- `src/aligner/allwave_impl.rs` - Allwave wrapper (NEW)
- `src/aligner/sweepga_impl.rs` - SweepGA wrapper (NEW)
- `src/seqrush.rs` - Added CLI flag and conditional compilation

### In sweepga:
- `src/lib.rs` - Exported `fastga_integration` module

## Testing Results

### Allwave (✅ Working)
```
Loaded 3 sequences
Building graph with 3 sequences (total length: 60)
Total sequence pairs: 9 (sparsification: None)
Built bidirected graph: 21 nodes, 21 edges, 3 paths
BidirectedGraph written to test_allwave.gfa: 4 nodes, 4 edges, 3 paths
```

### SweepGA (⚠️ Blocked by FastGA issue)
```
Using sweepga aligner
Running alignment with 3 sequences...
[sweepga] Aligning 3 sequences
[sweepga] Wrote 3 sequences to /tmp/.tmpRENTKJ
[sweepga] Creating GDB index...
Alignment failed: Failed to prepare GDB index
```

## Next Steps

1. **Immediate**: Use allwave for production work
2. **Short-term**: Debug fastga-rs integration
   - Add better error logging
   - Check FastGA binary stdout/stderr
   - Test with non-temp file paths
3. **Long-term**: Consider alternative FastGA bindings or direct FFI

## Branch Status

Branch: `sweepga`
- Ready to merge for allwave improvements
- SweepGA integration needs fastga-rs debugging
- All code compiles and passes type checking
- Architecture is sound and extensible

## Usage

```bash
# Default (allwave)
cargo build --release
./target/release/seqrush -s input.fa -o output.gfa

# Explicit allwave
./target/release/seqrush -s input.fa -o output.gfa --aligner allwave

# Sweepga (currently blocked)
./target/release/seqrush -s input.fa -o output.gfa --aligner sweepga
```
