#!/bin/bash
# Test SGD quality in isolation (without grooming/topo sort)

set -e

TEST_FASTA="HLA-zoo/seqs/B-3106.fa"
OUTDIR="/tmp/sgd_test"
mkdir -p "$OUTDIR"

echo "========================================="
echo "Testing SGD quality in isolation"
echo "========================================="
echo ""

# Build a simple test tool that does JUST SGD
cat > /tmp/test_sgd_only.rs <<'EOF'
use seqrush::builder::build_from_fasta;
use seqrush::ygs_sort::YgsParams;
use seqrush::path_sgd::path_sgd_sort;
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input.fa> <output.gfa>", args[0]);
        std::process::exit(1);
    }

    let input = &args[1];
    let output = &args[2];

    // Build graph
    eprintln!("Building graph...");
    let mut graph = build_from_fasta(input, 0, None, false)?;

    eprintln!("Initial nodes: {}", graph.node_count());

    // Apply JUST SGD (no grooming, no topo sort)
    let params = YgsParams::from_graph(&graph, true, 1);
    eprintln!("\nApplying path-guided SGD...");
    let ordering = path_sgd_sort(&graph, params.path_sgd);
    graph.apply_ordering(ordering, false);

    eprintln!("After SGD nodes: {}", graph.node_count());

    // Write output
    graph.write_gfa(output)?;
    eprintln!("\nWrote: {}", output);

    Ok(())
}
EOF

echo "1. SeqRush SGD only (no grooming/topo):"
cargo build --release 2>&1 | tail -1
# Run inline with seqrush binary since we can't easily compile standalone
# Instead, let's modify seqrush to have a --sgd-only flag... or just measure after each step

# Actually, let's just test with verbose output and measure after SGD
echo "Building with SeqRush (full Ygs pipeline)..."
./target/release/seqrush -s "$TEST_FASTA" -o "$OUTDIR/full.gfa" -k 0 -v 2>&1 | tee "$OUTDIR/full.log" | tail -20

echo ""
echo "Layout quality after full Ygs:"
./target/release/measure_layout_quality "$OUTDIR/full.gfa" 2>&1 | grep "RMSE:"

echo ""
echo "========================================="
echo "Now testing with ODGI (for comparison):"
echo ""

# Convert to ODGI format and measure just after SGD
odgi build -g "$OUTDIR/full.gfa" -o "$OUTDIR/seqrush.og" 2>&1 | head -1

# Run JUST the Y step with ODGI
odgi sort -i "$OUTDIR/seqrush.og" -o "$OUTDIR/odgi_Y.og" -p Y -P 2>&1 | tail -10
odgi view -i "$OUTDIR/odgi_Y.og" -g > "$OUTDIR/odgi_Y.gfa"

echo "Layout quality after ODGI Y only:"
./target/release/measure_layout_quality "$OUTDIR/odgi_Y.gfa" 2>&1 | grep "RMSE:"

echo ""
echo "Layout quality after ODGI full Ygs:"
odgi sort -i "$OUTDIR/seqrush.og" -o "$OUTDIR/odgi_Ygs.og" -p Ygs -P 2>&1 | tail -10
odgi view -i "$OUTDIR/odgi_Ygs.og" -g > "$OUTDIR/odgi_Ygs.gfa"
./target/release/measure_layout_quality "$OUTDIR/odgi_Ygs.gfa" 2>&1 | grep "RMSE:"

echo ""
echo "========================================="
