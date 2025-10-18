#!/bin/bash
# Compare layout quality between SeqRush and ODGI
# Tests different SGD iteration counts to find optimal parameters

set -e

SEQRUSH=./target/release/seqrush
MEASURE=./target/release/measure_layout_quality
RESULTS_DIR=/tmp/layout_benchmark
mkdir -p "$RESULTS_DIR"

# Test with a single HLA graph first
TEST_FASTA="HLA-zoo/seqs/B-3106.fa"
BASENAME=$(basename "$TEST_FASTA" .fa)

echo "================================================================================"
echo "Layout Quality Benchmark: SeqRush vs ODGI"
echo "================================================================================"
echo ""
echo "Test input: $TEST_FASTA"
echo ""

OUTDIR="$RESULTS_DIR/$BASENAME"
mkdir -p "$OUTDIR"

# Test 1: SeqRush with default parameters (iter_max=100)
echo "--- Test 1: SeqRush (default: 100 iterations) ---"
$SEQRUSH -s "$TEST_FASTA" -o "$OUTDIR/seqrush_100.gfa" -k 0 > "$OUTDIR/seqrush_100.log" 2>&1
echo "Layout quality:"
$MEASURE "$OUTDIR/seqrush_100.gfa" 2>&1 | grep -A20 "^Overall metrics:"
echo ""

# Test 2: ODGI with Ygs pipeline
echo "--- Test 2: ODGI (Ygs pipeline, default 100 iterations) ---"
# First build graph with allwave
allwave -t 8 -i "$TEST_FASTA" > "$OUTDIR/odgi.paf" 2>&1
# Build with seqwish (no sorting)
seqwish -s "$TEST_FASTA" -p "$OUTDIR/odgi.paf" -g "$OUTDIR/odgi_unsorted.gfa" 2>&1 | head -5
# Sort with ODGI Ygs
odgi build -g "$OUTDIR/odgi_unsorted.gfa" -o "$OUTDIR/odgi_unsorted.og" 2>&1 | head -3
odgi sort -i "$OUTDIR/odgi_unsorted.og" -o "$OUTDIR/odgi_sorted.og" -p Ygs -P > "$OUTDIR/odgi_sort.log" 2>&1
# Convert back to GFA for measurement
odgi view -i "$OUTDIR/odgi_sorted.og" -g > "$OUTDIR/odgi_sorted.gfa"
echo "Layout quality:"
$MEASURE "$OUTDIR/odgi_sorted.gfa" 2>&1 | grep -A20 "^Overall metrics:"
echo ""

echo "================================================================================"
echo "Comparison complete. Full results in: $OUTDIR"
echo "================================================================================"
