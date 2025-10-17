#!/bin/bash
set -e

# Test all HLA-zoo graphs: compare seqrush Ygs vs ODGI Ygs sorting

SEQRUSH=./target/release/seqrush
RESULTS_DIR=/tmp/hla_test_results
mkdir -p "$RESULTS_DIR"

# Summary file
SUMMARY="$RESULTS_DIR/summary.txt"
echo "HLA-Zoo Sorting Comparison: SeqRush vs ODGI" > "$SUMMARY"
echo "=============================================" >> "$SUMMARY"
echo "" >> "$SUMMARY"

SUCCESS_COUNT=0
TOTAL_COUNT=0

for FASTA in HLA-zoo/seqs/*.fa; do
    BASENAME=$(basename "$FASTA" .fa)
    echo "Testing $BASENAME..."
    TOTAL_COUNT=$((TOTAL_COUNT + 1))

    OUTDIR="$RESULTS_DIR/$BASENAME"
    mkdir -p "$OUTDIR"

    # 1. Build with seqrush (default Ygs sorting)
    echo "  Building with seqrush..."
    $SEQRUSH -s "$FASTA" -o "$OUTDIR/seqrush.gfa" -k 0 > "$OUTDIR/seqrush.log" 2>&1

    # 2. Build with seqwish for comparison baseline
    echo "  Building baseline graph..."
    allwave -i "$FASTA" -t 4 > "$OUTDIR/alignments.paf" 2>&1
    seqwish -s "$FASTA" -p "$OUTDIR/alignments.paf" -g "$OUTDIR/baseline.gfa" 2> "$OUTDIR/seqwish.log"

    # 3. Convert to ODGI format and sort with Ygs
    echo "  Sorting baseline with ODGI..."
    odgi build -g "$OUTDIR/baseline.gfa" -o "$OUTDIR/baseline.og" 2>&1 | head -5 > "$OUTDIR/odgi-build.log"
    odgi sort -i "$OUTDIR/baseline.og" -o "$OUTDIR/odgi-ygs.og" -p Ygs -P 2>&1 | head -20 >> "$OUTDIR/odgi-build.log"
    
    # 4. Get stats for both (after compaction)
    echo "  Collecting stats..."
    SEQRUSH_STATS=$(odgi stats -i "$OUTDIR/seqrush.gfa" -S 2>&1 | grep -v "^#" | grep -v "^\[" | head -1)
    ODGI_STATS=$(odgi stats -i "$OUTDIR/odgi-ygs.og" -S 2>&1 | grep -v "^#" | grep -v "^\[" | head -1)

    # Parse stats (format: length nodes edges paths steps)
    SR_NODES=$(echo "$SEQRUSH_STATS" | awk '{print $2}')
    SR_EDGES=$(echo "$SEQRUSH_STATS" | awk '{print $3}')
    ODGI_NODES=$(echo "$ODGI_STATS" | awk '{print $2}')
    ODGI_EDGES=$(echo "$ODGI_STATS" | awk '{print $3}')

    # Calculate differences
    NODE_DIFF=$((SR_NODES - ODGI_NODES))
    EDGE_DIFF=$((SR_EDGES - ODGI_EDGES))
    
    # Avoid division by zero
    if [ "$ODGI_NODES" -gt 0 ]; then
        NODE_DIFF_PCT=$(awk "BEGIN {printf \"%.1f\", ($NODE_DIFF * 100.0) / $ODGI_NODES}")
    else
        NODE_DIFF_PCT="N/A"
    fi

    # Record results
    echo "$BASENAME:" >> "$SUMMARY"
    echo "  SeqRush:  $SR_NODES nodes, $SR_EDGES edges" >> "$SUMMARY"
    echo "  ODGI:     $ODGI_NODES nodes, $ODGI_EDGES edges" >> "$SUMMARY"
    echo "  Diff:     $NODE_DIFF nodes ($NODE_DIFF_PCT%), $EDGE_DIFF edges" >> "$SUMMARY"

    # Check if within acceptable range (±5%)
    if [ "$NODE_DIFF_PCT" != "N/A" ]; then
        # Use awk for comparison since bc might not be available
        PASS=$(awk "BEGIN {print (($NODE_DIFF_PCT > -5.0) && ($NODE_DIFF_PCT < 5.0)) ? 1 : 0}")
        if [ "$PASS" -eq 1 ]; then
            echo "  Status:   ✓ PASS (within 5%)" >> "$SUMMARY"
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        else
            echo "  Status:   ✗ FAIL (difference > 5%)" >> "$SUMMARY"
        fi
    fi
    echo "" >> "$SUMMARY"

    # Create visualizations
    echo "  Creating visualizations..."
    odgi viz -i "$OUTDIR/seqrush.gfa" -o "$OUTDIR/seqrush.png" -x 1500 -y 500 -P 2>/dev/null || true
    odgi viz -i "$OUTDIR/odgi-ygs.og" -o "$OUTDIR/odgi.png" -x 1500 -y 500 -P 2>/dev/null || true
done

echo "========================================" >> "$SUMMARY"
echo "Results: $SUCCESS_COUNT/$TOTAL_COUNT graphs within 5% node count" >> "$SUMMARY"
echo "" >> "$SUMMARY"

cat "$SUMMARY"

echo ""
echo "Detailed results in: $RESULTS_DIR"
echo "Summary: $SUMMARY"
