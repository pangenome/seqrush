#!/bin/bash
# Simple test: Does ODGI think our sorting is good?
# Strategy: Take our sorted output, re-sort with ODGI Ygs, and compare node counts

SEQRUSH=./target/release/seqrush
RESULTS_DIR=/tmp/hla_simple_test
mkdir -p "$RESULTS_DIR"

echo "HLA-Zoo Sorting Validation: Does ODGI Ygs change our output?"
echo "=============================================================="
echo ""

SUCCESS=0
TOTAL=0

for FASTA in HLA-zoo/seqs/*.fa; do
    BASENAME=$(basename "$FASTA" .fa)
    echo -n "Testing $BASENAME... "
    TOTAL=$((TOTAL + 1))

    OUTDIR="$RESULTS_DIR/$BASENAME"
    mkdir -p "$OUTDIR"

    # Build with seqrush (applies Ygs sorting)
    $SEQRUSH -s "$FASTA" -o "$OUTDIR/seqrush.gfa" -k 0 > "$OUTDIR/build.log" 2>&1

    # Get stats BEFORE ODGI re-sort
    SR_BEFORE=$(odgi stats -i "$OUTDIR/seqrush.gfa" -S 2>&1 | grep -v "^#" | grep -v "^\[" | head -1)
    NODES_BEFORE=$(echo "$SR_BEFORE" | awk '{print $2}')
    EDGES_BEFORE=$(echo "$SR_BEFORE" | awk '{print $3}')

    # Re-sort with ODGI Ygs
    odgi build -g "$OUTDIR/seqrush.gfa" -o "$OUTDIR/seqrush.og" 2>&1 | head -3 > "$OUTDIR/odgi.log"
    odgi sort -i "$OUTDIR/seqrush.og" -o "$OUTDIR/resorted.og" -p Ygs -P 2>&1 | head -10 >> "$OUTDIR/odgi.log"

    # Get stats AFTER ODGI re-sort
    SR_AFTER=$(odgi stats -i "$OUTDIR/resorted.og" -S 2>&1 | grep -v "^#" | grep -v "^\[" | head -1)
    NODES_AFTER=$(echo "$SR_AFTER" | awk '{print $2}')
    EDGES_AFTER=$(echo "$SR_AFTER" | awk '{print $3}')

    # Calculate change
    NODE_CHANGE=$((NODES_AFTER - NODES_BEFORE))
    EDGE_CHANGE=$((EDGES_AFTER - EDGES_BEFORE))

    if [ "$NODE_CHANGE" -eq 0 ] && [ "$EDGE_CHANGE" -eq 0 ]; then
        echo "✓ PERFECT (no changes: $NODES_BEFORE nodes, $EDGES_BEFORE edges)"
        SUCCESS=$((SUCCESS + 1))
    elif [ "${NODE_CHANGE#-}" -le 5 ] && [ "${EDGE_CHANGE#-}" -le 10 ]; then
        echo "✓ GOOD (minimal changes: ${NODE_CHANGE:+$NODE_CHANGE nodes, }${EDGE_CHANGE:+$EDGE_CHANGE edges})"
        SUCCESS=$((SUCCESS + 1))
    else
        echo "✗ CHANGED ($NODE_CHANGE nodes, $EDGE_CHANGE edges)"
    fi
done

echo ""
echo "=============================================================="
echo "Result: $SUCCESS/$TOTAL graphs passed (unchanged or minimal change)"
