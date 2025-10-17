#!/bin/bash
# Compare SeqRush Ygs sorting vs ODGI Ygs sorting
# Strategy: Build graph, randomize it, sort with both tools, compare orderings

SORT_GFA=./target/release/sort_gfa
SEQRUSH=./target/release/seqrush
COMPARE=./compare_sorted_gfas.py
RESULTS_DIR=/tmp/sort_comparison
mkdir -p "$RESULTS_DIR"

echo "SeqRush vs ODGI Ygs Sorting Comparison"
echo "======================================"
echo ""

if [ ! -f "$SORT_GFA" ]; then
    echo "Error: sort_gfa binary not found at $SORT_GFA"
    exit 1
fi

if [ ! -f "$COMPARE" ]; then
    echo "Error: comparison script not found at $COMPARE"
    exit 1
fi

SUCCESS=0
TOTAL=0

for FASTA in HLA-zoo/seqs/*.fa; do
    BASENAME=$(basename "$FASTA" .fa)
    echo "Testing $BASENAME..."
    TOTAL=$((TOTAL + 1))

    OUTDIR="$RESULTS_DIR/$BASENAME"
    mkdir -p "$OUTDIR"

    # Step 1: Build initial graph with seqrush
    echo "  [1/5] Building graph from FASTA..."
    $SEQRUSH -s "$FASTA" -o "$OUTDIR/initial.gfa" -k 0 > "$OUTDIR/build.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "  ✗ FAILED: Graph construction failed"
        continue
    fi

    # Step 2: Randomize the graph
    echo "  [2/5] Randomizing graph..."
    odgi build -g "$OUTDIR/initial.gfa" -o "$OUTDIR/initial.og" 2>&1 | head -3 > "$OUTDIR/randomize.log"
    odgi sort -i "$OUTDIR/initial.og" -o "$OUTDIR/randomized.og" -R 2>&1 | head -5 >> "$OUTDIR/randomize.log"
    odgi view -i "$OUTDIR/randomized.og" -g > "$OUTDIR/randomized.gfa" 2>&1
    if [ $? -ne 0 ]; then
        echo "  ✗ FAILED: Randomization failed"
        continue
    fi

    # Step 3: Sort with SeqRush
    echo "  [3/5] Sorting with SeqRush Ygs..."
    $SORT_GFA -i "$OUTDIR/randomized.gfa" -o "$OUTDIR/seqrush_sorted.gfa" > "$OUTDIR/seqrush_sort.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "  ✗ FAILED: SeqRush sorting failed"
        continue
    fi

    # Step 4: Sort with ODGI
    echo "  [4/5] Sorting with ODGI Ygs..."
    odgi build -g "$OUTDIR/randomized.gfa" -o "$OUTDIR/random_for_odgi.og" 2>&1 | head -3 > "$OUTDIR/odgi_sort.log"
    odgi sort -i "$OUTDIR/random_for_odgi.og" -o "$OUTDIR/odgi_sorted.og" -p Ygs -P 2>&1 | head -10 >> "$OUTDIR/odgi_sort.log"
    odgi view -i "$OUTDIR/odgi_sorted.og" -g > "$OUTDIR/odgi_sorted.gfa" 2>&1
    if [ $? -ne 0 ]; then
        echo "  ✗ FAILED: ODGI sorting failed"
        continue
    fi

    # Step 5: Compare the orderings
    echo "  [5/5] Comparing orderings..."
    python3 $COMPARE "$OUTDIR/seqrush_sorted.gfa" "$OUTDIR/odgi_sorted.gfa" > "$OUTDIR/comparison.txt" 2>&1

    # Extract results
    MATCHED=$(grep "Matched:" "$OUTDIR/comparison.txt" | awk '{print $2}')
    IDENTICAL=$(grep "✓.*IDENTICAL" "$OUTDIR/comparison.txt" | wc -l)
    DIFFERENT=$(grep "✗.*differences" "$OUTDIR/comparison.txt" | wc -l)
    TOTAL_PATHS=$(grep "✓\|✗" "$OUTDIR/comparison.txt" | wc -l)

    if [ "$IDENTICAL" -eq "$TOTAL_PATHS" ]; then
        echo "  ✓ PERFECT: All $TOTAL_PATHS paths have identical ordering"
        SUCCESS=$((SUCCESS + 1))
    elif [ "$DIFFERENT" -le 2 ]; then
        echo "  ✓ GOOD: $IDENTICAL/$TOTAL_PATHS paths identical, minor differences in $DIFFERENT paths"
        SUCCESS=$((SUCCESS + 1))
    else
        echo "  ✗ DIFFERS: Only $IDENTICAL/$TOTAL_PATHS paths identical, $DIFFERENT paths differ"
        echo "     See $OUTDIR/comparison.txt for details"
    fi
    echo ""
done

echo "======================================"
echo "Result: $SUCCESS/$TOTAL graphs matched ODGI sorting"
echo "Details in $RESULTS_DIR"
