#!/bin/bash
# Proper test: Compare actual node ORDERING, not just counts

FASTA="HLA-zoo/seqs/B-3106.fa"
OUTDIR="/tmp/ordering_test"
mkdir -p "$OUTDIR"

echo "Testing if node ORDERING matches between seqrush and ODGI..."
echo ""

# Build with seqrush
./target/release/seqrush -s "$FASTA" -o "$OUTDIR/seqrush.gfa" -k 0 > /dev/null 2>&1

# Get seqrush node order from first path
SR_ORDER=$(grep "^P" "$OUTDIR/seqrush.gfa" | head -1 | cut -f3 | tr ',' '\n' | head -30)

# Re-sort with ODGI
odgi build -g "$OUTDIR/seqrush.gfa" -o "$OUTDIR/seqrush.og" 2>/dev/null
odgi sort -i "$OUTDIR/seqrush.og" -o "$OUTDIR/resorted.og" -p Ygs -P 2>/dev/null
odgi view -i "$OUTDIR/resorted.og" -g > "$OUTDIR/resorted.gfa"

# Get ODGI node order from first path
ODGI_ORDER=$(grep "^P" "$OUTDIR/resorted.gfa" | head -1 | cut -f3 | tr ',' '\n' | head -30)

echo "SeqRush path (first 30 steps):"
echo "$SR_ORDER"
echo ""
echo "ODGI re-sorted path (first 30 steps):"
echo "$ODGI_ORDER"
echo ""

# Compare
if [ "$SR_ORDER" = "$ODGI_ORDER" ]; then
    echo "✓ Path ordering is IDENTICAL"
else
    echo "✗ Path ordering DIFFERS"
    echo ""
    echo "Differences:"
    diff <(echo "$SR_ORDER") <(echo "$ODGI_ORDER") || true
fi
