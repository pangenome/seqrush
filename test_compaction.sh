#!/bin/bash

echo "=== Testing Node Compaction ==="

# Create a test FASTA with sequences that should create linear chains
cat > test_compaction.fa << 'EOF'
>seq1
ATCGATCGATCGATCG
>seq2
ATCGATCGATCGATCGTTTTGGGG
EOF

echo "Created test sequences that should form linear chains"
echo

# Test with compaction disabled
echo "1. Running with --no-compact (compaction disabled):"
cargo run --release -- \
    -s test_compaction.fa \
    -o test_no_compact.gfa \
    --no-compact \
    2>&1 | grep -E "(nodes|edges|Compacted)"

echo
echo "Graph without compaction:"
head -n 20 test_no_compact.gfa

echo
echo "2. Running with compaction enabled (default):"
cargo run --release -- \
    -s test_compaction.fa \
    -o test_with_compact.gfa \
    --verbose \
    2>&1 | grep -E "(nodes|edges|Compacted|compaction)"

echo
echo "Graph with compaction:"
head -n 20 test_with_compact.gfa

# Compare sizes
echo
echo "Comparing graph sizes:"
echo -n "Without compaction: "
grep "^S" test_no_compact.gfa | wc -l | xargs echo "nodes"
echo -n "With compaction: "
grep "^S" test_with_compact.gfa | wc -l | xargs echo "nodes"

# Clean up
rm -f test_compaction.fa test_no_compact.gfa test_with_compact.gfa

echo
echo "=== Compaction Test Complete ==="