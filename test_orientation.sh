#!/bin/bash

echo "=== Testing Orientation Check Feature ==="

# Create test sequences - one forward and one reverse complement
cat > test_orientation.fa << 'EOF'
>seq1
ATCGATCGATCGATCGATCG
>seq2_forward
ATCGATCGATCGATCGATCG
>seq3_reverse
CGATCGATCGATCGATCGAT
EOF

echo "Test sequences created:"
echo "  seq1: ATCGATCGATCGATCGATCG"
echo "  seq2_forward: ATCGATCGATCGATCGATCG (identical to seq1)"
echo "  seq3_reverse: CGATCGATCGATCGATCGAT (reverse complement of seq1)"
echo

# Run with verbose output to see orientation check scores
echo "Running seqrush with default orientation scores (0,1,1,1)..."
cargo run --release -- \
    -s test_orientation.fa \
    -o test_orientation_default.gfa \
    --verbose \
    2>&1 | grep -E "(Aligning|Orientation check|Full alignment score)"

echo
echo "Running seqrush with custom orientation scores (0,2,2,2)..."
cargo run --release -- \
    -s test_orientation.fa \
    -o test_orientation_custom.gfa \
    --orientation-scores 0,2,2,2 \
    --verbose \
    2>&1 | grep -E "(Aligning|Orientation check|Full alignment score)"

# Check the output
echo
echo "Checking output graphs..."
echo "Default orientation scores graph:"
head -n 5 test_orientation_default.gfa

echo
echo "Custom orientation scores graph:"
head -n 5 test_orientation_custom.gfa

# Clean up
rm -f test_orientation.fa test_orientation_default.gfa test_orientation_custom.gfa

echo
echo "=== Test Complete ==="