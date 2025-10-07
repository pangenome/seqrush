#!/bin/bash

echo "=== Testing Simple Compaction Case ==="

# Create a very simple test case
cat > test_simple.fa << 'EOF'
>seq1
ABCDE
EOF

echo "Created sequence: ABCDE"
echo "This should create nodes A->B->C->D->E"
echo

# Run with verbose to see what happens
cargo run --release -- \
    -s test_simple.fa \
    -o test_simple.gfa \
    --verbose \
    2>&1 | grep -E "(nodes|edges|Compacted|component)"

echo
echo "Output graph:"
cat test_simple.gfa

# Clean up
rm -f test_simple.fa test_simple.gfa

echo
echo "=== Test Complete ==="