#!/bin/bash

echo "=== Testing Orientation Check Performance ==="

# Create test dataset with many sequences that need orientation checking
cat > test_perf.fa << 'EOF'
>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
>seq3
TACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq4
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTA
>seq5
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
>seq6
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
>seq7
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq8
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
EOF

echo "Created test dataset with 8 sequences (28 pairwise alignments)"
echo

# Time the original approach (align in both orientations)
echo "1. Testing without orientation check (using test mode to simulate old behavior):"
time cargo run --release -- \
    -s test_perf.fa \
    -o test_perf_no_check.gfa \
    --test-mode \
    2>&1 | grep -E "(Building graph|Aligning)"

echo
echo "2. Testing with orientation check (new default behavior):"
time cargo run --release -- \
    -s test_perf.fa \
    -o test_perf_with_check.gfa \
    2>&1 | grep -E "(Building graph|Aligning)"

echo
echo "3. Testing with cheaper orientation scores (0,1,0,0 - only penalize mismatches):"
time cargo run --release -- \
    -s test_perf.fa \
    -o test_perf_cheap_check.gfa \
    --orientation-scores 0,1,0,0 \
    2>&1 | grep -E "(Building graph|Aligning)"

# Compare outputs
echo
echo "Comparing output sizes:"
wc -l test_perf_*.gfa

# Clean up
rm -f test_perf.fa test_perf_*.gfa

echo
echo "=== Performance Test Complete ==="