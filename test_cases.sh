#\!/bin/bash

echo "=== Test Case 1: Small identical sequences (should merge completely) ==="
cat > identical.fa << 'EOL'
>seq1
ATCGATCGATCGATCG
>seq2
ATCGATCGATCGATCG
>seq3
ATCGATCGATCGATCG
EOL

echo "Expected: Very few components (most characters merged)"
cargo run -- -s identical.fa -o identical.gfa --verbose --max-iterations 3
echo ""

echo "=== Test Case 2: Medium sequences with high similarity ==="
cat > medium.fa << 'EOL'
>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq3
ATCGATCGATCGATCGATCGATCGCTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOL

echo "Expected: Moderate component reduction"
cargo run -- -s medium.fa -o medium.gfa --verbose --max-iterations 3
echo ""

echo "=== Test Case 3: Large sequences (demonstrates the problem) ==="
# Take first 3 sequences from the problematic file
head -6 partition714.del.fa > large.fa

echo "Expected: PROBLEM - too many components due to full-sequence alignment"
timeout 30s cargo run -- -s large.fa -o large.gfa --verbose --max-iterations 2
echo ""
