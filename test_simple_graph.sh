#!/bin/bash
set -e

# Create a simple test case with 3 sequences
cat > test_simple.fa << EOF
>seq1
ACTGACTG
>seq2
ACTGACTG
>seq3
ACTGTTTG
EOF

echo "=== Test sequences ==="
cat test_simple.fa

echo -e "\n=== Generate alignments with allwave ==="
allwave -t 1 -i test_simple.fa > test_simple.paf
cat test_simple.paf

echo -e "\n=== Run seqwish ==="
seqwish -s test_simple.fa -p test_simple.paf -g test_simple_seqwish.gfa
odgi build -g test_simple_seqwish.gfa -o test_simple_seqwish.odgi
odgi stats -i test_simple_seqwish.odgi -S

echo -e "\n=== Run seqrush ==="
cargo run --release --bin seqrush -- -s test_simple.fa -o test_simple_seqrush.gfa -k 0 --no-compact -v
odgi build -g test_simple_seqrush.gfa -o test_simple_seqrush.odgi
odgi stats -i test_simple_seqrush.odgi -S

echo -e "\n=== Compare GFA files ==="
echo "Seqwish nodes:"
grep "^S" test_simple_seqwish.gfa | wc -l
echo "SeqRush nodes:"
grep "^S" test_simple_seqrush.gfa | wc -l