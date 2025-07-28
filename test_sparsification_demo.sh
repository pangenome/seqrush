#!/bin/bash

echo "=== SeqRush Sparsification Demo ==="
echo

# Create a test dataset with 20 sequences
cat > test_sparse_demo.fa << EOF
>seq1
ATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
ATCGATCGATCGATCGATCGATCGATCGATCG
>seq3
ATCGATCGATCGATCGATCGATCGATCGATCG
>seq4
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>seq5
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>seq6
TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGT
>seq7
TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGT
>seq8
CGATCGATCGATCGATCGATCGATCGATCGAT
>seq9
CGATCGATCGATCGATCGATCGATCGATCGAT
>seq10
ATATATATATATATATATATATATATATATAT
>seq11
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
>seq12
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
>seq13
CACACACACACACACACACACACACACACACA
>seq14
AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG
>seq15
TCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTC
>seq16
GATAGATAGATAGATAGATAGATAGATAGATA
>seq17
CTATCTATCTATCTATCTATCTATCTATCTAT
>seq18
ACGTACGTACGTACGTACGTACGTACGTACGT
>seq19
TGCATGCATGCATGCATGCATGCATGCATGCA
>seq20
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
EOF

echo "Created test dataset with 20 sequences"
echo "Total possible alignments: $(echo "20*19/2" | bc) = 190"
echo

# Run with different sparsification levels
echo "1. Running with no sparsification (100% of alignments)..."
./target/release/seqrush -s test_sparse_demo.fa -o test_100.gfa -t 1 -x 1.0 -v 2>&1 | grep -c "Aligning seq" | xargs -I {} echo "   Alignments performed: {}"

echo
echo "2. Running with 50% sparsification..."
./target/release/seqrush -s test_sparse_demo.fa -o test_50.gfa -t 1 -x 0.5 -v 2>&1 | grep -c "Aligning seq" | xargs -I {} echo "   Alignments performed: {}"

echo
echo "3. Running with 20% sparsification..."
./target/release/seqrush -s test_sparse_demo.fa -o test_20.gfa -t 1 -x 0.2 -v 2>&1 | grep -c "Aligning seq" | xargs -I {} echo "   Alignments performed: {}"

echo
echo "4. Running with auto sparsification (Erdős-Rényi model)..."
./target/release/seqrush -s test_sparse_demo.fa -o test_auto.gfa -t 1 -x auto -v 2>&1 | tee /tmp/auto_run.log | grep -E "Auto sparsification:|Aligning seq" | head -2
grep -c "Aligning seq" /tmp/auto_run.log | xargs -I {} echo "   Alignments performed: {}"

echo
echo "Comparing graph sizes:"
for f in test_100.gfa test_50.gfa test_20.gfa test_auto.gfa; do
    nodes=$(grep -c "^S" $f)
    edges=$(grep -c "^L" $f || echo 0)
    echo "   $f: $nodes nodes, $edges edges"
done

# Cleanup
rm -f test_sparse_demo.fa test_*.gfa /tmp/auto_run.log

echo
echo "=== Demo Complete ==="