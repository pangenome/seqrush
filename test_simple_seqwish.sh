#!/bin/bash

# Create simple test sequences
echo ">seq1" > simple_test.fa
echo "ACGT" >> simple_test.fa
echo ">seq2" >> simple_test.fa
echo "ACGT" >> simple_test.fa

# Generate alignments with allwave
allwave -i simple_test.fa > simple_test.paf

echo "PAF content:"
cat simple_test.paf

# Run seqwish
seqwish -s simple_test.fa -p simple_test.paf -g simple_seqwish.gfa

echo -e "\nSeqwish output:"
cat simple_seqwish.gfa

# Check stats
echo -e "\nSeqwish stats:"
odgi stats -i simple_seqwish.gfa -S