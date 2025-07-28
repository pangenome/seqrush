#!/bin/bash
set -e

# Create test sequences where seq2 is RC of seq1
cat > test_rc_debug.fa << EOF
>seq1
AAACCCGGGTTT
>seq2
AAACCCGGGTTT
EOF

# Create manual PAF with RC alignment
# seq1 forward matches seq2 reverse complement perfectly
cat > test_rc_debug.paf << EOF  
seq1	12	0	12	-	seq2	12	0	12	12	12	60	gi:f:1.0	cg:Z:12=
EOF

echo "=== Test FASTA ==="
cat test_rc_debug.fa

echo -e "\n=== Test PAF (seq1 RC matches seq2 forward) ==="
cat test_rc_debug.paf

echo -e "\n=== Running SeqRush ==="
RUST_LOG=debug cargo run --release --bin seqrush -- -s test_rc_debug.fa -o test_rc_debug_seqrush.gfa -k 0 -v --no-compact --output-alignments test_rc_debug_seqrush.paf 2>&1 | grep -E "(alignment|UNITE|CIGAR|nodes|Processing)"

echo -e "\n=== SeqRush GFA ==="
cat test_rc_debug_seqrush.gfa