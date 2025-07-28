#!/bin/bash
set -e

echo "Testing PAF output functionality..."

# Create test sequences
cat > test_sequences.fa << EOF
>seq1
ATCGATCGATCGATCG
>seq2
ATCGATGGATCGATCG
>seq3
CGATTGATCGATCGAT
EOF

# Run seqrush with PAF output
echo "Running seqrush with PAF output..."
cargo run -- -s test_sequences.fa -o test_output.gfa --output-alignments test_alignments.paf

# Check if PAF file was created
if [ -f test_alignments.paf ]; then
    echo "✓ PAF file created successfully"
    echo "PAF contents:"
    cat test_alignments.paf
else
    echo "✗ PAF file not created"
    exit 1
fi

# Check if GFA file was created
if [ -f test_output.gfa ]; then
    echo "✓ GFA file created successfully"
else
    echo "✗ GFA file not created"
    exit 1
fi

# Clean up
rm -f test_sequences.fa test_output.gfa test_alignments.paf

echo "Test completed successfully!"