#\!/bin/bash

# Test all HLA sequences with SeqRush
echo "Testing all HLA sequences with SeqRush..."

# Create output directory
mkdir -p hla_outputs

# Counter for successful builds
success=0
total=0

# Process each .fa file
for fa_file in HLA-zoo/seqs/*.fa; do
    if [ -f "$fa_file" ]; then
        base=$(basename "$fa_file" .fa)
        echo "Processing $base..."
        total=$((total + 1))
        
        # Run seqrush
        if cargo run --release --bin seqrush -- -s "$fa_file" -o "hla_outputs/${base}.gfa" -k 0 2>&1 | grep -q "Graph written"; then
            echo "✓ $base completed successfully"
            success=$((success + 1))
        else
            echo "✗ $base failed"
            # Show the error
            cargo run --release --bin seqrush -- -s "$fa_file" -o "hla_outputs/${base}.gfa" -k 0 2>&1 | tail -10
        fi
        echo "---"
    fi
done

echo "Summary: $success/$total sequences processed successfully"
EOF < /dev/null