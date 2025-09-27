#!/bin/bash

# Function to process a single HLA file
process_hla() {
    local fa_file=$1
    local basename=$(basename "$fa_file" .fa)

    echo "Processing $basename..."

    # Run SeqRush
    cargo run --release --bin seqrush -- \
        -s "$fa_file" \
        -o "HLA_graphs/${basename}.gfa" \
        --odgi-groom \
        -k 0 2>&1 | grep -E "(Graph written|ERROR)" || true

    if [ ! -f "HLA_graphs/${basename}.gfa" ]; then
        echo "  ERROR: Failed to create GFA for $basename"
        return 1
    fi

    # Convert to ODGI
    odgi build -g "HLA_graphs/${basename}.gfa" -o "HLA_graphs/${basename}.og" 2>/dev/null

    if [ ! -f "HLA_graphs/${basename}.og" ]; then
        echo "  ERROR: Failed to convert to ODGI for $basename"
        return 1
    fi

    # Optimize for visualization
    odgi sort -i "HLA_graphs/${basename}.og" -o "HLA_graphs/${basename}.sorted.og" -O 2>/dev/null

    if [ ! -f "HLA_graphs/${basename}.sorted.og" ]; then
        echo "  WARNING: Failed to sort, using unsorted"
        cp "HLA_graphs/${basename}.og" "HLA_graphs/${basename}.sorted.og"
    fi

    # Get stats
    stats=$(odgi stats -i "HLA_graphs/${basename}.sorted.og" -S 2>/dev/null | head -1)
    echo "  Stats: $stats"

    # Create visualization
    odgi viz -i "HLA_graphs/${basename}.sorted.og" \
        -o "HLA_visualizations/${basename}.png" \
        -x 1500 \
        -y 500 \
        -s 4 2>/dev/null

    if [ -f "HLA_visualizations/${basename}.png" ]; then
        echo "  ✓ Visualization created"
    else
        echo "  WARNING: Visualization failed"
    fi

    return 0
}

# Create directories
mkdir -p HLA_graphs HLA_visualizations

# Process all HLA files
echo "Processing all HLA sequences..."
echo "================================"

success_count=0
fail_count=0

for fa_file in HLA-zoo/seqs/*.fa; do
    if process_hla "$fa_file"; then
        ((success_count++))
    else
        ((fail_count++))
    fi
done

# Summary
echo ""
echo "================================"
echo "SUMMARY"
echo "================================"
echo "Successfully processed: $success_count"
echo "Failed: $fail_count"
echo ""
echo "Visualizations created:"
ls -la HLA_visualizations/*.png 2>/dev/null | wc -l
echo ""
echo "Graph files:"
ls -lh HLA_graphs/*.sorted.og 2>/dev/null | head -10

echo ""
echo "Done! Check HLA_visualizations/ for the PNG files."