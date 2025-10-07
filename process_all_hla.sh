#!/bin/bash

# Create output directories
mkdir -p HLA_graphs
mkdir -p HLA_visualizations

echo "Processing all HLA sequences..."

# Process each HLA file
for fa_file in HLA-zoo/seqs/*.fa; do
    basename=$(basename "$fa_file" .fa)
    echo ""
    echo "========================================="
    echo "Processing $basename"
    echo "========================================="

    # Run SeqRush with SGD sorting and ODGI-style grooming
    echo "Building graph with SeqRush..."
    cargo run --release --bin seqrush -- \
        -s "$fa_file" \
        -o "HLA_graphs/${basename}.gfa" \
        --odgi-groom \
        --sgd-sort \
        -k 0 2>&1 | tail -5

    if [ ! -f "HLA_graphs/${basename}.gfa" ]; then
        echo "ERROR: Failed to create GFA for $basename"
        continue
    fi

    # Convert to ODGI format
    echo "Converting to ODGI..."
    odgi build -g "HLA_graphs/${basename}.gfa" -o "HLA_graphs/${basename}.og" 2>/dev/null

    if [ ! -f "HLA_graphs/${basename}.og" ]; then
        echo "ERROR: Failed to convert to ODGI for $basename"
        continue
    fi

    # Get graph statistics
    echo "Graph statistics:"
    odgi stats -i "HLA_graphs/${basename}.og" -S

    # Create visualization
    echo "Creating visualization..."
    odgi viz -i "HLA_graphs/${basename}.og" \
        -o "HLA_visualizations/${basename}.png" \
        -x 1000 \
        -s \
        -P 10 2>/dev/null

    if [ -f "HLA_visualizations/${basename}.png" ]; then
        echo "✓ Visualization saved to HLA_visualizations/${basename}.png"
    else
        echo "WARNING: Failed to create visualization for $basename"
    fi
done

# Summary
echo ""
echo "========================================="
echo "SUMMARY"
echo "========================================="
echo "GFA files created:"
ls HLA_graphs/*.gfa 2>/dev/null | wc -l
echo ""
echo "ODGI files created:"
ls HLA_graphs/*.og 2>/dev/null | wc -l
echo ""
echo "Visualizations created:"
ls HLA_visualizations/*.png 2>/dev/null | wc -l

# Check for any failures
echo ""
echo "Checking for issues..."
for fa_file in HLA-zoo/seqs/*.fa; do
    basename=$(basename "$fa_file" .fa)
    if [ ! -f "HLA_graphs/${basename}.og" ]; then
        echo "Failed: $basename"
    fi
done

echo ""
echo "Done! Check HLA_visualizations/ for the PNG files."