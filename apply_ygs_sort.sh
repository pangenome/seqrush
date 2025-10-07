#!/bin/bash
# Script to apply path-guided SGD sort to all HLA graphs and regenerate visualizations

set -e

echo "============================================"
echo "Applying path-guided SGD sort to all HLA graphs"
echo "============================================"

mkdir -p HLA_graphs_sorted

for og_file in HLA_graphs/*.sorted.og; do
    basename=$(basename "$og_file" .sorted.og)

    if [ ! -f "$og_file" ]; then
        echo "Skipping $basename - no sorted.og file"
        continue
    fi

    echo ""
    echo "Processing $basename..."
    echo "--------------------------"

    # Step 1: Apply path-guided SGD sort
    echo "  1. Applying path-guided SGD sort..."
    odgi sort -p Ygs -i "$og_file" -o "HLA_graphs_sorted/${basename}.Ygs.og" -P 2>/dev/null

    # Step 2: Get stats
    stats=$(odgi stats -i "HLA_graphs_sorted/${basename}.Ygs.og" -S 2>/dev/null | tail -1)
    echo "  Stats: $stats"

    # Step 3: Create visualization
    echo "  2. Creating visualization..."
    odgi viz -i "HLA_graphs_sorted/${basename}.Ygs.og" \
        -o "HLA_graphs_sorted/${basename}.png" \
        -x 1500 -y 500 -P 2>/dev/null

    # Step 4: Convert to GFA for inspection
    echo "  3. Converting to GFA..."
    odgi view -i "HLA_graphs_sorted/${basename}.Ygs.og" -g > "HLA_graphs_sorted/${basename}.gfa" 2>/dev/null

    echo "  ✓ Complete: HLA_graphs_sorted/${basename}.png"
done

echo ""
echo "============================================"
echo "All HLA graphs sorted with path-guided SGD!"
echo "Sorted graphs and visualizations saved in HLA_graphs_sorted/"
echo "============================================"

# Create an HTML index to view all graphs
cat > HLA_graphs_sorted/index.html << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>HLA Graphs - Path-guided SGD Sorted</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #333; }
        .graph { margin: 20px 0; border: 1px solid #ddd; padding: 10px; }
        .graph h2 { margin-top: 0; color: #666; }
        img { max-width: 100%; height: auto; }
    </style>
</head>
<body>
    <h1>HLA Graphs - Path-guided SGD Sorted</h1>
EOF

for png in HLA_graphs_sorted/*.png; do
    basename=$(basename "$png" .png)
    echo "    <div class='graph'>" >> HLA_graphs_sorted/index.html
    echo "        <h2>$basename</h2>" >> HLA_graphs_sorted/index.html
    echo "        <img src='$(basename $png)' alt='$basename'>" >> HLA_graphs_sorted/index.html
    echo "    </div>" >> HLA_graphs_sorted/index.html
done

echo "</body></html>" >> HLA_graphs_sorted/index.html

echo ""
echo "HTML index created at HLA_graphs_sorted/index.html"