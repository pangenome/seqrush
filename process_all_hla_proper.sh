#!/bin/bash

# Function to process a single HLA file
process_hla() {
    local fa_file=$1
    local basename=$(basename "$fa_file" .fa)

    echo ""
    echo "Processing $basename..."
    echo "--------------------------"

    # Step 1: Run SeqRush to build graph
    echo "  1. Building graph with SeqRush..."
    cargo run --release --bin seqrush -- \
        -s "$fa_file" \
        -o "HLA_graphs/${basename}.gfa" \
        --odgi-groom \
        -k 0 2>&1 | grep -E "(Graph written|ERROR)" || true

    if [ ! -f "HLA_graphs/${basename}.gfa" ]; then
        echo "  ERROR: Failed to create GFA"
        return 1
    fi

    # Step 2: Convert to ODGI
    echo "  2. Converting to ODGI format..."
    odgi build -g "HLA_graphs/${basename}.gfa" -o "HLA_graphs/${basename}.og" 2>/dev/null

    if [ ! -f "HLA_graphs/${basename}.og" ]; then
        echo "  ERROR: Failed to convert to ODGI"
        return 1
    fi

    # Step 3: Optimize graph
    echo "  3. Optimizing graph..."
    odgi sort -i "HLA_graphs/${basename}.og" -o "HLA_graphs/${basename}.opt.og" -O 2>/dev/null

    # Step 4: Apply path-guided SGD sort (Ygs)
    echo "  4. Applying path-guided SGD sort..."
    odgi sort -p Ygs -i "HLA_graphs/${basename}.opt.og" -o "HLA_graphs/${basename}.Ygs.og" -P 2>/dev/null

    # Step 5: Get stats
    stats=$(odgi stats -i "HLA_graphs/${basename}.Ygs.og" -S 2>/dev/null | tail -1)
    echo "  Stats: $stats"

    # Step 6: Create visualizations
    echo "  5. Creating visualizations..."

    # SeqRush sorted version
    odgi viz -i "HLA_graphs/${basename}.opt.og" \
        -o "HLA_visualizations/${basename}.seqrush.png" \
        -x 2000 -y 600 -s 4 2>/dev/null

    # ODGI Ygs sorted version
    odgi viz -i "HLA_graphs/${basename}.Ygs.og" \
        -o "HLA_visualizations/${basename}.Ygs.png" \
        -x 2000 -y 600 -s 4 2>/dev/null

    # Color by path prefix version
    odgi viz -i "HLA_graphs/${basename}.Ygs.og" \
        -o "HLA_visualizations/${basename}.colored.png" \
        -x 2000 -y 600 -s : -P 2>/dev/null

    echo "  ✓ Complete"
    return 0
}

# Create directories
mkdir -p HLA_graphs HLA_visualizations

echo "============================================"
echo "Processing all HLA sequences with SeqRush"
echo "============================================"

success_count=0
fail_count=0

# Process all HLA files
for fa_file in HLA-zoo/seqs/*.fa; do
    if process_hla "$fa_file"; then
        ((success_count++))
    else
        ((fail_count++))
    fi
done

# Create an HTML index for easy viewing
cat > HLA_visualizations/index.html <<EOF
<!DOCTYPE html>
<html>
<head>
    <title>HLA Graph Visualizations</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #333; }
        .gene { margin-bottom: 30px; border: 1px solid #ddd; padding: 15px; }
        .gene h2 { color: #666; }
        img { max-width: 100%; height: auto; margin: 10px 0; }
        .comparison { display: grid; grid-template-columns: 1fr 1fr; gap: 10px; }
        .comparison img { width: 100%; }
    </style>
</head>
<body>
    <h1>HLA Graph Visualizations</h1>
    <p>Generated with SeqRush and ODGI</p>
EOF

for png in HLA_visualizations/*.Ygs.png; do
    basename=$(basename "$png" .Ygs.png)
    echo "    <div class='gene'>" >> HLA_visualizations/index.html
    echo "        <h2>$basename</h2>" >> HLA_visualizations/index.html
    echo "        <div class='comparison'>" >> HLA_visualizations/index.html
    echo "            <div>" >> HLA_visualizations/index.html
    echo "                <h3>SeqRush Sorted</h3>" >> HLA_visualizations/index.html
    echo "                <img src='${basename}.seqrush.png' alt='SeqRush sorted'>" >> HLA_visualizations/index.html
    echo "            </div>" >> HLA_visualizations/index.html
    echo "            <div>" >> HLA_visualizations/index.html
    echo "                <h3>ODGI Ygs Sorted</h3>" >> HLA_visualizations/index.html
    echo "                <img src='${basename}.Ygs.png' alt='ODGI Ygs sorted'>" >> HLA_visualizations/index.html
    echo "            </div>" >> HLA_visualizations/index.html
    echo "        </div>" >> HLA_visualizations/index.html
    echo "        <h3>Colored by Path</h3>" >> HLA_visualizations/index.html
    echo "        <img src='${basename}.colored.png' alt='Colored by path'>" >> HLA_visualizations/index.html
    echo "    </div>" >> HLA_visualizations/index.html
done

echo "</body></html>" >> HLA_visualizations/index.html

# Summary
echo ""
echo "============================================"
echo "SUMMARY"
echo "============================================"
echo "Successfully processed: $success_count"
echo "Failed: $fail_count"
echo ""
echo "Visualizations created:"
ls HLA_visualizations/*.png 2>/dev/null | wc -l
echo ""
echo "View all visualizations: open HLA_visualizations/index.html"