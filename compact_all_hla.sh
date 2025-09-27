#!/bin/bash
# Script to properly compact all HLA graphs using ODGI unchop

set -e

echo "============================================"
echo "Compacting all HLA graphs with ODGI unchop"
echo "============================================"

mkdir -p HLA_graphs_compacted

for fa_file in HLA-zoo/seqs/*.fa; do
    basename=$(basename "$fa_file" .fa)
    echo ""
    echo "Processing $basename..."
    echo "--------------------------"

    # Step 1: Build graph with SeqRush (no compact flag)
    echo "  1. Building graph with SeqRush..."
    cargo run --release --bin seqrush -- \
        -s "$fa_file" \
        -o "HLA_graphs_compacted/${basename}.raw.gfa" \
        -k 0 2>&1 | grep -E "(Graph written|ERROR)" || true

    if [ ! -f "HLA_graphs_compacted/${basename}.raw.gfa" ]; then
        echo "  ERROR: Failed to create GFA"
        continue
    fi

    # Step 2: Convert to ODGI
    echo "  2. Converting to ODGI format..."
    odgi build -g "HLA_graphs_compacted/${basename}.raw.gfa" -o "HLA_graphs_compacted/${basename}.og" 2>/dev/null

    # Step 3: Sort the graph
    echo "  3. Sorting graph..."
    odgi sort -i "HLA_graphs_compacted/${basename}.og" -o "HLA_graphs_compacted/${basename}.sorted.og" 2>/dev/null

    # Step 4: Compact with unchop
    echo "  4. Compacting with odgi unchop..."
    odgi unchop -i "HLA_graphs_compacted/${basename}.sorted.og" -o "HLA_graphs_compacted/${basename}.compacted.og" 2>/dev/null

    # Step 5: Apply path-guided SGD sort
    echo "  5. Applying path-guided SGD sort..."
    odgi sort -p Ygs -i "HLA_graphs_compacted/${basename}.compacted.og" -o "HLA_graphs_compacted/${basename}.final.og" -P 2>/dev/null

    # Step 6: Convert back to GFA
    echo "  6. Converting back to GFA..."
    odgi view -i "HLA_graphs_compacted/${basename}.final.og" -g > "HLA_graphs_compacted/${basename}.gfa" 2>/dev/null

    # Step 7: Get stats
    stats=$(odgi stats -i "HLA_graphs_compacted/${basename}.final.og" -S 2>/dev/null | tail -1)
    echo "  Stats: $stats"

    # Step 8: Create visualization
    echo "  7. Creating visualization..."
    odgi viz -i "HLA_graphs_compacted/${basename}.final.og" \
        -o "HLA_graphs_compacted/${basename}.png" \
        -x 1500 -y 500 -P 5 2>/dev/null

    echo "  ✓ Complete: HLA_graphs_compacted/${basename}.gfa"
done

echo ""
echo "============================================"
echo "All HLA graphs compacted successfully!"
echo "Compacted graphs saved in HLA_graphs_compacted/"
echo "============================================"