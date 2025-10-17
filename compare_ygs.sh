#!/bin/bash
set -e

# Compare ODGI vs seqrush Ygs sorting step-by-step

INPUT_FA="$1"
if [ -z "$INPUT_FA" ]; then
    echo "Usage: $0 <input.fa>"
    exit 1
fi

echo "=== Generating alignments ==="
sweepga --paf "$INPUT_FA" 2>/dev/null > test.paf
echo "Generated $(wc -l < test.paf) alignments"

echo ""
echo "=== Building initial graphs ==="

# ODGI: build base graph
seqwish -p test.paf -s "$INPUT_FA" -g test.base.gfa 2>&1 | grep -v "^\["
odgi build -g test.base.gfa -o test.base.og -P 2>&1 | tail -1

# Get stats
echo ""
echo "Base graph stats:"
odgi stats -i test.base.og -S 2>&1

echo ""
echo "=== Step 1: After SGD only ==="

# ODGI: Just Y (SGD) step
odgi sort -p Y -i test.base.og -o test.odgi.Y.og -P 2>&1 | grep "SGD"
echo "ODGI after Y:"
odgi stats -i test.odgi.Y.og -S 2>&1

# Seqrush: Build and get ordering after just SGD
echo ""
echo "Seqrush building graph..."
cargo run --release --bin seqrush -- -s "$INPUT_FA" -o test.seqrush.unsorted.gfa 2>&1 | grep -E "(nodes|edges|paths)" | head -3

echo ""
echo "=== Step 2: After SGD + Groom ==="

# ODGI: Y + g
odgi sort -p Yg -i test.base.og -o test.odgi.Yg.og -P 2>&1 | tail -1
echo "ODGI after Yg:"
odgi stats -i test.odgi.Yg.og -S 2>&1

echo ""
echo "=== Step 3: Full Ygs pipeline ==="

# ODGI: Full Ygs
odgi sort -p Ygs -i test.base.og -o test.odgi.Ygs.og -P 2>&1 | tail -1
echo "ODGI after Ygs:"
odgi stats -i test.odgi.Ygs.og -S 2>&1

# Seqrush: Full pipeline
cargo run --release --bin seqrush -- -s "$INPUT_FA" -o test.seqrush.Ygs.gfa --sgd-sort --verbose 2>&1 | grep -E "(Step|nodes|edges)" | head -10
odgi build -g test.seqrush.Ygs.gfa -o test.seqrush.Ygs.og -P 2>&1 | tail -1
echo "Seqrush after Ygs:"
odgi stats -i test.seqrush.Ygs.og -S 2>&1

echo ""
echo "=== Visualizations ==="
odgi viz -i test.odgi.Ygs.og -o test.odgi.Ygs.png 2>&1 | tail -1
odgi viz -i test.seqrush.Ygs.og -o test.seqrush.Ygs.png 2>&1 | tail -1
echo "Created:"
echo "  test.odgi.Ygs.png"
echo "  test.seqrush.Ygs.png"

echo ""
echo "=== Node order comparison (first 20 nodes) ==="
echo "ODGI node order:"
odgi view -i test.odgi.Ygs.og -g | grep "^S" | head -20 | cut -f2
echo ""
echo "Seqrush node order:"
grep "^S" test.seqrush.Ygs.gfa | head -20 | cut -f2
