#!/bin/bash
# Script to build and visualize HLA graph using SeqRush and odgi

set -e  # Exit on error

# Default values
INPUT="${1:-HLA-zoo/seqs/B-3106.fa}"
OUTPUT_PREFIX="${2:-b}"
MIN_MATCH="${3:-0}"

echo "Building graph with SeqRush..."
cargo run --release --bin seqrush -- \
    -s "$INPUT" \
    -o "${OUTPUT_PREFIX}.gfa" \
    -k "$MIN_MATCH" \
    --no-compact

echo "Converting to odgi format..."
odgi build -g "${OUTPUT_PREFIX}.gfa" -o "${OUTPUT_PREFIX}.og"

echo "Sorting graph..."
odgi sort -i "${OUTPUT_PREFIX}.og" -o "${OUTPUT_PREFIX}.sorted.og"

echo "Compacting with odgi unchop..."
odgi unchop -i "${OUTPUT_PREFIX}.sorted.og" -o "${OUTPUT_PREFIX}.unchopped.og"

echo "Converting back to GFA..."
odgi view -i "${OUTPUT_PREFIX}.unchopped.og" -g > "${OUTPUT_PREFIX}.s.gfa"

echo "Creating visualization..."
odgi viz -i "${OUTPUT_PREFIX}.unchopped.og" -o "${OUTPUT_PREFIX}.s.png"

echo "Graph statistics:"
odgi stats -i "${OUTPUT_PREFIX}.unchopped.og" -S

echo "Visualization saved to ${OUTPUT_PREFIX}.s.png"
echo "To view: eog ${OUTPUT_PREFIX}.s.png"