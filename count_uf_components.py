#!/usr/bin/env python3
import sys

# Simple script to count unique union-find representatives
# by analyzing a GFA and counting unique node sequences

nodes = {}
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('S'):
            parts = line.strip().split('\t')
            node_id = int(parts[1])
            seq = parts[2]
            nodes[node_id] = seq

print(f"Total nodes: {len(nodes)}")
print(f"Unique sequences: {len(set(nodes.values()))}")
print(f"Single-base nodes: {sum(1 for seq in nodes.values() if len(seq) == 1)}")
print(f"Multi-base nodes: {sum(1 for seq in nodes.values() if len(seq) > 1)}")

# Show distribution of node lengths
from collections import Counter
lengths = Counter(len(seq) for seq in nodes.values())
print(f"\nNode length distribution:")
for length in sorted(lengths.keys())[:20]:  # Show first 20
    print(f"  Length {length}: {lengths[length]} nodes")
