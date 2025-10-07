#!/usr/bin/env python3
"""Analyze edges in GFA file to find duplicates and complements."""

import sys
from collections import defaultdict

def flip_orientation(orient):
    return '-' if orient == '+' else '+'

def parse_edge(line):
    parts = line.strip().split('\t')
    if len(parts) >= 5 and parts[0] == 'L':
        from_node = int(parts[1])
        from_orient = parts[2]
        to_node = int(parts[3])
        to_orient = parts[4]
        return (from_node, from_orient, to_node, to_orient)
    return None

def edge_complement(edge):
    """Get the complement of an edge in a bidirected graph.
    If edge is A+ -> B+, complement is B- -> A-"""
    from_node, from_orient, to_node, to_orient = edge
    return (to_node, flip_orientation(to_orient),
            from_node, flip_orientation(from_orient))

def main():
    if len(sys.argv) != 2:
        print("Usage: python analyze_edges.py <gfa_file>")
        sys.exit(1)

    edges = []
    with open(sys.argv[1]) as f:
        for line in f:
            if line.startswith('L'):
                edge = parse_edge(line)
                if edge:
                    edges.append(edge)

    print(f"Total edges: {len(edges)}")

    # Find unique edges (considering complements)
    unique_edges = set()
    complement_pairs = 0

    for edge in edges:
        complement = edge_complement(edge)
        # Use canonical form (smaller edge first)
        canonical = min(edge, complement)
        if canonical in unique_edges:
            complement_pairs += 1
        unique_edges.add(canonical)

    print(f"Unique edges (treating complements as same): {len(unique_edges)}")
    print(f"Complement pairs found: {complement_pairs}")

    # Check for duplicate edges (exact same edge appearing multiple times)
    edge_counts = defaultdict(int)
    for edge in edges:
        edge_counts[edge] += 1

    duplicates = {e: c for e, c in edge_counts.items() if c > 1}
    if duplicates:
        print(f"\nDuplicate edges found: {len(duplicates)}")
        for edge, count in list(duplicates.items())[:5]:
            print(f"  {edge}: appears {count} times")
    else:
        print("\nNo exact duplicates found")

    # Analyze edge patterns
    print("\nEdge pattern analysis:")
    both_directions = 0
    edge_set = set(edges)

    for edge in edges:
        from_node, from_orient, to_node, to_orient = edge
        # Check if reverse edge exists (same nodes, same orientations, reversed direction)
        reverse = (to_node, to_orient, from_node, from_orient)
        if reverse in edge_set:
            both_directions += 1

    print(f"Edges with both directions stored: {both_directions // 2}")

    # Sample some edges to see patterns
    print("\nFirst 10 edges:")
    for edge in edges[:10]:
        from_node, from_orient, to_node, to_orient = edge
        print(f"  {from_node}{from_orient} -> {to_node}{to_orient}")
        complement = edge_complement(edge)
        if complement in edge_set:
            print(f"    (complement {complement[0]}{complement[1]} -> {complement[2]}{complement[3]} exists)")

if __name__ == "__main__":
    main()