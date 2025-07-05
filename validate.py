#!/usr/bin/env python3
import sys

def parse_gfa(filename):
    nodes = {}
    paths = {}
    
    with open(filename) as f:
        for line in f:
            if line.startswith('S\t'):
                parts = line.strip().split('\t')
                node_id = int(parts[1])
                sequence = parts[2]
                nodes[node_id] = sequence
            elif line.startswith('P\t'):
                parts = line.strip().split('\t')
                seq_id = parts[1]
                path = parts[2]
                paths[seq_id] = path
    
    return nodes, paths

def reconstruct_sequence(nodes, path_str):
    sequence = []
    for node_ref in path_str.split(','):
        if node_ref.endswith('+'):
            node_id = int(node_ref[:-1])
            sequence.append(nodes[node_id])
    return ''.join(sequence)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: validate.py <gfa_file>")
        sys.exit(1)
    
    nodes, paths = parse_gfa(sys.argv[1])
    
    print(f"Found {len(nodes)} nodes and {len(paths)} paths")
    print("\nNodes:")
    for node_id in sorted(nodes.keys()):
        print(f"  {node_id}: {nodes[node_id]}")
    
    print("\nPaths and reconstructed sequences:")
    for seq_id, path in paths.items():
        reconstructed = reconstruct_sequence(nodes, path)
        print(f"  {seq_id}: {reconstructed}")
        print(f"    Path: {path}")