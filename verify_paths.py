#!/usr/bin/env python3
import sys

def parse_gfa(filename):
    nodes = {}
    paths = {}
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                node_id = parts[1]
                sequence = parts[2]
                nodes[node_id] = sequence
            elif line.startswith('P'):
                parts = line.strip().split('\t')
                path_name = parts[1]
                steps = parts[2].split(',')
                paths[path_name] = steps
    
    return nodes, paths

def reconstruct_sequence(nodes, steps):
    sequence = []
    for step in steps:
        node_id = step[:-1]  # Remove orientation
        orientation = step[-1]
        if node_id in nodes:
            if orientation == '+':
                sequence.append(nodes[node_id])
            else:  # orientation == '-'
                # Reverse complement
                base = nodes[node_id]
                rc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
                sequence.append(rc.get(base, base))
    return ''.join(sequence)

def main():
    if len(sys.argv) != 2:
        print("Usage: verify_paths.py <gfa_file>")
        sys.exit(1)
    
    gfa_file = sys.argv[1]
    nodes, paths = parse_gfa(gfa_file)
    
    print(f"Loaded {len(nodes)} nodes and {len(paths)} paths")
    
    # Show first few paths
    for i, (path_name, steps) in enumerate(paths.items()):
        if i >= 3:
            break
        seq = reconstruct_sequence(nodes, steps)
        print(f"\nPath: {path_name}")
        print(f"Length: {len(seq)}")
        print(f"First 100bp: {seq[:100]}")
        print(f"Last 100bp: {seq[-100:]}")

if __name__ == "__main__":
    main()