#!/usr/bin/env python3
import sys

def analyze_edge_spans(gfa_file):
    edges = []
    node_positions = {}
    
    # Read GFA and extract edges
    with open(gfa_file, 'r') as f:
        position = 0
        for line in f:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                node_id = int(parts[1])
                node_positions[node_id] = position
                position += 1
            elif line.startswith('L'):
                parts = line.strip().split('\t')
                from_node = int(parts[1])
                to_node = int(parts[3])
                edges.append((from_node, to_node))
    
    # Calculate edge spans
    edge_spans = []
    for from_node, to_node in edges:
        if from_node in node_positions and to_node in node_positions:
            span = abs(node_positions[to_node] - node_positions[from_node])
            edge_spans.append((span, from_node, to_node))
    
    # Sort by span length
    edge_spans.sort(reverse=True)
    
    # Print statistics
    print(f"Total edges: {len(edge_spans)}")
    print(f"Average span: {sum(s[0] for s in edge_spans) / len(edge_spans):.2f}")
    print(f"Max span: {edge_spans[0][0] if edge_spans else 0}")
    print("\nTop 20 longest edge spans:")
    for i, (span, from_node, to_node) in enumerate(edge_spans[:20]):
        print(f"{i+1}. Span={span}: {from_node} -> {to_node}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        analyze_edge_spans(sys.argv[1])
    else:
        analyze_edge_spans("b.gfa")