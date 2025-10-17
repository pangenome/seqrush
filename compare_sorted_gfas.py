#!/usr/bin/env python3
"""
Compare two sorted GFAs by matching nodes via content/context, then comparing ordering.
Nodes are matched by: sequence, length, neighbors, and path membership.
"""

import sys
import hashlib
from collections import defaultdict

def parse_gfa(filename):
    """Parse GFA and extract node/edge/path information"""
    nodes = {}  # node_id -> {seq, length, checksum}
    edges = defaultdict(set)  # node_id -> set of neighbor node_ids
    paths = []  # list of (path_name, [node_steps])

    with open(filename) as f:
        for line in f:
            if line.startswith('S\t'):
                parts = line.strip().split('\t')
                node_id = parts[1]
                seq = parts[2]
                nodes[node_id] = {
                    'seq': seq,
                    'length': len(seq),
                    'checksum': hashlib.md5(seq.encode()).hexdigest()[:8]
                }
            elif line.startswith('L\t'):
                parts = line.strip().split('\t')
                from_id, to_id = parts[1], parts[3]
                edges[from_id].add(to_id)
                edges[to_id].add(from_id)
            elif line.startswith('P\t'):
                parts = line.strip().split('\t')
                path_name = parts[1]
                steps = parts[2].split(',')
                # Extract just node IDs (remove orientation)
                node_ids = [step[:-1] for step in steps]
                paths.append((path_name, node_ids))

    return nodes, edges, paths

def build_node_signature(node_id, nodes, edges, paths):
    """Build a signature for a node based on content and context"""
    node = nodes[node_id]

    # Sequence-based identity
    seq_sig = f"{node['length']}:{node['checksum']}"

    # Neighbor context (sorted for consistency)
    neighbors = sorted(edges.get(node_id, set()))
    neighbor_sigs = []
    for n in neighbors:
        if n in nodes:
            neighbor_sigs.append(f"{nodes[n]['length']}:{nodes[n]['checksum'][:4]}")
    neighbor_sig = '|'.join(neighbor_sigs[:5])  # First 5 neighbors

    # Path membership
    in_paths = set()
    for path_name, node_list in paths:
        if node_id in node_list:
            in_paths.add(path_name)
    path_sig = ','.join(sorted(in_paths))

    return f"{seq_sig}#{neighbor_sig}#{path_sig}"

def match_nodes(nodes1, edges1, paths1, nodes2, edges2, paths2):
    """Match nodes from graph1 to graph2 based on signatures"""
    # Build signatures for all nodes
    sigs1 = {nid: build_node_signature(nid, nodes1, edges1, paths1) for nid in nodes1}
    sigs2 = {nid: build_node_signature(nid, nodes2, edges2, paths2) for nid in nodes2}

    # Create reverse mapping: signature -> node_id
    sig_to_node2 = {}
    for nid, sig in sigs2.items():
        if sig in sig_to_node2:
            # Multiple nodes with same signature - add to list
            if isinstance(sig_to_node2[sig], list):
                sig_to_node2[sig].append(nid)
            else:
                sig_to_node2[sig] = [sig_to_node2[sig], nid]
        else:
            sig_to_node2[sig] = nid

    # Match nodes from graph1 to graph2
    matches = {}  # node_id1 -> node_id2
    ambiguous = []
    unmatched = []

    for nid1, sig in sigs1.items():
        if sig in sig_to_node2:
            match = sig_to_node2[sig]
            if isinstance(match, list):
                ambiguous.append((nid1, match))
            else:
                matches[nid1] = match
        else:
            unmatched.append(nid1)

    return matches, ambiguous, unmatched

def compare_path_order(path1, path2, node_mapping):
    """Compare ordering of nodes in paths"""
    # Map path1 nodes to path2 node IDs
    mapped_path1 = []
    for nid in path1:
        if nid in node_mapping:
            mapped_path1.append(node_mapping[nid])
        else:
            mapped_path1.append(f"UNMAPPED_{nid}")

    # Find differences in ordering
    diffs = []
    max_check = min(len(mapped_path1), len(path2))

    for i in range(max_check):
        if mapped_path1[i] != path2[i]:
            diffs.append((i, mapped_path1[i], path2[i]))

    return diffs

def main():
    if len(sys.argv) != 3:
        print("Usage: compare_sorted_gfas.py <gfa1> <gfa2>")
        print("  Compares ordering of two GFA files by matching nodes via content")
        sys.exit(1)

    gfa1, gfa2 = sys.argv[1], sys.argv[2]

    print(f"Comparing {gfa1} vs {gfa2}")
    print("=" * 70)

    # Parse both GFAs
    nodes1, edges1, paths1 = parse_gfa(gfa1)
    nodes2, edges2, paths2 = parse_gfa(gfa2)

    print(f"Graph 1: {len(nodes1)} nodes, {sum(len(e) for e in edges1.values())//2} edges, {len(paths1)} paths")
    print(f"Graph 2: {len(nodes2)} nodes, {sum(len(e) for e in edges2.values())//2} edges, {len(paths2)} paths")
    print()

    # Match nodes
    matches, ambiguous, unmatched = match_nodes(nodes1, edges1, paths1, nodes2, edges2, paths2)

    print(f"Node matching:")
    print(f"  Matched: {len(matches)}/{len(nodes1)} nodes")
    print(f"  Ambiguous: {len(ambiguous)} nodes")
    print(f"  Unmatched: {len(unmatched)} nodes")
    print()

    if unmatched:
        print(f"Unmatched nodes from graph1: {unmatched[:10]}")
        if len(unmatched) > 10:
            print(f"  ... and {len(unmatched) - 10} more")
        print()

    # Compare path orderings
    print("Path ordering comparison:")
    for (name1, path1_nodes), (name2, path2_nodes) in zip(paths1, paths2):
        diffs = compare_path_order(path1_nodes, path2_nodes, matches)

        if not diffs:
            print(f"  ✓ {name1}: IDENTICAL ordering ({len(path1_nodes)} nodes)")
        else:
            print(f"  ✗ {name1}: {len(diffs)} differences in first {min(len(path1_nodes), len(path2_nodes))} nodes")
            for i, (pos, id1, id2) in enumerate(diffs[:5]):
                print(f"      Position {pos}: {id1} vs {id2}")
            if len(diffs) > 5:
                print(f"      ... and {len(diffs) - 5} more differences")

    print()
    print("=" * 70)

if __name__ == '__main__':
    main()
