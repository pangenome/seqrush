# How Seqwish Works: A Precise Algorithm Description

This document provides a detailed, precise description of how the seqwish algorithm constructs variation graphs from sequences and pairwise alignments.

## Overview

Seqwish takes a set of sequences and pairwise alignments (in PAF format) and produces an unbiased variation graph where each base position is represented by a node, with nodes unified based on alignment matches.

## Input

1. **Sequences** (FASTA format): A collection of sequences S = {s₁, s₂, ..., sₙ}
2. **Alignments** (PAF format): Pairwise alignments A including self-alignments

## Core Data Structures

1. **Position Encoding**: Each position is encoded as a 64-bit integer:
   - Bit 0: Orientation (0 = forward, 1 = reverse)
   - Bits 1-63: Offset in the global sequence space
   - Formula: `pos = (offset << 1) | orientation`

2. **Union-Find (Disjoint Set)**: Lock-free parallel implementation
   - Stores parent-rank pairs in 128-bit atomic values
   - Supports concurrent find() and unite() operations
   - Uses path compression and union by rank

3. **Graph Sequence**: A linearized representation of all union components

## Algorithm Steps

### Step 1: Process Input Sequences

1. Concatenate all input sequences with separator characters
2. Assign global offsets to each sequence
3. Create bidirectional position mappings (forward and reverse complement)

### Step 2: Process Alignments

For each alignment in the PAF file:

1. **Parse Alignment**:
   - Extract query/target positions, orientations, and CIGAR string
   - Handle relative/absolute strand orientations correctly

2. **Process CIGAR Operations**:
   ```
   For each CIGAR operation:
     If operation is 'M' (match):
       For each matched position i in [0, match_length):
         pos1 = encode_position(seq1_offset + seq1_pos + i, seq1_orientation)
         pos2 = encode_position(seq2_offset + seq2_pos + i, seq2_orientation)
         union_find.unite(pos1, pos2)
   ```

3. **Handle Reverse Complements**:
   - When seq2 is reverse complement:
     - Map positions from 3' to 5' coordinates
     - Unite forward positions with reverse complement positions

### Step 3: Build Graph Sequence (Critical Difference)

Instead of creating one node per union component, seqwish builds a "graph sequence":

1. **Chunked Processing** (for memory efficiency):
   ```
   For each chunk of size batch_size (default 1M bases):
     1. Find all alignments overlapping this chunk
     2. Apply union-find to positions in chunk
     3. Sort union components by minimum position
     4. Emit components to graph sequence in sorted order
   ```

2. **Graph Sequence Construction**:
   ```
   graph_sequence = []
   component_to_graph_pos = {}
   
   For each union component C in sorted order:
     graph_pos = len(graph_sequence)
     For each position p in C:
       component_to_graph_pos[find(p)] = graph_pos
     graph_sequence.append(C)
   ```

### Step 4: Identify Node Boundaries

Seqwish doesn't use union components as nodes directly. Instead:

1. **Mark Sequence Boundaries**:
   ```
   node_boundaries = BitVector(len(graph_sequence))
   
   For each input sequence s:
     For each position i in s:
       graph_pos = component_to_graph_pos[find(encode_pos(s, i))]
       if i == 0 or i == len(s)-1:
         node_boundaries[graph_pos] = 1
       if previous_graph_pos != graph_pos - 1:
         node_boundaries[graph_pos] = 1
   ```

2. **Create Nodes**:
   - Each run of positions between boundaries becomes a single node
   - This naturally merges linear regions into single nodes

### Step 5: Build Final Graph

1. **Create Nodes**:
   ```
   nodes = []
   current_node = []
   
   For each position i in graph_sequence:
     current_node.append(graph_sequence[i])
     if node_boundaries[i]:
       nodes.append(Node(current_node))
       current_node = []
   ```

2. **Create Paths**:
   ```
   For each input sequence s:
     path = []
     For each position i in s:
       node_id = position_to_node[component_to_graph_pos[find(encode_pos(s, i))]]
       if path.empty() or path.back() != node_id:
         path.append(node_id)
     paths[s.name] = path
   ```

3. **Create Edges**:
   ```
   For each path p:
     For each consecutive pair (u, v) in p:
       edges.add(Edge(u, v))
   ```

## Key Insights

1. **Transitive Closure**: Union-find ensures if A matches B and B matches C, then A, B, C are unified
2. **Graph Sequence**: The linearization of union components by minimum position is crucial
3. **Node Boundaries**: Determined by sequence traversal patterns, not union components
4. **Self-Alignments**: Establish the backbone structure of each sequence
5. **Chunked Processing**: Enables handling of large inputs with bounded memory

## Complexity

- **Time**: O(|A| log |A|) dominated by sorting operations
- **Space**: O(|S|) for graph sequence, reduced from O(|S|²) by match compression
- **Memory**: Bounded by largest transitive closure component

## Why This Produces Fewer Nodes

The key is that seqwish doesn't create one node per union component. Instead:
1. Union components are sorted and linearized into a graph sequence
2. Node boundaries are determined by how input sequences traverse this graph sequence
3. Linear regions without branching naturally merge into single nodes
4. Result: Much larger nodes that better represent the underlying sequence structure

This approach creates nodes based on the graph topology and sequence traversal patterns, not just on union-find components.