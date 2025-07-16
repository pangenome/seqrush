# SeqRush Implementation Issues

This document describes precisely how the current SeqRush implementation differs from seqwish and why it produces incorrect results (too many nodes).

## Current SeqRush Implementation

### 1. Union-Find Implementation (CORRECT)

The union-find implementation in `src/bidirected_union_find.rs` is **fundamentally correct**:
- Properly handles bidirectional positions (forward/reverse)
- Correctly unites matching positions from alignments
- Implements transitive closure (if A→B and B→C, then A→C)
- All tests pass, confirming correct behavior

### 2. Alignment Processing (CORRECT)

The CIGAR parsing and position uniting logic is **correct**:
- Parses CIGAR strings properly
- Only unites positions when bases actually match
- Handles reverse complement alignments correctly
- Processes self-alignments to establish sequence backbones

### 3. Node Creation (INCORRECT - Core Issue)

**Current Implementation**:
```rust
// From src/seqrush.rs - build_initial_graph()
for seq in sequences {
    for position in seq {
        union_rep = union_find.find(position)
        if not seen(union_rep):
            create_new_node(union_rep)  // ONE NODE PER UNION COMPONENT
        add_node_to_path(node_id)
    }
}
```

**The Problem**: SeqRush creates **one node per union-find component**, resulting in:
- 2397 initial nodes (one per union component)
- 1785 nodes after compaction
- Seqwish: 471 nodes for the same input

### 4. Why This Is Wrong

The fundamental misunderstanding is treating union components as nodes directly. Consider this example:

**Input**: Three sequences with a shared region
```
Seq1: ATCGATCG
Seq2: ATCGATCG  
Seq3: ATCGATCG
```

**SeqRush Current Behavior**:
- Position 0: Union component {S1[0], S2[0], S3[0]} → Node 1
- Position 1: Union component {S1[1], S2[1], S3[1]} → Node 2
- Position 2: Union component {S1[2], S2[2], S3[2]} → Node 3
- ... and so on
- Result: 8 nodes (one per position)

**What Should Happen (Seqwish)**:
- All positions form a linear chain
- No branching points exist
- Result: 1 node containing "ATCGATCG"

### 5. The Missing Graph Sequence Step

SeqRush is missing the critical **graph sequence construction** step:

**What SeqRush Does**:
```
Alignments → Union-Find → Union Components → Nodes (1:1 mapping)
```

**What Seqwish Does**:
```
Alignments → Union-Find → Union Components → Graph Sequence → Node Boundaries → Nodes
```

The graph sequence step:
1. Sorts union components by minimum position
2. Creates a linear ordering of all components
3. Maps input sequences to this linear ordering
4. Identifies node boundaries based on sequence traversal patterns

### 6. Specific Implementation Errors

#### Error 1: Direct Component-to-Node Mapping
```rust
// WRONG: Creates one node per component
let node_id = match union_to_node.get(&union_rep) {
    Some(&id) => id,
    None => {
        let id = next_node_id;
        next_node_id += 1;
        union_to_node.insert(union_rep, id);
        // Creates a new node for EVERY union component
        graph.nodes.insert(id, Node { ... });
        id
    }
};
```

#### Error 2: No Graph Sequence Construction
SeqRush completely lacks the graph sequence construction phase that:
- Linearizes union components
- Establishes a total ordering
- Enables proper node boundary detection

#### Error 3: Inadequate Compaction
While SeqRush implements compaction, it operates on already-created nodes:
- Can only merge existing nodes
- Cannot create the large, multi-character nodes that seqwish produces
- Limited by the initial one-node-per-component creation

### 7. Why Compaction Isn't Enough

The compaction algorithm can only merge nodes that:
1. Appear consecutively in all paths
2. Have no branching between them

But with one node per union component:
- Most nodes represent single positions
- Limited merging opportunities
- Cannot achieve seqwish's level of node reduction

### 8. Concrete Example of the Problem

Given alignment: Seq1[0:8] matches Seq2[0:8]

**SeqRush**:
- Creates 8 union components (one per position)
- Creates 8 nodes initially
- Compaction might merge to 1 node (if no branching)

**Seqwish**:
- Creates 8 union components (same as SeqRush)
- Builds graph sequence: [C0, C1, C2, C3, C4, C5, C6, C7]
- Detects no boundaries within the match
- Creates 1 node directly

### 9. The Fix Required

To match seqwish's behavior, SeqRush needs to:

1. **Keep the union-find as-is** (it's correct)

2. **Add graph sequence construction**:
   ```rust
   struct GraphSequence {
       components: Vec<UnionComponent>,
       component_to_seq_pos: HashMap<Pos, usize>,
       seq_pos_to_component: Vec<Pos>,
   }
   ```

3. **Implement node boundary detection**:
   ```rust
   fn find_node_boundaries(graph_seq: &GraphSequence, sequences: &[Sequence]) -> BitVec {
       let mut boundaries = BitVec::new(graph_seq.len());
       for seq in sequences {
           for (i, pos) in seq.positions() {
               let seq_pos = graph_seq.component_to_seq_pos[union_find.find(pos)];
               // Mark boundaries at sequence starts/ends and discontinuities
               if i == 0 || i == seq.len()-1 || seq_pos != prev_seq_pos + 1 {
                   boundaries.set(seq_pos, true);
               }
           }
       }
       boundaries
   }
   ```

4. **Create nodes from boundaries**:
   ```rust
   fn create_nodes_from_boundaries(graph_seq: &GraphSequence, boundaries: &BitVec) -> Vec<Node> {
       let mut nodes = Vec::new();
       let mut current_node = Vec::new();
       
       for (i, component) in graph_seq.components.iter().enumerate() {
           current_node.extend(component.positions());
           if boundaries[i] {
               nodes.push(Node::from_positions(current_node));
               current_node = Vec::new();
           }
       }
       nodes
   }
   ```

## Summary

The SeqRush implementation is **mostly correct** - the union-find works properly and alignments are processed correctly. The fatal flaw is the **direct mapping from union components to nodes**, which creates far too many nodes. The solution requires implementing seqwish's graph sequence construction approach, which uses union components as an intermediate representation rather than the final node structure.