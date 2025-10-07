use seqrush::seqrush::{SeqRush, Sequence};
use seqrush::pos::make_pos;

fn main() {
    // Create a simple case that will produce duplicate edges
    // Two sequences with a path that creates both forward and reverse edges
    
    // Create sequences that will produce the problematic pattern
    let seq1 = Sequence {
        id: "seq1".to_string(),
        data: b"ABCDE".to_vec(),
        offset: 0,
    };
    
    let seq2 = Sequence {
        id: "seq2".to_string(),
        data: b"ABCDE".to_vec(), 
        offset: 5,
    };
    
    let mut seqrush = SeqRush::new(vec![seq1, seq2], 0);
    
    // Unite some positions to create shared nodes
    // Make it so node 2 and 3 are shared between sequences
    seqrush.union_find.unite(make_pos(1, false), make_pos(6, false)); // B positions
    seqrush.union_find.unite(make_pos(2, false), make_pos(7, false)); // C positions
    
    // Create a more complex uniting that will cause a path to visit a node multiple times
    // Unite position 3 with position 1 (both D and B map to same node)
    seqrush.union_find.unite(make_pos(3, false), make_pos(1, false));
    
    // Build graph
    println!("=== BUILDING GRAPH ===");
    let graph = seqrush.build_bidirected_graph(false).unwrap();
    
    // Print paths
    println!("\n=== PATHS ===");
    for path in &graph.paths {
        print!("{}: ", path.name);
        for (i, handle) in path.steps.iter().enumerate() {
            if i > 0 { print!(","); }
            print!("{}{}", handle.node_id(), if handle.is_reverse() { "-" } else { "+" });
        }
        println!();
    }
    
    // Print edges
    println!("\n=== EDGES (RAW) ===");
    let mut edges: Vec<_> = graph.edges.iter().collect();
    edges.sort_by_key(|e| (e.from.node_id(), e.from.is_reverse(), e.to.node_id(), e.to.is_reverse()));
    
    for edge in &edges {
        println!("L\t{}\t{}\t{}\t{}\t0M", 
            edge.from.node_id(),
            if edge.from.is_reverse() { "-" } else { "+" },
            edge.to.node_id(), 
            if edge.to.is_reverse() { "-" } else { "+" }
        );
    }
    
    // Show what happens with canonicalization
    println!("\n=== EDGES (CANONICAL FORM) ===");
    for edge in &edges {
        let canonical = edge.canonical();
        println!("L\t{}\t{}\t{}\t{}\t0M (original: {}{} -> {}{})", 
            canonical.from.node_id(),
            if canonical.from.is_reverse() { "-" } else { "+" },
            canonical.to.node_id(), 
            if canonical.to.is_reverse() { "-" } else { "+" },
            edge.from.node_id(),
            if edge.from.is_reverse() { "-" } else { "+" },
            edge.to.node_id(),
            if edge.to.is_reverse() { "-" } else { "+" }
        );
    }
    
    // Check for duplicate edge patterns
    println!("\n=== CHECKING FOR DUPLICATE PATTERNS ===");
    for edge1 in &edges {
        // The reverse complement of edge A+ -> B+ is B- -> A-
        // So if we have both A+ -> B+ and A- -> B-, that's a duplicate pattern
        let rc_edge_exists = edges.iter().any(|edge2| {
            edge1.from.node_id() == edge2.from.node_id() &&
            edge1.to.node_id() == edge2.to.node_id() &&
            edge1.from.is_reverse() != edge2.from.is_reverse() &&
            edge1.to.is_reverse() != edge2.to.is_reverse()
        });
        
        if rc_edge_exists && !edge1.from.is_reverse() {
            println!("Duplicate pattern found: {}+ -> {}+ AND {}- -> {}-",
                edge1.from.node_id(), edge1.to.node_id(),
                edge1.from.node_id(), edge1.to.node_id());
        }
    }
}