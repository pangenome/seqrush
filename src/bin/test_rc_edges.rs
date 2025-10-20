use seqrush::bidirected_union_find::BidirectedUnionFind;
use seqrush::pos::{make_pos};
use seqrush::seqrush::{SeqRush, Sequence};

fn main() {
    // Create two sequences that are reverse complements
    // seq1: ATCG
    // seq2: CGAT (reverse complement)
    
    let seq1 = Sequence {
        id: "seq1".to_string(),
        data: b"ATCG".to_vec(),
        offset: 0,
    };
    
    let seq2 = Sequence {
        id: "seq2".to_string(), 
        data: b"CGAT".to_vec(),
        offset: 4,
    };
    
    let mut seqrush = SeqRush::new(vec![seq1, seq2]);
    
    // Unite them via RC alignment
    // When seq1 is RC'd, it becomes CGAT which matches seq2
    println!("=== UNITING VIA RC ===");
    println!("seq1: ATCG at offset 0");
    println!("seq2: CGAT at offset 4");
    println!("When seq1 is RC'd: CGAT");
    println!("Calling unite_matching_region(0, 4, 0, 0, 4, true, 4)");
    seqrush.union_find.unite_matching_region(0, 4, 0, 0, 4, true, 4);
    
    // Manual verification
    println!("\nManual RC transformation:");
    for i in 0..4 {
        let rc_pos = i;
        let fwd_pos = 4 - 1 - rc_pos;
        println!("  RC[{}] -> forward[{}]", rc_pos, fwd_pos);
    }
    
    // Build graph
    println!("\n=== BUILDING GRAPH ===");
    let graph = seqrush.build_bidirected_graph(true).unwrap();
    
    // Debug: Print what's in union_to_node after building
    println!("\n=== FINAL NODE MAPPING ===");
    println!("Created {} nodes from unions", graph.nodes.len());
    
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
    println!("\n=== EDGES ===");
    let mut edges: Vec<_> = graph.edges.iter().collect();
    edges.sort_by_key(|e| (e.from.node_id(), e.from.is_reverse(), e.to.node_id(), e.to.is_reverse()));
    
    for edge in edges {
        println!("L\t{}\t{}\t{}\t{}\t0M", 
            edge.from.node_id(),
            if edge.from.is_reverse() { "-" } else { "+" },
            edge.to.node_id(), 
            if edge.to.is_reverse() { "-" } else { "+" }
        );
    }
    
    // Check what unions were created
    println!("\n=== UNION ANALYSIS ===");
    for seq_idx in 0..2 {
        let offset = seq_idx * 4;
        println!("Sequence {} (offset {}):", seq_idx + 1, offset);
        for pos in 0..4 {
            let global_pos = offset + pos;
            let fwd = make_pos(global_pos, false);
            let rev = make_pos(global_pos, true);
            let union_fwd = seqrush.union_find.find(fwd);
            let union_rev = seqrush.union_find.find(rev);
            println!("  pos {} (global {}): fwd union={}, rev union={}", 
                pos, global_pos, union_fwd, union_rev);
        }
    }
    
    // Check which unions are connected
    println!("\n=== UNION CONNECTIONS ===");
    println!("seq1 pos 0 fwd (union 0) same as seq2 pos 3 rev (union 14)? {}", 
        seqrush.union_find.same(make_pos(0, false), make_pos(7, true)));
    println!("seq1 pos 3 fwd (union 6) same as seq2 pos 0 rev (union 8)? {}", 
        seqrush.union_find.same(make_pos(3, false), make_pos(4, true)));
    
    // More detailed check
    println!("\nExpected connections from RC alignment:");
    println!("  seq1[0] fwd should connect to seq2[3] fwd (via RC)");
    println!("  seq1[1] fwd should connect to seq2[2] fwd (via RC)");
    println!("  seq1[2] fwd should connect to seq2[1] fwd (via RC)");
    println!("  seq1[3] fwd should connect to seq2[0] fwd (via RC)");
    
    println!("\nActual union connections:");
    for i in 0..4 {
        let seq1_pos = make_pos(i, false);
        let seq2_pos = make_pos(7-i, false); // RC position mapping
        println!("  seq1[{}] fwd (union {}) same as seq2[{}] fwd (union {})? {}",
            i, seqrush.union_find.find(seq1_pos),
            7-i, seqrush.union_find.find(seq2_pos),
            seqrush.union_find.same(seq1_pos, seq2_pos));
    }
    
    // Check what the unite_matching_region actually did
    println!("\nWhat unite_matching_region should have done:");
    println!("  Unite seq1[3] rev with seq2[0] fwd");
    println!("  Unite seq1[2] rev with seq2[1] fwd");
    println!("  Unite seq1[1] rev with seq2[2] fwd");
    println!("  Unite seq1[0] rev with seq2[3] fwd");
    
    println!("\nChecking those specific connections:");
    println!("  seq1[3] rev (union {}) same as seq2[0] fwd (union {})? {}",
        seqrush.union_find.find(make_pos(3, true)),
        seqrush.union_find.find(make_pos(4, false)),
        seqrush.union_find.same(make_pos(3, true), make_pos(4, false)));
    println!("  seq1[0] rev (union {}) same as seq2[3] fwd (union {})? {}",
        seqrush.union_find.find(make_pos(0, true)),
        seqrush.union_find.find(make_pos(7, false)),
        seqrush.union_find.same(make_pos(0, true), make_pos(7, false)));
}