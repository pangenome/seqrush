#[cfg(test)]
mod tests {
    use seqrush::bidirected_union_find::BidirectedUnionFind;
    use seqrush::pos::make_pos;

    #[test]
    fn test_rc_alignment_debug() {
    // Test case: two sequences, second is RC aligned
    // Seq1: ACGT at offset 0
    // Seq2: ACGT at offset 10 (but aligned as RC)
    
    let uf = BidirectedUnionFind::new(100);
    
    println!("Testing RC alignment:");
    println!("Seq1: ACGT (offset 0)");
    println!("Seq2: ACGT (offset 10, aligned as RC)");
    println!();
    
    // When seq2 is RC aligned, seq1[i] matches seq2[len-1-i]
    // So: seq1[0]=A matches seq2[3]=T in reverse (which is A)
    
    // Simulate what unite_matching_region does for RC
    for i in 0..4 {
        let pos1_fwd = make_pos(0 + i, false);
        
        // RC calculation
        let rc_local_pos = i;
        let forward_local_pos = 4 - 1 - rc_local_pos;
        let seq2_global_offset = 10 + forward_local_pos;
        
        let pos2_fwd = make_pos(seq2_global_offset, false);
        let pos2_rev = make_pos(seq2_global_offset, true);
        
        println!("Position {}: seq1[{}] matches seq2[{}] in RC", i, i, forward_local_pos);
        println!("  pos1_fwd = {}", pos1_fwd);
        println!("  pos2_fwd = {} (seq2 position {} forward)", pos2_fwd, forward_local_pos);
        println!("  pos2_rev = {} (seq2 position {} reverse)", pos2_rev, forward_local_pos);
        
        // The current code does:
        uf.unite(pos1_fwd, pos2_rev);  // This makes sense
        uf.unite(pos1_fwd, pos2_fwd);  // This seems wrong!
        
        println!("  United: {} <-> {} (correct)", pos1_fwd, pos2_rev);
        println!("  United: {} <-> {} (suspicious!)", pos1_fwd, pos2_fwd);
        println!();
    }
    
    // Check what got united
    println!("Checking union results:");
    for i in 0..4 {
        let pos1 = make_pos(i, false);
        let rep1 = uf.find(pos1);
        println!("Position {} has representative {}", i, rep1);
        
        // Check what it's united with
        for j in 10..14 {
            let pos2_fwd = make_pos(j, false);
            let pos2_rev = make_pos(j, true);
            if uf.same(pos1, pos2_fwd) {
                println!("  -> united with position {} forward", j);
            }
            if uf.same(pos1, pos2_rev) {
                println!("  -> united with position {} reverse", j);
            }
        }
    }
    }
}