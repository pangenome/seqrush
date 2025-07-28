#[cfg(test)]
mod tests {
    use seqrush::bidirected_union_find::BidirectedUnionFind;
    use seqrush::pos::make_pos;

    #[test] 
    fn test_seqrush_transitivity() {
        // Test that if we unite A-B and B-C, then A and C are connected
        let uf = BidirectedUnionFind::new(100);
        
        // Unite position 0 with position 10
        uf.unite(make_pos(0, false), make_pos(10, false));
        
        // Unite position 10 with position 20  
        uf.unite(make_pos(10, false), make_pos(20, false));
        
        // Check transitivity
        assert!(uf.same(make_pos(0, false), make_pos(20, false)));
        
        // Check that find gives same representative
        assert_eq!(uf.find(make_pos(0, false)), uf.find(make_pos(20, false)));
    }
    
    #[test]
    fn test_multiple_sequences_unite() {
        // Simulate what happens with multiple identical sequences
        let uf = BidirectedUnionFind::new(1000);
        
        // Seq1 at offset 0, seq2 at offset 100, seq3 at offset 200
        // All have same content "ACGT"
        
        // Unite seq1[0] with seq2[0]
        uf.unite(make_pos(0, false), make_pos(100, false));
        
        // Unite seq2[0] with seq3[0]
        uf.unite(make_pos(100, false), make_pos(200, false));
        
        // All should be connected
        assert!(uf.same(make_pos(0, false), make_pos(200, false)));
        
        // Now if we have 9 sequences all with 'A' at position 0
        let offsets = vec![0, 100, 200, 300, 400, 500, 600, 700, 800];
        
        // Unite them all in a chain
        for i in 0..offsets.len()-1 {
            uf.unite(
                make_pos(offsets[i], false),
                make_pos(offsets[i+1], false)
            );
        }
        
        // All should be connected
        let rep0 = uf.find(make_pos(0, false));
        for &offset in &offsets {
            assert_eq!(uf.find(make_pos(offset, false)), rep0);
        }
    }
}