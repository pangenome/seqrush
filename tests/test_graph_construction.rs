use seqrush::bidirected_union_find::BidirectedUnionFind;
use seqrush::pos::make_pos;

#[test]
fn test_unite_matching_region_forward() {
    // Test basic forward-forward alignment
    let uf = BidirectedUnionFind::new(100);
    
    // Unite positions 10-14 with positions 20-24
    uf.unite_matching_region(
        10,  // seq1_offset
        20,  // seq2_offset
        0,   // seq1_local_start
        0,   // seq2_local_start
        5,   // match_length
        false, // seq2_is_rc
        100,   // seq2_len (unused for forward)
    );
    
    // Check that corresponding positions are united
    assert!(uf.same(make_pos(10, false), make_pos(20, false)));
    assert!(uf.same(make_pos(11, false), make_pos(21, false)));
    assert!(uf.same(make_pos(12, false), make_pos(22, false)));
    assert!(uf.same(make_pos(13, false), make_pos(23, false)));
    assert!(uf.same(make_pos(14, false), make_pos(24, false)));
}

#[test]
fn test_unite_matching_region_with_offset() {
    // Test alignment with local offsets
    let uf = BidirectedUnionFind::new(100);
    
    // Unite with local offsets: seq1[5-7] matches seq2[10-12]
    uf.unite_matching_region(
        0,   // seq1_offset
        50,  // seq2_offset  
        5,   // seq1_local_start
        10,  // seq2_local_start
        3,   // match_length
        false, // seq2_is_rc
        100,   // seq2_len
    );
    
    // Check correct positions are united
    assert!(uf.same(make_pos(5, false), make_pos(60, false)));  // 0+5 with 50+10
    assert!(uf.same(make_pos(6, false), make_pos(61, false)));
    assert!(uf.same(make_pos(7, false), make_pos(62, false)));
    
    // Other positions should not be united
    assert!(!uf.same(make_pos(4, false), make_pos(59, false)));
    assert!(!uf.same(make_pos(8, false), make_pos(63, false)));
}

#[test]
fn test_unite_matching_region_reverse() {
    // Test reverse complement alignment
    let uf = BidirectedUnionFind::new(100);
    
    // Unite with RC: seq1[0-3] matches seq2[0-3] in reverse
    uf.unite_matching_region_seq2_rc(
        0,   // seq1_offset
        10,  // seq2_offset
        0,   // seq1_local_start
        0,   // seq2_local_start (in RC coordinates)
        4,   // match_length
        true, // seq2_is_rc
        4,    // seq2_len
    );
    
    // When seq2 is RC, seq1[i] matches seq2[len-1-i] in reverse orientation
    // seq1[0] should match seq2[3]- (position 13 reverse)
    // seq1[1] should match seq2[2]- (position 12 reverse)
    // seq1[2] should match seq2[1]- (position 11 reverse)
    // seq1[3] should match seq2[0]- (position 10 reverse)
    assert!(uf.same(make_pos(0, false), make_pos(13, true)));
    assert!(uf.same(make_pos(1, false), make_pos(12, true)));
    assert!(uf.same(make_pos(2, false), make_pos(11, true)));
    assert!(uf.same(make_pos(3, false), make_pos(10, true)));
}

#[test]
fn test_transitive_closure() {
    // Test that transitive relationships are captured
    let uf = BidirectedUnionFind::new(100);
    
    // Unite A-B
    uf.unite_matching_region(0, 10, 0, 0, 1, false, 100);
    // Unite B-C
    uf.unite_matching_region(10, 20, 0, 0, 1, false, 100);
    
    // Check transitivity: A should be connected to C
    assert!(uf.same(make_pos(0, false), make_pos(10, false)));
    assert!(uf.same(make_pos(10, false), make_pos(20, false)));
    assert!(uf.same(make_pos(0, false), make_pos(20, false))); // Transitive
}

#[test]
fn test_self_alignment() {
    // Test self-alignment behavior
    let uf = BidirectedUnionFind::new(100);
    
    // Self-alignment: seq[0-3] aligns to itself
    uf.unite_matching_region(
        10,  // seq_offset
        10,  // same seq_offset
        0,   // local_start
        0,   // local_start
        4,   // match_length
        false, // not RC
        4,    // seq_len
    );
    
    // Each position is united with itself (no-op)
    // This doesn't create connections between consecutive positions
    assert!(uf.same(make_pos(10, false), make_pos(10, false)));
    assert!(uf.same(make_pos(11, false), make_pos(11, false)));
    
    // Consecutive positions are NOT united by self-alignment
    assert!(!uf.same(make_pos(10, false), make_pos(11, false)));
}

#[test]
fn test_orientation_handling() {
    // Test that forward and reverse orientations are handled correctly
    let uf = BidirectedUnionFind::new(100);
    
    // Unite position 10+ with position 20+
    uf.unite(make_pos(10, false), make_pos(20, false));
    
    // This should also unite 10- with 20- (implicit RC relationship)
    // Actually, let me check the implementation...
    // The union-find just tracks which positions are equivalent
    assert!(uf.same(make_pos(10, false), make_pos(20, false)));
    
    // Different orientations of same position are NOT automatically united
    assert!(!uf.same(make_pos(10, false), make_pos(10, true)));
}

#[test]
fn test_multiple_sequences() {
    // Test alignment of multiple sequences
    let uf = BidirectedUnionFind::new(300);
    
    // Three sequences at offsets 0, 100, 200
    // Align seq1 to seq2
    uf.unite_matching_region(0, 100, 0, 0, 5, false, 100);
    // Align seq2 to seq3
    uf.unite_matching_region(100, 200, 0, 0, 5, false, 100);
    
    // Check transitivity across all three sequences
    for i in 0..5 {
        assert!(uf.same(make_pos(i, false), make_pos(100 + i, false)));
        assert!(uf.same(make_pos(100 + i, false), make_pos(200 + i, false)));
        assert!(uf.same(make_pos(i, false), make_pos(200 + i, false))); // Transitive
    }
}

#[test]
fn test_union_find_behavior() {
    // Test basic union-find behavior
    let uf = BidirectedUnionFind::new(100);
    
    // Initially, each position is its own component
    assert_ne!(uf.find(make_pos(0, false)), uf.find(make_pos(1, false)));
    
    // After uniting, they should have the same representative
    uf.unite(make_pos(0, false), make_pos(1, false));
    assert_eq!(uf.find(make_pos(0, false)), uf.find(make_pos(1, false)));
    
    // Unite with a third position
    uf.unite(make_pos(1, false), make_pos(2, false));
    
    // All three should now be connected
    assert_eq!(uf.find(make_pos(0, false)), uf.find(make_pos(2, false)));
}

#[test]
fn test_real_world_scenario() {
    // Test a realistic scenario with SNP
    let uf = BidirectedUnionFind::new(100);
    
    // Three sequences: ref, var1 (with SNP), var2 (identical to ref)
    // ref:  ATCGATCGATCG (offset 0)
    // var1: ATCGATGGATCG (offset 20, SNP at position 6: C->G)
    // var2: ATCGATCGATCG (offset 40)
    
    // Align ref to var1: 6=1X5=
    uf.unite_matching_region(0, 20, 0, 0, 6, false, 12);  // First 6 match
    // Skip position 6 (mismatch)
    uf.unite_matching_region(0, 20, 7, 7, 5, false, 12);  // Last 5 match
    
    // Align ref to var2: 12=
    uf.unite_matching_region(0, 40, 0, 0, 12, false, 12);
    
    // Align var1 to var2: 6=1X5=
    uf.unite_matching_region(20, 40, 0, 0, 6, false, 12);
    uf.unite_matching_region(20, 40, 7, 7, 5, false, 12);
    
    // Check results
    // Positions 0-5 should all be united across all three
    for i in 0..6 {
        assert!(uf.same(make_pos(i, false), make_pos(20 + i, false)));
        assert!(uf.same(make_pos(i, false), make_pos(40 + i, false)));
    }
    
    // Position 6: ref and var2 united, var1 separate
    assert!(uf.same(make_pos(6, false), make_pos(46, false)));
    assert!(!uf.same(make_pos(6, false), make_pos(26, false)));
    
    // Positions 7-11 all united
    for i in 7..12 {
        assert!(uf.same(make_pos(i, false), make_pos(20 + i, false)));
        assert!(uf.same(make_pos(i, false), make_pos(40 + i, false)));
    }
}