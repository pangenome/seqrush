use seqrush::pos::make_pos;
use seqrush::seqrush::{SeqRush, Sequence};
use std::collections::HashSet;

#[test]
fn test_forward_reverse_unification() {
    // Test that forward and reverse orientations of each position are unified
    let seq = Sequence {
        id: "test".to_string(),
        data: b"ACGT".to_vec(),
        offset: 0,
    };

    let seqrush = SeqRush::new(vec![seq], 0);

    // Check that forward and reverse of each position have the same representative
    for i in 0..4 {
        let pos_fwd = make_pos(i, false);
        let pos_rev = make_pos(i, true);

        let rep_fwd = seqrush.union_find.find(pos_fwd);
        let rep_rev = seqrush.union_find.find(pos_rev);

        assert_eq!(
            rep_fwd, rep_rev,
            "Forward and reverse orientations of position {} should have same representative",
            i
        );
    }
}

#[test]
fn test_transitive_closure() {
    // Test that union-find implements transitive closure correctly
    // If A aligns to B and B aligns to C, then A should be united with C

    let seq1 = Sequence {
        id: "seq1".to_string(),
        data: b"ACGT".to_vec(),
        offset: 0,
    };

    let seq2 = Sequence {
        id: "seq2".to_string(),
        data: b"ACGT".to_vec(),
        offset: 4,
    };

    let seq3 = Sequence {
        id: "seq3".to_string(),
        data: b"ACGT".to_vec(),
        offset: 8,
    };

    let seqrush = SeqRush::new(vec![seq1, seq2, seq3], 0);

    // Unite seq1[0] with seq2[0]
    seqrush
        .union_find
        .unite_matching_region(0, 4, 0, 0, 4, false, 4);

    // Unite seq2[0] with seq3[0]
    seqrush
        .union_find
        .unite_matching_region(4, 8, 0, 0, 4, false, 4);

    // Check transitive closure: seq1[0] should be united with seq3[0]
    let pos1 = make_pos(0, false);
    let pos3 = make_pos(8, false);

    assert!(
        seqrush.union_find.same(pos1, pos3),
        "Transitive closure failed: seq1[0] should be united with seq3[0]"
    );
}

#[test]
fn test_single_component_per_position() {
    // Test that each position belongs to exactly one component
    let sequences = vec![
        Sequence {
            id: "seq1".to_string(),
            data: b"ACGTACGT".to_vec(),
            offset: 0,
        },
        Sequence {
            id: "seq2".to_string(),
            data: b"ACGTACGT".to_vec(),
            offset: 8,
        },
    ];

    let seqrush = SeqRush::new(sequences, 0);

    // Unite matching positions
    seqrush
        .union_find
        .unite_matching_region(0, 8, 0, 0, 8, false, 8);

    // Collect all unique representatives
    let mut representatives = HashSet::new();
    for i in 0..16 {
        let pos = make_pos(i, false);
        let rep = seqrush.union_find.find(pos);
        representatives.insert(rep);
    }

    // Should have exactly 8 components (each position from both sequences united)
    assert_eq!(
        representatives.len(),
        8,
        "Expected 8 components, got {}",
        representatives.len()
    );
}

#[test]
fn test_reverse_complement_alignment() {
    // Test that RC alignments are handled correctly
    // seq1: ACGT
    // seq2: ACGT (reverse complement of seq1)

    let seq1 = Sequence {
        id: "seq1".to_string(),
        data: b"ACGT".to_vec(),
        offset: 0,
    };

    let seq2 = Sequence {
        id: "seq2".to_string(),
        data: b"ACGT".to_vec(),
        offset: 4,
    };

    let seqrush = SeqRush::new(vec![seq1, seq2], 0);

    // seq1 aligns to seq2 via reverse complement
    // seq1[0..4] RC aligns to seq2[0..4]
    seqrush
        .union_find
        .unite_matching_region(0, 4, 0, 0, 4, true, 4);

    // Check the alignment:
    // When seq1 is reverse complemented, the alignment works as follows:
    // seq1 RC coordinates: 0 1 2 3 map to forward coordinates 3 2 1 0
    // So seq1 RC[0] = seq1 forward[3] = T
    // seq1 RC[0] aligns with seq2[0] = A (complement match)

    // Due to RC transformation in unite_matching_region:
    // seq1 RC pos 0 = forward pos 3, so seq1[3] rev should be united with seq2[0] fwd
    let seq1_3_rev = make_pos(3, true);
    let seq2_0_fwd = make_pos(4, false);

    assert!(
        seqrush.union_find.same(seq1_3_rev, seq2_0_fwd),
        "seq1[3] reverse should be united with seq2[0] forward"
    );

    // Similarly, seq1 RC pos 3 = forward pos 0
    let seq1_0_rev = make_pos(0, true);
    let seq2_3_fwd = make_pos(7, false);
    assert!(
        seqrush.union_find.same(seq1_0_rev, seq2_3_fwd),
        "seq1[0] reverse should be united with seq2[3] forward"
    );
}

#[test]
fn test_identical_sequences_produce_minimal_components() {
    // Test that N identical sequences produce exactly L components (where L is sequence length)
    let sequences = vec![
        Sequence {
            id: "seq1".to_string(),
            data: b"ACGTACGT".to_vec(),
            offset: 0,
        },
        Sequence {
            id: "seq2".to_string(),
            data: b"ACGTACGT".to_vec(),
            offset: 8,
        },
        Sequence {
            id: "seq3".to_string(),
            data: b"ACGTACGT".to_vec(),
            offset: 16,
        },
    ];

    let seqrush = SeqRush::new(sequences, 0);

    // Unite all matching positions
    seqrush
        .union_find
        .unite_matching_region(0, 8, 0, 0, 8, false, 8);
    seqrush
        .union_find
        .unite_matching_region(0, 16, 0, 0, 8, false, 8);
    seqrush
        .union_find
        .unite_matching_region(8, 16, 0, 0, 8, false, 8);

    // Count unique components
    let mut components = HashSet::new();
    for i in 0..24 {
        let pos = make_pos(i, false);
        let rep = seqrush.union_find.find(pos);
        components.insert(rep);
    }

    assert_eq!(
        components.len(),
        8,
        "Three identical sequences of length 8 should produce exactly 8 components"
    );
}

#[test]
fn test_no_false_unifications() {
    // Test that positions that don't align are not unified
    let sequences = vec![
        Sequence {
            id: "seq1".to_string(),
            data: b"AAAAAAAA".to_vec(),
            offset: 0,
        },
        Sequence {
            id: "seq2".to_string(),
            data: b"TTTTTTTT".to_vec(),
            offset: 8,
        },
    ];

    let seqrush = SeqRush::new(sequences, 0);

    // Without any alignments, positions from different sequences should not be united
    let pos1 = make_pos(0, false);
    let pos2 = make_pos(8, false);

    assert!(
        !seqrush.union_find.same(pos1, pos2),
        "Positions from different sequences should not be united without alignments"
    );

    // But forward and reverse of same position should still be united
    let pos1_rev = make_pos(0, true);
    assert!(
        seqrush.union_find.same(pos1, pos1_rev),
        "Forward and reverse of same position should always be united"
    );
}

#[test]
fn test_partial_alignment() {
    // Test that partial alignments only unite the matching region
    let seq1 = Sequence {
        id: "seq1".to_string(),
        data: b"AAAAGGGGTTTT".to_vec(),
        offset: 0,
    };

    let seq2 = Sequence {
        id: "seq2".to_string(),
        data: b"CCCCGGGGAAAA".to_vec(),
        offset: 12,
    };

    let seqrush = SeqRush::new(vec![seq1, seq2], 0);

    // Only the GGGG region matches (positions 4-7 in seq1, 4-7 in seq2)
    seqrush
        .union_find
        .unite_matching_region(0, 12, 4, 4, 4, false, 12);

    // Check that only the GGGG positions are united
    for i in 0..4 {
        let seq1_pos = make_pos(4 + i, false);
        let seq2_pos = make_pos(16 + i, false);
        assert!(
            seqrush.union_find.same(seq1_pos, seq2_pos),
            "Position {} in GGGG region should be united",
            i
        );
    }

    // Check that other positions are not united
    let seq1_start = make_pos(0, false);
    let seq2_start = make_pos(12, false);
    assert!(
        !seqrush.union_find.same(seq1_start, seq2_start),
        "Non-matching positions should not be united"
    );
}
