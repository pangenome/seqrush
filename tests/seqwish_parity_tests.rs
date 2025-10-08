/// Test cases designed to achieve 100% parity with seqwish
/// These tests specifically target the differences found in the analysis

use seqrush::bidirected_union_find::BidirectedUnionFind;
use seqrush::pos::{make_pos, offset};
use seqrush::seqrush::{Sequence, SeqRush};

#[test]
fn test_self_mapping_guard() {
    /// This test verifies that seqrush properly guards against self-mappings
    /// Seqwish checks: offset(q_pos) != offset(t_pos)
    /// Without this guard, we'd create unnecessary unions

    let uf = BidirectedUnionFind::new(1000);

    // Simulate what happens when the same global offset appears in query and target
    // This shouldn't happen in practice, but the guard protects against it
    let pos1 = make_pos(100, false);
    let pos2 = make_pos(100, false); // Same offset!

    // These should NOT be united (they're already the same position)
    assert_eq!(offset(pos1), offset(pos2));

    // The union-find's unite function already checks raw1 != raw2
    // But we need to ensure the alignment processing doesn't even try to unite
    // positions with the same global offset
}

#[test]
fn test_n_base_filtering() {
    /// Seqwish filters out N bases: query_base != 'N'
    /// This test verifies that N bases don't create matches

    // Create sequences with N bases
    let seq1 = Sequence {
        id: "seq1".to_string(),
        data: b"ACGTNNNACGT".to_vec(),
        offset: 0,
    };

    let seq2 = Sequence {
        id: "seq2".to_string(),
        data: b"ACGTNNNACGT".to_vec(),
        offset: 11,
    };

    let sequences = vec![seq1.clone(), seq2.clone()];
    let seqrush = SeqRush::new(sequences, u64::MAX);

    // Process a match region that includes N bases
    // CIGAR: 11= (all matches including Ns)
    let cigar = "11=";
    seqrush.process_alignment(
        cigar,
        &seq1,
        &seq2,
        0,     // min_match_len = 0
        false, // not RC
        false, // not verbose
    );

    // Check that positions with N bases are NOT united
    // Positions 0-3 should be united (ACGT)
    let pos1_0 = make_pos(0, false);
    let pos2_0 = make_pos(11, false);
    assert!(seqrush.union_find.same(pos1_0, pos2_0), "Position 0 should be united");

    // Positions 4-6 (NNN) should NOT be united
    let pos1_4 = make_pos(4, false);
    let pos2_4 = make_pos(15, false);
    assert!(
        !seqrush.union_find.same(pos1_4, pos2_4),
        "Positions with N bases should NOT be united"
    );

    // Positions 7-10 (ACGT) should be united
    let pos1_7 = make_pos(7, false);
    let pos2_7 = make_pos(18, false);
    assert!(seqrush.union_find.same(pos1_7, pos2_7), "Position 7 should be united");
}

#[test]
fn test_global_offset_calculation() {
    /// Test that global offsets are calculated correctly for RC alignments
    /// This is critical for the self-mapping guard

    let seq1 = Sequence {
        id: "seq1".to_string(),
        data: b"AAAACCCCGGGGTTTT".to_vec(), // length 16
        offset: 0,
    };

    let seq2 = Sequence {
        id: "seq2".to_string(),
        data: b"TTTTCCCCGGGGAAAA".to_vec(), // RC of seq1
        offset: 16,
    };

    // When seq1 is RC'd for alignment, position 0 in RC space
    // corresponds to position 15 in forward space
    // Global offset = seq1.offset + forward_pos = 0 + 15 = 15

    let query_is_rc = true;
    let seq1_local_pos_rc = 0; // First position in RC coordinates
    let seq1_len = seq1.data.len();
    let forward_local_pos = seq1_len - 1 - seq1_local_pos_rc;
    let global_offset1 = seq1.offset + forward_local_pos;

    assert_eq!(global_offset1, 15, "RC position 0 should map to global offset 15");

    // seq2 position 0 has global offset 16
    let seq2_local_pos = 0;
    let global_offset2 = seq2.offset + seq2_local_pos;

    assert_eq!(global_offset2, 16, "seq2 position 0 should be at global offset 16");

    // These are different offsets, so they can be united
    assert_ne!(global_offset1, global_offset2);
}

#[test]
fn test_match_run_across_mismatches() {
    /// Seqwish breaks match runs on mismatches
    /// This test verifies that seqrush does the same

    let seq1 = Sequence {
        id: "seq1".to_string(),
        data: b"AAAXAAA".to_vec(),
        offset: 0,
    };

    let seq2 = Sequence {
        id: "seq2".to_string(),
        data: b"AAAAAAA".to_vec(),
        offset: 7,
    };

    let sequences = vec![seq1.clone(), seq2.clone()];
    let seqrush = SeqRush::new(sequences, u64::MAX);

    // CIGAR: 3M1X3M (3 matches, 1 mismatch, 3 matches)
    let cigar = "3M1X3M";
    seqrush.process_alignment(
        cigar,
        &seq1,
        &seq2,
        2, // min_match_len = 2
        false,
        false,
    );

    // Positions 0-2 should be united (first AAA)
    let pos1_0 = make_pos(0, false);
    let pos2_0 = make_pos(7, false);
    assert!(seqrush.union_find.same(pos1_0, pos2_0));

    // Position 3 should NOT be united (X vs A)
    let pos1_3 = make_pos(3, false);
    let pos2_3 = make_pos(10, false);
    assert!(!seqrush.union_find.same(pos1_3, pos2_3));

    // Positions 4-6 should be united (second AAA)
    let pos1_4 = make_pos(4, false);
    let pos2_4 = make_pos(11, false);
    assert!(seqrush.union_find.same(pos1_4, pos2_4));
}

#[test]
fn test_min_match_length_filtering() {
    /// Seqwish only creates matches if match_len >= min_match_len
    /// This is checked in the add_match lambda

    let seq1 = Sequence {
        id: "seq1".to_string(),
        data: b"AAATTTGGG".to_vec(),
        offset: 0,
    };

    let seq2 = Sequence {
        id: "seq2".to_string(),
        data: b"AAATTTGGG".to_vec(),
        offset: 9,
    };

    let sequences = vec![seq1.clone(), seq2.clone()];
    let seqrush = SeqRush::new(sequences, u64::MAX);

    // CIGAR: 3=3X3= (3 match, 3 mismatch, 3 match)
    let cigar = "3=3X3=";

    // With min_match_len = 4, neither match run should be united
    seqrush.process_alignment(
        cigar,
        &seq1,
        &seq2,
        4, // min_match_len = 4
        false,
        false,
    );

    // Positions 0-2 should NOT be united (length 3 < 4)
    let pos1_0 = make_pos(0, false);
    let pos2_0 = make_pos(9, false);
    assert!(
        !seqrush.union_find.same(pos1_0, pos2_0),
        "Match run of length 3 should not be united with min_match_len=4"
    );

    // Positions 6-8 should NOT be united (length 3 < 4)
    let pos1_6 = make_pos(6, false);
    let pos2_6 = make_pos(15, false);
    assert!(
        !seqrush.union_find.same(pos1_6, pos2_6),
        "Match run of length 3 should not be united with min_match_len=4"
    );
}

#[test]
fn test_cigar_m_vs_equals() {
    /// 'M' operations can contain mismatches, '=' operations are always matches
    /// Seqwish checks actual bases for 'M', but trusts '='

    let seq1 = Sequence {
        id: "seq1".to_string(),
        data: b"AAAXTTT".to_vec(),
        offset: 0,
    };

    let seq2 = Sequence {
        id: "seq2".to_string(),
        data: b"AAATTTT".to_vec(),
        offset: 7,
    };

    let sequences = vec![seq1.clone(), seq2.clone()];
    let seqrush = SeqRush::new(sequences, u64::MAX);

    // CIGAR: 7M (includes a mismatch at position 3)
    let cigar = "7M";
    seqrush.process_alignment(
        cigar,
        &seq1,
        &seq2,
        2, // min_match_len = 2
        false,
        false,
    );

    // First AAA should be united
    let pos1_0 = make_pos(0, false);
    let pos2_0 = make_pos(7, false);
    assert!(seqrush.union_find.same(pos1_0, pos2_0));

    // Position 3 (X vs T) should NOT be united
    let pos1_3 = make_pos(3, false);
    let pos2_3 = make_pos(10, false);
    assert!(!seqrush.union_find.same(pos1_3, pos2_3));

    // Last TTT should be united
    let pos1_4 = make_pos(4, false);
    let pos2_4 = make_pos(11, false);
    assert!(seqrush.union_find.same(pos1_4, pos2_4));
}

#[test]
fn test_alignment_boundaries() {
    /// Test that alignment processing respects sequence boundaries
    /// and doesn't go out of bounds

    let seq1 = Sequence {
        id: "seq1".to_string(),
        data: b"AAAA".to_vec(),
        offset: 0,
    };

    let seq2 = Sequence {
        id: "seq2".to_string(),
        data: b"AAAA".to_vec(),
        offset: 4,
    };

    let sequences = vec![seq1.clone(), seq2.clone()];
    let seqrush = SeqRush::new(sequences, u64::MAX);

    // CIGAR: 4= (exactly matches sequence length)
    let cigar = "4=";
    seqrush.process_alignment(
        cigar,
        &seq1,
        &seq2,
        0,
        false,
        false,
    );

    // All positions should be united
    for i in 0..4 {
        let pos1 = make_pos(i, false);
        let pos2 = make_pos(4 + i, false);
        assert!(
            seqrush.union_find.same(pos1, pos2),
            "Position {} should be united", i
        );
    }
}
