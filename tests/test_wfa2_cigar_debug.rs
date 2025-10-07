use lib_wfa2::affine_wavefront::{AffineWavefronts, MemoryMode};

#[test]
fn test_wfa2_cigar_parsing() {
    let aligner =
        AffineWavefronts::with_penalties_and_memory_mode(0, 5, 8, 2, MemoryMode::Ultralow);

    // Test case from test_simple.fa
    let seq1 = b"ATCGATCG"; // 8bp
    let seq2 = b"ATCGATCGATCG"; // 12bp

    aligner.align(seq1, seq2);

    println!(
        "\nQuery (seq1): {} ({}bp)",
        String::from_utf8_lossy(seq1),
        seq1.len()
    );
    println!(
        "Target (seq2): {} ({}bp)",
        String::from_utf8_lossy(seq2),
        seq2.len()
    );
    let cigar = String::from_utf8_lossy(aligner.cigar());
    println!("WFA2 CIGAR: {}", cigar);
    println!("Score: {}", aligner.score());

    // Check if this is what we expect
    assert_eq!(cigar, "MMMMMMMMIIII");
}
