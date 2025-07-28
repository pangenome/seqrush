use lib_wfa2::affine_wavefront::{AffineWavefronts, MemoryMode, AlignmentStatus};

#[test]
fn test_cigar_length_consistency() {
    // Test that CIGAR strings consume the correct amount of sequence
    let test_cases = vec![
        // (seq1, seq2, description)
        (b"ATCGATCGATCG".to_vec(), b"ATCGATCGATCG".to_vec(), "identical"),
        (b"ATCGATCGATCG".to_vec(), b"ATCGATCGAT".to_vec(), "deletion at end"),
        (b"ATCGATCGAT".to_vec(), b"ATCGATCGATCG".to_vec(), "insertion at end"),
        (b"ATCGATCGATCG".to_vec(), b"ATTGATCGATCG".to_vec(), "single mismatch"),
    ];
    
    for (seq1, seq2, desc) in test_cases {
        println!("Testing: {}", desc);
        
        let wf = AffineWavefronts::with_penalties_and_memory_mode(0, 5, 8, 2, MemoryMode::Ultralow);
        
        let status = wf.align(&seq1, &seq2);
        assert!(matches!(status, AlignmentStatus::Completed));
        
        let cigar_bytes = wf.cigar();
        let cigar = String::from_utf8_lossy(cigar_bytes);
        println!("  CIGAR: {}", cigar);
        
        // Parse CIGAR to verify lengths
        let (query_consumed, target_consumed) = parse_cigar_consumption(&cigar);
        
        println!("  Query: consumed {} of {} bases", query_consumed, seq1.len());
        println!("  Target: consumed {} of {} bases", target_consumed, seq2.len());
        
        // The CIGAR should consume exactly the sequence lengths
        assert_eq!(query_consumed, seq1.len(), 
            "CIGAR query consumption mismatch for {}", desc);
        assert_eq!(target_consumed, seq2.len(), 
            "CIGAR target consumption mismatch for {}", desc);
    }
}

fn parse_cigar_consumption(cigar: &str) -> (usize, usize) {
    let mut query_consumed = 0;
    let mut target_consumed = 0;
    let mut count = 0;
    
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            count = count * 10 + (ch as usize - '0' as usize);
        } else {
            if count == 0 { count = 1; }
            
            match ch {
                'M' | '=' | 'X' => {
                    query_consumed += count;
                    target_consumed += count;
                }
                'I' => {
                    // WFA2 convention: I means insertions needed to transform query to target
                    // So target has extra bases
                    target_consumed += count;
                }
                'D' => {
                    // WFA2 convention: D means deletions needed from query
                    // So query has extra bases
                    query_consumed += count;
                }
                _ => {}
            }
            count = 0;
        }
    }
    
    (query_consumed, target_consumed)
}

#[test] 
fn test_problematic_alignment() {
    // Test the specific alignment that's causing issues
    // This simulates aligning sequences of length 11068 and 11065
    let seq1 = vec![b'A'; 11068];
    let seq2 = vec![b'A'; 11065];
    
    let wf = AffineWavefronts::with_penalties_and_memory_mode(0, 5, 8, 2, MemoryMode::Ultralow);
    
    let status = wf.align(&seq1, &seq2);
    assert!(matches!(status, AlignmentStatus::Completed));
    
    let cigar_bytes = wf.cigar();
    let cigar = String::from_utf8_lossy(cigar_bytes);
    println!("Problematic CIGAR: {}", cigar);
    
    let (query_consumed, target_consumed) = parse_cigar_consumption(&cigar);
    println!("Query consumed: {} (expected: {})", query_consumed, seq1.len());
    println!("Target consumed: {} (expected: {})", target_consumed, seq2.len());
    
    // These should match the sequence lengths
    assert_eq!(query_consumed, seq1.len());
    assert_eq!(target_consumed, seq2.len());
}