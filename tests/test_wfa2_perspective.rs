use lib_wfa2::affine_wavefront::{AffineWavefronts, MemoryMode};

#[test]
fn test_wfa2_cigar_perspective() {
    // Test to understand WFA2's CIGAR perspective
    let aligner = AffineWavefronts::with_penalties_and_memory_mode(0, 5, 8, 2, MemoryMode::Ultralow);
    
    // Test 1: query is shorter than target
    let query = b"ATCG";
    let target = b"ATCGATCG";
    
    aligner.align(query, target);
    
    println!("\nTest 1:");
    println!("Query ({}bp): {}", query.len(), String::from_utf8_lossy(query));
    println!("Target ({}bp): {}", target.len(), String::from_utf8_lossy(target));
    println!("WFA2 CIGAR: {}", String::from_utf8_lossy(aligner.cigar()));
    println!("Score: {}", aligner.score());
    
    // Expected from standard perspective:
    // Query:  ATCG----
    // Target: ATCGATCG
    // Standard CIGAR would be: 4=4I (4 matches, 4 insertions in target)
    // But if WFA2 swaps I/D, it might say: 4=4D
    
    let cigar1 = String::from_utf8_lossy(aligner.cigar());
    
    // Test 2: query is longer than target  
    let query2 = b"ATCGATCG";
    let target2 = b"ATCG";
    
    aligner.align(query2, target2);
    
    println!("\nTest 2:");
    println!("Query ({}bp): {}", query2.len(), String::from_utf8_lossy(query2));
    println!("Target ({}bp): {}", target2.len(), String::from_utf8_lossy(target2));
    println!("WFA2 CIGAR: {}", String::from_utf8_lossy(aligner.cigar()));
    println!("Score: {}", aligner.score());
    
    // Expected from standard perspective:
    // Query:  ATCGATCG
    // Target: ATCG----
    // Standard CIGAR would be: 4=4D (4 matches, 4 deletions in target)
    // But if WFA2 swaps I/D, it might say: 4=4I
    
    let cigar2 = String::from_utf8_lossy(aligner.cigar());
    
    println!("\nAnalysis:");
    println!("When query < target: WFA2 says {}", cigar1);
    println!("When query > target: WFA2 says {}", cigar2);
    println!("\nThis confirms whether WFA2's CIGAR is from standard perspective or swapped.");
}