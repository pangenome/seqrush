use lib_wfa2::affine_wavefront::AffineWavefronts;

fn main() {
    let aligner = AffineWavefronts::default();
    
    // Test 1: Pattern shorter than text
    println!("Test 1: pattern < text");
    let pattern = b"ATCG";
    let text = b"ATCGATCG";
    
    aligner.align(pattern, text);
    println!("  Pattern ({}bp): {}", pattern.len(), String::from_utf8_lossy(pattern));
    println!("  Text ({}bp): {}", text.len(), String::from_utf8_lossy(text));
    println!("  CIGAR: {}", String::from_utf8_lossy(aligner.cigar()));
    println!("  Score: {}\n", aligner.score());
    
    // Test 2: Pattern longer than text
    println!("Test 2: pattern > text");
    let pattern2 = b"ATCGATCG";
    let text2 = b"ATCG";
    
    aligner.align(pattern2, text2);
    println!("  Pattern ({}bp): {}", pattern2.len(), String::from_utf8_lossy(pattern2));
    println!("  Text ({}bp): {}", text2.len(), String::from_utf8_lossy(text2));
    println!("  CIGAR: {}", String::from_utf8_lossy(aligner.cigar()));
    println!("  Score: {}\n", aligner.score());
    
    // Test 3: What we're doing in seqrush
    println!("Test 3: seqrush style (seq1=query, seq2=target)");
    let seq1 = b"ATCG";      // query
    let seq2 = b"ATCGATCG";  // target
    
    aligner.align(seq1, seq2);
    println!("  seq1/query ({}bp): {}", seq1.len(), String::from_utf8_lossy(seq1));
    println!("  seq2/target ({}bp): {}", seq2.len(), String::from_utf8_lossy(seq2));
    println!("  CIGAR: {}", String::from_utf8_lossy(aligner.cigar()));
    println!("  Score: {}", aligner.score());
    
    println!("\nCIGAR interpretation:");
    println!("  4M4I means: 4 matches, then 4 insertions");
    println!("  In WFA2: insertions are bases in text/target not in pattern/query");
    println!("  In PAF: insertions are bases in query not in target");
    println!("  So we need to swap I/D when converting to PAF!");
}