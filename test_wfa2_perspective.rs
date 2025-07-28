use lib_wfa2::affine_wavefront::AffineWavefronts;

fn main() {
    // Test to understand WFA2's CIGAR perspective
    let aligner = AffineWavefronts::default();
    
    // Simple test: query is shorter than target
    let query = b"ATCG";
    let target = b"ATCGATCG";
    
    aligner.align(query, target);
    
    println!("Query ({}bp): {}", query.len(), String::from_utf8_lossy(query));
    println!("Target ({}bp): {}", target.len(), String::from_utf8_lossy(target));
    println!("WFA2 CIGAR: {}", String::from_utf8_lossy(aligner.cigar()));
    println!("Score: {}", aligner.score());
    println!();
    
    // Another test: query is longer than target
    let query2 = b"ATCGATCG";
    let target2 = b"ATCG";
    
    aligner.align(query2, target2);
    
    println!("Query ({}bp): {}", query2.len(), String::from_utf8_lossy(query2));
    println!("Target ({}bp): {}", target2.len(), String::from_utf8_lossy(target2));
    println!("WFA2 CIGAR: {}", String::from_utf8_lossy(aligner.cigar()));
    println!("Score: {}", aligner.score());
}