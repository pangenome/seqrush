use lib_wfa2::affine_wavefront::AffineWavefronts;

fn main() {
    // Test sequences
    let query = b"ATCGATCGATCG";    // 12 bases
    let target = b"ATCGATTTGATCG";  // 13 bases - has TTT instead of just T
    
    let aligner = AffineWavefronts::default();
    aligner.align(query, target);
    
    println!("Query:  {}", String::from_utf8_lossy(query));
    println!("Target: {}", String::from_utf8_lossy(target));
    println!("Score: {}", aligner.score());
    println!("CIGAR: {}", String::from_utf8_lossy(aligner.cigar()));
    
    // Now check what this CIGAR means
    let cigar = String::from_utf8_lossy(aligner.cigar());
    let mut query_pos = 0;
    let mut target_pos = 0;
    let mut count = 0;
    
    println!("\nCIGAR interpretation:");
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            count = count * 10 + (ch as usize - '0' as usize);
        } else {
            if count == 0 { count = 1; }
            
            match ch {
                'M' => {
                    println!("  {}M: {} query bases match {} target bases", count, count, count);
                    query_pos += count;
                    target_pos += count;
                }
                'I' => {
                    println!("  {}I: {} bases inserted (present in target, not in query)", count, count);
                    target_pos += count;
                }
                'D' => {
                    println!("  {}D: {} bases deleted (present in query, not in target)", count, count);
                    query_pos += count;
                }
                'X' => {
                    println!("  {}X: {} mismatches", count, count);
                    query_pos += count;
                    target_pos += count;
                }
                _ => {}
            }
            count = 0;
        }
    }
    
    println!("\nFinal positions: query={}/{}, target={}/{}", query_pos, query.len(), target_pos, target.len());
}