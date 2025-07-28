use seqrush::bidirected_union_find::BidirectedUnionFind;
use seqrush::pos::make_pos;

fn main() {
    let uf = BidirectedUnionFind::new(1000);
    
    // Reverse-forward alignment (seq1/query was reverse complemented)
    uf.unite_matching_region(
        100, // seq1_offset
        200, // seq2_offset
        10,  // seq1_local_start (in RC coordinates)
        15,  // seq2_local_start
        3,   // match_length
        true, // seq1_is_rc
        50,   // seq1_len
    );
    
    // Test the first position
    let pos1_start = make_pos(139, true);   // 100 + (50-1-10) = 139, reverse
    let pos2_start = make_pos(215, false);  // 200 + 15 = 215, forward
    
    println!("pos1_start: {} (offset: {}, rev: {})", pos1_start, pos1_start >> 1, pos1_start & 1);
    println!("pos2_start: {} (offset: {}, rev: {})", pos2_start, pos2_start >> 1, pos2_start & 1);
    
    let rep1 = uf.find(pos1_start);
    let rep2 = uf.find(pos2_start);
    
    println!("rep1: {} (offset: {}, rev: {})", rep1, rep1 >> 1, rep1 & 1);
    println!("rep2: {} (offset: {}, rev: {})", rep2, rep2 >> 1, rep2 & 1);
    
    println!("Same? {}", uf.same(pos1_start, pos2_start));
    
    // Test the last position
    let pos1_end = make_pos(137, true);    // 100 + (50-1-12) = 137, reverse
    let pos2_end = make_pos(217, false);   // 200 + 15 + 2 = 217, forward
    
    println!("\npos1_end: {} (offset: {}, rev: {})", pos1_end, pos1_end >> 1, pos1_end & 1);
    println!("pos2_end: {} (offset: {}, rev: {})", pos2_end, pos2_end >> 1, pos2_end & 1);
    
    let rep1_end = uf.find(pos1_end);
    let rep2_end = uf.find(pos2_end);
    
    println!("rep1_end: {} (offset: {}, rev: {})", rep1_end, rep1_end >> 1, rep1_end & 1);
    println!("rep2_end: {} (offset: {}, rev: {})", rep2_end, rep2_end >> 1, rep2_end & 1);
    
    println!("Same? {}", uf.same(pos1_end, pos2_end));
}