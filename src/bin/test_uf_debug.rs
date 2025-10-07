use seqrush::bidirected_union_find::BidirectedUnionFind;
use seqrush::pos::make_pos;

fn main() {
    println!("Testing BidirectedUnionFind RC handling...\n");
    
    // Test direct unite
    println!("=== Direct Unite Test ===");
    let uf1 = BidirectedUnionFind::new(1000);
    let pos1 = make_pos(139, true);
    let pos2 = make_pos(215, false);
    println!("pos1 = make_pos(139, true) = {}", pos1);
    println!("pos2 = make_pos(215, false) = {}", pos2);
    
    uf1.unite(pos1, pos2);
    
    let rep1 = uf1.find(pos1);
    let rep2 = uf1.find(pos2);
    println!("After unite:");
    println!("  find(pos1) = {}", rep1);
    println!("  find(pos2) = {}", rep2);
    println!("  same(pos1, pos2) = {}", uf1.same(pos1, pos2));
    
    // Test unite_matching_region
    println!("\n=== Unite Matching Region Test ===");
    let uf2 = BidirectedUnionFind::new(1000);
    
    println!("Calling unite_matching_region with:");
    println!("  seq1_offset: 100");
    println!("  seq2_offset: 200");
    println!("  seq1_local_start: 10 (in RC coordinates)");
    println!("  seq2_local_start: 15");
    println!("  match_length: 3");
    println!("  seq1_is_rc: true");
    println!("  seq1_len: 50");
    
    uf2.unite_matching_region(
        100, // seq1_offset
        200, // seq2_offset
        10,  // seq1_local_start (in RC coordinates)
        15,  // seq2_local_start
        3,   // match_length
        true, // seq1_is_rc
        50,   // seq1_len
    );
    
    println!("\nExpected unions:");
    println!("  Position 0: RC pos 10 -> forward pos 39 -> global 139");
    println!("    Should unite: make_pos(139, true) with make_pos(215, false)");
    
    // Check first position
    let pos1_start = make_pos(139, true);
    let pos2_start = make_pos(215, false);
    println!("\nChecking first position:");
    println!("  pos1_start = {}", pos1_start);
    println!("  pos2_start = {}", pos2_start);
    println!("  same? {}", uf2.same(pos1_start, pos2_start));
    
    // Check last position
    let pos1_end = make_pos(137, true);
    let pos2_end = make_pos(217, false);
    println!("\nChecking last position:");
    println!("  pos1_end = {}", pos1_end);
    println!("  pos2_end = {}", pos2_end);
    println!("  same? {}", uf2.same(pos1_end, pos2_end));
}