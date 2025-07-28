// Test program to understand WFA2 CIGAR convention
fn main() {
    // Create test sequences
    let query = "ATCGATCGATCG";    // 12 bases
    let target = "ATCGATTTGATCG";  // 13 bases - has 2 extra T's
    
    println!("Query:  {}", query);
    println!("Target: {}", target);
    println!();
    
    // Manual alignment:
    // Query:  ATCGAT--CGATCG
    // Target: ATCGATTTGATCG
    //         ^^^^^^  ^^^^^^
    //         6 match, 2 insertion in target, 6 match
    
    println!("Expected alignment:");
    println!("Query:  ATCGAT--CGATCG");
    println!("Target: ATCGATTTGATCG");
    println!();
    
    println!("Standard CIGAR: 6M2I6M");
    println!("  - 6M: first 6 bases match");
    println!("  - 2I: 2 bases inserted in TARGET (TT)");
    println!("  - 6M: last 6 bases match");
    println!();
    
    println!("If WFA2 uses swapped convention:");
    println!("  - I means insertion TO TRANSFORM query to target");
    println!("  - So if query needs insertions to become target,");
    println!("  - that means target has extra bases");
    println!("  - This matches standard CIGAR!");
}