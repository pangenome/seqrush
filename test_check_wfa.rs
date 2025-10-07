fn main() {
    // Test if WFA2's CIGAR is really swapped
    
    // Test case 1: Query has insertion
    let query1 = "ATCGATTTCGATCG";  // 14 bases - has extra TT
    let target1 = "ATCGATCGATCG";   // 12 bases
    
    println!("Test 1: Query longer than target");
    println!("Query:  {} (len={})", query1, query1.len());
    println!("Target: {} (len={})", target1, target1.len());
    println!();
    println!("Expected alignment:");
    println!("Query:  ATCGATTTCGATCG");
    println!("Target: ATCGAT--CGATCG");
    println!();
    println!("Standard CIGAR would be: 6M2D6M");
    println!("  - 2D means 2 deletions from REFERENCE");
    println!("  - Since target is reference, this means query has extra bases");
    println!();
    
    // Test case 2: Target has insertion 
    let query2 = "ATCGATCGATCG";    // 12 bases
    let target2 = "ATCGATTTGATCG";  // 13 bases - has extra T
    
    println!("Test 2: Target longer than query");
    println!("Query:  {} (len={})", query2, query2.len());
    println!("Target: {} (len={})", target2, target2.len());
    println!();
    println!("Expected alignment:");
    println!("Query:  ATCGAT-CGATCG");
    println!("Target: ATCGATTTGATCG");
    println!();
    println!("Standard CIGAR would be: 6M1I6M");
    println!("  - 1I means 1 insertion in REFERENCE");
    println!("  - Since target is reference, this means target has extra base");
}