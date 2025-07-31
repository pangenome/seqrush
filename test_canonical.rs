fn main() {
    // Test canonical edge logic
    let h1_fwd = 559u64 << 1;      // Forward handle for node 559
    let h1_rev = (559u64 << 1) | 1; // Reverse handle for node 559
    let h2_fwd = 560u64 << 1;      // Forward handle for node 560
    let h2_rev = (560u64 << 1) | 1; // Reverse handle for node 560
    
    println!("h1_fwd = {}, h1_rev = {}", h1_fwd, h1_rev);
    println!("h2_fwd = {}, h2_rev = {}", h2_fwd, h2_rev);
    
    // Check ordering
    println!("h1_fwd < h2_fwd: {}", h1_fwd < h2_fwd);
    println!("h1_rev < h2_rev: {}", h1_rev < h2_rev);
    
    // Simulate canonical form logic
    // Edge 559+ -> 560+
    if h1_fwd <= h2_fwd {
        println!("Edge 559+ -> 560+ canonical: 559+ -> 560+");
    }
    
    // Edge 559- -> 560-
    if h1_rev <= h2_rev {
        println!("Edge 559- -> 560- canonical: 559- -> 560-");
    }
    
    // Edge 560- -> 559- (reverse direction)
    if h2_rev > h1_rev {
        // Need to flip both
        println!("Edge 560- -> 559- canonical: 559+ -> 560+ (flipped both)");
    }
}