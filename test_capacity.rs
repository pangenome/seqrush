fn main() {
    let max_offset = 1000;
    let old_capacity = (max_offset * 2) + 1000;
    let new_capacity = (max_offset << 1) + 2;
    
    println!("max_offset: {}", max_offset);
    println!("old capacity formula: ({} * 2) + 1000 = {}", max_offset, old_capacity);
    println!("new capacity formula: ({} << 1) + 2 = {}", max_offset, new_capacity);
    
    // Test positions from the failing test
    let p1 = (139 << 1) | 1;  // make_pos(139, true)
    let p2 = (215 << 1) | 0;  // make_pos(215, false)
    
    println!("\nTest positions:");
    println!("p1 = make_pos(139, true) = (139 << 1) | 1 = {}", p1);
    println!("p2 = make_pos(215, false) = (215 << 1) | 0 = {}", p2);
    
    println!("\nCapacity check:");
    println!("p1 ({}) < new_capacity ({})? {}", p1, new_capacity, p1 < new_capacity);
    println!("p2 ({}) < new_capacity ({})? {}", p2, new_capacity, p2 < new_capacity);
}