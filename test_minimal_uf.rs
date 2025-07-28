#[test]
fn test_uf_rush_directly() {
    use uf_rush::UFRush;
    
    let uf = UFRush::new(1000);
    
    // Test with the exact values from the failing test
    let p1 = 279usize;  // make_pos(139, true)
    let p2 = 430usize;  // make_pos(215, false)
    
    println!("Before unite: p1={}, p2={}", p1, p2);
    uf.unite(p1, p2);
    
    let r1 = uf.find(p1);
    let r2 = uf.find(p2);
    println!("After unite: find(p1)={}, find(p2)={}", r1, r2);
    
    assert_eq!(r1, r2, "Representatives should be equal after unite");
}

fn main() {
    test_uf_rush_directly();
    println!("Direct UFRush test passed!");
}