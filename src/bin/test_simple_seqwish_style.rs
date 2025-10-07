use seqrush::seqwish_style::SeqwishStyleBuilder;
use std::fs::File;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create simple test
    let mut builder = SeqwishStyleBuilder::new();
    builder.add_sequence("seq1".to_string(), b"ACGT".to_vec());
    builder.add_sequence("seq2".to_string(), b"ACGT".to_vec());
    
    println!("Building graph for simple test...");
    let graph = builder.build_from_paf("simple_test.paf", 0, true)?;
    
    // Write to GFA
    let mut file = File::create("simple_seqrush_style.gfa")?;
    graph.write_gfa(&mut file)?;
    
    println!("\nCreated simple_seqrush_style.gfa with {} nodes", graph.nodes.len());
    println!("Expected: 1 node (like seqwish)");
    
    Ok(())
}