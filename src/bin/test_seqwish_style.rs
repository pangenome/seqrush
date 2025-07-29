use seqrush::range_builder::{RangeBasedGraphBuilder, AlignmentRange};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Test with simple sequences
    let mut builder = RangeBasedGraphBuilder::new();
    
    // Add two identical sequences
    builder.add_sequence("seq1".to_string(), b"ACGT".to_vec());
    builder.add_sequence("seq2".to_string(), b"ACGT".to_vec());
    
    // Add alignment that says seq2 aligns to seq1 completely
    builder.add_alignment_range(AlignmentRange {
        seq1_start: 0,
        seq1_end: 4,
        seq2_start: 4,  // seq2 starts at offset 4 in the graph sequence
        seq2_end: 8,
        seq2_is_rc: false,
    });
    
    // Build the graph
    let graph = builder.build_graph(true)?;
    
    // Write to GFA
    let mut file = File::create("test_seqwish_style.gfa")?;
    graph.write_gfa(&mut file)?;
    
    println!("Created test_seqwish_style.gfa with {} nodes", graph.nodes.len());
    
    Ok(())
}