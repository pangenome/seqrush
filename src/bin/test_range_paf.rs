use seqrush::range_builder::{RangeBasedGraphBuilder, AlignmentRange};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::collections::HashMap;

fn parse_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
    let file = File::open(path).expect("Failed to open FASTA");
    let reader = BufReader::new(file);
    let mut sequences = Vec::new();
    let mut current_name = String::new();
    let mut current_seq = Vec::new();
    
    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('>') {
            if !current_name.is_empty() {
                sequences.push((current_name.clone(), current_seq.clone()));
            }
            current_name = line[1..].to_string();
            current_seq.clear();
        } else {
            current_seq.extend(line.as_bytes());
        }
    }
    
    if !current_name.is_empty() {
        sequences.push((current_name, current_seq));
    }
    
    sequences
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let fasta_path = "HLA-zoo/seqs/B-3106.fa";
    let paf_path = "b.paf";
    
    // Load sequences
    let sequences = parse_fasta(fasta_path);
    let mut builder = RangeBasedGraphBuilder::new();
    
    // Calculate offsets for sequences in the graph
    let mut seq_offsets = HashMap::new();
    let mut total_offset = 0;
    
    for (name, data) in &sequences {
        seq_offsets.insert(name.clone(), total_offset);
        builder.add_sequence(name.clone(), data.clone());
        total_offset += data.len();
    }
    
    // Parse PAF and add alignments
    let file = File::open(paf_path)?;
    let reader = BufReader::new(file);
    let mut alignment_count = 0;
    
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }
        
        let query_name = fields[0];
        let query_start = fields[2].parse::<usize>().unwrap();
        let query_end = fields[3].parse::<usize>().unwrap();
        let query_strand = fields[4];
        let target_name = fields[5];
        let target_start = fields[7].parse::<usize>().unwrap();
        let target_end = fields[8].parse::<usize>().unwrap();
        
        // Skip self-alignments for now (they're added automatically)
        if query_name == target_name && query_start == target_start && query_end == target_end {
            continue;
        }
        
        // Get offsets
        let query_offset = seq_offsets.get(query_name).cloned().unwrap_or(0);
        let target_offset = seq_offsets.get(target_name).cloned().unwrap_or(0);
        
        // Add alignment range
        builder.add_alignment_range(AlignmentRange {
            seq1_start: query_offset + query_start,
            seq1_end: query_offset + query_end,
            seq2_start: target_offset + target_start,
            seq2_end: target_offset + target_end,
            seq2_is_rc: query_strand == "-",
        });
        
        alignment_count += 1;
    }
    
    println!("Loaded {} sequences and {} alignments", sequences.len(), alignment_count);
    
    // Build the graph
    let graph = builder.build_graph(true)?;
    
    // Write to GFA
    let mut file = File::create("test_range_paf.gfa")?;
    graph.write_gfa(&mut file)?;
    
    println!("Created test_range_paf.gfa with {} nodes", graph.nodes.len());
    
    // Compare with seqwish
    println!("\nComparison:");
    println!("SeqRush (range-based): {} nodes", graph.nodes.len());
    println!("seqwish: 471 nodes (expected)");
    
    Ok(())
}