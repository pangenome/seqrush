use seqrush::seqwish_style::SeqwishStyleBuilder;
use std::fs::File;
use std::io::{BufRead, BufReader};

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
            // Only take the first part of the header (up to first space)
            current_name = line[1..].split_whitespace().next().unwrap_or("").to_string();
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
    let mut builder = SeqwishStyleBuilder::new();
    
    for (name, data) in sequences {
        println!("DEBUG: Adding sequence '{}'", name);
        builder.add_sequence(name, data);
    }
    
    println!("Building graph using seqwish-style algorithm...");
    
    // Build graph
    let graph = builder.build_from_paf(paf_path, 0, true)?;
    
    // Write to GFA
    let mut file = File::create("b_seqwish_style.gfa")?;
    graph.write_gfa(&mut file)?;
    
    println!("\nCreated b_seqwish_style.gfa with {} nodes", graph.nodes.len());
    println!("Expected: 471 nodes (seqwish output)");
    
    Ok(())
}