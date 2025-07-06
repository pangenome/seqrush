use seqrush::inversion_aware_seqrush::*;
use seqrush::seqrush_bidirected_simplified::*;
use seqrush::seqrush::Args;
use seqrush::bidirected_graph::reverse_complement;
use std::fs;
use std::io::Write;
use tempfile::NamedTempFile;
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

fn create_test_fasta(sequences: &[(&str, &str)]) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    for (id, seq) in sequences {
        writeln!(file, ">{}", id).unwrap();
        writeln!(file, "{}", seq).unwrap();
    }
    file.flush().unwrap();
    file
}

// Sequence generation utilities
struct SequenceGenerator {
    rng: StdRng,
}

impl SequenceGenerator {
    fn new(seed: u64) -> Self {
        Self {
            rng: StdRng::seed_from_u64(seed),
        }
    }
    
    fn random_dna(&mut self, length: usize) -> String {
        const BASES: [char; 4] = ['A', 'T', 'C', 'G'];
        (0..length)
            .map(|_| BASES[self.rng.gen_range(0..4)])
            .collect()
    }
    
    fn random_dna_weighted(&mut self, length: usize, gc_content: f64) -> String {
        (0..length)
            .map(|_| {
                if self.rng.gen::<f64>() < gc_content {
                    if self.rng.gen_bool(0.5) { 'G' } else { 'C' }
                } else {
                    if self.rng.gen_bool(0.5) { 'A' } else { 'T' }
                }
            })
            .collect()
    }
    
    // Generate sequence with specific motif density
    fn with_motifs(&mut self, length: usize, motif: &str, density: f64) -> String {
        let mut seq = self.random_dna(length);
        let motif_count = ((length as f64) * density / (motif.len() as f64)) as usize;
        
        for _ in 0..motif_count {
            let pos = self.rng.gen_range(0..length.saturating_sub(motif.len()));
            seq.replace_range(pos..pos + motif.len(), motif);
        }
        
        seq
    }
    
    // Apply random inversions
    fn apply_random_inversions(&mut self, seq: &str, count: usize, min_size: usize, max_size: usize) -> String {
        let mut result = seq.to_string();
        let seq_len = seq.len();
        
        for _ in 0..count {
            if seq_len <= max_size {
                break;
            }
            
            let size = self.rng.gen_range(min_size..=max_size.min(seq_len / 2));
            let start = self.rng.gen_range(0..seq_len.saturating_sub(size));
            let end = start + size;
            
            let prefix = &result[..start];
            let inverted = String::from_utf8(
                reverse_complement(result[start..end].as_bytes())
            ).unwrap();
            let suffix = &result[end..];
            
            result = format!("{}{}{}", prefix, inverted, suffix);
        }
        
        result
    }
    
    // Generate CNV with variable copy numbers
    fn generate_cnv(&mut self, unit: &str, min_copies: usize, max_copies: usize, samples: usize) -> Vec<String> {
        (0..samples)
            .map(|_| {
                let copies = self.rng.gen_range(min_copies..=max_copies);
                unit.repeat(copies)
            })
            .collect()
    }
    
    // Generate complex rearrangement
    fn complex_rearrangement(&mut self, segments: &[&str]) -> String {
        let mut indices: Vec<usize> = (0..segments.len()).collect();
        
        // Shuffle segments
        for i in (1..indices.len()).rev() {
            let j = self.rng.gen_range(0..=i);
            indices.swap(i, j);
        }
        
        // Apply inversions to some segments
        segments.iter().enumerate()
            .map(|(_i, &seg)| {
                if self.rng.gen_bool(0.3) { // 30% chance of inversion
                    String::from_utf8(reverse_complement(seg.as_bytes())).unwrap()
                } else {
                    seg.to_string()
                }
            })
            .collect::<Vec<_>>()
            .join("")
    }
}

#[test]
fn test_programmatic_inversions_scale() {
    println!("\n=== Testing Programmatic Inversions at Scale ===");
    
    let mut gen = SequenceGenerator::new(12345);
    let base_seq = gen.random_dna(1000);
    
    // Generate sequences with increasing number of inversions
    let inversion_counts = vec![1, 2, 5, 10];
    let mut sequences = vec![("original", base_seq.clone())];
    
    for &count in &inversion_counts {
        let seq_with_inv = gen.apply_random_inversions(&base_seq, count, 20, 100);
        sequences.push((
            Box::leak(format!("inv_{}", count).into_boxed_str()),
            seq_with_inv
        ));
    }
    
    let sequences_ref: Vec<(&str, &str)> = sequences.iter()
        .map(|(id, seq)| (*id as &str, seq.as_str()))
        .collect();
    
    let fasta = create_test_fasta(&sequences_ref);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 4,
        min_match_length: 10,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: Some(0.1),
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_inversion_aware_seqrush(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    
    // Verify all sequences are represented
    let path_count = gfa_content.lines()
        .filter(|l| l.starts_with('P'))
        .count();
    assert_eq!(path_count, inversion_counts.len() + 1);
}

#[test]
fn test_gc_content_variation() {
    println!("\n=== Testing GC Content Variation ===");
    
    let mut gen = SequenceGenerator::new(54321);
    let length = 500;
    
    // Generate sequences with different GC content
    let gc_contents = vec![0.2, 0.35, 0.5, 0.65, 0.8];
    let mut sequences = Vec::new();
    
    for &gc in &gc_contents {
        let seq = gen.random_dna_weighted(length, gc);
        
        // Apply inversion to GC-rich region
        let inv_seq = if gc > 0.5 {
            gen.apply_random_inversions(&seq, 1, 50, 100)
        } else {
            seq.clone()
        };
        
        sequences.push((
            Box::leak(format!("gc_{}", (gc * 100.0) as u32).into_boxed_str()),
            inv_seq
        ));
    }
    
    let sequences_ref: Vec<(&str, &str)> = sequences.iter()
        .map(|(id, seq)| (*id as &str, seq.as_str()))
        .collect();
    
    let fasta = create_test_fasta(&sequences_ref);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 8,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: Some(0.2),
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_seqrush_bidirected_simple(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    assert!(gfa_content.contains("VN:Z:1.0"));
}

#[test]
fn test_motif_based_variations() {
    println!("\n=== Testing Motif-Based Variations ===");
    
    let mut gen = SequenceGenerator::new(99999);
    
    // Common genomic motifs
    let motifs = vec![
        ("TATA", 0.05),     // TATA box
        ("GGGCGG", 0.02),   // GC box
        ("AATAAA", 0.03),   // PolyA signal
        ("CTCF", 0.04),     // CTCF binding
    ];
    
    let base_length = 800;
    let mut sequences = Vec::new();
    
    for (motif, density) in motifs {
        let seq = gen.with_motifs(base_length, motif, density);
        
        // Create variant with inverted motif regions
        let mut variant = seq.clone();
        let motif_positions: Vec<usize> = seq.match_indices(motif)
            .map(|(pos, _)| pos)
            .collect();
        
        // Invert around some motif occurrences
        for (i, &pos) in motif_positions.iter().enumerate() {
            if i % 2 == 0 && pos > 10 && pos + motif.len() + 10 < variant.len() {
                let start = pos.saturating_sub(5);
                let end = (pos + motif.len() + 5).min(variant.len());
                
                let prefix = variant[..start].to_string();
                let inverted = String::from_utf8(
                    reverse_complement(variant[start..end].as_bytes())
                ).unwrap();
                let suffix = variant[end..].to_string();
                
                variant = format!("{}{}{}", prefix, inverted, suffix);
            }
        }
        
        sequences.push((
            Box::leak(format!("{}_normal", motif).into_boxed_str()),
            seq
        ));
        sequences.push((
            Box::leak(format!("{}_variant", motif).into_boxed_str()),
            variant
        ));
    }
    
    let sequences_ref: Vec<(&str, &str)> = sequences.iter()
        .map(|(id, seq)| (*id as &str, seq.as_str()))
        .collect();
    
    let fasta = create_test_fasta(&sequences_ref);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 6,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_inversion_aware_seqrush(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let path_count = gfa_content.lines()
        .filter(|l| l.starts_with('P'))
        .count();
    assert_eq!(path_count, sequences.len());
}

#[test]
fn test_microsatellite_instability() {
    println!("\n=== Testing Microsatellite Instability ===");
    
    let mut gen = SequenceGenerator::new(11111);
    
    // Common microsatellite patterns
    let microsats = vec![
        ("A", 10, 20),      // PolyA
        ("CA", 8, 15),      // Dinucleotide
        ("CAG", 5, 25),     // Trinucleotide (Huntington's)
        ("GATA", 4, 12),    // Tetranucleotide
    ];
    
    let mut sequences = Vec::new();
    let flanking = gen.random_dna(50);
    
    for (unit, min_rep, max_rep) in microsats {
        let variants = gen.generate_cnv(unit, min_rep, max_rep, 3);
        
        for (i, var) in variants.iter().enumerate() {
            // Some variants have inversions near the microsatellite
            let full_seq = if i == 1 {
                let ms_with_flank = format!("{}{}{}", &flanking[..25], var, &flanking[25..]);
                gen.apply_random_inversions(&ms_with_flank, 1, 20, 40)
            } else {
                format!("{}{}{}", &flanking[..25], var, &flanking[25..])
            };
            
            sequences.push((
                Box::leak(format!("{}_var{}", unit, i).into_boxed_str()),
                full_seq
            ));
        }
    }
    
    let sequences_ref: Vec<(&str, &str)> = sequences.iter()
        .map(|(id, seq)| (*id as &str, seq.as_str()))
        .collect();
    
    let fasta = create_test_fasta(&sequences_ref);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 4,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: Some(0.1),
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_inversion_aware_seqrush(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    assert!(gfa_content.contains("VN:Z:1.0"));
}

#[test]
fn test_complex_structural_variation_patterns() {
    println!("\n=== Testing Complex Structural Variation Patterns ===");
    
    let mut gen = SequenceGenerator::new(33333);
    
    // Create segments that will be rearranged
    let segments = vec![
        gen.random_dna(100),
        gen.random_dna(150),
        gen.random_dna(80),
        gen.random_dna(120),
        gen.random_dna(90),
    ];
    
    let segment_refs: Vec<&str> = segments.iter()
        .map(|s| s.as_str())
        .collect();
    
    // Generate different rearrangement patterns
    let mut sequences = vec![
        ("reference", segments.join("")),
    ];
    
    // Simple inversion of middle segment
    let mut simple_inv = segments.clone();
    simple_inv[2] = String::from_utf8(reverse_complement(segments[2].as_bytes())).unwrap();
    sequences.push(("simple_inv", simple_inv.join("")));
    
    // Complex rearrangement
    let complex1 = gen.complex_rearrangement(&segment_refs);
    sequences.push(("complex1", complex1));
    
    // Another complex rearrangement
    let complex2 = gen.complex_rearrangement(&segment_refs);
    sequences.push(("complex2", complex2));
    
    // Deletion + inversion
    let mut del_inv = segments[..2].to_vec();
    del_inv.push(String::from_utf8(reverse_complement(segments[3].as_bytes())).unwrap());
    del_inv.push(segments[4].clone());
    sequences.push(("del_inv", del_inv.join("")));
    
    let sequences_ref: Vec<(&str, &str)> = sequences.iter()
        .map(|(id, seq)| (*id, seq.as_str()))
        .collect();
    
    let fasta = create_test_fasta(&sequences_ref);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 4,
        min_match_length: 15,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: Some(0.1),
        verbose: true,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_inversion_aware_seqrush(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    
    // Complex rearrangements should create a graph with branching
    let edges: Vec<&str> = gfa_content.lines()
        .filter(|l| l.starts_with('L'))
        .collect();
    
    let nodes: Vec<&str> = gfa_content.lines()
        .filter(|l| l.starts_with('S'))
        .collect();
    
    // Graph should represent all sequences
    println!("Graph complexity: {} nodes, {} edges", nodes.len(), edges.len());
    assert!(nodes.len() > 0);
    assert!(edges.len() > 0);
    
    let paths: Vec<&str> = gfa_content.lines()
        .filter(|l| l.starts_with('P'))
        .collect();
    assert_eq!(paths.len(), sequences.len());
}

#[test]
fn test_nested_tandem_inversions() {
    println!("\n=== Testing Nested Tandem Inversions ===");
    
    let _gen = SequenceGenerator::new(77777);
    
    // Create a pattern with nested structure
    let core = "ATCG";
    let level1 = format!("{}{}{}", core, "TTTT", core); // ATCGTTTTATCG
    let level2 = format!("{}{}{}", level1, "AAAA", level1); // (ATCGTTTTATCG)AAAA(ATCGTTTTATCG)
    
    // Create variations
    let mut sequences = vec![
        ("original", level2.clone()),
    ];
    
    // Invert inner units
    let inner_inv = level2.replace(
        &level1,
        &String::from_utf8(reverse_complement(level1.as_bytes())).unwrap()
    );
    sequences.push(("inner_inv", inner_inv));
    
    // Invert entire structure
    let full_inv = String::from_utf8(reverse_complement(level2.as_bytes())).unwrap();
    sequences.push(("full_inv", full_inv));
    
    // Mixed: some units normal, some inverted
    let mixed = format!(
        "{}{}{}",
        level1,
        "AAAA",
        String::from_utf8(reverse_complement(level1.as_bytes())).unwrap()
    );
    sequences.push(("mixed", mixed));
    
    let sequences_ref: Vec<(&str, &str)> = sequences.iter()
        .map(|(id, seq)| (*id, seq.as_str()))
        .collect();
    
    let fasta = create_test_fasta(&sequences_ref);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 4,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: None,
        verbose: true,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_inversion_aware_seqrush(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let paths: Vec<&str> = gfa_content.lines()
        .filter(|l| l.starts_with('P'))
        .collect();
    assert_eq!(paths.len(), sequences.len());
}