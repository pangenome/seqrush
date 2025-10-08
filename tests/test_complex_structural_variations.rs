use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use seqrush::bidirected_graph::reverse_complement;
use seqrush::inversion_aware_seqrush::*;
use seqrush::seqrush::Args;
use std::fs;
use std::io::Write;
use tempfile::NamedTempFile;

fn create_test_fasta(sequences: &[(&str, &str)]) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    for (id, seq) in sequences {
        writeln!(file, ">{}", id).unwrap();
        writeln!(file, "{}", seq).unwrap();
    }
    file.flush().unwrap();
    file
}

// Generate random DNA sequence
fn generate_random_dna(rng: &mut StdRng, length: usize) -> String {
    const BASES: [char; 4] = ['A', 'T', 'C', 'G'];
    (0..length).map(|_| BASES[rng.gen_range(0..4)]).collect()
}

// Generate tandem repeat
fn generate_tandem_repeat(unit: &str, copies: usize) -> String {
    unit.repeat(copies)
}

// Generate interspersed repeat
fn generate_interspersed_repeat(unit: &str, spacer: &str, copies: usize) -> String {
    (0..copies)
        .map(|i| {
            if i < copies - 1 {
                format!("{}{}", unit, spacer)
            } else {
                unit.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("")
}

// Apply inversion to a sequence
fn apply_inversion(seq: &str, start: usize, end: usize) -> String {
    let prefix = &seq[..start];
    let inverted = String::from_utf8(reverse_complement(seq[start..end].as_bytes())).unwrap();
    let suffix = &seq[end..];
    format!("{}{}{}", prefix, inverted, suffix)
}

// Apply duplication
fn apply_duplication(seq: &str, start: usize, end: usize, target: usize) -> String {
    let duplicated = &seq[start..end];
    let mut result = seq.to_string();
    result.insert_str(target, duplicated);
    result
}

// Generate complex CNV pattern
fn generate_cnv_pattern(unit: &str, copy_numbers: &[usize]) -> Vec<String> {
    copy_numbers.iter().map(|&n| unit.repeat(n)).collect()
}

#[test]
fn test_inverted_repeats() {
    println!("\n=== Testing Inverted Repeats ===");

    // Create inverted repeat structure: ABC...CBA
    let repeat_unit = "ATCGATCG";
    let spacer = "TTTTAAAA";
    let inv_repeat = String::from_utf8(reverse_complement(repeat_unit.as_bytes())).unwrap();

    let ref_seq = format!("{}{}{}", repeat_unit, spacer, repeat_unit);
    let inv_seq = format!("{}{}{}", repeat_unit, spacer, inv_repeat);

    let sequences = vec![("ref", ref_seq.as_str()), ("inv_repeat", inv_seq.as_str())];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 4,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: true,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: true,
        paf: None,
        seqwish_style: false,
        no_sort: false,
            groom: false,
        iterative_groom: None,
        odgi_style_groom: false,
        sort_groom_sort: false,
        sgd_sort: false,
        aligner: "allwave".to_string(),
    };

    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    assert!(gfa_content.contains("VN:Z:1.0"));

    // Should detect the inverted repeat structure
    let paths: Vec<&str> = gfa_content.lines().filter(|l| l.starts_with('P')).collect();
    assert_eq!(paths.len(), 2);
}

#[test]
fn test_complex_cnv_repeats() {
    println!("\n=== Testing Complex CNV Repeats ===");

    let unit = "ATCG";
    let cnv_patterns = generate_cnv_pattern(unit, &[1, 3, 5, 2, 4]);

    let sequences: Vec<(&str, String)> = cnv_patterns
        .iter()
        .enumerate()
        .map(|(i, seq)| {
            (
                Box::leak(format!("cnv_{}", i).into_boxed_str()) as &str,
                seq.clone(),
            )
        })
        .collect();

    let sequences_ref: Vec<(&str, &str)> = sequences
        .iter()
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
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: true,
        paf: None,
        seqwish_style: false,
        no_sort: false,
            groom: false,
        iterative_groom: None,
        odgi_style_groom: false,
        sort_groom_sort: false,
        sgd_sort: false,
        aligner: "allwave".to_string(),
    };

    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();

    // Should create a graph that represents the tandem repeat structure
    let _nodes: Vec<&str> = gfa_content.lines().filter(|l| l.starts_with('S')).collect();

    // Should have paths for each CNV pattern
    let paths: Vec<&str> = gfa_content.lines().filter(|l| l.starts_with('P')).collect();
    assert_eq!(paths.len(), cnv_patterns.len());
}

#[test]
fn test_y_shaped_breakpoint() {
    println!("\n=== Testing Y-shaped Graph at Breakpoint ===");

    // Create Y-shaped structure:
    // Common prefix -> Branch A
    //               -> Branch B (with inversion)
    let prefix = "ATCGATCGATCG";
    let branch_a = "AAAATTTTCCCC";
    let branch_b_original = "AAAATTTTCCCC";
    let branch_b = apply_inversion(branch_b_original, 4, 8); // Invert TTTT
    let suffix = "GGGGTTTTAAAA";

    let path_a = format!("{}{}{}", prefix, branch_a, suffix);
    let path_b = format!("{}{}{}", prefix, branch_b, suffix);

    let sequences = vec![("path_a", path_a.as_str()), ("path_b", path_b.as_str())];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 4,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: true,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: true,
        paf: None,
        seqwish_style: false,
        no_sort: false,
            groom: false,
        iterative_groom: None,
        odgi_style_groom: false,
        sort_groom_sort: false,
        sgd_sort: false,
        aligner: "allwave".to_string(),
    };

    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();

    // Should create a branching structure
    let edges: Vec<&str> = gfa_content.lines().filter(|l| l.starts_with('L')).collect();

    // With a Y-shape, we expect branching (more edges than a linear graph)
    assert!(edges.len() > 0);
}

#[test]
fn test_nested_inversions() {
    println!("\n=== Testing Nested Inversions ===");

    // Create nested inversion: large inversion containing a smaller inversion
    let seq = "AAAATTTTCCCCGGGGTTTTAAAA";
    // First apply small inversion
    let seq_small_inv = apply_inversion(&seq, 8, 12); // Invert CCCC
                                                      // Then apply large inversion that contains the small one
    let seq_nested = apply_inversion(&seq_small_inv, 4, 20); // Invert TTTT...TTTT

    let sequences = vec![("original", seq), ("nested_inv", seq_nested.as_str())];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 4,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: true,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: true,
        paf: None,
        seqwish_style: false,
        no_sort: false,
            groom: false,
        iterative_groom: None,
        odgi_style_groom: false,
        sort_groom_sort: false,
        sgd_sort: false,
        aligner: "allwave".to_string(),
    };

    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    assert!(gfa_content.contains("VN:Z:1.0"));
}

#[test]
fn test_random_rearrangements() {
    println!("\n=== Testing Random Rearrangements ===");

    let mut rng = StdRng::seed_from_u64(42);
    let base_length = 100;

    // Generate base sequence
    let base_seq = generate_random_dna(&mut rng, base_length);

    // Apply multiple random rearrangements
    let mut rearranged = base_seq.clone();

    // Random inversion
    let inv_start = rng.gen_range(10..40);
    let inv_end = rng.gen_range(inv_start + 10..60);
    rearranged = apply_inversion(&rearranged, inv_start, inv_end);

    // Random duplication
    let dup_start = rng.gen_range(60..80);
    let dup_end = rng.gen_range(dup_start + 5..90);
    let dup_target = rng.gen_range(0..20);
    rearranged = apply_duplication(&rearranged, dup_start, dup_end, dup_target);

    let sequences = vec![
        ("original", base_seq.as_str()),
        ("rearranged", rearranged.as_str()),
    ];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 5,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: true,
        paf: None,
        seqwish_style: false,
        no_sort: false,
            groom: false,
        iterative_groom: None,
        odgi_style_groom: false,
        sort_groom_sort: false,
        sgd_sort: false,
        aligner: "allwave".to_string(),
    };

    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    assert!(gfa_content.contains("VN:Z:1.0"));
}

#[test]
fn test_tandem_repeat_polymorphism() {
    println!("\n=== Testing Tandem Repeat Polymorphism ===");

    let repeat_unit = "CAG";
    let prefix = "ATCGATCG";
    let suffix = "TAGCTAGT";

    // Different copy numbers (like Huntington's disease CAG repeats)
    let copy_numbers = vec![10, 15, 20, 35, 50];

    let sequences: Vec<(String, String)> = copy_numbers
        .iter()
        .map(|&n| {
            let id = format!("repeat_{}", n);
            let seq = format!(
                "{}{}{}",
                prefix,
                generate_tandem_repeat(repeat_unit, n),
                suffix
            );
            (id, seq)
        })
        .collect();

    let sequences_ref: Vec<(&str, &str)> = sequences
        .iter()
        .map(|(id, seq)| (id.as_str(), seq.as_str()))
        .collect();

    let fasta = create_test_fasta(&sequences_ref);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 3,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: true,
        paf: None,
        seqwish_style: false,
        no_sort: false,
            groom: false,
        iterative_groom: None,
        odgi_style_groom: false,
        sort_groom_sort: false,
        sgd_sort: false,
        aligner: "allwave".to_string(),
    };

    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();

    // Should handle variable number tandem repeats
    let paths: Vec<&str> = gfa_content.lines().filter(|l| l.starts_with('P')).collect();
    assert_eq!(paths.len(), copy_numbers.len());
}

#[test]
fn test_complex_nested_variations() {
    println!("\n=== Testing Complex Nested Variations ===");

    // Create a complex sequence with multiple variation types
    let base = "AAAATTTTCCCCGGGGAAAATTTTCCCCGGGG";

    // Apply variations in layers
    // 1. Tandem duplication of CCCC
    let dup_seq = base.replace("CCCC", "CCCCCCCC");

    // 2. Inversion of part of the sequence including the duplication
    let inv_seq = apply_inversion(&dup_seq, 8, 24);

    // 3. Add some SNPs
    let mut snp_seq = inv_seq.chars().collect::<Vec<_>>();
    snp_seq[5] = 'G'; // A->G
    snp_seq[15] = 'A'; // ?->A
    let final_seq: String = snp_seq.into_iter().collect();

    let sequences = vec![("original", base), ("complex_var", final_seq.as_str())];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 4,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: Some(0.2),
        verbose: true,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: true,
        paf: None,
        seqwish_style: false,
        no_sort: false,
            groom: false,
        iterative_groom: None,
        odgi_style_groom: false,
        sort_groom_sort: false,
        sgd_sort: false,
        aligner: "allwave".to_string(),
    };

    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    assert!(gfa_content.contains("VN:Z:1.0"));
}

#[test]
fn test_interspersed_repeats_with_inversions() {
    println!("\n=== Testing Interspersed Repeats with Inversions ===");

    let repeat = "ATCGATCG";
    let spacer = "TTTT";

    // Create interspersed repeat
    let normal = generate_interspersed_repeat(repeat, spacer, 4);

    // Create version with some repeats inverted
    let mut parts: Vec<String> = Vec::new();
    for i in 0..4 {
        if i % 2 == 0 {
            parts.push(repeat.to_string());
        } else {
            parts.push(String::from_utf8(reverse_complement(repeat.as_bytes())).unwrap());
        }
        if i < 3 {
            parts.push(spacer.to_string());
        }
    }
    let inverted = parts.join("");

    let sequences = vec![
        ("normal_repeat", normal.as_str()),
        ("inv_repeat", inverted.as_str()),
    ];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 4,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: true,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: true,
        paf: None,
        seqwish_style: false,
        no_sort: false,
            groom: false,
        iterative_groom: None,
        odgi_style_groom: false,
        sort_groom_sort: false,
        sgd_sort: false,
        aligner: "allwave".to_string(),
    };

    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();

    // Should handle alternating orientations
    let paths: Vec<&str> = gfa_content.lines().filter(|l| l.starts_with('P')).collect();
    assert_eq!(paths.len(), 2);
}

#[test]
fn test_breakpoint_resolution() {
    println!("\n=== Testing Breakpoint Resolution ===");

    // Simulate a translocation breakpoint
    let chr1_left = "AAAATTTTCCCC";
    let chr1_right = "GGGGAAAATTTT";
    let chr2_left = "CCCCGGGGTTTT";
    let chr2_right = "AAAACCCCGGGG";

    // Normal sequences
    let normal1 = format!("{}{}", chr1_left, chr1_right);
    let normal2 = format!("{}{}", chr2_left, chr2_right);

    // Translocation: chr1_left + chr2_right and chr2_left + chr1_right
    let trans1 = format!("{}{}", chr1_left, chr2_right);
    let trans2 = format!("{}{}", chr2_left, chr1_right);

    // Translocation with inversion at breakpoint
    let trans_inv = format!(
        "{}{}",
        chr1_left,
        String::from_utf8(reverse_complement(chr2_right.as_bytes())).unwrap()
    );

    let sequences = vec![
        ("normal_chr1", normal1.as_str()),
        ("normal_chr2", normal2.as_str()),
        ("trans_der1", trans1.as_str()),
        ("trans_der2", trans2.as_str()),
        ("trans_inv", trans_inv.as_str()),
    ];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 4,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: true,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: true,
        paf: None,
        seqwish_style: false,
        no_sort: false,
            groom: false,
        iterative_groom: None,
        odgi_style_groom: false,
        sort_groom_sort: false,
        sgd_sort: false,
        aligner: "allwave".to_string(),
    };

    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();

    // Should create a complex graph representing all breakpoint combinations
    let nodes: Vec<&str> = gfa_content.lines().filter(|l| l.starts_with('S')).collect();

    let edges: Vec<&str> = gfa_content.lines().filter(|l| l.starts_with('L')).collect();

    let paths: Vec<&str> = gfa_content.lines().filter(|l| l.starts_with('P')).collect();

    // Should create a graph representing all breakpoint combinations
    assert_eq!(paths.len(), 5); // All 5 sequences should have paths
    assert!(nodes.len() > 0); // Should have nodes
    assert!(edges.len() > 0); // Should have edges
}
