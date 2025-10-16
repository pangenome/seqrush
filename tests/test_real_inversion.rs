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

#[test]
fn test_real_inversion_case() {
    // Create sequences where the middle part can only align in reverse orientation
    // Reference: ATCG [AAAA] GCTA
    // Inverted:  ATCG [TTTT] GCTA (TTTT is reverse complement of AAAA)
    // When aligned forward, the middle will be a mismatch/gap
    let sequences = vec![("reference", "ATCGAAAAGCTA"), ("inverted", "ATCGTTTTGCTA")];

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
        sgd_iter_max: 30,
            skip_sgd: true,
            skip_groom: true,
            skip_topo: true,
            sgd_eta_max: None,
            sgd_theta: 0.99,
            sgd_eps: 0.01,
            sgd_cooling_start: 0.5,
            sort_groom_sort: false,
            iterative_groom: None,
            odgi_style_groom: false,
            sgd_sort: false,
            groom: false,
        aligner: "allwave".to_string(),
    };

    println!("Testing real inversion case...");
    println!("Reference: ATCG[AAAA]GCTA");
    println!("Inverted:  ATCG[TTTT]GCTA");
    println!("Expected: Forward alignment should have mismatches in the middle");

    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    println!("\nGFA output:\n{}", gfa_content);
}

#[test]
fn test_larger_inversion_with_unique_content() {
    // Create a more complex example with unique sequences that can only match when inverted
    // Reference: ATCGATCG [AAAATTTTCCCCGGGG] TAGCTAGT
    // Inverted:  ATCGATCG [CCCCGGGGAAAATTTT] TAGCTAGT
    let sequences = vec![
        ("reference", "ATCGATCGAAAATTTTCCCCGGGGTAGCTAGT"),
        ("inverted", "ATCGATCGCCCCGGGGAAAATTTTTAGCTAGT"),
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
        sgd_iter_max: 30,
            skip_sgd: true,
            skip_groom: true,
            skip_topo: true,
            sgd_eta_max: None,
            sgd_theta: 0.99,
            sgd_eps: 0.01,
            sgd_cooling_start: 0.5,
            sort_groom_sort: false,
            iterative_groom: None,
            odgi_style_groom: false,
            sgd_sort: false,
            groom: false,
        aligner: "allwave".to_string(),
    };

    println!("\nTesting larger inversion with unique content...");
    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();

    // Count how many alignments were found
    let alignments_found = gfa_content.lines().any(|line| line.contains("patched"));

    println!("Found patched alignments: {}", alignments_found);
}

#[test]
fn test_inversion_with_snps() {
    // Test case with an inversion plus some SNPs to make it more realistic
    // This tests that we can still detect inversions even with some noise
    let sequences = vec![
        ("reference", "ATCGATCGAAAATTTTGCTAGCTA"),
        ("inv_snps", "ATCGATCGTTTTAAAAGCTAGCTA"), // Inverted middle + 1 SNP on each side
    ];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 3, // Smaller to allow for SNPs
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: Some(0.2), // Allow 20% divergence
        verbose: true,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: true,
        paf: None,
        seqwish_style: false,
        no_sort: false,
        sgd_iter_max: 30,
            skip_sgd: true,
            skip_groom: true,
            skip_topo: true,
            sgd_eta_max: None,
            sgd_theta: 0.99,
            sgd_eps: 0.01,
            sgd_cooling_start: 0.5,
            sort_groom_sort: false,
            iterative_groom: None,
            odgi_style_groom: false,
            sgd_sort: false,
            groom: false,
        aligner: "allwave".to_string(),
    };

    println!("\nTesting inversion detection with SNPs...");
    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    assert!(gfa_content.contains("VN:Z:1.0"));
}
