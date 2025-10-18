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
fn test_inversion_detection_simple() {
    // Create sequences with a clear inversion
    // Reference: AAAA TTTT GGGG
    // Inverted:  AAAA AAAA GGGG (middle segment is inverted)
    let sequences = vec![("reference", "AAAATTTTGGGG"), ("inverted", "AAAAAAAAGGGG")];

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
        frequency: None,
    };

    println!("Running inversion-aware alignment...");
    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    println!("\nGFA output:\n{}", gfa_content);

    // Check that we have the expected structure
    let mut node_count = 0;
    let mut edge_count = 0;
    let mut path_count = 0;

    for line in gfa_content.lines() {
        if line.starts_with('S') {
            node_count += 1;
        } else if line.starts_with('L') {
            edge_count += 1;
        } else if line.starts_with('P') {
            path_count += 1;
        }
    }

    println!(
        "Nodes: {}, Edges: {}, Paths: {}",
        node_count, edge_count, path_count
    );

    assert_eq!(path_count, 2);
    assert!(node_count > 0);
    assert!(edge_count > 0);
}

#[test]
fn test_inversion_detection_with_context() {
    // More realistic example with context around the inversion
    // Reference: ATCGATCG TTTTAAAA GCTAGCTA
    // Inverted:  ATCGATCG AAAATTTT GCTAGCTA (middle inverted)
    let sequences = vec![
        ("reference", "ATCGATCGTTTTAAAAGCTAGCTA"),
        ("inverted", "ATCGATCGAAAATTTTGCTAGCTA"),
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
        frequency: None,
    };

    println!("\nRunning inversion detection with context...");
    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();

    // Verify the graph was built
    assert!(gfa_content.contains("VN:Z:1.0"));

    // Count components
    let s_count = gfa_content.lines().filter(|l| l.starts_with('S')).count();
    let p_count = gfa_content.lines().filter(|l| l.starts_with('P')).count();

    println!("Graph has {} nodes and {} paths", s_count, p_count);
    assert_eq!(p_count, 2);
}

#[test]
fn test_no_inversion() {
    // Control test - sequences with no inversions
    let sequences = vec![("seq1", "ATCGATCGATCG"), ("seq2", "ATCGATCGATCG")];

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
        verbose: false,
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
        frequency: None,
    };

    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();

    // For identical sequences, we should get a simple graph
    let p_count = gfa_content.lines().filter(|l| l.starts_with('P')).count();
    assert_eq!(p_count, 2);
}

#[test]
fn test_multiple_inversions() {
    // Test with multiple small inversions
    // Reference: AAAA TTTT CCCC GGGG
    // Multi-inv: AAAA AAAA GGGG GGGG (two inversions)
    let sequences = vec![
        ("reference", "AAAATTTTCCCCGGGG"),
        ("multi_inv", "AAAAAAAAGGGGCCCC"),
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
        frequency: None,
    };

    println!("\nRunning multiple inversion detection...");
    run_inversion_aware_seqrush(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    println!(
        "Multiple inversions GFA:\n{}",
        &gfa_content[..500.min(gfa_content.len())]
    );

    assert!(gfa_content.contains("VN:Z:1.0"));
}
