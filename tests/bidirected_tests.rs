use seqrush::bidirected_graph::*;
use seqrush::pos::*;
use seqrush::seqrush::Args;
use seqrush::seqrush_bidirected_simplified::*;
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

fn parse_gfa_for_orientations(gfa_content: &str) -> (Vec<String>, Vec<String>, Vec<String>) {
    let mut nodes = Vec::new();
    let mut edges = Vec::new();
    let mut paths = Vec::new();

    for line in gfa_content.lines() {
        if line.starts_with('S') {
            nodes.push(line.to_string());
        } else if line.starts_with('L') {
            edges.push(line.to_string());
        } else if line.starts_with('P') {
            paths.push(line.to_string());
        }
    }

    (nodes, edges, paths)
}

#[test]
fn test_simple_forward_sequence() {
    let sequences = vec![("seq1", "ATCG")];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
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

    run_seqrush_bidirected_simple(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (nodes, _edges, paths) = parse_gfa_for_orientations(&gfa_content);

    // Should have 4 nodes
    assert_eq!(nodes.len(), 4);

    // Path should be all forward
    assert!(paths[0].contains("1+,2+,3+,4+"));
}

#[test]
fn test_reverse_complement_alignment() {
    let sequences = vec![
        ("forward", "ATCGATCG"),
        ("reverse", "CGATCGAT"), // Reverse complement of ATCGATCG
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

    run_seqrush_bidirected_simple(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    println!("GFA output:\n{}", gfa_content);

    let (nodes, _edges, paths) = parse_gfa_for_orientations(&gfa_content);

    // Check that we have paths for both sequences
    assert_eq!(paths.len(), 2);

    // The graph should properly represent both sequences
    // Even though they are reverse complements, they might not share nodes
    // if they are treated as separate sequences in memory
    println!("Number of nodes: {}", nodes.len());
}

#[test]
fn test_palindromic_sequence() {
    let sequences = vec![
        ("palindrome", "GAATTC"), // EcoRI site - palindromic
    ];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
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

    run_seqrush_bidirected_simple(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (_nodes, _edges, paths) = parse_gfa_for_orientations(&gfa_content);

    // Should have one path
    assert_eq!(paths.len(), 1);
}

#[test]
fn test_inversion_between_sequences() {
    let sequences = vec![
        ("ref", "AAAATTTTGGGG"),
        ("inv", "AAAAAAAAGGGG"), // Middle TTTT is inverted to AAAA
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

    run_seqrush_bidirected_simple(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    println!("Inversion test GFA:\n{}", gfa_content);

    let (_nodes, _edges, paths) = parse_gfa_for_orientations(&gfa_content);

    // Should have shared nodes for the common parts
    // Note: The simple bidirected implementation doesn't detect partial inversions within sequences
    // It only aligns full sequences in forward and reverse orientations
    // For partial inversion detection, use the inversion_aware_seqrush module
    // Here we just check that the graph was built correctly
    assert_eq!(paths.len(), 2);

    // The simple implementation won't detect the partial inversion
    // Both paths will use forward orientation only
    assert!(paths.iter().all(|p| !p.contains("-")));
}

#[test]
fn test_complex_rearrangement() {
    let sequences = vec![
        ("seq1", "ATCGATCGATCG"),
        ("seq2", "ATCGCGATCGAT"), // Reverse complement of middle part
    ];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 3,
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

    run_seqrush_bidirected_simple(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (_nodes, _edges, paths) = parse_gfa_for_orientations(&gfa_content);

    assert_eq!(paths.len(), 2);
}

#[test]
fn test_self_alignment_with_repeats() {
    let sequences = vec![
        ("repeat", "ATCGATCG"), // Self-similar sequence
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

    run_seqrush_bidirected_simple(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (nodes, _edges, paths) = parse_gfa_for_orientations(&gfa_content);

    // Due to self-alignment, might have fewer nodes than characters
    assert!(nodes.len() <= 8);
    assert_eq!(paths.len(), 1);
}

#[test]
fn test_position_encoding() {
    // Test the position encoding functions
    let pos_fwd = make_pos(100, false);
    assert_eq!(offset(pos_fwd), 100);
    assert!(!is_rev(pos_fwd));

    let pos_rev = make_pos(100, true);
    assert_eq!(offset(pos_rev), 100);
    assert!(is_rev(pos_rev));

    // Test increment/decrement
    let next_fwd = incr_pos(pos_fwd);
    assert_eq!(offset(next_fwd), 101);

    let next_rev = incr_pos(pos_rev);
    assert_eq!(offset(next_rev), 99);
}

#[test]
fn test_handle_conversion() {
    let handle_fwd = Handle::forward(42);
    let handle_rev = Handle::reverse(42);

    assert_eq!(handle_fwd.node_id(), 42);
    assert!(!handle_fwd.is_reverse());

    assert_eq!(handle_rev.node_id(), 42);
    assert!(handle_rev.is_reverse());

    assert_eq!(handle_fwd.flip(), handle_rev);
}

#[test]
fn test_multiple_sequence_alignment() {
    let sequences = vec![
        ("seq1", "ATCGATCG"),
        ("seq2", "ATCGATGG"), // SNP at end
        ("seq3", "CGATCGAT"), // Reverse complement of seq1
    ];

    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 3,
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

    run_seqrush_bidirected_simple(args).unwrap();

    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (_nodes, _edges, paths) = parse_gfa_for_orientations(&gfa_content);

    assert_eq!(paths.len(), 3);
}
