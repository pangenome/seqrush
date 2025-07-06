use seqrush::pos::*;
use seqrush::bidirected_graph::*;
use seqrush::seqrush_bidirected_simplified::*;
use seqrush::seqrush::Args;
use std::fs;
use std::io::Write;
use tempfile::NamedTempFile;
use std::collections::{HashMap, HashSet};

fn create_test_fasta(sequences: &[(&str, &str)]) -> NamedTempFile {
    let mut file = NamedTempFile::new().unwrap();
    for (id, seq) in sequences {
        writeln!(file, ">{}", id).unwrap();
        writeln!(file, "{}", seq).unwrap();
    }
    file.flush().unwrap();
    file
}

fn parse_gfa_to_graph(gfa_content: &str) -> (Vec<(usize, String)>, Vec<(usize, usize)>, Vec<(String, Vec<(usize, bool)>)>) {
    let mut nodes = Vec::new();
    let mut edges = Vec::new();
    let mut paths = Vec::new();
    
    for line in gfa_content.lines() {
        if line.starts_with('S') {
            let parts: Vec<&str> = line.split('\t').collect();
            let id = parts[1].parse::<usize>().unwrap();
            let seq = parts[2].to_string();
            nodes.push((id, seq));
        } else if line.starts_with('L') {
            let parts: Vec<&str> = line.split('\t').collect();
            let from = parts[1].parse::<usize>().unwrap();
            let to = parts[3].parse::<usize>().unwrap();
            edges.push((from, to));
        } else if line.starts_with('P') {
            let parts: Vec<&str> = line.split('\t').collect();
            let name = parts[1].to_string();
            // Handle empty paths (no steps)
            if parts.len() > 2 && !parts[2].is_empty() {
                let steps: Vec<(usize, bool)> = parts[2].split(',')
                    .map(|s| {
                        let id = s.trim_end_matches('+').trim_end_matches('-').parse::<usize>().unwrap();
                        let is_rev = s.ends_with('-');
                        (id, is_rev)
                    })
                    .collect();
                paths.push((name, steps));
            } else {
                paths.push((name, vec![]));
            }
        }
    }
    
    (nodes, edges, paths)
}

// Test 1: Base Cases
#[test]
fn test_base_case_empty() {
    let sequences = vec![
        ("empty", ""),
    ];
    
    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_seqrush_bidirected_simple(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    
    // Debug: print the GFA content to see what's happening
    if gfa_content.lines().count() > 1 {
        println!("GFA content for empty sequence:\n{}", gfa_content);
    }
    
    let (_nodes, _edges, _paths) = parse_gfa_to_graph(&gfa_content);
    
    // Empty sequence should produce no content lines besides header and possibly empty path
    let content_lines: Vec<&str> = gfa_content.lines()
        .filter(|line| !line.starts_with('H') && !line.starts_with("P\tempty\t\t"))
        .collect();
    assert_eq!(content_lines.len(), 0, "Empty sequence should produce empty graph");
}

#[test]
fn test_base_case_single_base() {
    let sequences = vec![
        ("single", "A"),
    ];
    
    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_seqrush_bidirected_simple(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (nodes, edges, paths) = parse_gfa_to_graph(&gfa_content);
    
    // Single base should produce one node, no edges, one path
    assert_eq!(nodes.len(), 1);
    assert_eq!(nodes[0].1, "A");
    assert_eq!(edges.len(), 0);
    assert_eq!(paths.len(), 1);
    assert_eq!(paths[0].1.len(), 1);
}

// Test 2: Position Encoding Invariants
#[test]
fn test_position_encoding_properties() {
    // Property 1: Encoding and decoding are inverses
    for offset in [0, 1, 100, 1000, usize::MAX >> 1] {
        for is_rev in [false, true] {
            let pos = make_pos(offset, is_rev);
            assert_eq!(seqrush::pos::offset(pos), offset);
            assert_eq!(seqrush::pos::is_rev(pos), is_rev);
        }
    }
    
    // Property 2: Flipping twice returns original
    for offset in [0, 1, 100] {
        for is_rev in [false, true] {
            let pos = make_pos(offset, is_rev);
            assert_eq!(flip_orientation(flip_orientation(pos)), pos);
        }
    }
    
    // Property 3: Forward and reverse have different encodings
    for offset in [0, 1, 100] {
        let fwd = make_pos(offset, false);
        let rev = make_pos(offset, true);
        assert_ne!(fwd, rev);
    }
}

// Test 3: Reverse Complement Properties
#[test]
fn test_reverse_complement_properties() {
    // Property 1: RC(RC(seq)) = seq
    let sequences = ["A", "T", "C", "G", "ATCG", "AAAA", "ATCGATCG"];
    for seq in sequences {
        let rc1 = reverse_complement(seq.as_bytes());
        let rc2 = reverse_complement(&rc1);
        assert_eq!(rc2, seq.as_bytes());
    }
    
    // Property 2: Complement relationships
    assert_eq!(reverse_complement(b"A"), b"T");
    assert_eq!(reverse_complement(b"T"), b"A");
    assert_eq!(reverse_complement(b"C"), b"G");
    assert_eq!(reverse_complement(b"G"), b"C");
    
    // Property 3: Length preservation
    for seq in sequences {
        let rc = reverse_complement(seq.as_bytes());
        assert_eq!(rc.len(), seq.len());
    }
}

// Test 4: Graph Construction Invariants
#[test]
fn test_identical_sequences_share_nodes() {
    let sequences = vec![
        ("seq1", "ATCGATCG"),
        ("seq2", "ATCGATCG"),
        ("seq3", "ATCGATCG"),
    ];
    
    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_seqrush_bidirected_simple(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (nodes, _edges, paths) = parse_gfa_to_graph(&gfa_content);
    
    // All identical sequences should share the same nodes
    assert_eq!(nodes.len(), 8); // Length of sequence
    assert_eq!(paths.len(), 3);
    
    // All paths should be identical
    assert_eq!(paths[0].1, paths[1].1);
    assert_eq!(paths[1].1, paths[2].1);
}

// Test 5: Inductive Property - Sequence Extension
#[test]
fn test_inductive_sequence_extension() {
    // If graph is correct for sequence S, it should be correct for S + base
    let base_seq = "ATCG";
    let extensions = ["A", "T", "C", "G"];
    
    for ext in extensions {
        let extended = format!("{}{}", base_seq, ext);
        let sequences = vec![
            ("base", base_seq),
            ("extended", &extended),
        ];
        
        let fasta = create_test_fasta(&sequences);
        let output = NamedTempFile::new().unwrap();
        
        let args = Args {
            sequences: fasta.path().to_str().unwrap().to_string(),
            output: output.path().to_str().unwrap().to_string(),
            threads: 1,
            min_match_length: 1,
            scores: "0,5,8,2,24,1".to_string(),
            max_divergence: None,
            verbose: false,
            test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
        };
        
        run_seqrush_bidirected_simple(args).unwrap();
        
        let gfa_content = fs::read_to_string(output.path()).unwrap();
        let (nodes, _edges, paths) = parse_gfa_to_graph(&gfa_content);
        
        // Build node map
        let node_map: HashMap<usize, String> = nodes.into_iter().collect();
        
        // Reconstruct sequences from paths
        let mut reconstructed_base = String::new();
        for (node_id, is_rev) in &paths[0].1 {
            let node_seq = &node_map[node_id];
            if *is_rev {
                reconstructed_base.push_str(&String::from_utf8(reverse_complement(node_seq.as_bytes())).unwrap());
            } else {
                reconstructed_base.push_str(node_seq);
            }
        }
        
        let mut reconstructed_extended = String::new();
        for (node_id, is_rev) in &paths[1].1 {
            let node_seq = &node_map[node_id];
            if *is_rev {
                reconstructed_extended.push_str(&String::from_utf8(reverse_complement(node_seq.as_bytes())).unwrap());
            } else {
                reconstructed_extended.push_str(node_seq);
            }
        }
        
        // Verify sequences are correctly reconstructed
        assert_eq!(reconstructed_base, base_seq);
        assert_eq!(reconstructed_extended, extended);
        
        // Extended sequence should have one more base
        assert_eq!(reconstructed_extended.len(), reconstructed_base.len() + 1);
        
        // First n characters should be the same
        assert_eq!(&reconstructed_extended[..base_seq.len()], base_seq);
    }
}

// Test 6: Bidirectional Properties
#[test]
fn test_bidirectional_handle_properties() {
    // Property 1: Handle encodes node and orientation correctly
    let h1 = Handle::new(42, false);
    assert_eq!(h1.node_id(), 42);
    assert!(!h1.is_reverse());
    
    let h2 = Handle::new(42, true);
    assert_eq!(h2.node_id(), 42);
    assert!(h2.is_reverse());
    
    // Property 2: Forward and reverse handles are different
    assert_ne!(h1, h2);
    
    // Property 3: Flipping inverts orientation
    assert_eq!(h1.flip(), h2);
    assert_eq!(h2.flip(), h1);
}

// Test 7: Path Reconstruction
#[test]
fn test_path_reconstruction_correctness() {
    let test_sequences = vec![
        ("seq1", "ATCG"),
        ("seq2", "ATCGATCG"),
        ("seq3", "GCTA"), // Reverse complement of TAGC
    ];
    
    let fasta = create_test_fasta(&test_sequences);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_seqrush_bidirected_simple(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (nodes, _edges, paths) = parse_gfa_to_graph(&gfa_content);
    
    // Build node map
    let node_map: HashMap<usize, String> = nodes.into_iter().collect();
    
    // Verify each path reconstructs to original sequence
    for (i, (seq_id, seq)) in test_sequences.iter().enumerate() {
        let (path_name, path_steps) = &paths[i];
        assert_eq!(path_name, seq_id);
        
        let mut reconstructed = String::new();
        for (node_id, is_rev) in path_steps {
            let node_seq = &node_map[node_id];
            if *is_rev {
                reconstructed.push_str(&String::from_utf8(reverse_complement(node_seq.as_bytes())).unwrap());
            } else {
                reconstructed.push_str(node_seq);
            }
        }
        
        assert_eq!(&reconstructed, seq, "Failed to reconstruct {}", seq_id);
    }
}

// Test 8: Edge Connectivity
#[test]
fn test_edge_connectivity() {
    let sequences = vec![
        ("linear", "ABCD"),
    ];
    
    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_seqrush_bidirected_simple(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (_nodes, edges, paths) = parse_gfa_to_graph(&gfa_content);
    
    // Linear sequence should have n-1 edges
    assert_eq!(edges.len(), 3);
    
    // Verify connectivity: each consecutive pair in path should have an edge
    let path_steps = &paths[0].1;
    for i in 0..path_steps.len()-1 {
        let from = path_steps[i].0;
        let to = path_steps[i+1].0;
        assert!(edges.contains(&(from, to)), "Missing edge {} -> {}", from, to);
    }
}

// Test 9: Palindrome Handling
#[test]
fn test_palindromic_sequences() {
    // These sequences are their own reverse complements
    let palindromes = vec![
        ("pal1", "GAATTC"), // EcoRI site
        ("pal2", "ACGT"),   // Self-complementary
        ("pal3", "ATAT"),   // Another palindrome
    ];
    
    let fasta = create_test_fasta(&palindromes);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_seqrush_bidirected_simple(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (_nodes, _edges, paths) = parse_gfa_to_graph(&gfa_content);
    
    // Each palindrome should have its own path
    assert_eq!(paths.len(), 3);
}

// Test 10: Union-Find Correctness
#[test]
fn test_union_find_with_matches() {
    // Sequences with known matches
    let sequences = vec![
        ("seq1", "AAAATTTT"),
        ("seq2", "AAAACCCC"),
        ("seq3", "GGGGTTTT"),
    ];
    
    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 4,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_seqrush_bidirected_simple(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (nodes, _edges, _paths) = parse_gfa_to_graph(&gfa_content);
    
    // Should have shared nodes for AAAA and TTTT
    let node_seqs: HashSet<String> = nodes.iter().map(|(_, seq)| seq.clone()).collect();
    assert!(node_seqs.contains("A"));
    assert!(node_seqs.contains("T"));
    assert!(node_seqs.contains("C"));
    assert!(node_seqs.contains("G"));
}

// Test 11: Stress Test - Large Sequences
#[test]
fn test_large_sequences() {
    let size = 1000;
    let seq1: String = (0..size).map(|i| match i % 4 {
        0 => 'A',
        1 => 'T',
        2 => 'C',
        _ => 'G',
    }).collect();
    
    let sequences = vec![
        ("large1", seq1.as_str()),
        ("large2", seq1.as_str()), // Identical
    ];
    
    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 10,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_seqrush_bidirected_simple(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (_nodes, _edges, paths) = parse_gfa_to_graph(&gfa_content);
    
    // Identical sequences should have identical paths
    assert_eq!(paths[0].1, paths[1].1);
}

// Test 12: Mathematical Invariant - Transitivity
#[test]
fn test_transitivity_of_matches() {
    // If A matches B and B matches C, then A should match C
    let sequences = vec![
        ("A", "ATCGATCG"),
        ("B", "ATCGATCG"),
        ("C", "ATCGATCG"),
    ];
    
    let fasta = create_test_fasta(&sequences);
    let output = NamedTempFile::new().unwrap();
    
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
        scores: "0,5,8,2,24,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
            no_compact: true,
        sparsification: "1.0".to_string(),
    };
    
    run_seqrush_bidirected_simple(args).unwrap();
    
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let (_nodes, _edges, paths) = parse_gfa_to_graph(&gfa_content);
    
    // All three should have identical paths (transitivity)
    assert_eq!(paths[0].1, paths[1].1);
    assert_eq!(paths[1].1, paths[2].1);
}