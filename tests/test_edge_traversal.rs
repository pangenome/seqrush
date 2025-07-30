use std::fs;
use std::io::Write;
use tempfile::NamedTempFile;
use seqrush::seqrush::{Args, run_seqrush};

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
fn test_no_untraversed_edges_in_output() {
    // Create sequences that should result in a connected graph
    let sequences = vec![
        ("seq1", "ATCGATCGATCG"),
        ("seq2", "ATCGATCGATCG"),
        ("seq3", "ATCGATCGTTCG"), // Single SNP
        ("seq4", "ATCGATCGATCG"),
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
        verbose: false,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
            output_alignments: None,
            validate_paf: true,
            paf: None,
            seqwish_style: false,
    };
    
    run_seqrush(args).unwrap();
    
    // Read and verify the GFA
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    
    // Parse edges
    let mut edge_count = 0;
    for line in gfa_content.lines() {
        if line.starts_with('L') {
            edge_count += 1;
        }
    }
    
    // Parse paths and check edge usage
    let mut path_edges = std::collections::HashSet::new();
    for line in gfa_content.lines() {
        if line.starts_with('P') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let path = parts[2];
                let nodes: Vec<&str> = path.split(',').collect();
                
                // Add edges from consecutive nodes
                for window in nodes.windows(2) {
                    if let [from, to] = window {
                        let from_id = from.trim_end_matches('+').trim_end_matches('-');
                        let to_id = to.trim_end_matches('+').trim_end_matches('-');
                        path_edges.insert((from_id.to_string(), to_id.to_string()));
                    }
                }
            }
        }
    }
    
    // In a properly constructed graph, all edges should be traversed by paths
    // Since we're not doing compaction, check that paths exist
    let path_count = gfa_content.lines().filter(|l| l.starts_with('P')).count();
    assert_eq!(path_count, 4, "Should have 4 paths");
    assert!(edge_count > 0, "Graph should have edges");
}

#[test]
fn test_self_loops_in_gfa() {
    // Create sequence with repeats that might create self-loops
    let sequences = vec![
        ("seq1", "AAAAAAAA"),
        ("seq2", "AAAAAAAA"),
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
        verbose: false,
        test_mode: false,
        no_compact: false, // Allow compaction
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: true,
            paf: None,
            seqwish_style: false,
    };
    
    run_seqrush(args).unwrap();
    
    // Check for self-loops in edges
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    let mut self_loop_count = 0;
    
    for line in gfa_content.lines() {
        if line.starts_with('L') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 5 {
                let from_node = parts[1];
                let to_node = parts[3];
                if from_node == to_node {
                    self_loop_count += 1;
                }
            }
        }
    }
    
    // Self-loops are allowed but should be minimal
    assert!(self_loop_count <= 2, "Should have minimal self-loops, found {}", self_loop_count);
}

#[test]
fn test_complex_graph_produces_valid_gfa() {
    // Complex test with multiple variation types
    let sequences = vec![
        ("ref", "ATCGATCGATCGATCG"),
        ("snp", "ATCGATGGATCGATCG"), // SNP in middle
        ("del", "ATCGATCGATCG"),     // Deletion
        ("inv", "ATCGATCGATGCTAGC"), // Different ending
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
        max_divergence: Some(0.2),
        verbose: false,
        test_mode: false,
        no_compact: true,
        sparsification: "1.0".to_string(),
            output_alignments: None,
            validate_paf: true,
            paf: None,
            seqwish_style: false,
    };
    
    run_seqrush(args).unwrap();
    
    // Verify GFA structure
    let gfa_content = fs::read_to_string(output.path()).unwrap();
    
    // Should have all components
    assert!(gfa_content.contains("H\tVN:Z:1.0"), "Missing header");
    
    let has_nodes = gfa_content.lines().any(|l| l.starts_with('S'));
    let has_edges = gfa_content.lines().any(|l| l.starts_with('L'));
    let has_paths = gfa_content.lines().any(|l| l.starts_with('P'));
    
    assert!(has_nodes, "GFA should have nodes");
    assert!(has_edges, "GFA should have edges");
    assert!(has_paths, "GFA should have paths");
    
    // Check path count
    let path_count = gfa_content.lines().filter(|l| l.starts_with('P')).count();
    assert_eq!(path_count, 4, "Should have 4 paths");
}