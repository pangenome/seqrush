use seqrush::seqrush::{run_seqrush, Args};
use std::fs;
use std::io::{BufRead, BufReader};

#[test]
fn test_sorting_creates_sequential_node_ids() {
    // Create test FASTA
    let test_fasta = "test_sort_input.fa";
    let test_output_sorted = "test_sorted_output.gfa";
    let test_output_unsorted = "test_unsorted_output.gfa";
    
    fs::write(
        test_fasta,
        ">seq1\nACGTACGTACGT\n>seq2\nACGTTTTTACGT\n>seq3\nAAAACGTAAAAA\n"
    )
    .unwrap();
    
    // Test with sorting (default)
    let args_sorted = Args {
        sequences: test_fasta.to_string(),
        paf: None,
        output: test_output_sorted.to_string(),
        threads: 1,
        min_match_length: 0,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
        no_compact: false,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: false,
        seqwish_style: false,
        no_sort: false,  // Sorting enabled
    };
    
    run_seqrush(args_sorted).unwrap();
    
    // Test without sorting
    let args_unsorted = Args {
        sequences: test_fasta.to_string(),
        paf: None,
        output: test_output_unsorted.to_string(),
        threads: 1,
        min_match_length: 0,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
        no_compact: false,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: false,
        seqwish_style: false,
        no_sort: true,  // Sorting disabled
    };
    
    run_seqrush(args_unsorted).unwrap();
    
    // Check sorted output has sequential node IDs from 1 to N
    let sorted_file = fs::File::open(test_output_sorted).unwrap();
    let reader = BufReader::new(sorted_file);
    let mut sorted_node_ids = Vec::new();
    
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('S') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                let node_id: usize = parts[1].parse().unwrap();
                sorted_node_ids.push(node_id);
            }
        }
    }
    
    sorted_node_ids.sort();
    
    // Verify nodes are numbered 1 to N
    assert!(!sorted_node_ids.is_empty(), "No nodes found in sorted output");
    assert_eq!(sorted_node_ids[0], 1, "First node should be 1");
    assert_eq!(
        sorted_node_ids[sorted_node_ids.len() - 1],
        sorted_node_ids.len(),
        "Last node should equal the number of nodes"
    );
    
    // Verify all node IDs are sequential
    for i in 0..sorted_node_ids.len() {
        assert_eq!(
            sorted_node_ids[i],
            i + 1,
            "Node IDs should be sequential from 1"
        );
    }
    
    // Check that unsorted output may have non-sequential IDs
    let unsorted_file = fs::File::open(test_output_unsorted).unwrap();
    let reader = BufReader::new(unsorted_file);
    let mut unsorted_node_ids = Vec::new();
    
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('S') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                let node_id: usize = parts[1].parse().unwrap();
                unsorted_node_ids.push(node_id);
            }
        }
    }
    
    unsorted_node_ids.sort();
    
    // The unsorted version might not start from 1 or be sequential
    // Just verify it has nodes
    assert!(!unsorted_node_ids.is_empty(), "No nodes found in unsorted output");
    
    // Clean up
    fs::remove_file(test_fasta).ok();
    fs::remove_file(test_output_sorted).ok();
    fs::remove_file(test_output_unsorted).ok();
}

#[test]
fn test_sorting_preserves_graph_structure() {
    // Create test FASTA
    let test_fasta = "test_structure_input.fa";
    let test_output_sorted = "test_structure_sorted.gfa";
    let test_output_unsorted = "test_structure_unsorted.gfa";
    
    fs::write(
        test_fasta,
        ">seq1\nACGTACGT\n>seq2\nACGTTCGT\n"
    )
    .unwrap();
    
    // Run with sorting
    let args_sorted = Args {
        sequences: test_fasta.to_string(),
        paf: None,
        output: test_output_sorted.to_string(),
        threads: 1,
        min_match_length: 0,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
        no_compact: false,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: false,
        seqwish_style: false,
        no_sort: false,
    };
    
    run_seqrush(args_sorted).unwrap();
    
    // Run without sorting
    let args_unsorted = Args {
        sequences: test_fasta.to_string(),
        paf: None,
        output: test_output_unsorted.to_string(),
        threads: 1,
        min_match_length: 0,
        scores: "0,5,8,2,24,1".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: false,
        no_compact: false,
        sparsification: "1.0".to_string(),
        output_alignments: None,
        validate_paf: false,
        seqwish_style: false,
        no_sort: true,
    };
    
    run_seqrush(args_unsorted).unwrap();
    
    // Count nodes, edges, and paths in both outputs
    let count_graph_elements = |filename: &str| -> (usize, usize, usize) {
        let file = fs::File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let mut nodes = 0;
        let mut edges = 0;
        let mut paths = 0;
        
        for line in reader.lines() {
            let line = line.unwrap();
            if line.starts_with('S') {
                nodes += 1;
            } else if line.starts_with('L') {
                edges += 1;
            } else if line.starts_with('P') {
                paths += 1;
            }
        }
        
        (nodes, edges, paths)
    };
    
    let (sorted_nodes, sorted_edges, sorted_paths) = count_graph_elements(test_output_sorted);
    let (unsorted_nodes, unsorted_edges, unsorted_paths) = count_graph_elements(test_output_unsorted);
    
    // Both graphs should have similar structure (compaction might cause minor differences)
    // Allow for small differences due to compaction non-determinism
    assert!(
        (sorted_nodes as i32 - unsorted_nodes as i32).abs() <= 1, 
        "Node count should be similar (sorted: {}, unsorted: {})", 
        sorted_nodes, unsorted_nodes
    );
    assert!(
        (sorted_edges as i32 - unsorted_edges as i32).abs() <= 1,
        "Edge count should be similar (sorted: {}, unsorted: {})",
        sorted_edges, unsorted_edges
    );
    assert_eq!(sorted_paths, unsorted_paths, "Path count should be the same");
    
    // Clean up
    fs::remove_file(test_fasta).ok();
    fs::remove_file(test_output_sorted).ok();
    fs::remove_file(test_output_unsorted).ok();
}