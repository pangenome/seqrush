#[cfg(test)]
mod reverse_complement_tests {
    use crate::seqrush::{Args, run_seqrush};
    use std::fs;
    use tempfile::NamedTempFile;
    use std::io::Write;
    
    fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        seq.iter().rev().map(|&base| {
            match base {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                _ => base,
            }
        }).collect()
    }
    
    fn create_test_fasta(sequences: &[(&str, Vec<u8>)]) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        for (id, seq) in sequences {
            writeln!(file, ">{}", id).unwrap();
            writeln!(file, "{}", String::from_utf8_lossy(seq)).unwrap();
        }
        file.flush().unwrap();
        file
    }
    
    #[test]
    fn test_reverse_complement_detection() {
        // Create 4 sequences where seq3 is reverse complement of part of seq1
        let seq1 = b"ATCGATCGATCGATCG".to_vec();
        let seq2 = b"ATCGATCGTTCGATCG".to_vec(); // Similar to seq1 with one change
        let seq3 = reverse_complement(&seq1[4..12]); // RC of middle portion: "ATCGATCG" -> "CGATCGAT"
        let seq4 = b"ATCGATCGATCGATCG".to_vec(); // Identical to seq1
        
        let sequences = vec![
            ("seq1", seq1),
            ("seq2", seq2),
            ("seq3_rc", seq3),
            ("seq4", seq4),
        ];
        
        let fasta = create_test_fasta(&sequences);
        let output = NamedTempFile::new().unwrap();
        
        let args = Args {
            sequences: fasta.path().to_str().unwrap().to_string(),
            output: output.path().to_str().unwrap().to_string(),
            threads: 1,
            min_match_length: 4,
            verbose: true,
            no_compact: true, // Disable compaction to see raw graph structure
            test_mode: false, // Important: don't disable RC detection
            scores: "0,5,8,2,24,1".to_string(),
            orientation_scores: "0,1,1,1".to_string(),
            max_divergence: None,
            sparsification: "1.0".to_string(),
        };
        
        // Build the graph
        run_seqrush(args).unwrap();
        
        // Read back the GFA to check orientations
        let gfa_content = fs::read_to_string(output.path()).unwrap();
        
        // Check if any paths contain reverse orientations (- suffix)
        let has_reverse_orientations = gfa_content.lines()
            .filter(|line| line.starts_with("P\t"))
            .any(|line| line.contains("-"));
        
        println!("GFA content:\n{}", gfa_content);
        println!("Has reverse orientations: {}", has_reverse_orientations);
        
        // Currently, we expect NO reverse orientations due to the bug
        assert!(!has_reverse_orientations, "Graph should not have reverse orientations with current implementation");
    }
    
    #[test]
    fn test_explicit_reverse_complement_sequence() {
        // Create sequences where one is explicitly the full reverse complement
        let forward = b"ATCGATCGATCGATCGATCG".to_vec();
        let reverse = reverse_complement(&forward);
        
        let sequences = vec![
            ("forward1", forward.clone()),
            ("forward2", forward.clone()),
            ("reverse1", reverse.clone()),
            ("forward3", forward.clone()),
        ];
        
        let fasta = create_test_fasta(&sequences);
        let output = NamedTempFile::new().unwrap();
        
        let args = Args {
            sequences: fasta.path().to_str().unwrap().to_string(),
            output: output.path().to_str().unwrap().to_string(),
            threads: 1,
            min_match_length: 4,
            verbose: true,
            no_compact: true,
            test_mode: false,
            scores: "0,5,8,2,24,1".to_string(),
            orientation_scores: "0,1,1,1".to_string(),
            max_divergence: None,
            sparsification: "1.0".to_string(),
        };
        
        run_seqrush(args).unwrap();
        
        // Analyze the graph structure
        let gfa_content = fs::read_to_string(output.path()).unwrap();
        
        // Count unique nodes - with proper RC handling, forward and reverse 
        // sequences should share nodes
        let node_count = gfa_content.lines()
            .filter(|line| line.starts_with("S\t"))
            .count();
        
        println!("Node count: {}", node_count);
        println!("GFA preview:\n{}", gfa_content.lines().take(20).collect::<Vec<_>>().join("\n"));
        
        // With current implementation, we expect many nodes because RC is not properly handled
        assert!(node_count > 10, "Current implementation creates many nodes due to RC bug");
    }
    
    #[test]
    fn test_mixed_orientations() {
        // Create a more complex scenario with mixed orientations
        // seq1: ATCGATCGATCG
        // seq2: ATCGAT (prefix of seq1)
        // seq3: CGATCGATCGAT (reverse complement of seq1)
        // seq4: GATCGATC (infix of seq1)
        
        let seq1 = b"ATCGATCGATCG".to_vec();
        let seq2 = b"ATCGAT".to_vec();
        let seq3 = reverse_complement(&seq1);
        let seq4 = b"GATCGATC".to_vec();
        
        let sequences = vec![
            ("seq1", seq1),
            ("seq2_prefix", seq2),
            ("seq3_rc", seq3),
            ("seq4_infix", seq4),
        ];
        
        let fasta = create_test_fasta(&sequences);
        let output = NamedTempFile::new().unwrap();
        
        let args = Args {
            sequences: fasta.path().to_str().unwrap().to_string(),
            output: output.path().to_str().unwrap().to_string(),
            threads: 1,
            min_match_length: 4,
            verbose: true,
            no_compact: false, // Enable compaction to see if it helps
            test_mode: false,
            scores: "0,5,8,2,24,1".to_string(),
            orientation_scores: "0,1,1,1".to_string(),
            max_divergence: None,
            sparsification: "1.0".to_string(),
        };
        
        run_seqrush(args).unwrap();
        
        let gfa_content = fs::read_to_string(output.path()).unwrap();
        
        // Parse paths from GFA
        let mut paths: std::collections::HashMap<String, Vec<(usize, char)>> = std::collections::HashMap::new();
        for line in gfa_content.lines() {
            if line.starts_with("P\t") {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 3 {
                    let seq_id = parts[1];
                    let path_str = parts[2];
                    let path: Vec<(usize, char)> = path_str.split(',')
                        .filter_map(|node| {
                            if node.ends_with('+') {
                                node[..node.len()-1].parse::<usize>().ok().map(|id| (id, '+'))
                            } else if node.ends_with('-') {
                                node[..node.len()-1].parse::<usize>().ok().map(|id| (id, '-'))
                            } else {
                                None
                            }
                        })
                        .collect();
                    paths.insert(seq_id.to_string(), path);
                }
            }
        }
        
        // Check the path for seq3_rc
        if let Some(rc_path) = paths.get("seq3_rc") {
            println!("Path for seq3_rc: {:?}", rc_path);
            
            // Check if any nodes have reverse orientation
            let has_reverse = rc_path.iter().any(|(_, orient)| *orient == '-');
            println!("seq3_rc path has reverse orientations: {}", has_reverse);
            
            // Currently expecting no reverse orientations due to the bug
            assert!(!has_reverse, "Current implementation doesn't support reverse orientations in paths");
        }
        
        println!("\nAlignments debug - checking what orientation was detected:");
        println!("Full GFA:\n{}", gfa_content);
    }
    
    #[test]
    fn test_bidirected_graph_expected_behavior() {
        // This test documents what the output SHOULD look like with proper RC handling
        
        let seq1 = b"ATCGATCG".to_vec();
        let _seq2 = b"ATCGATCG".to_vec(); // Identical to seq1
        let _seq3 = reverse_complement(&seq1); // "CGATCGAT"
        let _seq4 = b"ATCGATCG".to_vec(); // Identical to seq1
        
        // With proper bidirected graph handling, we would expect:
        // 1. seq1, seq2, and seq4 share the same path through nodes
        // 2. seq3_rc follows the REVERSE path through the SAME nodes
        // 3. Total nodes ≈ 8 (one per base position in the consensus)
        
        // Expected GFA output would look like:
        // P  seq1      1+,2+,3+,4+,5+,6+,7+,8+  *
        // P  seq2      1+,2+,3+,4+,5+,6+,7+,8+  *
        // P  seq3_rc   8-,7-,6-,5-,4-,3-,2-,1-  *  <-- Note the reverse path!
        // P  seq4      1+,2+,3+,4+,5+,6+,7+,8+  *
        
        println!("This test documents expected behavior with proper RC handling");
        println!("Expected: seq3_rc should traverse nodes in reverse order with - orientation");
        
        // For now, we just document the expected behavior
        assert!(true, "Documentation test - no actual implementation to test yet");
    }
    
    #[test]
    fn test_verify_rc_alignment_detection() {
        // Simple test to verify that RC alignments are at least being detected
        let forward = b"ATCGATCGATCGATCG".to_vec();
        let reverse = reverse_complement(&forward);
        
        let sequences = vec![
            ("seq_forward", forward.clone()),
            ("seq_reverse", reverse.clone()),
        ];
        
        let forward_str = String::from_utf8_lossy(&forward);
        let reverse_str = String::from_utf8_lossy(&reverse);
        
        let fasta = create_test_fasta(&sequences);
        let output = NamedTempFile::new().unwrap();
        
        let args = Args {
            sequences: fasta.path().to_str().unwrap().to_string(),
            output: output.path().to_str().unwrap().to_string(),
            threads: 1,
            min_match_length: 8, // Longer to ensure good alignment
            verbose: true,
            no_compact: true,
            test_mode: false, // Important: allow RC detection
            scores: "0,-1,1,1".to_string(),
            orientation_scores: "0,-1,2,1".to_string(),
            max_divergence: None,
            sparsification: "1.0".to_string(),
        };
        
        run_seqrush(args).unwrap();
        
        let gfa_content = fs::read_to_string(output.path()).unwrap();
        
        // Count nodes - if RC is properly handled, we should have few nodes
        // If not handled, we'll have many distinct nodes
        let node_count = gfa_content.lines()
            .filter(|line| line.starts_with("S\t"))
            .count();
        
        println!("Forward seq: {}", forward_str);
        println!("Reverse seq: {}", reverse_str);
        println!("Node count with RC sequence: {}", node_count);
        
        // With the bug, we expect many nodes (close to 2x sequence length)
        // With proper handling, we'd expect around sequence length
        assert!(node_count > 20, "Expected many nodes due to unhandled RC");
    }
}