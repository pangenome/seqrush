#[cfg(test)]
mod tests {
    use crate::bidirected_graph::*;
    use crate::bidirected_ops::*;
    
    #[test]
    fn test_simple_reverse_complement_path() {
        let mut graph = BidirectedGraph::new();
        
        // Create a simple sequence: ATCG
        graph.add_node(1, b"A".to_vec());
        graph.add_node(2, b"T".to_vec());
        graph.add_node(3, b"C".to_vec());
        graph.add_node(4, b"G".to_vec());
        
        // Add forward edges
        graph.add_edge(Handle::forward(1), Handle::forward(2));
        graph.add_edge(Handle::forward(2), Handle::forward(3));
        graph.add_edge(Handle::forward(3), Handle::forward(4));
        
        // Path 1: Forward traversal (ATCG)
        graph.build_path("forward".to_string(), vec![
            (1, false), (2, false), (3, false), (4, false)
        ]);
        
        // Path 2: Reverse complement traversal (CGAT)
        // We traverse the nodes in reverse order with reverse orientation
        graph.build_path("reverse_complement".to_string(), vec![
            (4, true), (3, true), (2, true), (1, true)
        ]);
        
        // Get sequences
        let get_node = |id: usize| graph.nodes.get(&id);
        
        let forward_seq = graph.paths[0].get_sequence(get_node);
        let rc_seq = graph.paths[1].get_sequence(get_node);
        
        assert_eq!(forward_seq, b"ATCG");
        assert_eq!(rc_seq, b"CGAT");
        assert_eq!(rc_seq, reverse_complement(&forward_seq));
    }
    
    #[test]
    fn test_palindromic_sequence() {
        let mut graph = BidirectedGraph::new();
        
        // Create a palindromic sequence: GAATTC (EcoRI site)
        graph.add_node(1, b"G".to_vec());
        graph.add_node(2, b"A".to_vec());
        graph.add_node(3, b"A".to_vec());
        graph.add_node(4, b"T".to_vec());
        graph.add_node(5, b"T".to_vec());
        graph.add_node(6, b"C".to_vec());
        
        // Add edges
        for i in 1..6 {
            graph.add_edge(Handle::forward(i), Handle::forward(i + 1));
        }
        
        // Forward path
        graph.build_path("forward".to_string(), vec![
            (1, false), (2, false), (3, false), (4, false), (5, false), (6, false)
        ]);
        
        // The reverse complement of GAATTC is also GAATTC
        // But traversed from the opposite direction
        graph.build_path("palindrome_rc".to_string(), vec![
            (6, true), (5, true), (4, true), (3, true), (2, true), (1, true)
        ]);
        
        let get_node = |id: usize| graph.nodes.get(&id);
        
        let forward_seq = graph.paths[0].get_sequence(get_node);
        let rc_seq = graph.paths[1].get_sequence(get_node);
        
        assert_eq!(forward_seq, b"GAATTC");
        assert_eq!(rc_seq, b"GAATTC"); // Palindrome property
    }
    
    #[test]
    fn test_mixed_orientation_path() {
        let mut graph = BidirectedGraph::new();
        
        // Create nodes
        graph.add_node(1, b"ATG".to_vec());
        graph.add_node(2, b"CCC".to_vec());
        graph.add_node(3, b"TAG".to_vec());
        
        // Add bidirectional edges
        graph.add_edge(Handle::forward(1), Handle::forward(2));
        graph.add_edge(Handle::forward(2), Handle::reverse(3)); // Forward to reverse
        
        // Path with mixed orientations: ATG -> CCC -> CTA (reverse of TAG)
        graph.build_path("mixed".to_string(), vec![
            (1, false), (2, false), (3, true)
        ]);
        
        let get_node = |id: usize| graph.nodes.get(&id);
        let seq = graph.paths[0].get_sequence(get_node);
        
        assert_eq!(seq, b"ATGCCCCTA");
    }
    
    #[test]
    fn test_self_reverse_complement() {
        let mut graph = BidirectedGraph::new();
        
        // Create a sequence that when reverse complemented maps to itself
        // Example: ACGT on forward = ACGT, on reverse = ACGT
        graph.add_node(1, b"ACGTACGT".to_vec());
        
        // Path that visits the same node twice in different orientations
        graph.build_path("self_rc".to_string(), vec![
            (1, false), // Forward: ACGTACGT
            (1, true),  // Reverse: ACGTACGT (happens to be self-RC)
        ]);
        
        let get_node = |id: usize| graph.nodes.get(&id);
        let seq = graph.paths[0].get_sequence(get_node);
        
        assert_eq!(seq, b"ACGTACGTACGTACGT");
    }
    
    #[test]
    fn test_inversion_representation() {
        let mut graph = BidirectedGraph::new();
        
        // Reference sequence: AAAA[TTTT]GGGG
        // Inverted sequence:  AAAA[AAAA]GGGG (TTTT -> AAAA by inversion)
        
        graph.add_node(1, b"AAAA".to_vec());
        graph.add_node(2, b"TTTT".to_vec());
        graph.add_node(3, b"GGGG".to_vec());
        
        // Add edges for both orientations
        graph.add_edge(Handle::forward(1), Handle::forward(2));
        graph.add_edge(Handle::forward(2), Handle::forward(3));
        graph.add_edge(Handle::forward(1), Handle::reverse(2)); // For inversion
        graph.add_edge(Handle::reverse(2), Handle::forward(3)); // For inversion
        
        // Reference path
        graph.build_path("reference".to_string(), vec![
            (1, false), (2, false), (3, false)
        ]);
        
        // Inverted path (middle segment is reverse complemented)
        graph.build_path("inverted".to_string(), vec![
            (1, false), (2, true), (3, false)
        ]);
        
        let get_node = |id: usize| graph.nodes.get(&id);
        
        let ref_seq = graph.paths[0].get_sequence(get_node);
        let inv_seq = graph.paths[1].get_sequence(get_node);
        
        assert_eq!(ref_seq, b"AAAATTTTGGGG");
        assert_eq!(inv_seq, b"AAAAAAAAGGGG"); // TTTT becomes AAAA when RC'd
    }
    
    #[test]
    fn test_complex_bidirected_graph() {
        let mut graph = BidirectedGraph::new();
        
        // Create a more complex graph with multiple paths
        graph.add_node(1, b"AT".to_vec());
        graph.add_node(2, b"GC".to_vec());
        graph.add_node(3, b"TA".to_vec());
        graph.add_node(4, b"CG".to_vec());
        
        // Add various edges
        graph.add_edge(Handle::forward(1), Handle::forward(2));
        graph.add_edge(Handle::forward(2), Handle::forward(3));
        graph.add_edge(Handle::forward(3), Handle::forward(4));
        graph.add_edge(Handle::forward(1), Handle::reverse(3)); // Alternative path
        graph.add_edge(Handle::reverse(3), Handle::forward(4));
        graph.add_edge(Handle::reverse(4), Handle::reverse(2)); // Back edge
        
        // Path 1: Simple forward
        graph.build_path("path1".to_string(), vec![
            (1, false), (2, false), (3, false), (4, false)
        ]);
        
        // Path 2: With inversion
        graph.build_path("path2".to_string(), vec![
            (1, false), (3, true), (4, false)
        ]);
        
        // Path 3: Complex with back-tracking
        graph.build_path("path3".to_string(), vec![
            (1, false), (2, false), (3, false), (4, true), (2, true)
        ]);
        
        let get_node = |id: usize| graph.nodes.get(&id);
        
        assert_eq!(graph.paths[0].get_sequence(get_node), b"ATGCTACG");
        assert_eq!(graph.paths[1].get_sequence(get_node), b"ATTACG"); // AT + TA (RC of TA) + CG
        assert_eq!(graph.paths[2].get_sequence(get_node), b"ATGCTACGGC"); // Forward then backward
    }
    
    #[test]
    fn test_gfa_with_orientations() {
        let mut graph = BidirectedGraph::new();
        
        graph.add_node(1, b"AAA".to_vec());
        graph.add_node(2, b"TTT".to_vec());
        
        // Add both orientations of edges
        graph.add_edge(Handle::forward(1), Handle::forward(2));
        graph.add_edge(Handle::forward(1), Handle::reverse(2));
        
        // Paths with different orientations
        graph.build_path("forward_path".to_string(), vec![(1, false), (2, false)]);
        graph.build_path("mixed_path".to_string(), vec![(1, false), (2, true)]);
        
        let mut output = Vec::new();
        graph.write_gfa(&mut output).unwrap();
        let gfa = String::from_utf8(output).unwrap();
        
        // Check nodes
        assert!(gfa.contains("S\t1\tAAA"));
        assert!(gfa.contains("S\t2\tTTT"));
        
        // Check edges with orientations
        assert!(gfa.contains("L\t1\t+\t2\t+\t0M"));
        assert!(gfa.contains("L\t1\t+\t2\t-\t0M"));
        
        // Check paths with orientations
        assert!(gfa.contains("P\tforward_path\t1+,2+\t*"));
        assert!(gfa.contains("P\tmixed_path\t1+,2-\t*"));
    }
    
    #[test]
    fn test_handle_edge_of_cases() {
        // Test maximum node ID
        let h1 = Handle::new(usize::MAX >> 1, false);
        assert_eq!(h1.node_id(), usize::MAX >> 1);
        assert!(!h1.is_reverse());
        
        // Test flipping
        let h2 = h1.flip();
        assert_eq!(h2.node_id(), usize::MAX >> 1);
        assert!(h2.is_reverse());
        
        // Test display
        assert_eq!(format!("{}", Handle::forward(42)), "42+");
        assert_eq!(format!("{}", Handle::reverse(42)), "42-");
    }
}