use std::fs;
use std::io::Write;
/// Integration tests for reverse complement alignment handling
/// These tests verify that sequences aligned in reverse complement
/// are properly integrated into the graph and not left as isolated nodes
use std::process::Command;
use tempfile::TempDir;

fn count_nodes_in_gfa(gfa_content: &str) -> usize {
    gfa_content
        .lines()
        .filter(|line| line.starts_with('S'))
        .count()
}

fn count_isolated_paths_in_gfa(gfa_content: &str) -> usize {
    // Count paths that consist of a single node (indicating isolation)
    gfa_content
        .lines()
        .filter(|line| line.starts_with('P'))
        .filter(|line| {
            // Extract the path part (third field)
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                // Count the number of nodes in the path by counting + or -
                let path_nodes = parts[2].split(',').count();
                path_nodes == 1
            } else {
                false
            }
        })
        .count()
}

#[test]
fn test_rc_alignment_creates_shared_nodes() {
    let temp_dir = TempDir::new().unwrap();

    // Create test sequences - forward and reverse complement
    let fasta_path = temp_dir.path().join("test_rc.fa");
    let mut fasta_file = fs::File::create(&fasta_path).unwrap();
    writeln!(fasta_file, ">seq1").unwrap();
    writeln!(fasta_file, "ATCGATCGATCG").unwrap();
    writeln!(fasta_file, ">seq2_rc").unwrap();
    writeln!(fasta_file, "CGATCGATCGAT").unwrap(); // RC of seq1

    // Create PAF with RC alignment
    let paf_path = temp_dir.path().join("test_rc.paf");
    let mut paf_file = fs::File::create(&paf_path).unwrap();
    // seq1 aligns to seq2_rc in reverse complement with perfect match
    writeln!(
        paf_file,
        "seq1\t12\t0\t12\t+\tseq2_rc\t12\t0\t12\t12\t12\t60\tcg:Z:12="
    )
    .unwrap();
    // Self-alignments
    writeln!(
        paf_file,
        "seq1\t12\t0\t12\t+\tseq1\t12\t0\t12\t12\t12\t60\tcg:Z:12="
    )
    .unwrap();
    writeln!(
        paf_file,
        "seq2_rc\t12\t0\t12\t+\tseq2_rc\t12\t0\t12\t12\t12\t60\tcg:Z:12="
    )
    .unwrap();
    // RC alignment - this is the key one
    writeln!(
        paf_file,
        "seq2_rc\t12\t0\t12\t-\tseq1\t12\t0\t12\t12\t12\t60\tcg:Z:12="
    )
    .unwrap();

    let gfa_path = temp_dir.path().join("test_rc.gfa");

    // Run seqrush (it will do all-vs-all alignment automatically)
    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "seqrush", "--"])
        .args(&["-s", fasta_path.to_str().unwrap()])
        .args(&["-o", gfa_path.to_str().unwrap()])
        .args(&["-k", "0"])
        .output()
        .expect("Failed to run seqrush");

    if !output.status.success() {
        eprintln!("STDOUT: {}", String::from_utf8_lossy(&output.stdout));
        eprintln!("STDERR: {}", String::from_utf8_lossy(&output.stderr));
        panic!("seqrush failed");
    }

    // Read and analyze the GFA
    let gfa_content = fs::read_to_string(&gfa_path).unwrap();
    let node_count = count_nodes_in_gfa(&gfa_content);
    let isolated_paths = count_isolated_paths_in_gfa(&gfa_content);

    println!("Node count: {}", node_count);
    println!("Isolated paths: {}", isolated_paths);

    // With proper RC handling, sequences should share nodes
    assert!(
        node_count < 24,
        "Sequences should share nodes, not have 24 separate nodes"
    );
    assert_eq!(isolated_paths, 0, "No sequence should be isolated");
}

#[test]
fn test_complex_rc_alignment_scenario() {
    let temp_dir = TempDir::new().unwrap();

    // Create test with 3 sequences: one forward, one RC, one partial match
    let fasta_path = temp_dir.path().join("test_complex_rc.fa");
    let mut fasta_file = fs::File::create(&fasta_path).unwrap();
    writeln!(fasta_file, ">forward").unwrap();
    writeln!(fasta_file, "AAAATTTTGGGGCCCC").unwrap();
    writeln!(fasta_file, ">reverse").unwrap();
    writeln!(fasta_file, "GGGGCCCCAAAATTTT").unwrap(); // RC of forward
    writeln!(fasta_file, ">partial").unwrap();
    writeln!(fasta_file, "AAAATTTTAAAATTTT").unwrap(); // Matches first half of forward

    let gfa_path = temp_dir.path().join("test_complex_rc.gfa");

    // Run seqrush without PAF (let it do all-vs-all alignment)
    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "seqrush", "--"])
        .args(&["-s", fasta_path.to_str().unwrap()])
        .args(&["-o", gfa_path.to_str().unwrap()])
        .args(&["-k", "0"])
        .output()
        .expect("Failed to run seqrush");

    if !output.status.success() {
        eprintln!("STDOUT: {}", String::from_utf8_lossy(&output.stdout));
        eprintln!("STDERR: {}", String::from_utf8_lossy(&output.stderr));
        panic!("seqrush failed");
    }

    // Read and analyze the GFA
    let gfa_content = fs::read_to_string(&gfa_path).unwrap();
    let isolated_paths = count_isolated_paths_in_gfa(&gfa_content);

    // No sequence should be completely isolated
    assert_eq!(
        isolated_paths, 0,
        "No sequence should be isolated, even RC ones"
    );

    // Verify all sequences have paths
    assert!(
        gfa_content.contains("P\tforward\t"),
        "Forward sequence should have a path"
    );
    assert!(
        gfa_content.contains("P\treverse\t"),
        "Reverse sequence should have a path"
    );
    assert!(
        gfa_content.contains("P\tpartial\t"),
        "Partial sequence should have a path"
    );
}

#[test]
fn test_palindromic_sequence_handling() {
    let temp_dir = TempDir::new().unwrap();

    // Create palindromic sequences (same forward and reverse)
    let fasta_path = temp_dir.path().join("test_palindrome.fa");
    let mut fasta_file = fs::File::create(&fasta_path).unwrap();
    writeln!(fasta_file, ">palindrome1").unwrap();
    writeln!(fasta_file, "AATTGGCCGGCCAATT").unwrap(); // Palindromic
    writeln!(fasta_file, ">palindrome2").unwrap();
    writeln!(fasta_file, "AATTGGCCGGCCAATT").unwrap(); // Same palindrome

    let gfa_path = temp_dir.path().join("test_palindrome.gfa");

    // Run seqrush
    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "seqrush", "--"])
        .args(&["-s", fasta_path.to_str().unwrap()])
        .args(&["-o", gfa_path.to_str().unwrap()])
        .args(&["-k", "0"])
        .output()
        .expect("Failed to run seqrush");

    assert!(
        output.status.success(),
        "seqrush should handle palindromic sequences"
    );

    // Read and verify the GFA
    let gfa_content = fs::read_to_string(&gfa_path).unwrap();
    let node_count = count_nodes_in_gfa(&gfa_content);

    // Identical palindromic sequences should collapse to the same nodes
    assert!(
        node_count <= 16,
        "Identical palindromes should share all nodes"
    );
}
