use seqrush::seqrush::{load_sequences, run_seqrush, Args};
use std::io::Write;
use std::time::Instant;
use tempfile::NamedTempFile;

#[test]
fn test_performance_scaling() {
    // Test that the algorithm scales reasonably with sequence size
    let sizes = vec![100, 500, 1000, 2000];
    let mut times = Vec::new();

    for size in sizes {
        // Generate random sequences
        let seq1: String = (0..size)
            .map(|i| match i % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            })
            .collect();

        let mut seq2 = seq1.clone();
        // Add some variations
        for i in (0..size).step_by(100) {
            if i < seq2.len() {
                seq2.replace_range(i..i + 1, "G");
            }
        }

        let mut fasta = NamedTempFile::new().unwrap();
        writeln!(fasta, ">seq1\n{}", seq1).unwrap();
        writeln!(fasta, ">seq2\n{}", seq2).unwrap();
        fasta.flush().unwrap();

        let output = NamedTempFile::new().unwrap();

        let args = Args {
            sequences: fasta.path().to_str().unwrap().to_string(),
            output: output.path().to_str().unwrap().to_string(),
            threads: 1,
            min_match_length: 15,
            scores: "0,5,8,2,24,1".to_string(),
            orientation_scores: "0,1,1,1".to_string(),
            max_divergence: None,
            verbose: false,
            test_mode: true,
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
        sgd_iter_max: 100,
        aligner: "allwave".to_string(),
        };

        let start = Instant::now();
        run_seqrush(args).unwrap();
        let duration = start.elapsed();

        times.push((size, duration.as_millis()));
        println!("Size {}: {} ms", size, duration.as_millis());
    }

    // Check that time doesn't explode (should be roughly O(n²) for pairwise alignment)
    // Allow for some variance but ensure it's not exponential
    if times.len() >= 2 {
        // Find first non-zero time measurement to use as baseline
        let first_nonzero = times.iter().find(|(_, t)| *t > 0);

        if let Some((base_size, base_time)) = first_nonzero {
            // Find the last measurement that is also non-zero and different from baseline
            let last_nonzero = times
                .iter()
                .rev()
                .find(|(size, t)| *t > 0 && *size != *base_size);

            if let Some((last_size, last_time)) = last_nonzero {
                // Only check scaling if we have meaningful measurements and different sizes
                if *last_time > 0 && *base_time > 0 && *last_size != *base_size {
                    let ratio = *last_time as f64 / *base_time as f64;
                    let size_ratio = *last_size as f64 / *base_size as f64;

                    // Time should grow no worse than cubic with size
                    assert!(
                        ratio <= size_ratio.powi(3) * 2.0,
                        "Performance scaling is worse than cubic: {} vs {} (with 2x tolerance)",
                        ratio,
                        size_ratio.powi(3)
                    );
                }
            }
        }
    }
}

#[test]
fn test_real_world_example() {
    // Test with a more realistic example: multiple sequences with various mutations
    let base = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";

    let mut fasta = NamedTempFile::new().unwrap();

    // Reference
    writeln!(fasta, ">ref\n{}", base).unwrap();

    // SNP variant
    let mut snp_variant = base.to_string();
    snp_variant.replace_range(20..21, "G");
    writeln!(fasta, ">snp\n{}", snp_variant).unwrap();

    // Deletion variant
    let del_variant = format!("{}{}", &base[..15], &base[20..]);
    writeln!(fasta, ">del\n{}", del_variant).unwrap();

    // Insertion variant
    let ins_variant = format!("{}AAAAA{}", &base[..30], &base[30..]);
    writeln!(fasta, ">ins\n{}", ins_variant).unwrap();

    fasta.flush().unwrap();

    let output = NamedTempFile::new().unwrap();

    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 2,
        min_match_length: 5,
        scores: "0,5,8,2".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: true,
        test_mode: true,
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
        sgd_iter_max: 100,
        aligner: "allwave".to_string(),
    };

    run_seqrush(args).unwrap();

    // Verify the output exists and is valid GFA
    let content = std::fs::read_to_string(output.path()).unwrap();
    assert!(content.starts_with("H\tVN:Z:1.0"));
    assert!(content.contains("S\t")); // Has segments
    assert!(content.contains("P\t")); // Has paths
    assert!(content.contains("L\t")); // Has links
}

#[test]
fn test_command_line_interface() {
    // Test that the basic CLI works
    let mut fasta = NamedTempFile::new().unwrap();
    writeln!(fasta, ">test\nACGTACGT").unwrap();
    fasta.flush().unwrap();

    let output = NamedTempFile::new().unwrap();

    // Test with minimal arguments
    let args = Args {
        sequences: fasta.path().to_str().unwrap().to_string(),
        output: output.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 15,
        scores: "0,5,8,2".to_string(),
        orientation_scores: "0,1,1,1".to_string(),
        max_divergence: None,
        verbose: false,
        test_mode: true,
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
        sgd_iter_max: 100,
        aligner: "allwave".to_string(),
    };

    run_seqrush(args).unwrap();

    // Test that sequences were loaded correctly
    let sequences = load_sequences(fasta.path().to_str().unwrap()).unwrap();
    assert_eq!(sequences.len(), 1);
    assert_eq!(sequences[0].id, "test");
    assert_eq!(sequences[0].data, b"ACGTACGT");
}

#[test]
fn test_parallel_consistency() {
    // Test that results are consistent regardless of thread count
    let mut fasta = NamedTempFile::new().unwrap();
    for i in 0..10 {
        writeln!(fasta, ">seq{}", i).unwrap();
        let seq: String = (0..50)
            .map(|j| match (i + j) % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            })
            .collect();
        writeln!(fasta, "{}", seq).unwrap();
    }
    fasta.flush().unwrap();

    let mut results = Vec::new();

    for threads in [1, 2, 4] {
        let output = NamedTempFile::new().unwrap();

        let args = Args {
            sequences: fasta.path().to_str().unwrap().to_string(),
            output: output.path().to_str().unwrap().to_string(),
            threads,
            min_match_length: 10,
            scores: "0,5,8,2,24,1".to_string(),
            orientation_scores: "0,1,1,1".to_string(),
            max_divergence: None,
            verbose: false,
            test_mode: true,
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
        sgd_iter_max: 100,
        aligner: "allwave".to_string(),
        };

        run_seqrush(args).unwrap();

        let content = std::fs::read_to_string(output.path()).unwrap();
        let node_count = content.lines().filter(|l| l.starts_with("S\t")).count();
        results.push(node_count);
    }

    // All thread counts should produce the same number of nodes
    assert_eq!(results[0], results[1]);
    assert_eq!(results[1], results[2]);
}
