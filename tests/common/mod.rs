use seqrush::seqrush::Args;

/// Create default Args for testing
pub fn default_test_args(sequences: String, output: String) -> Args {
    Args {
        sequences,
        output,
        threads: 1,
        min_match_length: 0,
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
}