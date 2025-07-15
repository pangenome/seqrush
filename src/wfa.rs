use lib_wfa2::affine_wavefront::{
    AffineWavefronts, AlignmentScope, AlignmentSpan, MemoryMode, HeuristicStrategy
};

/// Convert WFA2 CIGAR bytes to standard CIGAR string
/// WFA2 uses opposite convention for I/D operations:
/// - WFA2: I = consume reference, D = consume query  
/// - Standard: I = consume query, D = consume reference
pub fn cigar_bytes_to_string(cigar_bytes: &[u8]) -> String {
    let mut cigar_str = String::new();
    let mut i = 0;

    while i < cigar_bytes.len() {
        let op = cigar_bytes[i];
        let mut count = 1;
        let mut j = i + 1;

        // Count consecutive same operations
        while j < cigar_bytes.len() && cigar_bytes[j] == op {
            count += 1;
            j += 1;
        }

        // Convert to CIGAR format with I/D swap
        let op_char = match op {
            b'M' => '=', // Match (WFA2 uses M for exact match)
            b'X' => 'X', // Mismatch
            b'I' => 'D', // WFA2 'I' means standard 'D'
            b'D' => 'I', // WFA2 'D' means standard 'I'
            _ => '?',
        };

        cigar_str.push_str(&format!("{}{}", count, op_char));
        i = j;
    }

    cigar_str
}

/// Create and configure WFA2 aligner with proper settings
pub fn create_aligner(
    match_score: i32,
    mismatch_penalty: i32,
    gap1_open: i32,
    gap1_extend: i32,
    gap2_open: Option<i32>,
    gap2_extend: Option<i32>,
) -> AffineWavefronts {
    let mut wf = if let (Some(gap2_open), Some(gap2_extend)) = (gap2_open, gap2_extend) {
        AffineWavefronts::with_penalties_affine2p_and_memory_mode(
            match_score,
            mismatch_penalty,
            gap1_open,
            gap1_extend,
            gap2_open,
            gap2_extend,
            MemoryMode::Ultralow,
        )
    } else {
        AffineWavefronts::with_penalties_and_memory_mode(
            match_score,
            mismatch_penalty,
            gap1_open,
            gap1_extend,
            MemoryMode::Ultralow,
        )
    };

    // Critical configuration settings for reliable global alignment
    wf.set_alignment_scope(AlignmentScope::Alignment);
    wf.set_alignment_span(AlignmentSpan::End2End);
    wf.set_heuristic(&HeuristicStrategy::None);

    wf
}