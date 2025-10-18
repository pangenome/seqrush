/// Allwave aligner implementation
use super::{Aligner, AlignmentRecord, AlignmentSequence};
use allwave::{AlignmentParams, AllPairIterator, SparsificationStrategy};
use std::error::Error;

pub struct AllwaveAligner {
    threads: usize,
    verbose: bool,
    params: AlignmentParams,
}

impl AllwaveAligner {
    pub fn new(threads: usize, verbose: bool) -> Result<Self, Box<dyn Error>> {
        // Default alignment parameters (matching seqrush defaults)
        let params = AlignmentParams {
            match_score: 2,
            mismatch_penalty: 4,
            gap_open: 4,
            gap_extend: 2,
            gap2_open: Some(24),
            gap2_extend: Some(1),
            max_divergence: None,
        };

        Ok(AllwaveAligner {
            threads,
            verbose,
            params,
        })
    }

    pub fn with_params(
        threads: usize,
        verbose: bool,
        params: AlignmentParams,
    ) -> Result<Self, Box<dyn Error>> {
        Ok(AllwaveAligner {
            threads,
            verbose,
            params,
        })
    }
}

impl AllwaveAligner {
    /// Get prioritized pair list using AllWave's sparsification
    pub fn get_prioritized_pairs(
        sequences: &[AlignmentSequence],
        sparsification: SparsificationStrategy,
        verbose: bool,
    ) -> Result<Vec<(usize, usize)>, Box<dyn Error>> {
        if verbose {
            eprintln!("[allwave] Computing prioritized pair list for {} sequences", sequences.len());
        }

        // Convert to allwave format
        let allwave_sequences: Vec<allwave::Sequence> = sequences
            .iter()
            .map(|s| allwave::Sequence {
                id: s.id.clone(),
                seq: s.seq.clone(),
            })
            .collect();

        // Create iterator with sparsification to get prioritized pairs
        let params = AlignmentParams {
            match_score: 2,
            mismatch_penalty: 4,
            gap_open: 4,
            gap_extend: 2,
            gap2_open: Some(24),
            gap2_extend: Some(1),
            max_divergence: None,
        };

        let aligner = AllPairIterator::with_options(
            &allwave_sequences,
            params,
            false, // exclude self-alignments
            false, // don't use mash orientation
            sparsification,
        );

        // Extract the pair list without performing alignments
        let pairs = aligner.get_pairs();

        if verbose {
            eprintln!("[allwave] Generated {} prioritized pairs", pairs.len());
        }

        Ok(pairs)
    }
}

impl Aligner for AllwaveAligner {
    fn align_sequences(
        &self,
        sequences: &[AlignmentSequence],
    ) -> Result<Vec<AlignmentRecord>, Box<dyn Error>> {
        if self.verbose {
            eprintln!("[allwave] Aligning {} sequences", sequences.len());
        }

        // Convert to allwave format
        let allwave_sequences: Vec<allwave::Sequence> = sequences
            .iter()
            .map(|s| allwave::Sequence {
                id: s.id.clone(),
                seq: s.seq.clone(),
            })
            .collect();

        // Create iterator over all pairs
        let aligner = AllPairIterator::with_options(
            &allwave_sequences,
            self.params.clone(),
            false, // exclude self-alignments
            false, // don't use mash orientation
            SparsificationStrategy::None,
        );

        // Collect alignments
        let mut records = Vec::new();
        for alignment in aligner {
            let cigar = allwave::cigar_bytes_to_string(&alignment.cigar_bytes);

            let record = AlignmentRecord {
                query_name: allwave_sequences[alignment.query_idx].id.clone(),
                query_len: allwave_sequences[alignment.query_idx].seq.len(),
                query_start: alignment.query_start,
                query_end: alignment.query_end,
                strand: if alignment.is_reverse { '-' } else { '+' },
                target_name: allwave_sequences[alignment.target_idx].id.clone(),
                target_len: allwave_sequences[alignment.target_idx].seq.len(),
                target_start: alignment.target_start,
                target_end: alignment.target_end,
                cigar,
            };

            records.push(record);
        }

        if self.verbose {
            eprintln!("[allwave] Alignment complete: {} records", records.len());
        }

        Ok(records)
    }
}
