/// Abstraction over different alignment backends (allwave, sweepga, etc.)
use std::error::Error;

/// Represents a single sequence
#[derive(Clone, Debug)]
pub struct AlignmentSequence {
    pub id: String,
    pub seq: Vec<u8>,
}

/// Represents a PAF-like alignment record
#[derive(Clone, Debug)]
pub struct AlignmentRecord {
    pub query_name: String,
    pub query_len: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: char,
    pub target_name: String,
    pub target_len: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub cigar: String,
}

/// Trait for alignment backends
pub trait Aligner {
    /// Align sequences and return alignment records
    fn align_sequences(
        &self,
        sequences: &[AlignmentSequence],
    ) -> Result<Vec<AlignmentRecord>, Box<dyn Error>>;
}

/// Factory function to create the appropriate aligner based on compile-time features
pub fn create_aligner(
    threads: usize,
    verbose: bool,
    #[allow(unused_variables)] frequency: Option<usize>,
) -> Result<Box<dyn Aligner>, Box<dyn Error>> {
    #[cfg(feature = "use-allwave")]
    {
        Ok(Box::new(allwave_impl::AllwaveAligner::new(
            threads, verbose,
        )?))
    }

    #[cfg(not(feature = "use-allwave"))]
    {
        Err("No aligner feature enabled. Enable 'use-allwave'".into())
    }
}

// Import implementations based on features
#[cfg(feature = "use-allwave")]
pub mod allwave_impl;

// Sweepga support temporarily disabled for CI
// #[cfg(feature = "use-sweepga")]
// pub mod sweepga_impl;
