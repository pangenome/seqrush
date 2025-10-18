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

/// Enum to select which aligner to use at runtime
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlignerBackend {
    AllWave,
    SweepGA,
}

impl std::str::FromStr for AlignerBackend {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "allwave" => Ok(AlignerBackend::AllWave),
            "sweepga" => Ok(AlignerBackend::SweepGA),
            _ => Err(format!("Unknown aligner: {}", s)),
        }
    }
}

impl std::fmt::Display for AlignerBackend {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlignerBackend::AllWave => write!(f, "allwave"),
            AlignerBackend::SweepGA => write!(f, "sweepga"),
        }
    }
}

/// Factory function to create the appropriate aligner based on runtime selection
pub fn create_aligner(
    backend: AlignerBackend,
    threads: usize,
    verbose: bool,
    #[allow(unused_variables)] frequency: Option<usize>,
) -> Result<Box<dyn Aligner>, Box<dyn Error>> {
    match backend {
        AlignerBackend::AllWave => {
            #[cfg(feature = "use-allwave")]
            {
                Ok(Box::new(allwave_impl::AllwaveAligner::new(
                    threads, verbose,
                )?))
            }
            #[cfg(not(feature = "use-allwave"))]
            {
                Err("AllWave aligner not available. Rebuild with --features use-allwave".into())
            }
        }
        AlignerBackend::SweepGA => {
            #[cfg(feature = "use-sweepga")]
            {
                Ok(Box::new(sweepga_impl::SweepgaAligner::new(
                    frequency, threads, verbose,
                )?))
            }
            #[cfg(not(feature = "use-sweepga"))]
            {
                Err("SweepGA aligner not available. Rebuild with --features use-sweepga".into())
            }
        }
    }
}

// Import implementations based on features
#[cfg(feature = "use-allwave")]
pub mod allwave_impl;

#[cfg(feature = "use-sweepga")]
pub mod sweepga_impl;
