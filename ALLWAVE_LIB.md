# Allwave Library Design

## Overview

This document outlines the design for refactoring allwave from a standalone CLI tool into a library that can be used by seqrush and other tools for all-to-all sequence alignment.

## Core Design Principles

1. **Lazy Iteration**: Generate alignments on-demand rather than computing all upfront
2. **Parallelizable**: Support efficient parallel processing with Rayon
3. **Interruptible**: Allow consumers to stop alignment generation at any point
4. **Minimal Overhead**: Avoid unnecessary allocations or string formatting
5. **Resource Reuse**: Reuse WFA2 aligners across multiple alignments

## Library Structure

```
allwave/
├── Cargo.toml
├── src/
│   ├── lib.rs          # Public API
│   ├── main.rs         # CLI tool (uses library)
│   ├── alignment.rs    # Core alignment logic
│   ├── iterator.rs     # Alignment iterator implementation
│   └── types.rs        # Shared types
```

## Public API

### Core Types

```rust
// In types.rs
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    // Indices into the original sequence array
    pub query_idx: usize,
    pub target_idx: usize,
    
    // Alignment coordinates (0-based)
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    
    // Orientation
    pub is_reverse: bool,  // true if target was reverse-complemented
    
    // Alignment details
    pub cigar_bytes: Vec<u8>,  // Raw CIGAR from WFA2
    pub score: i32,            // WFA2 score (lower is better)
    pub num_matches: usize,    // Number of matching bases
    pub alignment_length: usize, // Total alignment length
}

#[derive(Debug, Clone)]
pub struct AlignmentParams {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub gap2_open: Option<i32>,
    pub gap2_extend: Option<i32>,
    pub max_divergence: Option<f64>,
}
```

### Iterator Interface

```rust
// In iterator.rs
pub struct AllPairIterator<'a> {
    sequences: &'a [Sequence],
    params: AlignmentParams,
    orientation_params: AlignmentParams,
    exclude_self: bool,
    sparsification: SparsificationStrategy,
    pair_iter: Box<dyn Iterator<Item = (usize, usize)> + Send + 'a>,
}

#[derive(Debug, Clone)]
pub enum SparsificationStrategy {
    /// Align all pairs
    None,
    /// Random sampling with given probability
    Random(f64),
    /// Automatic sparsification based on sequence count
    Auto,
}

impl<'a> AllPairIterator<'a> {
    /// Create iterator for all-vs-all alignments
    pub fn new(sequences: &'a [Sequence], params: AlignmentParams) -> Self {
        Self::with_options(sequences, params, true, SparsificationStrategy::None)
    }
    
    /// Create iterator with custom options
    pub fn with_options(
        sequences: &'a [Sequence], 
        params: AlignmentParams,
        exclude_self: bool,
        sparsification: SparsificationStrategy,
    ) -> Self {
        let n = sequences.len();
        
        // Generate all pairs
        let mut pairs: Vec<(usize, usize)> = if exclude_self {
            (0..n).flat_map(|i| (0..n).filter(move |&j| i != j).map(move |j| (i, j))).collect()
        } else {
            (0..n).flat_map(|i| (0..n).map(move |j| (i, j))).collect()
        };
        
        // Apply sparsification
        match sparsification {
            SparsificationStrategy::None => {},
            SparsificationStrategy::Random(keep_fraction) => {
                pairs = apply_random_sparsification(pairs, keep_fraction);
            },
            SparsificationStrategy::Auto => {
                let keep_fraction = compute_auto_sparsification(n);
                pairs = apply_random_sparsification(pairs, keep_fraction);
            },
        }
        
        Self {
            sequences,
            params,
            orientation_params: AlignmentParams {
                match_score: 0,
                mismatch_penalty: 1,
                gap_open: 1,
                gap_extend: 1,
                gap2_open: None,
                gap2_extend: None,
                max_divergence: None,
            },
            exclude_self,
            sparsification,
            pair_iter: Box::new(pairs.into_iter()),
        }
    }
    
    /// Set custom orientation detection parameters
    pub fn with_orientation_params(mut self, params: AlignmentParams) -> Self {
        self.orientation_params = params;
        self
    }
    
    /// Set sparsification strategy
    pub fn with_sparsification(mut self, strategy: SparsificationStrategy) -> Self {
        // Need to regenerate pairs with new sparsification
        *self = Self::with_options(self.sequences, self.params.clone(), self.exclude_self, strategy);
        self
    }
}

impl<'a> Iterator for AllPairIterator<'a> {
    type Item = AlignmentResult;
    
    fn next(&mut self) -> Option<Self::Item> {
        // Get next pair
        let (i, j) = self.pair_iter.next()?;
        
        // Perform alignment (implementation in alignment.rs)
        let result = align_pair(
            &self.sequences[i],
            &self.sequences[j],
            i,
            j,
            &self.params,
            &self.orientation_params,
        );
        
        Some(result)
    }
}

// Support for Rayon parallel iteration
impl<'a> ParallelIterator for AllPairIterator<'a> {
    type Item = AlignmentResult;
    
    fn drive_unindexed<C>(self, consumer: C) -> C::Result
    where
        C: UnindexedConsumer<Self::Item>,
    {
        // Create pairs upfront for parallel processing
        let pairs: Vec<_> = self.pair_iter.collect();
        
        pairs
            .into_par_iter()
            .map(|(i, j)| {
                align_pair(
                    &self.sequences[i],
                    &self.sequences[j],
                    i,
                    j,
                    &self.params,
                    &self.orientation_params,
                )
            })
            .drive_unindexed(consumer)
    }
}
```

### Convenience Functions

```rust
// In lib.rs

/// Convert CIGAR bytes to standard string representation
pub fn cigar_bytes_to_string(cigar: &[u8]) -> String {
    // Implementation here
}

/// Format alignment result as PAF record
pub fn alignment_to_paf(
    result: &AlignmentResult,
    sequences: &[Sequence],
) -> String {
    // Implementation here
}

/// Apply random sparsification to pairs using deterministic hashing
fn apply_random_sparsification(
    mut pairs: Vec<(usize, usize)>, 
    keep_fraction: f64
) -> Vec<(usize, usize)> {
    use std::hash::{Hash, Hasher};
    use std::collections::hash_map::DefaultHasher;
    
    pairs.retain(|(i, j)| {
        let mut hasher = DefaultHasher::new();
        i.hash(&mut hasher);
        j.hash(&mut hasher);
        let hash = hasher.finish();
        
        // Convert hash to a value between 0 and 1
        let normalized = (hash as f64) / (u64::MAX as f64);
        normalized < keep_fraction
    });
    
    pairs
}

/// Compute automatic sparsification factor based on sequence count
fn compute_auto_sparsification(n: usize) -> f64 {
    match n {
        0..=50 => 1.0,          // No sparsification for small datasets
        51..=200 => 0.8,        // Keep 80% of pairs
        201..=500 => 0.5,       // Keep 50% of pairs
        501..=1000 => 0.3,      // Keep 30% of pairs
        1001..=5000 => 0.1,     // Keep 10% of pairs
        _ => 0.05,              // Keep 5% for very large datasets
    }
}
```

## Usage Examples

### Sequential Processing
```rust
use allwave::{AllPairIterator, AlignmentParams};

let sequences = load_sequences("input.fa");
let params = AlignmentParams::default();

let aligner = AllPairIterator::new(&sequences, params);

for alignment in aligner {
    println!("Aligned {} to {}", alignment.query_idx, alignment.target_idx);
    // Process alignment...
    
    // Can stop early if needed
    if some_condition {
        break;
    }
}
```

### Parallel Processing with Rayon
```rust
use rayon::prelude::*;

let aligner = AllPairIterator::new(&sequences, params);

aligner
    .par_bridge()
    .for_each(|alignment| {
        // Process each alignment in parallel
        process_alignment(alignment);
    });
```

### Using Sparsification
```rust
// Random sparsification - keep 30% of pairs
let aligner = AllPairIterator::new(&sequences, params)
    .with_sparsification(SparsificationStrategy::Random(0.3));

// Automatic sparsification based on dataset size
let aligner = AllPairIterator::new(&sequences, params)
    .with_sparsification(SparsificationStrategy::Auto);

for alignment in aligner {
    // Process sparse subset of alignments
}
```

### Collecting Results
```rust
// Get all alignments meeting some criteria
let good_alignments: Vec<_> = AllPairIterator::new(&sequences, params)
    .filter(|a| a.score < threshold)
    .collect();
```

### Writing PAF Output (in main.rs)
```rust
use std::io::{BufWriter, Write};
use rayon::prelude::*;

let aligner = AllPairIterator::new(&sequences, params);
let paf_records: Vec<String> = aligner
    .par_bridge()
    .map(|alignment| allwave::alignment_to_paf(&alignment, &sequences))
    .collect();

let mut writer = BufWriter::new(output_file);
for record in paf_records {
    writeln!(writer, "{}", record)?;
}
```

## Implementation Notes

1. **WFA2 Aligner Reuse**: The `align_pair` function should use thread-local storage to reuse WFA2 aligners:
   ```rust
   thread_local! {
       static WFA_ALIGNER: RefCell<Option<AffineWavefronts>> = RefCell::new(None);
       static WFA_ORIENTATION: RefCell<Option<AffineWavefronts>> = RefCell::new(None);
   }
   ```

2. **Memory Efficiency**: Avoid cloning sequences unnecessarily. The iterator holds references to the original sequence array.

3. **Error Handling**: The iterator should handle alignment failures gracefully by skipping failed pairs rather than panicking.

4. **Backward Compatibility**: The existing CLI interface should continue to work using the library implementation.

## Migration Path for seqrush

1. Add allwave as a dependency:
   ```toml
   [dependencies]
   allwave = { path = "../allwave" }
   ```

2. Replace alignment logic with allwave iterator:
   ```rust
   use allwave::{AllPairIterator, AlignmentParams, SparsificationStrategy};
   
   let params = AlignmentParams {
       match_score: scores.match_score,
       mismatch_penalty: scores.mismatch_penalty,
       // ...
   };
   
   // Parse sparsification from command line
   let sparsification = match args.sparsification.as_str() {
       "1.0" => SparsificationStrategy::None,
       "auto" => SparsificationStrategy::Auto,
       s => {
           let factor = s.parse::<f64>().expect("Invalid sparsification factor");
           SparsificationStrategy::Random(factor)
       }
   };
   
   let aligner = AllPairIterator::new(&sequences, params)
       .with_orientation_params(orientation_params)
       .with_sparsification(sparsification);
   
   // For PAF output
   if let Some(paf_writer) = paf_writer {
       let paf_aligner = AllPairIterator::new(&sequences, params)
           .with_options(&sequences, params, true, SparsificationStrategy::None); // No sparsification for PAF
   }
   
   aligner
       .par_bridge()
       .for_each(|alignment| {
           // Process alignment for graph building
           let cigar = allwave::cigar_bytes_to_string(&alignment.cigar_bytes);
           self.process_alignment(
               &cigar,
               &sequences[alignment.query_idx],
               &sequences[alignment.target_idx],
               args.min_match_length,
               alignment.is_reverse,
               args.verbose
           );
       });
   ```

## Benefits

1. **Single Source of Truth**: One implementation of all-to-all alignment logic
2. **Flexible**: Can be used for PAF output, graph building, or other purposes
3. **Efficient**: Lazy evaluation and parallel processing support
4. **Interruptible**: Natural support for early termination
5. **Testable**: Library can be unit tested independently of CLI