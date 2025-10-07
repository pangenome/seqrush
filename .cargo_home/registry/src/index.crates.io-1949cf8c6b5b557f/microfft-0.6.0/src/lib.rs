//! Library for computing fast fourier transforms on embedded systems.
//!
//! microfft provides an in-place implementation of the Radix-2 FFT algorithm.
//! All computations are performed directly on the input buffer and require no
//! additional allocations. This makes microfft suitable for `no_std`
//! environments.
//!
//! This crate provides three FFT implementations:
//!  * [`complex`]: FFT on [`Complex32`] input values (CFFT).
//!  * [`real`]: FFT on real (`f32`) input values (RFFT). An `N`-point RFFT
//!    internally computes an `N/2`-point CFFT, making it roughly twice as fast
//!    as the complex variant.
//!  * [`inverse`]: Inverse FFT (IFFT), implemented in terms of a CFFT.
//!
//! # Example
//!
//! ```
//! use std::convert::TryInto;
//! use std::f32::consts::PI;
//!
//! // generate 16 samples of a sine wave at frequency 3
//! let sample_count = 16;
//! let signal_freq = 3.;
//! let sample_interval = 1. / sample_count as f32;
//! let mut samples: Vec<_> = (0..sample_count)
//!     .map(|i| (2. * PI * signal_freq * sample_interval * i as f32).sin())
//!     .collect();
//!
//! // compute the RFFT of the samples
//! let mut samples: [_; 16] = samples.try_into().unwrap();
//! let spectrum = microfft::real::rfft_16(&mut samples);
//! // since the real-valued coefficient at the Nyquist frequency is packed into the
//! // imaginary part of the DC bin, it must be cleared before computing the amplitudes
//! spectrum[0].im = 0.0;
//!
//! // the spectrum has a spike at index `signal_freq`
//! let amplitudes: Vec<_> = spectrum.iter().map(|c| c.norm() as u32).collect();
//! assert_eq!(&amplitudes, &[0, 0, 0, 8, 0, 0, 0, 0]);
//! ```
//!
//! [`complex`]: complex/index.html
//! [`inverse`]: inverse/index.html
//! [`real`]: real/index.html
//! [`Complex32`]: type.Complex32.html

#![no_std]
#![deny(missing_docs)]
#![warn(rust_2018_idioms)]

pub mod complex;
pub mod inverse;
pub mod real;

pub use num_complex::Complex32;

mod impls {
    pub(crate) mod cfft;
    pub(crate) mod ifft;
    pub(crate) mod rfft;
}
mod tables;

use static_assertions::assert_cfg;

assert_cfg!(
    any(
        feature = "size-4",
        feature = "size-8",
        feature = "size-16",
        feature = "size-32",
        feature = "size-64",
        feature = "size-128",
        feature = "size-256",
        feature = "size-512",
        feature = "size-1024",
        feature = "size-2048",
        feature = "size-4096",
        feature = "size-8192",
        feature = "size-16384",
        feature = "size-32768",
    ),
    "At least one of the `size-*` features of this crate must be set."
);
