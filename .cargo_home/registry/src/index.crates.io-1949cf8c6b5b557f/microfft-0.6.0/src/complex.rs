//! FFT on complex inputs (CFFT)

use crate::impls::cfft::*;
use crate::Complex32;

macro_rules! cfft_impls {
    ( $( $N:expr => ($cfft_N:ident, $CFftN:ident $(, $feature:expr)?), )* ) => {
        $(
            #[doc = concat!("Perform an in-place ", stringify!($N), "-point CFFT.")]
            #[doc = ""]
            #[doc = "# Example"]
            #[doc = ""]
            #[doc = "```"]
            #[doc = concat!("use microfft::{Complex32, complex::", stringify!($cfft_N), "};")]
            #[doc = ""]
            #[doc = concat!("let mut input = [Complex32::default(); ", stringify!($N), "];")]
            #[doc = concat!("let result = ", stringify!($cfft_N), "(&mut input);")]
            #[doc = "```"]
            $( #[cfg(feature = $feature)] )?
            #[inline]
            #[must_use]
            pub fn $cfft_N(input: &mut [Complex32; $N]) -> &mut [Complex32; $N] {
                $CFftN::transform(input);
                input
            }
        )*
    };
}

cfft_impls! {
    2 => (cfft_2, CFftN2),
    4 => (cfft_4, CFftN4, "size-4"),
    8 => (cfft_8, CFftN8, "size-8"),
    16 => (cfft_16, CFftN16, "size-16"),
    32 => (cfft_32, CFftN32, "size-32"),
    64 => (cfft_64, CFftN64, "size-64"),
    128 => (cfft_128, CFftN128, "size-128"),
    256 => (cfft_256, CFftN256, "size-256"),
    512 => (cfft_512, CFftN512, "size-512"),
    1024 => (cfft_1024, CFftN1024, "size-1024"),
    2048 => (cfft_2048, CFftN2048, "size-2048"),
    4096 => (cfft_4096, CFftN4096, "size-4096"),
    8192 => (cfft_8192, CFftN8192, "size-8192"),
    16384 => (cfft_16384, CFftN16384, "size-16384"),
    32768 => (cfft_32768, CFftN32768, "size-32768"),
}
