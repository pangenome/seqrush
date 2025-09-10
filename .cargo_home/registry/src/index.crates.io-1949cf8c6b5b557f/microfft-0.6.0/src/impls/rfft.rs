use core::slice;

use static_assertions::{assert_eq_align, assert_eq_size};

use crate::impls::cfft::*;
use crate::{tables, Complex32};

pub(crate) trait RFft {
    type CFft: CFft;

    const N: usize = Self::CFft::N * 2;

    #[inline]
    fn transform(x: &mut [f32]) -> &mut [Complex32] {
        debug_assert_eq!(x.len(), Self::N);

        let x = Self::pack_complex(x);

        Self::CFft::transform(x);
        Self::recombine(x);
        x
    }

    #[inline]
    fn pack_complex(x: &mut [f32]) -> &mut [Complex32] {
        assert_eq_size!(Complex32, [f32; 2]);
        assert_eq_align!(Complex32, f32);
        assert_eq!(x.len(), Self::N);

        let len = Self::N / 2;
        let data = x.as_mut_ptr().cast::<Complex32>();
        unsafe { slice::from_raw_parts_mut(data, len) }
    }

    #[inline]
    fn recombine(x: &mut [Complex32]) {
        let m = Self::CFft::N;
        debug_assert_eq!(x.len(), m);

        let table_len = tables::SINE.len();
        let table_stride = (table_len + 1) * 4 / Self::N;

        // The real part of the first element is the DC value.
        // Additionally, the real-valued coefficient at the Nyquist frequency
        // is stored in the imaginary part.
        let x0 = x[0];
        x[0] = Complex32::new(x0.re + x0.im, x0.re - x0.im);

        let u = m / 2;
        for k in 1..u {
            let s = k * table_stride;
            let twiddle_re = tables::SINE[table_len - s] * -1.;
            let twiddle_im = tables::SINE[s - 1];

            let (x_k, x_nk) = (x[k], x[m - k]);
            let sum = (x_k + x_nk) / 2.;
            let diff = (x_k - x_nk) / 2.;

            x[k] = Complex32::new(
                sum.re + twiddle_re * sum.im + twiddle_im * diff.re,
                diff.im + twiddle_im * sum.im - twiddle_re * diff.re,
            );
            x[m - k] = Complex32::new(
                sum.re - twiddle_re * sum.im - twiddle_im * diff.re,
                -diff.im + twiddle_im * sum.im - twiddle_re * diff.re,
            );
        }

        let xu = x[u];
        x[u] = Complex32::new(xu.re, -xu.im);
    }
}

pub(crate) struct RFftN2;

impl RFft for RFftN2 {
    type CFft = CFftN1;

    #[inline]
    fn recombine(x: &mut [Complex32]) {
        debug_assert_eq!(x.len(), 1);

        // The real part of the first element is the DC value.
        // Additionally, the real-valued coefficient at the Nyquist frequency
        // is stored in the imaginary part.
        let x0 = x[0];
        x[0] = Complex32::new(x0.re + x0.im, x0.re - x0.im);
    }
}

macro_rules! rfft_impls {
    ( $( ($RFftN:ident, $CFftN:ident), )* ) => {
        $(
            #[allow(dead_code)]
            pub(crate) struct $RFftN;

            impl RFft for $RFftN {
                type CFft = $CFftN;
            }
        )*
    };
}

rfft_impls! {
    (RFftN4, CFftN2),
    (RFftN8, CFftN4),
    (RFftN16, CFftN8),
    (RFftN32, CFftN16),
    (RFftN64, CFftN32),
    (RFftN128, CFftN64),
    (RFftN256, CFftN128),
    (RFftN512, CFftN256),
    (RFftN1024, CFftN512),
    (RFftN2048, CFftN1024),
    (RFftN4096, CFftN2048),
    (RFftN8192, CFftN4096),
    (RFftN16384, CFftN8192),
    (RFftN32768, CFftN16384),
}
