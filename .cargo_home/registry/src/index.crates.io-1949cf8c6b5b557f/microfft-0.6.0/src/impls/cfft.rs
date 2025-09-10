use crate::{tables, Complex32};

pub(crate) trait CFft {
    type Half: CFft;

    const N: usize;
    const LOG2_N: usize;

    #[cfg(feature = "bitrev-tables")]
    const BITREV_TABLE: &'static [u16] = tables::BITREV[Self::LOG2_N];

    #[inline]
    fn transform(x: &mut [Complex32]) -> &mut [Complex32] {
        debug_assert_eq!(x.len(), Self::N);

        Self::bit_reverse_reorder(x);
        Self::compute_butterflies(x);
        x
    }

    #[cfg(feature = "bitrev-tables")]
    #[inline]
    fn bit_reverse_reorder(x: &mut [Complex32]) {
        debug_assert_eq!(x.len(), Self::N);

        for i in 0..Self::N {
            let j = Self::BITREV_TABLE[i] as usize;
            x.swap(i, j);
        }
    }

    #[cfg(not(feature = "bitrev-tables"))]
    #[inline]
    fn bit_reverse_reorder(x: &mut [Complex32]) {
        debug_assert_eq!(x.len(), Self::N);

        let shift = core::mem::size_of::<usize>() * 8 - Self::LOG2_N;
        for i in 0..Self::N {
            let rev = i.reverse_bits();
            let j = rev >> shift;
            if j > i {
                x.swap(i, j);
            }
        }
    }

    #[inline]
    fn compute_butterflies(x: &mut [Complex32]) {
        debug_assert_eq!(x.len(), Self::N);

        let m = Self::N / 2;
        let u = m / 2;

        let table_len = tables::SINE.len();
        let table_stride = (table_len + 1) * 4 / Self::N;

        Self::Half::compute_butterflies(&mut x[..m]);
        Self::Half::compute_butterflies(&mut x[m..]);

        // [k = 0] twiddle factor: `1 + 0i`
        let (x_0, x_m) = (x[0], x[m]);
        x[0] = x_0 + x_m;
        x[m] = x_0 - x_m;

        // [k in [1, m/2)] twiddle factor:
        //   - re from SINE table backwards and negative
        //   - im from SINE table directly
        for k in 1..u {
            let s = k * table_stride;
            let re = tables::SINE[table_len - s] * -1.;
            let im = tables::SINE[s - 1];
            let twiddle = Complex32::new(re, im);

            let (x_k, x_km) = (x[k], x[k + m]);
            let y = twiddle * x_km;
            x[k] = x_k + y;
            x[k + m] = x_k - y;
        }

        // [k = m/2] twiddle factor: `0 - 1i`
        let (x_u, x_um) = (x[u], x[u + m]);
        let y = x_um * Complex32::new(0., -1.);
        x[u] = x_u + y;
        x[u + m] = x_u - y;

        // [k in (m/2, m)] twiddle factor:
        //   - re from SINE table directly
        //   - im from SINE table backwards
        for k in (u + 1)..m {
            let s = (k - u) * table_stride;
            let re = tables::SINE[s - 1];
            let im = tables::SINE[table_len - s];
            let twiddle = Complex32::new(re, im);

            let (x_k, x_km) = (x[k], x[k + m]);
            let y = twiddle * x_km;
            x[k] = x_k + y;
            x[k + m] = x_k - y;
        }
    }
}

pub(crate) struct CFftN1;

impl CFft for CFftN1 {
    type Half = Self;

    const N: usize = 1;
    const LOG2_N: usize = 0;

    #[inline]
    fn bit_reverse_reorder(x: &mut [Complex32]) {
        debug_assert_eq!(x.len(), 1);
    }

    #[inline]
    fn compute_butterflies(x: &mut [Complex32]) {
        debug_assert_eq!(x.len(), 1);
    }
}

pub(crate) struct CFftN2;

impl CFft for CFftN2 {
    type Half = CFftN1;

    const N: usize = 2;
    const LOG2_N: usize = 1;

    #[inline]
    fn compute_butterflies(x: &mut [Complex32]) {
        debug_assert_eq!(x.len(), 2);

        let (x_0, x_1) = (x[0], x[1]);
        x[0] = x_0 + x_1;
        x[1] = x_0 - x_1;
    }
}

macro_rules! cfft_impls {
    ( $( $I:expr => ($N:expr, $CFftN:ident, $Half:ident), )* ) => {
        $(
            #[allow(dead_code)]
            pub(crate) struct $CFftN;

            impl CFft for $CFftN {
                type Half = $Half;

                const N: usize = $N;
                const LOG2_N: usize = $I;
            }
        )*
    };
}

cfft_impls! {
     2 => (4, CFftN4, CFftN2),
     3 => (8, CFftN8, CFftN4),
     4 => (16, CFftN16, CFftN8),
     5 => (32, CFftN32, CFftN16),
     6 => (64, CFftN64, CFftN32),
     7 => (128, CFftN128, CFftN64),
     8 => (256, CFftN256, CFftN128),
     9 => (512, CFftN512, CFftN256),
    10 => (1024, CFftN1024, CFftN512),
    11 => (2048, CFftN2048, CFftN1024),
    12 => (4096, CFftN4096, CFftN2048),
    13 => (8192, CFftN8192, CFftN4096),
    14 => (16384, CFftN16384, CFftN8192),
    15 => (32768, CFftN32768, CFftN16384),
}
