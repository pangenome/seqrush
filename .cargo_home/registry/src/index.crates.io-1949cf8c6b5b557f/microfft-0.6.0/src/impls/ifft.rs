use crate::impls::cfft::*;
use crate::Complex32;

pub(crate) trait IFft {
    type CFft: CFft;

    const N: usize = Self::CFft::N;

    #[inline]
    fn transform(x: &mut [Complex32]) {
        debug_assert_eq!(x.len(), Self::N);

        Self::reorder(x);
        Self::CFft::transform(x);
        Self::normalize(x);
    }

    #[inline]
    fn reorder(x: &mut [Complex32]) {
        debug_assert_eq!(x.len(), Self::N);

        let m = Self::N / 2;
        for i in 1..m {
            x.swap(i, Self::N - i);
        }
    }

    #[inline]
    fn normalize(x: &mut [Complex32]) {
        for c in x {
            *c /= Self::N as f32;
        }
    }
}

macro_rules! ifft_impls {
    ( $( ($IFftN:ident, $CFftN:ident), )* ) => {
        $(
            pub(crate) struct $IFftN;

            impl IFft for $IFftN {
                type CFft = $CFftN;
            }
        )*
    };
}

ifft_impls! {
    (IFftN2, CFftN2),
    (IFftN4, CFftN4),
    (IFftN8, CFftN8),
    (IFftN16, CFftN16),
    (IFftN32, CFftN32),
    (IFftN64, CFftN64),
    (IFftN128, CFftN128),
    (IFftN256, CFftN256),
    (IFftN512, CFftN512),
    (IFftN1024, CFftN1024),
    (IFftN2048, CFftN2048),
    (IFftN4096, CFftN4096),
    (IFftN8192, CFftN8192),
    (IFftN16384, CFftN16384),
    (IFftN32768, CFftN32768),
}
