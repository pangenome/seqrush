//! A simple bloom filter implementation for `u64` hash values only.

#[cfg(feature = "no_std")]
use std::prelude::v1::*;
use std::{
    mem,
    ops::{Add, AddAssign},
};

/// Number of bits for a `usize`.
const USIZE_BITS: usize = mem::size_of::<usize>() * 8;

/// Number of `usize` values required for 256 bits.
const SIZE: usize = 256 / USIZE_BITS;

/// _(internals)_ A simple bloom filter implementation for `u64` hash values only - i.e. all 64 bits are assumed
/// to be relatively random.
/// Exported under the `internals` feature only.
///
/// For this reason, the implementation is simplistic - it just looks at the least significant byte
/// of the `u64` hash value and sets the corresponding bit in a 256-long bit vector.
///
/// The rationale of this type is to avoid pulling in another dependent crate.
#[derive(Debug, Clone, Eq, PartialEq, Hash, Default)]
pub struct BloomFilterU64([usize; SIZE]);

impl BloomFilterU64 {
    /// Get the bit position of a `u64` hash value.
    #[inline(always)]
    #[must_use]
    const fn calc_hash(value: u64) -> (usize, usize) {
        let hash = (value & 0x00ff) as usize;
        (hash / USIZE_BITS, 0x01 << (hash % USIZE_BITS))
    }
    /// Create a new [`BloomFilterU64`].
    #[inline(always)]
    #[must_use]
    pub const fn new() -> Self {
        Self([0; SIZE])
    }
    /// Is this [`BloomFilterU64`] empty?
    #[inline(always)]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0 == [0; SIZE]
    }
    /// Clear this [`BloomFilterU64`].
    #[inline(always)]
    pub fn clear(&mut self) -> &mut Self {
        self.0 = [0; SIZE];
        self
    }
    /// Mark a `u64` hash into this [`BloomFilterU64`].
    #[inline]
    pub fn mark(&mut self, hash: u64) -> &mut Self {
        let (offset, mask) = Self::calc_hash(hash);
        self.0[offset] |= mask;
        self
    }
    /// Is a `u64` hash definitely absent from this [`BloomFilterU64`]?
    #[inline]
    #[must_use]
    pub const fn is_absent(&self, hash: u64) -> bool {
        let (offset, mask) = Self::calc_hash(hash);
        (self.0[offset] & mask) == 0
    }
    /// If a `u64` hash is absent from this [`BloomFilterU64`], return `true` and then mark it.
    /// Otherwise return `false`.
    #[inline]
    #[must_use]
    pub fn is_absent_and_set(&mut self, hash: u64) -> bool {
        let (offset, mask) = Self::calc_hash(hash);
        let result = (self.0[offset] & mask) == 0;
        self.0[offset] |= mask;
        result
    }
}

impl Add for &BloomFilterU64 {
    type Output = BloomFilterU64;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        let mut buf = [0; SIZE];

        self.0
            .iter()
            .zip(rhs.0.iter())
            .map(|(&a, &b)| a | b)
            .zip(buf.iter_mut())
            .for_each(|(v, x)| *x = v);

        BloomFilterU64(buf)
    }
}

impl Add<BloomFilterU64> for &BloomFilterU64 {
    type Output = BloomFilterU64;

    #[inline(always)]
    fn add(self, rhs: BloomFilterU64) -> Self::Output {
        self + &rhs
    }
}

impl AddAssign<Self> for BloomFilterU64 {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl AddAssign<&Self> for BloomFilterU64 {
    #[inline]
    fn add_assign(&mut self, rhs: &Self) {
        self.0
            .iter_mut()
            .zip(rhs.0.iter())
            .for_each(|(x, &v)| *x |= v);
    }
}
