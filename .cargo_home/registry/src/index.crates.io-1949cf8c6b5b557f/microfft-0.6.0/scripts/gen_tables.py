#! /usr/bin/env python3

"""
Script for generating the pre-computed tables used by microfft:
  - radix-2 FFT twiddle tables
  - bit reversal tables

Used to create the file `src/tables.rs`.
"""

import argparse
import math


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("N", type=int, help="Max FFT size")
    return parser.parse_args()


def emit_header():
    print("//! This file was generated with gen_tables.py.")
    print()
    print(f"#![allow(clippy::excessive_precision)]")
    print(f"#![allow(clippy::unreadable_literal)]")
    print()


def emit_sine(max_n):
    print("cfg_if::cfg_if! {")

    n = max_n
    while n > 2:
        kw = "if" if n == max_n else "else if"
        print(f'{kw} #[cfg(feature = "size-{n}")] {{')
        emit_sine_table(n)
        print("}", end=" ")
        n //= 2

    print()
    print("}")
    print()


def emit_sine_table(n):
    print("pub(crate) const SINE: &[f32] = &[")
    for k in range(1, n // 4):
        sine = math.sin(-2 * math.pi * k / n)
        print(f"    {sine},")
    print("];")
    print()


def emit_bitrev(max_n):
    print('#[cfg(feature = "bitrev-tables")]')
    print("pub(crate) const BITREV: &[&[u16]] = &[")

    n = 1
    while n <= max_n:
        emit_bitrev_table(n)
        n *= 2

    print("];")
    print()


def emit_bitrev_table(n):
    print("    &[")
    nbits = int(math.log2(n))
    for i in range(n):
        rev = reverse_bits(i, nbits)
        entry = rev if rev > i else i
        print(f"        {entry},")
    print("    ],")


def reverse_bits(num, nbits):
    reverse = 0
    for i in range(nbits):
        if num & (1 << i):
            reverse |= 1 << (nbits - 1 - i)
    return reverse


def main():
    args = parse_args()

    emit_header()
    emit_sine(args.N)
    emit_bitrev(args.N)


if __name__ == "__main__":
    main()
