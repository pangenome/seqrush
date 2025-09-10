# Changelog

All notable changes to this project will be documented in this file.

## 0.6.0 (2024-04-14)

### Changed

- **Breaking:** The MSRV has been increased to 1.61.0.

### Added

- Support for FFT size 32768.


## 0.5.1 (2023-05-18)

### Added

- Support for inverse FFTs.


## 0.5.0 (2022-06-19)

### Changed

- **Breaking:** The MSRV has been increased to 1.56.0.
- **Breaking:** The former `maxn-*` features of this crate have been renamed
  to `size-*`. They still work basically the same way but are now additive,
  which means it is allowed that more than one of them is enabled at the same
  time. This makes microfft work as expected in the face of Cargo's feature
  unification.

### Added

- A `std` feature that is disabled by default. It transitively enables the
  `std` feature on the `num-complex` dependency, enabling more methods on the
  `Complex32` values returned by the FFT functions.
- Support for FFT sizes 8192 and 16384.


## 0.4.0 (2021-04-03)

### Changed

- **Breaking:** The MSRV has been increased to 1.51.0.
- **Breaking:** All API functions for computing FFTs now take references to
  arrays instead of references to slices. This has the benefit of moving the
  length checking to compile time, which makes those functions panic-proof.
- **Breaking:** When computing RFFTs, the real-valued coefficient at the
  Nyquist frequency is now packed into the imaginary part of the DC bin
  (bin 0). Previously this value was simply dropped from the RFFT output.
- Added `#[must_use]` annotations to all FFT API functions.

### Fixed

- Thanks to Cargo's new resolver it is now possible to build microfft as a
  stand-alone library.


## 0.3.1 (2020-11-09)

### Fixed

- Fixed a bug during the RFFT recombination calculation that caused a wrong
  output value in bin `N / 4`.


## 0.3.0 (2020-03-08)

### Changed

- Store only the largest sine table, instead of one for each FFT size. By
  default this is the table for the 4096-point FFT (the largest supported one),
  but smaller ones can be selected via the new `maxn-*` crate features.


## 0.2.0 (2020-03-08)

### Changed

- Bitrev tables are not used anymore by default, instead the bit-reversed
  indices are computed directly at runtime. This significantly reduces the
  memory usage of microfft. On architectures that provide a dedicated
  bit-reversal instruction (like `RBIT` on ARMv7), speed is also increased.
  The `bitrev-tables` feature can be enabled to still opt into using bitrev
  tables.


## 0.1.2 (2020-03-07)

### Changed

- Store pre-computed sine values instead of full twiddles, reducing the size
  of the twiddle tables to one fourth the prior size.

## 0.1.1 (2020-03-05)

### Added

- Support for FFT sizes 2048 and 4096.


## 0.1.0 (2020-03-04)

### Added

- Support for complex and real FFTs up to size 1024.
