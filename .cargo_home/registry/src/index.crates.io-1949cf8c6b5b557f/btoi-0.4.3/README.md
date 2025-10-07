btoi
====

Parse integers from ASCII byte slices.

[![crates.io](https://img.shields.io/crates/v/btoi.svg)](https://crates.io/crates/btoi)
[![docs.rs](https://docs.rs/btoi/badge.svg)](https://docs.rs/btoi)

Introduction
------------

Provides functions similar to [`from_str_radix`](https://doc.rust-lang.org/std/primitive.u32.html#method.from_str_radix),
but is faster when parsing directly from byte slices instead of strings.

Supports `#![no_std]`.

```rust
use btoi::btoi;

assert_eq!(Ok(42), btoi(b"42"));
assert_eq!(Ok(-1000), btoi(b"-1000"));
```

Documentation
-------------

[Read the documentation](https://docs.rs/btoi)

Changelog
---------

* 0.4.3
  - Use `#[track_caller]`.
* 0.4.2
  - No longer `!#[deny(warnings)]`, which is is a forwards compability hazard
    in libraries.
  - Explicit `!#[forbid(unsafe_code)]`.
* 0.4.1
  - `-` was parsed as zero, but should have errored. Thanks @wayslog.
* 0.4.0
  - Change type of radix to `u32` (from `u8`) to mirror the standard library.
  - No need to `#[inline]` generic functions.
* 0.3.0
  - New default feature `std`. Disable for `#![no_std]` support.
  - Mark functions as `#[inline]`.
* 0.2.0
  - No longer reexport num-traits.
* 0.1.3
  - Update to num-traits 0.2 (semver compatible).
* 0.1.2
  - Fix documentation warnings.
  - Update dependencies.
* 0.1.1
  - Documentation fixes.
* 0.1.0
  - Initial release.

License
-------

btoi is dual licensed under the [Apache 2.0](http://www.apache.org/licenses/LICENSE-2.0)
and [MIT](http://opensource.org/licenses/MIT) license, at your option.
