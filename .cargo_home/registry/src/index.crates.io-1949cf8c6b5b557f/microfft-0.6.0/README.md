# microfft

microfft is a library for computing fast fourier transforms that targets
embedded systems. It provides an in-place implementation of the Radix-2 FFT
algorithm. All computations are performed directly on the input buffer and
require no additional allocations. This makes microfft suitable for `no_std`
environments, like microcontrollers.

Speed is achieved mainly by maintaining a pre-computed sine table that is used
to look up the necessary twiddle factors. By replacing arithmetic operations
with simple memory lookups, we reduce the number of CPU cycles spent.
Unfortunately, the pre-computed table also claims a considerable amount of
memory, which might be a deal-breaker for some embedded projects (see
[Memory Usage](#memory-usage)).

microfft also implements a specialized algorithm for FFTs on real (instead
of complex) values. Naively one would calculate a real FFT simply by converting
the input to complex values (leaving the imaginary part empty) and running a
CFFT. microfft's RFFT algorithm instead packs pairs of real values into
a single complex one each, then computes a CFFT of half the original input
size, followed by some recombination magic. This has the effect of roughly
halving the number of CPU cycles required, as can be seen in the
[benchmark results][bench/README.md].

## Example

The following example demonstrates computing a 16-point RFFT on a set of
samples generated from a sine signal:

```rust
use std::convert::TryInto;
use std::f32::consts::PI;

// generate 16 samples of a sine wave at frequency 3
let sample_count = 16;
let signal_freq = 3.;
let sample_interval = 1. / sample_count as f32;
let mut samples: Vec<_> = (0..sample_count)
    .map(|i| (2. * PI * signal_freq * sample_interval * i as f32).sin())
    .collect();

// compute the RFFT of the samples
let mut samples: [_; 16] = samples.try_into().unwrap();
let spectrum = microfft::real::rfft_16(&mut samples);
// since the real-valued coefficient at the Nyquist frequency is packed into the
// imaginary part of the DC bin, it must be cleared before computing the amplitudes
spectrum[0].im = 0.0;

// the spectrum has a spike at index `signal_freq`
let amplitudes: Vec<_> = spectrum.iter().map(|c| c.norm() as u32).collect();
assert_eq!(&amplitudes, &[0, 0, 0, 8, 0, 0, 0, 0]);
```

## Requirements

Requires Rust version **1.61.0** or newer.

## Sine Tables

microfft keeps a single sine table to calculate the twiddle factors for all
FFT sizes. This removes some memory overhead compared to keeping a separate
table for each FFT size, as there would be duplication between those tables.

The default sine table supports FFTs up to size 4096. If you only want to
compute FFTs of smaller sizes, it is recommended to select the appropriate
`size-*` feature, to not waste memory. For example, if your maximum FFT size is
1024, add this to your `Cargo.toml`:

```toml
[dependencies.microfft]
default-features = false
features = ["size-1024"]
```

This tells microfft to not provide functions for computing FFTs of sizes larger
than 1024 and to keep only the 1024-point sine table.

If you want to compute FFTs with more than 4096 points, you also need to enable
the respective feature. In this case, disabling the default features is not
required as microfft always determines the sine table size based on the largest
size requested. So this works as expected:

```toml
[dependencies.microfft]
features = ["size-8192"]
```

## Bit-reversal Tables

The optional feature `bitrev-tables` enables the use of pre-computed tables of
bit-reversed indices required for the reordering of input values performed at
the start of each FFT. If this feature is disabled (the default), the
bit-reversals are computed at runtime instead.

Note that enabling bitrev tables significantly increases the memory usage of
microfft. While it can speed up FFT computation on some systems, there are also
architectures that provide dedicated bit-reversal instructions (like `RBIT` on
ARMv7). On such architectures, switching on bitrev tables is usually
detrimental to performance.

## `std` Usage

microfft provides a `std` feature meant to make the library more useful for
applications that can make use of the Rust standard library. Currently the only
thing it does is transitively enabling the `std` feature of the `num-complex`
crate, thereby making more methods available on the `Complex32` values returned
by the FFT functions.

As embedded applications usually run on targets that don't have a Rust standard
library, the `std` feature is disabled by default. You can enable it in your
`Cargo.toml`:

```toml
[dependencies.microfft]
features = ["std"]
```

## Limitations

microfft has a few limitations, mostly due to its focus on speed, that might
make it unsuitable for some embedded projects. You should know about these
if you consider using this library:

### Memory Usage

The use of pre-computed sine and bitrev tables means that microfft has
considerable requirements on read-only memory. If your chip doesn't have much
flash to begin with, this can be an issue.

The amount of memory required for tables depends on the the configuration of
the [`size-*`](#sine-tables) and [`bitrev-tables`](#bit-reversal-tables)
features:

| `size-*`     | without `bitrev-tables` | with `bitrev-tables` |
| ------------ | ----------------------: | -------------------: |
| `size-4`     |                       0 |                    8 |
| `size-8`     |                       4 |                   20 |
| `size-16`    |                      12 |                   44 |
| `size-32`    |                      28 |                   92 |
| `size-64`    |                      60 |                  188 |
| `size-128`   |                     124 |                  380 |
| `size-256`   |                     252 |                  764 |
| `size-512`   |                     508 |                1,532 |
| `size-1024`  |                   1,020 |                3,068 |
| `size-2048`  |                   2,044 |                6,140 |
| `size-4096`  |                   4,092 |               12,284 |
| `size-8192`  |                   8,188 |               24,572 |
| `size-16384` |                  16,380 |               49,148 |
| `size-32768` |                  32,764 |               98,300 |

In addition, the code size also increases with FFT size.

### Supported FFT Sizes

microfft only supports FFT point-sizes that are powers of two, a limitation of
the Radix-2 algorithm. Additionally, the maximum supported size is currently
32768, although this limit can be increased in the future as necessary.

### `f64` Support

This library currently only supports single-precision floating-point inputs.
Similarly to the FFT size limit, this is a restriction that might be lifted
in the future, should the need arise.

## License

This project is licensed under the MIT license ([LICENSE](LICENSE) or
http://opensource.org/licenses/MIT).

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in microfft by you, shall be licensed as above, without any
additional terms or conditions.
