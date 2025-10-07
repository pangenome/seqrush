# fasteval
Fast and safe evaluation of algebraic expressions

`fasteval` is a library for parsing, compiling, and evaluating algebraic expressions.
It can be used directly as a calculator language (much like `python`), and it is
an excellent foundation for building higher-level-languages.

Documentation:

* [API Reference (docs.rs)](https://docs.rs/fasteval/)


## Usage

Add this to your Cargo.toml:

    [dependencies]
    fasteval = "0.2.4"


You should **always** build with `RUSTFLAGS="--emit=asm"` because [it greatly improves LLVM's compile-time optimizations](http://likebike.com/posts/How_To_Write_Fast_Rust_Code.html#emit-asm).

If you are using a 'nightly' Rust compiler, you can build with `--features nightly` to enable optimizations that aren't yet available in 'stable' Rust.

You can build with `--no-default-features` to disable alphabetical keywords like `and`, `or`, `NaN`, `inf`.  (These words might be important to your applications.)

You can build with `--features unsafe-vars` to enable [Unsafe Variables](https://docs.rs/fasteval/#unsafe-variables).


## Features
* No dependencies.
* Safe execution of untrusted expressions.
* Works with stable Rust.
* Supports interpretation (i.e. parse & eval) as well as compiled execution (i.e. parse, compile, eval).
* Supports Variables and Custom Functions.
* `fasteval` is a good base for building higher-level languages.
* Supports many built-in functions and constants.
* Supports all the standard algebraic unary and binary operators (+ - * / ^ %),
  as well as comparisons (< <= == != >= >) and logical operators (&& ||) with
  short-circuit support.
* Easy integration into many different types of applications, including scoped evaluation.
* Very fast performance.


## Easy Example

Here is one simple example.  See the [API Reference](https://docs.rs/fasteval/#examples) for many more!

The `ez_eval()` function performs the entire allocation-parse-eval process
for you.  It is slightly inefficient because it always allocates a
fresh [`Slab`](https://docs.rs/fasteval/latest/fasteval/slab/index.html), but it is very simple to use:

```rust
fn main() -> Result<(), fasteval::Error> {
    // This example doesn't use any variables, so just use an EmptyNamespace:
    let mut ns = fasteval::EmptyNamespace;

    let val = fasteval::ez_eval(
        "1+2*3/4^5%6 + log(100K) + log(e(),100) + [3*(3-3)/3] + (2<3) && 1.23",    &mut ns)?;
    //    |            |      |    |   |          |               |   |
    //    |            |      |    |   |          |               |   boolean logic with short-circuit support
    //    |            |      |    |   |          |               comparisons
    //    |            |      |    |   |          square-brackets act like parenthesis
    //    |            |      |    |   built-in constants: e(), pi()
    //    |            |      |    'log' can take an optional first 'base' argument, defaults to 10
    //    |            |      numeric literal with suffix: p, n, µ, m, K, M, G, T
    //    |            many built-in functions: print, int, ceil, floor, abs, sign, log, round, min, max, sin, asin, ...
    //    standard binary operators

    assert_eq!(val, 1.23);

    Ok(())
}
```


## REPL Demo
```text
github.com/likebike/fasteval$ rlwrap cargo run --release --example repl
    Finished release [optimized] target(s) in 0.01s
     Running `target/release/examples/repl`
>>> print("Hello fasteval", 1, 2, 3)
Hello fasteval 1 2 3
3
>>> _ + 1
4
>>> _ + 1
5
>>> _ * 2
10
>>> _ ^ 0.5
3.1622776601683795
>>> let a = 1
1
>>> let b = a + 1
2
>>> let c = a + b * 3
7
>>> a + b + c
10
>>> push
Entered scope[1]
>>> let b = b + 10
12
>>> a + b + c
20
>>> pop
Exited scope[1]
>>> a + b + c
10
>>> 1+2*3/4^5%6 + log(100K) + log(e(),100) + [3*(3-3)/3] + (2<3) && 1.23
1.23
>>> 1+2*3/4^5%6 + print("log(100K) =",log(100K)) + log(e(),100) + [3*(3-3)/3] + (2<3) && 1.23
log(100K) = 5
1.23
```

## Safety

`fasteval` is designed to evaluate untrusted expressions safely.  By
default, an expression can only perform math operations; there is no way
for it to access other types of operations (like network or filesystem or
external commands).  Additionally, we guard against malicious expressions:

* Expressions that are too large (greater than 4KB).
* Expressions that are too-deeply nested (greater than 32 levels).
* Expressions with too many values (greater than 64).
* Expressions with too many sub-expressions (greater than 64).

All limits can be customized at parse time.  If any limits are exceeded,
[`parse()`](https://docs.rs/fasteval/latest/fasteval/parser/struct.Parser.html#method.parse) will return an
[Error](https://docs.rs/fasteval/latest/fasteval/error/enum.Error.html).

Note that it *is* possible for you (the developer) to define custom functions
which might perform dangerous operations.  It is your responsibility to make
sure that all custom functionality is safe.


## Performance Benchmarks

Here is a short summary of the performance benchmarks.  For a more complete report and anlysis, see the [official documentation](https://docs.rs/fasteval/#performance-benchmarks).

### Charts
Note that the following charts use logarithmic scales.  Therefore, tiny
visual differences actually represent very significant performance
differences.

**Performance of evaluation of a compiled expression:**  
![Compiled Eval Performance](https://raw.githubusercontent.com/likebike/fasteval/master/benches/results/20191225/fasteval-compiled.png)

**Performance of one-time interpretation (parse and eval):**  
![Interpretation Performance](https://raw.githubusercontent.com/likebike/fasteval/master/benches/results/20191225/fasteval-interp.png)

**Performance of compiled Unsafe Variables, compared to the tinyexpr C library (the
only other library in our test set that supports this mode):**  
![Unsafe Compiled Eval Performance](https://raw.githubusercontent.com/likebike/fasteval/master/benches/results/20191225/fasteval-compiled-unsafe.png)

**Performance of interpreted Unsafe Variables, compared to the tinyexpr C library (the
only other library in our test set that supports this mode):**  
![Unsafe Interpretation Performance](https://raw.githubusercontent.com/likebike/fasteval/master/benches/results/20191225/fasteval-interp-unsafe.png)

### Summary

The impressive thing about these results is that `fasteval` consistently
achieves the fastest times across every benchmark and in every mode of
operation (interpreted, compiled, and unsafe).  It's easy to create a
design to claim the #1 spot in any one of these metrics by sacrificing
performance in another, but it is difficult to create a design that can be
#1 across-the-board.

Because of the broad and robust performance advantages, `fasteval` is very
likely to be an excellent choice for your dynamic evaluation needs.


## License
`fasteval` is distributed under the terms the MIT license.

See [LICENSE](https://github.com/likebike/fasteval/blob/master/LICENSE) for details.

