Rhai Tools
==========

Tools for working with Rhai scripts.

| Tool                                                                             | Required feature(s) | Description                                           |
| -------------------------------------------------------------------------------- | :-----------------: | ----------------------------------------------------- |
| [`rhai-run`](https://github.com/rhaiscript/rhai/blob/main/src/bin/rhai-run.rs)   |                     | runs each filename passed to it as a Rhai script      |
| [`rhai-repl`](https://github.com/rhaiscript/rhai/blob/main/src/bin/rhai-repl.rs) |     `rustyline`     | a simple REPL that interactively evaluates statements |
| [`rhai-dbg`](https://github.com/rhaiscript/rhai/blob/main/src/bin/rhai-dbg.rs)   |     `debugging`     | the _Rhai Debugger_                                   |

For convenience, a feature named `bin-features` is available which is a combination of the following:

* `decimal` &ndash; support for decimal numbers
* `metadata` &ndash; access functions metadata
* `serde` &ndash; export functions metadata to JSON
* `debugging` &ndash; required by `rhai-dbg`
* `rustyline` &ndash; required by `rhai-repl`


How to Run
----------

```sh
cargo run --features bin-features --bin sample_app_to_run
```


How to Install
--------------

To install these all tools (with full features), use the following command:

```sh
cargo install --path . --bins  --features bin-features
```

or specifically:

```sh
cargo install --path . --bin sample_app_to_run  --features bin-features
```
