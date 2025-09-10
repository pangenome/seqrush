Source Structure
================

Root Sources
------------

| Source file    | Description                                                                     |
| -------------- | ------------------------------------------------------------------------------- |
| `lib.rs`       | Crate root                                                                      |
| `engine.rs`    | The scripting engine, defines the `Engine` type                                 |
| `tokenizer.rs` | Script tokenizer/lexer                                                          |
| `parser.rs`    | Script parser                                                                   |
| `optimizer.rs` | Script optimizer                                                                |
| `defer.rs`     | Utilities for deferred clean-up of resources                                    |
| `reify.rs`     | Utilities for making generic types concrete                                     |
| `tests.rs`     | Unit tests (not integration tests, which are in the main `tests` sub-directory) |


Sub-Directories
---------------

| Sub-directory | Description                                                        |
| ------------- | ------------------------------------------------------------------ |
| `config`      | Configuration                                                      |
| `types`       | Common data types (e.g. `Dynamic`, errors)                         |
| `api`         | Public API for the scripting engine                                |
| `ast`         | AST definition                                                     |
| `module`      | Support for modules                                                |
| `packages`    | Pre-defined packages                                               |
| `func`        | Registering and calling functions (native Rust and script-defined) |
| `eval`        | AST evaluation                                                     |
| `serde`       | Support for [`serde`](https://crates.io/crates/serde) and metadata |
| `bin`         | Pre-built CLI binaries                                             |
