# Issue Draft for edlib-rs

## Title: Update needletail dependency to avoid deprecated buf_redux

## Description

The current version of edlib-rs uses needletail v0.4.1, which depends on buf_redux v0.8.4. The buf_redux crate has been unmaintained since 2019 and now generates deprecation warnings:

```
warning: the following packages contain code that will be rejected by a future version of Rust: buf_redux v0.8.4
```

This warning indicates that buf_redux contains code that will not compile with future versions of Rust.

## Proposed Solution

Update the needletail dependency from v0.4 to v0.6.3 (the latest version). The newer version of needletail likely doesn't depend on buf_redux anymore.

In `Cargo.toml`, change:
```toml
needletail = "0.4"
```

to:
```toml
needletail = "0.6"
```

## Impact

This change would:
- Remove the deprecation warning for downstream users
- Ensure compatibility with future Rust versions
- Potentially bring performance improvements from newer needletail versions

## Additional Context

We're using edlib-rs through the allwave library in our project seqrush, and this deprecation warning is the only remaining warning in our build.

Would you be open to accepting a PR for this update? I'd be happy to test the changes and ensure compatibility.

---
*Note: This is a draft. You can submit this issue at https://github.com/jean-pierreBoth/edlib-rs/issues/new*