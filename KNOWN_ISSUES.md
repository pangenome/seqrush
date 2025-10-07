# Known Issues

## Deprecation Warning: buf_redux

There is a deprecation warning about the `buf_redux` crate:
```
warning: the following packages contain code that will be rejected by a future version of Rust: buf_redux v0.8.4
```

### Dependency Chain
The warning comes from a transitive dependency:
- `seqrush` → `allwave` → `edlib_rs` → `needletail v0.4.1` → `buf_redux v0.8.4`

### Issue Details
- `buf_redux` is unmaintained since 2019
- `needletail` has a newer version (0.6.3) that likely doesn't use `buf_redux`
- `edlib_rs` needs to update its dependency on `needletail` from v0.4 to v0.6

### Potential Solutions
1. **Upstream Fix**: File an issue with [edlib-rs](https://github.com/jean-pierreBoth/edlib-rs) to update needletail dependency
2. **Fork**: Fork edlib-rs and update the dependency ourselves
3. **Alternative**: Find an alternative to edlib-rs that doesn't have this issue

### Current Status
This warning doesn't affect functionality but may become an error in future Rust versions. We should monitor this and plan to address it before it becomes blocking.