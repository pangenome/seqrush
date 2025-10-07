# lcs-diff-rs
Compute differences between two slices using LCS algorithm.

[![CircleCI](https://circleci.com/gh/bokuweb/lcs-diff-rs.svg?style=svg)](https://circleci.com/gh/bokuweb/lcs-diff-rs)

## Example

``` rust
extern crate lcs_diff;

use self::lcs_diff::*;

fn main() {
    let old = vec!["foo", "bar", "baz"];
    let new = vec!["foo", "baz", "hoge"];

    for diff in lcs_diff::diff(&old, &new) {
        match diff {
            DiffResult::Added(a) => println!("+{} new index = {}", a.data, a.new_index.unwrap()),
            DiffResult::Common(c) => {
                println!(" {} old index = {}, new index = {}",
                         c.data,
                         c.old_index.unwrap(),
                         c.new_index.unwrap())
            }
            DiffResult::Removed(r) => println!("-{} old index = {}", r.data, r.old_index.unwrap()),
        }
    }
}
```

You can also run example as shown below.

``` bash
rustup run nightly cargo run --example example
```


## LICENSE

The MIT License (MIT)

Copyright (c) 2018 @bokuweb

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
