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