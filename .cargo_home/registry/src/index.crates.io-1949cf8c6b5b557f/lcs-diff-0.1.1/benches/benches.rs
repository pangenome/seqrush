#![feature(plugin, test)]
#![plugin(speculate)]

extern crate lcs_diff;
extern crate test;

speculate! {
    describe "diff" {
        bench "1000 equal items" |b| {
            let slice = [0u8; 1000];
            b.iter(|| ::lcs_diff::diff(&slice, &slice));
        }

        bench "1000 non-equal items" |b| {
            let (left, right) = ([0u8; 1000], [1u8; 1000]);
            b.iter(|| ::lcs_diff::diff(&left, &right));
        }
    }
}