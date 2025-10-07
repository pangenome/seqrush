#[macro_use]
extern crate bencher;
extern crate btoi;

use bencher::{black_box, Bencher};

fn bench_btou(b: &mut Bencher) {
    b.iter(|| assert_eq!(Ok(123_456_789u32), btoi::btou(black_box(b"123456789"))));
}

fn bench_btou_saturating(b: &mut Bencher) {
    b.iter(|| {
        assert_eq!(
            Ok(123_456_789u32),
            btoi::btou_saturating(black_box(b"123456789"))
        )
    });
}

fn bench_small_btou(b: &mut Bencher) {
    b.iter(|| assert_eq!(Ok(42u8), btoi::btou(black_box(b"42"))));
}

fn bench_small_btou_saturating(b: &mut Bencher) {
    b.iter(|| assert_eq!(Ok(255u8), btoi::btou_saturating(black_box(b"256"))));
}

fn bench_btoi(b: &mut Bencher) {
    b.iter(|| assert_eq!(Ok(-123_456_789i32), btoi::btoi(black_box(b"-123456789"))));
}

fn bench_from_str(b: &mut Bencher) {
    fn btou_from_str(s: &[u8]) -> Option<u32> {
        ::std::str::from_utf8(s).ok().and_then(|s| s.parse().ok())
    }

    b.iter(|| assert_eq!(Some(123_456_789), btou_from_str(black_box(b"123456789"))));
}

fn bench_from_str_unchecked(b: &mut Bencher) {
    unsafe fn btou_from_str_unchecked(s: &[u8]) -> Option<u32> {
        ::std::str::from_utf8_unchecked(s).parse().ok()
    }

    b.iter(|| {
        assert_eq!(Some(123_456_789), unsafe {
            btou_from_str_unchecked(black_box(b"123456789"))
        })
    });
}

fn bench_small_from_str(b: &mut Bencher) {
    fn btou_from_str(s: &[u8]) -> Option<u8> {
        ::std::str::from_utf8(s).ok().and_then(|s| s.parse().ok())
    }

    b.iter(|| assert_eq!(Some(42u8), btou_from_str(black_box(b"42"))));
}

fn bench_small_from_str_unchecked(b: &mut Bencher) {
    unsafe fn btou_from_str_unchecked(s: &[u8]) -> Option<u8> {
        ::std::str::from_utf8_unchecked(s).parse().ok()
    }

    b.iter(|| {
        assert_eq!(Some(42u8), unsafe {
            btou_from_str_unchecked(black_box(b"42"))
        })
    });
}

benchmark_group!(
    benches,
    bench_btou,
    bench_btou_saturating,
    bench_btoi,
    bench_from_str,
    bench_from_str_unchecked,
    bench_small_btou,
    bench_small_btou_saturating,
    bench_small_from_str,
    bench_small_from_str_unchecked
);

benchmark_main!(benches);
