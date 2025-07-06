#[cfg(feature = "cli")]
use seqrush::Args;

#[cfg(feature = "cli")]
/// Parse command line arguments into `Args`.
pub fn parse() -> Args {
    let mut sequences = None;
    let mut output = None;
    let mut threads = 1_usize;
    let mut min_match_length = 15_usize;

    let mut iter = std::env::args().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "-s" | "--sequences" => sequences = iter.next(),
            "-o" | "--output" => output = iter.next(),
            "-t" | "--threads" => {
                if let Some(t) = iter.next() {
                    if let Ok(v) = t.parse() {
                        threads = v;
                    }
                }
            }
            "-k" | "--min-match-length" => {
                if let Some(m) = iter.next() {
                    if let Ok(v) = m.parse() {
                        min_match_length = v;
                    }
                }
            }
            _ => {}
        }
    }
    let sequences = sequences.expect("input FASTA required");
    let output = output.expect("output file required");
    Args {
        sequences,
        output,
        threads,
        min_match_length,
    }
}

