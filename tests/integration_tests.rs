use seqrush::{Args, run_seqrush, load_sequences};
use std::fs::{self, File};
use std::io::Write;
use std::env::temp_dir;

fn temp_file(name: &str) -> std::path::PathBuf {
    let mut path = temp_dir();
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_micros();
    path.push(format!("{}_{}", name, now));
    path
}

#[test]
fn load_sequences_parses_fasta() {
    let path = temp_file("seqs");
    let mut f = File::create(&path).unwrap();
    writeln!(f, ">a\nACGT\n>b\nTTTT").unwrap();
    f.sync_all().unwrap();
    let seqs = load_sequences(path.to_str().unwrap()).unwrap();
    assert_eq!(seqs.len(), 2);
    assert_eq!(seqs[0].id, "a");
}

#[test]
fn run_seqrush_writes_output() {
    let in_path = temp_file("in");
    let mut f = File::create(&in_path).unwrap();
    writeln!(f, ">x\nAAAA\n>y\nGGGG").unwrap();
    f.sync_all().unwrap();
    let out_path = temp_file("out");
    let args = Args {
        sequences: in_path.to_str().unwrap().to_string(),
        output: out_path.to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
    };
    run_seqrush(args).unwrap();
    let content = fs::read_to_string(out_path).unwrap();
    assert!(content.starts_with("H\tVN:Z:1.0"));
}

#[cfg(feature = "cli")]
#[test]
fn cli_parses_flags() {
    use std::process::Command;
    let in_path = temp_file("cli_in");
    let mut f = File::create(&in_path).unwrap();
    writeln!(f, ">z\nAAAA").unwrap();
    f.sync_all().unwrap();
    let out_path = temp_file("cli_out");
    let status = Command::new(env!("CARGO_BIN_EXE_seqrush"))
        .args([
            "-s",
            in_path.to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
            "-t",
            "2",
            "-k",
            "5",
        ])
        .status()
        .unwrap();
    assert!(status.success());
}
