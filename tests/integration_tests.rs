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
    fs::remove_file(path).unwrap();
}

#[test]
fn load_sequences_multiline_sequence() {
    let path = temp_file("multi");
    let mut f = File::create(&path).unwrap();
    writeln!(f, ">id\nACG\nTGA").unwrap();
    f.sync_all().unwrap();
    let seqs = load_sequences(path.to_str().unwrap()).unwrap();
    assert_eq!(seqs.len(), 1);
    assert_eq!(seqs[0].data, b"ACGTGA".to_vec());
    fs::remove_file(path).unwrap();
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
    let content = fs::read_to_string(&out_path).unwrap();
    assert!(content.starts_with("H\tVN:Z:1.0"));
    assert!(content.contains("P\tp1\tx+,y+\t*"));
    let links: Vec<_> = content
        .lines()
        .filter(|l| l.starts_with('L'))
        .collect();
    assert_eq!(links.len(), 1);
    assert_eq!(links[0], "L\tx\t+\ty\t+\t0M");
    fs::remove_file(&in_path).unwrap();
    fs::remove_file(&out_path).unwrap();
}

#[test]
fn run_seqrush_multi_sequence_links() {
    let in_path = temp_file("multi_in");
    let mut f = File::create(&in_path).unwrap();
    writeln!(f, ">a\nAAAA\n>b\nCCCC\n>c\nGGGG").unwrap();
    f.sync_all().unwrap();
    let out_path = temp_file("multi_out");
    let args = Args {
        sequences: in_path.to_str().unwrap().to_string(),
        output: out_path.to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
    };
    run_seqrush(args).unwrap();
    let content = fs::read_to_string(&out_path).unwrap();
    assert!(content.contains("P\tp1\ta+,b+,c+\t*"));
    let links: Vec<_> = content
        .lines()
        .filter(|l| l.starts_with('L'))
        .collect();
    assert_eq!(links.len(), 2);
    assert!(links.contains(&"L\ta\t+\tb\t+\t0M"));
    assert!(links.contains(&"L\tb\t+\tc\t+\t0M"));
    fs::remove_file(in_path).unwrap();
    fs::remove_file(out_path).unwrap();
}
#[test]
fn load_sequences_empty_input() {
    let path = temp_file("empty");
    File::create(&path).unwrap();
    let seqs = load_sequences(path.to_str().unwrap()).unwrap();
    assert!(seqs.is_empty());
}

#[test]
fn load_sequences_missing_file() {
    let path = temp_file("missing");
    let result = load_sequences(path.to_str().unwrap());
    assert!(result.is_err());
}

#[test]
fn load_sequences_large_input() {
    let path = temp_file("large");
    let mut f = File::create(&path).unwrap();
    for i in 0..1000 {
        writeln!(f, ">{}\nACGTACGTACGTACGTACGTACGTACGTACGT", i).unwrap();
    }
    f.sync_all().unwrap();
    let seqs = load_sequences(path.to_str().unwrap()).unwrap();
    assert_eq!(seqs.len(), 1000);
    fs::remove_file(path).unwrap();
}

#[test]
fn run_seqrush_missing_input() {
    let out_path = temp_file("out_missing");
    let args = Args {
        sequences: temp_file("noexist").to_str().unwrap().to_string(),
        output: out_path.to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
    };
    let result = run_seqrush(args);
    assert!(result.is_err());
}

use std::process::Command;

#[cfg(feature = "cli")]
#[test]
fn cli_no_arguments() {
    let exe = env!("CARGO_BIN_EXE_seqrush");
    let output = Command::new(exe).output().unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("required arguments"));
}

#[cfg(feature = "cli")]
#[test]
fn cli_missing_output() {
    let exe = env!("CARGO_BIN_EXE_seqrush");
    let output = Command::new(exe)
        .args(["-s", "somefile"]) 
        .output()
        .unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("required arguments"));
}

#[cfg(feature = "cli")]
#[test]
fn cli_parses_flags() {
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

#[test]
fn run_seqrush_single_sequence_no_links() {
    let fasta_path = temp_file("single");
    let mut f = File::create(&fasta_path).unwrap();
    writeln!(f, ">id\nACGT").unwrap();
    f.sync_all().unwrap();

    let gfa_path = temp_file("single_out");
    let args = Args {
        sequences: fasta_path.to_str().unwrap().to_string(),
        output: gfa_path.to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
    };
    run_seqrush(args).unwrap();

    let content = fs::read_to_string(&gfa_path).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines[0], "H\tVN:Z:1.0");
    let s_lines: Vec<&&str> = lines.iter().filter(|l| l.starts_with("S\t")).collect();
    assert_eq!(s_lines.len(), 1);
    assert_eq!(s_lines[0], &"S\tid\tACGT");
    let p_lines: Vec<&&str> = lines.iter().filter(|l| l.starts_with("P\t")).collect();
    assert_eq!(p_lines.len(), 1);
    assert_eq!(p_lines[0], &"P\tp1\tid+\t*");
    assert!(lines.iter().all(|l| !l.starts_with("L\t")));

    fs::remove_file(&fasta_path).unwrap();
    fs::remove_file(&gfa_path).unwrap();
}
