use seqrush::{load_sequences, run_seqrush, Args};
use std::fs::{self, File};
use std::io::Write;
use tempfile::NamedTempFile;

fn temp_file() -> NamedTempFile {
    NamedTempFile::new().unwrap()
}

#[test]
fn load_sequences_parses_fasta() {
    let mut file = temp_file();
    writeln!(file, ">a\nACGT\n>b\nTTTT").unwrap();
    file.as_file_mut().sync_all().unwrap();
    let seqs = load_sequences(file.path().to_str().unwrap()).unwrap();
    assert_eq!(seqs.len(), 2);
    assert_eq!(seqs[0].id, "a");
}

#[test]
fn load_sequences_multiline_sequence() {
    let mut file = temp_file();
    writeln!(file, ">id\nACG\nTGA").unwrap();
    file.as_file_mut().sync_all().unwrap();
    let seqs = load_sequences(file.path().to_str().unwrap()).unwrap();
    assert_eq!(seqs.len(), 1);
    assert_eq!(seqs[0].data, b"ACGTGA".to_vec());
}

#[test]
fn load_sequences_windows_newlines() {
    let mut file = temp_file();
    write!(file, ">id\r\nACGT\r\n").unwrap();
    file.as_file_mut().sync_all().unwrap();
    let seqs = load_sequences(file.path().to_str().unwrap()).unwrap();
    assert_eq!(seqs.len(), 1);
    assert_eq!(seqs[0].data, b"ACGT".to_vec());
}

#[test]
fn run_seqrush_writes_output() {
    let mut in_file = temp_file();
    writeln!(in_file, ">x\nAAAA\n>y\nGGGG").unwrap();
    in_file.as_file_mut().sync_all().unwrap();
    let out_file = temp_file();
    let args = Args {
        sequences: in_file.path().to_str().unwrap().to_string(),
        output: out_file.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
    };
    run_seqrush(args).unwrap();
    let content = fs::read_to_string(out_file.path()).unwrap();
    assert!(content.starts_with("H\tVN:Z:1.0"));
    assert!(content.contains("P\tp1\tx+,y+\t*"));
    let links: Vec<_> = content.lines().filter(|l| l.starts_with('L')).collect();
    assert_eq!(links.len(), 1);
    assert_eq!(links[0], "L\tx\t+\ty\t+\t0M");
}

#[test]
fn run_seqrush_multi_sequence_links() {
    let mut in_file = temp_file();
    writeln!(in_file, ">a\nAAAA\n>b\nCCCC\n>c\nGGGG").unwrap();
    in_file.as_file_mut().sync_all().unwrap();
    let out_file = temp_file();
    let args = Args {
        sequences: in_file.path().to_str().unwrap().to_string(),
        output: out_file.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
    };
    run_seqrush(args).unwrap();
    let content = fs::read_to_string(out_file.path()).unwrap();
    assert!(content.contains("P\tp1\ta+,b+,c+\t*"));
    let links: Vec<_> = content.lines().filter(|l| l.starts_with('L')).collect();
    assert_eq!(links.len(), 2);
    assert!(links.contains(&"L\ta\t+\tb\t+\t0M"));
    assert!(links.contains(&"L\tb\t+\tc\t+\t0M"));
}
#[test]
fn load_sequences_empty_input() {
    let file = temp_file();
    let seqs = load_sequences(file.path().to_str().unwrap()).unwrap();
    assert!(seqs.is_empty());
}

#[test]
fn load_sequences_missing_file() {
    let tmp = temp_file();
    let path = tmp.path().to_path_buf();
    drop(tmp);
    let result = load_sequences(path.to_str().unwrap());
    assert!(result.is_err());
}

#[test]
fn load_sequences_large_input() {
    let mut file = temp_file();
    for i in 0..1000 {
        writeln!(file, ">{}\nACGTACGTACGTACGTACGTACGTACGTACGT", i).unwrap();
    }
    file.as_file_mut().sync_all().unwrap();
    let seqs = load_sequences(file.path().to_str().unwrap()).unwrap();
    assert_eq!(seqs.len(), 1000);
}

#[test]
fn run_seqrush_missing_input() {
    let out_file = temp_file();
    let missing_path = {
        let tmp = temp_file();
        let path = tmp.path().to_path_buf();
        drop(tmp);
        path
    };
    let args = Args {
        sequences: missing_path.to_str().unwrap().to_string(),
        output: out_file.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
    };
    let result = run_seqrush(args);
    assert!(result.is_err());
    let out_path = out_file.path();
    if out_path.exists() {
        fs::remove_file(out_path).unwrap();
    }
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
    let output = Command::new(exe).args(["-s", "somefile"]).output().unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("required arguments"));
}

#[cfg(feature = "cli")]
#[test]
fn cli_parses_flags() {
    let mut in_file = temp_file();
    writeln!(in_file, ">z\nAAAA").unwrap();
    in_file.as_file_mut().sync_all().unwrap();
    let out_file = temp_file();
    let status = Command::new(env!("CARGO_BIN_EXE_seqrush"))
        .args([
            "-s",
            in_file.path().to_str().unwrap(),
            "-o",
            out_file.path().to_str().unwrap(),
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
    let mut fasta_file = temp_file();
    writeln!(fasta_file, ">id\nACGT").unwrap();
    fasta_file.as_file_mut().sync_all().unwrap();

    let gfa_file = temp_file();
    let args = Args {
        sequences: fasta_file.path().to_str().unwrap().to_string(),
        output: gfa_file.path().to_str().unwrap().to_string(),
        threads: 1,
        min_match_length: 1,
    };
    run_seqrush(args).unwrap();

    let content = fs::read_to_string(gfa_file.path()).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines[0], "H\tVN:Z:1.0");
    let s_lines: Vec<&&str> = lines.iter().filter(|l| l.starts_with("S\t")).collect();
    assert_eq!(s_lines.len(), 1);
    assert_eq!(s_lines[0], &"S\tid\tACGT");
    let p_lines: Vec<&&str> = lines.iter().filter(|l| l.starts_with("P\t")).collect();
    assert_eq!(p_lines.len(), 1);
    assert_eq!(p_lines[0], &"P\tp1\tid+\t*");
    assert!(lines.iter().all(|l| !l.starts_with("L\t")));
}
