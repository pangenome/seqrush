/// This file provides examples and tests similar to the file edaligner.cpp in  apps directory
/// of edlib. It loads simple fasta files of one sequence in directory test_data of edlib.
///

extern crate edlib_rs;
use edlib_rs::edlibrs::*;

use clap::{App, Arg};
use std::path::Path;
use std::process;
use log::*;

use ::cpu_time::ProcessTime;
use std::time::Duration;
/// example
/// edaligner --dirdata  "/Soft/edlib/test_data/Enterobacteria_Phage_1"
///           --tf "Enterobacteria_phage_1.fasta"
///           --qf "mutated_60_perc.fasta"
fn main() {
    // initialize logger from env variable RUST_LOG
    env_logger::Builder::from_default_env().init();

    let dirdata : String;
    let qfile : String;
    let tfile : String;

    let matches = App::new("edaligner")
        .arg(Arg::with_name("dirdata")
            .long("dirdata")
            .required(true)
            .takes_value(true)
            .help("expection directory of data files"))
        .arg(Arg::with_name("qfile")
            .long("qf")
            .takes_value(true)
            .help("expection query file for seq"))
        .arg(Arg::with_name("tfile")
            .long("tf")
            .takes_value(true)
            .help("expection target file for seq"))
        .get_matches();

    // get data directory
    if matches.is_present("dirdata") {
        println!("dirdata");
        dirdata = matches.value_of("dirdata").ok_or("bad value").unwrap().parse::<String>().unwrap();
        println!("got dirdata , {}", dirdata);
    }
    else {
        println!("dirdata is mandatory");
        process::exit(1);
    }

    // get query file
    if matches.is_present("qfile") {
        qfile = matches.value_of("qfile").ok_or("bad value").unwrap().parse::<String>().unwrap();
        println!("got qfile , {}", qfile);
    }
    else {
        println!("query file is mandatory");
        process::exit(1);
    }

    // get target file
    if matches.is_present("tfile") {
        tfile = matches.value_of("tfile").ok_or("bad value").unwrap().parse::<String>().unwrap();
        println!("got target file , {}", tfile);
    }
    else {
        println!("target file is mandatory");
        process::exit(1);
    }

    let qfname = Path::new(&dirdata).join(qfile);
    info!(" got query file : {:?} ", qfname);
    let tfname = Path::new(&dirdata).join(tfile);
    info!(" got target file : {:?} ",tfname);
    // use logger
    let qseq : Vec<u8>;
    // get sequences
    let mut reader = needletail::parse_fastx_file(&qfname).expect("expecting valid query filename");
    if let Some(record) = reader.next() {
        let qrec = record.expect("invalid record");
        let n_bases = qrec.num_bases();
        println!(" query nb bases : {}", n_bases);
        qseq = qrec.seq().into_owned();
    } else {
        std::process::exit(1);
    } // end for query seq
    //
    let tseq: Vec<u8>;
    // get sequences
    let mut reader = needletail::parse_fastx_file(&tfname).expect("expecting valid target filename");
    if let Some(record) = reader.next() {
        let trec = record.expect("invalid record");
        let n_bases = trec.num_bases();
        println!(" target nb bases : {}", n_bases);
        tseq = trec.seq().into_owned();
    } else {
        std::process::exit(1);
    } // end for query seq
    let mut config = EdlibAlignConfigRs::default();
    let mut mod_str;
    //
    mod_str = "EDLIB_MODE_NW";
    config.mode = EdlibAlignModeRs::EDLIB_MODE_NW;
    let start = ProcessTime::try_now().unwrap();
    let align_res = edlibAlignRs(&qseq, &tseq, &config);
    assert_eq!(align_res.status, EDLIB_STATUS_OK);
    let cpu_time: Duration = start.try_elapsed().unwrap();
    println!("\n mode : {}, cpu time (ms) {} distance : {} ", mod_str , cpu_time.as_millis(), align_res.editDistance);
    //
    mod_str = "EDLIB_MODE_SHW";
    config.mode = EdlibAlignModeRs::EDLIB_MODE_SHW;
    let start = ProcessTime::try_now().unwrap();
    let align_res = edlibAlignRs(&qseq, &tseq, &config);
    assert_eq!(align_res.status, EDLIB_STATUS_OK);
    let cpu_time: Duration = start.try_elapsed().unwrap();
    println!("\n mode : {}, cpu time (ms) {} distance : {} ", mod_str , cpu_time.as_millis(), align_res.editDistance);
    //
    mod_str = "EDLIB_MODE_HW";
    config.mode = EdlibAlignModeRs::EDLIB_MODE_HW;
    let start = ProcessTime::try_now().unwrap();
    let align_res = edlibAlignRs(&qseq, &tseq, &config);
    assert_eq!(align_res.status, EDLIB_STATUS_OK);
    let cpu_time: Duration = start.try_elapsed().unwrap();
    println!("\n mode : {}, cpu time (ms) {} distance : {} ", mod_str , cpu_time.as_millis(), align_res.editDistance);
}