#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
//include!(concat!(env!("OUT_DIR"), "/bindings.rs"));



use env_logger::{Builder};

#[macro_use]
extern crate lazy_static;

pub mod edlibrs;


lazy_static! {
    #[allow(dead_code)]
    pub static ref LOG: u64 = {
        let res = init_log();
        res
    };
}
// install a logger facility
fn init_log() -> u64 {
    let mut builder = Builder::from_default_env();
    builder.init();
    println!("\n ************** initializing logger from env *****************\n");    
    return 1;
}

