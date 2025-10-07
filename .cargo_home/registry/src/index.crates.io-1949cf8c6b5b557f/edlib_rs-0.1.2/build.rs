extern crate bindgen;


use std::env;
use std::path::PathBuf;
use std::fs::File;
use std::io::prelude::*;

// We have path to edlib, generate wrapper.h with  path to correct include 
fn write_wrapper(wrapper_path : &PathBuf) {
    let mut inclpath = PathBuf::new();
    inclpath.push("include");
    inclpath.push("edlib.h");
    let mut incl = String::from("#include\"");
    incl.push_str(& inclpath.to_str().unwrap());
    incl.push_str("\"");
    // now we write file wrapper.h
    let mut f = File::create(wrapper_path).unwrap();
    f.write_all(incl.as_bytes()).expect("cannot write wrapper.h");
}



fn main() {
    //  try to build edlib in build.rs with edlib cloned in edlib-c (from which we removed )
    let dst = cmake::Config::new("edlib-c")
                .cflag("-D CMAKE_BUILD_TYPE=Release")
                .build();
    println!("cargo:rustc-link-search=native={}", dst.display());
    //
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    let wrapper_path = out_path.join("wrapper.h");
    //
    write_wrapper(&wrapper_path);
    // adapted fix from sam217pa. libstdc++ name varies...
    let lib_std : &str;
    if cfg!(target_os = "macos") {
        lib_std  = "cargo:rustc-link-lib=c++";
    }
    else {
        lib_std = "cargo:rustc-link-lib=stdc++";
    }
    
    let mut libdir = String::from("cargo:rustc-link-search=");
    libdir.push_str(out_path.to_str().unwrap());
    libdir.push_str("/lib");
    println!("{}",libdir);
    println!("cargo:rustc-link-lib=edlib");
    println!("{}", lib_std);

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header(wrapper_path.to_str().unwrap())
        //Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
         // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    //
    println!("cargo:rerun-if-changed=binding.rs");
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
    //
}
