// From vergen documentation
extern crate vergen;
extern crate rustc_version;
use std::process::Command;
use rustc_version::version;
use vergen::{ConstantsFlags, generate_cargo_keys};

fn main() {
    let flags = ConstantsFlags::all();
    generate_cargo_keys(flags).expect("Must be able to generate cargo keys");
    // Based on 
    // https://vallentin.io/2019/06/06/versioning
    // https://stackoverflow.com/questions/43753491/include-git-commit-hash-as-string-into-rust-program
    // https://unix.stackexchange.com/questions/155046/determine-if-git-working-directory-is-clean-from-a-script
    let out = Command::new("git")
        .arg("status")
        .arg("--porcelain")
        .output()
        .expect("Can run git command successfully");
    if out.stdout.is_empty() {
        println!("cargo:rustc-env=WD_IS_CLEAN=true");
    } else {
        println!("cargo:rustc-env=WD_IS_CLEAN=false");
    }
    let cargo_out = Command::new(env!("CARGO"))
        .arg("--version")
        .output()
        .expect("Can run cargo");
    println!("cargo:rustc-env=C_VER={}", std::str::from_utf8(&cargo_out.stdout).unwrap());
    // https://stackoverflow.com/questions/35806568/is-there-a-way-to-detect-the-compiler-version-from-within-a-rust-program
    println!("cargo:rustc-env=V_RUSTC={}", version().unwrap());
    println!("cargo:rustc-env=S_RUSTFLAGS={:?}", option_env!("RUSTFLAGS"));
}
