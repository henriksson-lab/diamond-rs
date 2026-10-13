fn main() {
    println!("cargo::rustc-check-cfg=cfg(diamond_ffi_available)");
    build();
}

#[cfg(all(feature = "ffi", not(windows)))]
fn build() {
    use cmake::Config;
    use std::path::Path;

    // The C++ adapter is a repository-only conformance diagnostic.  crates.io
    // packages intentionally omit the large adjacent upstream source tree, so
    // keep `--all-features` usable there without pretending that the adapter
    // was linked.
    if !Path::new("diamond/CMakeLists.txt").is_file() {
        println!(
            "cargo:warning=the `ffi` feature requires the adjacent `diamond/` source tree; building without the C++ diagnostic adapter"
        );
        println!("cargo:rerun-if-changed=build.rs");
        return;
    }

    let dst = Config::new("diamond")
        .build_target("diamond_core")
        .define("CMAKE_BUILD_TYPE", "Release")
        .build();

    println!("cargo:rustc-link-search=native={}/build", dst.display());
    println!("cargo:rustc-link-lib=static=diamond_core");

    // Link C++ standard library
    #[cfg(target_os = "linux")]
    println!("cargo:rustc-link-lib=stdc++");
    #[cfg(target_os = "macos")]
    println!("cargo:rustc-link-lib=c++");

    // System dependencies that DIAMOND requires
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=pthread");
    println!("cargo:rustc-link-lib=dl");
    println!("cargo:rustc-cfg=diamond_ffi_available");

    // Rerun if C++ sources change
    println!("cargo:rerun-if-changed=diamond/src/");
    println!("cargo:rerun-if-changed=diamond/CMakeLists.txt");
}

#[cfg(any(not(feature = "ffi"), windows))]
fn build() {
    println!("cargo:rerun-if-changed=build.rs");
}
