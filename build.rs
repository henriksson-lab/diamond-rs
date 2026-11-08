fn main() {
    println!("cargo::rustc-check-cfg=cfg(diamond_ffi_available)");
    println!("cargo:rerun-if-changed=build.rs");
    build();
}

#[cfg(all(feature = "ffi", not(windows)))]
fn build() {
    use cmake::Config;
    use std::path::PathBuf;

    println!("cargo:rerun-if-env-changed=DIAMOND_FFI_SOURCE_DIR");
    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    if target_os == "windows" {
        println!("cargo:warning=the C++ diagnostic adapter is unavailable for Windows targets");
        return;
    }
    // An isolated oracle checkout can be used without building against or
    // altering the user's adjacent source tree.
    let source = std::env::var_os("DIAMOND_FFI_SOURCE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("diamond"));

    // The C++ adapter is a repository-only conformance diagnostic.  crates.io
    // packages intentionally omit the large adjacent upstream source tree, so
    // keep `--all-features` usable there without pretending that the adapter
    // was linked.
    if !source.join("CMakeLists.txt").is_file() {
        println!(
            "cargo:warning=the `ffi` feature requires C++ sources at {}; building without the C++ diagnostic adapter", source.display()
        );
        println!("cargo:rerun-if-changed=build.rs");
        return;
    }

    let on_off = |enabled| if enabled { "ON" } else { "OFF" };
    let dst = Config::new(&source)
        .build_target("diamond_core")
        .define("CMAKE_BUILD_TYPE", "Release")
        .define("EXTRA", on_off(cfg!(feature = "extra")))
        .define("WITH_DNA", on_off(cfg!(feature = "dna")))
        .define("WITH_MCL", on_off(cfg!(feature = "mcl")))
        .define("WITH_ZSTD", on_off(cfg!(feature = "zstd")))
        .build();

    println!("cargo:rustc-link-search=native={}/build", dst.display());
    println!("cargo:rustc-link-lib=static=diamond_core");

    // Link C++ standard library
    match target_os.as_str() {
        "linux" => println!("cargo:rustc-link-lib=stdc++"),
        "macos" => println!("cargo:rustc-link-lib=c++"),
        _ => {}
    }

    // System dependencies that DIAMOND requires
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=sqlite3");
    if cfg!(feature = "zstd") {
        println!("cargo:rustc-link-lib=zstd");
    }
    println!("cargo:rustc-link-lib=pthread");
    println!("cargo:rustc-link-lib=dl");
    println!("cargo:rustc-cfg=diamond_ffi_available");

    // Rerun if C++ sources change
    println!("cargo:rerun-if-changed={}", source.join("src").display());
    println!(
        "cargo:rerun-if-changed={}",
        source.join("CMakeLists.txt").display()
    );
}

#[cfg(any(not(feature = "ffi"), windows))]
fn build() {
    println!("cargo:rerun-if-changed=build.rs");
}
