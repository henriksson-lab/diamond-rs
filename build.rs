use cmake::Config;

fn main() {
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
    println!("cargo:rustc-link-lib=sqlite3");
    println!("cargo:rustc-link-lib=dl");

    // Rerun if C++ sources change
    println!("cargo:rerun-if-changed=diamond/src/");
    println!("cargo:rerun-if-changed=diamond/CMakeLists.txt");
}
