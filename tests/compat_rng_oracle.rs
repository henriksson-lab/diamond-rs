//! Live cross-platform comparison with the target C/C++ runtime.
//!
//! This test is ignored by default because it invokes a C++ compiler and needs
//! the pinned upstream checkout at `diamond/`. Dedicated CI runs it on each
//! supported OS; normal Rust builds remain independent of a C++ toolchain.

#[path = "../src/util/compat_rng.rs"]
mod compat_rng;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use diamond::basic::value::AMINO_ACID_COUNT;
use diamond::stats::matrices::BLOSUM62;
use diamond::stats::pvalues::pvalues;
use diamond::stats::sls_alignment_evaluer::AlignmentEvaluer;
use diamond::stats::sls_basic::AlignmentEvaluerParametersWithErrors;
use diamond::tantan::LambdaCalculator;

#[test]
#[ignore = "opt-in test compiles the pinned target-native C++ RNG oracle"]
fn target_c_runtime_and_downstream_paths_match_pinned_cpp() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let upstream = root.join("diamond/src");
    assert!(
        upstream.join("lib/tantan/LambdaCalculator.cc").is_file(),
        "pinned upstream checkout missing at {}",
        upstream.display()
    );

    let output_dir = std::env::temp_dir().join(format!(
        "diamond-rng-oracle-{}-{}",
        std::process::id(),
        std::env::consts::OS
    ));
    std::fs::create_dir_all(&output_dir).unwrap();
    let executable = compile_probe(&root, &upstream, &output_dir);
    let output = Command::new(&executable).output().unwrap_or_else(|error| {
        panic!(
            "failed to execute RNG oracle {}: {error}",
            executable.display()
        )
    });
    assert!(
        output.status.success(),
        "RNG oracle failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let oracle_stdout = String::from_utf8(output.stdout).unwrap();
    let vectors = parse_vectors(&oracle_stdout);

    compare_rand(&vectors, 1);
    compare_rand(&vectors, 12_345);

    compat_rng::c_srand(12_345);
    let rust_normals = (0..16)
        .map(|_| pvalues::standard_normal())
        .collect::<Vec<_>>();
    compare_float_vectors(
        "standard-normal",
        &rust_normals,
        vector(&vectors, "standard-normal"),
        3.0e-14,
    );

    compat_rng::c_srand(1);
    let mut calculator = LambdaCalculator::new();
    calculator.calculate_flat_i8(&BLOSUM62.scores, AMINO_ACID_COUNT, 20);
    let mut rust_tantan = Vec::with_capacity(41);
    rust_tantan.push(calculator.lambda());
    rust_tantan.extend_from_slice(calculator.letter_probs1().unwrap());
    rust_tantan.extend_from_slice(calculator.letter_probs2().unwrap());
    compare_float_vectors("tantan", &rust_tantan, vector(&vectors, "tantan"), 3.0e-13);

    let mut evaluator = AlignmentEvaluer::new();
    evaluator
        .init_parameters_with_errors(&parameters_with_errors())
        .unwrap();
    let mut rust_evaluator = Vec::with_capacity(40);
    rust_evaluator.extend_from_slice(&evaluator.parameters().m_LambdaSbs);
    rust_evaluator.extend_from_slice(&evaluator.parameters().m_TauSbs);
    compare_float_vectors(
        "evaluator",
        &rust_evaluator,
        vector(&vectors, "evaluator"),
        3.0e-13,
    );

    // Write retained evidence only after every comparison succeeds. A partial
    // oracle transcript from a failed test must not look like a conformance
    // artifact. Dedicated CI sets an artifact path; local runs retain the same
    // validated capture beside the temporary oracle executable.
    let capture_path = std::env::var_os("DIAMOND_RNG_ORACLE_CAPTURE")
        .map(PathBuf::from)
        .unwrap_or_else(|| output_dir.join("compat-rng-oracle.txt"));
    write_capture(&capture_path, &oracle_stdout);
}

fn write_capture(path: &Path, cpp_output: &str) {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).unwrap();
    }

    let mut capture = format!(
        "target-os {}\ntarget-arch {}\ntarget-env {}\nstatus passed\n\n# Target-native C++ oracle\n{}",
        std::env::consts::OS,
        std::env::consts::ARCH,
        target_env(),
        cpp_output
    );
    for seed in [1_u32, 12_345] {
        compat_rng::c_srand(seed);
        capture.push_str(&format!("rust-rand-{seed}"));
        for _ in 0..100 {
            capture.push_str(&format!(" {}", compat_rng::c_rand()));
        }
        capture.push('\n');
    }
    std::fs::write(path, capture).unwrap();
    eprintln!("RNG conformance vectors captured at {}", path.display());
}

fn target_env() -> &'static str {
    if cfg!(target_env = "gnu") {
        "gnu"
    } else if cfg!(target_env = "musl") {
        "musl"
    } else if cfg!(target_env = "msvc") {
        "msvc"
    } else {
        "none"
    }
}

fn compile_probe(root: &Path, upstream: &Path, output_dir: &Path) -> PathBuf {
    let source = root.join("tests/fixtures/compat_rng_oracle.cpp");
    let lambda = upstream.join("lib/tantan/LambdaCalculator.cc");
    let executable = output_dir.join(if cfg!(windows) {
        "compat_rng_oracle.exe"
    } else {
        "compat_rng_oracle"
    });
    let compiler = std::env::var_os("CXX").unwrap_or_else(|| {
        if cfg!(target_env = "msvc") {
            "cl.exe".into()
        } else {
            "c++".into()
        }
    });
    let mut command = Command::new(&compiler);
    if cfg!(target_env = "msvc") {
        command
            .current_dir(output_dir)
            .args(["/nologo", "/EHsc", "/std:c++17", "/O2"])
            .arg(format!("/I{}", upstream.display()))
            .arg(&source)
            .arg(&lambda)
            .arg(format!("/Fe:{}", executable.display()));
    } else {
        command
            .args(["-std=c++17", "-O2", "-I"])
            .arg(upstream)
            .arg(&source)
            .arg(&lambda)
            .arg("-o")
            .arg(&executable);
    }
    let output = command.output().unwrap_or_else(|error| {
        panic!(
            "failed to invoke target C++ compiler {:?}: {error}",
            compiler
        )
    });
    assert!(
        output.status.success(),
        "C++ RNG oracle compilation failed\ncommand: {command:?}\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    executable
}

fn parse_vectors(output: &str) -> BTreeMap<String, Vec<String>> {
    output
        .lines()
        .map(|line| {
            let mut fields = line.split_ascii_whitespace();
            let label = fields.next().expect("oracle emitted an empty line");
            (label.to_owned(), fields.map(str::to_owned).collect())
        })
        .collect()
}

fn vector<'a>(vectors: &'a BTreeMap<String, Vec<String>>, label: &str) -> &'a [String] {
    vectors
        .get(label)
        .unwrap_or_else(|| panic!("C++ oracle did not emit {label}"))
}

fn compare_rand(vectors: &BTreeMap<String, Vec<String>>, seed: u32) {
    let label = format!("rand-{seed}");
    let cpp = vector(vectors, &label)
        .iter()
        .map(|value| value.parse::<i32>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(cpp.len(), 100, "{label} C++ vector length");
    compat_rng::c_srand(seed);
    let rust = (0..100).map(|_| compat_rng::c_rand()).collect::<Vec<_>>();
    assert_eq!(rust, cpp, "{label} differs from the target C runtime");
}

fn compare_float_vectors(label: &str, rust: &[f64], cpp: &[String], tolerance: f64) {
    let cpp = cpp
        .iter()
        .map(|value| value.parse::<f64>().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(rust.len(), cpp.len(), "{label} vector length");
    for (index, (&rust, &cpp)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (rust - cpp).abs() <= tolerance,
            "{label}[{index}] differs: Rust {rust:?}, C++ {cpp:?}"
        );
    }
}

fn parameters_with_errors() -> AlignmentEvaluerParametersWithErrors {
    AlignmentEvaluerParametersWithErrors {
        d_lambda: 0.267,
        d_lambda_error: 0.001,
        d_k: 0.041,
        d_k_error: 0.001,
        d_a1: 1.9,
        d_a1_error: 0.01,
        d_b1: 4.0,
        d_b1_error: 0.01,
        d_a2: 2.1,
        d_a2_error: 0.01,
        d_b2: 5.0,
        d_b2_error: 0.01,
        d_alpha1: 1.7,
        d_alpha1_error: 0.01,
        d_beta1: 6.0,
        d_beta1_error: 0.01,
        d_alpha2: 1.8,
        d_alpha2_error: 0.01,
        d_beta2: 7.0,
        d_beta2_error: 0.01,
        d_sigma: 43.0,
        d_sigma_error: 0.01,
        d_tau: 8.0,
        d_tau_error: 0.01,
    }
}
