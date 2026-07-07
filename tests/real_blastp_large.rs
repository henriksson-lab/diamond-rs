//! Larger real-data native blastp conformance tests.
//!
//! These tests are ignored by default because they require the local benchmark
//! fixtures generated during translation work:
//!   - `/tmp/bench_queries.faa`
//!   - `/tmp/bench_rust.dmnd`
//!   - `/tmp/bench_cpp.dmnd` for dynamic C++ comparison tests
//!   - `tests/data/cpp_blastp_mouse_500q.tsv`
//!
//! For faster runs against the optimized binary:
//!
//! ```text
//! DIAMOND_RS_BIN=target/release/diamond cargo test --test real_blastp_large -- --ignored --nocapture
//! ```

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const QUERIES_500: &str = "/tmp/bench_queries.faa";
const DB: &str = "/tmp/bench_rust.dmnd";
const CPP_DB: &str = "/tmp/bench_cpp.dmnd";
const CPP_BIN: &str = "diamond/build/diamond";
const CPP_500_SNAPSHOT: &str = "tests/data/cpp_blastp_mouse_500q.tsv";

fn snapshot_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(CPP_500_SNAPSHOT)
}

fn require_large_inputs() {
    for path in [QUERIES_500, DB] {
        assert!(
            Path::new(path).exists(),
            "missing required large real-data fixture: {path}"
        );
    }
    assert!(
        snapshot_path().exists(),
        "missing C++ snapshot: {}",
        snapshot_path().display()
    );
}

fn diamond_bin() -> String {
    std::env::var("DIAMOND_RS_BIN").unwrap_or_else(|_| env!("CARGO_BIN_EXE_diamond").to_string())
}

fn cpp_diamond_bin() -> PathBuf {
    std::env::var("DIAMOND_CPP_BIN")
        .map(PathBuf::from)
        .unwrap_or_else(|_| Path::new(env!("CARGO_MANIFEST_DIR")).join(CPP_BIN))
}

fn require_cpp_comparison_inputs() {
    require_large_inputs();
    assert!(
        Path::new(CPP_DB).exists(),
        "missing required C++ comparison database: {CPP_DB}"
    );
    assert!(
        cpp_diamond_bin().exists(),
        "missing C++ diamond binary: {}",
        cpp_diamond_bin().display()
    );
}

fn first_n_fasta_records(input: &Path, n: usize, output: &Path) -> BTreeSet<String> {
    let fasta = fs::read_to_string(input).expect("read FASTA");
    let mut ids = BTreeSet::new();
    let mut out = String::new();
    let mut records = 0usize;
    let mut keep = false;

    for line in fasta.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if records == n {
                break;
            }
            records += 1;
            keep = true;
            let id = header.split_whitespace().next().unwrap_or(header);
            ids.insert(id.to_string());
        }
        if keep {
            out.push_str(line);
            out.push('\n');
        }
    }

    assert_eq!(records, n, "expected at least {n} FASTA records");
    fs::write(output, out).expect("write FASTA subset");
    ids
}

fn fasta_ids(input: &Path) -> BTreeSet<String> {
    fs::read_to_string(input)
        .expect("read FASTA")
        .lines()
        .filter_map(|line| line.strip_prefix('>'))
        .map(|header| {
            header
                .split_whitespace()
                .next()
                .unwrap_or(header)
                .to_string()
        })
        .collect()
}

fn expected_snapshot_rows(ids: &BTreeSet<String>) -> String {
    fs::read_to_string(snapshot_path())
        .expect("read C++ snapshot")
        .lines()
        .filter(|line| {
            line.split('\t')
                .next()
                .is_some_and(|query_id| ids.contains(query_id))
        })
        .map(|line| {
            let mut s = String::with_capacity(line.len() + 1);
            s.push_str(line);
            s.push('\n');
            s
        })
        .collect()
}

fn run_native_blastp(query: &Path, output: &Path) {
    run_blastp_with_args(
        &diamond_bin(),
        query,
        DB,
        output,
        &[],
        "native diamond blastp",
    );
}

fn run_blastp_with_args(
    bin: impl AsRef<std::ffi::OsStr>,
    query: &Path,
    db: &str,
    output: &Path,
    extra_args: &[&str],
    label: &str,
) {
    let status = Command::new(bin)
        .args([
            "blastp",
            "-q",
            query.to_str().expect("query path"),
            "-d",
            db,
            "-p1",
            "-o",
            output.to_str().expect("output path"),
        ])
        .args(extra_args)
        .status()
        .unwrap_or_else(|e| panic!("run {label}: {e}"));
    assert!(status.success(), "{label} exited with {status}");
}

fn assert_exact_snapshot_match(query: &Path, expected_ids: &BTreeSet<String>, label: &str) {
    let out = std::env::temp_dir().join(format!("diamond_rs_{label}.tsv"));
    run_native_blastp(query, &out);

    let expected = expected_snapshot_rows(expected_ids);
    let actual = fs::read_to_string(&out).expect("read native output");
    assert_eq!(
        actual.lines().count(),
        expected.lines().count(),
        "{label}: row count changed"
    );
    assert_eq!(actual, expected, "{label}: native output differs from C++");

    let _ = fs::remove_file(out);
}

fn assert_dynamic_cpp_match(query: &Path, extra_args: &[&str], label: &str) {
    let cpp_out = std::env::temp_dir().join(format!("diamond_cpp_{label}.tsv"));
    let rust_out = std::env::temp_dir().join(format!("diamond_rs_{label}.tsv"));

    run_blastp_with_args(
        cpp_diamond_bin(),
        query,
        CPP_DB,
        &cpp_out,
        extra_args,
        "C++ diamond blastp",
    );
    run_blastp_with_args(
        diamond_bin(),
        query,
        DB,
        &rust_out,
        extra_args,
        "native diamond blastp",
    );

    let expected = fs::read_to_string(&cpp_out).expect("read C++ output");
    let actual = fs::read_to_string(&rust_out).expect("read native output");
    assert_eq!(
        actual.lines().count(),
        expected.lines().count(),
        "{label}: row count changed"
    );
    assert_eq!(actual, expected, "{label}: native output differs from C++");

    let _ = fs::remove_file(cpp_out);
    let _ = fs::remove_file(rust_out);
}

fn assert_first_n_dynamic_cpp_match(n: usize, extra_args: &[&str], label: &str) {
    require_cpp_comparison_inputs();

    let query = std::env::temp_dir().join(format!("diamond_rs_{label}.faa"));
    first_n_fasta_records(Path::new(QUERIES_500), n, &query);
    assert_dynamic_cpp_match(&query, extra_args, label);
    let _ = fs::remove_file(query);
}

#[test]
#[ignore = "requires /tmp/bench_queries.faa and /tmp/bench_rust.dmnd; large real-data run"]
fn first_100_complete_mouse_queries_match_cpp_snapshot() {
    require_large_inputs();

    let query = std::env::temp_dir().join("diamond_rs_first100_complete.faa");
    let ids = first_n_fasta_records(Path::new(QUERIES_500), 100, &query);
    assert_exact_snapshot_match(&query, &ids, "first100_complete_mouse");
    let _ = fs::remove_file(query);
}

#[test]
#[ignore = "requires /tmp/bench_queries.faa and /tmp/bench_rust.dmnd; medium real-data run"]
fn first_250_complete_mouse_queries_match_cpp_snapshot() {
    require_large_inputs();

    let query = std::env::temp_dir().join("diamond_rs_first250_complete.faa");
    let ids = first_n_fasta_records(Path::new(QUERIES_500), 250, &query);
    assert_exact_snapshot_match(&query, &ids, "first250_complete_mouse");
    let _ = fs::remove_file(query);
}

#[test]
#[ignore = "requires /tmp/bench_queries.faa and /tmp/bench_rust.dmnd; full 500-query real-data run"]
fn full_500_mouse_queries_match_cpp_snapshot() {
    require_large_inputs();

    let query = Path::new(QUERIES_500);
    let ids = fasta_ids(query);
    assert_eq!(ids.len(), 500, "expected the 500-query mouse benchmark");
    assert_exact_snapshot_match(query, &ids, "full500_mouse");
}

#[test]
#[ignore = "requires /tmp/bench_* databases and C++ diamond; dynamic real-data comparison"]
fn first_100_comp_based_stats_0_matches_cpp() {
    assert_first_n_dynamic_cpp_match(100, &["--comp-based-stats", "0"], "first100_comp0");
}

#[test]
#[ignore = "requires /tmp/bench_* databases and C++ diamond; larger dynamic real-data comparison"]
fn first_250_comp_based_stats_0_matches_cpp() {
    assert_first_n_dynamic_cpp_match(250, &["--comp-based-stats", "0"], "first250_comp0");
}

#[test]
#[ignore = "requires /tmp/bench_* databases and C++ diamond; dynamic real-data comparison"]
fn first_25_blosum45_gap_14_2_matches_cpp() {
    assert_first_n_dynamic_cpp_match(
        25,
        &[
            "--matrix",
            "BLOSUM45",
            "--gapopen",
            "14",
            "--gapextend",
            "2",
        ],
        "first25_blosum45_gap14_2",
    );
}

#[test]
#[ignore = "requires /tmp/bench_* databases and C++ diamond; dynamic real-data comparison"]
fn first_100_fast_sensitivity_matches_cpp() {
    assert_first_n_dynamic_cpp_match(100, &["--fast"], "first100_fast");
}

#[test]
#[ignore = "requires /tmp/bench_* databases and C++ diamond; larger dynamic real-data comparison"]
fn first_250_fast_sensitivity_matches_cpp() {
    assert_first_n_dynamic_cpp_match(250, &["--fast"], "first250_fast");
}

#[test]
#[ignore = "requires /tmp/bench_* databases and C++ diamond; dynamic real-data comparison"]
fn first_25_mid_sensitive_matches_cpp() {
    require_cpp_comparison_inputs();

    let query = std::env::temp_dir().join("diamond_rs_first25_mid_sensitive.faa");
    first_n_fasta_records(Path::new(QUERIES_500), 25, &query);
    assert_dynamic_cpp_match(&query, &["--mid-sensitive"], "first25_mid_sensitive");
    let _ = fs::remove_file(query);
}

#[test]
#[ignore = "requires /tmp/bench_* databases and C++ diamond; dynamic real-data comparison"]
fn first_100_max_target_seqs_5_matches_cpp() {
    assert_first_n_dynamic_cpp_match(100, &["-k", "5"], "first100_k5");
}

#[test]
#[ignore = "requires /tmp/bench_* databases and C++ diamond; larger dynamic real-data comparison"]
fn first_250_max_target_seqs_5_matches_cpp() {
    assert_first_n_dynamic_cpp_match(250, &["-k", "5"], "first250_k5");
}

#[test]
#[ignore = "requires /tmp/bench_* databases and C++ diamond; dynamic real-data comparison"]
fn first_100_strict_evalue_matches_cpp() {
    require_cpp_comparison_inputs();

    let query = std::env::temp_dir().join("diamond_rs_first100_evalue_1e20.faa");
    first_n_fasta_records(Path::new(QUERIES_500), 100, &query);
    assert_dynamic_cpp_match(&query, &["-e", "1e-20"], "first100_evalue_1e20");
    let _ = fs::remove_file(query);
}

#[test]
#[ignore = "requires /tmp/bench_* databases and C++ diamond; dynamic real-data comparison"]
fn first_100_explicit_default_tabular_outfmt_matches_cpp() {
    assert_first_n_dynamic_cpp_match(
        100,
        &[
            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        ],
        "first100_explicit_default_outfmt",
    );
}

#[test]
#[ignore = "requires /tmp/bench_* databases and C++ diamond; larger dynamic real-data comparison"]
fn first_250_explicit_default_tabular_outfmt_matches_cpp() {
    assert_first_n_dynamic_cpp_match(
        250,
        &[
            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        ],
        "first250_explicit_default_outfmt",
    );
}
