//! Tests for blastp rows where Rust output diverges from C++.
//!
//! Inputs (kept outside the repo because they're large):
//!   - `/tmp/bench_queries.faa`  — 500 mouse proteins from `uniprot.fa`
//!   - `/tmp/bench_rust.dmnd`   — diamond DB built by Rust `makedb`
//!
//! Reference output (snapshot of C++ DIAMOND on the same inputs):
//!   - `tests/data/cpp_blastp_mouse_500q.tsv`
//!
//! As of the snapshot date, Rust produced 6242 alignment rows matching C++
//! on 6105 of them (97.8%). The remaining 137 rows share the same query/
//! subject coordinates and identity counts but report different
//! evalue/bitscore values — suggesting a scoring (not alignment-path)
//! divergence, likely in CBS bias accumulation.
//!
//! These tests are `#[ignore]` by default because they require the
//! benchmark `.dmnd` and query FASTA files at fixed `/tmp/` paths. Run
//! them with:
//!
//! ```
//! cargo test --release --test divergent_blastp_rows -- --ignored --nocapture
//! ```
//!
//! Each test asserts that Rust's bitscore/evalue for one specific
//! (qseqid, sseqid) pair matches the C++ snapshot. Tests will fail today
//! and pass once the underlying scoring divergence is resolved.

use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::process::Command;

const RUST_BIN: &str = env!("CARGO_BIN_EXE_diamond");
const QUERIES: &str = "/tmp/bench_queries.faa";
const DB: &str = "/tmp/bench_rust.dmnd";
const CPP_REF: &str = "tests/data/cpp_blastp_mouse_500q.tsv";

fn snapshot_path() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(CPP_REF)
}

fn bench_inputs_present() -> bool {
    Path::new(QUERIES).exists() && Path::new(DB).exists() && snapshot_path().exists()
}

/// Run Rust blastp on the benchmark and return parsed rows keyed by
/// `(qseqid, sseqid)`. Each value is the full row of tab-separated fields.
fn run_rust_blastp() -> HashMap<(String, String), Vec<String>> {
    let out = std::env::temp_dir().join("divergent_blastp_rust.tsv");
    let status = Command::new(RUST_BIN)
        .args([
            "blastp",
            "-q",
            QUERIES,
            "-d",
            DB,
            "-o",
            out.to_str().unwrap(),
            "-f",
            "6",
            "--threads",
            "8",
        ])
        .status()
        .expect("failed to spawn rust diamond");
    assert!(status.success(), "rust diamond exited {:?}", status);
    parse_tsv(&out)
}

fn parse_tsv(path: &Path) -> HashMap<(String, String), Vec<String>> {
    let content = fs::read_to_string(path).expect("read tsv");
    let mut out = HashMap::new();
    for line in content.lines() {
        if line.is_empty() {
            continue;
        }
        let cols: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if cols.len() < 12 {
            continue;
        }
        out.insert((cols[0].clone(), cols[1].clone()), cols);
    }
    out
}

/// Assert that Rust's row for `(qseqid, sseqid)` matches the C++ snapshot
/// field-for-field. Fields 11 (evalue) and 12 (bitscore) are the typical
/// disagreement sites; the rest match.
fn assert_row_matches(
    rust: &HashMap<(String, String), Vec<String>>,
    cpp: &HashMap<(String, String), Vec<String>>,
    qseqid: &str,
    sseqid: &str,
) {
    let key = (qseqid.to_string(), sseqid.to_string());
    let r = rust
        .get(&key)
        .unwrap_or_else(|| panic!("Rust missing row {}/{}", qseqid, sseqid));
    let c = cpp
        .get(&key)
        .unwrap_or_else(|| panic!("C++ snapshot missing row {}/{}", qseqid, sseqid));
    assert_eq!(
        r,
        c,
        "Row {}/{} diverges from C++:\n  rust: {}\n  cpp:  {}",
        qseqid,
        sseqid,
        r.join("\t"),
        c.join("\t")
    );
}

/// Each #[test] below pins one (query, subject) pair where Rust currently
/// diverges from C++. The tests share a one-shot run of Rust blastp; we
/// use `std::sync::OnceLock` so it executes at most once even when
/// `--test-threads=8`.
static RUST_OUTPUT: std::sync::OnceLock<HashMap<(String, String), Vec<String>>> =
    std::sync::OnceLock::new();
static CPP_OUTPUT: std::sync::OnceLock<HashMap<(String, String), Vec<String>>> =
    std::sync::OnceLock::new();

fn rust_output() -> &'static HashMap<(String, String), Vec<String>> {
    RUST_OUTPUT.get_or_init(run_rust_blastp)
}

fn cpp_output() -> &'static HashMap<(String, String), Vec<String>> {
    CPP_OUTPUT.get_or_init(|| parse_tsv(&snapshot_path()))
}

macro_rules! divergent_row_test {
    ($name:ident, $q:expr, $s:expr) => {
        #[test]
        #[ignore = "requires bench data at /tmp/bench_*; run with --ignored"]
        fn $name() {
            if !bench_inputs_present() {
                panic!(
                    "missing /tmp/bench_queries.faa, /tmp/bench_rust.dmnd, or {}",
                    CPP_REF
                );
            }
            assert_row_matches(rust_output(), cpp_output(), $q, $s);
        }
    };
}

// Self-blast Chd5: identical sequence, no mismatches, no gaps — yet the
// bitscore diverges (Rust 3572 vs C++ 3624, ~52 bits). Almost certainly a
// CBS accumulation difference since the raw alignment path is trivial.
divergent_row_test!(chd5_self_blast, "Chd5", "Chd5");

// Chd5 vs Chd1/Chd1l/Chd2/Chd3/Chd4 family — all share the same coords
// and identity counts between Rust and C++ but report different bitscores.
divergent_row_test!(chd5_vs_chd1, "Chd5", "Chd1");
divergent_row_test!(chd5_vs_chd1l, "Chd5", "Chd1l");
divergent_row_test!(chd5_vs_chd2, "Chd5", "Chd2");
divergent_row_test!(chd5_vs_chd3, "Chd5", "Chd3");
divergent_row_test!(chd5_vs_chd4, "Chd5", "Chd4");

// Btaf1 — cross-family hit; same path, different score.
divergent_row_test!(chd5_vs_btaf1, "Chd5", "Btaf1");

// Olfactory receptors — short sequences exposing the short-query path.
divergent_row_test!(olfr1024_self_or_family, "Olfr1024", "Olfr1024");
divergent_row_test!(olfr1490_family, "Olfr1490", "Olfr1490");

// Long proteins: titin (Ttn), obscurin (Obscn) — stress test for
// chunked / banded scoring on multi-kb queries.
divergent_row_test!(ttn_self, "Ttn", "Ttn");
divergent_row_test!(obscn_self, "Obscn", "Obscn");

/// Summary test: enumerates ALL pairs where Rust diverges from C++ and
/// prints them. Always asserts the divergence count is <= 137 (the
/// known-bad count at the snapshot date) — so this test will start
/// FAILING in the GOOD direction if a fix reduces divergence below that
/// threshold (prompting an update to lower the bound) or in the BAD
/// direction if a regression widens it.
#[test]
#[ignore = "requires bench data at /tmp/bench_*; run with --ignored"]
fn count_divergent_rows_against_cpp_snapshot() {
    if !bench_inputs_present() {
        panic!(
            "missing /tmp/bench_queries.faa, /tmp/bench_rust.dmnd, or {}",
            CPP_REF
        );
    }
    let rust = rust_output();
    let cpp = cpp_output();

    let mut divergent: Vec<((String, String), Vec<String>, Vec<String>)> = Vec::new();
    for (key, r) in rust {
        if let Some(c) = cpp.get(key) {
            if r != c {
                divergent.push((key.clone(), r.clone(), c.clone()));
            }
        }
    }
    divergent.sort_by(|a, b| a.0.cmp(&b.0));

    eprintln!(
        "{} divergent rows out of {} shared (rust)",
        divergent.len(),
        rust.len()
    );
    for ((q, s), r, c) in divergent.iter().take(10) {
        eprintln!("  {}/{}", q, s);
        eprintln!("    rust: {}", r.join("\t"));
        eprintln!("    cpp:  {}", c.join("\t"));
    }
    if divergent.len() > 10 {
        eprintln!("  ... ({} more)", divergent.len() - 10);
    }

    assert!(
        divergent.len() <= 137,
        "Divergent row count grew to {} (snapshot baseline: 137). \
         Investigate the regression.",
        divergent.len()
    );
    assert_eq!(
        divergent.len(),
        0,
        "Rust output still diverges from C++ on {} rows. Drop this assertion \
         and lower `<= 137` once the underlying scoring divergence is fixed.",
        divergent.len()
    );
}
