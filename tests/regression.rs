use std::fs;
use std::path::{Path, PathBuf};

fn test_data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("diamond")
        .join("src")
        .join("test")
}

fn run_diamond(args: &[&str]) -> i32 {
    diamond::ffi::run(args)
}

/// Run the built-in regression test suite (20 test cases with hash comparison).
/// This is the most comprehensive test — it exercises blastp with various
/// sensitivity levels, output formats, comp-based-stats modes, etc.
#[test]
fn builtin_regression_tests() {
    let code = run_diamond(&["diamond", "test"]);
    assert_eq!(code, 0, "Built-in regression tests failed (diamond test returned {})", code);
}

/// Helper: run a diamond command, write output to a temp file, compare with reference.
fn run_and_compare(name: &str, args: &str) {
    let td = test_data_dir();
    let out_dir = std::env::temp_dir().join(format!("diamond_test_{}", name));
    fs::create_dir_all(&out_dir).unwrap();
    let out_file = out_dir.join(format!("{}.out", name));

    let mut full_args: Vec<&str> = vec!["diamond"];
    let arg_parts: Vec<String> = args
        .split_whitespace()
        .map(|s| s.replace("${TD}", &td.to_string_lossy()))
        .collect();
    let arg_refs: Vec<&str> = arg_parts.iter().map(|s| s.as_str()).collect();
    full_args.extend_from_slice(&arg_refs);
    full_args.push("-o");
    let out_str = out_file.to_string_lossy().to_string();
    full_args.push(&out_str);

    let code = run_diamond(&full_args);
    assert_eq!(code, 0, "{}: diamond exited with code {}", name, code);

    let expected_path = td.join(format!("{}.out", name));
    let expected = fs::read(&expected_path)
        .unwrap_or_else(|e| panic!("{}: failed to read reference file {:?}: {}", name, expected_path, e));
    let actual = fs::read(&out_file)
        .unwrap_or_else(|e| panic!("{}: failed to read output file {:?}: {}", name, out_file, e));

    assert_eq!(
        actual, expected,
        "{}: output differs from reference.\nExpected {} bytes, got {} bytes",
        name,
        expected.len(),
        actual.len()
    );

    // Clean up
    let _ = fs::remove_dir_all(&out_dir);
}

#[test]
fn ctest_blastp() {
    run_and_compare(
        "blastp",
        "blastp -q ${TD}/1.faa -d ${TD}/2.faa -p1",
    );
}

#[test]
fn ctest_blastp_mid_sens() {
    run_and_compare(
        "blastp-mid-sens",
        "blastp -q ${TD}/3.faa -d ${TD}/4.faa --mid-sensitive -p1",
    );
}

#[test]
fn ctest_blastp_f0() {
    run_and_compare(
        "blastp-f0",
        "blastp -q ${TD}/1.faa -d ${TD}/2.faa -f0 -p1",
    );
}

#[test]
fn ctest_galaxy_7() {
    run_and_compare(
        "galaxy_7",
        "blastx --threads 2 --db ${TD}/galaxy/db.dmnd --query ${TD}/galaxy/nucleotide.fasta --query-gencode 1 --strand both --min-orf 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --header 0 --compress 0 --matrix BLOSUM62 --comp-based-stats 1 --masking tantan --max-target-seqs 25 --evalue 0.001 --id 0.0 --approx-id 0.0 --query-cover 0.0 --subject-cover 0.0 --block-size 2.0 --motif-masking 0 --soft-masking 0 --swipe --algo 0 --index-chunks 4 --file-buffer-size 67108864",
    );
}

#[test]
fn ctest_galaxy_9() {
    run_and_compare(
        "galaxy_9",
        "blastx --threads 2 --db ${TD}/galaxy/db.dmnd --query ${TD}/galaxy/nucleotide.fasta --query-gencode 1 --strand both --min-orf 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --header 0 --compress 0 --matrix BLOSUM62 --comp-based-stats 1 --masking tantan --max-target-seqs 25 --evalue 0.001 --id 0.0 --approx-id 0.0 --query-cover 0.0 --subject-cover 0.0 --block-size 2.0 --motif-masking 0 --soft-masking 0 --algo 0 --global-ranking 10 --index-chunks 4 --file-buffer-size 67108864",
    );
}
