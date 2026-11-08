//! Standard upstream composition modes, without the optional C++ adapter.
use std::path::Path;
use std::process::Command;

const FIELDS: &[&str] = &[
    "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
    "send", "evalue", "bitscore", "score", "nident", "positive", "gaps",
];

#[test]
fn composition_modes_match_upstream_scores_statistics_and_archives() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let root = std::env::temp_dir().join(format!("diamond-native-cbs-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    let db = root.join("db");
    let result = Command::new(binary)
        .args(["makedb", "--in"])
        .arg(fixtures.join("multi_hsp/reference.faa"))
        .arg("-d")
        .arg(&db)
        .arg("-p1")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_ne!(
        std::fs::read(fixtures.join("cbs/compact.cbs0.tsv")).unwrap(),
        std::fs::read(fixtures.join("cbs/compact.cbs4.tsv")).unwrap(),
        "the fixture must exercise adjusted scoring"
    );
    for mode in 0..6 {
        for maximum in ["1", "0"] {
            let label = if maximum == "0" {
                format!("compact.max0.cbs{mode}")
            } else {
                format!("compact.cbs{mode}")
            };
            for (format, memory, threads) in
                [("tsv", "16G", "1"), ("tsv", "0G", "4"), ("daa", "16G", "1")]
            {
                let output = root.join(format!("output.{format}"));
                let mut process = Command::new(binary);
                process
                    .args(["blastp", "-q"])
                    .arg(fixtures.join("multi_hsp/query.faa"))
                    .arg("-d")
                    .arg(&db)
                    .args([
                        "--masking",
                        "0",
                        "--comp-based-stats",
                        &mode.to_string(),
                        "--max-hsps",
                        maximum,
                        "--memory-limit",
                        memory,
                        "-p",
                        threads,
                    ])
                    .arg("-o")
                    .arg(&output)
                    .arg("--outfmt");
                if format == "daa" {
                    process.arg("100");
                } else {
                    process.args(FIELDS);
                }
                let result = process.output().unwrap();
                assert!(
                    result.status.success(),
                    "{label} {format} {memory}: {}",
                    String::from_utf8_lossy(&result.stderr)
                );
                let expected =
                    std::fs::read(fixtures.join(format!("cbs/{label}.{format}"))).unwrap();
                assert_eq!(
                    std::fs::read(&output).unwrap(),
                    expected,
                    "{label} {format} {memory}"
                );
            }
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn unsupported_composition_modes_fail_before_input_io() {
    let binary = env!("CARGO_BIN_EXE_diamond");
    let invalid = "Invalid value for --comp-based-stats. Permitted values: 0, 1, 2, 3, 4, 5.";
    for (command, mode, extra, message) in [
        ("blastp", "6", None, invalid),
        ("blastp", "-1", None, invalid),
        (
            "blastx",
            "2",
            None,
            "This mode of composition based stats is not supported for translated searches.",
        ),
        (
            "blastx",
            "5",
            None,
            "This mode of composition based stats is not supported for translated searches.",
        ),
        (
            "blastp",
            "4",
            Some("1"),
            "Global ranking is not supported in this mode.",
        ),
    ] {
        let mut process = Command::new(binary);
        process.args([command, "--comp-based-stats", mode]);
        if let Some(ranking) = extra {
            process.args(["--global-ranking", ranking]);
        }
        let result = process.output().unwrap();
        assert!(!result.status.success());
        assert_eq!(
            String::from_utf8_lossy(&result.stderr).trim(),
            format!("Error: {message}")
        );
    }
}

#[test]
fn frequency_ratio_starting_matrix_matches_natural_upstream_regression() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cbs");
    let root =
        std::env::temp_dir().join(format!("diamond-native-cbs-scaled-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    let db = root.join("db");
    let result = Command::new(binary)
        .args(["makedb", "--in"])
        .arg(fixtures.join("reference.faa"))
        .arg("-d")
        .arg(&db)
        .arg("-p1")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    for mode in 0..6 {
        for (format, memory) in [("tsv", "16G"), ("tsv", "0G"), ("daa", "16G")] {
            let output = root.join(format!("output.{format}"));
            let mut process = Command::new(binary);
            process
                .args(["blastp", "-q"])
                .arg(fixtures.join("query.faa"))
                .arg("-d")
                .arg(&db)
                .args([
                    "--more-sensitive",
                    "--masking",
                    "0",
                    "--comp-based-stats",
                    &mode.to_string(),
                    "--memory-limit",
                    memory,
                    "-p4",
                ])
                .arg("-o")
                .arg(&output)
                .arg("--outfmt");
            if format == "daa" {
                process.arg("100");
            } else {
                process.args(FIELDS);
            }
            let result = process.output().unwrap();
            assert!(
                result.status.success(),
                "mode {mode} {format} {memory}: {}",
                String::from_utf8_lossy(&result.stderr)
            );
            let expected =
                std::fs::read(fixtures.join(format!("scaled.cbs{mode}.{format}"))).unwrap();
            assert!(
                std::fs::read(&output).unwrap() == expected,
                "mode {mode} {format} {memory}: bytes differ"
            );
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}
