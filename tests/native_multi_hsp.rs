//! Upstream snapshots with three disjoint reordered domains per target.

use std::path::Path;
use std::process::Command;

const FIELDS: &[&str] = &[
    "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
    "send", "evalue", "bitscore", "score", "qlen", "slen", "nident", "positive", "gaps", "qframe",
];

#[test]
fn multiple_hsps_match_upstream_limits_ordering_frames_and_archive_bytes() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/multi_hsp");
    let root = std::env::temp_dir().join(format!("diamond-native-hsps-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    let db = root.join("db");
    let status = Command::new(binary)
        .arg("makedb")
        .arg("--in")
        .arg(fixtures.join("reference.faa"))
        .arg("-d")
        .arg(&db)
        .arg("-p1")
        .output()
        .unwrap();
    assert!(
        status.status.success(),
        "{}",
        String::from_utf8_lossy(&status.stderr)
    );

    for (command, extension) in [("blastp", "faa"), ("blastx", "fna")] {
        for cbs in ["0", "1"] {
            for maximum in ["0", "1", "2", "3"] {
                let label = format!("{command}.hsps{maximum}.cbs{cbs}");
                let expected = std::fs::read(fixtures.join(format!("{label}.tsv"))).unwrap();
                let rows = expected.iter().filter(|&&byte| byte == b'\n').count();
                assert_eq!(
                    rows,
                    (if command == "blastx" { 2 } else { 1 })
                        * match maximum {
                            "1" => 2,
                            "2" => 3,
                            _ => 4,
                        },
                    "{label}"
                );
                for (format, memory, threads) in
                    [("tsv", "16G", "1"), ("tsv", "0G", "4"), ("daa", "16G", "1")]
                {
                    let output = root.join(format!("output.{format}"));
                    let mut process = Command::new(binary);
                    process
                        .arg(command)
                        .arg("-q")
                        .arg(fixtures.join(format!("query.{extension}")))
                        .arg("-d")
                        .arg(&db)
                        .args([
                            "--masking",
                            "0",
                            "--comp-based-stats",
                            cbs,
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
                    let status = process.output().unwrap();
                    assert!(
                        status.status.success(),
                        "{label} {format} {memory}: {}",
                        String::from_utf8_lossy(&status.stderr)
                    );
                    let expected = if format == "daa" {
                        std::fs::read(fixtures.join(format!("{label}.daa"))).unwrap()
                    } else {
                        expected.clone()
                    };
                    let actual = std::fs::read(output).unwrap();
                    assert!(
                        actual == expected,
                        "{label} {format} {memory}: wire bytes differ ({} vs {} bytes)",
                        actual.len(),
                        expected.len()
                    );
                }
            }
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn malformed_max_hsps_is_rejected_without_cpp_or_input_files() {
    for command in ["blastp", "blastx"] {
        for (arguments, expected) in [
            (
                vec!["--max-hsps"],
                "Invalid parameter count for option '--max-hsps'",
            ),
            (
                vec!["--max-hsps", "1", "2"],
                "Invalid parameter count for option '--max-hsps'",
            ),
            (vec!["--max-hsps=2"], "Invalid option: max-hsps=2"),
        ] {
            let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
                .arg(command)
                .args(arguments)
                .output()
                .unwrap();
            assert!(!result.status.success());
            let error = String::from_utf8_lossy(&result.stderr);
            assert!(error.contains(expected), "{error}");
            assert!(!error.contains("FFI"), "{error}");
        }
    }
}

#[test]
fn full_matrix_alternatives_preserve_banded_ties_and_transcript_statistics() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/multi_hsp");
    let root = std::env::temp_dir().join(format!("diamond-hsp-ties-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    let db = root.join("db");
    let build = Command::new(binary)
        .arg("makedb")
        .arg("--in")
        .arg(fixtures.join("reduced-reference.faa"))
        .arg("-d")
        .arg(&db)
        .arg("-p1")
        .output()
        .unwrap();
    assert!(
        build.status.success(),
        "{}",
        String::from_utf8_lossy(&build.stderr)
    );
    for format in ["tsv", "daa"] {
        let output = root.join(format!("output.{format}"));
        let mut command = Command::new(binary);
        command
            .arg("blastp")
            .arg("-q")
            .arg(fixtures.join("reduced-query.faa"))
            .arg("-d")
            .arg(&db)
            .args([
                "--more-sensitive",
                "--comp-based-stats",
                "1",
                "--max-hsps",
                "0",
                "-p1",
                "--outfmt",
            ]);
        if format == "daa" {
            command.arg("100");
        } else {
            command.args(FIELDS);
        }
        let result = command.arg("-o").arg(&output).output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(
            std::fs::read(output).unwrap()
                == std::fs::read(fixtures.join(format!("reduced.{format}"))).unwrap(),
            "{format} differs from upstream"
        );
    }
    std::fs::remove_dir_all(root).unwrap();
}
