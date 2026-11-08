#[path = "support/source_startup.rs"]
mod source_startup;

use std::path::Path;
use std::process::Command;

#[test]
fn native_transcript_formats_match_upstream_for_both_search_modes() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let fixtures = data.join("search_output");
    let root = std::env::temp_dir().join(format!("diamond-search-output-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    let variants: &[(&str, &[&str])] = &[
        ("cbs0-all", &["--comp-based-stats", "0", "--max-hsps", "0"]),
        ("cbs1-all", &["--comp-based-stats", "1", "--max-hsps", "0"]),
        (
            "coverage-unal0",
            &[
                "--comp-based-stats",
                "0",
                "--query-cover",
                "99",
                "--subject-cover",
                "99",
                "--unal",
                "0",
            ],
        ),
        (
            "coverage-unal1",
            &[
                "--comp-based-stats",
                "0",
                "--query-cover",
                "99",
                "--subject-cover",
                "99",
                "--unal",
                "1",
            ],
        ),
    ];
    for case in ["compact", "natural"] {
        let reference = if case == "compact" {
            fixtures.join("reference.faa")
        } else {
            data.join("reassign/natural.faa")
        };
        let db = root.join(format!("{case}.dmnd"));
        let result = Command::new(binary)
            .args(["makedb", "--in"])
            .arg(reference)
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
        for mode in ["blastp", "blastx"] {
            let query = fixtures.join(format!(
                "{}query.{}",
                if case == "natural" { "natural-" } else { "" },
                if mode == "blastp" { "faa" } else { "fna" }
            ));
            for format in ["0", "xml", "paf"] {
                for (variant, options) in variants {
                    if case == "natural" && *variant != "cbs1-all" {
                        continue;
                    }
                    let expected = std::fs::read(
                        fixtures.join(format!("{case}.{mode}.{format}.{variant}.out")),
                    )
                    .unwrap();
                    for (threads, memory) in [("1", "16G"), ("4", "0G")] {
                        let alias = match (format, threads) {
                            ("xml", "4") => "5",
                            ("paf", "4") => "103",
                            _ => format,
                        };
                        let output = root.join("output");
                        let result = Command::new(binary)
                            .arg(mode)
                            .arg("-q")
                            .arg(&query)
                            .arg("-d")
                            .arg(&db)
                            .args([
                                "-p",
                                threads,
                                "--memory-limit",
                                memory,
                                "--masking",
                                "0",
                                "-f",
                                alias,
                            ])
                            .args(*options)
                            .arg("-o")
                            .arg(&output)
                            .output()
                            .unwrap();
                        let label = format!("{case} {mode} {format} {variant} p{threads} {memory}");
                        assert!(
                            result.status.success(),
                            "{label}: {}",
                            String::from_utf8_lossy(&result.stderr)
                        );
                        let actual = std::fs::read_to_string(output).unwrap();
                        let actual = actual.replace(
                            &format!("<BlastOutput_db>{}</BlastOutput_db>", db.display()),
                            "<BlastOutput_db>@DATABASE@</BlastOutput_db>",
                        );
                        assert!(actual.as_bytes() == expected, "{label}: bytes differ");
                    }
                }
            }
        }
        if case == "compact" {
            let result = Command::new(binary)
                .arg("blastp")
                .arg("-q")
                .arg(fixtures.join("query.faa"))
                .arg("-d")
                .arg(&db)
                .args([
                    "-p1",
                    "-f",
                    "xml",
                    "--masking",
                    "0",
                    "--comp-based-stats",
                    "0",
                    "-e",
                    "1e-10",
                ])
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            let actual = String::from_utf8(result.stdout).unwrap().replace(
                &format!("<BlastOutput_db>{}</BlastOutput_db>", db.display()),
                "<BlastOutput_db>@DATABASE@</BlastOutput_db>",
            );
            assert_eq!(
                actual,
                std::fs::read_to_string(fixtures.join("small-evalue.xml")).unwrap()
            );
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn unaligned_option_cardinality_matches_the_upstream_parser() {
    for (args, expected) in [
        (
            vec!["--unal"],
            "Invalid parameter count for option '--unal'",
        ),
        (
            vec!["--unal", "1", "0"],
            "Invalid parameter count for option '--unal'",
        ),
        (vec!["--unal=1"], "Invalid option: unal=1"),
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .args(["blastp", "-f", "xml"])
            .args(args)
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stderr).ends_with(&format!("Error: {expected}\n")));
    }
}

#[test]
fn unseeded_queries_are_skipped_inside_nonempty_weighted_bins() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/search_output");
    let root = std::env::temp_dir().join(format!("diamond-output-bins-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    let db = root.join("db.dmnd");
    let result = Command::new(binary)
        .args(["makedb", "--in"])
        .arg(fixtures.join("reference.faa"))
        .arg("-d")
        .arg(&db)
        .arg("-p1")
        .output()
        .unwrap();
    assert!(result.status.success());
    for mode in ["blastp", "blastx"] {
        let query = fixtures.join(if mode == "blastp" {
            "bins-query.faa"
        } else {
            "bins-query.fna"
        });
        for format in ["0", "xml", "paf"] {
            let expected =
                std::fs::read(fixtures.join(format!("bins.{mode}.{format}.out"))).unwrap();
            for (threads, memory) in [("1", "16G"), ("4", "0G")] {
                let result = Command::new(binary)
                    .arg(mode)
                    .arg("-q")
                    .arg(&query)
                    .arg("-d")
                    .arg(&db)
                    .args([
                        "-p",
                        threads,
                        "--memory-limit",
                        memory,
                        "--masking",
                        "0",
                        "--comp-based-stats",
                        "0",
                        "--max-hsps",
                        "0",
                        "-f",
                        format,
                    ])
                    .output()
                    .unwrap();
                assert!(
                    result.status.success(),
                    "{}",
                    String::from_utf8_lossy(&result.stderr)
                );
                let actual = String::from_utf8(result.stdout).unwrap().replace(
                    &format!("<BlastOutput_db>{}</BlastOutput_db>", db.display()),
                    "<BlastOutput_db>@DATABASE@</BlastOutput_db>",
                );
                assert!(
                    actual.as_bytes() == expected,
                    "{mode} {format} p{threads}: bytes differ"
                );
            }
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn deprecated_workflows_use_native_validation_and_status() {
    for command in ["roc", "prepdb"] {
        for prefix in ["", "--"] {
            let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
                .arg(format!("{prefix}{command}"))
                .output()
                .unwrap();
            assert_eq!(
                result.status.code(),
                Some(if command == "roc" { 1 } else { 0 })
            );
            assert_eq!(
                source_startup::stdout_body_before_suffix(&result),
                b"Options:\n\n"
            );
            if command == "roc" {
                assert!(String::from_utf8_lossy(&result.stderr)
                    .ends_with("Error: Deprecated command: roc\n"));
            }
        }
        for (args, message, extra_accepts) in [
            (
                vec!["--db", "missing"],
                "Option is not permitted for this workflow: db",
                true,
            ),
            (
                vec!["-qmissing"],
                "Option is not permitted for this workflow: q",
                true,
            ),
            (
                vec!["--family-map", "missing"],
                "Option is not permitted for this workflow: family-map",
                true,
            ),
            (vec!["--output-hits"], "Invalid option: output-hits", true),
            (
                vec!["--definitely-unknown"],
                "Invalid option: definitely-unknown",
                false,
            ),
        ] {
            let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
                .arg(command)
                .args(args)
                .output()
                .unwrap();
            if cfg!(feature = "extra") && extra_accepts {
                // EXTRA registers these options and disables workflow membership
                // checks. Their values do not change either deprecated command.
                assert_eq!(
                    result.status.code(),
                    Some(if command == "roc" { 1 } else { 0 })
                );
                assert_eq!(source_startup::stdout_body_before_suffix(&result), b"");
                if command == "roc" {
                    assert!(String::from_utf8_lossy(&result.stderr)
                        .ends_with("Error: Deprecated command: roc\n"));
                }
            } else {
                assert_eq!(result.status.code(), Some(1));
                assert!(String::from_utf8_lossy(&result.stderr)
                    .ends_with(&format!("Error: {message}\n")));
            }
        }
    }
}

#[test]
fn translated_coverage_does_not_length_sort_equal_score_targets() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/search_output");
    let root = std::env::temp_dir().join(format!("diamond-output-ties-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    let db = root.join("tie.dmnd");
    let result = Command::new(binary)
        .args(["makedb", "--in"])
        .arg(fixtures.join("tie-reference.faa"))
        .arg("-d")
        .arg(&db)
        .arg("-p1")
        .output()
        .unwrap();
    assert!(result.status.success());
    for format in ["xml", "paf"] {
        let expected = std::fs::read(fixtures.join(format!("tie.{format}.out"))).unwrap();
        for (threads, memory) in [("1", "16G"), ("4", "0G")] {
            let result = Command::new(binary)
                .args(["blastx", "-q"])
                .arg(fixtures.join("tie-query.fna"))
                .arg("-d")
                .arg(&db)
                .args([
                    "-p",
                    threads,
                    "--memory-limit",
                    memory,
                    "--masking",
                    "0",
                    "--comp-based-stats",
                    "0",
                    "--query-cover",
                    "50",
                    "--subject-cover",
                    "50",
                    "--strand",
                    "plus",
                    "-f",
                    format,
                ])
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            let actual = String::from_utf8(result.stdout).unwrap().replace(
                &format!("<BlastOutput_db>{}</BlastOutput_db>", db.display()),
                "<BlastOutput_db>@DATABASE@</BlastOutput_db>",
            );
            assert!(
                actual.as_bytes() == expected,
                "{format} p{threads}: target order differs"
            );
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}
