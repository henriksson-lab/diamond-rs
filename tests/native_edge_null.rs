//! Source binary edge and null output through the actual native search pipeline.
#[path = "support/source_startup.rs"]
mod source_startup;

use std::path::Path;
use std::process::Command;

#[test]
fn native_edge_and_null_match_source_across_blocks_algorithms_and_filters() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/edge_null");
    let inputs = root.join("tests/data/tabular_context");
    let directory = std::env::temp_dir().join(format!("diamond-edge-null-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let database = directory.join("reference.dmnd");
    let binary = env!("CARGO_BIN_EXE_diamond");
    let result = Command::new(binary)
        .args(["makedb", "--in"])
        .arg(inputs.join("reference.faa"))
        .arg("-d")
        .arg(&database)
        .arg("-p1")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let cases = std::fs::read_to_string(fixtures.join("cases.tsv")).unwrap();
    let mut checked = 0;
    for row in cases.lines() {
        let c: Vec<_> = row.split('\t').collect();
        assert_eq!(c.len(), 10);
        let expected = std::fs::read(fixtures.join(format!("{}.out", c[0]))).unwrap();
        for (threads, memory) in [("1", "16G"), ("4", "0G")] {
            let output = directory.join("output");
            let mut command = Command::new(binary);
            command
                .arg(c[1])
                .arg("-q")
                .arg(inputs.join(c[3]))
                .arg("-d")
                .arg(&database)
                .args([
                    "-f",
                    c[2],
                    "--algo",
                    c[4],
                    "-b",
                    c[5],
                    "--masking",
                    c[7],
                    "--comp-based-stats",
                    c[8],
                    "--max-hsps",
                    "0",
                    "-p",
                    threads,
                    "--memory-limit",
                    memory,
                ]);
            if c[6] != "default" {
                command.args(["--unal", c[6]]);
            }
            match c[9] {
                "none" => (),
                "identity" => {
                    command.args(["--id", "95"]);
                }
                "coverage" => {
                    command.args(["--query-cover", "25", "--subject-cover", "50"]);
                }
                "top" => {
                    command.args(["--top", "10"]);
                }
                "minscore" => {
                    command.args(["--min-score", "170"]);
                }
                _ => panic!("unknown fixture filter"),
            }
            command.arg("-o").arg(&output);
            let result = command.output().unwrap();
            assert!(
                result.status.success(),
                "{}/p{threads}: {}",
                c[0],
                String::from_utf8_lossy(&result.stderr)
            );
            let actual = std::fs::read(&output).unwrap();
            let offset = actual
                .iter()
                .zip(&expected)
                .position(|(a, b)| a != b)
                .unwrap_or(actual.len().min(expected.len()));
            assert!(
                actual == expected,
                "{}/p{threads}: first differing byte {offset}; lengths {} and {}",
                c[0],
                actual.len(),
                expected.len()
            );
            if c[2] == "null" {
                assert!(
                    String::from_utf8_lossy(&result.stderr).contains("Seed matches:"),
                    "null must execute search"
                );
            }
            if threads == "1" && ["case-000", "case-012", "case-083"].contains(&c[0]) {
                let args = std::iter::once(command.get_program())
                    .chain(command.get_args())
                    .map(|s| s.to_string_lossy().into_owned())
                    .collect::<Vec<_>>();
                assert_eq!(
                    diamond::run(&args.iter().map(String::as_str).collect::<Vec<_>>()),
                    0
                );
                assert_eq!(std::fs::read(&output).unwrap(), expected);
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 168);
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn corrected_score_binary_precision_matches_source_for_varied_scores_and_frames() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/edge_null");
    let directory =
        std::env::temp_dir().join(format!("diamond-edge-precision-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let database = directory.join("reference.dmnd");
    let binary = env!("CARGO_BIN_EXE_diamond");
    let result = Command::new(binary)
        .args(["makedb", "--in"])
        .arg(root.join("tests/data/tabular_context/reference.faa"))
        .arg("-d")
        .arg(&database)
        .arg("-p1")
        .output()
        .unwrap();
    assert!(result.status.success());
    for (command, ext) in [("blastp", "faa"), ("blastx", "fna")] {
        for algorithm in ["0", "1", "mutual"] {
            for block in ["1", "0.0000001"] {
                let label = if algorithm == "mutual" {
                    format!("precision.{command}.mutual.b{block}")
                } else {
                    format!("precision.{command}.a{algorithm}.b{block}")
                };
                let expected = std::fs::read(fixtures.join(format!("{label}.out"))).unwrap();
                assert_eq!(
                    expected.len(),
                    (if algorithm == "mutual" { 120 } else { 180 }) * 32
                );
                for (threads, memory) in [("1", "16G"), ("4", "0G")] {
                    let output = directory.join("output");
                    let result = Command::new(binary)
                        .arg(command)
                        .arg("-q")
                        .arg(fixtures.join(format!("precision.{ext}")))
                        .arg("-d")
                        .arg(&database)
                        .args([
                            "-f",
                            "edge",
                            "--algo",
                            if algorithm == "mutual" {
                                "0"
                            } else {
                                algorithm
                            },
                            "-b",
                            block,
                            "--masking",
                            "0",
                            "--max-hsps",
                            "0",
                            "-p",
                            threads,
                            "--memory-limit",
                            memory,
                        ])
                        .args(if algorithm == "mutual" {
                            vec!["--query-cover", "50", "--subject-cover", "50"]
                        } else {
                            Vec::new()
                        })
                        .arg("-o")
                        .arg(&output)
                        .output()
                        .unwrap();
                    assert!(
                        result.status.success(),
                        "{label}: {}",
                        String::from_utf8_lossy(&result.stderr)
                    );
                    let actual = std::fs::read(output).unwrap();
                    let offset = actual
                        .iter()
                        .zip(&expected)
                        .position(|(a, b)| a != b)
                        .unwrap_or(actual.len().min(expected.len()));
                    assert!(
                        actual == expected,
                        "{label}/p{threads}: first differing byte {offset}"
                    );
                }
            }
        }
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn null_and_edge_still_validate_actual_query_records() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    for format in ["edge", "null"] {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("blastp")
            .arg("-q")
            .arg(root.join("tests/data/tabular_context/multiline-sequence.fq"))
            .arg("-d")
            .arg(root.join("tests/data/tabular_context/reference.faa"))
            .args(["-f", format])
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(String::from_utf8_lossy(&result.stderr).ends_with(
            &source_startup::source_console_text("Error: Malformed FASTQ record\n")
        ));
        assert!(!String::from_utf8_lossy(&result.stderr).contains("C++ FFI"));
    }
}
