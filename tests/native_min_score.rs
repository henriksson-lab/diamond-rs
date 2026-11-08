//! Frozen C++ cutoffs exercise score-only, coordinate/statistics and transcript
//! result paths, including merged reference/query blocks and adjusted matrices.
use std::path::Path;
use std::process::Command;

const FIELDS: &[&str] = &[
    "6", "qseqid", "sseqid", "score", "bitscore", "evalue", "qstart", "qend", "sstart", "send",
    "positive",
];

#[test]
fn bit_score_cutoffs_override_evalues_in_native_search() {
    let binary = env!("CARGO_BIN_EXE_diamond");
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let inputs = root.join("tests/data/multi_hsp");
    let fixtures = root.join("tests/data/min_score");
    let temporary = std::env::temp_dir().join(format!("diamond-min-score-{}", std::process::id()));
    std::fs::create_dir_all(&temporary).unwrap();
    for (command, extension, modes) in [("blastp", "faa", 6), ("blastx", "fna", 2)] {
        for cbs in 0..modes {
            for cutoff in ["0", "100", "170", "174.5", "175", "180"] {
                let label = format!("{command}.cbs{cbs}.score{cutoff}");
                let expected = std::fs::read(fixtures.join(format!("{label}.tsv"))).unwrap();
                if cutoff == "0" || cutoff == "180" {
                    assert!(expected.is_empty(), "{label}");
                }
                if cutoff == "100" {
                    assert!(!expected.is_empty(), "override must retain {label}");
                }
                for (block, memory, threads, format) in [
                    ("1", "16G", "1", "transcript"),
                    ("0.0000001", "0G", "4", "score"),
                    ("0.0000001", "16G", "1", "coordinates"),
                ] {
                    let output = temporary.join("output.tsv");
                    let mut process = Command::new(binary);
                    process
                        .arg(command)
                        .arg("-q")
                        .arg(inputs.join(format!("query.{extension}")))
                        .arg("-d")
                        .arg(inputs.join("reference.faa"))
                        .args([
                            "--masking",
                            "0",
                            "--comp-based-stats",
                            &cbs.to_string(),
                            "--max-hsps",
                            "0",
                            "--min-score",
                            cutoff,
                            "-e",
                            "1e-200",
                            "-b",
                            block,
                            "--memory-limit",
                            memory,
                            "-p",
                            threads,
                        ])
                        .arg("-o")
                        .arg(&output)
                        .arg("-f");
                    let fields = match format {
                        "score" => &FIELDS[..4],
                        "coordinates" => &FIELDS[..10],
                        _ => FIELDS,
                    };
                    let result = process.args(fields).output().unwrap();
                    assert!(
                        result.status.success(),
                        "{label} {format}: {}",
                        String::from_utf8_lossy(&result.stderr)
                    );
                    let expected = if format == "transcript" {
                        expected.clone()
                    } else {
                        String::from_utf8(expected.clone())
                            .unwrap()
                            .lines()
                            .map(|row| {
                                format!(
                                    "{}\n",
                                    row.split('\t')
                                        .take(fields.len() - 1)
                                        .collect::<Vec<_>>()
                                        .join("\t")
                                )
                            })
                            .collect::<String>()
                            .into_bytes()
                    };
                    assert_eq!(std::fs::read(output).unwrap(), expected, "{label} {format}");
                }
                if cbs == 1 && matches!(cutoff, "100" | "170") {
                    let output = temporary.join("output.daa");
                    let result = Command::new(binary)
                        .arg(command)
                        .arg("-q")
                        .arg(inputs.join(format!("query.{extension}")))
                        .arg("-d")
                        .arg(inputs.join("reference.faa"))
                        .args([
                            "--masking",
                            "0",
                            "--comp-based-stats",
                            "1",
                            "--max-hsps",
                            "0",
                            "--min-score",
                            cutoff,
                            "-e",
                            "1e-200",
                            "-p1",
                            "-f100",
                        ])
                        .arg("-o")
                        .arg(&output)
                        .output()
                        .unwrap();
                    assert!(
                        result.status.success(),
                        "{label}: {}",
                        String::from_utf8_lossy(&result.stderr)
                    );
                    assert_eq!(
                        std::fs::read(output).unwrap(),
                        std::fs::read(fixtures.join(format!("{label}.daa"))).unwrap(),
                        "{label} DAA"
                    );
                }
            }
        }
    }
    std::fs::remove_dir_all(temporary).unwrap();
}

#[test]
fn invalid_cutoff_arguments_do_not_require_cpp() {
    for command in ["blastp", "blastx"] {
        for (options, message) in [
            (
                vec!["--min-score"],
                "Invalid parameter count for option '--min-score'",
            ),
            (
                vec!["--min-score", "1", "2"],
                "Invalid parameter count for option '--min-score'",
            ),
            (vec!["--min-score=1"], "Invalid option: min-score=1"),
        ] {
            let output = Command::new(env!("CARGO_BIN_EXE_diamond"))
                .arg(command)
                .args(options)
                .output()
                .unwrap();
            assert_eq!(output.status.code(), Some(1));
            assert!(output.stdout.is_empty());
            assert_eq!(
                String::from_utf8(output.stderr).unwrap(),
                format!("Error: {message}\n")
            );
        }
    }
}

#[test]
fn cutoff_parser_preserves_source_numeric_conversion_and_last_value() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let inputs = root.join("tests/data/multi_hsp");
    let temporary =
        std::env::temp_dir().join(format!("diamond-min-score-parser-{}", std::process::id()));
    std::fs::create_dir_all(&temporary).unwrap();
    for (command, extension) in [("blastp", "faa"), ("blastx", "fna")] {
        for (options, expected_cutoff) in [
            (vec!["--min-score", "100tail"], "100"),
            (vec!["--min-score", "0x1.54p7"], "170"),
            (vec!["--min-score", "-3"], "100"),
            (vec!["--min-score", "junk"], "0"),
            (vec!["--min-score", "nan"], "0"),
            (vec!["--min-score", "inf"], "180"),
            (vec!["--min-score", "180", "--min-score", "100tail"], "100"),
        ] {
            let output = temporary.join("output.tsv");
            let mut args = vec![
                command.to_owned(),
                "-q".into(),
                inputs
                    .join(format!("query.{extension}"))
                    .to_string_lossy()
                    .into_owned(),
                "-d".into(),
                inputs.join("reference.faa").to_string_lossy().into_owned(),
            ];
            args.extend(
                [
                    "--masking",
                    "0",
                    "--comp-based-stats",
                    "1",
                    "--max-hsps",
                    "0",
                    "-e",
                    "1e-200",
                    "-p1",
                ]
                .into_iter()
                .chain(options)
                .map(str::to_owned),
            );
            args.extend([
                "-o".into(),
                output.to_string_lossy().into_owned(),
                "-f".into(),
            ]);
            args.extend(FIELDS.iter().map(|value| (*value).to_owned()));
            // Exercise the public in-process entry point too.
            let mut argv = vec!["diamond"];
            argv.extend(args.iter().map(String::as_str));
            assert_eq!(diamond::run(&argv), 0);
            let expected = std::fs::read(root.join(format!(
                "tests/data/min_score/{command}.cbs1.score{expected_cutoff}.tsv"
            )))
            .unwrap();
            assert_eq!(std::fs::read(output).unwrap(), expected, "{args:?}");
        }
    }
    std::fs::remove_dir_all(temporary).unwrap();
}
