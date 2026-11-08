#[path = "support/source_startup.rs"]
mod source_startup;

use std::path::Path;
use std::process::Command;

#[test]
fn reassignment_matches_the_explicitly_repaired_cpp_diagnostic() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/reassign");
    let root = std::env::temp_dir().join(format!("diamond-native-reassign-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    let variants: &[(&str, &[&str])] = &[
        ("default", &[]),
        ("mutual80", &["--mutual-cover", "80"]),
        ("member90", &["--member-cover", "90"]),
        ("approx90", &["--approx-id", "90"]),
        ("cbs0", &["--comp-based-stats", "0"]),
        ("cbs4", &["--comp-based-stats", "4"]),
        ("seg", &["--masking", "seg"]),
        ("soft-seg", &["--soft-masking", "seg"]),
        ("pam250", &["--matrix", "PAM250"]),
        ("full", &["--ext", "full"]),
        ("oid", &["--oid-output"]),
    ];
    assert_ne!(
        std::fs::read(fixtures.join("reference.default.diagnostic.tsv")).unwrap(),
        std::fs::read(fixtures.join("reference.mutual80.diagnostic.tsv")).unwrap(),
        "the fixture must exercise coverage-dependent reassignment"
    );
    for case in ["reference", "natural"] {
        let fasta = fixtures.join(format!("{case}.faa"));
        let db = root.join(format!("{case}.dmnd"));
        let result = Command::new(binary)
            .args(["makedb", "--in"])
            .arg(&fasta)
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
        for (variant, options) in variants {
            if case == "natural" && !["default", "mutual80", "approx90"].contains(variant) {
                continue;
            }
            let expected = std::fs::read(fixtures.join(if *variant == "soft-seg" {
                "reference.soft-seg.normalized-diagnostic.tsv".to_owned()
            } else {
                format!("{case}.{variant}.diagnostic.tsv")
            }))
            .unwrap();
            for database in [&fasta, &db] {
                for (memory, threads) in [("16G", "1"), ("0G", "4")] {
                    let output = root.join("output.tsv");
                    let result = Command::new(binary)
                        .args(["reassign", "-d"])
                        .arg(database)
                        .arg("--clusters")
                        .arg(fixtures.join(if case == "reference" {
                            "clusters.tsv"
                        } else {
                            "natural.clusters.tsv"
                        }))
                        .args(["--header", "-M", memory, "-p", threads])
                        .args(*options)
                        .arg("-o")
                        .arg(&output)
                        .output()
                        .unwrap();
                    let label = format!("{case} {variant} {database:?} {memory} p{threads}");
                    assert!(
                        result.status.success(),
                        "{label}: {}",
                        String::from_utf8_lossy(&result.stderr)
                    );
                    assert!(
                        std::fs::read(output).unwrap() == expected,
                        "{label}: bytes differ"
                    );
                }
            }
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn singleton_reassignment_matches_unmodified_upstream() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/reassign");
    let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args(["reassign", "-d"])
        .arg(fixtures.join("reference.faa"))
        .arg("--clusters")
        .arg(fixtures.join("singletons.tsv"))
        .args(["--header", "-p4"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        source_startup::stdout_body(&result),
        std::fs::read(fixtures.join("singletons.upstream.tsv")).unwrap()
    );
}

#[test]
fn reassignment_validates_configuration_before_input_io() {
    let cases: &[(&[&str], &str)] = &[
        (&["-d", "missing"], "Missing parameter: --clusters/-"),
        (
            &[
                "-d",
                "missing",
                "--clusters",
                "missing",
                "--member-cover",
                "80",
                "--mutual-cover",
                "80",
            ],
            "--member-cover and --mutual-cover are mutually exclusive.",
        ),
        (
            &["--query", "missing"],
            if cfg!(feature = "extra") {
                "Missing parameter: --db/-d"
            } else {
                "Option is not permitted for this workflow: query"
            },
        ),
        (
            &["--diag-filter-id", "50"],
            if cfg!(feature = "extra") {
                "Missing parameter: --db/-d"
            } else {
                "Invalid option: diag-filter-id"
            },
        ),
        (
            &["--comp-based-stats", "6"],
            "Invalid value for --comp-based-stats. Permitted values: 0, 1, 2, 3, 4, 5.",
        ),
    ];
    for (args, message) in cases {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("reassign")
            .args(*args)
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert!(
            String::from_utf8_lossy(&result.stderr).ends_with(&format!("Error: {message}\n")),
            "{args:?}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
}
