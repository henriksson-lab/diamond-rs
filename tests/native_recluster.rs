#[path = "support/source_startup.rs"]
mod source_startup;

use std::path::Path;
use std::process::Command;

#[test]
fn recursive_clustering_matches_the_explicitly_repaired_cpp_diagnostic() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let fixtures = data.join("recluster");
    let root =
        std::env::temp_dir().join(format!("diamond-native-recluster-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    let variants: &[(&str, &[&str])] = &[
        ("default", &[]),
        ("approx90", &["--approx-id", "90"]),
        ("mutual80", &["--mutual-cover", "80"]),
        ("cbs4", &["--comp-based-stats", "4"]),
        ("one-round", &["--cluster-steps", "default"]),
        ("full", &["--ext", "full"]),
        ("seg", &["--masking", "seg"]),
        ("pam250", &["--matrix", "PAM250"]),
        ("oid", &["--oid-output"]),
    ];
    assert_ne!(
        std::fs::read(fixtures.join("reference.default.diagnostic.tsv")).unwrap(),
        std::fs::read(fixtures.join("reference.approx90.diagnostic.tsv")).unwrap(),
        "the recursive fixture must distinguish identity thresholds"
    );
    // Recluster's empty-but-present iterate setting chooses earlier hits.
    // A direct final-sensitivity search gives different natural assignments.
    assert_ne!(
        std::fs::read(fixtures.join("natural.default.diagnostic.tsv")).unwrap(),
        std::fs::read(data.join("reassign/natural.default.diagnostic.tsv")).unwrap()
    );
    for case in ["reference", "natural"] {
        let fasta = if case == "reference" {
            fixtures.join("reference.faa")
        } else {
            data.join("reassign/natural.faa")
        };
        let clusters = if case == "reference" {
            fixtures.join("clusters.tsv")
        } else {
            data.join("reassign/natural.clusters.tsv")
        };
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
            if case == "natural" && !["default", "approx90", "mutual80"].contains(variant) {
                continue;
            }
            let expected =
                std::fs::read(fixtures.join(format!("{case}.{variant}.diagnostic.tsv"))).unwrap();
            for database in [&fasta, &db] {
                for (threads, memory) in [("1", "16G"), ("4", "32K")] {
                    // PAM250's debug-build optimizer is expensive; the source
                    // conformance sweep covers all four storage/thread pairs.
                    if *variant == "pam250" && (database != &fasta || threads != "1") {
                        continue;
                    }
                    let output = root.join("output.tsv");
                    let result = Command::new(binary)
                        .args(["recluster", "-d"])
                        .arg(database)
                        .arg("--clusters")
                        .arg(&clusters)
                        .args(["--header", "-p", threads, "-M", memory])
                        .args(*options)
                        .arg("-o")
                        .arg(&output)
                        .output()
                        .unwrap();
                    let label = format!("{case} {variant} {database:?} p{threads} {memory}");
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
fn already_valid_clusters_match_unmodified_upstream() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    for extra in [
        vec![],
        vec!["--comp-based-stats", "0"],
        vec!["--approx-id", "0"],
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .args(["recluster", "-d"])
            .arg(data.join("reassign/reference.faa"))
            .arg("--clusters")
            .arg(data.join("recluster/valid.clusters.tsv"))
            .args(["--header", "-p4", "-M", "32K"])
            .args(extra)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            source_startup::stdout_body(&result),
            source_startup::source_console_text(
                &std::fs::read_to_string(data.join("recluster/valid.upstream.tsv")).unwrap()
            )
            .as_bytes()
        );
    }
}

#[test]
fn reclustering_validates_options_before_input_io() {
    let cases: &[(&[&str], &str)] = &[
        (&["-d", "missing"], "Missing parameter: --clusters/-"),
        (
            &["--query", "missing"],
            if cfg!(feature = "extra") {
                "Missing parameter: --db/-d"
            } else {
                "Option is not permitted for this workflow: query"
            },
        ),
        (
            &["--file-buffer-size", "1"],
            if cfg!(feature = "extra") {
                "Missing parameter: --db/-d"
            } else {
                "Option is not permitted for this workflow: file-buffer-size"
            },
        ),
        (
            &["--comp-based-stats", "6"],
            "Invalid value for --comp-based-stats. Permitted values: 0, 1, 2, 3, 4, 5.",
        ),
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
    ];
    for (args, message) in cases {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("recluster")
            .args(*args)
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert!(
            String::from_utf8_lossy(&result.stderr).ends_with(
                &source_startup::source_console_text(&format!("Error: {message}\n"))
            ),
            "{args:?}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
}
