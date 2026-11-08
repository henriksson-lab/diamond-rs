#[path = "support/source_startup.rs"]
mod source_startup;

use std::path::Path;
use std::process::Command;

const STATS: &[&str] = &[
    "6",
    "qseqid",
    "sseqid",
    "approx_pident",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "score",
    "qlen",
    "slen",
    "nident",
    "positive",
    "gaps",
];

#[test]
fn native_realign_matches_upstream_blocks_modes_fields_and_lane_refills() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/realign");
    let root = std::env::temp_dir().join(format!("diamond-native-realign-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    for case in ["compact", "dynamic"] {
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
        for database in [&fasta, &db] {
            for cbs in 0..if case == "compact" { 6 } else { 2 } {
                for (memory, threads) in [
                    ("16G", "1"),
                    (if case == "compact" { "1K" } else { "16K" }, "4"),
                ] {
                    for format in ["default", "stats", "ids"] {
                        let output = root.join("output.tsv");
                        let mut process = Command::new(binary);
                        process
                            .args(["realign", "-d"])
                            .arg(database)
                            .arg("--clusters")
                            .arg(fixtures.join(format!("{case}.clusters.tsv")))
                            .args([
                                "--header",
                                "--comp-based-stats",
                                &cbs.to_string(),
                                "-M",
                                memory,
                                "-p",
                                threads,
                            ])
                            .arg("-o")
                            .arg(&output);
                        if format == "stats" {
                            process.arg("--outfmt").args(STATS);
                        }
                        if format == "ids" {
                            process.args(["--outfmt", "6", "qseqid", "sseqid", "score"]);
                        }
                        let result = process.output().unwrap();
                        let label = format!("{case} {database:?} cbs{cbs} {memory} {format}");
                        assert!(
                            result.status.success(),
                            "{label}: {}",
                            String::from_utf8_lossy(&result.stderr)
                        );
                        // Upstream realign uses Hauser correction only at mode 1;
                        // the other accepted modes do not adjust these matrices.
                        let golden_mode = if cbs == 1 { 1 } else { 0 };
                        let expected = std::fs::read(
                            fixtures.join(format!("{case}.cbs{golden_mode}.{memory}.{format}.tsv")),
                        )
                        .unwrap();
                        let actual = std::fs::read(&output).unwrap();
                        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
                        let avx2 = std::arch::is_x86_feature_detected!("avx2");
                        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
                        let avx2 = false;
                        if avx2 {
                            assert!(actual == expected, "{label}: bytes differ");
                        } else {
                            // Upstream's dynamic lane completion order depends
                            // on SIMD width. These snapshots use AVX2; other
                            // hosts still check every complete output row.
                            let mut actual = actual.split(|b| *b == b'\n').collect::<Vec<_>>();
                            let mut expected = expected.split(|b| *b == b'\n').collect::<Vec<_>>();
                            actual.sort_unstable();
                            expected.sort_unstable();
                            assert_eq!(actual, expected, "{label}: rows differ");
                        }
                    }
                }
            }
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn native_realign_rejects_invalid_workflows_and_output_before_input_io() {
    let cases: &[(&[&str], &str)] = &[
        (&["-d", "missing"], "Missing parameter: --clusters/-"),
        (
            &["--matrix", "PAM30"],
            if cfg!(feature = "extra") {
                "Missing parameter: --db/-d"
            } else {
                "Option is not permitted for this workflow: matrix"
            },
        ),
        (
            &["--cutoff-score-8bit", "100"],
            if cfg!(feature = "extra") {
                "Missing parameter: --db/-d"
            } else {
                "Invalid option: cutoff-score-8bit"
            },
        ),
        (
            &["--comp-based-stats", "6"],
            "Invalid value for --comp-based-stats. Permitted values: 0, 1, 2, 3, 4, 5.",
        ),
        (
            &["-d", "missing", "--clusters", "missing", "--outfmt", "100"],
            "The realign workflow only supports tabular output format.",
        ),
        (
            &[
                "-d",
                "missing",
                "--clusters",
                "missing",
                "--outfmt",
                "6",
                "qseq",
            ],
            "Unsupported output field for the realign workflow: qseq",
        ),
    ];
    for (args, message) in cases {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("realign")
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

#[test]
fn native_realign_extra_fields_match_the_intermediate_record_contract() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/realign");
    let root = std::env::temp_dir().join(format!(
        "diamond-native-realign-extra-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&root).unwrap();
    let clusters = root.join("clusters.tsv");
    let records = std::fs::read_to_string(fixtures.join("compact.clusters.tsv")).unwrap();
    std::fs::write(&clusters, records.split_once('\n').unwrap().1).unwrap();
    for cbs in ["0", "1"] {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .args(["realign", "-d"])
            .arg(fixtures.join("compact.faa"))
            .arg("--clusters")
            .arg(&clusters)
            .args([
                "--header",
                "0",
                "--comp-based-stats",
                cbs,
                "-p4",
                "--qnum-offset",
                "100",
                "--snum-offset",
                "-50",
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "qnum",
                "snum",
                "qtitle",
                "stitle",
                "sallseqid",
                "salltitles",
                "qcovhsp",
                "scovhsp",
                "corrected_bitscore",
                "normalized_nident",
                "hspnum",
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let actual = std::str::from_utf8(source_startup::stdout_body(&result)).unwrap();
        let mut defined_fields = String::new();
        for row in actual.lines() {
            let (fields, hspnum) = row.rsplit_once('\t').unwrap();
            // Upstream leaves hsp_num uninitialized after deserialization.
            // Rust deliberately retains its deterministic default of zero.
            assert_eq!(hspnum, "0");
            defined_fields.push_str(fields);
            defined_fields.push('\n');
        }
        let expected =
            std::fs::read_to_string(fixtures.join(format!("compact.extras.cbs{cbs}.tsv"))).unwrap();
        assert_eq!(defined_fields, expected);
    }
    std::fs::remove_dir_all(root).unwrap();
}
