//! Real taxonomy classification contexts against the unchanged source executable.
use std::path::Path;
use std::process::Command;

fn databases(root: &Path, directory: &Path) {
    let source = root.join("tests/data/search_blocks/taxonomy");
    let fixtures = root.join("tests/data/taxon_search");
    for variant in [
        "full", "none", "map", "no_names", "leaf", "unmapped", "multi", "root", "weak",
    ] {
        let reference = if ["leaf", "unmapped", "multi", "weak"].contains(&variant) {
            fixtures.join(format!("{variant}.faa"))
        } else {
            source.join("reference.faa")
        };
        let mut command = Command::new(env!("CARGO_BIN_EXE_diamond"));
        command
            .args(["makedb", "--in"])
            .arg(reference)
            .arg("-d")
            .arg(directory.join(format!("{variant}.dmnd")))
            .arg("-p1");
        if variant != "none" {
            command.arg("--taxonmap").arg(if variant == "weak" {
                fixtures.join("weak-mapping.tsv")
            } else if variant == "root" {
                fixtures.join("root-mapping.tsv")
            } else {
                source.join("mapping.tsv")
            });
        }
        if !["none", "map"].contains(&variant) {
            command
                .arg("--taxonnodes")
                .arg(if ["root", "weak"].contains(&variant) {
                    fixtures.join("root-nodes.dmp")
                } else {
                    source.join("nodes.dmp")
                });
        }
        if !["none", "map", "no_names"].contains(&variant) {
            command
                .arg("--taxonnames")
                .arg(if ["root", "weak"].contains(&variant) {
                    fixtures.join("root-names.dmp")
                } else {
                    source.join("names.dmp")
                });
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{variant}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
}
fn assert_wire(actual: Vec<u8>, expected: &[u8], label: &str) {
    let offset = actual
        .iter()
        .zip(expected)
        .position(|(a, b)| a != b)
        .unwrap_or(actual.len().min(expected.len()));
    assert!(
        actual == expected,
        "{label}: first differing byte {offset}; native: {} source: {}",
        String::from_utf8_lossy(&actual),
        String::from_utf8_lossy(expected)
    );
}
#[test]
fn native_taxonomy_uses_real_mappings_nodes_and_optional_names() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/taxon_search");
    let inputs = root.join("tests/data/tabular_context");
    let directory =
        std::env::temp_dir().join(format!("diamond-taxon-matrix-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    databases(root, &directory);
    let mut checked = 0;
    // The safe coverage reference is an explicitly documented source-UB exception.
    for (manifest, suffix) in [
        ("cases.tsv", ".tsv"),
        ("safe-coverage-cases.tsv", ".safe.tsv"),
    ] {
        for row in std::fs::read_to_string(fixtures.join(manifest))
            .unwrap()
            .lines()
        {
            let c: Vec<_> = row.split('\t').collect();
            assert_eq!(c.len(), 10);
            let expected = std::fs::read(fixtures.join(format!("{}{suffix}", c[0]))).unwrap();
            for (threads, memory) in [("1", "16G"), ("4", "0G")] {
                let output = directory.join("output.tsv");
                let mut command = Command::new(env!("CARGO_BIN_EXE_diamond"));
                command
                    .arg(c[1])
                    .arg("-q")
                    .arg(inputs.join(c[3]))
                    .arg("-d")
                    .arg(directory.join(format!("{}.dmnd", c[2])))
                    .args([
                        "--masking",
                        "0",
                        "--comp-based-stats",
                        "1",
                        "--algo",
                        c[4],
                        "-b",
                        c[5],
                        "--max-hsps",
                        c[8],
                        "-p",
                        threads,
                        "--memory-limit",
                        memory,
                        "-f",
                        "102",
                    ]);
                if c[6] == "1" {
                    command.arg("--include-lineage");
                }
                if c[7] != "default" {
                    command.args(["--unal", c[7]]);
                }
                match c[9] {
                    "none" => (),
                    "best" => {
                        command.args(["-k", "1"]);
                    }
                    "top" => {
                        command.args(["--top", "0"]);
                    }
                    "score" => {
                        command.args(["--min-score", "180"]);
                    }
                    "evalue" => {
                        command.args(["-e", "1e-58"]);
                    }
                    "coverage" => {
                        command.args(["--query-cover", "50", "--subject-cover", "50"]);
                    }
                    _ => panic!("unknown fixture filter"),
                }
                command.arg("-o").arg(&output);
                let result = command.output().unwrap();
                assert!(
                    result.status.success(),
                    "{}: {}",
                    c[0],
                    String::from_utf8_lossy(&result.stderr)
                );
                assert_wire(std::fs::read(&output).unwrap(), &expected, c[0]);
                if threads == "1" && ["case-000", "case-033", "case-086"].contains(&c[0]) {
                    let args = std::iter::once(command.get_program())
                        .chain(command.get_args())
                        .map(|s| s.to_string_lossy().into_owned())
                        .collect::<Vec<_>>();
                    assert_eq!(
                        diamond::run(&args.iter().map(String::as_str).collect::<Vec<_>>()),
                        0
                    );
                    assert_wire(std::fs::read(&output).unwrap(), &expected, c[0]);
                }
                checked += 1;
            }
        }
    }
    assert_eq!(checked, 256);
    std::fs::remove_dir_all(directory).unwrap();
}
#[test]
fn classifier_top_default_and_explicit_overrides_match_source() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/taxon_search");
    let directory =
        std::env::temp_dir().join(format!("diamond-taxon-defaults-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    databases(root, &directory);
    for (command, ext) in [("blastp", "faa"), ("blastx", "fna")] {
        for block in ["1", "0.0000001"] {
            for mode in ["default", "all", "score"] {
                let expected =
                    std::fs::read(fixtures.join(format!("weak.{command}.b{block}.{mode}.tsv")))
                        .unwrap();
                for (threads, memory) in [("1", "16G"), ("4", "0G")] {
                    let output = directory.join("output.tsv");
                    let mut process = Command::new(env!("CARGO_BIN_EXE_diamond"));
                    process
                        .arg(command)
                        .arg("-q")
                        .arg(root.join(format!("tests/data/tabular_context/query.{ext}")))
                        .arg("-d")
                        .arg(directory.join("weak.dmnd"))
                        .args([
                            "-f",
                            "102",
                            "--include-lineage",
                            "--masking",
                            "0",
                            "-b",
                            block,
                            "-p",
                            threads,
                            "--memory-limit",
                            memory,
                        ]);
                    if mode == "all" {
                        process.args(["--top", "100"]);
                    }
                    if mode == "score" {
                        process.args(["--min-score", "80"]);
                    }
                    let result = process.arg("-o").arg(&output).output().unwrap();
                    assert!(
                        result.status.success(),
                        "{}",
                        String::from_utf8_lossy(&result.stderr)
                    );
                    assert_wire(std::fs::read(output).unwrap(), &expected, mode);
                }
            }
        }
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn taxonomy_metadata_and_option_errors_match_source() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let directory =
        std::env::temp_dir().join(format!("diamond-taxon-errors-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    databases(root, &directory);
    for (variant,lineage,missing) in [
        ("none",false,"taxonomy mapping information (--taxonmap option), taxonomy nodes information (--taxonnodes option)"),
        ("none",true,"taxonomy mapping information (--taxonmap option), taxonomy nodes information (--taxonnodes option), taxonomy names information (--taxonnames option)"),
        ("map",false,"taxonomy nodes information (--taxonnodes option)"),
        ("map",true,"taxonomy nodes information (--taxonnodes option), taxonomy names information (--taxonnames option)"),
        ("no_names",true,"taxonomy names information (--taxonnames option)"),
    ] {
        let mut command=Command::new(env!("CARGO_BIN_EXE_diamond"));command.arg("blastp").arg("-q").arg(root.join("tests/data/tabular_context/query.faa")).arg("-d").arg(directory.join(format!("{variant}.dmnd"))).args(["-f","102"]);if lineage {command.arg("--include-lineage");}
        let result=command.output().unwrap();assert_eq!(result.status.code(),Some(1));
        let expected=format!("Error: Options require taxonomy information included in the database. Please use the respective options to build this information into the database when running diamond makedb: {missing}\n");assert!(String::from_utf8_lossy(&result.stderr).ends_with(&expected),"{}",String::from_utf8_lossy(&result.stderr));
    }
    for (options, error) in [
        (
            vec!["--include-lineage", "1"],
            "Invalid parameter count for option '--include-lineage'",
        ),
        (
            vec!["--include-lineage=1"],
            "Invalid option: include-lineage=1",
        ),
        (
            vec!["--top", "101"],
            "Allowed value range for --top is between 0.0 and 100.0",
        ),
        (
            vec!["--top", "-1"],
            "Allowed value range for --top is between 0.0 and 100.0",
        ),
        (
            vec!["--top", "10", "-k1"],
            "--top and -k/--max-target-seqs are mutually exclusive.",
        ),
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("blastp")
            .arg("-q")
            .arg(root.join("tests/data/tabular_context/query.faa"))
            .arg("-d")
            .arg(directory.join("full.dmnd"))
            .args(["-f", "102"])
            .args(options)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(
            String::from_utf8_lossy(&result.stderr).ends_with(&format!("Error: {error}\n")),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(!String::from_utf8_lossy(&result.stderr).contains("C++ FFI"));
    }
    std::fs::remove_dir_all(directory).unwrap();
}
