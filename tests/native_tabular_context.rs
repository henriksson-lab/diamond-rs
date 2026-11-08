//! Concrete tabular callbacks against the immutable source's output files.
use std::path::Path;
use std::process::Command;

fn assert_wire(actual: Vec<u8>, expected: &[u8], label: &str) {
    let offset = actual
        .iter()
        .zip(expected)
        .position(|(a, b)| a != b)
        .unwrap_or(actual.len().min(expected.len()));
    assert!(
        actual == expected,
        "{label}: wire bytes first differ at {offset}; lengths {} and {}",
        actual.len(),
        expected.len()
    );
}

#[test]
fn tabular_fields_retain_real_context_across_native_block_merging() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/tabular_context");
    let temporary =
        std::env::temp_dir().join(format!("diamond-tabular-context-{}", std::process::id()));
    std::fs::create_dir_all(&temporary).unwrap();
    let database = temporary.join("reference.dmnd");
    let binary = env!("CARGO_BIN_EXE_diamond");
    let result = Command::new(binary)
        .arg("makedb")
        .arg("--in")
        .arg(fixtures.join("reference.faa"))
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
        let columns = row.split('\t').collect::<Vec<_>>();
        assert_eq!(columns.len(), 7);
        let expected = std::fs::read(fixtures.join(format!("{}.tsv", columns[0]))).unwrap();
        assert!(!expected.is_empty(), "{}", columns[0]);
        for (threads, memory) in [("1", "16G"), ("4", "0G")] {
            let output = temporary.join("output.tsv");
            let mut process = Command::new(binary);
            process.arg(columns[1]).arg("-q");
            for query in columns[2].split(',') {
                process.arg(fixtures.join(query));
            }
            process.arg("-d").arg(&database).args([
                "--masking",
                "0",
                "--comp-based-stats",
                columns[3],
                "--max-hsps",
                "0",
                "-b",
                columns[4],
                "-p",
                threads,
                "--memory-limit",
                memory,
            ]);
            if !columns[6].is_empty() {
                process.args(columns[6].split(','));
            }
            let result = process
                .arg("-o")
                .arg(&output)
                .args(["-f", "6"])
                .args(columns[5].split(','))
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{}/p{threads}: {}",
                columns[0],
                String::from_utf8_lossy(&result.stderr)
            );
            assert_eq!(
                std::fs::read(&output).unwrap(),
                expected,
                "{}/p{threads}",
                columns[0]
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 132);
    std::fs::remove_dir_all(temporary).unwrap();
}

#[test]
fn quality_and_field_validation_match_source_without_ffi() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/tabular_context");
    let temporary =
        std::env::temp_dir().join(format!("diamond-tabular-errors-{}", std::process::id()));
    std::fs::create_dir_all(&temporary).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    for name in ["short-quality", "multiline-sequence"] {
        let output = temporary.join("output.tsv");
        let result = Command::new(binary)
            .arg("blastp")
            .arg("-q")
            .arg(fixtures.join(format!("{name}.fq")))
            .arg("-d")
            .arg(fixtures.join("reference.faa"))
            .args([
                "--masking",
                "0",
                "--max-hsps",
                "0",
                "-p1",
                "-f",
                "6",
                "qseqid",
                "full_qqual",
            ])
            .arg("-o")
            .arg(&output)
            .output()
            .unwrap();
        assert_eq!(
            result.status.code(),
            Some(if name == "multiline-sequence" { 1 } else { 0 })
        );
        assert_eq!(
            std::fs::read(output).unwrap(),
            std::fs::read(fixtures.join(format!("{name}.tsv"))).unwrap()
        );
        if name == "multiline-sequence" {
            assert!(String::from_utf8_lossy(&result.stderr)
                .ends_with("Error: Malformed FASTQ record\n"));
        }
        assert!(!String::from_utf8_lossy(&result.stderr).contains("C++ FFI"));
    }
    for (fields, extra, expected) in [
        ("qacc", vec![], "Invalid output field: qacc"),
        (
            "full_qseq_mate",
            vec![],
            "Output field only supported for translated search.",
        ),
        (
            "qseq_translated",
            vec![],
            "Output field only supported for translated search.",
        ),
        (
            "approx_pident",
            vec!["--matrix", "PAM30"],
            "Approximate identity is only supported for the BLOSUM62 scoring matrix.",
        ),
    ] {
        let result = Command::new(binary)
            .arg("blastp")
            .arg("-q")
            .arg(fixtures.join("query.faa"))
            .arg("-d")
            .arg(fixtures.join("reference.faa"))
            .args(extra)
            .args(["-f", "6", fields])
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(String::from_utf8_lossy(&result.stderr).ends_with(&format!("Error: {expected}\n")));
        assert!(!String::from_utf8_lossy(&result.stderr).contains("C++ FFI"));
    }
    for field in [
        "hspnum",
        "normalized_bitscore",
        "normalized_nident",
        "corrected_bitscore",
    ] {
        let output = temporary.join("error-output.tsv");
        let result = Command::new(binary)
            .arg("blastp")
            .arg("-q")
            .arg(fixtures.join("query.faa"))
            .arg("-d")
            .arg(fixtures.join("reference.faa"))
            .args([
                "--masking",
                "0",
                "--max-hsps",
                "0",
                "--unal",
                "1",
                "-p1",
                "-b0.0000001",
                "-f",
                "6",
                "qseqid",
                "qtitle",
                field,
            ])
            .arg("-o")
            .arg(&output)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(String::from_utf8_lossy(&result.stderr)
            .ends_with(&format!("Error: Invalid output field: {field}\n")));
        assert_eq!(
            std::fs::read(output).unwrap(),
            std::fs::read(fixtures.join(format!("unaligned-error.{field}.tsv"))).unwrap()
        );
    }
    std::fs::remove_dir_all(temporary).unwrap();
}

#[test]
fn all_dynamic_taxonomic_rank_fields_use_database_metadata() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/tabular_context");
    let taxonomy = root.join("tests/data/search_blocks/taxonomy");
    let temporary =
        std::env::temp_dir().join(format!("diamond-tabular-ranks-{}", std::process::id()));
    std::fs::create_dir_all(&temporary).unwrap();
    let database = temporary.join("taxonomy.dmnd");
    let binary = env!("CARGO_BIN_EXE_diamond");
    let result = Command::new(binary)
        .args(["makedb", "--in"])
        .arg(taxonomy.join("reference.faa"))
        .arg("-d")
        .arg(&database)
        .arg("--taxonmap")
        .arg(taxonomy.join("mapping.tsv"))
        .arg("--taxonnodes")
        .arg(taxonomy.join("nodes.dmp"))
        .arg("--taxonnames")
        .arg(taxonomy.join("names.dmp"))
        .arg("-p1")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let ranks = diamond::data::taxonomy_nodes::Rank::NAMES
        .iter()
        .skip(1)
        .map(|name| format!("s{}", name.replace(' ', "_")))
        .collect::<Vec<_>>();
    assert_eq!(ranks.len(), 48);
    for (command, extension) in [("blastp", "faa"), ("blastx", "fna")] {
        for block in ["1", "0.0000001"] {
            let expected =
                std::fs::read(fixtures.join(format!("{command}.dynamic-ranks.b{block}.tsv")))
                    .unwrap();
            assert!(!expected.is_empty());
            for (threads, memory) in [("1", "16G"), ("4", "0G")] {
                let output = temporary.join("output.tsv");
                let result = Command::new(binary)
                    .arg(command)
                    .arg("-q")
                    .arg(fixtures.join(format!("query.{extension}")))
                    .arg("-d")
                    .arg(&database)
                    .args([
                        "--masking",
                        "0",
                        "-p",
                        threads,
                        "-b",
                        block,
                        "--memory-limit",
                        memory,
                        "-f",
                        "6",
                        "qseqid",
                        "sseqid",
                    ])
                    .args(&ranks)
                    .arg("-o")
                    .arg(&output)
                    .output()
                    .unwrap();
                assert!(
                    result.status.success(),
                    "{}",
                    String::from_utf8_lossy(&result.stderr)
                );
                assert_eq!(
                    std::fs::read(output).unwrap(),
                    expected,
                    "{command}/b{block}/p{threads}"
                );
            }
        }
    }
    std::fs::remove_dir_all(temporary).unwrap();
}

#[test]
fn flat_json_preserves_source_callbacks_and_separator_quirks() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/tabular_context");
    let taxonomy = root.join("tests/data/search_blocks/taxonomy");
    let temporary =
        std::env::temp_dir().join(format!("diamond-json-context-{}", std::process::id()));
    std::fs::create_dir_all(&temporary).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    for (name, input) in [
        ("plain", fixtures.join("reference.faa")),
        ("single", fixtures.join("single-target.faa")),
        ("taxonomy", taxonomy.join("reference.faa")),
    ] {
        let mut process = Command::new(binary);
        process
            .args(["makedb", "--in"])
            .arg(input)
            .arg("-d")
            .arg(temporary.join(format!("{name}.dmnd")))
            .arg("-p1");
        if name == "taxonomy" {
            process
                .arg("--taxonmap")
                .arg(taxonomy.join("mapping.tsv"))
                .arg("--taxonnodes")
                .arg(taxonomy.join("nodes.dmp"))
                .arg("--taxonnames")
                .arg(taxonomy.join("names.dmp"));
        }
        let result = process.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    let cases = std::fs::read_to_string(fixtures.join("json-cases.tsv")).unwrap();
    let mut checked = 0;
    for row in cases.lines() {
        let columns = row.split('\t').collect::<Vec<_>>();
        assert_eq!(columns.len(), 8);
        let expected = std::fs::read(fixtures.join(format!("{}.out", columns[0]))).unwrap();
        for (threads, memory) in [("1", "16G"), ("4", "0G")] {
            let output = temporary.join("output.json");
            let mut process = Command::new(binary);
            process.arg(columns[1]).arg("-q");
            for query in columns[2].split(',') {
                process.arg(fixtures.join(query));
            }
            let result = process
                .arg("-d")
                .arg(temporary.join(format!("{}.dmnd", columns[7])))
                .args([
                    "--masking",
                    "0",
                    "--comp-based-stats",
                    columns[3],
                    "--max-hsps",
                    "0",
                    "-b",
                    columns[4],
                    "-p",
                    threads,
                    "--memory-limit",
                    memory,
                    "-f",
                    columns[5],
                ])
                .args(columns[6].split(','))
                .arg("-o")
                .arg(&output)
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{}/p{threads}: {}",
                columns[0],
                String::from_utf8_lossy(&result.stderr)
            );
            assert_wire(
                std::fs::read(output).unwrap(),
                &expected,
                &format!("{}/p{threads}", columns[0]),
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 44);
    std::fs::remove_dir_all(temporary).unwrap();
}

#[test]
fn masked_search_reports_original_subject_sequences() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/tabular_context");
    let temporary =
        std::env::temp_dir().join(format!("diamond-tabular-mask-{}", std::process::id()));
    std::fs::create_dir_all(&temporary).unwrap();
    let database = temporary.join("reference.dmnd");
    let binary = env!("CARGO_BIN_EXE_diamond");
    let result = Command::new(binary)
        .args(["makedb", "--in"])
        .arg(fixtures.join("masked-reference.faa"))
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
    for (command, extension) in [("blastp", "faa"), ("blastx", "fna")] {
        for mask in ["1", "seg"] {
            for block in ["1", "0.0000001"] {
                let label = format!("{command}.mask{mask}.b{block}");
                for algorithm in ["0", "1"] {
                    let suffix = if algorithm == "1" { ".algo1" } else { "" };
                    let expected =
                        std::fs::read(fixtures.join(format!("{label}{suffix}.tsv"))).unwrap();
                    for (threads, memory) in [("1", "16G"), ("4", "0G")] {
                        let output = temporary.join("output.tsv");
                        let result = Command::new(binary)
                            .arg(command)
                            .arg("-q")
                            .arg(fixtures.join(format!("query.{extension}")))
                            .arg("-d")
                            .arg(&database)
                            .args([
                                "--masking",
                                mask,
                                "--max-hsps",
                                "0",
                                "--algo",
                                algorithm,
                                "-b",
                                block,
                                "-p",
                                threads,
                                "--memory-limit",
                                memory,
                                "-f",
                                "6",
                                "qseqid",
                                "sseqid",
                                "qseq",
                                "sseq",
                                "full_sseq",
                                "full_qseq",
                                "cigar",
                                "normalized_bitscore",
                                "approx_pident",
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
                        assert_wire(
                            std::fs::read(output).unwrap(),
                            &expected,
                            &format!("{label}/a{algorithm}/p{threads}"),
                        );
                    }
                }
            }
        }
    }
    std::fs::remove_dir_all(temporary).unwrap();
}
