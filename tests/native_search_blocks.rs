//! Real query/reference blocks compared with frozen, unmodified upstream output.
use std::path::{Path, PathBuf};
use std::process::Command;

fn compare(root: &Path, label: &str, mode: &str, format: &str, budget: &str, coverage: bool) {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let input = data.join("search_output");
    let reference = input.join("reference.faa");
    let query = input.join(if coverage {
        "bins-query.faa"
    } else if mode == "blastp" {
        "query.faa"
    } else {
        "query.fna"
    });
    let suffix = if format == "100" { "daa" } else { "out" };
    let expected =
        std::fs::read(data.join("search_blocks").join(format!("{label}.{suffix}"))).unwrap();
    for (threads, memory) in [("1", "16G"), ("4", "0G")] {
        let output = root.join(format!("{label}.{threads}.{suffix}"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_diamond"));
        command
            .arg(mode)
            .arg("-q")
            .arg(&query)
            .arg("-d")
            .arg(&reference)
            .args([
                "-b",
                budget,
                "-p",
                threads,
                "--memory-limit",
                memory,
                "--masking",
                "0",
                "-f",
                format,
                "-o",
            ])
            .arg(&output);
        if coverage {
            command.args(["--query-cover", "99", "--subject-cover", "99"]);
        } else {
            command.args(["--max-hsps", "0"]);
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{label}/p{threads}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let mut actual = std::fs::read(&output).unwrap();
        if format == "xml" {
            actual = String::from_utf8(actual)
                .unwrap()
                .replace(
                    &format!("<BlastOutput_db>{}</BlastOutput_db>", reference.display()),
                    "<BlastOutput_db>@DATABASE@</BlastOutput_db>",
                )
                .into_bytes();
        }
        assert!(
            actual == expected,
            "{label}/p{threads}/{memory}: bytes differ"
        );
    }
}

fn directory(name: &str) -> PathBuf {
    let root = std::env::temp_dir().join(format!("diamond-blocks-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    root
}

#[test]
fn native_blocks_preserve_all_five_output_formats_and_reverse_dna_coordinates() {
    let root = directory("formats");
    for mode in ["blastp", "blastx"] {
        for format in ["6", "0", "xml", "paf", "100"] {
            compare(
                &root,
                &format!("{mode}.{format}"),
                mode,
                format,
                "0.0000001",
                false,
            );
        }
        compare(
            &root,
            &format!("{mode}.direct.100"),
            mode,
            "100",
            "0.00002",
            false,
        );
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn query_block_sorting_retains_original_ids_and_mutual_coverage_seed_filtering() {
    let root = directory("coverage");
    for format in ["6", "xml", "paf"] {
        compare(
            &root,
            &format!("coverage.{format}"),
            "blastp",
            format,
            "0.0000006",
            true,
        );
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn paired_query_blocks_reject_unsynchronized_file_streams() {
    let root = directory("unequal-pair");
    let input = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/search_output");
    let query = input.join("query.fna");
    let records = std::fs::read_to_string(&query).unwrap();
    let first = records.split("\n>").next().unwrap();
    let mate = root.join("mate.fna");
    std::fs::write(&mate, format!("{first}\n")).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .arg("blastx")
        .arg("-q")
        .arg(query)
        .arg(mate)
        .arg("-d")
        .arg(input.join("reference.faa"))
        .args(["-b0.0000001", "-p1", "--masking", "0", "-o"])
        .arg(root.join("out.tsv"))
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr)
        .contains("Unequal number of sequences in paired read files."));
    std::fs::remove_dir_all(root).unwrap();
}
