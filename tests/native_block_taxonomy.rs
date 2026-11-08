//! Block-local target IDs must retain real database taxonomy after merging.
use std::path::Path;
use std::process::Command;

#[test]
fn blocked_taxonomy_matches_upstream_for_mapped_unmapped_and_multi_accession_targets() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixture = root.join("tests/data/search_blocks/taxonomy");
    let data = root.join("tests/data/search_output");
    let directory = std::env::temp_dir().join(format!("diamond-tax-blocks-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let database = directory.join("database.dmnd");
    let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args(["makedb", "--in"])
        .arg(fixture.join("reference.faa"))
        .arg("-d")
        .arg(&database)
        .arg("--taxonmap")
        .arg(fixture.join("mapping.tsv"))
        .arg("--taxonnodes")
        .arg(fixture.join("nodes.dmp"))
        .arg("--taxonnames")
        .arg(fixture.join("names.dmp"))
        .arg("-p1")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    for mode in ["blastp", "blastx"] {
        for budget in ["0.0000001", "0.00002"] {
            let expected = std::fs::read(fixture.join(format!("{mode}.b{budget}.tsv"))).unwrap();
            for (threads, memory) in [("1", "16G"), ("4", "0G")] {
                let output = directory.join("output.tsv");
                let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
                    .arg(mode)
                    .arg("-q")
                    .arg(data.join(if mode == "blastp" {
                        "query.faa"
                    } else {
                        "query.fna"
                    }))
                    .arg("-d")
                    .arg(&database)
                    .args([
                        "-p",
                        threads,
                        "--memory-limit",
                        memory,
                        "-b",
                        budget,
                        "--masking",
                        "0",
                        "--max-hsps",
                        "0",
                        "--algo",
                        "0",
                        "-f",
                        "6",
                        "qseqid",
                        "sseqid",
                        "qstart",
                        "qend",
                        "evalue",
                        "bitscore",
                        "staxids",
                        "sscinames",
                        "sskingdoms",
                        "skingdoms",
                        "sphylums",
                        "slineages",
                        "-o",
                    ])
                    .arg(&output)
                    .output()
                    .unwrap();
                assert!(
                    result.status.success(),
                    "{mode}/b{budget}/p{threads}: {}",
                    String::from_utf8_lossy(&result.stderr)
                );
                assert_eq!(
                    std::fs::read(output).unwrap(),
                    expected,
                    "{mode}/b{budget}/p{threads}"
                );
            }
        }
    }
    std::fs::remove_dir_all(directory).unwrap();
}
