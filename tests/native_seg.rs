//! Upstream snapshots exercise SEG on targets while leaving queries unmasked.

use std::path::Path;
use std::process::Command;

#[test]
fn seg_search_matches_upstream_in_protein_and_translated_modes() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/seg");
    let root = std::env::temp_dir().join(format!("diamond-native-seg-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    let db = root.join("db");
    let result = Command::new(binary)
        .arg("makedb")
        .arg("--in")
        .arg(fixtures.join("reference.faa"))
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
    for (command, extension) in [("blastp", "faa"), ("blastx", "fna")] {
        for cbs in ["0", "1"] {
            for masking in ["seg", "0"] {
                let expected =
                    std::fs::read(fixtures.join(format!("{command}.{masking}.cbs{cbs}.tsv")))
                        .unwrap();
                assert!(!expected.is_empty());
                for memory in ["16G", "0G"] {
                    let output = root.join("output.tsv");
                    let result = Command::new(binary)
                        .arg(command)
                        .arg("-q")
                        .arg(fixtures.join(format!("query.{extension}")))
                        .arg("-d")
                        .arg(&db)
                        .args([
                            "-p1",
                            "--masking",
                            masking,
                            "--comp-based-stats",
                            cbs,
                            "--memory-limit",
                            memory,
                        ])
                        .arg("-o")
                        .arg(&output)
                        .args([
                            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch",
                            "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                            "score", "qlen", "slen", "nident", "positive", "gaps", "qframe",
                        ])
                        .output()
                        .unwrap();
                    assert!(
                        result.status.success(),
                        "{command} {masking} CBS {cbs} {memory}: {}",
                        String::from_utf8_lossy(&result.stderr)
                    );
                    assert_eq!(
                        std::fs::read(&output).unwrap(),
                        expected,
                        "{command} {masking} CBS {cbs} {memory}"
                    );
                }
                let archive = root.join("output.daa");
                let result = Command::new(binary)
                    .arg(command)
                    .arg("-q")
                    .arg(fixtures.join(format!("query.{extension}")))
                    .arg("-d")
                    .arg(&db)
                    .args([
                        "-p1",
                        "--masking",
                        masking,
                        "--comp-based-stats",
                        cbs,
                        "--outfmt",
                        "100",
                    ])
                    .arg("-o")
                    .arg(&archive)
                    .output()
                    .unwrap();
                assert!(
                    result.status.success(),
                    "{}",
                    String::from_utf8_lossy(&result.stderr)
                );
                assert_eq!(
                    std::fs::read(&archive).unwrap(),
                    std::fs::read(fixtures.join(format!("{command}.{masking}.cbs{cbs}.daa")))
                        .unwrap(),
                    "{command} {masking} CBS {cbs} DAA wire parity"
                );
            }
            // Both the mask-sensitive internal interval and the completely
            // masked target must have an observable effect on the output.
            assert_ne!(
                std::fs::read(fixtures.join(format!("{command}.seg.cbs{cbs}.tsv"))).unwrap(),
                std::fs::read(fixtures.join(format!("{command}.0.cbs{cbs}.tsv"))).unwrap()
            );
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}
