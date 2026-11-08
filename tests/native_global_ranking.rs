//! Source-frozen two-phase ranking, final full extension and ordered output.
#[path = "support/source_startup.rs"]
mod source_startup;

use std::path::{Path, PathBuf};
use std::process::Command;

const FIELDS: &[&str] = &[
    "qseqid", "sseqid", "bitscore", "score", "qframe", "qstart", "qend", "sstart", "send", "cigar",
    "hspnum",
];

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}
fn directory(kind: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "diamond-global-ranking-{kind}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&path).unwrap();
    path
}
fn database(directory: &Path) -> PathBuf {
    let path = directory.join("reference.dmnd");
    let output = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args(["makedb", "--in"])
        .arg(root().join("tests/data/search_controls/reference.faa"))
        .arg("-d")
        .arg(&path)
        .arg("-p1")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    path
}

fn expected(mut bytes: Vec<u8>, format: &str, invocation: &str, database: &str) -> Vec<u8> {
    if ["edge", "null"].contains(&format) {
        return bytes;
    }
    let mut text = String::from_utf8(std::mem::take(&mut bytes)).unwrap();
    if format == "101" {
        let marker = "\tCL:";
        let begin = text.find(marker).unwrap() + marker.len();
        let end = begin + text[begin..].find('\n').unwrap();
        text.replace_range(begin..end, invocation);
    }
    if format == "5" {
        let marker = "<BlastOutput_db>";
        let begin = text.find(marker).unwrap() + marker.len();
        let end = begin + text[begin..].find("</BlastOutput_db>").unwrap();
        text.replace_range(begin..end, database);
    }
    text.into_bytes()
}

#[test]
fn global_ranking_matches_source_across_reference_blocks_and_final_formats() {
    let directory = directory("matrix");
    let db = database(&directory);
    let mut comparisons = 0;
    for row in std::fs::read_to_string(root().join("tests/data/global_ranking/cases.tsv"))
        .unwrap()
        .lines()
    {
        let c: Vec<_> = row.split('\t').collect();
        assert_eq!(c.len(), 8);
        for (threads, memory) in [("1", "16G"), ("4", "0G")] {
            let path = directory.join("output");
            let mut command = Command::new(env!("CARGO_BIN_EXE_diamond"));
            command
                .arg(c[1])
                .arg("-q")
                .arg(
                    root()
                        .join("tests/data/search_controls")
                        .join(if c[1] == "blastx" {
                            "query.fna"
                        } else {
                            "query.faa"
                        }),
                )
                .arg("-d")
                .arg(&db)
                .args([
                    "-p",
                    threads,
                    "--memory-limit",
                    memory,
                    "--comp-based-stats",
                    c[6],
                    "--masking",
                    "0",
                    "--soft-masking",
                    "none",
                    "--global-ranking",
                    c[2],
                    "--algo",
                    c[4],
                    "-c1",
                    "-b",
                    c[3],
                    "-f",
                    c[5],
                ]);
            if ["6", "104"].contains(&c[5]) {
                command.args(FIELDS);
            }
            command.args(c[7].split_whitespace()).arg("-o").arg(&path);
            let invocation = std::iter::once(command.get_program())
                .chain(command.get_args())
                .map(|value| value.to_string_lossy())
                .collect::<Vec<_>>()
                .join(" ");
            let output = command.output().unwrap();
            assert!(
                output.status.success(),
                "{} t{threads}: {}",
                c[0],
                String::from_utf8_lossy(&output.stderr)
            );
            let actual = std::fs::read(&path).unwrap();
            let source = expected(
                std::fs::read(root().join(format!("tests/data/global_ranking/{}.out", c[0])))
                    .unwrap(),
                c[5],
                &invocation,
                &db.to_string_lossy(),
            );
            let first = actual
                .iter()
                .zip(&source)
                .position(|(a, b)| a != b)
                .unwrap_or(actual.len().min(source.len()));
            assert!(
                actual == source,
                "{} t{threads}: first difference {first}, native/source lengths {}/{}",
                c[0],
                actual.len(),
                source.len()
            );
            comparisons += 1;
        }
    }
    assert_eq!(comparisons, 96);
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn global_ranking_option_errors_match_source() {
    let directory = directory("errors");
    let db = database(&directory);
    for (flags, diagnostic) in [
        (
            vec!["--global-ranking"],
            "Invalid parameter count for option '-g/--global-ranking'",
        ),
        (
            vec!["--global-ranking", "1", "2"],
            "Invalid parameter count for option '-g/--global-ranking'",
        ),
        (
            vec!["--global-ranking=2"],
            "Invalid option: global-ranking=2",
        ),
        (
            vec!["-g"],
            "Invalid parameter count for option '-g/--global-ranking'",
        ),
        (
            vec!["-g2", "--comp-based-stats", "2"],
            "Global ranking is not supported in this mode.",
        ),
        (
            vec!["-g2", "-f", "100"],
            "The DAA format is not supported in global ranking mode.",
        ),
        (
            vec!["-g2", "--frameshift", "15"],
            "Global ranking is not supported in this mode.",
        ),
        (
            vec!["-g2", "--range-culling"],
            "Global ranking is not supported in this mode.",
        ),
        (
            vec!["-g2", "--multiprocessing"],
            "Global ranking is not supported in this mode.",
        ),
        (
            vec!["--global-ranking", "-1"],
            "cannot create std::vector larger than max_size()",
        ),
        (
            vec!["--global-ranking", "9223372036854775807"],
            "cannot create std::vector larger than max_size()",
        ),
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("blastx")
            .arg("-q")
            .arg(root().join("tests/data/search_controls/query.fna"))
            .arg("-d")
            .arg(&db)
            .arg("-p1")
            .args(&flags)
            .arg("-o")
            .arg(directory.join("error"))
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(1), "{flags:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr).ends_with(
                &source_startup::source_console_text(&format!("Error: {diagnostic}\n"))
            ),
            "{flags:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    std::fs::remove_dir_all(directory).unwrap();
}
