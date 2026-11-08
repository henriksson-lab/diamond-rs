//! SAM search callbacks compared with frozen output from the unchanged oracle.
use std::path::Path;
use std::process::Command;

fn with_invocation(bytes: Vec<u8>, invocation: &str) -> Vec<u8> {
    let text = String::from_utf8(bytes).unwrap();
    let (before, after) = text.split_once("\tCL:").unwrap();
    let (_, rest) = after.split_once('\n').unwrap();
    format!("{before}\tCL:{invocation}\n{rest}").into_bytes()
}

#[test]
fn native_sam_search_matches_source_across_blocks_algorithms_and_options() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/sam_search");
    let inputs = root.join("tests/data/tabular_context");
    let directory = std::env::temp_dir().join(format!("diamond-sam-search-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let database = directory.join("reference.dmnd");
    let binary = env!("CARGO_BIN_EXE_diamond");
    let result = Command::new(binary)
        .args(["makedb", "--in"])
        .arg(inputs.join("reference.faa"))
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
        let columns: Vec<_> = row.split('\t').collect();
        assert_eq!(columns.len(), 10);
        let expected = std::fs::read(fixtures.join(format!("{}.sam", columns[0]))).unwrap();
        for (threads, memory) in [("1", "16G"), ("4", "0G")] {
            let output = directory.join("output.sam");
            let mut command = Command::new(binary);
            command
                .arg(columns[1])
                .arg("-q")
                .arg(inputs.join(columns[2]))
                .arg("-d")
                .arg(&database)
                .args([
                    "--algo",
                    columns[3],
                    "-b",
                    columns[4],
                    "--masking",
                    columns[7],
                    "--comp-based-stats",
                    columns[8],
                    "--max-hsps",
                    "0",
                    "-p",
                    threads,
                    "--memory-limit",
                    memory,
                    "-f",
                    columns[9],
                ]);
            if columns[5] != "default" {
                command.args(["--unal", columns[5]]);
            }
            if columns[6] == "both" {
                command.args(["--salltitles", "--sam-query-len"]);
            }
            command.arg("-o").arg(&output);
            let invocation = std::iter::once(command.get_program())
                .chain(command.get_args())
                .map(|s| s.to_string_lossy())
                .collect::<Vec<_>>()
                .join(" ");
            let result = command.output().unwrap();
            assert!(
                result.status.success(),
                "{}/p{threads}: {}",
                columns[0],
                String::from_utf8_lossy(&result.stderr)
            );
            let actual = std::fs::read(&output).unwrap();
            let expected = with_invocation(expected.clone(), &invocation);
            let offset = actual
                .iter()
                .zip(&expected)
                .position(|(a, b)| a != b)
                .unwrap_or(actual.len().min(expected.len()));
            assert!(
                actual == expected,
                "{}/p{threads}: first differing byte {offset}; lengths {} and {}",
                columns[0],
                actual.len(),
                expected.len()
            );
            // The public in-process dispatcher must also select this native path.
            if threads == "1" && ["case-000", "case-041", "case-083"].contains(&columns[0]) {
                let args = std::iter::once(command.get_program())
                    .chain(command.get_args())
                    .map(|s| s.to_string_lossy().into_owned())
                    .collect::<Vec<_>>();
                assert_eq!(
                    diamond::run(&args.iter().map(String::as_str).collect::<Vec<_>>()),
                    0
                );
                assert_eq!(std::fs::read(&output).unwrap(), expected);
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 168);
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn sam_boolean_options_reject_values_like_source() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    for (flags, error) in [
        (
            vec!["--sam-query-len", "1"],
            "Invalid parameter count for option '--sam-query-len'",
        ),
        (vec!["--sam-query-len=1"], "Invalid option: sam-query-len=1"),
        (
            vec!["--salltitles", "1"],
            "Invalid parameter count for option '--salltitles'",
        ),
        (vec!["--salltitles=1"], "Invalid option: salltitles=1"),
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("blastp")
            .arg("-q")
            .arg(root.join("tests/data/tabular_context/query.faa"))
            .arg("-d")
            .arg(root.join("tests/data/tabular_context/reference.faa"))
            .args(["-f", "101"])
            .args(flags)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(String::from_utf8_lossy(&result.stderr).ends_with(&format!("Error: {error}\n")));
        assert!(!String::from_utf8_lossy(&result.stderr).contains("C++ FFI"));
    }
}
