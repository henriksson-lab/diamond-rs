//! Exercise the public runner in one process, without a C++ backend.
use std::path::Path;
use std::process::Command;

#[test]
fn native_library_errors_return_to_the_caller() {
    const CHILD: &str = "DIAMOND_NATIVE_RUN_RETURN_CHECK";
    if std::env::var_os(CHILD).is_some() {
        for args in [
            vec!["diamond", "blastp", "-p1", "--masking", "invalid"],
            vec!["diamond", "blastp", "-p1", "--memory-limit", "invalid"],
            vec!["diamond", "blastp", "-p1", "--target-indexed", "-c2"],
            vec!["diamond", "blastp", "-p1", "--fast", "--sensitive"],
            vec!["diamond", "blastp", "--max-hsps"],
            vec!["diamond", "getseq", "--reverse"],
            // EXTRA removes workflow restrictions: this becomes a valid,
            // complete regression run rather than an error-return case.
            #[cfg(not(feature = "extra"))]
            vec!["diamond", "test", "--log"],
            vec!["diamond", "unknown-command"],
            vec!["diamond", "blastp", "--legacy", "-p1"],
        ] {
            assert_eq!(diamond::run(&args), 1, "{args:?}");
        }
        assert_eq!(diamond::run(&["diamond", "version"]), 0);
        println!("native-library-returned-after-errors");
        return;
    }
    let result = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "native_library_errors_return_to_the_caller",
            "--nocapture",
        ])
        .env(CHILD, "1")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        String::from_utf8_lossy(&result.stdout).contains("native-library-returned-after-errors")
    );
    assert!(!String::from_utf8_lossy(&result.stderr).contains("C++ FFI"));
}

fn run(args: Vec<String>) {
    let borrowed: Vec<&str> = args.iter().map(String::as_str).collect();
    assert_eq!(diamond::run(&borrowed), 0, "{args:?}");
}

#[test]
fn public_runner_builds_searches_indexes_and_views_without_ffi() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let data = root.join("tests/data/search_output");
    let directory = std::env::temp_dir().join(format!("diamond-native-run-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let database = directory.join("database.dmnd");
    run(vec![
        "diamond".into(),
        "makedb".into(),
        "--in".into(),
        data.join("reference.faa").to_string_lossy().into_owned(),
        "-d".into(),
        database.to_string_lossy().into_owned(),
        "-p1".into(),
    ]);
    run(vec![
        "diamond".into(),
        "dbinfo".into(),
        "-d".into(),
        database.to_string_lossy().into_owned(),
    ]);
    run(vec![
        "diamond".into(),
        "makeidx".into(),
        "-d".into(),
        database.to_string_lossy().into_owned(),
        "--fast".into(),
        "-p4".into(),
    ]);
    assert!(directory.join("database.dmnd.seed_idx").exists());
    for mode in ["blastp", "blastx"] {
        let expected =
            std::fs::read(root.join(format!("tests/data/search_blocks/{mode}.6.out"))).unwrap();
        for (threads, memory) in [("1", "16G"), ("4", "0G")] {
            let output = directory.join("output.tsv");
            run(vec![
                "diamond".into(),
                mode.into(),
                "-q".into(),
                data.join(if mode == "blastp" {
                    "query.faa"
                } else {
                    "query.fna"
                })
                .to_string_lossy()
                .into_owned(),
                "-d".into(),
                database.to_string_lossy().into_owned(),
                "-p".into(),
                threads.into(),
                "--memory-limit".into(),
                memory.into(),
                "--masking".into(),
                "0".into(),
                "--max-hsps".into(),
                "0".into(),
                "-o".into(),
                output.to_string_lossy().into_owned(),
            ]);
            assert_eq!(
                std::fs::read(output).unwrap(),
                expected,
                "{mode}/p{threads}"
            );
        }
    }
    let archive = directory.join("output.daa");
    run(vec![
        "diamond".into(),
        "blastp".into(),
        "-q".into(),
        data.join("query.faa").to_string_lossy().into_owned(),
        "-d".into(),
        database.to_string_lossy().into_owned(),
        "-f100".into(),
        "--masking".into(),
        "0".into(),
        "-p1".into(),
        "-o".into(),
        archive.to_string_lossy().into_owned(),
    ]);
    let viewed = directory.join("viewed.tsv");
    run(vec![
        "diamond".into(),
        "view".into(),
        "-a".into(),
        archive.to_string_lossy().into_owned(),
        "-o".into(),
        viewed.to_string_lossy().into_owned(),
        "-f6".into(),
        "-p1".into(),
    ]);
    assert!(!std::fs::read(viewed).unwrap().is_empty());
    std::fs::remove_dir_all(directory).unwrap();
}
