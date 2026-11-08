//! Frozen upstream regression hashes against the real native search pipeline.
//! No C++ checkout or executable is read by these tests.
use diamond::commands::test_cmd::NativeRegressionBackend;
use diamond::upstream_test::{data, test, test_cases};
use test::RegressionBackend;

#[test]
fn native_search_matches_twenty_upstream_builtin_hashes() {
    let mut backend = NativeRegressionBackend::default();
    let mut query = backend.create_sequence_file("query").unwrap();
    let mut db = backend.create_sequence_file("db").unwrap();
    test::load_seqs(&mut query, data::seqs()).unwrap();
    test::load_seqs(&mut db, data::seqs()).unwrap();
    let mut passed = 0;
    for (index, case) in test_cases::TEST_CASES.iter().enumerate() {
        let mut output = Vec::new();
        let count = test::run_testcase(
            case,
            test_cases::REF_HASHES[index],
            &mut db,
            &mut query,
            0,
            test::TestRunConfig::default(),
            &mut backend,
            &mut output,
        )
        .unwrap();
        assert_eq!(
            count,
            1,
            "{}: {}",
            case.description,
            String::from_utf8_lossy(&output)
        );
        passed += count;
    }
    assert_eq!(passed, 20);
}

#[test]
fn regression_debug_logging_uses_connected_native_search() {
    const CHILD: &str = "DIAMOND_REGRESSION_LOG_CHILD";
    if std::env::var_os(CHILD).is_some() {
        let mut backend = NativeRegressionBackend::default();
        let mut query = backend.create_sequence_file("query").unwrap();
        let mut db = backend.create_sequence_file("db").unwrap();
        test::load_seqs(&mut query, &data::seqs()[..1]).unwrap();
        test::load_seqs(&mut db, &data::seqs()[..1]).unwrap();
        let args = ["diamond", "blastp", "--log", "-p1"].map(str::to_string);
        let output = backend.run_search(&args, &mut db, &mut query).unwrap();
        assert!(!output.is_empty(), "real self-search must produce rows");
        drop(backend);
        let log = std::fs::read_to_string("diamond.log").unwrap();
        assert!(log.contains("Scoring parameters:"), "{log}");
        assert!(log.contains("Current RSS:"), "{log}");
        return;
    }
    // Source's logger uses a fixed relative filename. Isolate the real
    // resource in a child process without changing the parallel suite's cwd.
    let root = std::env::temp_dir().join(format!("diamond-regression-log-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    let result = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "regression_debug_logging_uses_connected_native_search",
            "--nocapture",
        ])
        .env(CHILD, "1")
        .current_dir(&root)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "stdout: {}\nstderr: {}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn native_blocked_search_matches_upstream_builtin_hash() {
    let mut backend = NativeRegressionBackend::default();
    let mut query = backend.create_sequence_file("query").unwrap();
    let mut db = backend.create_sequence_file("db").unwrap();
    test::load_seqs(&mut query, data::seqs()).unwrap();
    test::load_seqs(&mut db, data::seqs()).unwrap();
    let mut output = Vec::new();
    assert_eq!(
        test::run_testcase(
            &test_cases::TEST_CASES[2],
            test_cases::REF_HASHES[2],
            &mut db,
            &mut query,
            0,
            test::TestRunConfig::default(),
            &mut backend,
            &mut output
        )
        .unwrap(),
        1,
        "{}",
        String::from_utf8_lossy(&output)
    );
}

#[test]
fn native_query_indexed_cli_matches_upstream_for_both_query_types_and_masks() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let input = root.join("tests/data/search_output");
    let expected = root.join("tests/data/query_indexed");
    let output =
        std::env::temp_dir().join(format!("diamond-query-indexed-{}.tsv", std::process::id()));
    for mode in ["blastp", "blastx"] {
        for mask in ["0", "tantan", "seg"] {
            let snapshot = std::fs::read(expected.join(format!("{mode}.{mask}.tsv"))).unwrap();
            for (chunks, threads, memory) in [("1", "1", "16G"), ("2", "4", "0G"), ("4", "4", "0G")]
            {
                let result = std::process::Command::new(env!("CARGO_BIN_EXE_diamond"))
                    .arg(mode)
                    .arg("-q")
                    .arg(input.join(if mode == "blastp" {
                        "query.faa"
                    } else {
                        "query.fna"
                    }))
                    .arg("-d")
                    .arg(input.join("reference.faa"))
                    .args([
                        "--algo",
                        "1",
                        "--masking",
                        mask,
                        "--index-chunks",
                        chunks,
                        "-p",
                        threads,
                        "--memory-limit",
                        memory,
                        "-o",
                    ])
                    .arg(&output)
                    .output()
                    .unwrap();
                assert!(
                    result.status.success(),
                    "{mode}/{mask}/c{chunks}: {}",
                    String::from_utf8_lossy(&result.stderr)
                );
                assert_eq!(
                    std::fs::read(&output).unwrap(),
                    snapshot,
                    "{mode}/{mask}/c{chunks}"
                );
            }
        }
    }
    std::fs::remove_file(output).unwrap();
}

#[test]
fn native_index_generation_honors_all_previously_ignored_sensitivity_flags() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let temp =
        std::env::temp_dir().join(format!("diamond-index-sensitivity-{}", std::process::id()));
    std::fs::create_dir_all(&temp).unwrap();
    let binary = env!("CARGO_BIN_EXE_diamond");
    let database = temp.join("reference.dmnd");
    let result = std::process::Command::new(binary)
        .args(["makedb", "--in"])
        .arg(root.join("tests/data/search_output/reference.faa"))
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
    for (mode, flag) in [
        ("lin20", "--linclust-20"),
        ("shapes6", "--shapes-6x10"),
        ("shapes30", "--shapes-30x10"),
    ] {
        let result = std::process::Command::new(binary)
            .arg("makeidx")
            .arg("-d")
            .arg(&database)
            .args(["-p1", flag])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{flag}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let actual = std::fs::read(database.with_file_name("reference.dmnd.seed_idx")).unwrap();
        let expected =
            std::fs::read(root.join(format!("tests/data/query_indexed/index-{mode}.seed_idx")))
                .unwrap();
        assert_eq!(actual, expected, "{flag} did not produce the source index");
    }
    std::fs::remove_dir_all(temp).unwrap();
}
