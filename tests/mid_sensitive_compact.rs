use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn run(args: &[&str]) {
    let output = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args(args)
        .output()
        .expect("run native diamond");
    assert!(
        output.status.success(),
        "diamond {:?} failed:\n{}",
        args,
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn compact_mid_sensitive_preserves_reciprocal_hits_across_index_chunks() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let input = root.join("diamond/src/test/5.faa");
    let work =
        std::env::temp_dir().join(format!("diamond-rs-mid-sensitive-{}", std::process::id()));
    fs::create_dir_all(&work).expect("create temporary test directory");
    let db = work.join("compact");
    let out = work.join("mid-sensitive.tsv");

    run(&[
        "makedb",
        "--in",
        input.to_str().unwrap(),
        "--db",
        db.to_str().unwrap(),
        "--threads",
        "1",
    ]);
    run(&[
        "blastp",
        "--query",
        input.to_str().unwrap(),
        "--db",
        db.to_str().unwrap(),
        "--out",
        out.to_str().unwrap(),
        "--threads",
        "1",
        "--mid-sensitive",
    ]);

    let rows = fs::read_to_string(&out).expect("read blastp output");
    assert_eq!(rows.lines().count(), 1338);
    assert!(rows.lines().any(|line| {
        line == "d3t6ka_\td3crna1\t30.2\t116\t79\t1\t5\t120\t3\t116\t5.74e-17\t65.5"
    }));
    assert!(rows.lines().any(|line| {
        line == "d3crna1\td3t6ka_\t30.1\t113\t77\t1\t3\t113\t5\t117\t1.26e-15\t62.0"
    }));

    let _ = fs::remove_dir_all(work);
}
