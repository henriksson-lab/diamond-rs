use std::io;
use std::time::Instant;

/// Run the built-in regression test suite.
///
/// Runs blastp on test data via FFI and verifies against reference output files.
pub fn run() -> io::Result<()> {
    let start = Instant::now();
    let td = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test");

    eprintln!("Running regression tests...");
    eprintln!();

    let mut passed = 0u32;
    let mut total = 0u32;

    // Helper closure to run a test
    let mut run_test = |name: &str, args: &[&str], ref_file: &str| {
        total += 1;
        let ref_output = match std::fs::read(ref_file) {
            Ok(data) => data,
            Err(e) => {
                eprintln!("{:<28} [ \x1b[31mError: {}\x1b[0m ]", name, e);
                return;
            }
        };
        let out_file = std::env::temp_dir().join(format!("diamond_test_{}.out", total));
        let mut full_args: Vec<&str> = vec!["diamond"];
        full_args.extend_from_slice(args);
        full_args.push("-o");
        let out_str = out_file.to_string_lossy().to_string();
        full_args.push(&out_str);

        let code = crate::ffi::run(&full_args);
        if code == 0 {
            let actual = std::fs::read(&out_file).unwrap_or_default();
            if actual == ref_output {
                eprintln!("{:<28} [ \x1b[32mPassed\x1b[0m ]", name);
                passed += 1;
            } else {
                eprintln!("{:<28} [ \x1b[31mFailed\x1b[0m ]", name);
            }
        } else {
            eprintln!("{:<28} [ \x1b[31mError ({})\x1b[0m ]", name, code);
        }
        let _ = std::fs::remove_file(&out_file);
    };

    let q1 = format!("{}/1.faa", td);
    let d2 = format!("{}/2.faa", td);
    let q3 = format!("{}/3.faa", td);
    let d4 = format!("{}/4.faa", td);

    run_test(
        "blastp (default)",
        &["blastp", "-q", &q1, "-d", &d2, "-p1"],
        &format!("{}/blastp.out", td),
    );

    run_test(
        "blastp (mid-sensitive)",
        &["blastp", "-q", &q3, "-d", &d4, "--mid-sensitive", "-p1"],
        &format!("{}/blastp-mid-sens.out", td),
    );

    run_test(
        "blastp (pairwise format)",
        &["blastp", "-q", &q1, "-d", &d2, "-f0", "-p1"],
        &format!("{}/blastp-f0.out", td),
    );

    eprintln!();
    eprintln!(
        "#Test cases passed: {}/{} ({:.1}s)",
        passed, total,
        start.elapsed().as_secs_f64()
    );

    if passed == total {
        Ok(())
    } else {
        Err(io::Error::other(format!("{} test(s) failed", total - passed)))
    }
}
