//! Exact Source Config release-age diagnostics (config.cpp:924-930).
//! This validates a known protocol prefix, not arbitrary ANSI stripping.

#![allow(dead_code)]

const WARNING: &[u8] = b"Warning: This version of DIAMOND is more than 180 days old. It is recommended to always use the latest version.\n";

fn age_warning() -> bool {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_secs()
        > 1_772_912_406 + 180 * 24 * 60 * 60
}

pub fn age_warning_stderr() -> Vec<u8> {
    if !age_warning() {
        return Vec::new();
    }
    // Source uses console attributes on Windows, escape bytes elsewhere.
    let mut expected = if cfg!(windows) {
        Vec::new()
    } else {
        b"\x1b[1;33m".to_vec()
    };
    expected.extend_from_slice(WARNING);
    expected
}

pub fn assert_age_warning(output: &std::process::Output) {
    let warning = age_warning_stderr();
    if warning.is_empty() {
        assert!(
            !output
                .stderr
                .windows(WARNING.len())
                .any(|bytes| bytes == WARNING),
            "unexpected Source release-age warning"
        );
    } else {
        assert_eq!(
            output
                .stderr
                .windows(warning.len())
                .filter(|bytes| *bytes == warning.as_slice())
                .count(),
            1,
            "missing, repeated, or changed Source release-age warning"
        );
    }
}

fn stdout_reset() -> &'static [u8] {
    if age_warning() && !cfg!(windows) {
        b"\x1b[0;39m"
    } else {
        b""
    }
}

pub fn stdout_body(output: &std::process::Output) -> &[u8] {
    assert_age_warning(output);
    // reset_color() defaults to stdout even though set_color(YELLOW,true)
    // and the warning go to stderr. Windows changes console attributes only.
    let prefix = stdout_reset();
    output
        .stdout
        .strip_prefix(prefix)
        .expect("missing Source startup stdout prefix")
}

/// Source argc==2 documentation is flushed before Config emits its warning;
/// with no data OutputFile owner, the reset remains at the end of stdout.
pub fn stdout_body_before_suffix(output: &std::process::Output) -> &[u8] {
    assert_age_warning(output);
    output
        .stdout
        .strip_suffix(stdout_reset())
        .expect("missing Source startup stdout suffix")
}
