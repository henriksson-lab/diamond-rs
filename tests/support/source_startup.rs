//! Exact Source Config release-age diagnostics (config.cpp:924-930).
//! This validates a known protocol prefix, not arbitrary ANSI stripping.

#![allow(dead_code)]

const WARNING: &[u8] = b"Warning: This version of DIAMOND is more than 180 days old. It is recommended to always use the latest version.\n";

/// Expected Source CLI cout/cerr text is written through CRT text streams.
/// Project fixture text only; never normalize captured stdout, stderr, or files.
pub fn source_console_text(expected: &str) -> String {
    if cfg!(windows) {
        expected.replace('\n', "\r\n")
    } else {
        expected.to_owned()
    }
}

/// Source FileSink uses the actual CRT stdout on MSVC, retaining text mode
/// even when its requested fopen mode is "wb". Other targets fdopen stdout
/// with that binary mode. Named files and in-process Vec writers are separate.
pub fn source_file_stdout_text(expected: &str) -> String {
    if cfg!(all(windows, target_env = "msvc")) {
        expected.replace('\n', "\r\n")
    } else {
        expected.to_owned()
    }
}

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
    expected
        .extend_from_slice(source_console_text(std::str::from_utf8(WARNING).unwrap()).as_bytes());
    expected
}

pub fn assert_age_warning(output: &std::process::Output) {
    let warning = age_warning_stderr();
    let warning_text = source_console_text(std::str::from_utf8(WARNING).unwrap());
    if warning.is_empty() {
        assert!(
            !output
                .stderr
                .windows(warning_text.len())
                .any(|bytes| bytes == warning_text.as_bytes()),
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
