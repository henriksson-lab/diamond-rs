use std::io;

/// The former built-in test command depended on an adjacent checkout of the
/// upstream C++ source tree. Repository tests are now self-contained; parity
/// comparisons belong in the external benchmark harness instead.
pub fn run() -> io::Result<()> {
    Err(io::Error::other(
        "the source-dependent built-in regression command has been removed; run the external parity harness instead",
    ))
}
