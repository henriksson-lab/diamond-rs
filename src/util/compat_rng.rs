//! Compatibility access to the platform C runtime pseudo-random generator.
//!
//! DIAMOND's translated statistical and tantan code relies on the exact
//! process-global `rand`/`srand` state supplied by the active C runtime. The
//! sequence, range, and implementation are platform-dependent, and calls from
//! any thread or library code in the process can interleave with these calls.
//! Seeding also replaces that shared process-wide state. These wrappers exist
//! only to preserve upstream compatibility; they are not suitable for general
//! randomness or security-sensitive use.

/// Advances and returns the platform C runtime's process-global RNG state.
#[inline]
pub(crate) fn c_rand() -> i32 {
    // SAFETY: C `rand` accepts no arguments and has no caller-side safety
    // preconditions. `libc` supplies the target-specific declaration; its
    // global-state semantics are documented above.
    unsafe { libc::rand() }
}

/// Replaces the platform C runtime's process-global RNG state.
#[inline]
pub(crate) fn c_srand(seed: u32) {
    // SAFETY: C `srand` accepts every `u32` seed and has no caller-side safety
    // preconditions. `libc` supplies the target-specific declaration; its
    // global-state effects are documented above.
    unsafe { libc::srand(seed) }
}
