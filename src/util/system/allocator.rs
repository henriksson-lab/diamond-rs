//! Allocator-specific compatibility operations.
//!
//! DIAMOND releases several very large transient buffers between pipeline
//! stages. On glibc, explicitly trimming freed arenas materially lowers RSS
//! and makes the process-level memory limit useful. There is no portable
//! standard-library equivalent, so the native call is deliberately isolated
//! here and is a no-op on other allocators and platforms.

#[inline]
pub(crate) fn trim_freed_heap_pages() {
    #[cfg(all(
        target_os = "linux",
        target_env = "gnu",
        not(feature = "disable-malloc-trim")
    ))]
    unsafe {
        // libc owns the target-specific declaration; this calls the same
        // glibc extension as upstream without maintaining a local ABI block.
        let _ = libc::malloc_trim(0);
    }
}
