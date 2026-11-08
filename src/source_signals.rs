//! Preserve the executable's inherited SIGPIPE disposition across Rust's
//! runtime initialization. Source's C++ main leaves that disposition intact.
//!
//! This module belongs only to the binary; linking the library installs no
//! initializer and changes no host-process signal policy.

use std::mem::MaybeUninit;
use std::sync::atomic::{AtomicBool, Ordering};

static mut INHERITED_SIGPIPE: MaybeUninit<libc::sigaction> = MaybeUninit::uninit();
static CAPTURED: AtomicBool = AtomicBool::new(false);

unsafe extern "C" fn capture_sigpipe() {
    // Loader initializers run before Rust's runtime sets SIGPIPE to SIG_IGN.
    // Query only: do not install a handler, allocate, or invoke Rust I/O here.
    let mut action: libc::sigaction = unsafe { std::mem::zeroed() };
    if unsafe { libc::sigaction(libc::SIGPIPE, std::ptr::null(), &mut action) } == 0 {
        // The loader is the only writer. Publish the complete action so its
        // mask and flags are retained along with the handler disposition.
        unsafe {
            std::ptr::addr_of_mut!(INHERITED_SIGPIPE)
                .cast::<libc::sigaction>()
                .write(action);
        }
        CAPTURED.store(true, Ordering::Release);
    }
}

// Both supported Unix object formats invoke these C-ABI initializer entries
// before main. #[used] retains the entry through code generation and LTO.
#[used]
#[cfg_attr(
    target_vendor = "apple",
    unsafe(link_section = "__DATA,__mod_init_func")
)]
#[cfg_attr(not(target_vendor = "apple"), unsafe(link_section = ".init_array"))]
static CAPTURE_SIGPIPE: unsafe extern "C" fn() = capture_sigpipe;

pub(crate) fn restore() {
    if CAPTURED.load(Ordering::Acquire) {
        // The initializer completed before main and never mutates this action
        // again. Borrow only its raw pointer; no mutable static reference escapes.
        unsafe {
            libc::sigaction(
                libc::SIGPIPE,
                std::ptr::addr_of!(INHERITED_SIGPIPE).cast::<libc::sigaction>(),
                std::ptr::null_mut(),
            );
        }
    } else {
        // Preserve the normal Source C-runtime default if querying the action
        // was denied during loader initialization.
        unsafe { libc::signal(libc::SIGPIPE, libc::SIG_DFL) };
    }
}
