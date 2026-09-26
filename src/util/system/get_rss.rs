//! Resident-set-size helpers translated from `util/system/getRSS.cpp`.

#[cfg(unix)]
use std::ffi::{c_int, c_long, c_void};

/// C++ `getCurrentRSS`.
pub fn get_current_rss() -> usize {
    #[cfg(target_os = "linux")]
    {
        return std::fs::read_to_string("/proc/self/statm")
            .ok()
            .and_then(|statm| current_rss_from_statm(&statm, page_size()))
            .unwrap_or(0);
    }

    #[cfg(target_os = "macos")]
    unsafe {
        let mut info = MachTaskBasicInfo::default();
        let mut count = MACH_TASK_BASIC_INFO_COUNT;
        if task_info(
            mach_task_self_,
            MACH_TASK_BASIC_INFO,
            (&mut info as *mut MachTaskBasicInfo).cast(),
            &mut count,
        ) != KERN_SUCCESS
        {
            return 0;
        }
        return info.resident_size;
    }

    #[cfg(windows)]
    unsafe {
        let mut info = ProcessMemoryCounters::default();
        if get_process_memory_info(
            get_current_process(),
            &mut info,
            std::mem::size_of::<ProcessMemoryCounters>() as u32,
        ) == 0
        {
            return 0;
        }
        return info.working_set_size;
    }

    #[cfg(not(any(target_os = "linux", target_os = "macos", windows)))]
    {
        0
    }
}

/// C++ `getPeakRSS`.
pub fn get_peak_rss() -> usize {
    #[cfg(unix)]
    unsafe {
        // `ru_maxrss` follows two `timeval` values on the supported Unix
        // layouts. The oversized scratch storage also gives `getrusage` room
        // for platform-specific trailing counters.
        let mut usage = [0 as c_long; 64];
        if getrusage(RUSAGE_SELF, usage.as_mut_ptr().cast()) != 0 {
            return 0;
        }
        let max_rss = usage[4].max(0) as usize;
        #[cfg(target_os = "macos")]
        return max_rss;
        #[cfg(not(target_os = "macos"))]
        return max_rss.saturating_mul(1024);
    }

    #[cfg(windows)]
    unsafe {
        let mut info = ProcessMemoryCounters::default();
        if get_process_memory_info(
            get_current_process(),
            &mut info,
            std::mem::size_of::<ProcessMemoryCounters>() as u32,
        ) == 0
        {
            return 0;
        }
        return info.peak_working_set_size;
    }

    #[cfg(not(any(unix, windows)))]
    {
        0
    }
}

#[cfg(target_os = "linux")]
fn current_rss_from_statm(statm: &str, page_size: usize) -> Option<usize> {
    statm
        .split_whitespace()
        .nth(1)?
        .parse::<usize>()
        .ok()?
        .checked_mul(page_size)
}

#[cfg(target_os = "linux")]
fn page_size() -> usize {
    let size = unsafe { sysconf(SC_PAGESIZE) };
    usize::try_from(size)
        .ok()
        .filter(|&size| size > 0)
        .unwrap_or(0)
}

#[cfg(unix)]
const RUSAGE_SELF: c_int = 0;
#[cfg(target_os = "linux")]
const SC_PAGESIZE: c_int = 30;

#[cfg(unix)]
unsafe extern "C" {
    fn getrusage(who: c_int, usage: *mut c_void) -> c_int;
}

#[cfg(target_os = "linux")]
unsafe extern "C" {
    fn sysconf(name: c_int) -> c_long;
}

#[cfg(target_os = "macos")]
const KERN_SUCCESS: c_int = 0;
#[cfg(target_os = "macos")]
const MACH_TASK_BASIC_INFO: c_int = 20;
#[cfg(target_os = "macos")]
const MACH_TASK_BASIC_INFO_COUNT: u32 =
    (std::mem::size_of::<MachTaskBasicInfo>() / std::mem::size_of::<u32>()) as u32;

#[cfg(target_os = "macos")]
#[repr(C)]
#[derive(Default)]
struct MachTaskBasicInfo {
    virtual_size: usize,
    resident_size: usize,
    resident_size_max: usize,
    user_time: [i32; 2],
    system_time: [i32; 2],
    policy: i32,
    suspend_count: i32,
}

#[cfg(target_os = "macos")]
unsafe extern "C" {
    static mach_task_self_: u32;
    fn task_info(target: u32, flavor: c_int, info: *mut c_int, count: *mut u32) -> c_int;
}

#[cfg(windows)]
#[repr(C)]
#[derive(Default)]
struct ProcessMemoryCounters {
    cb: u32,
    page_fault_count: u32,
    peak_working_set_size: usize,
    working_set_size: usize,
    quota_peak_paged_pool_usage: usize,
    quota_paged_pool_usage: usize,
    quota_peak_non_paged_pool_usage: usize,
    quota_non_paged_pool_usage: usize,
    pagefile_usage: usize,
    peak_pagefile_usage: usize,
}

#[cfg(windows)]
#[link(name = "kernel32")]
unsafe extern "system" {
    #[link_name = "GetCurrentProcess"]
    fn get_current_process() -> *mut std::ffi::c_void;
}

#[cfg(windows)]
#[link(name = "psapi")]
unsafe extern "system" {
    #[link_name = "GetProcessMemoryInfo"]
    fn get_process_memory_info(
        process: *mut std::ffi::c_void,
        counters: *mut ProcessMemoryCounters,
        size: u32,
    ) -> i32;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(target_os = "linux")]
    #[test]
    fn parses_linux_statm_in_bytes() {
        assert_eq!(current_rss_from_statm("100 7 3 2\n", 4096), Some(28_672));
        assert_eq!(current_rss_from_statm("100", 4096), None);
        assert_eq!(current_rss_from_statm("100 invalid", 4096), None);
    }

    #[test]
    fn rss_functions_are_callable() {
        let current = get_current_rss();
        let peak = get_peak_rss();
        #[cfg(any(target_os = "linux", target_os = "macos", windows))]
        {
            assert!(current > 0);
            assert!(peak > 0);
        }
    }
}
