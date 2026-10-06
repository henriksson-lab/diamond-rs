//! Resident-set-size helpers translated from `util/system/getRSS.cpp`.

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
        use mach2::task::task_info;
        use mach2::task_info::{
            mach_task_basic_info, MACH_TASK_BASIC_INFO, MACH_TASK_BASIC_INFO_COUNT,
        };
        use mach2::traps::mach_task_self;

        let mut info = mach_task_basic_info::default();
        let mut count = MACH_TASK_BASIC_INFO_COUNT;
        if task_info(
            mach_task_self(),
            MACH_TASK_BASIC_INFO,
            (&mut info as *mut mach_task_basic_info).cast(),
            &mut count,
        ) != mach2::kern_return::KERN_SUCCESS
        {
            return 0;
        }
        return info.resident_size;
    }

    #[cfg(windows)]
    unsafe {
        use windows_sys::Win32::System::ProcessStatus::{
            K32GetProcessMemoryInfo, PROCESS_MEMORY_COUNTERS,
        };
        use windows_sys::Win32::System::Threading::GetCurrentProcess;

        let mut info: PROCESS_MEMORY_COUNTERS = std::mem::zeroed();
        info.cb = std::mem::size_of::<PROCESS_MEMORY_COUNTERS>() as u32;
        if K32GetProcessMemoryInfo(GetCurrentProcess(), &mut info, info.cb) == 0 {
            return 0;
        }
        return info.WorkingSetSize;
    }

    #[cfg(not(any(target_os = "linux", target_os = "macos", windows)))]
    {
        0
    }
}

/// C++ `getPeakRSS`.
pub fn get_peak_rss() -> usize {
    #[cfg(unix)]
    {
        use nix::sys::resource::{getrusage, UsageWho};

        let Ok(usage) = getrusage(UsageWho::RUSAGE_SELF) else {
            return 0;
        };
        let max_rss = usage.max_rss().max(0) as usize;
        #[cfg(target_vendor = "apple")]
        return max_rss;
        #[cfg(not(target_vendor = "apple"))]
        return max_rss.saturating_mul(1024);
    }

    #[cfg(windows)]
    unsafe {
        use windows_sys::Win32::System::ProcessStatus::{
            K32GetProcessMemoryInfo, PROCESS_MEMORY_COUNTERS,
        };
        use windows_sys::Win32::System::Threading::GetCurrentProcess;

        let mut info: PROCESS_MEMORY_COUNTERS = std::mem::zeroed();
        info.cb = std::mem::size_of::<PROCESS_MEMORY_COUNTERS>() as u32;
        if K32GetProcessMemoryInfo(GetCurrentProcess(), &mut info, info.cb) == 0 {
            return 0;
        }
        return info.PeakWorkingSetSize;
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
    rustix::param::page_size()
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

    #[cfg(any(target_os = "linux", target_os = "macos", windows))]
    #[test]
    fn current_rss_observes_touched_memory_in_an_isolated_process() {
        const CHILD_ENV: &str = "DIAMOND_RSS_TEST_CHILD";
        if std::env::var_os(CHILD_ENV).is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .args([
                    "--exact",
                    "util::system::get_rss::tests::current_rss_observes_touched_memory_in_an_isolated_process",
                    "--nocapture",
                ])
                .env(CHILD_ENV, "1")
                .status()
                .expect("launch isolated RSS test");
            assert!(status.success(), "isolated RSS test failed: {status}");
            return;
        }

        let before = get_current_rss();
        let mut allocation = vec![0u8; 32 * 1024 * 1024];
        for byte in allocation.iter_mut().step_by(4096) {
            // A volatile store ensures that every page is physically touched.
            unsafe { std::ptr::write_volatile(byte, 1) };
        }
        std::hint::black_box(&allocation);
        let after = get_current_rss();
        assert!(before > 0 && after > 0);
        assert!(
            after.saturating_sub(before) >= 16 * 1024 * 1024,
            "current RSS did not observe touched pages: before={before}, after={after}"
        );
    }
}
