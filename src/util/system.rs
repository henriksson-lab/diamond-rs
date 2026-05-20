use std::fs;
use std::path::Path;
#[cfg(windows)]
use std::path::PathBuf;

#[cfg(not(windows))]
use std::ffi::{c_void, CString};
#[cfg(target_os = "linux")]
use std::os::raw::{c_int, c_long};

#[cfg(windows)]
pub const PATH_SEPARATOR: char = '\\';
#[cfg(not(windows))]
pub const PATH_SEPARATOR: char = '/';

#[cfg(windows)]
pub const DEFAULT_LINE_DELIMITER: &str = "\r\n";
#[cfg(not(windows))]
pub const DEFAULT_LINE_DELIMITER: &str = "\n";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Color {
    Red,
    Green,
    Yellow,
}

pub fn set_color(color: Color, err: bool) {
    let code = match color {
        Color::Red => "31",
        Color::Green => "32",
        Color::Yellow => "1;33",
    };
    if err {
        eprint!("\x1b[{code}m");
    } else {
        print!("\x1b[{code}m");
    }
}

pub fn reset_color(err: bool) {
    if err {
        eprint!("\x1b[0;39m");
    } else {
        print!("\x1b[0;39m");
    }
}

pub fn executable_path() -> Result<String, String> {
    std::env::current_exe()
        .map(|p| p.to_string_lossy().into_owned())
        .map_err(|e| e.to_string())
}

pub fn exists(file_name: &str) -> bool {
    Path::new(file_name).exists()
}

pub fn auto_append_extension(str_: &mut String, ext: &str) {
    if !str_.is_empty() && !str_.ends_with(ext) {
        str_.push_str(ext);
    }
}

pub fn auto_append_extension_if_exists(str_: &str, ext: &str) -> String {
    let candidate = format!("{str_}{ext}");
    if !str_.ends_with(ext) && exists(&candidate) {
        candidate
    } else {
        str_.to_string()
    }
}

pub fn get_current_rss() -> usize {
    #[cfg(target_os = "linux")]
    {
        let statm = match fs::read_to_string("/proc/self/statm") {
            Ok(s) => s,
            Err(_) => return 0,
        };
        let rss_pages = match statm
            .split_whitespace()
            .nth(1)
            .and_then(|s| s.parse::<usize>().ok())
        {
            Some(x) => x,
            None => return 0,
        };
        rss_pages * page_size()
    }
    #[cfg(not(target_os = "linux"))]
    {
        0
    }
}

pub fn get_peak_rss() -> usize {
    #[cfg(target_os = "linux")]
    {
        let mut rusage = RUsage::default();
        unsafe {
            getrusage(RUSAGE_SELF, &mut rusage);
        }
        rusage.ru_maxrss as usize * 1024
    }
    #[cfg(not(target_os = "linux"))]
    {
        0
    }
}

#[cfg(target_os = "linux")]
const RUSAGE_SELF: c_int = 0;
#[cfg(target_os = "linux")]
const _SC_PAGESIZE: c_int = 30;
#[cfg(target_os = "linux")]
const _SC_LEVEL3_CACHE_SIZE: c_int = 194;

#[cfg(target_os = "linux")]
#[repr(C)]
#[derive(Default)]
struct TimeVal {
    tv_sec: c_long,
    tv_usec: c_long,
}

#[cfg(target_os = "linux")]
#[repr(C)]
#[derive(Default)]
struct RUsage {
    ru_utime: TimeVal,
    ru_stime: TimeVal,
    ru_maxrss: c_long,
    ru_ixrss: c_long,
    ru_idrss: c_long,
    ru_isrss: c_long,
    ru_minflt: c_long,
    ru_majflt: c_long,
    ru_nswap: c_long,
    ru_inblock: c_long,
    ru_oublock: c_long,
    ru_msgsnd: c_long,
    ru_msgrcv: c_long,
    ru_nsignals: c_long,
    ru_nvcsw: c_long,
    ru_nivcsw: c_long,
}

#[cfg(target_os = "linux")]
unsafe extern "C" {
    fn getrusage(who: c_int, usage: *mut RUsage) -> c_int;
    fn sysconf(name: c_int) -> c_long;
}

pub fn log_rss() -> std::io::Result<()> {
    use crate::util::log_stream::log_stream;
    use crate::util::string::convert_size;

    let mut out = log_stream().lock().unwrap();
    out.write("Current RSS: ")?
        .write(convert_size(get_current_rss()))?
        .write(", Peak RSS: ")?
        .write(convert_size(get_peak_rss()))?
        .endl()?;
    Ok(())
}

pub fn file_size(name: &str) -> usize {
    fs::metadata(name)
        .map(|m| m.len() as usize)
        .unwrap_or(usize::MAX)
}

pub fn total_ram() -> f64 {
    #[cfg(target_os = "linux")]
    {
        let meminfo = match fs::read_to_string("/proc/meminfo") {
            Ok(s) => s,
            Err(_) => return 0.0,
        };
        for line in meminfo.lines() {
            if let Some(rest) = line.strip_prefix("MemTotal:") {
                let kb = rest
                    .split_whitespace()
                    .next()
                    .and_then(|s| s.parse::<f64>().ok())
                    .unwrap_or(0.0);
                return kb * 1024.0 / 1e9;
            }
        }
        0.0
    }
    #[cfg(not(target_os = "linux"))]
    {
        0.0
    }
}

#[cfg(not(windows))]
unsafe extern "C" {
    fn open(pathname: *const i8, flags: i32) -> i32;
    fn mmap(
        addr: *mut c_void,
        length: usize,
        prot: i32,
        flags: i32,
        fd: i32,
        offset: isize,
    ) -> *mut c_void;
    fn munmap(addr: *mut c_void, length: usize) -> i32;
    fn close(fd: i32) -> i32;
}

pub fn mmap_file(filename: &str) -> Result<(*mut u8, usize, i32), String> {
    #[cfg(windows)]
    {
        let _ = filename;
        Err("Memory mapping not supported on Windows.".to_string())
    }
    #[cfg(not(windows))]
    {
        const O_RDONLY: i32 = 0;
        const PROT_READ: i32 = 1;
        const MAP_SHARED: i32 = 1;
        let c_filename =
            CString::new(filename).map_err(|_| format!("Error opening file: {filename}"))?;
        let fd = unsafe { open(c_filename.as_ptr(), O_RDONLY) };
        if fd == -1 {
            return Err(format!("Error opening file: {filename}"));
        }
        let length = fs::metadata(filename)
            .map_err(|_| format!("Error calling fstat on file: {filename}"))?
            .len() as usize;
        let addr = unsafe { mmap(std::ptr::null_mut(), length, PROT_READ, MAP_SHARED, fd, 0) };
        if addr as isize == -1 {
            return Err(format!("Error calling mmap on file: {filename}"));
        }
        Ok((addr as *mut u8, length, fd))
    }
}

pub fn unmap_file(ptr: *mut u8, size: usize, fd: i32) {
    #[cfg(windows)]
    {
        let _ = (ptr, size, fd);
    }
    #[cfg(not(windows))]
    unsafe {
        munmap(ptr as *mut c_void, size);
        close(fd);
    }
}

pub fn l3_cache_size() -> usize {
    #[cfg(target_os = "linux")]
    {
        let s = unsafe { sysconf(_SC_LEVEL3_CACHE_SIZE) };
        if s == -1 {
            0
        } else {
            s as usize
        }
    }
    #[cfg(not(target_os = "linux"))]
    {
        0
    }
}

pub fn mkdir(dir: &str) -> Result<(), String> {
    match fs::create_dir(dir) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => Ok(()),
        Err(_) => Err(format!("could not create temporary directory {dir}")),
    }
}

pub fn rmdir(dir: &str) -> Result<(), String> {
    fs::remove_dir(dir).map_err(|e| e.to_string())
}

pub fn absolute_path(file_path: &str) -> (String, String) {
    let fp = if file_path.is_empty() { "." } else { file_path };
    let base = last_component(fp);
    let treat_as_dir = ends_with_sep(fp) || base == "." || base == "..";

    #[cfg(not(windows))]
    let joined = if is_abs_posix(fp) {
        fp.to_string()
    } else {
        let cwd = get_cwd_posix();
        if cwd.is_empty() {
            return (String::new(), String::new());
        }
        format!("{cwd}/{fp}")
    };

    #[cfg(windows)]
    let joined = if is_absolute_path_for_absolute_path(fp) {
        PathBuf::from(fp).to_string_lossy().into_owned()
    } else {
        std::env::current_dir()
            .unwrap_or_else(|_| PathBuf::from("."))
            .join(fp)
            .to_string_lossy()
            .into_owned()
    };

    let normalized = lex_normalize_posix(&joined);
    if treat_as_dir {
        (normalized, String::new())
    } else {
        (parent_dir_posix(&normalized), base)
    }
}

pub fn is_absolute_path(path: &str) -> bool {
    if path.is_empty() {
        return false;
    }
    let bytes = path.as_bytes();
    bytes[0] == b'/'
        || bytes[0] == b'\\'
        || (bytes.len() >= 3
            && bytes[0].is_ascii_alphabetic()
            && bytes[1] == b':'
            && (bytes[2] == b'/' || bytes[2] == b'\\'))
}

pub fn stdout_is_a_tty() -> bool {
    use std::io::IsTerminal;
    std::io::stdout().is_terminal()
}

pub fn containing_directory(file_name: &str) -> String {
    match file_name.rfind(PATH_SEPARATOR) {
        Some(pos) => file_name[..pos].to_string(),
        None => String::new(),
    }
}

fn page_size() -> usize {
    #[cfg(target_os = "linux")]
    {
        unsafe { sysconf(_SC_PAGESIZE) as usize }
    }
    #[cfg(not(target_os = "linux"))]
    {
        4096
    }
}

fn is_sep_char(c: char) -> bool {
    #[cfg(windows)]
    {
        c == '\\' || c == '/'
    }
    #[cfg(not(windows))]
    {
        c == '/'
    }
}

fn ends_with_sep(s: &str) -> bool {
    s.chars().last().is_some_and(is_sep_char)
}

#[cfg(windows)]
fn is_absolute_path_for_absolute_path(path: &str) -> bool {
    is_absolute_path(path)
}

#[cfg(not(windows))]
pub fn is_abs_posix(path: &str) -> bool {
    !path.is_empty() && path.starts_with('/')
}

#[cfg(not(windows))]
pub fn get_cwd_posix() -> String {
    std::env::current_dir()
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_default()
}

fn last_component(path: &str) -> String {
    if path.is_empty() {
        return String::new();
    }
    let trimmed = path.trim_end_matches(is_sep_char);
    trimmed
        .rsplit(is_sep_char)
        .next()
        .unwrap_or_default()
        .to_string()
}

#[cfg(not(windows))]
pub fn lex_normalize_posix(path: &str) -> String {
    let abs = is_abs_posix(path);
    let mut stack = Vec::new();
    for token in path.split('/') {
        if token.is_empty() || token == "." {
            continue;
        }
        if token == ".." {
            if !stack.is_empty() {
                stack.pop();
            } else if !abs {
                stack.push("..");
            }
        } else {
            stack.push(token);
        }
    }

    let mut out = String::new();
    if abs {
        out.push('/');
    }
    out.push_str(&stack.join("/"));
    if out.is_empty() {
        if abs {
            "/".to_string()
        } else {
            ".".to_string()
        }
    } else {
        out
    }
}

#[cfg(windows)]
fn lex_normalize_posix(path: &str) -> String {
    let abs = path.starts_with('/');
    let mut stack = Vec::new();
    for token in path.split('/') {
        if token.is_empty() || token == "." {
            continue;
        }
        if token == ".." {
            if !stack.is_empty() {
                stack.pop();
            } else if !abs {
                stack.push("..");
            }
        } else {
            stack.push(token);
        }
    }

    let mut out = String::new();
    if abs {
        out.push('/');
    }
    out.push_str(&stack.join("/"));
    if out.is_empty() {
        if abs {
            "/".to_string()
        } else {
            ".".to_string()
        }
    } else {
        out
    }
}

#[cfg(not(windows))]
pub fn parent_dir_posix(abs_path: &str) -> String {
    if abs_path == "/" {
        return "/".to_string();
    }
    match abs_path.rfind('/') {
        Some(0) => "/".to_string(),
        Some(pos) => abs_path[..pos].to_string(),
        None => ".".to_string(),
    }
}

#[cfg(windows)]
fn parent_dir_posix(abs_path: &str) -> String {
    if abs_path == "/" {
        return "/".to_string();
    }
    match abs_path.rfind('/') {
        Some(0) => "/".to_string(),
        Some(pos) => abs_path[..pos].to_string(),
        None => ".".to_string(),
    }
}

#[cfg(windows)]
fn widen_utf8(s: &str) -> Vec<u16> {
    if s.is_empty() {
        Vec::new()
    } else {
        s.encode_utf16().collect()
    }
}

#[cfg(windows)]
fn narrow_utf8(w: &[u16]) -> String {
    if w.is_empty() {
        String::new()
    } else {
        String::from_utf16_lossy(w)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extension_helpers() {
        let mut s = "db".to_string();
        auto_append_extension(&mut s, ".dmnd");
        assert_eq!(s, "db.dmnd");
        auto_append_extension(&mut s, ".dmnd");
        assert_eq!(s, "db.dmnd");
        assert_eq!(
            auto_append_extension_if_exists("definitely_missing", ".dmnd"),
            "definitely_missing"
        );
    }

    #[test]
    fn test_paths() {
        assert!(is_absolute_path("/tmp/x"));
        assert!(is_absolute_path("\\tmp\\x"));
        assert!(is_absolute_path("C:\\tmp\\x"));
        assert!(!is_absolute_path("tmp/x"));
        assert_eq!(containing_directory("/tmp/x"), "/tmp");
        assert_eq!(last_component("/tmp/x/"), "x");
        #[cfg(not(windows))]
        {
            assert!(is_abs_posix("/tmp/x"));
            assert!(!is_abs_posix("tmp/x"));
            assert!(!get_cwd_posix().is_empty());
            assert_eq!(lex_normalize_posix("/tmp/./a/../b"), "/tmp/b");
            assert_eq!(parent_dir_posix("/tmp/b"), "/tmp");
            let (dir, file) = absolute_path("C:\\tmp\\x");
            assert_eq!(file, "C:\\tmp\\x");
            assert!(dir.starts_with('/'));
        }
    }

    #[test]
    fn test_file_and_sizes() {
        assert!(exists("Cargo.toml"));
        assert!(file_size("Cargo.toml") > 0);
        assert_eq!(file_size("definitely_missing"), usize::MAX);
        let _ = l3_cache_size();
        log_rss().unwrap();
    }

    #[test]
    fn test_mmap_file_and_unmap_file() {
        #[cfg(not(windows))]
        {
            let path = std::env::temp_dir().join(format!(
                "diamond-rs-mmap-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ));
            fs::write(&path, b"abcd").unwrap();
            let path_str = path.to_str().unwrap();
            let (ptr, size, fd) = mmap_file(path_str).unwrap();
            assert_eq!(size, 4);
            unsafe {
                assert_eq!(std::slice::from_raw_parts(ptr, size), b"abcd");
            }
            unmap_file(ptr, size, fd);
            fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn test_stdout_is_a_tty_is_callable() {
        let _ = stdout_is_a_tty();
    }
}
