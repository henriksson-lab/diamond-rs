use std::fs;
use std::path::Path;

use super::{get_current_rss, get_peak_rss};

#[cfg(windows)]
use std::ffi::c_void;
#[cfg(not(windows))]
use std::ffi::{c_char, c_void, CString};
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
    #[cfg(windows)]
    unsafe {
        let attribute = match color {
            Color::Red => FOREGROUND_RED,
            Color::Green => FOREGROUND_GREEN,
            Color::Yellow => FOREGROUND_RED | FOREGROUND_GREEN,
        };
        // The upstream Windows implementation always addresses stdout and
        // ignores `err`; preserve that behavior.
        let _ = err;
        set_console_text_attribute(get_std_handle(STD_OUTPUT_HANDLE), attribute);
    }
    #[cfg(not(windows))]
    {
        let sequence = color_escape_sequence(color);
        if err {
            eprint!("{sequence}");
        } else {
            print!("{sequence}");
        }
    }
}

pub fn reset_color(err: bool) {
    #[cfg(windows)]
    unsafe {
        let _ = err;
        set_console_text_attribute(
            get_std_handle(STD_OUTPUT_HANDLE),
            FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE,
        );
    }
    #[cfg(not(windows))]
    {
        if err {
            eprint!("{}", reset_color_escape_sequence());
        } else {
            print!("{}", reset_color_escape_sequence());
        }
    }
}

#[cfg(not(windows))]
fn color_escape_sequence(color: Color) -> &'static str {
    match color {
        Color::Red => "\x1b[31m",
        Color::Green => "\x1b[32m",
        Color::Yellow => "\x1b[1;33m",
    }
}

#[cfg(not(windows))]
fn reset_color_escape_sequence() -> &'static str {
    "\x1b[0;39m"
}

#[cfg(windows)]
const STD_OUTPUT_HANDLE: u32 = -11i32 as u32;
#[cfg(windows)]
const FOREGROUND_BLUE: u16 = 0x0001;
#[cfg(windows)]
const FOREGROUND_GREEN: u16 = 0x0002;
#[cfg(windows)]
const FOREGROUND_RED: u16 = 0x0004;

#[cfg(windows)]
#[link(name = "kernel32")]
unsafe extern "system" {
    #[link_name = "GetStdHandle"]
    fn get_std_handle(which: u32) -> *mut c_void;
    #[link_name = "SetConsoleTextAttribute"]
    fn set_console_text_attribute(console: *mut c_void, attributes: u16) -> i32;
}

#[cfg(windows)]
fn widen_utf8(value: &str) -> Vec<u16> {
    value.encode_utf16().collect()
}

#[cfg(windows)]
fn narrow_utf8(value: &[u16]) -> String {
    String::from_utf16_lossy(value)
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

#[cfg(target_os = "linux")]
const _SC_LEVEL3_CACHE_SIZE: c_int = 194;

#[cfg(target_os = "linux")]
unsafe extern "C" {
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

// POSIX declares `open` as `int open(const char *, int, ...)`. Rust's standard
// library also treats this runtime symbol as variadic on Unix, even when
// O_RDONLY means that no optional mode argument is passed.
#[cfg(not(windows))]
unsafe extern "C" {
    fn open(pathname: *const c_char, flags: i32, ...) -> i32;
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
            .map_err(|_| {
                unsafe { close(fd) };
                format!("Error calling fstat on file: {filename}")
            })?
            .len() as usize;
        let addr = unsafe { mmap(std::ptr::null_mut(), length, PROT_READ, MAP_SHARED, fd, 0) };
        if addr as isize == -1 {
            unsafe { close(fd) };
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
    #[cfg(unix)]
    let result = {
        use std::os::unix::fs::DirBuilderExt;
        let mut builder = fs::DirBuilder::new();
        builder.mode(0o755).create(dir)
    };
    #[cfg(not(unix))]
    let result = fs::create_dir(dir);

    match result {
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

    #[cfg(not(windows))]
    {
        let normalized = lex_normalize_posix(&joined);
        if treat_as_dir {
            (normalized, String::new())
        } else {
            (parent_dir_posix(&normalized), base)
        }
    }

    #[cfg(windows)]
    {
        let full = match std::path::absolute(fp) {
            Ok(path) => path,
            Err(_) => return (String::new(), String::new()),
        };
        let full_string = narrow_utf8(&widen_utf8(&full.to_string_lossy()));
        if treat_as_dir {
            return (full_string, String::new());
        }
        let parent = full
            .parent()
            .map(|path| narrow_utf8(&widen_utf8(&path.to_string_lossy())))
            .unwrap_or_default();
        (parent, base)
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
        // C++ `substr(0, npos)` returns the entire string.
        None => file_name.to_string(),
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
        assert_eq!(containing_directory("file.dmnd"), "file.dmnd");
        assert_eq!(last_component("/tmp/x/"), "x");
        assert_eq!(absolute_path("").1, "");
        assert_eq!(absolute_path(".").1, "");
        assert_eq!(absolute_path("dir/").1, "");
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
        let executable = executable_path().unwrap();
        assert!(exists(&executable));
        assert!(exists("Cargo.toml"));
        assert!(file_size("Cargo.toml") > 0);
        assert_eq!(file_size("definitely_missing"), usize::MAX);
        let _ = l3_cache_size();
        log_rss().unwrap();
    }

    #[test]
    fn test_directory_and_existing_extension_helpers() {
        let root = std::env::temp_dir().join(format!(
            "diamond-rs-system-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let root_string = root.to_string_lossy();
        mkdir(&root_string).unwrap();
        mkdir(&root_string).unwrap();
        assert!(exists(&root_string));
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            assert_eq!(
                fs::metadata(&root).unwrap().permissions().mode() & 0o777,
                0o755
            );
        }

        let base = root.join("archive");
        let extended = root.join("archive.daa");
        fs::write(&extended, b"DAA").unwrap();
        assert_eq!(file_size(extended.to_str().unwrap()), 3);
        assert_eq!(
            auto_append_extension_if_exists(base.to_str().unwrap(), ".daa"),
            extended.to_string_lossy()
        );

        fs::remove_file(extended).unwrap();
        rmdir(&root_string).unwrap();
        assert!(!exists(&root_string));
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

    #[test]
    #[cfg(not(windows))]
    fn test_terminal_escape_sequences_match_cpp() {
        assert_eq!(color_escape_sequence(Color::Red), "\x1b[31m");
        assert_eq!(color_escape_sequence(Color::Green), "\x1b[32m");
        assert_eq!(color_escape_sequence(Color::Yellow), "\x1b[1;33m");
        assert_eq!(reset_color_escape_sequence(), "\x1b[0;39m");
    }
}
