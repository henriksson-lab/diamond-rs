use std::fs;
use std::path::Path;

use memmap2::{Mmap, MmapOptions};

use super::{get_current_rss, get_peak_rss};

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
        use windows_sys::Win32::System::Console::{
            GetStdHandle, SetConsoleTextAttribute, FOREGROUND_GREEN, FOREGROUND_RED,
            STD_OUTPUT_HANDLE,
        };

        let attribute = match color {
            Color::Red => FOREGROUND_RED,
            Color::Green => FOREGROUND_GREEN,
            Color::Yellow => FOREGROUND_RED | FOREGROUND_GREEN,
        };
        // The upstream Windows implementation always addresses stdout and
        // ignores `err`; preserve that behavior.
        let _ = err;
        SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), attribute);
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
        use windows_sys::Win32::System::Console::{
            GetStdHandle, SetConsoleTextAttribute, FOREGROUND_BLUE, FOREGROUND_GREEN,
            FOREGROUND_RED, STD_OUTPUT_HANDLE,
        };

        let _ = err;
        SetConsoleTextAttribute(
            GetStdHandle(STD_OUTPUT_HANDLE),
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

/// Open `filename` as an owning, read-only memory map.
///
/// The map owns its platform mapping and releases it on drop.
/// Unlike the former raw Unix-only helper, this intentionally uses memmap2's
/// native Windows implementation as well.
///
/// # Safety
///
/// The caller must ensure that no process truncates or otherwise mutates the
/// file while the returned map is alive. This is the safety condition required
/// by all file-backed memory maps.
pub unsafe fn mmap_file(filename: &str) -> Result<Mmap, String> {
    let file = fs::File::open(filename).map_err(|_| format!("Error opening file: {filename}"))?;
    if file.metadata().map(|metadata| metadata.len()).unwrap_or(0) == 0 {
        return Err(format!(
            "Error mapping file {filename}: empty files cannot be memory-mapped"
        ));
    }
    // SAFETY: this function creates a read-only map and never mutates the
    // backing file itself. The external-mutation requirement is documented on
    // the returned public API above.
    unsafe { MmapOptions::new().map(&file) }
        .map_err(|error| format!("Error mapping file {filename}: {error}"))
}

pub fn l3_cache_size() -> usize {
    #[cfg(target_os = "linux")]
    {
        linux_l3_cache_size_at(Path::new("/sys/devices/system/cpu"))
    }
    #[cfg(not(target_os = "linux"))]
    {
        0
    }
}

#[cfg(target_os = "linux")]
fn linux_l3_cache_size_at(cpu_root: &Path) -> usize {
    let Ok(cpus) = fs::read_dir(cpu_root) else {
        return 0;
    };
    let mut largest = 0;
    for cpu in cpus.flatten() {
        let cpu_name = cpu.file_name();
        let cpu_name = cpu_name.to_string_lossy();
        if !cpu_name.strip_prefix("cpu").is_some_and(|suffix| {
            !suffix.is_empty() && suffix.bytes().all(|byte| byte.is_ascii_digit())
        }) {
            continue;
        }
        let Ok(indices) = fs::read_dir(cpu.path().join("cache")) else {
            continue;
        };
        for index in indices.flatten() {
            if !index.file_name().to_string_lossy().starts_with("index") {
                continue;
            }
            let path = index.path();
            let level = fs::read_to_string(path.join("level")).unwrap_or_default();
            let kind = fs::read_to_string(path.join("type")).unwrap_or_default();
            if level.trim() != "3" || !matches!(kind.trim(), "Data" | "Unified") {
                continue;
            }
            if let Ok(size) = fs::read_to_string(path.join("size")) {
                largest = largest.max(parse_linux_cache_size(&size).unwrap_or(0));
            }
        }
    }
    largest
}

#[cfg(target_os = "linux")]
fn parse_linux_cache_size(value: &str) -> Option<usize> {
    let value = value.trim();
    let digits = value.bytes().take_while(u8::is_ascii_digit).count();
    if digits == 0 {
        return None;
    }
    let amount = value[..digits].parse::<usize>().ok()?;
    let suffix = value[digits..].trim().to_ascii_lowercase();
    let multiplier = match suffix.as_str() {
        "" | "b" => 1,
        "k" | "kb" | "kib" => 1024,
        "m" | "mb" | "mib" => 1024 * 1024,
        "g" | "gb" | "gib" => 1024 * 1024 * 1024,
        _ => return None,
    };
    amount.checked_mul(multiplier)
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

    #[cfg(target_os = "linux")]
    #[test]
    fn parses_linux_cache_sizes_without_overflow() {
        assert_eq!(parse_linux_cache_size("32K\n"), Some(32 * 1024));
        assert_eq!(parse_linux_cache_size("4 MiB"), Some(4 * 1024 * 1024));
        assert_eq!(parse_linux_cache_size("1024"), Some(1024));
        assert_eq!(parse_linux_cache_size(""), None);
        assert_eq!(parse_linux_cache_size("K"), None);
        assert_eq!(parse_linux_cache_size("12XB"), None);
        assert_eq!(parse_linux_cache_size(&format!("{}G", usize::MAX)), None);
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn discovers_largest_shared_l3_cache_from_sysfs_layout() {
        let root = std::env::temp_dir().join(format!(
            "diamond-rs-cache-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let caches = [
            ("cpu0/cache/index0", "1", "Data", "32K"),
            ("cpu0/cache/index3", "3", "Unified", "8M"),
            ("cpu1/cache/index3", "3", "Unified", "8M"),
            ("cpu8/cache/index3", "3", "Data", "32M"),
            ("cpu9/cache/index3", "3", "Instruction", "64M"),
        ];
        for (relative, level, kind, size) in caches {
            let directory = root.join(relative);
            fs::create_dir_all(&directory).unwrap();
            fs::write(directory.join("level"), level).unwrap();
            fs::write(directory.join("type"), kind).unwrap();
            fs::write(directory.join("size"), size).unwrap();
        }
        // Non-CPU directories must not affect discovery.
        fs::create_dir_all(root.join("cpufreq/cache/index3")).unwrap();

        assert_eq!(linux_l3_cache_size_at(&root), 32 * 1024 * 1024);
        fs::remove_dir_all(root).unwrap();
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
    fn test_mmap_file_is_an_owning_slice() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-mmap-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::write(&path, b"abcd").unwrap();
        // SAFETY: this test exclusively owns the file until the map is dropped.
        let map = unsafe { mmap_file(path.to_str().unwrap()) }.unwrap();
        assert_eq!(&map[..], b"abcd");
        drop(map);
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_mmap_file_preserves_empty_file_failure() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-empty-mmap-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::write(&path, []).unwrap();
        // SAFETY: this test exclusively owns the file during the attempted map.
        assert!(unsafe { mmap_file(path.to_str().unwrap()) }.is_err());
        fs::remove_file(path).unwrap();
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
