//! Translation of `diamond/src/util/parallel/filestack.{h,cpp}`.

use std::fs::{File, OpenOptions};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::time::Duration;

pub const DEFAULT_FILE_NAME: &str = "default_stack.idx";
pub const DEFAULT_MAX_LINE_LENGTH: usize = 4096;
const MINIMUM_LINE_LENGTH: usize = 8;

/// A newline-delimited, file-backed LIFO stack.
///
/// The process-local mutex corresponds to the C++ `mtx_`; `File::lock` adds
/// the operating-system advisory lock used by the original `fcntl`/
/// `LockFileEx` implementation, so distinct `FileStack` instances and
/// processes serialize compound operations such as `fetch_add`.
#[derive(Debug)]
pub struct FileStack {
    file_name: PathBuf,
    file: Mutex<File>,
    max_line_length: usize,
    mtx: Mutex<()>,
}

impl Default for FileStack {
    /// C++ `FileStack::FileStack()`.
    fn default() -> Self {
        eprintln!("FileStack: Using default file name {DEFAULT_FILE_NAME}");
        Self::new(DEFAULT_FILE_NAME)
    }
}

impl FileStack {
    pub const DEFAULT_MAX_LINE_LENGTH: usize = DEFAULT_MAX_LINE_LENGTH;

    /// C++ `FileStack::FileStack(const string&)`.
    pub fn new(file_name: impl AsRef<Path>) -> Self {
        Self::with_max_line_length(file_name, DEFAULT_MAX_LINE_LENGTH)
    }

    /// C++ `FileStack::FileStack(const string&, int)`.
    pub fn with_max_line_length(file_name: impl AsRef<Path>, maximum_line_length: usize) -> Self {
        let file_name = file_name.as_ref().to_path_buf();
        let file = OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .open(&file_name)
            .unwrap_or_else(|_| panic!("could not open file {}", file_name.display()));
        Self {
            file_name,
            file: Mutex::new(file),
            max_line_length: maximum_line_length.max(MINIMUM_LINE_LENGTH),
            mtx: Mutex::new(()),
        }
    }

    /// C++ `FileStack::size()`.
    pub fn size(&self) -> Result<usize, String> {
        self.with_lock(|| self.size_non_locked())
    }

    /// C++ `FileStack::pop(int64_t&)`, returning the popped value directly.
    pub fn pop_i64(&self) -> Result<i64, String> {
        let mut value = 0;
        self.pop_i64_into(&mut value)
    }

    /// C++ `FileStack::pop(int64_t&)` compatibility form.
    pub fn pop_i64_into(&self, value: &mut i64) -> Result<i64, String> {
        self.with_lock(|| self.pop_non_locked_i64(value))
    }

    /// C++ `FileStack::pop(string&)`.
    pub fn pop_string(&self) -> Result<Option<String>, String> {
        self.with_lock(|| self.pop_non_locked_string(false))
    }

    /// C++ `FileStack::pop(string&, size_t&)`.
    pub fn pop_string_size(&self) -> Result<(Option<String>, usize), String> {
        self.with_lock(|| {
            let value = self.pop_non_locked_string(false)?;
            Ok((value, self.size_non_locked()?))
        })
    }

    /// C++ `FileStack::top(int64_t&)`, returning the value directly.
    pub fn top_i64(&self) -> Result<i64, String> {
        self.with_lock(|| match self.pop_non_locked_string(true)? {
            Some(value) => value.parse::<i64>().map_err(|e| e.to_string()),
            None => Ok(-1),
        })
    }

    /// C++ `FileStack::top(int64_t&)` compatibility form.
    pub fn top_i64_into(&self, value: &mut i64) -> Result<i64, String> {
        *value = self.top_i64()?;
        Ok(*value)
    }

    /// C++ `FileStack::top(string&)`.
    pub fn top_string(&self) -> Result<Option<String>, String> {
        self.with_lock(|| self.pop_non_locked_string(true))
    }

    /// C++ `FileStack::remove(const string&)`.
    pub fn remove(&self, line: &str) -> Result<(), String> {
        self.with_lock(|| {
            let data = std::fs::read_to_string(&self.file_name).map_err(|e| e.to_string())?;
            let tokens = super::multiprocessing::split(&data, '\n');
            let mut output = String::new();
            for token in tokens.into_iter().filter(|token| token != line) {
                output.push_str(&token);
                output.push('\n');
            }
            std::fs::write(&self.file_name, output).map_err(|e| e.to_string())
        })
    }

    /// C++ `FileStack::push(int64_t)`.
    pub fn push_i64(&self, value: i64) -> Result<i64, String> {
        self.with_lock(|| self.push_non_locked_i64(value))
    }

    /// C++ `FileStack::push(const string&)`.
    pub fn push_string(&self, value: &str) -> Result<i64, String> {
        self.with_lock(|| self.push_exclusive(value))
    }

    /// C++ `FileStack::push(const string&, size_t&)`.
    pub fn push_string_size(&self, value: &str) -> Result<(i64, usize), String> {
        self.with_lock(|| {
            let written = self.push_exclusive(value)?;
            Ok((written, self.size_non_locked()?))
        })
    }

    /// C++ `FileStack::fetch_add(int64_t)`.
    pub fn fetch_add(&self, n: i64) -> Result<i64, String> {
        self.with_lock(|| {
            let mut value = -1;
            self.pop_non_locked_i64(&mut value)?;
            if value == -1 {
                value = 0;
            }
            self.push_non_locked_i64(value + n)?;
            Ok(value)
        })
    }

    /// The default-argument form `fetch_add()`.
    pub fn fetch_add_one(&self) -> Result<i64, String> {
        self.fetch_add(1)
    }

    /// C++ `FileStack::get_max_line_length()`.
    pub fn get_max_line_length(&self) -> usize {
        self.max_line_length
    }

    /// C++ `FileStack::set_max_line_length(int)`.
    pub fn set_max_line_length(&mut self, n: usize) -> usize {
        self.max_line_length = n.max(MINIMUM_LINE_LENGTH);
        self.max_line_length
    }

    /// C++ `FileStack::clear()`.
    pub fn clear(&self) -> Result<(), String> {
        self.with_lock(|| std::fs::write(&self.file_name, []).map_err(|e| e.to_string()))
    }

    /// C++ `FileStack::seek(int64_t, int)`.
    pub fn seek(&self, offset: i64, mode: SeekFrom) -> Result<i64, String> {
        let target = match mode {
            SeekFrom::Start(_) => SeekFrom::Start(offset as u64),
            SeekFrom::End(_) => SeekFrom::End(offset),
            SeekFrom::Current(_) => SeekFrom::Current(offset),
        };
        self.file
            .lock()
            .map_err(|e| e.to_string())?
            .seek(target)
            .map(|position| position as i64)
            .map_err(|e| e.to_string())
    }

    /// C++ `FileStack::read(char*, size_t)`.
    pub fn read(&self, buf: &mut [u8]) -> Result<usize, String> {
        self.file
            .lock()
            .map_err(|e| e.to_string())?
            .read(buf)
            .map_err(|e| format!("Error reading file {}: {e}", self.file_name.display()))
    }

    /// C++ `FileStack::write(const char*, size_t)`.
    pub fn write(&self, buf: &[u8]) -> Result<i64, String> {
        self.file
            .lock()
            .map_err(|e| e.to_string())?
            .write(buf)
            .map(|n| n as i64)
            .map_err(|e| format!("Error writing file {}: {e}", self.file_name.display()))
    }

    /// C++ `FileStack::truncate(size_t)`.
    pub fn truncate(&self, size: usize) -> Result<(), String> {
        self.file
            .lock()
            .map_err(|e| e.to_string())?
            .set_len(size as u64)
            .map_err(|e| e.to_string())
    }

    /// C++ `FileStack::poll_query(...)`.
    pub fn poll_query(&self, query: &str, sleep_s: f64, max_iter: usize) -> Result<bool, String> {
        for _ in 0..max_iter {
            if let Some(value) = self.top_string()? {
                if value.contains(query) {
                    return Ok(true);
                }
                if value.contains("STOP") {
                    return Err(format!("STOP on FileStack {}", self.file_name.display()));
                }
            }
            std::thread::sleep(Duration::from_secs_f64(sleep_s));
        }
        Err(format!(
            "Could not discover keyword {} on FileStack {} within {} seconds.",
            query,
            self.file_name.display(),
            max_iter as f64 * sleep_s
        ))
    }

    /// C++ `FileStack::poll_size(...)`.
    pub fn poll_size(&self, size: usize, sleep_s: f64, max_iter: usize) -> Result<bool, String> {
        for _ in 0..max_iter {
            if self.size()? == size {
                return Ok(true);
            }
            std::thread::sleep(Duration::from_secs_f64(sleep_s));
        }
        Err(format!(
            "Could not detect size {} of FileStack {} within {} seconds.",
            size,
            self.file_name.display(),
            max_iter as f64 * sleep_s
        ))
    }

    /// C++ `FileStack::file_name()`.
    pub fn file_name(&self) -> &Path {
        &self.file_name
    }

    /// C++ `FileStack::pop_exclusive(string&)`.
    ///
    /// Like the original, this assumes the caller already owns any required
    /// compound-operation lock.
    pub fn pop_exclusive(&self) -> Result<Option<String>, String> {
        self.pop_non_locked_string(false)
    }

    /// C++ `FileStack::push_exclusive(const string&)`.
    pub fn push_exclusive(&self, value: &str) -> Result<i64, String> {
        let mut file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.file_name)
            .map_err(|e| e.to_string())?;
        file.write_all(value.as_bytes())
            .map_err(|e| e.to_string())?;
        let mut written = value.len() as i64;
        if !value.ends_with('\n') {
            file.write_all(b"\n").map_err(|e| e.to_string())?;
            written += 1;
        }
        Ok(written)
    }

    /// C++ `FileStack::pop_non_locked(string&, bool, size_t&)`.
    fn pop_non_locked_string(&self, keep: bool) -> Result<Option<String>, String> {
        let data = std::fs::read(&self.file_name).map_err(|e| e.to_string())?;
        if data.is_empty() {
            return Ok(None);
        }

        let jump = data.len().saturating_sub(self.max_line_length);
        let chunk = &data[jump..];
        let mut begin = 0;
        let mut end = chunk.iter().rposition(|&byte| byte == b'\n').unwrap_or(0);
        if end > 0 {
            if let Some(found) = chunk[..end].iter().rposition(|&byte| byte == b'\n') {
                begin = found + 1;
            }
        }
        let line_size = end - begin + 1;
        let value_end = if chunk[begin + line_size - 1] == b'\n' {
            begin + line_size - 1
        } else {
            begin + line_size
        };
        let value =
            String::from_utf8(chunk[begin..value_end].to_vec()).map_err(|e| e.to_string())?;
        if !keep {
            end = data.len() - line_size;
            std::fs::OpenOptions::new()
                .write(true)
                .open(&self.file_name)
                .map_err(|e| e.to_string())?
                .set_len(end as u64)
                .map_err(|e| e.to_string())?;
        }
        Ok(Some(value))
    }

    /// C++ `FileStack::pop_non_locked(int64_t&)`.
    fn pop_non_locked_i64(&self, value: &mut i64) -> Result<i64, String> {
        match self.pop_non_locked_string(false)? {
            Some(text) => {
                *value = text.parse::<i64>().map_err(|e| e.to_string())?;
            }
            None => *value = -1,
        }
        Ok(*value)
    }

    /// C++ `FileStack::push_non_locked(int64_t)`.
    fn push_non_locked_i64(&self, value: i64) -> Result<i64, String> {
        self.push_exclusive(&value.to_string())
    }

    fn size_non_locked(&self) -> Result<usize, String> {
        let data = std::fs::read(&self.file_name).map_err(|e| e.to_string())?;
        Ok(data.iter().filter(|&&byte| byte == b'\n').count())
    }

    /// C++ `lock`/`unlock`, expressed as a scoped operation so error paths
    /// cannot leave either lock held.
    fn with_lock<T>(&self, operation: impl FnOnce() -> Result<T, String>) -> Result<T, String> {
        let _thread_lock = self.mtx.lock().map_err(|e| e.to_string())?;
        let file = self.file.lock().map_err(|e| e.to_string())?;
        file.lock().map_err(|e| {
            format!(
                "could not put lock on file {}: {e}",
                self.file_name.display()
            )
        })?;
        let result = operation();
        let unlock = file
            .unlock()
            .map_err(|e| format!("could not unlock file {}: {e}", self.file_name.display()));
        match (result, unlock) {
            (Err(error), _) => Err(error),
            (Ok(_), Err(error)) => Err(error),
            (Ok(value), Ok(())) => Ok(value),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    fn path(tag: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-filestack-mirror-{}-{tag}.idx",
            std::process::id()
        ))
    }

    #[test]
    fn remove_preserves_empty_records_like_cpp_split_and_rewrite() {
        let path = path("empty-records");
        let _ = std::fs::remove_file(&path);
        std::fs::write(&path, b"alpha\n\nbeta\n").unwrap();
        let stack = FileStack::new(&path);
        stack.remove("alpha").unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"\nbeta\n");
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn independent_instances_serialize_fetch_add_with_file_lock() {
        let path = path("file-lock");
        let _ = std::fs::remove_file(&path);
        let first = Arc::new(FileStack::new(&path));
        let second = Arc::new(FileStack::new(&path));
        let mut threads = Vec::new();
        for stack in [first.clone(), second.clone()] {
            threads.push(std::thread::spawn(move || {
                for _ in 0..100 {
                    stack.fetch_add_one().unwrap();
                }
            }));
        }
        for thread in threads {
            thread.join().unwrap();
        }
        assert_eq!(first.top_i64().unwrap(), 200);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn maximum_line_length_controls_cpp_tail_window() {
        let path = path("tail-window");
        let _ = std::fs::remove_file(&path);
        let stack = FileStack::with_max_line_length(&path, 8);
        stack.push_string("123456789").unwrap();
        assert_eq!(stack.pop_string().unwrap().as_deref(), Some("3456789"));
        assert_eq!(std::fs::read(&path).unwrap(), b"12");
        let _ = std::fs::remove_file(path);
    }
}
