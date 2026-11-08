//! Translation of `diamond/src/util/parallel/filestack.{h,cpp}`.

use crate::util::io::source_c_path;
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
    errors: crate::util::io::IoErrorOwner,
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

    /// Fallible construction for translated callers which propagate the
    /// Source constructor's open failure rather than unwinding through panic.
    pub fn try_new(file_name: impl AsRef<Path>) -> Result<Self, String> {
        Self::try_with_max_line_length(file_name, DEFAULT_MAX_LINE_LENGTH)
    }

    /// C++ `FileStack::FileStack(const string&, int)`.
    pub fn with_max_line_length(file_name: impl AsRef<Path>, maximum_line_length: usize) -> Self {
        Self::try_with_max_line_length(file_name, maximum_line_length)
            .unwrap_or_else(|error| panic!("{error}"))
    }

    fn try_with_max_line_length(
        file_name: impl AsRef<Path>,
        maximum_line_length: usize,
    ) -> Result<Self, String> {
        Self::try_with_owner(file_name.as_ref(), maximum_line_length, Default::default())
    }

    pub fn try_new_owned(
        file_name: impl AsRef<Path>,
        errors: crate::util::io::IoErrorOwner,
    ) -> Result<Self, String> {
        Self::try_with_owner(file_name.as_ref(), DEFAULT_MAX_LINE_LENGTH, errors)
    }

    fn try_with_owner(
        file_name: &Path,
        maximum_line_length: usize,
        errors: crate::util::io::IoErrorOwner,
    ) -> Result<Self, String> {
        let file_name = file_name.to_path_buf();
        let mut options = OpenOptions::new();
        options.create(true).read(true).write(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            options.mode(0o664);
        }
        #[cfg(windows)]
        {
            use std::os::windows::fs::OpenOptionsExt;
            use windows_sys::Win32::Storage::FileSystem::{FILE_SHARE_READ, FILE_SHARE_WRITE};
            options.share_mode(FILE_SHARE_READ | FILE_SHARE_WRITE);
        }
        let file = options
            .open(source_c_path(&file_name).as_ref())
            .map_err(|_| {
                errors
                    .capture::<()>(Err(std::io::Error::other(
                        crate::util::io::IoError::path_message(
                            if cfg!(windows) {
                                b"Error opening file "
                            } else {
                                b"could not open file "
                            },
                            &file_name,
                        ),
                    )))
                    .unwrap_err()
            })?;
        Ok(Self {
            file_name,
            file: Mutex::new(file),
            max_line_length: maximum_line_length.max(MINIMUM_LINE_LENGTH),
            mtx: Mutex::new(()),
            errors,
        })
    }

    /// C++ `FileStack::size()`.
    pub fn size(&self) -> Result<usize, String> {
        self.with_lock(|file| self.size_non_locked(file))
    }

    /// C++ `FileStack::pop(int64_t&)`, returning the popped value directly.
    pub fn pop_i64(&self) -> Result<i64, String> {
        let mut value = 0;
        self.pop_i64_into(&mut value)
    }

    /// C++ `FileStack::pop(int64_t&)` compatibility form.
    pub fn pop_i64_into(&self, value: &mut i64) -> Result<i64, String> {
        self.with_lock(|file| self.pop_non_locked_i64(file, value))
    }

    /// C++ `FileStack::pop(string&)`.
    pub fn pop_string(&self) -> Result<Option<String>, String> {
        self.with_lock(|file| self.pop_non_locked_string(file, false))
    }

    /// C++ `FileStack::pop(string&, size_t&)`.
    pub fn pop_string_size(&self) -> Result<(Option<String>, usize), String> {
        self.with_lock(|file| {
            let value = self.pop_non_locked_string(file, false)?;
            Ok((value, self.size_non_locked(file)?))
        })
    }

    /// C++ `FileStack::top(int64_t&)`, returning the value directly.
    pub fn top_i64(&self) -> Result<i64, String> {
        self.with_lock(|file| match self.pop_non_locked_string(file, true)? {
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
        self.with_lock(|file| self.pop_non_locked_string(file, true))
    }

    /// C++ `FileStack::remove(const string&)`.
    pub fn remove(&self, line: &str) -> Result<(), String> {
        self.remove_bytes(line.as_bytes())
    }

    pub fn remove_bytes(&self, line: &[u8]) -> Result<(), String> {
        self.with_lock(|file| {
            file.seek(SeekFrom::Start(0)).map_err(|e| e.to_string())?;
            let mut data = Vec::new();
            file.read_to_end(&mut data).map_err(|e| e.to_string())?;
            let mut output = Vec::new();
            let length = data.len() - usize::from(data.last() == Some(&b'\n'));
            if !data.is_empty() {
                for token in data[..length]
                    .split(|&byte| byte == b'\n')
                    .filter(|token| *token != line)
                {
                    output.extend_from_slice(token);
                    output.push(b'\n');
                }
            }
            file.set_len(0).map_err(|e| e.to_string())?;
            file.seek(SeekFrom::Start(0)).map_err(|e| e.to_string())?;
            file.write_all(&output).map_err(|e| e.to_string())
        })
    }

    pub fn pop_bytes(&self) -> Result<Option<Vec<u8>>, String> {
        self.with_lock(|file| self.pop_non_locked_bytes(file, false))
    }
    pub fn top_bytes(&self) -> Result<Option<Vec<u8>>, String> {
        self.with_lock(|file| self.pop_non_locked_bytes(file, true))
    }
    pub fn pop_bytes_size(&self) -> Result<(Option<Vec<u8>>, usize), String> {
        self.with_lock(|file| {
            Ok((
                self.pop_non_locked_bytes(file, false)?,
                self.size_non_locked(file)?,
            ))
        })
    }
    pub fn push_bytes(&self, value: &[u8]) -> Result<i64, String> {
        self.with_lock(|file| self.push_non_locked_bytes(file, value))
    }
    pub fn push_bytes_size(&self, value: &[u8]) -> Result<(i64, usize), String> {
        self.with_lock(|file| {
            Ok((
                self.push_non_locked_bytes(file, value)?,
                self.size_non_locked(file)?,
            ))
        })
    }
    pub fn pop_exclusive_bytes(&self) -> Result<Option<Vec<u8>>, String> {
        let mut file = self.file.lock().map_err(|e| e.to_string())?;
        self.pop_non_locked_bytes(&mut file, false)
    }
    pub fn push_exclusive_bytes(&self, value: &[u8]) -> Result<i64, String> {
        let mut file = self.file.lock().map_err(|e| e.to_string())?;
        self.push_non_locked_bytes(&mut file, value)
    }

    /// C++ `FileStack::push(int64_t)`.
    pub fn push_i64(&self, value: i64) -> Result<i64, String> {
        self.with_lock(|file| self.push_non_locked_i64(file, value))
    }

    /// C++ `FileStack::push(const string&)`.
    pub fn push_string(&self, value: &str) -> Result<i64, String> {
        self.with_lock(|file| self.push_non_locked_string(file, value))
    }

    /// C++ `FileStack::push(const string&, size_t&)`.
    pub fn push_string_size(&self, value: &str) -> Result<(i64, usize), String> {
        self.with_lock(|file| {
            let written = self.push_non_locked_string(file, value)?;
            Ok((written, self.size_non_locked(file)?))
        })
    }

    /// C++ `FileStack::fetch_add(int64_t)`.
    pub fn fetch_add(&self, n: i64) -> Result<i64, String> {
        self.with_lock(|file| {
            let mut value = -1;
            self.pop_non_locked_i64(file, &mut value)?;
            if value == -1 {
                value = 0;
            }
            self.push_non_locked_i64(file, value + n)?;
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
        self.with_lock(|file| {
            file.set_len(0).map_err(|e| e.to_string())?;
            file.seek(SeekFrom::Start(0))
                .map(|_| ())
                .map_err(|e| e.to_string())
        })
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
            .map_err(|e| self.message_error(b"Error reading file ", format!(": {e}").as_bytes()))
    }

    /// C++ `FileStack::write(const char*, size_t)`.
    pub fn write(&self, buf: &[u8]) -> Result<i64, String> {
        self.file
            .lock()
            .map_err(|e| e.to_string())?
            .write(buf)
            .map(|n| n as i64)
            .map_err(|e| self.message_error(b"Error writing file ", format!(": {e}").as_bytes()))
    }

    /// C++ `FileStack::truncate(size_t)`.
    pub fn truncate(&self, size: usize) -> Result<(), String> {
        self.file
            .lock()
            .map_err(|e| e.to_string())?
            .set_len(size as u64)
            .map_err(|e| e.to_string())
    }

    fn message_error(&self, prefix: &[u8], suffix: &[u8]) -> String {
        let bytes = crate::util::io::path_diagnostic_bytes(&self.file_name).map(|path| {
            let mut bytes = prefix.to_vec();
            bytes.extend_from_slice(&path);
            bytes.extend_from_slice(suffix);
            bytes
        });
        self.errors
            .capture::<()>(Err(match bytes {
                Ok(bytes) => std::io::Error::other(crate::util::io::RawDiagnosticError(bytes)),
                Err(error) => std::io::Error::other(error),
            }))
            .unwrap_err()
    }

    /// C++ string searches operate on bytes, including filenames in queue lines.
    pub fn poll_query(&self, query: &str, sleep_s: f64, max_iter: usize) -> Result<bool, String> {
        self.poll_query_bytes(query.as_bytes(), sleep_s, max_iter)
    }
    pub fn poll_query_bytes(
        &self,
        query: &[u8],
        sleep_s: f64,
        max_iter: usize,
    ) -> Result<bool, String> {
        for _ in 0..max_iter {
            if let Some(value) = self.top_bytes()? {
                if query.is_empty() || value.windows(query.len()).any(|part| part == query) {
                    return Ok(true);
                }
                if value.windows(4).any(|part| part == b"STOP") {
                    return Err(self.message_error(b"STOP on FileStack ", b""));
                }
            }
            std::thread::sleep(Duration::from_secs_f64(sleep_s));
        }
        let mut prefix = b"Could not discover keyword ".to_vec();
        prefix.extend_from_slice(query);
        prefix.extend_from_slice(b" on FileStack ");
        Err(self.message_error(
            &prefix,
            format!(" within {} seconds.", max_iter as f64 * sleep_s).as_bytes(),
        ))
    }

    pub fn poll_size(&self, size: usize, sleep_s: f64, max_iter: usize) -> Result<bool, String> {
        for _ in 0..max_iter {
            if self.size()? == size {
                return Ok(true);
            }
            std::thread::sleep(Duration::from_secs_f64(sleep_s));
        }
        Err(self.message_error(
            format!("Could not detect size {size} of FileStack ").as_bytes(),
            format!(" within {} seconds.", max_iter as f64 * sleep_s).as_bytes(),
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
        let mut file = self.file.lock().map_err(|e| e.to_string())?;
        self.pop_non_locked_string(&mut file, false)
    }

    /// C++ `FileStack::push_exclusive(const string&)`.
    pub fn push_exclusive(&self, value: &str) -> Result<i64, String> {
        let mut file = self.file.lock().map_err(|e| e.to_string())?;
        self.push_non_locked_string(&mut file, value)
    }

    fn push_non_locked_string(&self, file: &mut File, value: &str) -> Result<i64, String> {
        self.push_non_locked_bytes(file, value.as_bytes())
    }

    fn push_non_locked_bytes(&self, file: &mut File, value: &[u8]) -> Result<i64, String> {
        file.seek(SeekFrom::End(0)).map_err(|e| e.to_string())?;
        file.write_all(value).map_err(|e| e.to_string())?;
        let mut written = value.len() as i64;
        if !value.ends_with(b"\n") {
            file.write_all(b"\n").map_err(|e| e.to_string())?;
            written += 1;
        }
        Ok(written)
    }

    /// C++ `FileStack::pop_non_locked(string&, bool, size_t&)`.
    fn pop_non_locked_string(&self, file: &mut File, keep: bool) -> Result<Option<String>, String> {
        // The checked compatibility adapter validates before consuming the row.
        let value = self
            .pop_non_locked_bytes(file, true)?
            .map(|bytes| String::from_utf8(bytes).map_err(|error| error.to_string()))
            .transpose()?;
        if !keep && value.is_some() {
            self.pop_non_locked_bytes(file, false)?;
        }
        Ok(value)
    }

    fn pop_non_locked_bytes(&self, file: &mut File, keep: bool) -> Result<Option<Vec<u8>>, String> {
        let size = file.seek(SeekFrom::End(0)).map_err(|e| e.to_string())? as usize;
        if size == 0 {
            return Ok(None);
        }

        let jump = size.saturating_sub(self.max_line_length);
        file.seek(SeekFrom::Start(jump as u64))
            .map_err(|e| e.to_string())?;
        let mut chunk = Vec::with_capacity(size - jump);
        file.read_to_end(&mut chunk).map_err(|e| e.to_string())?;
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
        let value = chunk[begin..value_end].to_vec();
        if !keep {
            end = size - line_size;
            file.set_len(end as u64).map_err(|e| e.to_string())?;
            file.seek(SeekFrom::Start(end as u64))
                .map_err(|e| e.to_string())?;
        }
        Ok(Some(value))
    }

    /// C++ `FileStack::pop_non_locked(int64_t&)`.
    fn pop_non_locked_i64(&self, file: &mut File, value: &mut i64) -> Result<i64, String> {
        match self.pop_non_locked_string(file, false)? {
            Some(text) => {
                *value = text.parse::<i64>().map_err(|e| e.to_string())?;
            }
            None => *value = -1,
        }
        Ok(*value)
    }

    /// C++ `FileStack::push_non_locked(int64_t)`.
    fn push_non_locked_i64(&self, file: &mut File, value: i64) -> Result<i64, String> {
        self.push_non_locked_string(file, &value.to_string())
    }

    fn size_non_locked(&self, file: &mut File) -> Result<usize, String> {
        file.seek(SeekFrom::Start(0)).map_err(|e| e.to_string())?;
        let mut data = Vec::new();
        file.read_to_end(&mut data).map_err(|e| e.to_string())?;
        Ok(data.iter().filter(|&&byte| byte == b'\n').count())
    }

    /// C++ `lock`/`unlock`, expressed as a scoped operation so error paths
    /// cannot leave either lock held.
    fn with_lock<T>(
        &self,
        operation: impl FnOnce(&mut File) -> Result<T, String>,
    ) -> Result<T, String> {
        let _thread_lock = self.mtx.lock().map_err(|e| e.to_string())?;
        let mut file = self.file.lock().map_err(|e| e.to_string())?;
        lock_file(&file).map_err(|e| {
            self.message_error(b"could not put lock on file ", format!(": {e}").as_bytes())
        })?;
        let result = operation(&mut file);
        let unlock = unlock_file(&file).map_err(|e| {
            self.message_error(b"could not unlock file ", format!(": {e}").as_bytes())
        });
        match (result, unlock) {
            (Err(error), _) => Err(error),
            (Ok(_), Err(error)) => Err(error),
            (Ok(value), Ok(())) => Ok(value),
        }
    }
}

/// Exclusive advisory lock on the whole file (equivalent of `File::lock`,
/// which is only stable since Rust 1.89).
#[cfg(unix)]
fn lock_file(file: &File) -> std::io::Result<()> {
    use std::os::unix::io::AsRawFd;
    loop {
        if unsafe { libc::flock(file.as_raw_fd(), libc::LOCK_EX) } == 0 {
            return Ok(());
        }
        let err = std::io::Error::last_os_error();
        if err.kind() != std::io::ErrorKind::Interrupted {
            return Err(err);
        }
    }
}

#[cfg(unix)]
fn unlock_file(file: &File) -> std::io::Result<()> {
    use std::os::unix::io::AsRawFd;
    if unsafe { libc::flock(file.as_raw_fd(), libc::LOCK_UN) } == 0 {
        Ok(())
    } else {
        Err(std::io::Error::last_os_error())
    }
}

#[cfg(windows)]
fn lock_file(file: &File) -> std::io::Result<()> {
    use std::os::windows::io::AsRawHandle;
    use windows_sys::Win32::Storage::FileSystem::{LockFileEx, LOCKFILE_EXCLUSIVE_LOCK};
    use windows_sys::Win32::System::IO::OVERLAPPED;
    let mut overlapped: OVERLAPPED = unsafe { std::mem::zeroed() };
    let ok = unsafe {
        LockFileEx(
            file.as_raw_handle() as _,
            LOCKFILE_EXCLUSIVE_LOCK,
            0,
            u32::MAX,
            u32::MAX,
            &mut overlapped,
        )
    };
    if ok != 0 {
        Ok(())
    } else {
        Err(std::io::Error::last_os_error())
    }
}

#[cfg(windows)]
fn unlock_file(file: &File) -> std::io::Result<()> {
    use std::os::windows::io::AsRawHandle;
    use windows_sys::Win32::Storage::FileSystem::UnlockFileEx;
    use windows_sys::Win32::System::IO::OVERLAPPED;
    let mut overlapped: OVERLAPPED = unsafe { std::mem::zeroed() };
    let ok = unsafe {
        UnlockFileEx(
            file.as_raw_handle() as _,
            0,
            u32::MAX,
            u32::MAX,
            &mut overlapped,
        )
    };
    if ok != 0 {
        Ok(())
    } else {
        Err(std::io::Error::last_os_error())
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn raw_stack_roundtrip_keeps_non_utf8_lines() {
        let path = std::env::temp_dir().join(format!("diamond-raw-stack-{}", std::process::id()));
        let _ = std::fs::remove_file(&path);
        let stack = super::FileStack::try_new(&path).unwrap();
        stack.push_bytes(b"volume_\xff\t2").unwrap();
        stack.push_bytes(b"second_\x80").unwrap();
        assert_eq!(stack.top_bytes().unwrap(), Some(b"second_\x80".to_vec()));
        stack.remove_bytes(b"second_\x80").unwrap();
        assert_eq!(stack.pop_bytes().unwrap(), Some(b"volume_\xff\t2".to_vec()));
        assert_eq!(stack.pop_bytes().unwrap(), None);
        drop(stack);
        std::fs::remove_file(path).unwrap();
    }

    use super::*;
    use std::sync::Arc;

    fn path(tag: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-filestack-mirror-{}-{tag}.idx",
            std::process::id()
        ))
    }

    #[cfg(windows)]
    #[test]
    fn source_stack_allows_peer_access_but_denies_unlink_until_close() {
        let path = path("sharing");
        let first = FileStack::try_new(&path).unwrap();
        let second = FileStack::try_new(&path).unwrap();
        assert_eq!(
            std::fs::remove_file(&path).unwrap_err().raw_os_error(),
            Some(32)
        );
        drop(first);
        assert_eq!(
            std::fs::remove_file(&path).unwrap_err().raw_os_error(),
            Some(32)
        );
        drop(second);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn remove_preserves_empty_records_like_cpp_split_and_rewrite() {
        let path = path("empty-records");
        let _ = std::fs::remove_file(&path);
        std::fs::write(&path, b"alpha\n\nbeta\n").unwrap();
        let stack = FileStack::new(&path);
        stack.remove("alpha").unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"\nbeta\n");
        drop(stack);
        std::fs::remove_file(path).unwrap();
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
        drop((first, second));
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn maximum_line_length_controls_cpp_tail_window() {
        let path = path("tail-window");
        let _ = std::fs::remove_file(&path);
        let stack = FileStack::with_max_line_length(&path, 8);
        stack.push_string("123456789").unwrap();
        assert_eq!(stack.pop_string().unwrap().as_deref(), Some("3456789"));
        assert_eq!(std::fs::read(&path).unwrap(), b"12");
        drop(stack);
        std::fs::remove_file(path).unwrap();
    }
}
