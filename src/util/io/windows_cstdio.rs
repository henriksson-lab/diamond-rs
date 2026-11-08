//! Windows CRT stream owner; FILE remains the opaque maintained libc type.
//!
//! Extra CRT declarations are generated from maintained headers, with no
//! runtime DLL override: all functions use the CRT selected by Rust/libc.
//!
//! Named Rust Unicode paths use equivalent wide CRT entry points. This
//! prevents UTF-8 path bytes being interpreted as the Windows ANSI codepage.
//! Paths outside the old Source narrow argv codepage round-trip domain are a
//! Rust Unicode extension. Argument/diagnostic conversion uses actual CP_ACP
//! flags0; Windows Source runtime/argv equivalence still requires host proof.
use super::cstdio_file::{StdioError, StdioOperation};
use super::windows_crt_bindings as crt;
use super::{FileOpenError, IoError, IoResult};
use std::ffi::{CStr, CString, OsStr, OsString};
use std::fs::File;
use std::io::SeekFrom;
use std::os::windows::ffi::{OsStrExt, OsStringExt};
use std::os::windows::io::{FromRawHandle, IntoRawHandle};
use std::path::{Path, PathBuf};
use std::ptr::NonNull;
use windows_sys::Win32::Foundation::{DuplicateHandle, DUPLICATE_SAME_ACCESS};
use windows_sys::Win32::Globalization::{MultiByteToWideChar, WideCharToMultiByte, CP_ACP};
use windows_sys::Win32::Storage::FileSystem::{SetFilePointerEx, FILE_CURRENT};
use windows_sys::Win32::System::Threading::GetCurrentProcess;

fn errno() -> libc::c_int {
    #[cfg(target_env = "msvc")]
    {
        let mut value = 0;
        // A valid pointer cannot trigger the invalid-parameter branch.
        if unsafe { crt::_get_errno(&mut value) } == 0 {
            return value;
        }
    }
    // GNU uses its matching CRT's thread-local errno, not GetLastError.
    unsafe { *crt::_errno() }
}

fn reason(code: libc::c_int) -> Vec<u8> {
    let pointer = unsafe { libc::strerror(code) };
    if pointer.is_null() {
        return format!("CRT error {code}").into_bytes();
    }
    unsafe { CStr::from_ptr(pointer) }.to_bytes().to_vec()
}

fn open_error(name: &Path, code: libc::c_int) -> IoError {
    // FileOpenError::new expects a Win32 OS error on Windows. This operation
    // instead failed inside the CRT, so retain its actual saved errno/reason.
    IoError::FileOpenDetails(FileOpenError {
        path: name.to_path_buf(),
        reason: reason(code),
        raw_os_error: Some(code),
    })
}

fn mode_string(mode: &str) -> IoResult<CString> {
    if !matches!(mode, "rb" | "wb" | "r+b" | "w+b") {
        return Err(IoError::Other("Invalid fopen mode.".into()));
    }
    Ok(CString::new(mode).expect("validated stdio mode"))
}

/// Source narrow argv/diagnostics use the process ANSI codepage. This is
/// deliberately separate from the original UTF16 Path retained by wide opens.
/// Flags0/default-character behavior follows the actual Win32 ACP conversion;
/// non-ACP Unicode path support is a qualified Native extension, not byte parity.
pub(crate) fn source_os_bytes(value: &OsStr) -> IoResult<Vec<u8>> {
    let wide: Vec<u16> = value.encode_wide().collect();
    if wide.is_empty() {
        return Ok(Vec::new());
    }
    let length = i32::try_from(wide.len())
        .map_err(|_| IoError::Other("OS string exceeds Windows conversion length.".into()))?;
    let count = unsafe {
        WideCharToMultiByte(
            CP_ACP,
            0,
            wide.as_ptr(),
            length,
            std::ptr::null_mut(),
            0,
            std::ptr::null(),
            std::ptr::null_mut(),
        )
    };
    if count == 0 {
        return Err(conversion_error("WideCharToMultiByte"));
    }
    let mut bytes = vec![0; count as usize];
    let written = unsafe {
        WideCharToMultiByte(
            CP_ACP,
            0,
            wide.as_ptr(),
            length,
            bytes.as_mut_ptr(),
            count,
            std::ptr::null(),
            std::ptr::null_mut(),
        )
    };
    if written == 0 {
        return Err(conversion_error("WideCharToMultiByte"));
    }
    bytes.truncate(written as usize);
    Ok(bytes)
}
pub(crate) fn source_bytes_os(bytes: &[u8]) -> IoResult<OsString> {
    if bytes.is_empty() {
        return Ok(OsString::new());
    }
    let length = i32::try_from(bytes.len()).map_err(|_| {
        IoError::Other("Source byte string exceeds Windows conversion length.".into())
    })?;
    let count =
        unsafe { MultiByteToWideChar(CP_ACP, 0, bytes.as_ptr(), length, std::ptr::null_mut(), 0) };
    if count == 0 {
        return Err(conversion_error("MultiByteToWideChar"));
    }
    let mut wide = vec![0; count as usize];
    let written =
        unsafe { MultiByteToWideChar(CP_ACP, 0, bytes.as_ptr(), length, wide.as_mut_ptr(), count) };
    if written == 0 {
        return Err(conversion_error("MultiByteToWideChar"));
    }
    wide.truncate(written as usize);
    Ok(OsString::from_wide(&wide))
}
fn conversion_error(operation: &str) -> IoError {
    // These are Win32 conversion failures, not CRT errno/file-open failures.
    let error = std::io::Error::last_os_error();
    IoError::Other(format!("{operation} failed: {error}"))
}
fn wide_path(name: &Path) -> IoResult<Vec<libc::wchar_t>> {
    let mut path: Vec<_> = super::source_c_path(name)
        .as_os_str()
        .encode_wide()
        .collect();
    path.push(0);
    Ok(path)
}

/// One uniquely owned FILE, except MSVC's borrowed process stdout. All stream
/// operations require mutable access; no FILE pointer escapes the wrapper.
#[derive(Debug)]
pub(crate) struct OwnedCstdio {
    stream: Option<NonNull<libc::FILE>>,
    name: PathBuf,
    close_stream: bool,
    file_view: Option<File>,
    view_dirty: bool,
}

// SAFETY: moving the owner transfers exclusive access to FILE. Immutable
// access exposes neither FILE nor mutating operations; CRT stdio also locks
// its streams. Shared consumers use SharedCstdio's mutex.
unsafe impl Send for OwnedCstdio {}
unsafe impl Sync for OwnedCstdio {}

impl OwnedCstdio {
    /// Source's distinct File facade uses fopen directly (including r/a
    /// modes), rather than FileSink's GNU open+fdopen branch.
    pub(crate) fn open_file_role_path(name: &Path, mode: &str) -> IoResult<Self> {
        let path = wide_path(name)?;
        let mut mode: Vec<u16> = mode.encode_utf16().collect();
        if mode.contains(&0) {
            return Err(open_error(name, libc::EINVAL));
        }
        mode.push(0);
        let pointer = unsafe { crt::_wfopen(path.as_ptr(), mode.as_ptr()) };
        let stream = NonNull::new(pointer).ok_or_else(|| open_error(name, errno()))?;
        Self::from_stream(name, stream, true, 0)
    }
    /// One actual fread; Source File::read inspects only the returned count.
    /// Save errno immediately, without strerror or feof/ferror interpretation.
    pub(crate) fn read_file_role(&mut self, bytes: &mut [u8]) -> IoResult<(usize, i32)> {
        self.sync_view()?;
        let pointer = self.pointer()?;
        let count =
            unsafe { libc::fread(bytes.as_mut_ptr().cast(), 1, bytes.len(), pointer.as_ptr()) };
        let saved = errno();
        Ok((count, saved))
    }
    pub(crate) fn set_file_role_buffer(&mut self, size: usize) -> IoResult<()> {
        let pointer = self.pointer()?;
        if unsafe { libc::setvbuf(pointer.as_ptr(), std::ptr::null_mut(), libc::_IOFBF, size) } != 0
        {
            let saved = errno();
            let mut message = b"Error setting buffer size for temporary file ".to_vec();
            message.extend_from_slice(&super::path_diagnostic_bytes(&self.name)?);
            message.extend_from_slice(b". ");
            message.extend_from_slice(&reason(saved));
            return Err(IoError::RawDetails(super::RawDiagnosticError(message)));
        }
        Ok(())
    }
    /// File::close ignores fclose's return; do not materialize FileSink
    /// diagnostics/strerror or alter its stale-errno behavior before ignoring it.
    pub(crate) fn close_file_role(&mut self) -> IoResult<()> {
        self.file_view.take();
        self.view_dirty = false;
        if self.close_stream {
            if let Some(pointer) = self.stream.take() {
                unsafe { libc::fclose(pointer.as_ptr()) };
            }
        }
        Ok(())
    }

    pub(crate) fn open(name: &str, mode: &str, buffer_size: usize) -> IoResult<Self> {
        Self::open_path(Path::new(name), mode, buffer_size)
    }
    pub(crate) fn open_path(name: &Path, mode: &str, buffer_size: usize) -> IoResult<Self> {
        let mode_c = mode_string(mode)?;
        let path = wide_path(name)?;
        #[cfg(target_env = "msvc")]
        {
            let _ = mode_c;
            // _wfopen changes only filename encoding, not stream translation,
            // sharing, ownership or the Source fopen constructor boundary.
            let mode_wide: Vec<libc::wchar_t> = mode.encode_utf16().chain(Some(0)).collect();
            let pointer = unsafe { crt::_wfopen(path.as_ptr(), mode_wide.as_ptr()) };
            let stream = NonNull::new(pointer).ok_or_else(|| open_error(name, errno()))?;
            Self::from_stream(name, stream, true, buffer_size)
        }
        #[cfg(not(target_env = "msvc"))]
        {
            // Source's non-MSVC branch uses open followed by fdopen, without
            // adding a separate O_BINARY flag to its posix_flags result.
            let flags = match mode {
                "rb" => libc::O_RDONLY,
                "wb" => libc::O_WRONLY | libc::O_CREAT | libc::O_TRUNC,
                "r+b" => libc::O_RDWR,
                "w+b" => libc::O_RDWR | libc::O_CREAT | libc::O_TRUNC,
                _ => unreachable!(),
            };
            let descriptor = unsafe { libc::wopen(path.as_ptr(), flags, 0o666 as libc::c_int) };
            if descriptor < 0 {
                return Err(open_error(name, errno()));
            }
            Self::from_descriptor(name, descriptor, &mode_c, buffer_size)
        }
    }

    pub(crate) fn from_file(
        name: &str,
        file: File,
        mode: &str,
        buffer_size: usize,
    ) -> IoResult<Self> {
        Self::from_file_path(Path::new(name), file, mode, buffer_size)
    }
    pub(crate) fn from_file_path(
        name: &Path,
        file: File,
        mode: &str,
        buffer_size: usize,
    ) -> IoResult<Self> {
        let mode_c = mode_string(mode)?;
        let flags = if mode.contains('+') {
            libc::O_RDWR
        } else if mode.starts_with('r') {
            libc::O_RDONLY
        } else {
            libc::O_WRONLY
        } | libc::O_BINARY;
        let handle = file.into_raw_handle();
        // The CRT takes ownership of the actual HANDLE on success. No
        // pathname is reopened and no duplicate is made at this boundary.
        let descriptor = unsafe { libc::open_osfhandle(handle as libc::intptr_t, flags) };
        if descriptor < 0 {
            let saved = errno();
            // Failure did not transfer ownership to the CRT.
            drop(unsafe { File::from_raw_handle(handle) });
            return Err(open_error(name, saved));
        }
        Self::from_descriptor(name, descriptor, &mode_c, buffer_size)
    }

    fn from_descriptor(
        name: &Path,
        descriptor: libc::c_int,
        mode: &CString,
        buffer_size: usize,
    ) -> IoResult<Self> {
        let pointer = unsafe { libc::fdopen(descriptor, mode.as_ptr()) };
        let Some(stream) = NonNull::new(pointer) else {
            let saved = errno();
            // fdopen failed: the descriptor (and any adopted HANDLE) remains
            // ours. Closing it exactly once also closes that HANDLE.
            unsafe { libc::close(descriptor) };
            return Err(open_error(name, saved));
        };
        Self::from_stream(name, stream, true, buffer_size)
    }

    fn from_stream(
        name: &Path,
        stream: NonNull<libc::FILE>,
        close_stream: bool,
        buffer_size: usize,
    ) -> IoResult<Self> {
        let result = Self {
            stream: Some(stream),
            name: name.to_path_buf(),
            close_stream,
            file_view: None,
            view_dirty: false,
        };
        if buffer_size != 0
            && unsafe {
                libc::setvbuf(
                    stream.as_ptr(),
                    std::ptr::null_mut(),
                    libc::_IOFBF,
                    buffer_size,
                )
            } != 0
        {
            return Err(IoError::Other("Error calling setvbuf.".into()));
        }
        Ok(result)
    }

    pub(crate) fn standard_output(name: &str, mode: &str, buffer_size: usize) -> IoResult<Self> {
        let name = Path::new(name);
        let mode = mode_string(mode)?;
        #[cfg(target_env = "msvc")]
        {
            let _ = mode;
            let stream = NonNull::new(unsafe { crt::__acrt_iob_func(1) })
                .ok_or_else(|| open_error(name, errno()))?;
            // Source FileSink::close explicitly skips the actual stdout FILE.
            Self::from_stream(name, stream, false, buffer_size)
        }
        #[cfg(not(target_env = "msvc"))]
        Self::from_descriptor(name, 1, &mode, buffer_size)
    }

    pub(crate) fn standard_input(name: &str) -> IoResult<Self> {
        let name = Path::new(name);
        #[cfg(target_env = "msvc")]
        {
            let stream = NonNull::new(unsafe { crt::__acrt_iob_func(0) })
                .ok_or_else(|| open_error(name, errno()))?;
            // FileSource::close closes stdin too, unlike FileSink's stdout.
            Self::from_stream(name, stream, true, 0)
        }
        #[cfg(not(target_env = "msvc"))]
        Self::from_descriptor(name, 0, &mode_string("rb")?, 0)
    }

    pub(crate) fn source_seekable(name: &str) -> IoResult<bool> {
        Self::source_seekable_path(Path::new(name))
    }
    pub(crate) fn source_seekable_path(name: &Path) -> IoResult<bool> {
        let path = wide_path(name)?;
        let mut metadata = std::mem::MaybeUninit::<libc::stat>::uninit();
        if unsafe { libc::wstat(path.as_ptr(), metadata.as_mut_ptr()) } < 0 {
            return Err(open_error(name, errno()));
        }
        let metadata = unsafe { metadata.assume_init() };
        Ok((metadata.st_mode as libc::c_int & libc::S_IFMT) == libc::S_IFREG)
    }

    fn pointer(&self) -> IoResult<NonNull<libc::FILE>> {
        self.stream
            .ok_or_else(|| IoError::path_message(b"Error writing file ", &self.name))
    }

    fn failure(&self, operation: StdioOperation, code: libc::c_int) -> IoError {
        let (prefix, named): (&[u8], bool) = match operation {
            StdioOperation::Write => (b"Error writing file ", true),
            StdioOperation::Read => (b"Error reading file ", true),
            StdioOperation::Close => (b"Error closing file ", true),
            StdioOperation::Seek => (b"Error calling fseek.", false),
            StdioOperation::Tell if cfg!(target_env = "msvc") => {
                (b"Error executing ftell on stream ", true)
            }
            StdioOperation::Tell => (b"Error calling ftell.", false),
        };
        let mut message = prefix.to_vec();
        if named {
            match super::path_diagnostic_bytes(&self.name) {
                Ok(path) => message.extend_from_slice(&path),
                Err(error) => return error,
            }
        }
        IoError::StdioDetails(StdioError {
            operation,
            path: self.name.clone(),
            raw_os_error: Some(code),
            perror: (!(cfg!(target_env = "msvc") && operation == StdioOperation::Tell))
                .then(|| reason(code)),
            message,
        })
    }

    fn os_handle(&self) -> IoResult<*mut std::ffi::c_void> {
        let descriptor = unsafe { libc::fileno(self.pointer()?.as_ptr()) };
        if descriptor < 0 {
            return Err(open_error(&self.name, errno()));
        }
        let handle = unsafe { libc::get_osfhandle(descriptor) };
        if handle == -1 {
            return Err(open_error(&self.name, errno()));
        }
        Ok(handle as *mut std::ffi::c_void)
    }

    fn sync_view(&mut self) -> IoResult<()> {
        if !self.view_dirty {
            return Ok(());
        }
        let mut position = 0;
        if unsafe { SetFilePointerEx(self.os_handle()?, 0, &mut position, FILE_CURRENT) } == 0 {
            // This is the Rust legacy-handle interop boundary, not a failed CRT
            // operation: retain its real Win32 error rather than sampling errno.
            return Err(IoError::file_open(
                &self.name,
                std::io::Error::last_os_error(),
            ));
        }
        self.seek_raw(position, SeekFrom::Start(0), false)?;
        self.view_dirty = false;
        Ok(())
    }

    /// Exception-disabled Source ofstream transport. The caller owns its
    /// sticky failbit; short writes do not materialize errno diagnostics.
    pub(crate) fn write_text_output(&mut self, bytes: &[u8]) -> IoResult<bool> {
        self.sync_view()?;
        let pointer = self.pointer()?;
        Ok(
            unsafe { libc::fwrite(bytes.as_ptr().cast(), 1, bytes.len(), pointer.as_ptr()) }
                == bytes.len(),
        )
    }
    /// Flush exactly this text stream at a Source endl boundary.
    pub(crate) fn flush_text_output(&mut self) -> IoResult<bool> {
        self.sync_view()?;
        let pointer = self.pointer()?;
        Ok(unsafe { libc::fflush(pointer.as_ptr()) } == 0)
    }

    pub(crate) fn write_bytes(&mut self, bytes: &[u8]) -> IoResult<()> {
        self.sync_view()?;
        let pointer = self.pointer()?;
        if unsafe { libc::fwrite(bytes.as_ptr().cast(), 1, bytes.len(), pointer.as_ptr()) }
            != bytes.len()
        {
            return Err(self.failure(StdioOperation::Write, errno()));
        }
        Ok(())
    }

    pub(crate) fn read_bytes(&mut self, bytes: &mut [u8]) -> IoResult<usize> {
        self.sync_view()?;
        let pointer = self
            .stream
            .ok_or_else(|| IoError::path_message(b"Error reading file ", &self.name))?;
        let count =
            unsafe { libc::fread(bytes.as_mut_ptr().cast(), 1, bytes.len(), pointer.as_ptr()) };
        let saved = errno();
        if count != bytes.len() && unsafe { libc::feof(pointer.as_ptr()) } == 0 {
            return Err(self.failure(StdioOperation::Read, saved));
        }
        Ok(count)
    }

    pub(crate) fn eof(&mut self) -> IoResult<bool> {
        let pointer = self
            .stream
            .ok_or_else(|| IoError::path_message(b"Error reading file ", &self.name))?;
        Ok(unsafe { libc::feof(pointer.as_ptr()) } != 0)
    }

    fn seek_raw(&mut self, position: i64, origin: SeekFrom, input: bool) -> IoResult<()> {
        #[cfg(target_env = "msvc")]
        let _ = input;
        let pointer = self.pointer()?;
        let whence = match origin {
            SeekFrom::Start(_) => libc::SEEK_SET,
            SeekFrom::Current(_) => libc::SEEK_CUR,
            SeekFrom::End(_) => libc::SEEK_END,
        };
        #[cfg(target_env = "msvc")]
        let failed = unsafe { crt::_fseeki64(pointer.as_ptr(), position, whence) } != 0;
        #[cfg(not(target_env = "msvc"))]
        let failed = {
            let result = unsafe { libc::fseek(pointer.as_ptr(), position as libc::c_long, whence) };
            if input {
                result < 0
            } else {
                result != 0
            }
        };
        if failed {
            return Err(self.failure(StdioOperation::Seek, errno()));
        }
        self.view_dirty = false;
        Ok(())
    }

    pub(crate) fn seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        if matches!(origin, SeekFrom::Current(_)) {
            self.sync_view()?;
        }
        self.seek_raw(position, origin, false)
    }

    pub(crate) fn position(&mut self) -> IoResult<i64> {
        self.sync_view()?;
        let pointer = self.pointer()?;
        #[cfg(target_env = "msvc")]
        let position = unsafe { crt::_ftelli64(pointer.as_ptr()) };
        #[cfg(not(target_env = "msvc"))]
        let position = unsafe { libc::ftell(pointer.as_ptr()) } as i64;
        if position < 0 {
            return Err(self.failure(StdioOperation::Tell, errno()));
        }
        Ok(position)
    }

    fn input_error(&self, error: IoError, operation: StdioOperation) -> IoError {
        match error {
            IoError::StdioDetails(mut error) => {
                let prefix: &[u8] = match operation {
                    StdioOperation::Seek => b"Error executing seek on file ",
                    StdioOperation::Tell => b"Error executing ftell on stream ",
                    _ => return IoError::StdioDetails(error),
                };
                let path = match super::path_diagnostic_bytes(&self.name) {
                    Ok(bytes) => bytes,
                    Err(error) => return error,
                };
                error.message = prefix.to_vec();
                error.message.extend_from_slice(&path);
                IoError::StdioDetails(error)
            }
            error => error,
        }
    }

    pub(crate) fn input_seek(&mut self, position: i64, origin: SeekFrom) -> IoResult<()> {
        if matches!(origin, SeekFrom::Current(_)) {
            self.sync_view()?;
        }
        self.seek_raw(position, origin, true)
            .map_err(|error| self.input_error(error, StdioOperation::Seek))
    }

    pub(crate) fn input_position(&mut self) -> IoResult<i64> {
        self.position()
            .map_err(|error| self.input_error(error, StdioOperation::Tell))
    }

    pub(crate) fn rewind(&mut self) -> IoResult<()> {
        unsafe { libc::rewind(self.pointer()?.as_ptr()) };
        self.view_dirty = false;
        Ok(())
    }

    pub(crate) fn close(&mut self) -> IoResult<()> {
        self.file_view.take();
        self.view_dirty = false;
        if !self.close_stream {
            return Ok(());
        }
        if let Some(pointer) = self.stream.take() {
            if unsafe { libc::fclose(pointer.as_ptr()) } != 0 {
                return Err(self.failure(StdioOperation::Close, errno()));
            }
        }
        Ok(())
    }

    pub(crate) fn clone_file(&mut self) -> IoResult<File> {
        self.sync_view()?;
        let process = unsafe { GetCurrentProcess() };
        let mut duplicate = std::ptr::null_mut();
        if unsafe {
            DuplicateHandle(
                process,
                self.os_handle()?,
                process,
                &mut duplicate,
                0,
                0,
                DUPLICATE_SAME_ACCESS,
            )
        } == 0
        {
            return Err(IoError::file_open(
                &self.name,
                std::io::Error::last_os_error(),
            ));
        }
        self.view_dirty = true;
        Ok(unsafe { File::from_raw_handle(duplicate) })
    }

    pub(crate) fn prepare_legacy(&mut self) -> IoResult<()> {
        self.sync_view()?;
        if unsafe { libc::fflush(self.pointer()?.as_ptr()) } != 0 {
            return Err(self.failure(StdioOperation::Write, errno()));
        }
        self.view_dirty = true;
        Ok(())
    }

    pub(crate) fn legacy_file(&mut self) -> IoResult<&mut File> {
        self.prepare_legacy()?;
        if self.file_view.is_none() {
            self.file_view = Some(self.clone_file()?);
        }
        self.view_dirty = true;
        Ok(self.file_view.as_mut().unwrap())
    }
}

impl Drop for OwnedCstdio {
    fn drop(&mut self) {
        self.file_view.take();
        if self.close_stream {
            if let Some(pointer) = self.stream.take() {
                // Explicit close owns diagnostics/fatal policy; abandonment is
                // resource cleanup and does not retry a failed fclose.
                unsafe { libc::fclose(pointer.as_ptr()) };
            }
        }
    }
}

pub(crate) fn source_external_stdio_limit() {
    // Source ignores this return value at External::external's actual boundary.
    unsafe { crt::_setmaxstdio(8192) };
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Seek, Write};

    fn path(label: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-windows-cstdio-{label}-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ))
    }

    #[test]
    fn source_ansi_codec_roundtrips_ascii_and_wide_path_preserves_non_unicode_units() {
        let input = OsStr::new("source ASCII\0tail");
        let bytes = source_os_bytes(input).unwrap();
        assert_eq!(bytes, b"source ASCII\0tail");
        assert_eq!(source_bytes_os(&bytes).unwrap().as_os_str(), input);
        assert!(source_os_bytes(OsStr::new("")).unwrap().is_empty());
        assert_eq!(source_bytes_os(b"").unwrap(), OsString::new());
        let raw = OsString::from_wide(&[b'a' as u16, 0xd800, b'z' as u16]);
        assert_eq!(
            wide_path(Path::new(&raw)).unwrap(),
            [b'a' as u16, 0xd800, b'z' as u16, 0]
        );
    }

    #[test]
    fn source_file_role_modes_retain_the_actual_path_and_unchecked_close() {
        let path = path("file-role-café-λ");
        let mut writer = OwnedCstdio::open_file_role_path(&path, "wb").unwrap();
        writer.write_bytes(b"first").unwrap();
        writer.close_file_role().unwrap();
        let mut append = OwnedCstdio::open_file_role_path(&path, "ab").unwrap();
        append.write_bytes(b" second").unwrap();
        append.close_file_role().unwrap();
        let mut reader = OwnedCstdio::open_file_role_path(&path, "rb").unwrap();
        let mut bytes = [0; 12];
        let (count, _) = reader.read_file_role(&mut bytes).unwrap();
        assert_eq!(count, bytes.len());
        assert_eq!(&bytes, b"first second");
        reader.close_file_role().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), bytes);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn transferred_handle_and_legacy_duplicate_keep_one_stream_position() {
        let path = path("transfer");
        let file = std::fs::OpenOptions::new()
            .create_new(true)
            .read(true)
            .write(true)
            .open(&path)
            .unwrap();
        let mut owner = OwnedCstdio::from_file(path.to_str().unwrap(), file, "w+b", 4096).unwrap();
        owner.write_bytes(b"abcdef").unwrap();
        assert_eq!(owner.position().unwrap(), 6);
        let view = owner.legacy_file().unwrap();
        view.seek(SeekFrom::Start(2)).unwrap();
        view.write_all(b"X").unwrap();
        owner.write_bytes(b"Y").unwrap();
        assert_eq!(owner.position().unwrap(), 4);
        owner.rewind().unwrap();
        let mut bytes = [0; 6];
        assert_eq!(owner.read_bytes(&mut bytes).unwrap(), 6);
        assert_eq!(&bytes, b"abXYef");
        assert!(!owner.eof().unwrap());
        assert_eq!(owner.read_bytes(&mut [0]).unwrap(), 0);
        assert!(owner.eof().unwrap());
        owner.close().unwrap();
        owner.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"abXYef");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn named_unicode_paths_retain_the_same_real_file_for_read_write_and_stat() {
        // Non-ACP names are an intentional Unicode library extension. This
        // also exercises names representable in common ACPs (é), for which
        // UTF-8 bytes passed to narrow fopen would identify the wrong path.
        let path = path("café-λ-漢字-🧬");
        let name = path.to_str().unwrap();
        let mut output = OwnedCstdio::open(name, "wb", 4096).unwrap();
        output.write_bytes(b"same physical file\n").unwrap();
        output.close().unwrap();
        assert!(OwnedCstdio::source_seekable(name).unwrap());
        assert_eq!(std::fs::read(&path).unwrap(), b"same physical file\n");
        let mut input = OwnedCstdio::open(name, "rb", 4096).unwrap();
        let mut bytes = [0; 19];
        assert_eq!(input.read_bytes(&mut bytes).unwrap(), bytes.len());
        assert_eq!(&bytes, b"same physical file\n");
        input.close().unwrap();
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn named_path_embedded_nul_preserves_invalid_argument_error() {
        let IoError::FileOpenDetails(error) =
            OwnedCstdio::open("invalid\0path", "rb", 0).unwrap_err()
        else {
            panic!("embedded NUL must remain a concrete open error");
        };
        assert_eq!(error.raw_os_error, Some(libc::EINVAL));
    }

    #[test]
    fn missing_input_preserves_the_actual_crt_errno() {
        let path = path("missing");
        let IoError::FileOpenDetails(error) =
            OwnedCstdio::open(path.to_str().unwrap(), "rb", 0).unwrap_err()
        else {
            panic!("fopen must retain its CRT open error");
        };
        assert_eq!(error.raw_os_error, Some(libc::ENOENT));
        assert_eq!(error.path, path);
        assert_eq!(error.reason, reason(libc::ENOENT));
    }

    #[cfg(target_env = "msvc")]
    #[test]
    fn msvc_seek_and_tell_preserve_positions_above_signed32_without_allocating_data() {
        let path = path("wide-position");
        let mut owner = OwnedCstdio::open(path.to_str().unwrap(), "w+b", 0).unwrap();
        let position = (3_i64 << 30) + 17;
        owner.seek(position, SeekFrom::Start(0)).unwrap();
        assert_eq!(owner.position().unwrap(), position);
        owner.rewind().unwrap();
        assert_eq!(owner.position().unwrap(), 0);
        owner.close().unwrap();
        assert_eq!(std::fs::metadata(&path).unwrap().len(), 0);
        std::fs::remove_file(path).unwrap();
    }
}
