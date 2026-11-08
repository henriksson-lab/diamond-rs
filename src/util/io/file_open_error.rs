//! Source FileOpenException metadata captured at the failing OS operation.
//! The C++ FileSource/FileSink constructors call perror before throwing; Rust
//! retains that diagnostic for the CLI boundary instead of printing in a
//! reusable library opener or inspecting a later, unrelated errno value.

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileOpenError {
    pub path: std::path::PathBuf,
    pub reason: Vec<u8>,
    pub raw_os_error: Option<i32>,
}

impl FileOpenError {
    pub fn new(path: impl AsRef<std::path::Path>, error: std::io::Error) -> Self {
        let raw_os_error = error.raw_os_error();
        #[cfg(windows)]
        let errno = match error.kind() {
            std::io::ErrorKind::NotFound => Some(libc::ENOENT),
            std::io::ErrorKind::PermissionDenied => Some(libc::EACCES),
            std::io::ErrorKind::AlreadyExists => Some(libc::EEXIST),
            std::io::ErrorKind::InvalidInput => Some(libc::EINVAL),
            _ => None,
        };
        #[cfg(not(windows))]
        let errno = raw_os_error;
        let reason = errno
            .and_then(|code| {
                // strerror describes the saved error code; it does not sample errno.
                let ptr = unsafe { libc::strerror(code) };
                (!ptr.is_null())
                    .then(|| unsafe { std::ffi::CStr::from_ptr(ptr) }.to_bytes().to_vec())
            })
            .unwrap_or_else(|| error.to_string().into_bytes());
        Self {
            path: path.as_ref().to_path_buf(),
            reason,
            raw_os_error,
        }
    }

    pub fn try_message_bytes(&self) -> super::IoResult<Vec<u8>> {
        let mut message = b"Error opening file ".to_vec();
        message.extend_from_slice(&path_diagnostic_bytes(&self.path)?);
        Ok(message)
    }
    pub fn try_perror_bytes(&self) -> super::IoResult<Vec<u8>> {
        let mut message = b"\n".to_vec();
        let label = self.try_message_bytes()?;
        message.extend_from_slice(label.split(|&byte| byte == 0).next().unwrap_or_default());
        message.extend_from_slice(b": ");
        message.extend_from_slice(&self.reason);
        message.push(b'\n');
        Ok(message)
    }

    pub fn into_io(path: impl AsRef<std::path::Path>, error: std::io::Error) -> std::io::Error {
        let kind = error.kind();
        std::io::Error::new(kind, Self::new(path, error))
    }
}

impl std::fmt::Display for FileOpenError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Keep the translated typed I/O error's historical display contract.
        write!(f, "Error opening file {}", self.path.display())
    }
}
impl std::error::Error for FileOpenError {}

pub fn find<'a>(error: &'a (dyn std::error::Error + 'static)) -> Option<&'a FileOpenError> {
    let mut current = Some(error);
    while let Some(error) = current {
        if let Some(open) = error.downcast_ref::<FileOpenError>() {
            return Some(open);
        }
        if let Some(io) = error.downcast_ref::<std::io::Error>() {
            if let Some(inner) = io.get_ref() {
                current = Some(inner);
                continue;
            }
        }
        current = error.source();
    }
    None
}

/// Source POSIX names are opaque bytes. Windows opens stay wide and diagnostics
/// use the CRT/ACP compatibility domain, independently of Rust's UTF-8 strings.
pub fn path_diagnostic_bytes(path: &std::path::Path) -> super::IoResult<Vec<u8>> {
    super::source_os_bytes(path.as_os_str()).map(|bytes| bytes.into_owned())
}
