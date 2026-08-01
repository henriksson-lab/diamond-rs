use std::fs::{File as StdFile, OpenOptions};
use std::io::{Read, Seek, SeekFrom, Write};
#[cfg(unix)]
use std::os::fd::{FromRawFd, IntoRawFd};
use std::sync::Mutex as StdMutex;

#[cfg(unix)]
unsafe extern "C" {
    fn dup(oldfd: i32) -> i32;
}

const DEFAULT_FILE_BUFFER_SIZE: usize = 1 << 20;
pub const MEGABYTES: usize = 1 << 20;
pub const GIGABYTES: usize = 1 << 30;
pub const KILOBYTES: usize = 1 << 10;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compressor {
    None,
    Zlib,
    Zstd,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TempFileData {
    pub name: String,
    pub fd: i32,
    pub unlinked: bool,
}

impl TempFileData {
    pub fn init(unlink: bool) -> IoResult<Self> {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "diamond-tmp-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
            .map_err(|_| IoError::Other(format!("Error opening temporary file {name}")))?;
        #[cfg(unix)]
        let fd = file.into_raw_fd();
        #[cfg(not(unix))]
        let fd = {
            let _ = &file;
            -1
        };
        let unlinked = if unlink {
            std::fs::remove_file(&path).is_ok()
        } else {
            false
        };
        Ok(Self { name, fd, unlinked })
    }
}

pub fn detect_compressor(b: &[u8]) -> Compressor {
    if b.len() >= 2
        && ((b[0] == 0x1f && b[1] == 0x8b)
            || (b[0] == 0x78 && (b[1] == 0x01 || b[1] == 0x9c || b[1] == 0xda)))
    {
        return Compressor::Zlib;
    }
    if b.len() >= 4 && b[0] == 0x28 && b[1] == 0xb5 && b[2] == 0x2f && b[3] == 0xfd {
        return Compressor::Zstd;
    }
    Compressor::None
}

pub fn decompress(src: &[u8], dst: &mut [u8], compressor: Compressor) -> IoResult<usize> {
    match compressor {
        Compressor::Zlib => zlib_decompress(src, dst),
        Compressor::Zstd => zstd_decompress(src, dst),
        Compressor::None => Err(IoError::Other(
            "Invalid compressor in decompress".to_string(),
        )),
    }
}

pub fn zlib_decompress(src: &[u8], dst: &mut [u8]) -> IoResult<usize> {
    let mut out = Vec::new();
    if src.len() >= 2 && src[0] == 0x1f && src[1] == 0x8b {
        let mut decoder = flate2::read::MultiGzDecoder::new(src);
        decoder
            .read_to_end(&mut out)
            .map_err(|e| IoError::Other(format!("Error during zlib decompression: {e}")))?;
    } else {
        let mut remaining = src;
        while !remaining.is_empty() {
            let mut decoder = flate2::read::ZlibDecoder::new(remaining);
            decoder
                .read_to_end(&mut out)
                .map_err(|e| IoError::Other(format!("Error during zlib decompression: {e}")))?;
            let consumed = decoder.total_in() as usize;
            if consumed == 0 {
                return Err(IoError::Other("Unexpected end of zlib stream".to_string()));
            }
            remaining = &remaining[consumed..];
        }
    }
    if out.len() > dst.len() {
        return Err(IoError::Other(
            "zlib_decompress: output buffer too small".to_string(),
        ));
    }
    dst[..out.len()].copy_from_slice(&out);
    Ok(out.len())
}

pub fn zstd_decompress(src: &[u8], dst: &mut [u8]) -> IoResult<usize> {
    let out = zstd::stream::decode_all(src)
        .map_err(|e| IoError::Other(format!("Failed decompressing zstd stream: {e}")))?;
    if out.len() > dst.len() {
        return Err(IoError::Other(
            "Failed decompressing zstd stream: output buffer too small".to_string(),
        ));
    }
    dst[..out.len()].copy_from_slice(&out);
    Ok(out.len())
}

#[cfg(unix)]
pub fn posix_flags(mode: &str) -> IoResult<i32> {
    const O_WRONLY: i32 = 1;
    const O_RDWR: i32 = 2;
    const O_CREAT: i32 = 64;
    const O_TRUNC: i32 = 512;
    match mode {
        "wb" => Ok(O_WRONLY | O_CREAT | O_TRUNC),
        "r+b" => Ok(O_RDWR),
        "w+b" => Ok(O_RDWR | O_CREAT | O_TRUNC),
        _ => Err(IoError::Other("Invalid fopen mode.".to_string())),
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IoError {
    UnsupportedOperation,
    FileOpen(String),
    FileRead(String),
    FileWrite(String),
    EndOfStream,
    StreamRead { line_count: usize, msg: String },
    Other(String),
}

impl std::fmt::Display for IoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IoError::UnsupportedOperation => f.write_str("Unsupported I/O operation."),
            IoError::FileOpen(file_name) => write!(f, "Error opening file {file_name}"),
            IoError::FileRead(file_name) => write!(f, "Error reading file {file_name}"),
            IoError::FileWrite(file_name) => write!(f, "Error writing file {file_name}"),
            IoError::EndOfStream => f.write_str("Unexpected end of input."),
            IoError::StreamRead { line_count, msg } => {
                write!(f, "Error reading input stream at line {line_count}: {msg}")
            }
            IoError::Other(msg) => f.write_str(msg),
        }
    }
}

impl std::error::Error for IoError {}

pub type IoResult<T> = Result<T, IoError>;

pub trait Consumer {
    fn consume(&mut self, ptr: &[u8]) -> IoResult<()>;
    fn finalize(&mut self) -> IoResult<()> {
        Ok(())
    }
}

pub trait StreamEntity {
    fn rewind(&mut self) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn seek(&mut self, _p: i64, _origin: SeekFrom) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn tell(&mut self) -> IoResult<i64> {
        Err(IoError::UnsupportedOperation)
    }

    fn read(&mut self, _ptr: &mut [u8]) -> IoResult<usize> {
        Err(IoError::UnsupportedOperation)
    }

    fn fetch(&mut self) -> IoResult<bool> {
        Err(IoError::UnsupportedOperation)
    }

    fn close(&mut self) -> IoResult<()> {
        Ok(())
    }

    fn file_name(&self) -> &str {
        ""
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        Err(IoError::UnsupportedOperation)
    }

    fn data(&self) -> &[u8] {
        &[]
    }

    fn write(&mut self, _ptr: &[u8]) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn flush(&mut self) -> IoResult<()> {
        Err(IoError::UnsupportedOperation)
    }

    fn eof(&mut self) -> IoResult<bool> {
        Err(IoError::UnsupportedOperation)
    }

    fn file_size(&mut self) -> IoResult<i64> {
        Err(IoError::UnsupportedOperation)
    }

    fn seekable(&self) -> bool {
        false
    }
}

#[derive(Debug)]
pub struct File {
    file: Option<StdFile>,
    auto_delete: bool,
    unlinked: bool,
    file_name: String,
    scratch: Vec<u8>,
}

impl File {
    pub fn new(name: &str, mode: &str) -> IoResult<Self> {
        Self::open(name, mode)
    }

    pub fn new_temporary(_temporary: Temporary) -> IoResult<Self> {
        Self::temporary()
    }

    pub fn open(name: &str, mode: &str) -> IoResult<Self> {
        let mut options = OpenOptions::new();
        if mode.contains('r') {
            options.read(true);
        }
        if mode.contains('w') {
            options.write(true).create(true).truncate(true);
        }
        if mode.contains('a') {
            options.append(true).create(true);
        }
        if mode.contains('+') {
            options.read(true).write(true);
        }
        let file = options
            .open(name)
            .map_err(|e| IoError::Other(format!("Error opening file {name}. {e}")))?;
        Ok(Self {
            file: Some(file),
            auto_delete: false,
            unlinked: false,
            file_name: name.to_string(),
            scratch: Vec::new(),
        })
    }

    pub fn temporary() -> IoResult<Self> {
        #[cfg(unix)]
        {
            let data = TempFileData::init(true)?;
            return Ok(Self {
                file: Some(unsafe { StdFile::from_raw_fd(data.fd) }),
                auto_delete: true,
                unlinked: data.unlinked,
                file_name: data.name,
                scratch: Vec::new(),
            });
        }
        #[cfg(not(unix))]
        {
            let mut path = std::env::temp_dir();
            path.push(format!(
                "diamond-rs-{}-{}.tmp",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ));
            let name = path.to_string_lossy().into_owned();
            let file = OpenOptions::new()
                .read(true)
                .write(true)
                .create_new(true)
                .open(&path)
                .map_err(|e| IoError::Other(format!("Error opening temporary file {name}. {e}")))?;
            Ok(Self {
                file: Some(file),
                auto_delete: true,
                unlinked: false,
                file_name: name,
                scratch: Vec::new(),
            })
        }
    }

    pub fn close(&mut self) -> IoResult<()> {
        self.file.take();
        if self.auto_delete && !self.unlinked {
            let _ = std::fs::remove_file(&self.file_name);
        }
        Ok(())
    }

    pub fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::FileWrite(self.file_name.clone()))?
            .write_all(ptr)
            .map_err(|_| IoError::FileWrite(self.file_name.clone()))
    }

    pub fn write_value<T: FilePrimitive>(&mut self, x: T) -> IoResult<()> {
        self.write(&x.to_ne_bytes_vec())
    }

    pub fn read_value<T: FilePrimitive>(&mut self) -> IoResult<T> {
        let mut buf = vec![0u8; T::SIZE];
        self.read_exact(&mut buf)?;
        T::from_ne_bytes_slice(&buf)
    }

    pub fn seek(&mut self, p: i64, origin: SeekFrom) -> IoResult<()> {
        let target = match origin {
            SeekFrom::Start(_) => SeekFrom::Start(p as u64),
            SeekFrom::End(_) => SeekFrom::End(p),
            SeekFrom::Current(_) => SeekFrom::Current(p),
        };
        self.file
            .as_mut()
            .ok_or_else(|| IoError::Other("Error calling fseek.".to_string()))?
            .seek(target)
            .map(|_| ())
            .map_err(|_| IoError::Other("Error calling fseek.".to_string()))
    }

    pub fn tell(&mut self) -> IoResult<i64> {
        self.file
            .as_mut()
            .ok_or_else(|| {
                IoError::Other(format!(
                    "Error executing ftell on stream {}",
                    self.file_name
                ))
            })?
            .stream_position()
            .map(|p| p as i64)
            .map_err(|_| IoError::Other("Error calling ftell.".to_string()))
    }

    pub fn size(&mut self) -> IoResult<i64> {
        let pos = self.tell()?;
        self.seek(0, SeekFrom::End(0))?;
        let size = self.tell()?;
        self.seek(pos, SeekFrom::Start(0))?;
        Ok(size)
    }

    pub fn read_exact(&mut self, ptr: &mut [u8]) -> IoResult<()> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::FileRead(self.file_name.clone()))?
            .read_exact(ptr)
            .map_err(|_| IoError::FileRead(self.file_name.clone()))
    }

    pub fn read_max(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::FileRead(self.file_name.clone()))?
            .read(ptr)
            .map_err(|_| IoError::FileRead(self.file_name.clone()))
    }

    pub fn read(&mut self, n: usize) -> IoResult<&[u8]> {
        self.scratch.resize(n, 0);
        let mut buf = std::mem::take(&mut self.scratch);
        self.read_exact(&mut buf)?;
        self.scratch = buf;
        Ok(&self.scratch)
    }

    pub fn eof(&mut self) -> IoResult<bool> {
        let file = self
            .file
            .as_mut()
            .ok_or_else(|| IoError::FileRead(self.file_name.clone()))?;
        let pos = file
            .stream_position()
            .map_err(|_| IoError::Other("Error calling ftell.".to_string()))?;
        let len = file
            .metadata()
            .map_err(|_| IoError::FileRead(self.file_name.clone()))?
            .len();
        Ok(pos >= len)
    }

    pub fn file_name(&self) -> &str {
        &self.file_name
    }

    pub fn file(&mut self) -> IoResult<&mut StdFile> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::Other(format!("Closed file {}", self.file_name)))
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Temporary;

impl Drop for File {
    fn drop(&mut self) {
        let _ = self.close();
    }
}

pub trait FilePrimitive: Copy {
    const SIZE: usize;
    fn to_ne_bytes_vec(self) -> Vec<u8>;
    fn from_ne_bytes_slice(bytes: &[u8]) -> IoResult<Self>;
}

macro_rules! impl_file_primitive {
    ($($t:ty),* $(,)?) => {
        $(
            impl FilePrimitive for $t {
                const SIZE: usize = std::mem::size_of::<$t>();

                fn to_ne_bytes_vec(self) -> Vec<u8> {
                    self.to_ne_bytes().to_vec()
                }

                fn from_ne_bytes_slice(bytes: &[u8]) -> IoResult<Self> {
                    Ok(<$t>::from_ne_bytes(bytes.try_into().unwrap()))
                }
            }
        )*
    };
}

impl_file_primitive!(u8, i8, u16, i16, u32, i32, u64, i64, f32, f64);

pub struct CompressedBuffer {
    encoder: Option<flate2::write::ZlibEncoder<Vec<u8>>>,
    buf: Vec<u8>,
}

impl Default for CompressedBuffer {
    fn default() -> Self {
        Self::new()
    }
}

impl CompressedBuffer {
    pub const BUF_SIZE: usize = 32768;

    pub fn new() -> Self {
        let mut out = Self {
            encoder: None,
            buf: Vec::with_capacity(Self::BUF_SIZE),
        };
        out.clear();
        out
    }

    pub fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.encoder
            .as_mut()
            .ok_or_else(|| IoError::Other("CompressedBuffer stream is closed.".to_string()))?
            .write_all(ptr)
            .map_err(|e| IoError::Other(e.to_string()))
    }

    pub fn write_value<T: FilePrimitive>(&mut self, x: T) -> IoResult<()> {
        self.write(&x.to_ne_bytes_vec())
    }

    pub fn finish(&mut self) -> IoResult<()> {
        if let Some(encoder) = self.encoder.take() {
            self.buf = encoder
                .finish()
                .map_err(|e| IoError::Other(e.to_string()))?;
        }
        Ok(())
    }

    pub fn clear(&mut self) {
        self.buf.clear();
        self.encoder = Some(flate2::write::ZlibEncoder::new(
            Vec::with_capacity(Self::BUF_SIZE),
            flate2::Compression::default(),
        ));
    }

    pub fn data(&self) -> &[u8] {
        if let Some(encoder) = &self.encoder {
            encoder.get_ref()
        } else {
            &self.buf
        }
    }

    pub fn size(&self) -> usize {
        self.data().len()
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct VecStream {
    data: Vec<u8>,
    pos: usize,
    name: String,
}

impl VecStream {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_vec(data: Vec<u8>) -> Self {
        Self {
            data,
            pos: 0,
            name: String::new(),
        }
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }
}

impl StreamEntity for VecStream {
    fn rewind(&mut self) -> IoResult<()> {
        self.pos = 0;
        Ok(())
    }

    fn seek(&mut self, p: i64, origin: SeekFrom) -> IoResult<()> {
        let pos = match origin {
            SeekFrom::Start(_) => p,
            SeekFrom::End(_) => self.data.len() as i64 + p,
            SeekFrom::Current(_) => self.pos as i64 + p,
        };
        if pos < 0 {
            return Err(IoError::UnsupportedOperation);
        }
        self.pos = pos as usize;
        if self.pos > self.data.len() {
            self.data.resize(self.pos, 0);
        }
        Ok(())
    }

    fn tell(&mut self) -> IoResult<i64> {
        Ok(self.pos as i64)
    }

    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        let n = ptr.len().min(self.data.len().saturating_sub(self.pos));
        ptr[..n].copy_from_slice(&self.data[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }

    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        if self.pos + ptr.len() > self.data.len() {
            self.data.resize(self.pos + ptr.len(), 0);
        }
        self.data[self.pos..self.pos + ptr.len()].copy_from_slice(ptr);
        self.pos += ptr.len();
        Ok(())
    }

    fn flush(&mut self) -> IoResult<()> {
        Ok(())
    }

    fn eof(&mut self) -> IoResult<bool> {
        Ok(self.pos >= self.data.len())
    }

    fn file_size(&mut self) -> IoResult<i64> {
        Ok(self.data.len() as i64)
    }

    fn seekable(&self) -> bool {
        true
    }

    fn file_name(&self) -> &str {
        &self.name
    }
}

#[derive(Debug)]
pub struct FileSource {
    file: StdFile,
    file_name: String,
    seekable: bool,
    eof: bool,
}

impl FileSource {
    pub fn new(file_name: &str) -> IoResult<Self> {
        let is_stdin = file_name.is_empty() || file_name == "-";
        let open_name = if is_stdin { "/dev/stdin" } else { file_name };
        let file =
            StdFile::open(open_name).map_err(|_| IoError::FileOpen(file_name.to_string()))?;
        let seekable = !is_stdin && file.metadata().map(|m| m.is_file()).unwrap_or(false);
        Ok(Self {
            file,
            file_name: file_name.to_string(),
            seekable,
            eof: false,
        })
    }

    pub fn from_file(file_name: &str, file: StdFile) -> Self {
        Self {
            file,
            file_name: file_name.to_string(),
            seekable: false,
            eof: false,
        }
    }
}

impl StreamEntity for FileSource {
    fn rewind(&mut self) -> IoResult<()> {
        self.file
            .seek(SeekFrom::Start(0))
            .map(|_| {
                self.eof = false;
            })
            .map_err(|_| IoError::Other(format!("Error executing seek on file {}", self.file_name)))
    }

    fn seek(&mut self, p: i64, origin: SeekFrom) -> IoResult<()> {
        let target = match origin {
            SeekFrom::Start(_) => SeekFrom::Start(p as u64),
            SeekFrom::End(_) => SeekFrom::End(p),
            SeekFrom::Current(_) => SeekFrom::Current(p),
        };
        self.file
            .seek(target)
            .map(|_| {
                self.eof = false;
            })
            .map_err(|_| IoError::Other(format!("Error executing seek on file {}", self.file_name)))
    }

    fn tell(&mut self) -> IoResult<i64> {
        self.file.stream_position().map(|p| p as i64).map_err(|_| {
            IoError::Other(format!(
                "Error executing ftell on stream {}",
                self.file_name
            ))
        })
    }

    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        match self.file.read(ptr) {
            Ok(n) => {
                self.eof = n == 0 || n < ptr.len();
                Ok(n)
            }
            Err(_) => Err(IoError::FileRead(self.file_name.clone())),
        }
    }

    fn close(&mut self) -> IoResult<()> {
        Ok(())
    }

    fn file_size(&mut self) -> IoResult<i64> {
        let pos = self.tell()?;
        self.seek(0, SeekFrom::End(0))?;
        let size = self.tell()?;
        self.seek(pos, SeekFrom::Start(0))?;
        Ok(size)
    }

    fn eof(&mut self) -> IoResult<bool> {
        Ok(self.eof)
    }

    fn file_name(&self) -> &str {
        &self.file_name
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        Ok(&mut self.file)
    }

    fn seekable(&self) -> bool {
        self.seekable
    }
}

#[derive(Debug)]
pub struct FileSink {
    file: Option<StdFile>,
    file_name: String,
    async_: bool,
    mtx: StdMutex<()>,
}

impl FileSink {
    pub fn new(file_name: &str, mode: &str, async_: bool, _buffer_size: usize) -> IoResult<Self> {
        #[cfg(unix)]
        let _ = posix_flags(mode)?;
        let mut options = OpenOptions::new();
        let is_stdout = file_name.is_empty();
        if is_stdout {
            options.write(true);
            if mode.contains('+') {
                options.read(true);
            }
        } else {
            match mode {
                "wb" => {
                    options.write(true).create(true).truncate(true);
                }
                "r+b" => {
                    options.read(true).write(true);
                }
                "w+b" => {
                    options.read(true).write(true).create(true).truncate(true);
                }
                _ => return Err(IoError::Other("Invalid fopen mode.".to_string())),
            }
        }
        let open_name = if is_stdout { "/dev/stdout" } else { file_name };
        let file = options
            .open(open_name)
            .map_err(|_| IoError::FileOpen(file_name.to_string()))?;
        Ok(Self {
            file: Some(file),
            file_name: file_name.to_string(),
            async_,
            mtx: StdMutex::new(()),
        })
    }

    #[cfg(unix)]
    pub unsafe fn from_fd(
        file_name: &str,
        fd: i32,
        mode: &str,
        async_: bool,
        _buffer_size: usize,
    ) -> IoResult<Self> {
        let _ = posix_flags(mode)?;
        Ok(Self {
            file: Some(unsafe { StdFile::from_raw_fd(fd) }),
            file_name: file_name.to_string(),
            async_,
            mtx: StdMutex::new(()),
        })
    }
}

impl StreamEntity for FileSink {
    fn close(&mut self) -> IoResult<()> {
        self.file.take();
        Ok(())
    }

    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        if self.async_ {
            let _guard = self.mtx.lock().unwrap();
            self.file
                .as_mut()
                .ok_or_else(|| IoError::FileWrite(self.file_name.clone()))?
                .write_all(ptr)
                .map_err(|_| IoError::FileWrite(self.file_name.clone()))
        } else {
            self.file
                .as_mut()
                .ok_or_else(|| IoError::FileWrite(self.file_name.clone()))?
                .write_all(ptr)
                .map_err(|_| IoError::FileWrite(self.file_name.clone()))
        }
    }

    fn seek(&mut self, p: i64, origin: SeekFrom) -> IoResult<()> {
        let target = match origin {
            SeekFrom::Start(_) => SeekFrom::Start(p as u64),
            SeekFrom::End(_) => SeekFrom::End(p),
            SeekFrom::Current(_) => SeekFrom::Current(p),
        };
        self.file
            .as_mut()
            .ok_or_else(|| IoError::Other("Error calling fseek.".to_string()))?
            .seek(target)
            .map(|_| ())
            .map_err(|_| IoError::Other("Error calling fseek.".to_string()))
    }

    fn rewind(&mut self) -> IoResult<()> {
        self.seek(0, SeekFrom::Start(0))
    }

    fn tell(&mut self) -> IoResult<i64> {
        self.file
            .as_mut()
            .ok_or_else(|| {
                IoError::Other(format!(
                    "Error executing ftell on stream {}",
                    self.file_name
                ))
            })?
            .stream_position()
            .map(|p| p as i64)
            .map_err(|_| IoError::Other("Error calling ftell.".to_string()))
    }

    fn file_name(&self) -> &str {
        &self.file_name
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        self.file
            .as_mut()
            .ok_or_else(|| IoError::FileWrite(self.file_name.clone()))
    }

    fn flush(&mut self) -> IoResult<()> {
        self.file
            .as_mut()
            .map(|f| {
                f.flush()
                    .map_err(|_| IoError::FileWrite(self.file_name.clone()))
            })
            .unwrap_or(Ok(()))
    }

    fn file_size(&mut self) -> IoResult<i64> {
        let pos = self.tell()?;
        self.seek(0, SeekFrom::End(0))?;
        let size = self.tell()?;
        self.seek(pos, SeekFrom::Start(0))?;
        Ok(size)
    }

    fn seekable(&self) -> bool {
        true
    }
}

#[derive(Debug, Clone)]
pub struct InputStreamBuffer<S: StreamEntity> {
    prev: S,
    buf_size: usize,
    buf: Vec<u8>,
    begin: usize,
    end: usize,
    file_offset: usize,
    async_: bool,
}

impl<S: StreamEntity> InputStreamBuffer<S> {
    pub const ASYNC: i32 = 4;

    pub fn new(prev: S, flags: i32) -> Self {
        Self {
            prev,
            buf_size: DEFAULT_FILE_BUFFER_SIZE,
            buf: vec![0; DEFAULT_FILE_BUFFER_SIZE],
            begin: 0,
            end: 0,
            file_offset: 0,
            async_: (flags & Self::ASYNC) != 0,
        }
    }

    pub fn data(&self) -> &[u8] {
        &self.buf[self.begin..self.end]
    }

    pub fn pop(&mut self, dst: &mut [u8]) -> usize {
        let n = dst.len().min(self.end.saturating_sub(self.begin));
        dst[..n].copy_from_slice(&self.buf[self.begin..self.begin + n]);
        self.begin += n;
        n
    }
}

impl<S: StreamEntity> StreamEntity for InputStreamBuffer<S> {
    fn rewind(&mut self) -> IoResult<()> {
        self.prev.rewind()?;
        self.file_offset = 0;
        self.begin = 0;
        self.end = 0;
        Ok(())
    }

    fn seek(&mut self, p: i64, origin: SeekFrom) -> IoResult<()> {
        self.prev.seek(p, origin)?;
        self.file_offset = 0;
        self.begin = 0;
        self.end = 0;
        Ok(())
    }

    fn fetch(&mut self) -> IoResult<bool> {
        let n = self.prev.read(&mut self.buf[..self.buf_size])?;
        if self.prev.seekable() {
            self.file_offset = self.prev.tell()? as usize;
        }
        self.begin = 0;
        self.end = n;
        let _ = self.async_;
        Ok(self.end > self.begin)
    }

    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        let mut total = 0usize;
        while total < ptr.len() {
            if self.begin == self.end && !self.fetch()? {
                break;
            }
            total += self.pop(&mut ptr[total..]);
        }
        Ok(total)
    }

    fn close(&mut self) -> IoResult<()> {
        self.prev.close()
    }

    fn tell(&mut self) -> IoResult<i64> {
        if !self.seekable() {
            return Err(IoError::Other(
                "Calling tell on non seekable stream.".to_string(),
            ));
        }
        Ok(self.file_offset as i64 - (self.end - self.begin) as i64)
    }

    fn eof(&mut self) -> IoResult<bool> {
        self.prev.eof()
    }

    fn file_size(&mut self) -> IoResult<i64> {
        self.prev.file_size()
    }

    fn seekable(&self) -> bool {
        self.prev.seekable()
    }

    fn file_name(&self) -> &str {
        self.prev.file_name()
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        self.prev.file()
    }

    fn data(&self) -> &[u8] {
        &self.buf[self.begin..self.end]
    }
}

struct PrevReader<S: StreamEntity>(InputStreamBuffer<S>);

impl<S: StreamEntity> std::io::Read for PrevReader<S> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.0
            .read(buf)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, format!("{e}")))
    }
}

type PeekChain<S> = std::io::Chain<std::io::Cursor<Vec<u8>>, PrevReader<S>>;

enum ZlibInner<S: StreamEntity> {
    Gz(flate2::read::MultiGzDecoder<PeekChain<S>>),
    Zl(flate2::read::ZlibDecoder<PeekChain<S>>),
}

impl<S: StreamEntity> ZlibInner<S> {
    fn into_prev(self) -> InputStreamBuffer<S> {
        let chain = match self {
            ZlibInner::Gz(d) => d.into_inner(),
            ZlibInner::Zl(d) => d.into_inner(),
        };
        let (_cursor, prev_reader) = chain.into_inner();
        prev_reader.0
    }

    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        use std::io::Read as _;
        match self {
            ZlibInner::Gz(d) => d.read(buf),
            ZlibInner::Zl(d) => d.read(buf),
        }
    }
}

fn zlib_wrap<S: StreamEntity>(mut prev: InputStreamBuffer<S>) -> IoResult<ZlibInner<S>> {
    use std::io::Read as _;
    let mut peek = [0u8; 2];
    let n = prev.read(&mut peek)?;
    let peeked = peek[..n].to_vec();
    let is_gz = n >= 2 && peeked[0] == 0x1f && peeked[1] == 0x8b;
    let chain = std::io::Cursor::new(peeked).chain(PrevReader(prev));
    if is_gz {
        Ok(ZlibInner::Gz(flate2::read::MultiGzDecoder::new(chain)))
    } else {
        Ok(ZlibInner::Zl(flate2::read::ZlibDecoder::new(chain)))
    }
}

pub struct ZlibSource<S: StreamEntity> {
    inner: Option<ZlibInner<S>>,
    eos: bool,
    file_name_cache: String,
}

impl<S: StreamEntity> ZlibSource<S> {
    pub fn new(prev: InputStreamBuffer<S>) -> IoResult<Self> {
        let file_name_cache = prev.file_name().to_string();
        let inner = zlib_wrap(prev)?;
        Ok(Self {
            inner: Some(inner),
            eos: false,
            file_name_cache,
        })
    }
}

impl<S: StreamEntity> StreamEntity for ZlibSource<S> {
    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        if self.eos {
            return Ok(0);
        }
        let inner = self
            .inner
            .as_mut()
            .ok_or_else(|| IoError::Other("ZlibSource closed".to_string()))?;
        let mut total = 0;
        while total < ptr.len() {
            let n = inner.read(&mut ptr[total..]).map_err(|_| {
                IoError::Other(format!(
                    "Error reading gzip-compressed input file. The file may be corrupted: {}",
                    self.file_name_cache
                ))
            })?;
            if n == 0 {
                self.eos = true;
                break;
            }
            total += n;
        }
        Ok(total)
    }

    fn close(&mut self) -> IoResult<()> {
        if let Some(inner) = self.inner.take() {
            let mut prev = inner.into_prev();
            prev.close()?;
        }
        Ok(())
    }

    fn rewind(&mut self) -> IoResult<()> {
        let inner = self
            .inner
            .take()
            .ok_or_else(|| IoError::Other("ZlibSource closed".to_string()))?;
        let mut prev = inner.into_prev();
        prev.rewind()?;
        self.inner = Some(zlib_wrap(prev)?);
        self.eos = false;
        Ok(())
    }

    fn eof(&mut self) -> IoResult<bool> {
        Ok(self.eos)
    }

    fn file_name(&self) -> &str {
        &self.file_name_cache
    }
}

pub struct ZstdSource<S: StreamEntity> {
    inner: Option<zstd::stream::Decoder<'static, std::io::BufReader<PrevReader<S>>>>,
    eos: bool,
    file_name_cache: String,
}

impl<S: StreamEntity> ZstdSource<S> {
    pub fn new(prev: InputStreamBuffer<S>) -> IoResult<Self> {
        let file_name_cache = prev.file_name().to_string();
        let dec = zstd::stream::Decoder::new(PrevReader(prev))
            .map_err(|e| IoError::Other(format!("ZSTD_decompressStream: {e}")))?;
        Ok(Self {
            inner: Some(dec),
            eos: false,
            file_name_cache,
        })
    }
}

impl<S: StreamEntity> StreamEntity for ZstdSource<S> {
    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        if self.eos {
            return Ok(0);
        }
        use std::io::Read as _;
        let dec = self
            .inner
            .as_mut()
            .ok_or_else(|| IoError::Other("ZstdSource closed".to_string()))?;
        let mut total = 0;
        while total < ptr.len() {
            let n = dec
                .read(&mut ptr[total..])
                .map_err(|e| IoError::Other(format!("ZSTD_decompressStream: {e}")))?;
            if n == 0 {
                self.eos = true;
                break;
            }
            total += n;
        }
        Ok(total)
    }

    fn close(&mut self) -> IoResult<()> {
        if let Some(dec) = self.inner.take() {
            let mut prev = dec.finish().into_inner().0;
            prev.close()?;
        }
        Ok(())
    }

    fn rewind(&mut self) -> IoResult<()> {
        let dec = self
            .inner
            .take()
            .ok_or_else(|| IoError::Other("ZstdSource closed".to_string()))?;
        let mut prev = dec.finish().into_inner().0;
        prev.rewind()?;
        let new_dec = zstd::stream::Decoder::new(PrevReader(prev))
            .map_err(|e| IoError::Other(format!("ZSTD_decompressStream: {e}")))?;
        self.inner = Some(new_dec);
        self.eos = false;
        Ok(())
    }

    fn eof(&mut self) -> IoResult<bool> {
        Ok(self.eos)
    }

    fn file_name(&self) -> &str {
        &self.file_name_cache
    }
}

#[derive(Debug, Clone)]
pub struct OutputStreamBuffer<S: StreamEntity> {
    prev: S,
    buf_size: usize,
    buf: Vec<u8>,
}

impl<S: StreamEntity> OutputStreamBuffer<S> {
    pub const STDOUT_BUF_SIZE: usize = 4096;

    pub fn new(prev: S) -> Self {
        let buf_size = if prev.file_name().is_empty() {
            Self::STDOUT_BUF_SIZE
        } else {
            DEFAULT_FILE_BUFFER_SIZE
        };
        Self {
            prev,
            buf_size,
            buf: Vec::with_capacity(buf_size),
        }
    }

    pub fn write_buffer(&mut self) -> &mut Vec<u8> {
        &mut self.buf
    }

    pub fn write_buffer_range(&mut self) -> (&mut [u8], usize) {
        if self.buf.len() != self.buf_size {
            self.buf.resize(self.buf_size, 0);
        }
        (self.buf.as_mut_slice(), self.buf_size)
    }

    pub fn flush_count(&mut self, count: usize) -> IoResult<()> {
        if self.buf.len() < count {
            return Err(IoError::Other(
                "OutputStreamBuffer::flush count exceeds buffer size.".to_string(),
            ));
        }
        self.prev.write(&self.buf[..count])?;
        self.buf.clear();
        Ok(())
    }
}

impl<S: StreamEntity> StreamEntity for OutputStreamBuffer<S> {
    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.buf.extend_from_slice(ptr);
        if self.buf.len() >= self.buf_size {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> IoResult<()> {
        if !self.buf.is_empty() {
            self.prev.write(&self.buf)?;
            self.buf.clear();
        }
        self.prev.flush()
    }

    fn seek(&mut self, p: i64, origin: SeekFrom) -> IoResult<()> {
        self.flush()?;
        self.prev.seek(p, origin)
    }

    fn rewind(&mut self) -> IoResult<()> {
        self.flush()?;
        self.prev.rewind()
    }

    fn tell(&mut self) -> IoResult<i64> {
        self.flush()?;
        self.prev.tell()
    }

    fn close(&mut self) -> IoResult<()> {
        self.flush()?;
        self.prev.close()
    }

    fn file_name(&self) -> &str {
        self.prev.file_name()
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        self.flush()?;
        self.prev.file()
    }

    fn file_size(&mut self) -> IoResult<i64> {
        self.flush()?;
        self.prev.file_size()
    }

    fn seekable(&self) -> bool {
        self.prev.seekable()
    }
}

#[derive(Debug, Clone)]
pub struct ZlibSink<S: StreamEntity> {
    prev: OutputStreamBuffer<S>,
    pending: Vec<u8>,
    closed: bool,
}

impl<S: StreamEntity> ZlibSink<S> {
    pub const CHUNK_SIZE: usize = 1 << 20;

    pub fn new(prev: OutputStreamBuffer<S>) -> Self {
        Self {
            prev,
            pending: Vec::new(),
            closed: false,
        }
    }

    pub fn deflate_loop(&mut self, ptr: &[u8], finish: bool) -> IoResult<()> {
        if self.closed {
            return Ok(());
        }
        self.pending.extend_from_slice(ptr);
        if !finish {
            return Ok(());
        }
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder
            .write_all(&self.pending)
            .map_err(|e| IoError::Other(format!("deflate error: {e}")))?;
        let compressed = encoder
            .finish()
            .map_err(|e| IoError::Other(format!("deflate error: {e}")))?;
        self.prev.write(&compressed)?;
        self.pending.clear();
        self.closed = true;
        Ok(())
    }
}

impl<S: StreamEntity> StreamEntity for ZlibSink<S> {
    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.deflate_loop(ptr, false)
    }

    fn close(&mut self) -> IoResult<()> {
        self.deflate_loop(&[], true)?;
        self.prev.close()
    }

    fn flush(&mut self) -> IoResult<()> {
        self.prev.flush()
    }

    fn file_name(&self) -> &str {
        self.prev.file_name()
    }
}

#[derive(Debug, Clone)]
pub struct ZstdSink<S: StreamEntity> {
    prev: OutputStreamBuffer<S>,
    pending: Vec<u8>,
    closed: bool,
}

impl<S: StreamEntity> ZstdSink<S> {
    pub fn new(prev: OutputStreamBuffer<S>) -> Self {
        Self {
            prev,
            pending: Vec::new(),
            closed: false,
        }
    }
}

impl<S: StreamEntity> StreamEntity for ZstdSink<S> {
    fn write(&mut self, ptr: &[u8]) -> IoResult<()> {
        if !self.closed {
            self.pending.extend_from_slice(ptr);
        }
        Ok(())
    }

    fn close(&mut self) -> IoResult<()> {
        if !self.closed {
            let compressed = zstd::stream::encode_all(&self.pending[..], 0)
                .map_err(|e| IoError::Other(format!("ZSTD_endStream: {e}")))?;
            self.prev.write(&compressed)?;
            self.pending.clear();
            self.closed = true;
        }
        self.prev.close()
    }

    fn flush(&mut self) -> IoResult<()> {
        self.prev.flush()
    }

    fn file_name(&self) -> &str {
        self.prev.file_name()
    }
}

#[derive(Debug, Clone)]
pub struct Serializer<S: StreamEntity> {
    buffer: S,
    pending: Vec<u8>,
    pending_pos: usize,
}

impl<S: StreamEntity> Serializer<S> {
    pub fn new(buffer: S) -> Self {
        Self {
            buffer,
            pending: Vec::new(),
            pending_pos: 0,
        }
    }

    // The C++ `Serializer::operator<<` writes scalars in NATIVE endianness
    // (`util/io/serializer.h:89-102`). The `big_endian_byteswap` helper there
    // is misnamed — it actually maps to `psnip_endian_le16/32/64` which is a
    // no-op on x86_64 LE. We must emit native bytes too, or any `.dmnd`/.daa
    // written via `Serializer` (e.g. `taxonomy.save`) is non-interoperable
    // with C++.
    pub fn write_i32(&mut self, x: i32) -> IoResult<&mut Self> {
        self.write_raw(&x.to_ne_bytes())?;
        Ok(self)
    }

    pub fn write_i16(&mut self, x: i16) -> IoResult<&mut Self> {
        self.write_raw(&x.to_ne_bytes())?;
        Ok(self)
    }

    pub fn write_u16(&mut self, x: u16) -> IoResult<&mut Self> {
        self.write_raw(&x.to_ne_bytes())?;
        Ok(self)
    }

    pub fn write_i64(&mut self, x: i64) -> IoResult<&mut Self> {
        self.write_raw(&x.to_ne_bytes())?;
        Ok(self)
    }

    pub fn write_u32(&mut self, x: u32) -> IoResult<&mut Self> {
        self.write_raw(&x.to_ne_bytes())?;
        Ok(self)
    }

    pub fn write_u64(&mut self, x: u64) -> IoResult<&mut Self> {
        self.write_raw(&x.to_ne_bytes())?;
        Ok(self)
    }

    pub fn write_f64(&mut self, x: f64) -> IoResult<&mut Self> {
        self.write_raw(&x.to_ne_bytes())?;
        Ok(self)
    }

    pub fn write_value<T: FilePrimitive>(&mut self, x: T) -> IoResult<&mut Self> {
        self.write_raw(&x.to_ne_bytes_vec())?;
        Ok(self)
    }

    pub fn write_values<T: FilePrimitive>(&mut self, ptr: &[T]) -> IoResult<&mut Self> {
        for &x in ptr {
            self.write_raw(&x.to_ne_bytes_vec())?;
        }
        Ok(self)
    }

    pub fn write_string(&mut self, s: &str) -> IoResult<&mut Self> {
        self.write_raw(s.as_bytes())?;
        self.write_raw(&[0])?;
        Ok(self)
    }

    pub fn write_raw(&mut self, ptr: &[u8]) -> IoResult<()> {
        let mut offset = 0usize;
        while offset < ptr.len() {
            if self.pending.is_empty() || self.pending_pos == self.pending.len() {
                self.reset_buffer();
            }
            let n = (self.pending.len() - self.pending_pos).min(ptr.len() - offset);
            self.pending[self.pending_pos..self.pending_pos + n]
                .copy_from_slice(&ptr[offset..offset + n]);
            self.pending_pos += n;
            offset += n;
            if self.pending_pos == self.pending.len() {
                self.flush()?;
                self.reset_buffer();
            }
        }
        Ok(())
    }

    pub fn file_size(&mut self) -> IoResult<i64> {
        self.flush()?;
        self.buffer.file_size()
    }

    pub fn seek(&mut self, p: i64, origin: SeekFrom) -> IoResult<()> {
        self.flush()?;
        self.buffer.seek(p, origin)?;
        self.reset_buffer();
        Ok(())
    }

    pub fn rewind(&mut self) -> IoResult<()> {
        self.flush()?;
        self.buffer.rewind()?;
        self.reset_buffer();
        Ok(())
    }

    pub fn tell(&mut self) -> IoResult<usize> {
        self.flush()?;
        self.reset_buffer();
        self.buffer.tell().map(|p| p as usize)
    }

    pub fn close(&mut self) -> IoResult<()> {
        self.flush()?;
        self.buffer.close()
    }

    pub fn file_name(&self) -> String {
        self.buffer.file_name().to_string()
    }

    pub fn file(&mut self) -> IoResult<&mut StdFile> {
        self.flush()?;
        self.buffer.file()
    }

    pub fn consume(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.write_raw(ptr)
    }

    pub fn finalize(&mut self) -> IoResult<()> {
        self.close()
    }

    pub fn flush(&mut self) -> IoResult<()> {
        if self.pending_pos != 0 {
            self.buffer.write(&self.pending[..self.pending_pos])?;
            self.pending_pos = 0;
        }
        self.buffer.flush()
    }

    pub fn reset_buffer(&mut self) {
        if self.pending.is_empty() {
            self.pending = vec![0; DEFAULT_FILE_BUFFER_SIZE];
        }
        self.pending_pos = 0;
    }

    pub fn into_inner(mut self) -> IoResult<S> {
        self.flush()?;
        Ok(self.buffer)
    }
}

impl<S: StreamEntity> Consumer for Serializer<S> {
    fn consume(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.write_raw(ptr)
    }

    fn finalize(&mut self) -> IoResult<()> {
        self.close()
    }
}

#[derive(Debug, Clone)]
pub struct Deserializer<S: StreamEntity> {
    buffer: S,
}

impl<S: StreamEntity> Deserializer<S> {
    pub fn new(buffer: S) -> Self {
        Self { buffer }
    }

    pub fn rewind(&mut self) -> IoResult<()> {
        self.buffer.rewind()
    }

    pub fn seek(&mut self, pos: i64) -> IoResult<&mut Self> {
        self.buffer.seek(pos, SeekFrom::Start(0))?;
        Ok(self)
    }

    pub fn seek_forward(&mut self, n: usize) -> IoResult<()> {
        self.buffer.seek(n as i64, SeekFrom::Current(0))
    }

    pub fn seek_forward_delim(&mut self, delimiter: u8) -> IoResult<bool> {
        let mut b = [0u8; 1];
        loop {
            match self.buffer.read(&mut b)? {
                0 => return Ok(false),
                1 if b[0] == delimiter => return Ok(true),
                _ => {}
            }
        }
    }

    pub fn close(&mut self) -> IoResult<()> {
        self.buffer.close()
    }

    // C++ `Deserializer::operator>>` reads NATIVE endianness (matches the
    // native-write side in `Serializer`). The previous big-endian read here
    // would silently flip bytes on every scalar, producing garbage when
    // reading any non-Rust-written file. See `write_u32` above.
    pub fn read_u32(&mut self) -> IoResult<u32> {
        let mut b = [0u8; 4];
        self.read_exact(&mut b)?;
        Ok(u32::from_ne_bytes(b))
    }

    pub fn read_i32(&mut self) -> IoResult<i32> {
        let mut b = [0u8; 4];
        self.read_exact(&mut b)?;
        Ok(i32::from_ne_bytes(b))
    }

    pub fn read_i16(&mut self) -> IoResult<i16> {
        let mut b = [0u8; 2];
        self.read_exact(&mut b)?;
        Ok(i16::from_ne_bytes(b))
    }

    pub fn read_u16(&mut self) -> IoResult<u16> {
        let mut b = [0u8; 2];
        self.read_exact(&mut b)?;
        Ok(u16::from_ne_bytes(b))
    }

    pub fn read_i64(&mut self) -> IoResult<i64> {
        let mut b = [0u8; 8];
        self.read_exact(&mut b)?;
        Ok(i64::from_ne_bytes(b))
    }

    pub fn read_u64(&mut self) -> IoResult<u64> {
        let mut b = [0u8; 8];
        self.read_exact(&mut b)?;
        Ok(u64::from_ne_bytes(b))
    }

    pub fn read_f64(&mut self) -> IoResult<f64> {
        let mut b = [0u8; 8];
        self.read_exact(&mut b)?;
        Ok(f64::from_ne_bytes(b))
    }

    pub fn read_value<T: FilePrimitive>(&mut self) -> IoResult<T> {
        let mut b = vec![0u8; T::SIZE];
        self.read_exact(&mut b)?;
        T::from_ne_bytes_slice(&b)
    }

    pub fn read_values<T: FilePrimitive>(&mut self, ptr: &mut [T]) -> IoResult<()> {
        for x in ptr {
            *x = self.read_value()?;
        }
        Ok(())
    }

    pub fn read_string(&mut self) -> IoResult<String> {
        let mut out = Vec::new();
        if !self.read_to(&mut out, 0)? {
            return Err(IoError::EndOfStream);
        }
        String::from_utf8(out).map_err(|e| IoError::Other(e.to_string()))
    }

    pub fn read_exact(&mut self, ptr: &mut [u8]) -> IoResult<()> {
        let n = self.read_raw(ptr)?;
        if n != ptr.len() {
            Err(IoError::EndOfStream)
        } else {
            Ok(())
        }
    }

    pub fn read_raw(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        let mut total = 0usize;
        while total < ptr.len() {
            let n = self.buffer.read(&mut ptr[total..])?;
            if n == 0 {
                break;
            }
            total += n;
        }
        Ok(total)
    }

    pub fn pop(&mut self, dst: &mut [u8]) -> IoResult<()> {
        self.read_exact(dst)
    }

    pub fn read_to(&mut self, dst: &mut Vec<u8>, delimiter: u8) -> IoResult<bool> {
        let mut b = [0u8; 1];
        loop {
            match self.buffer.read(&mut b)? {
                0 => return Ok(false),
                1 if b[0] == delimiter => return Ok(true),
                _ => dst.push(b[0]),
            }
        }
    }

    pub fn data(&self) -> &[u8] {
        self.buffer.data()
    }

    pub fn read_to_record_start(
        &mut self,
        dst: &mut Vec<u8>,
        line_delimiter: u8,
        record_start: u8,
    ) -> IoResult<(bool, i64)> {
        let mut nl = 0i64;
        let mut b = [0u8; 1];
        loop {
            match self.buffer.read(&mut b)? {
                0 => return Ok((false, nl)),
                _ if b[0] == line_delimiter => {
                    nl += 1;
                    let pos = self.buffer.tell()?;
                    let mut next = [0u8; 1];
                    match self.buffer.read(&mut next)? {
                        0 => return Ok((false, nl)),
                        _ if next[0] == record_start => {
                            self.buffer.seek(pos, SeekFrom::Start(0))?;
                            return Ok((true, nl));
                        }
                        _ => {
                            dst.push(b[0]);
                            dst.push(next[0]);
                        }
                    }
                }
                _ => dst.push(b[0]),
            }
        }
    }

    pub fn peek(&mut self, n: usize) -> IoResult<String> {
        let data = self.buffer.data();
        if data.len() >= n {
            return String::from_utf8(data[..n].to_vec())
                .map_err(|e| IoError::Other(e.to_string()));
        }
        let pos = self.buffer.tell()?;
        let mut buf = vec![0u8; n];
        let read = self.buffer.read(&mut buf)?;
        self.buffer.seek(pos, SeekFrom::Start(0))?;
        buf.truncate(read);
        String::from_utf8(buf).map_err(|e| IoError::Other(e.to_string()))
    }

    pub fn file_size(&mut self) -> IoResult<i64> {
        self.buffer.file_size()
    }

    pub fn file(&mut self) -> IoResult<&mut StdFile> {
        self.buffer.file()
    }

    pub fn into_inner(self) -> S {
        self.buffer
    }
}

#[derive(Debug)]
pub struct OutputFile {
    serializer: Serializer<OutputStreamBuffer<FileSink>>,
    file_name: String,
    compressor: Compressor,
}

impl OutputFile {
    pub fn new(file_name: &str, compressor: Compressor, mode: &str) -> IoResult<Self> {
        Ok(Self {
            serializer: Serializer::new(OutputStreamBuffer::new(FileSink::new(
                file_name, mode, false, 0,
            )?)),
            file_name: file_name.to_string(),
            compressor,
        })
    }

    pub fn from_temp_file_data(
        data: &TempFileData,
        compressor: Compressor,
        mode: &str,
    ) -> IoResult<Self> {
        #[cfg(unix)]
        {
            let fd = unsafe { dup(data.fd) };
            if fd < 0 {
                return Err(IoError::Other(format!(
                    "Error opening temporary file {}",
                    data.name
                )));
            }
            return Ok(Self {
                serializer: Serializer::new(OutputStreamBuffer::new(unsafe {
                    FileSink::from_fd(&data.name, fd, mode, false, 0)?
                })),
                file_name: data.name.clone(),
                compressor,
            });
        }
        #[cfg(not(unix))]
        {
            Self::new(&data.name, compressor, mode)
        }
    }

    pub fn remove(&mut self) {
        if std::fs::remove_file(&self.file_name).is_err() {
            eprintln!("Warning: Failed to delete file {}", self.file_name);
        }
    }

    pub fn advise_need(&mut self) {}

    pub fn file_name(&self) -> String {
        self.file_name.clone()
    }

    pub fn write_raw(&mut self, ptr: &[u8]) -> IoResult<()> {
        self.serializer.write_raw(ptr)
    }

    pub fn write_i32(&mut self, x: i32) -> IoResult<&mut Self> {
        self.serializer.write_i32(x)?;
        Ok(self)
    }

    pub fn write_i16(&mut self, x: i16) -> IoResult<&mut Self> {
        self.serializer.write_i16(x)?;
        Ok(self)
    }

    pub fn write_u16(&mut self, x: u16) -> IoResult<&mut Self> {
        self.serializer.write_u16(x)?;
        Ok(self)
    }

    pub fn write_i64(&mut self, x: i64) -> IoResult<&mut Self> {
        self.serializer.write_i64(x)?;
        Ok(self)
    }

    pub fn write_u32(&mut self, x: u32) -> IoResult<&mut Self> {
        self.serializer.write_u32(x)?;
        Ok(self)
    }

    pub fn write_u64(&mut self, x: u64) -> IoResult<&mut Self> {
        self.serializer.write_u64(x)?;
        Ok(self)
    }

    pub fn write_f64(&mut self, x: f64) -> IoResult<&mut Self> {
        self.serializer.write_f64(x)?;
        Ok(self)
    }

    pub fn write_value<T: FilePrimitive>(&mut self, x: T) -> IoResult<&mut Self> {
        self.serializer.write_value(x)?;
        Ok(self)
    }

    pub fn write_values<T: FilePrimitive>(&mut self, ptr: &[T]) -> IoResult<&mut Self> {
        self.serializer.write_values(ptr)?;
        Ok(self)
    }

    pub fn write_string(&mut self, s: &str) -> IoResult<&mut Self> {
        self.serializer.write_string(s)?;
        Ok(self)
    }

    pub fn file_size(&mut self) -> IoResult<i64> {
        self.serializer.file_size()
    }

    pub fn seek(&mut self, p: i64, origin: SeekFrom) -> IoResult<()> {
        self.serializer.seek(p, origin)
    }

    pub fn rewind(&mut self) -> IoResult<()> {
        self.serializer.rewind()
    }

    pub fn tell(&mut self) -> IoResult<usize> {
        self.serializer.tell()
    }

    pub fn close(&mut self) -> IoResult<()> {
        self.compress_pending_output()?;
        self.serializer.close()
    }

    pub fn flush(&mut self) -> IoResult<()> {
        self.serializer.flush()
    }

    pub fn file(&mut self) -> IoResult<&mut StdFile> {
        self.serializer.file()
    }

    fn compress_pending_output(&mut self) -> IoResult<()> {
        if self.compressor == Compressor::None {
            return Ok(());
        }
        self.serializer.flush()?;
        let file = self.serializer.file()?;
        file.seek(SeekFrom::Start(0))
            .map_err(|_| IoError::Other("Error calling fseek.".to_string()))?;
        let mut plain = Vec::new();
        file.read_to_end(&mut plain)
            .map_err(|_| IoError::FileRead(self.file_name.clone()))?;
        let compressed = match self.compressor {
            Compressor::None => plain,
            Compressor::Zlib => {
                let mut encoder =
                    flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
                encoder
                    .write_all(&plain)
                    .map_err(|e| IoError::Other(format!("deflate error: {e}")))?;
                encoder
                    .finish()
                    .map_err(|e| IoError::Other(format!("deflate error: {e}")))?
            }
            Compressor::Zstd => zstd::stream::encode_all(&plain[..], 0)
                .map_err(|e| IoError::Other(format!("ZSTD_endStream: {e}")))?,
        };
        file.set_len(0)
            .map_err(|_| IoError::FileWrite(self.file_name.clone()))?;
        file.seek(SeekFrom::Start(0))
            .map_err(|_| IoError::Other("Error calling fseek.".to_string()))?;
        file.write_all(&compressed)
            .map_err(|_| IoError::FileWrite(self.file_name.clone()))?;
        file.flush()
            .map_err(|_| IoError::FileWrite(self.file_name.clone()))?;
        self.compressor = Compressor::None;
        Ok(())
    }
}

#[derive(Debug)]
pub struct TempFile {
    output: OutputFile,
    unlinked: bool,
}

impl TempFile {
    pub fn new(unlink: bool) -> IoResult<Self> {
        Self::from_temp_file_data(&TempFileData::init(unlink)?)
    }

    pub fn from_file_name(file_name: &str) -> IoResult<Self> {
        Ok(Self {
            output: OutputFile::new(file_name, Compressor::None, "w+b")?,
            unlinked: false,
        })
    }

    pub fn from_temp_file_data(data: &TempFileData) -> IoResult<Self> {
        Ok(Self {
            output: OutputFile::from_temp_file_data(data, Compressor::None, "w+b")?,
            unlinked: data.unlinked,
        })
    }

    pub fn get_temp_dir() -> IoResult<String> {
        let tmp = Self::new(false)?;
        let dir = std::path::Path::new(&tmp.output.file_name)
            .parent()
            .map(|p| p.to_string_lossy().into_owned())
            .unwrap_or_default();
        let _ = std::fs::remove_file(&tmp.output.file_name);
        Ok(dir)
    }

    pub fn output_file(&mut self) -> &mut OutputFile {
        &mut self.output
    }

    pub fn file_name(&self) -> String {
        self.output.file_name()
    }

    pub fn unlinked(&self) -> bool {
        self.unlinked
    }
}

impl std::ops::Deref for TempFile {
    type Target = OutputFile;

    fn deref(&self) -> &Self::Target {
        &self.output
    }
}

impl std::ops::DerefMut for TempFile {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.output
    }
}

#[derive(Debug)]
pub struct InputFile {
    deserializer: Deserializer<InputFileStream>,
    pub file_name: String,
    pub unlinked: bool,
    pub temp_file: bool,
}

#[derive(Debug)]
pub enum InputFileStream {
    File(InputStreamBuffer<FileSource>),
    Memory(VecStream),
}

impl StreamEntity for InputFileStream {
    fn rewind(&mut self) -> IoResult<()> {
        match self {
            Self::File(stream) => stream.rewind(),
            Self::Memory(stream) => stream.rewind(),
        }
    }

    fn seek(&mut self, p: i64, origin: SeekFrom) -> IoResult<()> {
        match self {
            Self::File(stream) => stream.seek(p, origin),
            Self::Memory(stream) => stream.seek(p, origin),
        }
    }

    fn tell(&mut self) -> IoResult<i64> {
        match self {
            Self::File(stream) => stream.tell(),
            Self::Memory(stream) => stream.tell(),
        }
    }

    fn read(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        match self {
            Self::File(stream) => stream.read(ptr),
            Self::Memory(stream) => stream.read(ptr),
        }
    }

    fn fetch(&mut self) -> IoResult<bool> {
        match self {
            Self::File(stream) => stream.fetch(),
            Self::Memory(stream) => stream.fetch(),
        }
    }

    fn close(&mut self) -> IoResult<()> {
        match self {
            Self::File(stream) => stream.close(),
            Self::Memory(stream) => stream.close(),
        }
    }

    fn file_name(&self) -> &str {
        match self {
            Self::File(stream) => stream.file_name(),
            Self::Memory(stream) => stream.file_name(),
        }
    }

    fn file(&mut self) -> IoResult<&mut StdFile> {
        match self {
            Self::File(stream) => stream.file(),
            Self::Memory(stream) => stream.file(),
        }
    }

    fn data(&self) -> &[u8] {
        match self {
            Self::File(stream) => stream.data(),
            Self::Memory(stream) => stream.data(),
        }
    }

    fn eof(&mut self) -> IoResult<bool> {
        match self {
            Self::File(stream) => stream.eof(),
            Self::Memory(stream) => stream.eof(),
        }
    }

    fn file_size(&mut self) -> IoResult<i64> {
        match self {
            Self::File(stream) => stream.file_size(),
            Self::Memory(stream) => stream.file_size(),
        }
    }

    fn seekable(&self) -> bool {
        match self {
            Self::File(stream) => stream.seekable(),
            Self::Memory(stream) => stream.seekable(),
        }
    }
}

impl InputFile {
    pub const BUFFERED: i32 = 1;
    pub const NO_AUTODETECT: i32 = 2;

    pub fn new(file_name: &str, flags: i32) -> IoResult<Self> {
        let mut stream = InputStreamBuffer::new(FileSource::new(file_name)?, flags);
        if !file_name.is_empty()
            && file_name != "-"
            && (flags & Self::NO_AUTODETECT) == 0
            && stream.seekable()
        {
            let mut begin = [0u8; 4];
            let read = stream.read(&mut begin)?;
            stream.seek(0, SeekFrom::Start(0))?;
            let compressor = if read < 4 {
                Compressor::None
            } else {
                detect_compressor(&begin)
            };
            if compressor != Compressor::None {
                let compressed = std::fs::read(file_name)
                    .map_err(|e| IoError::Other(format!("Error reading file {file_name}. {e}")))?;
                let mut decoded = Vec::new();
                match compressor {
                    Compressor::Zlib => {
                        let mut capacity = compressed.len().saturating_mul(4).max(1024);
                        loop {
                            decoded.resize(capacity, 0);
                            match zlib_decompress(&compressed, &mut decoded) {
                                Ok(n) => {
                                    decoded.truncate(n);
                                    break;
                                }
                                Err(IoError::Other(msg))
                                    if msg.contains("output buffer too small") =>
                                {
                                    capacity = capacity.saturating_mul(2);
                                }
                                Err(e) => return Err(e),
                            }
                        }
                    }
                    Compressor::Zstd => {
                        decoded = zstd::stream::decode_all(std::io::Cursor::new(&compressed))
                            .map_err(|e| {
                                IoError::Other(format!("Failed decompressing zstd stream: {e}"))
                            })?;
                    }
                    Compressor::None => {}
                }
                return Ok(Self {
                    deserializer: Deserializer::new(InputFileStream::Memory(VecStream::from_vec(
                        decoded,
                    ))),
                    file_name: file_name.to_string(),
                    unlinked: false,
                    temp_file: false,
                });
            }
        }
        Ok(Self {
            deserializer: Deserializer::new(InputFileStream::File(stream)),
            file_name: file_name.to_string(),
            unlinked: false,
            temp_file: false,
        })
    }

    pub fn from_temp_file_data(
        data: &TempFileData,
        flags: i32,
        compressor: Compressor,
    ) -> IoResult<Self> {
        #[cfg(unix)]
        if compressor != Compressor::None {
            let fd = unsafe { dup(data.fd) };
            if fd < 0 {
                return Err(IoError::Other(format!(
                    "Error opening temporary file {}",
                    data.name
                )));
            }
            let mut file = unsafe { StdFile::from_raw_fd(fd) };
            file.seek(SeekFrom::Start(0))
                .map_err(|_| IoError::Other("Error calling fseek.".to_string()))?;
            let mut compressed = Vec::new();
            file.read_to_end(&mut compressed)
                .map_err(|_| IoError::FileRead(data.name.clone()))?;
            let mut decoded = Vec::new();
            match compressor {
                Compressor::Zlib => {
                    let mut capacity = compressed.len().saturating_mul(4).max(1024);
                    loop {
                        decoded.resize(capacity, 0);
                        match zlib_decompress(&compressed, &mut decoded) {
                            Ok(n) => {
                                decoded.truncate(n);
                                break;
                            }
                            Err(IoError::Other(msg)) if msg.contains("output buffer too small") => {
                                capacity = capacity.saturating_mul(2);
                            }
                            Err(e) => return Err(e),
                        }
                    }
                }
                Compressor::Zstd => {
                    decoded = zstd::stream::decode_all(std::io::Cursor::new(&compressed)).map_err(
                        |e| IoError::Other(format!("Failed decompressing zstd stream: {e}")),
                    )?;
                }
                Compressor::None => {}
            }
            return Ok(Self {
                deserializer: Deserializer::new(InputFileStream::Memory(VecStream::from_vec(
                    decoded,
                ))),
                file_name: data.name.clone(),
                unlinked: data.unlinked,
                temp_file: true,
            });
        }

        #[cfg(not(unix))]
        if compressor != Compressor::None {
            let mut compressed = std::fs::read(&data.name)
                .map_err(|e| IoError::Other(format!("Error reading file {}. {e}", data.name)))?;
            let mut decoded = Vec::new();
            match compressor {
                Compressor::Zlib => {
                    let mut capacity = compressed.len().saturating_mul(4).max(1024);
                    loop {
                        decoded.resize(capacity, 0);
                        match zlib_decompress(&compressed, &mut decoded) {
                            Ok(n) => {
                                decoded.truncate(n);
                                break;
                            }
                            Err(IoError::Other(msg)) if msg.contains("output buffer too small") => {
                                capacity = capacity.saturating_mul(2);
                            }
                            Err(e) => return Err(e),
                        }
                    }
                }
                Compressor::Zstd => {
                    decoded = zstd::stream::decode_all(std::io::Cursor::new(&compressed)).map_err(
                        |e| IoError::Other(format!("Failed decompressing zstd stream: {e}")),
                    )?;
                }
                Compressor::None => {}
            }
            compressed.clear();
            return Ok(Self {
                deserializer: Deserializer::new(InputFileStream::Memory(VecStream::from_vec(
                    decoded,
                ))),
                file_name: data.name.clone(),
                unlinked: data.unlinked,
                temp_file: true,
            });
        }
        Ok(Self {
            deserializer: Deserializer::new(InputFileStream::File(InputStreamBuffer::new(
                FileSource::new(&data.name)?,
                flags,
            ))),
            file_name: data.name.clone(),
            unlinked: data.unlinked,
            temp_file: true,
        })
    }

    pub fn from_temp_file(
        tmp_file: &mut TempFile,
        flags: i32,
        compressor: Compressor,
    ) -> IoResult<Self> {
        tmp_file.rewind()?;
        let file_name = tmp_file.file_name();
        let file = tmp_file
            .file()?
            .try_clone()
            .map_err(|_| IoError::FileOpen(file_name.clone()))?;
        if compressor != Compressor::None {
            let mut file = file;
            file.seek(SeekFrom::Start(0))
                .map_err(|_| IoError::Other("Error calling fseek.".to_string()))?;
            let mut compressed = Vec::new();
            file.read_to_end(&mut compressed)
                .map_err(|_| IoError::FileRead(file_name.clone()))?;
            let mut decoded = Vec::new();
            match compressor {
                Compressor::Zlib => {
                    let mut capacity = compressed.len().saturating_mul(4).max(1024);
                    loop {
                        decoded.resize(capacity, 0);
                        match zlib_decompress(&compressed, &mut decoded) {
                            Ok(n) => {
                                decoded.truncate(n);
                                break;
                            }
                            Err(IoError::Other(msg)) if msg.contains("output buffer too small") => {
                                capacity = capacity.saturating_mul(2);
                            }
                            Err(e) => return Err(e),
                        }
                    }
                }
                Compressor::Zstd => {
                    decoded = zstd::stream::decode_all(std::io::Cursor::new(&compressed)).map_err(
                        |e| IoError::Other(format!("Failed decompressing zstd stream: {e}")),
                    )?;
                }
                Compressor::None => {}
            }
            return Ok(Self {
                deserializer: Deserializer::new(InputFileStream::Memory(VecStream::from_vec(
                    decoded,
                ))),
                file_name,
                unlinked: tmp_file.unlinked(),
                temp_file: true,
            });
        }
        Ok(Self {
            deserializer: Deserializer::new(InputFileStream::File(InputStreamBuffer::new(
                FileSource::from_file(&file_name, file),
                flags,
            ))),
            file_name,
            unlinked: tmp_file.unlinked(),
            temp_file: true,
        })
    }

    pub fn from_output_file(tmp_file: &mut OutputFile, flags: i32) -> IoResult<Self> {
        tmp_file.rewind()?;
        let file_name = tmp_file.file_name();
        let file = tmp_file
            .file()?
            .try_clone()
            .map_err(|_| IoError::FileOpen(file_name.clone()))?;
        Ok(Self {
            deserializer: Deserializer::new(InputFileStream::File(InputStreamBuffer::new(
                FileSource::from_file(&file_name, file),
                flags,
            ))),
            file_name,
            unlinked: false,
            temp_file: true,
        })
    }

    pub fn close_and_delete(&mut self) -> IoResult<()> {
        self.close()?;
        if !self.unlinked && std::fs::remove_file(&self.file_name).is_err() {
            eprintln!(
                "Warning: Failed to delete temporary file {}",
                self.file_name
            );
        }
        Ok(())
    }

    pub fn hash(&mut self) -> IoResult<u64> {
        let mut seed = [0u8; 16];
        let mut buf = [0u8; 4096];
        loop {
            let n = self.read_raw(&mut buf)?;
            if n == 0 {
                break;
            }
            seed = crate::util::hash::murmurhash3_x64_128(&buf[..n], &seed);
        }
        Ok(u64::from_le_bytes(seed[0..8].try_into().unwrap()))
    }

    pub fn rewind(&mut self) -> IoResult<()> {
        self.deserializer.rewind()
    }

    pub fn seek(&mut self, pos: i64) -> IoResult<&mut Self> {
        self.deserializer.seek(pos)?;
        Ok(self)
    }

    pub fn seek_forward(&mut self, n: usize) -> IoResult<()> {
        self.deserializer.seek_forward(n)
    }

    pub fn seek_forward_delim(&mut self, delimiter: u8) -> IoResult<bool> {
        self.deserializer.seek_forward_delim(delimiter)
    }

    pub fn close(&mut self) -> IoResult<()> {
        self.deserializer.close()
    }

    pub fn read_u32(&mut self) -> IoResult<u32> {
        self.deserializer.read_u32()
    }

    pub fn read_i32(&mut self) -> IoResult<i32> {
        self.deserializer.read_i32()
    }

    pub fn read_i16(&mut self) -> IoResult<i16> {
        self.deserializer.read_i16()
    }

    pub fn read_u16(&mut self) -> IoResult<u16> {
        self.deserializer.read_u16()
    }

    pub fn read_i64(&mut self) -> IoResult<i64> {
        self.deserializer.read_i64()
    }

    pub fn read_u64(&mut self) -> IoResult<u64> {
        self.deserializer.read_u64()
    }

    pub fn read_f64(&mut self) -> IoResult<f64> {
        self.deserializer.read_f64()
    }

    pub fn read_value<T: FilePrimitive>(&mut self) -> IoResult<T> {
        self.deserializer.read_value()
    }

    pub fn read_values<T: FilePrimitive>(&mut self, ptr: &mut [T]) -> IoResult<()> {
        self.deserializer.read_values(ptr)
    }

    pub fn read_string(&mut self) -> IoResult<String> {
        self.deserializer.read_string()
    }

    pub fn read_raw(&mut self, ptr: &mut [u8]) -> IoResult<usize> {
        self.deserializer.read_raw(ptr)
    }

    pub fn read_to(&mut self, dst: &mut Vec<u8>, delimiter: u8) -> IoResult<bool> {
        self.deserializer.read_to(dst, delimiter)
    }

    pub fn peek(&mut self, n: usize) -> IoResult<String> {
        self.deserializer.peek(n)
    }

    pub fn file_size(&mut self) -> IoResult<i64> {
        self.deserializer.file_size()
    }

    pub fn file(&mut self) -> IoResult<&mut StdFile> {
        self.deserializer.file()
    }
}

#[derive(Debug, Clone)]
pub struct TextInputFile<S: StreamEntity> {
    input: Deserializer<S>,
    pub line: String,
    pub line_count: usize,
    putback_line: bool,
    eof: bool,
    line_separator: u8,
}

impl<S: StreamEntity> TextInputFile<S> {
    pub fn new(input: Deserializer<S>, line_separator: u8) -> Self {
        Self {
            input,
            line: String::new(),
            line_count: 0,
            putback_line: false,
            eof: false,
            line_separator,
        }
    }

    pub fn rewind(&mut self) -> IoResult<()> {
        self.input.rewind()?;
        self.line_count = 0;
        self.putback_line = false;
        self.eof = false;
        self.line.clear();
        Ok(())
    }

    pub fn eof(&self) -> bool {
        self.eof
    }

    pub fn getline(&mut self) -> IoResult<()> {
        if self.putback_line {
            self.putback_line = false;
            self.line_count += 1;
            return Ok(());
        }
        let mut line = Vec::new();
        self.eof = !self.input.read_to(&mut line, self.line_separator)?;
        self.line_count += 1;
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        self.line = String::from_utf8(line).map_err(|e| IoError::Other(e.to_string()))?;
        Ok(())
    }

    pub fn putback_line(&mut self) {
        self.putback_line = true;
        self.line_count = self.line_count.saturating_sub(1);
    }

    pub fn is_open(&self) -> bool {
        !self.eof()
    }

    pub fn into_inner(self) -> Deserializer<S> {
        self.input
    }
}

impl TextInputFile<InputFileStream> {
    pub fn from_file_name(file_name: &str, line_separator: u8) -> IoResult<Self> {
        let input = InputFile::new(file_name, 0)?;
        Ok(Self::new(input.deserializer, line_separator))
    }

    pub fn from_temp_file(tmp_file: &mut TempFile, line_separator: u8) -> IoResult<Self> {
        let input = InputFile::from_temp_file(tmp_file, 0, Compressor::None)?;
        Ok(Self::new(input.deserializer, line_separator))
    }

    pub fn from_output_file(out_file: &mut OutputFile, line_separator: u8) -> IoResult<Self> {
        let input = InputFile::from_output_file(out_file, 0)?;
        Ok(Self::new(input.deserializer, line_separator))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_serializer_deserializer_roundtrip() {
        let stream = VecStream::new();
        let mut ser = Serializer::new(stream);
        ser.write_i16(-2).unwrap();
        ser.write_u16(0x0102).unwrap();
        ser.write_i32(-7).unwrap();
        ser.write_u64(0x0102_0304_0506_0708).unwrap();
        ser.write_f64(1.5).unwrap();
        ser.write_value(-9i32).unwrap();
        ser.write_values(&[11u16, 12u16, 13u16]).unwrap();
        ser.write_string("abc").unwrap();
        let mut stream = ser.into_inner().unwrap();
        stream.rewind().unwrap();

        let mut de = Deserializer::new(stream);
        assert_eq!(de.read_i16().unwrap(), -2);
        assert_eq!(de.read_u16().unwrap(), 0x0102);
        assert_eq!(de.read_i32().unwrap(), -7);
        assert_eq!(de.read_u64().unwrap(), 0x0102_0304_0506_0708);
        assert_eq!(de.read_f64().unwrap(), 1.5);
        assert_eq!(de.read_value::<i32>().unwrap(), -9);
        let mut values = [0u16; 3];
        de.read_values(&mut values).unwrap();
        assert_eq!(values, [11, 12, 13]);
        assert_eq!(de.read_string().unwrap(), "abc");
    }

    #[test]
    fn test_serializer_buffered_write_seek_and_reset() {
        let stream = VecStream::new();
        let mut ser = Serializer::new(stream);
        let data = vec![b'a'; DEFAULT_FILE_BUFFER_SIZE + 17];
        ser.write_raw(&data).unwrap();
        assert_eq!(ser.tell().unwrap(), DEFAULT_FILE_BUFFER_SIZE + 17);
        ser.seek(4, SeekFrom::Start(0)).unwrap();
        ser.write_raw(b"BC").unwrap();
        let stream = ser.into_inner().unwrap();
        assert_eq!(&stream.data()[..4], b"aaaa");
        assert_eq!(&stream.data()[4..6], b"BC");
        assert_eq!(stream.data().len(), DEFAULT_FILE_BUFFER_SIZE + 17);

        let stream = VecStream::new();
        let mut ser = Serializer::new(stream);
        ser.reset_buffer();
        ser.write_raw(b"xy").unwrap();
        assert_eq!(ser.into_inner().unwrap().data(), b"xy");
    }

    #[test]
    fn test_deserializer_read_to_seek_and_peek() {
        let stream = VecStream::from_vec(b"abc,def\n>rec\nx".to_vec());
        let mut de = Deserializer::new(stream);
        assert!(de.data().is_empty());
        assert_eq!(de.peek(3).unwrap(), "abc");
        let mut out = Vec::new();
        assert!(de.read_to(&mut out, b',').unwrap());
        assert_eq!(out, b"abc");
        de.seek_forward(3).unwrap();
        let mut out = Vec::new();
        assert_eq!(
            de.read_to_record_start(&mut out, b'\n', b'>').unwrap(),
            (true, 1)
        );
        assert!(out.is_empty());
        assert_eq!(de.read_string().unwrap_err(), IoError::EndOfStream);

        let mut input = InputStreamBuffer::new(VecStream::from_vec(b"peek-pop".to_vec()), 0);
        input.fetch().unwrap();
        let mut de = Deserializer::new(input);
        assert_eq!(de.peek(4).unwrap(), "peek");
        let mut popped = [0u8; 4];
        de.pop(&mut popped).unwrap();
        assert_eq!(&popped, b"peek");
    }

    #[test]
    fn test_file_temporary_write_read_size() {
        let mut file = File::new_temporary(Temporary).unwrap();
        file.write(b"abcd").unwrap();
        file.write_value(7u16).unwrap();
        assert_eq!(file.size().unwrap(), 6);
        file.seek(0, SeekFrom::Start(0)).unwrap();
        let mut buf = [0u8; 4];
        file.read_exact(&mut buf).unwrap();
        assert_eq!(&buf, b"abcd");
        assert!(!file.eof().unwrap());
        assert_eq!(file.read(2).unwrap(), &7u16.to_ne_bytes());
        assert!(file.eof().unwrap());
        file.seek(4, SeekFrom::Start(0)).unwrap();
        assert_eq!(file.read_value::<u16>().unwrap(), 7);
        assert!(!file.file_name().is_empty());
        #[cfg(unix)]
        assert!(!std::path::Path::new(file.file_name()).exists());
        assert!(file.file().is_ok());
        file.close().unwrap();

        let path = std::env::temp_dir().join(format!(
            "diamond-rs-file-new-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut named = File::new(&name, "w+b").unwrap();
            named.write(b"x").unwrap();
            assert_eq!(named.file_name(), name);
        }
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_stream_entity_file_delegation() {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "diamond-rs-stream-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut serializer = Serializer::new(OutputStreamBuffer::new(
                FileSink::new(&name, "w+b", false, 0).unwrap(),
            ));
            serializer.write_raw(b"abcdef").unwrap();
            assert!(serializer.file().is_ok());
            serializer.close().unwrap();
        }
        {
            let mut de =
                Deserializer::new(InputStreamBuffer::new(FileSource::new(&name).unwrap(), 0));
            assert!(de.file().is_ok());
            let mut b = [0u8; 2];
            de.read_exact(&mut b).unwrap();
            assert_eq!(&b, b"ab");
            de.close().unwrap();
        }
        let _ = std::fs::remove_file(&name);
    }

    #[test]
    fn test_errors_display() {
        assert_eq!(
            IoError::UnsupportedOperation.to_string(),
            "Unsupported I/O operation."
        );
        assert_eq!(IoError::EndOfStream.to_string(), "Unexpected end of input.");
        assert!(IoError::StreamRead {
            line_count: 3,
            msg: "bad".to_string()
        }
        .to_string()
        .contains("line 3"));
    }

    #[test]
    fn test_detect_compressor_and_decompress() {
        assert_eq!(detect_compressor(&[0x1f, 0x8b, 0, 0]), Compressor::Zlib);
        assert_eq!(detect_compressor(&[0x78, 0x9c, 0, 0]), Compressor::Zlib);
        assert_eq!(
            detect_compressor(&[0x28, 0xb5, 0x2f, 0xfd]),
            Compressor::Zstd
        );
        assert_eq!(detect_compressor(b"plain"), Compressor::None);

        let mut zlib_encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        zlib_encoder.write_all(b"abcabc").unwrap();
        let compressed = zlib_encoder.finish().unwrap();
        let mut out = [0u8; 16];
        let n = decompress(&compressed, &mut out, Compressor::Zlib).unwrap();
        assert_eq!(&out[..n], b"abcabc");

        let mut gzip_encoder =
            flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip_encoder.write_all(b"xyz").unwrap();
        let compressed = gzip_encoder.finish().unwrap();
        let n = decompress(&compressed, &mut out, Compressor::Zlib).unwrap();
        assert_eq!(&out[..n], b"xyz");
        let n = zlib_decompress(&compressed, &mut out).unwrap();
        assert_eq!(&out[..n], b"xyz");

        let compressed = zstd::stream::encode_all(&b"zstd-data"[..], 0).unwrap();
        let n = decompress(&compressed, &mut out, Compressor::Zstd).unwrap();
        assert_eq!(&out[..n], b"zstd-data");
        let n = zstd_decompress(&compressed, &mut out).unwrap();
        assert_eq!(&out[..n], b"zstd-data");
        assert!(zstd_decompress(&compressed, &mut [0u8; 2]).is_err());

        assert!(decompress(&compressed, &mut [0u8; 2], Compressor::Zlib).is_err());

        let mut first = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        first.write_all(b"one").unwrap();
        let mut concatenated = first.finish().unwrap();
        let mut second =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        second.write_all(b"two").unwrap();
        concatenated.extend_from_slice(&second.finish().unwrap());
        let mut out = [0u8; 6];
        let n = zlib_decompress(&concatenated, &mut out).unwrap();
        assert_eq!(&out[..n], b"onetwo");
    }

    #[test]
    fn test_compressed_buffer_roundtrip_and_clear() {
        let mut compressed = CompressedBuffer::new();
        compressed.write(b"abc").unwrap();
        compressed.write_value(7u16).unwrap();
        compressed.write(b"xyzxyzxyz").unwrap();
        compressed.finish().unwrap();
        assert!(compressed.size() > 0);
        assert_eq!(detect_compressor(compressed.data()), Compressor::Zlib);

        let mut out = [0u8; 32];
        let n = decompress(compressed.data(), &mut out, Compressor::Zlib).unwrap();
        let mut expected = Vec::new();
        expected.extend_from_slice(b"abc");
        expected.extend_from_slice(&7u16.to_ne_bytes());
        expected.extend_from_slice(b"xyzxyzxyz");
        assert_eq!(&out[..n], expected.as_slice());

        compressed.clear();
        assert_eq!(compressed.size(), 0);
        compressed.write(b"second").unwrap();
        compressed.finish().unwrap();
        let n = decompress(compressed.data(), &mut out, Compressor::Zlib).unwrap();
        assert_eq!(&out[..n], b"second");
    }

    #[test]
    fn test_file_source_sink_and_stream_buffers() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let mut out = OutputStreamBuffer::new(sink);
            out.write(b"abc").unwrap();
            out.write(b"def").unwrap();
            assert_eq!(out.tell().unwrap(), 6);
            out.close().unwrap();
        }
        {
            let source = FileSource::new(&name).unwrap();
            assert!(source.seekable());
            let mut input = InputStreamBuffer::new(source, 0);
            let mut buf = [0u8; 6];
            assert_eq!(input.read(&mut buf).unwrap(), 6);
            assert_eq!(&buf, b"abcdef");
            assert_eq!(input.file_size().unwrap(), 6);
        }
        let _ = std::fs::remove_file(path);
    }

    #[cfg(unix)]
    #[test]
    fn test_file_source_sink_standard_stream_names() {
        let stdin_empty = FileSource::new("").unwrap();
        assert_eq!(stdin_empty.file_name(), "");
        assert!(!stdin_empty.seekable());

        let stdin_dash = FileSource::new("-").unwrap();
        assert_eq!(stdin_dash.file_name(), "-");
        assert!(!stdin_dash.seekable());

        let mut stdout_sink = FileSink::new("", "wb", false, 0).unwrap();
        assert_eq!(stdout_sink.file_name(), "");
        stdout_sink.close().unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn test_file_sink_posix_flags_and_from_fd() {
        use std::os::fd::IntoRawFd;

        assert_eq!(posix_flags("wb").unwrap(), 1 | 64 | 512);
        assert_eq!(posix_flags("r+b").unwrap(), 2);
        assert_eq!(posix_flags("w+b").unwrap(), 2 | 64 | 512);
        assert_eq!(
            posix_flags("bad").unwrap_err().to_string(),
            "Invalid fopen mode."
        );

        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-fd-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        let fd = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&path)
            .unwrap()
            .into_raw_fd();
        let mut sink = unsafe { FileSink::from_fd(&name, fd, "w+b", false, 0).unwrap() };
        sink.write(b"fd").unwrap();
        sink.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"fd");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_output_stream_buffer_write_buffer_and_flush_count() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-output-buffer-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let mut out = OutputStreamBuffer::new(sink);
            let (buf, size) = out.write_buffer_range();
            assert_eq!(size, DEFAULT_FILE_BUFFER_SIZE);
            buf[..4].copy_from_slice(b"abcd");
            out.flush_count(4).unwrap();
            out.write(b"ef").unwrap();
            out.close().unwrap();
        }
        assert_eq!(std::fs::read(&path).unwrap(), b"abcdef");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_compressed_stream_entities_roundtrip() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-compressed-stream-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();

        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let out = OutputStreamBuffer::new(sink);
            let mut zlib = ZlibSink::new(out);
            zlib.write(b"stream-zlib").unwrap();
            zlib.close().unwrap();
        }
        assert_eq!(
            detect_compressor(&std::fs::read(&path).unwrap()),
            Compressor::Zlib
        );
        {
            let source = FileSource::new(&name).unwrap();
            let input = InputStreamBuffer::new(source, 0);
            let mut zlib = ZlibSource::new(input).unwrap();
            let mut buf = [0u8; 11];
            assert_eq!(zlib.read(&mut buf).unwrap(), 11);
            assert_eq!(&buf, b"stream-zlib");
            let mut tail = [0u8; 1];
            assert_eq!(zlib.read(&mut tail).unwrap(), 0);
            assert!(zlib.eof().unwrap());
            zlib.rewind().unwrap();
            let mut head = [0u8; 6];
            assert_eq!(zlib.read(&mut head).unwrap(), 6);
            assert_eq!(&head, b"stream");
        }

        {
            let sink = FileSink::new(&name, "w+b", false, 0).unwrap();
            let out = OutputStreamBuffer::new(sink);
            let mut zstd = ZstdSink::new(out);
            zstd.write(b"stream-zstd").unwrap();
            zstd.close().unwrap();
        }
        assert_eq!(
            detect_compressor(&std::fs::read(&path).unwrap()),
            Compressor::Zstd
        );
        {
            let source = FileSource::new(&name).unwrap();
            let input = InputStreamBuffer::new(source, 0);
            let mut zstd = ZstdSource::new(input).unwrap();
            let mut buf = [0u8; 11];
            assert_eq!(zstd.read(&mut buf).unwrap(), 11);
            assert_eq!(&buf, b"stream-zstd");
            let mut tail = [0u8; 1];
            assert_eq!(zstd.read(&mut tail).unwrap(), 0);
            assert!(zstd.eof().unwrap());
        }

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_temp_file_data_and_temp_file() {
        let data = TempFileData::init(false).unwrap();
        {
            let mut out = OutputFile::from_temp_file_data(&data, Compressor::None, "w+b").unwrap();
            out.write_raw(b"temp").unwrap();
            out.close().unwrap();
        }
        assert_eq!(std::fs::read(&data.name).unwrap(), b"temp");
        std::fs::remove_file(&data.name).unwrap();

        let dir = TempFile::get_temp_dir().unwrap();
        assert!(!dir.is_empty());

        let mut tmp = TempFile::new(true).unwrap();
        assert!(tmp.unlinked());
        tmp.write_raw(b"x").unwrap();
        tmp.close().unwrap();
    }

    #[test]
    fn test_output_file_input_file_and_hash() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-file-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut out = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
            out.write_i16(-123).unwrap();
            out.write_u16(0xCAFE).unwrap();
            out.write_i32(-7).unwrap();
            out.write_u64(0x0102_0304_0506_0708).unwrap();
            out.write_value(99u32).unwrap();
            out.write_values(&[3i16, 4i16]).unwrap();
            out.write_string("abc").unwrap();
            assert!(out.file().is_ok());
            assert!(out.tell().unwrap() > 0);
            out.close().unwrap();
        }
        {
            let mut input = InputFile::new(&name, 0).unwrap();
            assert_eq!(input.read_i16().unwrap(), -123);
            assert_eq!(input.read_u16().unwrap(), 0xCAFE);
            assert_eq!(input.read_i32().unwrap(), -7);
            assert_eq!(input.read_u64().unwrap(), 0x0102_0304_0506_0708);
            assert_eq!(input.read_value::<u32>().unwrap(), 99);
            let mut values = [0i16; 2];
            input.read_values(&mut values).unwrap();
            assert_eq!(values, [3, 4]);
            assert_eq!(input.read_string().unwrap(), "abc");
            assert!(input.file().is_ok());
            input.rewind().unwrap();
            let h = input.hash().unwrap();
            let bytes = std::fs::read(&name).unwrap();
            assert_eq!(h, crate::util::hash::file_hash(&bytes));
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_output_file_compressors_roundtrip() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-compressor-out-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();

        {
            let mut out = OutputFile::new(&name, Compressor::Zlib, "w+b").unwrap();
            out.write_raw(b"gzip-output").unwrap();
            out.close().unwrap();
        }
        let compressed = std::fs::read(&path).unwrap();
        assert_eq!(detect_compressor(&compressed), Compressor::Zlib);
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 11];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"gzip-output");

        {
            let mut out = OutputFile::new(&name, Compressor::Zstd, "w+b").unwrap();
            out.write_raw(b"zstd-output").unwrap();
            out.close().unwrap();
        }
        let compressed = std::fs::read(&path).unwrap();
        assert_eq!(detect_compressor(&compressed), Compressor::Zstd);
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 11];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zstd-output");

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_input_file_temp_constructors_with_compressors() {
        let mut tmp = TempFile::new(true).unwrap();
        let mut zlib_encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        zlib_encoder.write_all(b"temp-zlib").unwrap();
        tmp.write_raw(&zlib_encoder.finish().unwrap()).unwrap();
        let mut input = InputFile::from_temp_file(&mut tmp, 0, Compressor::Zlib).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"temp-zlib");
        input.close().unwrap();
        tmp.close().unwrap();

        let data = TempFileData::init(false).unwrap();
        {
            let mut out = OutputFile::from_temp_file_data(&data, Compressor::Zstd, "w+b").unwrap();
            out.write_raw(b"temp-zstd").unwrap();
            out.close().unwrap();
        }
        let mut input = InputFile::from_temp_file_data(&data, 0, Compressor::Zstd).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"temp-zstd");
        input.close().unwrap();
        std::fs::remove_file(&data.name).unwrap();
    }

    #[test]
    fn test_input_file_autodetects_zlib_and_gzip() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-compressed-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();

        let mut zlib_encoder =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        zlib_encoder.write_all(b"zlib-data").unwrap();
        std::fs::write(&path, zlib_encoder.finish().unwrap()).unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zlib-data");
        input.rewind().unwrap();
        let mut buf = [0u8; 4];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zlib");

        let mut first = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        first.write_all(b"one").unwrap();
        let mut concatenated = first.finish().unwrap();
        let mut second =
            flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
        second.write_all(b"two").unwrap();
        concatenated.extend_from_slice(&second.finish().unwrap());
        std::fs::write(&path, &concatenated).unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 6];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"onetwo");

        let mut gzip_encoder =
            flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        gzip_encoder.write_all(b"gzip-data").unwrap();
        let gzip_compressed = gzip_encoder.finish().unwrap();
        std::fs::write(&path, &gzip_compressed).unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"gzip-data");

        std::fs::write(
            &path,
            zstd::stream::encode_all(&b"zstd-data"[..], 0).unwrap(),
        )
        .unwrap();
        let mut input = InputFile::new(&name, 0).unwrap();
        let mut buf = [0u8; 9];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"zstd-data");

        std::fs::write(&path, &gzip_compressed).unwrap();
        let mut input = InputFile::new(&name, InputFile::NO_AUTODETECT).unwrap();
        let mut header = [0u8; 2];
        input.read_raw(&mut header).unwrap();
        assert_eq!(header, [0x1f, 0x8b]);

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_input_file_from_output_file() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-io-from-output-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        {
            let mut out = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
            out.write_raw(b"from-output").unwrap();
            let mut input = InputFile::from_output_file(&mut out, 0).unwrap();
            let mut buf = [0u8; 11];
            input.read_raw(&mut buf).unwrap();
            assert_eq!(&buf, b"from-output");
            assert!(input.temp_file);
            input.close().unwrap();
            out.close().unwrap();
        }
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_input_file_from_unlinked_temp_file() {
        let mut tmp = TempFile::new(true).unwrap();
        assert!(tmp.unlinked());
        tmp.write_raw(b"hidden-temp").unwrap();
        let mut input = InputFile::from_temp_file(&mut tmp, 0, Compressor::None).unwrap();
        assert!(input.temp_file);
        assert!(input.unlinked);
        #[cfg(unix)]
        assert!(!std::path::Path::new(&input.file_name).exists());
        let mut buf = [0u8; 11];
        input.read_raw(&mut buf).unwrap();
        assert_eq!(&buf, b"hidden-temp");
        input.close().unwrap();
        tmp.close().unwrap();
    }

    #[test]
    fn test_text_input_file_getline_putback() {
        let stream = VecStream::from_vec(b"a\r\nb\n".to_vec());
        let de = Deserializer::new(stream);
        let mut text = TextInputFile::new(de, b'\n');
        text.getline().unwrap();
        assert_eq!(text.line, "a");
        assert_eq!(text.line_count, 1);
        text.putback_line();
        text.getline().unwrap();
        assert_eq!(text.line, "a");
        text.getline().unwrap();
        assert_eq!(text.line, "b");
        assert!(!text.eof());
        text.getline().unwrap();
        assert_eq!(text.line, "");
        assert!(text.eof());
        text.rewind().unwrap();
        assert_eq!(text.line_count, 0);
    }

    #[test]
    fn test_text_input_file_constructors() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-text-input-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        std::fs::write(&path, b"one\r\ntwo\n").unwrap();

        let mut text = TextInputFile::from_file_name(&name, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "one");
        text.getline().unwrap();
        assert_eq!(text.line, "two");

        let mut out = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
        out.write_raw(b"out\n").unwrap();
        let mut text = TextInputFile::from_output_file(&mut out, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "out");
        out.close().unwrap();

        let data = TempFileData::init(false).unwrap();
        let mut tmp = TempFile::from_temp_file_data(&data).unwrap();
        tmp.write_raw(b"tmp\n").unwrap();
        let mut text = TextInputFile::from_temp_file(&mut tmp, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "tmp");
        tmp.close().unwrap();
        std::fs::remove_file(&data.name).unwrap();

        let mut unlinked = TempFile::new(true).unwrap();
        unlinked.write_raw(b"unlinked\n").unwrap();
        let mut text = TextInputFile::from_temp_file(&mut unlinked, b'\n').unwrap();
        text.getline().unwrap();
        assert_eq!(text.line, "unlinked");
        unlinked.close().unwrap();

        std::fs::remove_file(path).unwrap();
    }
}
