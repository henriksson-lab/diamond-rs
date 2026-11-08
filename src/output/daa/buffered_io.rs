//! Source InputFile/OutputFile adapters for native DAA consumers.
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::path::Path;

pub(crate) const DEFAULT_BUFFER_SIZE: usize = 67_108_864;

pub(crate) enum DaaInput {
    Source(crate::util::io::InputFile),
}
impl DaaInput {
    pub(crate) fn open(path: &Path, capacity: usize) -> io::Result<Self> {
        crate::util::io::InputFile::with_source_profile_path(path, 0, capacity)
            .map(Self::Source)
            .map_err(io::Error::other)
    }

    pub(crate) fn seek_absolute(&mut self, position: u64) -> io::Result<()> {
        match self {
            Self::Source(input) => {
                input
                    .seek(i64::try_from(position).map_err(io::Error::other)?)
                    .map_err(io::Error::other)?;
            }
        }
        Ok(())
    }
}
impl Read for DaaInput {
    fn read_exact(&mut self, mut buffer: &mut [u8]) -> io::Result<()> {
        while !buffer.is_empty() {
            match self.read(buffer) {
                Ok(0) => {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        crate::util::io::IoError::EndOfStream,
                    ))
                }
                Ok(count) => buffer = &mut buffer[count..],
                Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(error),
            }
        }
        Ok(())
    }
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Source(input) => input.read_raw(buffer).map_err(io::Error::other),
        }
    }
}
// Source DAAFile/InputFile destructors do not execute an explicit close
// operation (merge even retains raw DAAFile pointers). Rust-owned primitives
// still release OS resources best-effort, without introducing close diagnostics.

pub(crate) enum DaaOutput {
    Source(crate::util::io::OutputFile),
    // Source ViewWriter explicitly closes in its noexcept destructor. Plain
    // merge OutputFile destruction instead abandons Serializer staging.
    SourceView(crate::util::io::OutputFile),
}
impl DaaOutput {
    pub(crate) fn create(path: &Path, capacity: usize) -> io::Result<Self> {
        Self::create_with_compressor(path, capacity, crate::util::io::Compressor::None)
    }
    pub(crate) fn create_with_compressor(
        path: &Path,
        capacity: usize,
        compressor: crate::util::io::Compressor,
    ) -> io::Result<Self> {
        Self::create_with_policy(path, capacity, compressor, false)
    }
    fn create_with_policy(
        path: &Path,
        capacity: usize,
        compressor: crate::util::io::Compressor,
        view: bool,
    ) -> io::Result<Self> {
        // Source OutputFile owns all pathname encodings and the stdout policy.
        crate::util::io::OutputFile::with_source_profile_path(path, compressor, "wb", capacity)
            .map(|output| {
                if view {
                    Self::SourceView(output)
                } else {
                    Self::Source(output)
                }
            })
            .map_err(io::Error::other)
    }

    pub(crate) fn create_view_with_compressor(
        path: &Path,
        capacity: usize,
        compressor: crate::util::io::Compressor,
    ) -> io::Result<Self> {
        Self::create_with_policy(path, capacity, compressor, true)
    }
    pub(crate) fn close(&mut self) -> io::Result<()> {
        match self {
            Self::Source(output) => output.close().map_err(io::Error::other),
            Self::SourceView(output) => {
                if let Err(error) = output.close() {
                    if let Err(error) = crate::util::io::emit_source_perror(&error) {
                        crate::util::io::abort_source_diagnostic_error(&error);
                    }
                    crate::util::io::abort_source_process();
                }
                Ok(())
            }
        }
    }
}
impl Write for DaaOutput {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        match self {
            Self::Source(output) | Self::SourceView(output) => {
                output.write_raw(bytes).map_err(io::Error::other)?;
                Ok(bytes.len())
            }
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Source(output) | Self::SourceView(output) => {
                output.flush().map_err(io::Error::other)
            }
        }
    }
}
impl Seek for DaaOutput {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        match self {
            Self::Source(output) | Self::SourceView(output) => {
                let (offset, origin) = match position {
                    SeekFrom::Start(offset) => (
                        i64::try_from(offset).map_err(io::Error::other)?,
                        SeekFrom::Start(0),
                    ),
                    SeekFrom::Current(offset) => (offset, SeekFrom::Current(0)),
                    SeekFrom::End(offset) => (offset, SeekFrom::End(0)),
                };
                output.seek(offset, origin).map_err(io::Error::other)?;
                output.tell().map(|p| p as u64).map_err(io::Error::other)
            }
        }
    }
}
impl Drop for DaaOutput {
    fn drop(&mut self) {
        if matches!(self, Self::SourceView(_)) {
            // ViewWriter's exception-disabled destructor closes on success
            // and unwind. A close failure terminates, not a generic catch.
            let _ = self.close();
        }
        // Source merge's plain OutputFile dtor discards pending outer data.
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn path(label: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-daa-{label}-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ))
    }
    #[cfg(unix)]
    #[test]
    fn raw_path_uses_source_file_owner_and_preserves_open_error() {
        let root = path("raw-path");
        std::fs::create_dir(&root).unwrap();
        let (name, _) =
            crate::util::test_os_path_or_unicode(&root, b"archive_\xff.daa", "archive-é.daa");
        {
            let mut output = DaaOutput::create(&name, 8192).unwrap();
            output.write_all(b"raw").unwrap();
            output.close().unwrap();
        }
        let mut input = DaaInput::open(&name, 8192).unwrap();
        let mut bytes = [0; 3];
        input.read_exact(&mut bytes).unwrap();
        assert_eq!(&bytes, b"raw");
        drop(input);
        std::fs::remove_file(&name).unwrap();
        let error = match DaaInput::open(&name, 8192) {
            Err(error) => error,
            Ok(_) => panic!("missing archive opened"),
        };
        let details = crate::util::io::file_open_error::find(&error).unwrap();
        assert_eq!(details.path, name);
        std::fs::remove_dir(root).unwrap();
    }

    #[test]
    fn merge_output_abandons_pending_header_but_view_writer_closes_it() {
        let merge = path("merge-discard");
        {
            let mut output = DaaOutput::create(&merge, 8192).unwrap();
            crate::data::daa::init_daa(&mut output).unwrap();
        }
        assert_eq!(std::fs::read(&merge).unwrap(), b"");
        let view = path("view-close");
        {
            let mut output = DaaOutput::create_view_with_compressor(
                &view,
                8192,
                crate::util::io::Compressor::None,
            )
            .unwrap();
            crate::data::daa::init_daa(&mut output).unwrap();
        }
        assert_eq!(
            std::fs::metadata(&view).unwrap().len(),
            (crate::data::daa::DAA_HEADER1_SIZE + crate::data::daa::DAA_HEADER2_SIZE) as u64
        );
        std::fs::remove_file(merge).unwrap();
        std::fs::remove_file(view).unwrap();
    }
    #[test]
    fn physical_daa_eof_preserves_source_deserializer_diagnostic() {
        let path = path("eof");
        std::fs::write(&path, []).unwrap();
        let mut input = DaaInput::open(&path, 8192).unwrap();
        let error = input.read_exact(&mut [0; 16]).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
        assert_eq!(error.to_string(), "Unexpected end of input.");
        drop(input);
        std::fs::remove_file(path).unwrap();
    }
}
