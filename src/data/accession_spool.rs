//! Source SequenceFile's one-String-column temporary accession TSV owner.
//! Rows stay on disk during search. Source materializes the actual typed table
//! only at write_accession_list, after creating its output file.
use std::io;
use std::path::Path;
use std::sync::Mutex;

use crate::util::io::{Compressor, InputFile, TempFile};
use crate::util::tsv::{Table, Type};

pub(crate) struct AccessionSpool {
    file: Mutex<SpoolFile>,
    errors: crate::util::io::IoErrorOwner,
}
struct SpoolFile {
    writer: TempFile,
    input: InputFile,
}
impl std::fmt::Debug for AccessionSpool {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("AccessionSpool(disk TSV)")
    }
}
impl AccessionSpool {
    pub(crate) fn new(directory: &Path, no_unlink: bool, buffer_size: usize) -> io::Result<Self> {
        Self::new_owned(directory, no_unlink, buffer_size, Default::default())
    }
    pub(crate) fn new_owned(
        directory: &Path,
        no_unlink: bool,
        buffer_size: usize,
        errors: crate::util::io::IoErrorOwner,
    ) -> io::Result<Self> {
        let data = TempFile::init_in(true, directory, no_unlink).map_err(io::Error::other)?;
        let mut file = TempFile::from_temp_file_data_with_buffer_size(&data, buffer_size)
            .map_err(io::Error::other)?;
        let input =
            InputFile::from_temp_file_with_buffer_size(&mut file, 0, Compressor::None, buffer_size)
                .map_err(io::Error::other)?;
        Ok(Self {
            file: Mutex::new(SpoolFile {
                writer: file,
                input,
            }),
            errors,
        })
    }
    pub(crate) fn io_errors(&self) -> crate::util::io::IoErrorOwner {
        self.errors.clone()
    }
    /// Source File::write_record(string): raw text followed by LF, not a
    /// private binary or escaped format. Parsing is deferred until final read.
    pub(crate) fn push(&self, accession: &str) -> io::Result<()> {
        self.push_bytes(accession.as_bytes())
    }
    pub(crate) fn push_bytes(&self, accession: &[u8]) -> io::Result<()> {
        self.errors.capture_io(self.push_inner(accession))
    }
    fn push_inner(&self, accession: &[u8]) -> io::Result<()> {
        let mut file = self
            .file
            .lock()
            .map_err(|_| io::Error::other("Accession spool lock poisoned."))?;
        file.writer.write_raw(accession).map_err(io::Error::other)?;
        file.writer.write_raw(b"\n").map_err(io::Error::other)?;
        Ok(())
    }
    pub(crate) fn read_table(&self, threads: usize) -> io::Result<Table> {
        self.errors.capture_io(self.read_table_inner(threads))
    }
    fn read_table_inner(&self, threads: usize) -> io::Result<Table> {
        let mut file = self
            .file
            .lock()
            .map_err(|_| io::Error::other("Accession spool lock poisoned."))?;
        // Source TSV File::rewind flushes the writer before rewinding the
        // already-existing TextInputFile alias of that same FILE.
        file.writer.rewind().map_err(io::Error::other)?;
        file.input.rewind().map_err(io::Error::other)?;
        crate::util::tsv::streaming::StreamingFile::read_input_io(
            &mut file.input,
            vec![Type::String],
            i64::MAX,
            threads,
        )
    }
}

/// Keep the real final owner alive until its operation's primary I/O error is
/// saved. This is a lexical destructor boundary, not a second mapping/cache.
pub(crate) fn source_scope<T>(
    spool: Option<std::sync::Arc<AccessionSpool>>,
    operation: impl FnOnce() -> io::Result<T>,
) -> io::Result<T> {
    let errors = spool
        .as_ref()
        .map_or_else(Default::default, |spool| spool.io_errors());
    let result = operation();
    errors.capture_io(result)
}

impl Drop for AccessionSpool {
    fn drop(&mut self) {
        let file = match self.file.get_mut() {
            Ok(file) => file,
            Err(poisoned) => poisoned.into_inner(),
        };
        let name = file.writer.file_path().to_path_buf();
        let unlinked = file.writer.unlinked();
        if let Err(error) = file.input.close() {
            if self.errors.source_destructors() {
                if let Err(error) = self.errors.emit_pending_perror() {
                    crate::util::io::abort_source_diagnostic_error(&error);
                }
                if let Err(conversion) =
                    crate::util::io::emit_source_perror(&io::Error::other(error))
                {
                    crate::util::io::abort_source_diagnostic_error(&conversion);
                }
                crate::util::io::abort_source_process();
            }
            // The reusable Rust extension remains nonfatal; Source's close
            // failed before its explicit remove, so do not delete this file.
            return;
        }
        // Source TSV File::close explicitly calls InputFile::close_and_delete;
        // --no-unlink only controls initial unlinking, not this owner cleanup.
        if !unlinked && std::fs::remove_file(&name).is_err() {
            crate::util::log_stream::flush_tied_stdout();
            let path = match crate::util::io::path_diagnostic_bytes(&name) {
                Ok(path) => path,
                Err(error) => {
                    if self.errors.source_destructors() {
                        crate::util::io::abort_source_diagnostic_error(&error);
                    }
                    let _ = self.errors.capture::<()>(Err(io::Error::other(error)));
                    return;
                }
            };
            {
                let mut message = b"Warning: Failed to delete temporary file ".to_vec();
                message.extend_from_slice(&path);
                message.push(b'\n');
                let _ =
                    std::io::Write::write_all(&mut crate::util::io::console::stderr(), &message);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn retained_input_reads_rewound_writer_and_drop_discards_pending_serializer() {
        let spool = AccessionSpool::new(&std::env::temp_dir(), true, 64).unwrap();
        let name = spool.file.lock().unwrap().writer.file_name();
        spool.push("one").unwrap();
        assert_eq!(std::fs::metadata(&name).unwrap().len(), 0);
        let table = spool.read_table(1).unwrap();
        assert_eq!(table.size(), 1);
        assert_eq!(table.record(0).string_bytes(0), b"one");
        // The same retained alias is rewound, not replaced/closed by read_table.
        assert_eq!(spool.read_table(1).unwrap().size(), 1);
        spool.push("pending").unwrap();
        assert_eq!(std::fs::metadata(&name).unwrap().len(), 4);
        drop(spool);
        assert!(!Path::new(&name).exists());
    }
}
