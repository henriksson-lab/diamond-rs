//! Concrete `SequenceFile::db_info` command for DMND and NCBI BLAST databases.
use std::ffi::OsString;
use std::io::{self, Write};
use std::path::Path;

use crate::data::blastdb::blastdb::BlastDB;
use crate::data::dmnd::{DatabaseOpenOptions, SourceDatabaseFile};
use crate::data::sequence_file::{
    detect_type_source_with_buffer, SequenceFileConfig, SequenceFileFlags, SequenceFileType,
};

pub fn run(database_path: &str) -> io::Result<()> {
    run_path(Path::new(database_path))
}
pub fn run_path(database_path: &Path) -> io::Result<()> {
    run_path_with_buffer_size(database_path, 67_108_864)
}
pub(crate) fn run_os(args: &[OsString]) -> io::Result<()> {
    run_with_startup_os(args, None)
}
pub(crate) fn run_with_startup_os(
    args: &[OsString],
    startup: Option<std::sync::Arc<super::source_config::Startup>>,
) -> io::Result<()> {
    use crate::util::command_line_parser::ReadOption;
    let startup = match startup {
        Some(startup) => startup,
        None => super::source_config::prepare_os(args)?,
    };
    let values = super::argument_registry::values_os(args)?;
    let mut size = 67_108_864i64;
    if let Some(value) = values.get("file-buffer-size") {
        i64::read_option_os(&mut size, value)?;
    }
    run_path_with_buffer_size(&startup.config.database, size as usize)
}

pub fn run_with_writer(database_path: &str, out: &mut impl Write) -> io::Result<()> {
    run_with_writer_and_buffer_size(database_path, out, 67_108_864)
}

pub fn run_with_buffer_size(database_path: &str, file_buffer_size: usize) -> io::Result<()> {
    run_with_writer_and_buffer_size(
        database_path,
        &mut crate::util::io::console::stdout_endl(),
        file_buffer_size,
    )
}

pub fn run_with_writer_and_buffer_size(
    database_path: &str,
    out: &mut impl Write,
    file_buffer_size: usize,
) -> io::Result<()> {
    run_path_with_writer_and_buffer_size(Path::new(database_path), out, file_buffer_size)
}
pub fn run_path_with_buffer_size(database_path: &Path, file_buffer_size: usize) -> io::Result<()> {
    run_path_with_writer_and_buffer_size(
        database_path,
        &mut crate::util::io::console::stdout_endl(),
        file_buffer_size,
    )
}
pub fn run_path_with_writer_and_buffer_size(
    database_path: &Path,
    out: &mut impl Write,
    file_buffer_size: usize,
) -> io::Result<()> {
    if database_path.as_os_str().is_empty() {
        return Err(io::Error::other(
            "Missing option for database file: --db/-d.",
        ));
    }
    let flags = SequenceFileFlags::NO_FASTA | SequenceFileFlags::NO_COMPATIBILITY_CHECK;
    let (kind, path) = detect_type_source_with_buffer(
        &[database_path.to_owned()],
        flags,
        SequenceFileConfig::default(),
        file_buffer_size,
    )?;
    let (version, build, sequences, letters) = match kind {
        SequenceFileType::Dmnd => {
            // Inspect unsupported versions while retaining the real database
            // constructor's header/record reads and malformed-file errors.
            let database = SourceDatabaseFile::open(
                &path,
                DatabaseOpenOptions {
                    no_compatibility_check: true,
                    file_buffer_size: Some(file_buffer_size),
                    ..Default::default()
                },
            )?;
            (
                database.ref_header.db_version as i64,
                Some(database.ref_header.build),
                database.ref_header.sequences,
                database.ref_header.letters,
            )
        }
        SequenceFileType::Blast => {
            let database = BlastDB::new_path(&path, flags)?;
            (
                database.db_version() as i64,
                None,
                database.sequence_count(),
                database.letters(),
            )
        }
        _ => {
            return Err(io::Error::other(
                "Sequence file does not have a supported format.",
            ))
        }
    };
    writeln!(out, "{:>23}  {kind}", "Database type")?;
    writeln!(out, "{:>23}  {version}", "Database format version")?;
    if let Some(build) = build {
        writeln!(out, "{:>23}  {build}", "Diamond build")?;
    }
    writeln!(out, "{:>23}  {sequences}", "Sequences")?;
    writeln!(out, "{:>23}  {letters}", "Letters")?;
    Ok(())
}
