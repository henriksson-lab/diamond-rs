use std::io;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::basic::value::SequenceType;
use crate::data::db_builder;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MakedbOptions {
    pub taxonmap: Option<PathBuf>,
    pub taxonnodes: Option<PathBuf>,
    pub taxonnames: Option<PathBuf>,
    pub no_parse_seqids: bool,
}

pub(crate) struct RuntimeOptions {
    pub temporary_directory: std::path::PathBuf,
    pub file_buffer_size: Option<usize>,
    pub sequence_type: SequenceType,
    pub no_auto_append: bool,
    pub mask: bool,
    pub matrix: Vec<u8>,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub cbs: u32,
    pub tuning: super::blastp::SearchTuning,
}
impl Default for RuntimeOptions {
    fn default() -> Self {
        Self {
            temporary_directory: std::env::temp_dir(),
            file_buffer_size: None,
            sequence_type: SequenceType::AminoAcid,
            no_auto_append: false,
            mask: true,
            matrix: b"blosum62".to_vec(),
            gap_open: -1,
            gap_extend: -1,
            cbs: 1,
            tuning: Default::default(),
        }
    }
}

/// Run the makedb command — build a DIAMOND database from FASTA file(s).
pub fn run(input_files: &[String], database: &str, threads: i32) -> io::Result<()> {
    run_with_options(input_files, database, threads, &MakedbOptions::default())
}

pub fn run_with_options(
    input_files: &[String],
    database: &str,
    threads: i32,
    options: &MakedbOptions,
) -> io::Result<()> {
    run_with_buffer(input_files, database, threads, options, None)
}

pub(crate) fn run_with_buffer(
    input_files: &[String],
    database: &str,
    threads: i32,
    options: &MakedbOptions,
    file_buffer_size: Option<usize>,
) -> io::Result<()> {
    run_with_runtime(
        input_files,
        database,
        threads,
        options,
        RuntimeOptions {
            file_buffer_size,
            ..Default::default()
        },
    )
}

pub(crate) fn run_with_runtime(
    input_files: &[String],
    database: &str,
    threads: i32,
    options: &MakedbOptions,
    runtime: RuntimeOptions,
) -> io::Result<()> {
    run_with_runtime_paths(
        &input_files.iter().map(PathBuf::from).collect::<Vec<_>>(),
        Path::new(database),
        threads,
        options,
        runtime,
    )
}

pub(crate) fn run_with_runtime_paths(
    input_files: &[PathBuf],
    database: &Path,
    threads: i32,
    options: &MakedbOptions,
    mut runtime: RuntimeOptions,
) -> io::Result<()> {
    // DatabaseFile::make_db resets this global before creating its FASTA and
    // output streams, regardless of the user-supplied command capacity.
    runtime.file_buffer_size = Some(4 * 1024 * 1024);

    // C++ rejects multiple values in the effective `--in` occurrence with
    // "Too many arguments provided for option --in."
    // (`diamond/src/data/dmnd/dmnd.cpp:261-262`). The CLI parser separately
    // mirrors upstream assignment semantics where a repeated occurrence
    // replaces the earlier one.
    if input_files.len() > 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Too many arguments provided for option --in.",
        ));
    }
    if threads > 0 && rayon::current_thread_index().is_none() {
        return rayon::ThreadPoolBuilder::new()
            .num_threads(threads as usize)
            .build()
            .map_err(io::Error::other)?
            .install(|| run_with_runtime_paths(input_files, database, threads, options, runtime));
    }
    if runtime.cbs >= if cfg!(feature = "extra") { 7 } else { 6 } {
        return Err(io::Error::other(
            "Invalid value for --comp-based-stats. Permitted values: 0, 1, 2, 3, 4, 5.",
        ));
    }
    if !runtime.tuning.matrix_file.as_os_str().is_empty() && runtime.cbs > 1 {
        return Err(io::Error::other("This value for --comp-based-stats is not supported when using a custom scoring matrix."));
    }
    let matrix = super::blastp::native_score_matrix_with_tuning_bytes(
        &runtime.matrix,
        runtime.gap_open,
        runtime.gap_extend,
        0,
        &runtime.tuning,
        0,
    )?;
    let masker = runtime.tuning.runtime.tantan().unwrap_or_else(|| {
        std::sync::Arc::new(crate::masking::tantan::TantanMasker::from_score_matrix(
            &matrix,
            runtime.tuning.tantan_min_mask_prob as f32,
        ))
    });
    let input_refs: Vec<&Path> = if input_files.is_empty() {
        vec![Path::new("")]
    } else {
        input_files.iter().map(PathBuf::as_path).collect()
    };

    // Source auto_append_extension compares the literal ASCII suffix while
    // retaining all original pathname units, rather than replacing extensions.
    let mut db_path = database.to_path_buf();
    if !runtime.no_auto_append
        && !db_path.as_os_str().is_empty()
        && !db_path.as_os_str().as_encoded_bytes().ends_with(b".dmnd")
    {
        db_path.as_mut_os_string().push(".dmnd");
    }

    let mut database_message = b"Database file: ".to_vec();
    database_message.extend_from_slice(
        &crate::util::io::source_os_bytes(db_path.as_os_str()).map_err(io::Error::other)?,
    );
    database_message.push(b'\n');
    crate::util::log_stream::message_stream()
        .lock()
        .unwrap()
        .write_bytes(&database_message)?;
    let mut inputs_message = b"Input files: ".to_vec();
    for (index, path) in input_files.iter().enumerate() {
        if index != 0 {
            inputs_message.extend_from_slice(b", ");
        }
        inputs_message.extend_from_slice(
            &crate::util::io::source_os_bytes(path.as_os_str()).map_err(io::Error::other)?,
        );
    }
    inputs_message.push(b'\n');
    crate::util::log_stream::message_stream()
        .lock()
        .unwrap()
        .write_bytes(&inputs_message)?;

    let start = Instant::now();
    db_builder::build_db_with_runtime_paths(
        &input_refs,
        &db_path,
        runtime.sequence_type,
        &db_builder::BuildOptions {
            taxonmap: options.taxonmap.clone(),
            taxonnodes: options.taxonnodes.clone(),
            taxonnames: options.taxonnames.clone(),
            no_parse_seqids: options.no_parse_seqids,
            file_buffer_size: runtime.file_buffer_size,
        },
        db_builder::BuildRuntime {
            temporary_directory: runtime.temporary_directory,
            mask: runtime.mask,
            masker: Some(masker),
            ignore_warnings: runtime.tuning.ignore_warnings,
            no_unlink: runtime.tuning.no_unlink,
            source_stats_start: Some(start),
        },
    )?;

    Ok(())
}
