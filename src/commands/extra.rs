//! Optional Source EXTRA tools connected to real native dependencies.
//!
//! Source comments out its record reads in mask/fastq2fasta/reverse, leaving
//! nonterminating loops. These commands implement the intended record reads
//! with defined EOF termination. This is not literal nontermination parity.
//! Merely enabling the EXTRA registry does not implement its other commands.

#[cfg(feature = "extra")]
use std::io;

#[cfg(feature = "extra")]
pub(crate) struct PreparedExtra {
    args: Vec<std::ffi::OsString>,
    original_args: Vec<std::ffi::OsString>,
    command_name: &'static str,
    startup: std::sync::Arc<super::source_config::Startup>,
}

#[cfg(feature = "extra")]
impl PreparedExtra {
    pub(crate) fn command(&self) -> crate::run::main::RunCommand {
        use crate::run::main::RunCommand::*;
        match self.command_name {
            "random-seqs" => RandomSeqs,
            "mask" => Mask,
            "fastq2fasta" => Fastq2Fasta,
            "info" => Info,
            "smith-waterman" => SmithWaterman,
            "split" => Split,
            "reverse" => ReverseSeqs,
            "roc-id" => RocId,
            "find-shapes" => FindShapes,
            "hashseqs" => HashSeqs,
            "listseeds" => ListSeeds,
            "blastn" => BlastN,
            "wc" => WordCount,
            "cut" => Cut,
            "profile-recluster" => ProfileRecluster,
            _ => Unknown,
        }
    }
    pub(crate) fn execute(self) -> io::Result<i32> {
        if super::extra_tools::source_fallthrough_os(&self.args)? {
            return Ok(1);
        }
        if enabled::run(&self.args, &self.original_args, &self.startup)? {
            Ok(0)
        } else {
            Err(io::Error::other(
                "Native execution of this command or option is not yet connected.",
            ))
        }
    }
}

#[cfg(feature = "extra")]
pub(crate) fn prepare(
    args: &[String],
    original_args: &[String],
    startup: std::sync::Arc<super::source_config::Startup>,
) -> io::Result<Option<PreparedExtra>> {
    let args = args
        .iter()
        .map(std::ffi::OsString::from)
        .collect::<Vec<_>>();
    let original_args = original_args
        .iter()
        .map(std::ffi::OsString::from)
        .collect::<Vec<_>>();
    prepare_os(&args, &original_args, startup)
}
#[cfg(feature = "extra")]
pub(crate) fn prepare_os(
    args: &[std::ffi::OsString],
    original_args: &[std::ffi::OsString],
    startup: std::sync::Arc<super::source_config::Startup>,
) -> io::Result<Option<PreparedExtra>> {
    let command_name = super::argument_registry::command_name_os(args)?.unwrap_or("");
    if !super::argument_registry::is_extra_command(command_name) {
        return Ok(None);
    }
    Ok(Some(PreparedExtra {
        args: args.to_vec(),
        original_args: original_args.to_vec(),
        command_name,
        startup,
    }))
}

#[cfg(feature = "extra")]
pub(crate) fn run_os(args: &[std::ffi::OsString]) -> io::Result<i32> {
    run_with_startup_os(args, None)
}
#[cfg(feature = "extra")]
pub(crate) fn run_with_startup_os(
    args: &[std::ffi::OsString],
    startup: Option<std::sync::Arc<super::source_config::Startup>>,
) -> io::Result<i32> {
    let startup = match startup {
        Some(startup) => startup,
        None => super::source_config::prepare_os(args)?,
    };
    prepare_os(args, args, startup)?
        .ok_or_else(|| io::Error::other("Not an EXTRA command"))?
        .execute()
}

#[cfg(feature = "extra")]
mod enabled {
    use std::collections::BTreeMap;
    use std::ffi::OsString;
    use std::io::{self, Write};
    use std::path::{Path, PathBuf};

    use crate::basic::packed_transcript::EditOperation;
    use crate::basic::value::{Letter, SequenceType, NUCLEOTIDE_ALPHABET};
    use crate::data::dmnd::DatabaseFile;
    use crate::data::fasta::{self, ByteFastaRecord, FastaRecord};
    use crate::run::tools::{
        self, PairwiseColumn, PairwiseOperation, PairwiseSettings, ToolSequenceType, ToolsBackend,
        ToolsConfig,
    };
    use crate::util::command_line_parser::ReadOption;
    use crate::util::compat_rng;
    use crate::util::io::{Compressor, OutputFile};

    /// Concrete file-backed backend for the translated generic tools.
    pub(crate) struct NativeToolsBackend<'a> {
        startup: &'a super::super::source_config::Startup,
        file_buffer_size: usize,
        // The translated trait returns String. Retain a failed concrete I/O
        // operation separately and restore its type at this command boundary.
        errors: crate::util::io::IoErrorOwner,
    }

    impl ToolsBackend for NativeToolsBackend<'_> {
        fn mask_sequence(
            &mut self,
            sequence: &mut [Letter],
            _probability: f32,
        ) -> Result<(), String> {
            self.startup
                .masking
                .as_ref()
                .ok_or("Masking was not initialized")?
                .apply_masking(
                    sequence,
                    crate::masking::MaskingAlgo::Tantan as u32,
                    0,
                    None,
                )?;
            Ok(())
        }
        fn get_seq(&mut self, config: &ToolsConfig) -> Result<(), String> {
            let database = config
                .database
                .as_deref()
                .ok_or_else(|| "Missing required parameter: database".to_owned())?;
            self.errors
                .capture(super::super::getseq::run_path_with_output_numbers(
                    database,
                    &config.seq_no,
                    config.output_file.as_deref(),
                ))
        }

        fn database_sequences(&mut self, database: &str) -> Result<Vec<FastaRecord>, String> {
            self.database_sequences_bytes(database)?
                .into_iter()
                .map(|record| {
                    record
                        .try_into()
                        .map_err(|error: std::string::FromUtf8Error| error.to_string())
                })
                .collect()
        }
        fn database_sequences_bytes(
            &mut self,
            database: &str,
        ) -> Result<Vec<ByteFastaRecord>, String> {
            self.database_sequences_path_bytes(Path::new(database))
        }
        fn database_sequences_path_bytes(
            &mut self,
            database: &Path,
        ) -> Result<Vec<ByteFastaRecord>, String> {
            // Source random-seqs specifically opens DatabaseFile, rather than
            // auto_create; FASTA/BLAST are not silently substituted here.
            super::super::blastp::preflight_file_buffer_size(self.file_buffer_size)
                .map_err(|error| error.to_string())?;
            let mut database = self.errors.capture(DatabaseFile::open_with_options_io(
                database,
                crate::data::dmnd::DatabaseOpenOptions {
                    file_buffer_size: Some(self.file_buffer_size),
                    ..Default::default()
                },
            ))?;
            let count = usize::try_from(database.sequence_count())
                .map_err(|_| "Database sequence count exceeds address space.".to_owned())?;
            let mut records = Vec::new();
            records
                .try_reserve(count)
                .map_err(|error| error.to_string())?;
            for oid in 0..count {
                let (sequence, id) = self.errors.capture(database.record_bytes_io(oid as u64))?;
                records.push(ByteFastaRecord {
                    id,
                    sequence,
                    quality: None,
                });
            }
            Ok(records)
        }

        fn read_sequences(
            &mut self,
            path: &str,
            kind: ToolSequenceType,
        ) -> Result<Vec<FastaRecord>, String> {
            self.read_sequences_bytes(path, kind)?
                .into_iter()
                .map(|record| {
                    record
                        .try_into()
                        .map_err(|error: std::string::FromUtf8Error| error.to_string())
                })
                .collect()
        }
        fn read_sequences_bytes(
            &mut self,
            path: &str,
            kind: ToolSequenceType,
        ) -> Result<Vec<ByteFastaRecord>, String> {
            self.read_sequences_path_bytes(Path::new(path), kind)
        }
        fn read_sequences_path_bytes(
            &mut self,
            path: &Path,
            kind: ToolSequenceType,
        ) -> Result<Vec<ByteFastaRecord>, String> {
            let kind = match kind {
                ToolSequenceType::AminoAcid => SequenceType::AminoAcid,
                ToolSequenceType::Nucleotide => SequenceType::Nucleotide,
            };
            super::super::blastp::preflight_file_buffer_size(self.file_buffer_size)
                .map_err(|error| error.to_string())?;
            let mut records = Vec::new();
            let emit = |record| {
                records.push(record);
                Ok(())
            };
            let result = fasta::for_each_fasta_file_with_buffer_size_bytes(
                path,
                kind,
                self.file_buffer_size,
                emit,
            );
            self.errors.capture(result)?;
            Ok(records)
        }

        fn random_index(&mut self, count: usize) -> usize {
            // Source rand()*RAND_MAX+rand() overflows signed int on GNU.
            // Specify the observed wrapping arithmetic safely, followed by
            // its size_t conversion, using the translated GNU RNG stream.
            let first = compat_rng::next_i31();
            let second = compat_rng::next_i31();
            let index = first
                .wrapping_mul(compat_rng::RAND_MAX)
                .wrapping_add(second) as usize
                % count;
            // The source set stores unsigned, even when r.size() is size_t.
            index as u32 as usize
        }

        fn pairwise_alignment(
            &mut self,
            reference: &[Letter],
            query: &[Letter],
            settings: PairwiseSettings,
        ) -> Result<Vec<PairwiseColumn>, String> {
            // This dependency is real DP, not a fabricated transcript. The
            // checked-in Source has no DNA entry in StandardMatrix::get:
            // pairwise's hardcoded "DNA" fails before reading records.
            // Requesting that matrix preserves its explicit native error.
            let matrix = crate::stats::score_matrix::ScoreMatrix::new(
                settings.matrix_name,
                settings.gap_open,
                settings.gap_extend,
                settings.reward,
                settings.penalty,
                0,
            )?;
            let result = crate::dp::smith_waterman::smith_waterman(query, reference, &matrix);
            let mut query_pos = result.query_begin;
            let mut subject_pos = result.subject_begin;
            let mut columns = Vec::new();
            for (operation, count) in result.operations {
                for _ in 0..count {
                    match operation {
                        EditOperation::Match | EditOperation::Substitution => {
                            let index = (query[query_pos as usize] as u8 & 127) as usize;
                            columns.push(PairwiseColumn {
                                operation: PairwiseOperation::Substitution,
                                subject_pos,
                                query_pos,
                                query_char: NUCLEOTIDE_ALPHABET[index] as char,
                            });
                            subject_pos += 1;
                            query_pos += 1;
                        }
                        EditOperation::Deletion => {
                            columns.push(PairwiseColumn {
                                operation: PairwiseOperation::Deletion,
                                subject_pos,
                                query_pos: -1,
                                query_char: '-',
                            });
                            subject_pos += 1;
                        }
                        EditOperation::Insertion => query_pos += 1,
                        _ => {
                            return Err("Unsupported frameshift in pairwise tool alignment.".into())
                        }
                    }
                }
            }
            Ok(columns)
        }

        fn prepare_pairwise(&mut self, settings: PairwiseSettings) -> Result<(), String> {
            // tools.cpp initializes this before opening input or inspecting
            // the worker count. StandardMatrix::get in the checked-in Source
            // has no DNA entry, so smith-waterman is a broken runtime mode
            // with a defined error, even for empty/missing input or -p0.
            crate::stats::score_matrix::ScoreMatrix::new(
                settings.matrix_name,
                settings.gap_open,
                settings.gap_extend,
                settings.reward,
                settings.penalty,
                0,
            )
            .map(|_| ())
        }
    }

    struct ToolOutput(OutputFile, crate::util::io::IoErrorOwner);
    impl Write for ToolOutput {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.1
                .capture_io(self.0.write_raw(bytes).map_err(io::Error::other))?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            self.1.capture_io(self.0.flush().map_err(io::Error::other))
        }
    }

    fn scalar(values: &BTreeMap<&str, Vec<OsString>>, name: &str) -> Option<PathBuf> {
        values
            .get(name)
            .and_then(|items| items.first())
            .cloned()
            .map(Into::into)
    }

    fn config(args: &[OsString]) -> io::Result<ToolsConfig> {
        let values = super::super::argument_registry::extra_values_os(args)?;
        let mut config = ToolsConfig {
            database: scalar(&values, "db"),
            // Source single_query_file defaults to stdin and rejects >1 file.
            query_file: Some(PathBuf::new()),
            output_file: scalar(&values, "out"),
            seq_no: values
                .get("seq")
                .into_iter()
                .flatten()
                .map(|s| {
                    crate::util::io::source_os_bytes(s)
                        .map(|s| s.into_owned())
                        .map_err(io::Error::other)
                })
                .collect::<io::Result<_>>()?,
            reverse: values.contains_key("reverse"),
            hardmasked: values.contains_key("hard-masked"),
            ..ToolsConfig::default()
        };
        if let Some(paths) = values.get("query") {
            if paths.len() > 1 {
                return Err(io::Error::other("--query/-q has more than one argument."));
            }
            config.query_file = paths.first().cloned().map(Into::into);
        }
        if let Some(items) = values.get("threads") {
            i32::read_option_os(&mut config.threads, items)?;
        }
        if let Some(items) = values.get("block-size") {
            f64::read_option_os(&mut config.chunk_size, items)?;
        }
        if let Some(items) = values.get("tantan-minMaskProb") {
            let mut probability = f64::from(config.tantan_min_mask_prob);
            f64::read_option_os(&mut probability, items)?;
            config.tantan_min_mask_prob = probability as f32;
        }
        Ok(config)
    }

    pub(super) fn run(
        args: &[OsString],
        original_args: &[OsString],
        startup: &super::super::source_config::Startup,
    ) -> io::Result<bool> {
        let command = super::super::argument_registry::command_name_os(args)?;
        #[cfg(not(feature = "dna"))]
        let _ = original_args;
        #[cfg(feature = "dna")]
        if command == Some("blastn") {
            super::super::blastn::run_with_startup_os(original_args, Some(startup))?;
            return Ok(true);
        }
        #[cfg(not(feature = "dna"))]
        if command == Some("blastn") {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "Nucleotide search requires the `dna` Cargo feature.",
            ));
        }
        if super::super::extra_tools::run_os(args, startup)? {
            return Ok(true);
        }
        if !matches!(
            command,
            Some("random-seqs" | "mask" | "fastq2fasta" | "info" | "reverse" | "smith-waterman")
        ) {
            return Ok(false);
        }
        let mut config = config(args)?;
        config.output_file = Some(startup.config.output_file.clone());
        let values = super::super::argument_registry::extra_values_os(args)?;
        let mut file_buffer_size = 67_108_864usize;
        if let Some(items) = values.get("file-buffer-size") {
            usize::read_option_os(&mut file_buffer_size, items)?;
        }
        let mut backend = NativeToolsBackend {
            startup,
            file_buffer_size,
            errors: crate::util::io::IoErrorOwner::default(),
        };
        // Mask/FASTQ conversion/status each use Source endl; reverse and
        // pairwise insert complete buffers, while info ends with one endl.
        let mut out = if matches!(command, Some("mask" | "fastq2fasta" | "random-seqs")) {
            crate::util::io::console::stdout_endl()
        } else {
            crate::util::io::console::stdout()
        };
        let result = match command.unwrap() {
            "info" => tools::info(&mut out),
            "mask" => tools::run_masker(
                &mut backend,
                &config,
                &mut out,
                &mut crate::util::io::console::stderr(),
            )
            .map(|_| ()),
            "fastq2fasta" => tools::fastq2fasta(&mut backend, &config, &mut out),
            "reverse" => tools::reverse(&mut backend, &config, &mut out),
            "smith-waterman" => tools::pairwise(&mut backend, &config, &mut out),
            "random-seqs" => {
                // Source opens its output after reading/selecting the DB.
                // Delay output creation until the first output byte so a
                // database/read/count error does not truncate the output.
                let mut output = LazyOutput {
                    path: config.output_file.as_deref().unwrap_or(Path::new("")),
                    errors: backend.errors.clone(),
                    file: None,
                    buffer_size: file_buffer_size,
                };
                let result = tools::random_seqs(&mut backend, &config, &mut output, &mut out);
                if let Err(message) = result {
                    Err(message)
                } else {
                    backend.errors.capture(output.finish())
                }
            }
            _ => unreachable!(),
        };
        if command == Some("info") {
            // Source info's prefix/flag fragments do not flush; its final endl
            // does. Keep this after the complete generic renderer.
            use std::io::Write;
            let _ = out.flush();
        }
        backend.errors.finish(result)?;
        Ok(true)
    }

    struct LazyOutput<'a> {
        path: &'a Path,
        errors: crate::util::io::IoErrorOwner,
        file: Option<ToolOutput>,
        buffer_size: usize,
    }
    impl LazyOutput<'_> {
        fn output(&mut self) -> io::Result<&mut ToolOutput> {
            if self.file.is_none() {
                let output = self.errors.capture_io((|| {
                    super::super::blastp::preflight_file_buffer_size(self.buffer_size)?;
                    OutputFile::with_source_profile_path(
                        self.path,
                        Compressor::None,
                        "wb",
                        self.buffer_size,
                    )
                    .map_err(io::Error::other)
                })())?;
                self.file = Some(ToolOutput(output, self.errors.clone()));
            }
            Ok(self.file.as_mut().unwrap())
        }
        fn finish(mut self) -> io::Result<()> {
            self.output()?.0.close().map_err(io::Error::other)
        }
    }
    impl Write for LazyOutput<'_> {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.output()?.write(bytes)
        }
        fn flush(&mut self) -> io::Result<()> {
            self.output()?.flush()
        }
    }
}
