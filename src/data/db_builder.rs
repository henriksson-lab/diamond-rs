use std::io::{self, BufWriter, Seek, SeekFrom, Write};
use std::path::Path;

use rayon::prelude::*;

use super::dmnd::{ReferenceHeader, ReferenceHeader2, CURRENT_DB_VERSION_PROT};
use super::fasta::{self, FastaRecord};
use crate::basic::value::OId;
use crate::basic::value::SequenceType;
use crate::data::taxon_list::TaxonList;
use crate::data::taxonomy::ByteTaxonomy;
use crate::data::taxonomy_nodes::TaxonomyNodes;
use crate::masking::tantan::mask_tantan;
use crate::util::algo::{ExternalSortIoOptions, ExternalSorter};
use crate::util::hash::murmurhash3_x64_128;
use crate::util::sequence::{accession_from_title_bytes, AccessionParsing};

/// Physical size of both reference headers: 40 + 8 + 48 (baseline),
/// or 40 + 8 + 52 for Source EXTRA. The latter's declared record size also
/// includes C ABI tail padding, which the source does not physically write.
/// Matches `out->tell()` after `*out << header << header2` in diamond/src/data/dmnd/dmnd.cpp.
#[cfg(test)]
const HEADER_TOTAL: u64 = ReferenceHeader::SIZE as u64 + 8 + ReferenceHeader2::PAYLOAD_SIZE as u64;

/// `DatabaseFile::make_db` asks `FastaFile::load_seqs` for one billion
/// letters per block.  The loader includes the record that reaches/crosses
/// this limit and then returns the owned block for masking, writing and
/// hashing.
const MAKEDB_BLOCK_LETTERS: u64 = 1_000_000_000;

/// Build a DIAMOND database from FASTA input files.
///
/// This is the Rust implementation of the `makedb` command.
pub fn build_db_paths(
    input_files: &[&Path],
    output_path: &Path,
    seq_type: SequenceType,
) -> io::Result<BuildStats> {
    build_db_with_options_paths(input_files, output_path, seq_type, &BuildOptions::default())
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BuildOptions {
    pub taxonmap: Option<std::path::PathBuf>,
    pub taxonnodes: Option<std::path::PathBuf>,
    pub taxonnames: Option<std::path::PathBuf>,
    pub no_parse_seqids: bool,
    pub file_buffer_size: Option<usize>,
}

/// Concrete command controls kept separate from the stable builder options.
pub(crate) struct BuildRuntime {
    pub temporary_directory: std::path::PathBuf,
    pub mask: bool,
    pub masker: Option<std::sync::Arc<crate::masking::tantan::TantanMasker>>,
    pub ignore_warnings: bool,
    pub no_unlink: bool,
    /// Source command reporting; library builders retain their quiet API.
    pub source_stats_start: Option<std::time::Instant>,
}
impl Default for BuildRuntime {
    fn default() -> Self {
        Self {
            temporary_directory: std::env::temp_dir(),
            mask: true,
            masker: None,
            ignore_warnings: true,
            no_unlink: false,
            source_stats_start: None,
        }
    }
}

pub(crate) fn build_db_with_runtime_paths(
    input_files: &[&Path],
    output_path: &Path,
    seq_type: SequenceType,
    options: &BuildOptions,
    runtime: BuildRuntime,
) -> io::Result<BuildStats> {
    build_db_with_block_limit_runtime_paths(
        input_files,
        output_path,
        seq_type,
        MAKEDB_BLOCK_LETTERS,
        options,
        runtime,
    )
}

pub fn build_db_with_options_paths(
    input_files: &[&Path],
    output_path: &Path,
    seq_type: SequenceType,
    options: &BuildOptions,
) -> io::Result<BuildStats> {
    build_db_with_block_limit_paths(
        input_files,
        output_path,
        seq_type,
        MAKEDB_BLOCK_LETTERS,
        options,
    )
}

fn build_db_with_block_limit_paths(
    input_files: &[&Path],
    output_path: &Path,
    seq_type: SequenceType,
    block_letter_limit: u64,
    options: &BuildOptions,
) -> io::Result<BuildStats> {
    build_db_with_block_limit_runtime_paths(
        input_files,
        output_path,
        seq_type,
        block_letter_limit,
        options,
        BuildRuntime::default(),
    )
}
fn build_db_with_block_limit_runtime_paths(
    input_files: &[&Path],
    output_path: &Path,
    seq_type: SequenceType,
    block_letter_limit: u64,
    options: &BuildOptions,
    runtime: BuildRuntime,
) -> io::Result<BuildStats> {
    debug_assert!(block_letter_limit > 0);
    // FastaFile construction opens and guesses the input format before
    // make_db creates/truncates the destination. Keep those readers alive so
    // pipes and stdin are never reopened after the constructor peek.
    let buffer_size = options.file_buffer_size.unwrap_or(8192);
    let mut inputs = Vec::with_capacity(input_files.len());
    for input_file in input_files {
        inputs.push(fasta::prepare_source_sequence_reader(
            input_file,
            buffer_size,
        )?);
    }
    let mut output =
        DmndBuildWriter::create_with_runtime(output_path, seq_type, options.clone(), runtime)?;
    let load_result = (|| {
        let mut block = Vec::new();
        let mut block_letters = 0u64;

        for input in inputs {
            let consume = |mut record: fasta::ByteFastaRecord| {
                // DMND stores no quality column. Source's database block does
                // not retain FASTQ qualities, so release them before batching.
                record.quality = None;
                if output.sequences == 0 {
                    output.probe_streamed = true;
                    output
                        .protein_probe
                        .observe(&record.sequence)
                        .map_err(io::Error::other)?;
                }
                block_letters = block_letters.saturating_add(record.sequence.len() as u64);
                block.push(record);
                // Upstream's one-pass loader is a do/while: the record that
                // reaches the bound belongs to this block.
                if block_letters >= block_letter_limit {
                    output.write_block_bytes(&mut block)?;
                    block_letters = 0;
                }
                Ok(())
            };
            fasta::for_each_prepared_sequence_reader_bytes(input, seq_type, buffer_size, consume)?;
        }
        output.write_block_bytes(&mut block)
    })();

    // Source's deletion handler covers the sequence loop only. Trailer and
    // taxonomy failures leave the actual, still-incomplete database behind.
    if let Err(error) = load_result {
        // Source closes before removing; a close/flush failure supersedes the
        // original loop error and skips removal. Do not retry a failed flush
        // from BufWriter's destructor.
        let close_result = output.writer.close();
        output.writer.discard_pending();
        close_result?;
        let _ = std::fs::remove_file(output_path);
        return Err(error);
    }
    output.finish()
}

/// Build a DIAMOND database from in-memory sequences.
pub fn build_db_from_records_path(
    records: &[FastaRecord],
    output_path: &Path,
    seq_type: SequenceType,
) -> io::Result<BuildStats> {
    // Keep this convenience API non-mutating.  The production `makedb` path
    // owns its records and masks them in place below, just like C++ mutates
    // each loaded Block before writing it.
    let mut records = records.to_vec();
    write_dmnd(output_path, &mut records, seq_type)
}

/// Statistics from database building.
#[derive(Debug, Clone)]
pub struct BuildStats {
    pub sequences: u64,
    pub letters: u64,
}

fn write_dmnd(
    output_path: &Path,
    records: &mut Vec<FastaRecord>,
    seq_type: SequenceType,
) -> io::Result<BuildStats> {
    let mut output = DmndBuildWriter::create(output_path, seq_type, BuildOptions::default())?;
    output.write_block(records)?;
    output.finish()
}

/// Source Serializer writes fill the remaining staging space before committing
/// a full chunk. Std BufWriter alone instead flushes early for a write that
/// does not fit, changing which bytes survive later taxonomy exceptions.
struct DmndOutputFile {
    #[cfg(any(unix, windows))]
    inner: crate::util::io::cstdio_file::OwnedCstdio,
    // Targets outside the supported POSIX/Windows C runtimes retain the
    // existing raw-file library extension.
    #[cfg(not(any(unix, windows)))]
    inner: Option<std::fs::File>,
}

impl DmndOutputFile {
    fn open(path: &Path) -> io::Result<Self> {
        #[cfg(any(unix, windows))]
        // Source OutputFile uses FileSink's default C buffering. Its separate
        // configured Serializer staging does not call setvbuf on FILE.
        let inner = crate::util::io::cstdio_file::OwnedCstdio::open_path(path, "wb", 0)
            .map_err(io::Error::other)?;
        #[cfg(not(any(unix, windows)))]
        let inner = Some(
            std::fs::File::create(path)
                .map_err(|error| crate::util::io::FileOpenError::into_io(path, error))?,
        );
        Ok(Self { inner })
    }

    fn close(&mut self) -> io::Result<()> {
        #[cfg(any(unix, windows))]
        return self.inner.close().map_err(io::Error::other);
        #[cfg(not(any(unix, windows)))]
        {
            if let Some(file) = self.inner.as_mut() {
                file.flush()?;
            }
            self.inner.take();
            Ok(())
        }
    }
}

impl Write for DmndOutputFile {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        #[cfg(any(unix, windows))]
        {
            self.inner.write_bytes(bytes).map_err(io::Error::other)?;
            Ok(bytes.len())
        }
        #[cfg(not(any(unix, windows)))]
        self.inner
            .as_mut()
            .ok_or_else(|| io::Error::other("Database output is closed"))?
            .write(bytes)
    }

    fn flush(&mut self) -> io::Result<()> {
        // Source Serializer::flush commits to FILE via fwrite; it does not
        // fflush the C stream. Actual fseek/fclose own that effect.
        #[cfg(any(unix, windows))]
        return Ok(());
        #[cfg(not(any(unix, windows)))]
        self.inner
            .as_mut()
            .ok_or_else(|| io::Error::other("Database output is closed"))?
            .flush()
    }
}

impl Seek for DmndOutputFile {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        #[cfg(any(unix, windows))]
        {
            let offset = match position {
                SeekFrom::Start(offset) => offset as i64,
                SeekFrom::Current(offset) | SeekFrom::End(offset) => offset,
            };
            self.inner
                .seek(offset, position)
                .map_err(io::Error::other)?;
            self.stream_position()
        }
        #[cfg(not(any(unix, windows)))]
        self.inner
            .as_mut()
            .ok_or_else(|| io::Error::other("Database output is closed"))?
            .seek(position)
    }

    fn stream_position(&mut self) -> io::Result<u64> {
        #[cfg(any(unix, windows))]
        return self
            .inner
            .position()
            .map(|position| position as u64)
            .map_err(io::Error::other);
        #[cfg(not(any(unix, windows)))]
        self.inner
            .as_mut()
            .ok_or_else(|| io::Error::other("Database output is closed"))?
            .stream_position()
    }
}

struct DmndBufferedFile {
    inner: BufWriter<DmndOutputFile>,
}

impl DmndBufferedFile {
    fn close(&mut self) -> io::Result<()> {
        self.inner.flush()?;
        self.inner.get_mut().close()
    }
    fn discard_pending(self) {
        let (file, _) = self.inner.into_parts();
        drop(file);
    }

    fn write_scalar(&mut self, bytes: &[u8]) -> io::Result<()> {
        // Source Serializer::write<T> does not flush when a scalar fits,
        // including one that exactly fills the staging buffer.
        if bytes.len() <= self.inner.capacity() - self.inner.buffer().len() {
            self.inner.write_all(bytes)
        } else {
            self.write_all(bytes)
        }
    }
}

impl Write for DmndBufferedFile {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.write_all(bytes)?;
        Ok(bytes.len())
    }

    fn write_all(&mut self, mut bytes: &[u8]) -> io::Result<()> {
        if self.inner.capacity() == 0 {
            return self.inner.write_all(bytes);
        }
        while !bytes.is_empty() {
            if self.inner.buffer().len() == self.inner.capacity() {
                self.inner.flush()?;
            }
            let remaining = self.inner.capacity() - self.inner.buffer().len();
            let count = remaining.min(bytes.len());
            self.inner.write_all(&bytes[..count])?;
            bytes = &bytes[count..];
            if count == remaining {
                self.inner.flush()?;
            }
        }
        Ok(())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

impl Seek for DmndBufferedFile {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        self.inner.seek(position)
    }

    fn stream_position(&mut self) -> io::Result<u64> {
        self.inner.flush()?;
        self.inner.get_mut().stream_position()
    }
}

struct DmndBuildWriter {
    writer: DmndBufferedFile,
    seq_type: SequenceType,
    pos_array: Vec<(u64, u32)>,
    current_pos: u64,
    hash: [u8; 16],
    sequences: u64,
    letters: u64,
    options: BuildOptions,
    runtime: BuildRuntime,
    protein_probe: crate::data::sequence_file::ProteinInputProbe,
    probe_streamed: bool,
    accessions: ExternalSorter<(Vec<u8>, OId)>,
    accession_stats: AccessionParsing,
}

impl DmndBuildWriter {
    fn create(
        output_path: &Path,
        seq_type: SequenceType,
        options: BuildOptions,
    ) -> io::Result<Self> {
        Self::create_with_runtime(output_path, seq_type, options, BuildRuntime::default())
    }
    fn create_with_runtime(
        output_path: &Path,
        seq_type: SequenceType,
        options: BuildOptions,
        runtime: BuildRuntime,
    ) -> io::Result<Self> {
        let file = DmndOutputFile::open(output_path)?;
        let mut writer = DmndBufferedFile {
            inner: match options.file_buffer_size {
                Some(size) => BufWriter::with_capacity(size, file),
                None => BufWriter::new(file),
            },
        };
        let mut header = ReferenceHeader::new();
        header.db_version = match seq_type {
            SequenceType::AminoAcid => CURRENT_DB_VERSION_PROT,
            SequenceType::Nucleotide => 4,
        };
        header.write_to(&mut writer)?;
        ReferenceHeader2::new().write_record_to(&mut writer)?;
        // Source initializes its record offset with out->tell(). Serializer's
        // tell flushes the initial headers before any sequence is buffered.
        let current_pos = match writer.stream_position() {
            Ok(position) => position,
            Err(error) => {
                writer.discard_pending();
                return Err(error);
            }
        };
        let accessions = ExternalSorter::in_directory_with_io(
            &runtime.temporary_directory,
            ExternalSortIoOptions {
                no_unlink: runtime.no_unlink,
                file_buffer_size: options.file_buffer_size.unwrap_or(1 << 20),
                ..Default::default()
            },
        );
        Ok(Self {
            writer,
            seq_type,
            pos_array: Vec::new(),
            current_pos,
            hash: [0; 16],
            sequences: 0,
            letters: 0,
            options,
            probe_streamed: false,
            protein_probe: crate::data::sequence_file::ProteinInputProbe::new(
                runtime.ignore_warnings || seq_type != SequenceType::AminoAcid,
            ),
            runtime,
            accessions,
            accession_stats: AccessionParsing::default(),
        })
    }

    fn write_block(&mut self, records: &mut Vec<FastaRecord>) -> io::Result<()> {
        self.write_block_impl(records)
    }
    fn write_block_bytes(&mut self, records: &mut Vec<fasta::ByteFastaRecord>) -> io::Result<()> {
        self.write_block_impl(records)
    }
    fn write_block_impl<T: AsRef<[u8]> + Send>(
        &mut self,
        records: &mut Vec<FastaRecord<T>>,
    ) -> io::Result<()> {
        if records.is_empty() {
            return Ok(());
        }
        // The upstream block is masked in place and dropped immediately after
        // this method.  No second complete encoded copy is retained.
        if self.seq_type == SequenceType::AminoAcid {
            if self.sequences == 0 {
                // In-memory convenience construction has no streaming callback.
                if !self.probe_streamed && !self.runtime.ignore_warnings {
                    for record in records.iter() {
                        self.protein_probe
                            .observe(&record.sequence)
                            .map_err(io::Error::other)?;
                    }
                }
                self.protein_probe
                    .finish_block()
                    .map_err(io::Error::other)?;
            }
            if self.runtime.mask {
                records.par_iter_mut().for_each(|record| {
                    if let Some(masker) = &self.runtime.masker {
                        masker.mask_bit(&mut record.sequence);
                    } else {
                        mask_tantan(&mut record.sequence);
                    }
                });
            }
        }
        for record in records.iter() {
            // `Letter` is i8 and DMND stores its bit pattern verbatim.  Both
            // element types have size/alignment one, and the borrow cannot outlive
            // `record.sequence`.
            let seq_bytes = unsafe {
                std::slice::from_raw_parts(
                    record.sequence.as_ptr().cast::<u8>(),
                    record.sequence.len(),
                )
            };
            // C++ `dmnd.cpp:313-314` throws "File format error: sequence of
            // length 0" — accepting empty sequences here writes a SeqInfo with
            // `seq_len == 0`, which the C++ partition logic in `dmnd.cpp:483`
            // mistakes for the sentinel (`r_next.seq_len == 0` terminates
            // partitioning), corrupting downstream reads of the same DB.
            if seq_bytes.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "File format error: sequence of length 0",
                ));
            }
            self.pos_array
                .push((self.current_pos, seq_bytes.len() as u32));

            self.writer.write_all(&[0xFF])?;
            self.writer.write_all(seq_bytes)?;
            self.writer.write_all(&[0xFF])?;
            self.writer.write_all(record.id.as_ref())?;
            self.writer.write_all(&[0])?;
            self.current_pos += seq_bytes.len() as u64 + record.id.as_ref().len() as u64 + 3;
            self.sequences += 1;
            self.letters += seq_bytes.len() as u64;

            self.hash = murmurhash3_x64_128(seq_bytes, &self.hash);
            self.hash = murmurhash3_x64_128(record.id.as_ref(), &self.hash);

            if self.options.taxonmap.is_some() {
                for accession in accession_from_title_bytes(
                    record.id.as_ref(),
                    !self.options.no_parse_seqids,
                    &mut self.accession_stats,
                ) {
                    self.accessions
                        .try_push((accession, self.sequences - 1))
                        .map_err(io::Error::other)?;
                }
            }
        }
        records.clear();
        Ok(())
    }

    fn finish(mut self) -> io::Result<BuildStats> {
        let finalized = self.finish_buffered();
        let closed = if finalized.is_ok() {
            self.writer.close()
        } else {
            Ok(())
        };
        // Serializer's destructor closes its owner without committing pending
        // bytes on an exception. BufWriter::drop would instead flush them.
        // Keep actual earlier tell/seek/full-buffer writes, discard only the
        // still-pending buffer, and close the same descriptor on either path.
        self.writer.discard_pending();
        closed?;
        let (stats, mut source_stats, hash) = finalized?;
        if let Some(start) = self.runtime.source_stats_start {
            source_stats
                .add_string("Database hash", &crate::util::misc::hex_print(&hash))
                .add_f64(
                    "Total time",
                    start.elapsed().as_millis() as f64 / 1000.0,
                    "s",
                );
            crate::util::log_stream::message_stream()
                .lock()
                .unwrap()
                .write(format!("\n{source_stats}\n"))?;
        }
        Ok(stats)
    }

    fn finish_buffered(&mut self) -> io::Result<(BuildStats, crate::util::table::Table, [u8; 16])> {
        let pos_array_offset = self.current_pos;

        // Position array (u64 pos + u32 seq_len + u32 padding per entry).
        for &(pos, seq_len) in &self.pos_array {
            self.writer.write_scalar(&pos.to_le_bytes())?;
            self.writer.write_scalar(&seq_len.to_le_bytes())?;
            self.writer.write_scalar(&0u32.to_le_bytes())?;
        }
        // Sentinel: (current end-of-seqdata, 0, 0).
        self.writer.write_scalar(&self.current_pos.to_le_bytes())?;
        self.writer.write_scalar(&0u32.to_le_bytes())?;
        self.writer.write_scalar(&0u32.to_le_bytes())?;

        let mut header2 = ReferenceHeader2::new();
        header2.hash = self.hash;
        header2.db_type = self.seq_type as i32;

        // Upstream calls the process-wide `taxonomy.init()` before building
        // the accession mapping. Preserve both lifetime and error order: the
        // sparse name vector remains resident through mapping/node output.
        // Source owns this Table across taxonomy construction and prints it
        // only after the actual database output has been closed.
        let mut source_stats = crate::util::table::Table::new();
        source_stats
            .add_u64("Database sequences", self.sequences, "")
            .add_u64("Database letters", self.letters, "");
        let mut taxonomy = ByteTaxonomy::default();
        if let Some(names_path) = self.options.taxonnames.as_deref() {
            taxonomy.init_source_bytes_io_path(
                names_path,
                self.options.file_buffer_size.unwrap_or(8192),
            )?;
        }

        if let Some(mapping_path) = self.options.taxonmap.as_deref() {
            header2.taxon_array_offset = self.writer.stream_position()?;
            let errors = self.accessions.io_options().errors;
            let build = errors.finish_with_kind(
                TaxonList::write_from_source_file_bytes_with_buffer_size_path(
                    &mut self.accessions,
                    mapping_path,
                    self.sequences,
                    self.options.no_parse_seqids,
                    &mut self.writer,
                    self.options.file_buffer_size.unwrap_or(8192),
                ),
                io::ErrorKind::InvalidData,
            )?;
            // Source measures this through out->tell(), committing mapping
            // bytes before a subsequent nodes-file constructor can fail.
            header2.taxon_array_size = self.writer.stream_position()? - header2.taxon_array_offset;
            source_stats
                .add_u64(
                    "Accessions in database",
                    build.stats.accessions_in_database as u64,
                    "",
                )
                .add_u64(
                    "Entries in accession to taxid file",
                    build.stats.entries_in_accession_to_taxid_file as u64,
                    "",
                )
                .add_u64(
                    "Database accessions mapped to taxid",
                    build.stats.database_accessions_mapped_to_taxid as u64,
                    "",
                )
                .add_u64(
                    "Database sequences mapped to taxid",
                    build.stats.database_sequences_mapped_to_taxid as u64,
                    "",
                );
        }

        if let Some(nodes_path) = self.options.taxonnodes.as_deref() {
            let nodes = TaxonomyNodes::from_nodes_dmp_source_io_path(
                nodes_path,
                self.options.file_buffer_size.unwrap_or(8192),
            )?;
            header2.taxon_nodes_offset = self.writer.stream_position()?;
            nodes.write_to(&mut self.writer)?;
        }

        if self.options.taxonnames.is_some() {
            header2.taxon_names_offset = self.writer.stream_position()?;
            write_string_vector(&mut self.writer, &taxonomy.name)?;
        }

        self.writer.flush()?;

        // Source seeks the same output owner to zero. Both seek and close
        // flush at their real successful boundaries; no pathname is reopened.
        let mut header = ReferenceHeader::new();
        header.db_version = match self.seq_type {
            SequenceType::AminoAcid => CURRENT_DB_VERSION_PROT,
            SequenceType::Nucleotide => 4,
        };
        header.sequences = self.sequences;
        header.letters = self.letters;
        header.pos_array_offset = pos_array_offset;
        self.writer.seek(SeekFrom::Start(0))?;
        header.write_to(&mut self.writer)?;
        header2.write_record_to(&mut self.writer)?;
        self.writer.flush()?;

        Ok((
            BuildStats {
                sequences: self.sequences,
                letters: self.letters,
            },
            source_stats,
            header2.hash,
        ))
    }
}

fn write_string_vector<W: Write, T: AsRef<[u8]>>(writer: &mut W, values: &[T]) -> io::Result<()> {
    writer.write_all(&(values.len() as u32).to_le_bytes())?;
    for value in values {
        writer.write_all(value.as_ref())?;
        writer.write_all(&[0])?;
    }
    Ok(())
}

pub fn build_db(
    input_files: &[&str],
    output_path: &str,
    seq_type: SequenceType,
) -> io::Result<BuildStats> {
    let paths: Vec<&Path> = input_files.iter().map(|path| Path::new(path)).collect();
    build_db_paths(&paths, Path::new(output_path), seq_type)
}

pub(crate) fn build_db_with_runtime(
    input_files: &[&str],
    output_path: &str,
    seq_type: SequenceType,
    options: &BuildOptions,
    runtime: BuildRuntime,
) -> io::Result<BuildStats> {
    let paths: Vec<&Path> = input_files.iter().map(|path| Path::new(path)).collect();
    build_db_with_runtime_paths(&paths, Path::new(output_path), seq_type, options, runtime)
}

pub fn build_db_with_options(
    input_files: &[&str],
    output_path: &str,
    seq_type: SequenceType,
    options: &BuildOptions,
) -> io::Result<BuildStats> {
    let paths: Vec<&Path> = input_files.iter().map(|path| Path::new(path)).collect();
    build_db_with_options_paths(&paths, Path::new(output_path), seq_type, options)
}

fn build_db_with_block_limit(
    input_files: &[&str],
    output_path: &str,
    seq_type: SequenceType,
    block_letter_limit: u64,
    options: &BuildOptions,
) -> io::Result<BuildStats> {
    let paths: Vec<&Path> = input_files.iter().map(|path| Path::new(path)).collect();
    build_db_with_block_limit_paths(
        &paths,
        Path::new(output_path),
        seq_type,
        block_letter_limit,
        options,
    )
}

fn build_db_with_block_limit_runtime(
    input_files: &[&str],
    output_path: &str,
    seq_type: SequenceType,
    block_letter_limit: u64,
    options: &BuildOptions,
    runtime: BuildRuntime,
) -> io::Result<BuildStats> {
    let paths: Vec<&Path> = input_files.iter().map(|path| Path::new(path)).collect();
    build_db_with_block_limit_runtime_paths(
        &paths,
        Path::new(output_path),
        seq_type,
        block_letter_limit,
        options,
        runtime,
    )
}

pub fn build_db_from_records(
    records: &[FastaRecord],
    output_path: &str,
    seq_type: SequenceType,
) -> io::Result<BuildStats> {
    build_db_from_records_path(records, Path::new(output_path), seq_type)
}
#[cfg(test)]
mod tests {
    #[cfg(unix)]
    #[test]
    fn source_builder_retains_os_paths_across_taxonomy_and_dmnd_owners() {
        use std::os::unix::ffi::OsStringExt;
        let root = std::env::temp_dir().join(format!(
            "diamond-data-path-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ));
        std::fs::create_dir(&root).unwrap();
        let (directory, raw_names) =
            crate::util::test_os_path_or_unicode(&root, b"files\xff", "files-é");
        std::fs::create_dir(&directory).unwrap();
        let component = |raw: &[u8], unicode: &str| {
            if raw_names {
                directory.join(std::ffi::OsString::from_vec(raw.to_vec()))
            } else {
                directory.join(unicode)
            }
        };
        let input = component(b"input\x80.faa", "input-λ.faa");
        let output = component(b"output\xfe.dmnd", "output-é.dmnd");
        let mapping = component(b"map\xff.tsv", "map-λ.tsv");
        let nodes = component(b"nodes\x80.dmp", "nodes-é.dmp");
        let names = component(b"names\xfe.dmp", "names-λ.dmp");
        std::fs::write(
            &input,
            b">ACC raw\xff\nMKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYT\n",
        )
        .unwrap();
        std::fs::write(&mapping, b"accession.version\ttaxid\nACC\t10\n").unwrap();
        std::fs::write(
            &nodes,
            b"1\t|\t1\t|\tno rank\t|\t\n10\t|\t1\t|\tspecies\t|\t\n",
        )
        .unwrap();
        std::fs::write(&names, b"10\t|\tName\xff\t|\t\t|\tscientific name\t|\n").unwrap();
        let options = super::BuildOptions {
            taxonmap: Some(mapping),
            taxonnodes: Some(nodes),
            taxonnames: Some(names),
            file_buffer_size: Some(64),
            ..Default::default()
        };
        let runtime = super::BuildRuntime {
            temporary_directory: directory.clone(),
            mask: false,
            ..Default::default()
        };
        let stats = super::build_db_with_runtime_paths(
            &[input.as_path()],
            &output,
            crate::basic::value::SequenceType::AminoAcid,
            &options,
            runtime,
        )
        .unwrap();
        assert_eq!(stats.sequences, 1);
        let mut db = crate::data::dmnd::DatabaseFile::open_source_with_options_io(
            &output,
            crate::data::dmnd::DatabaseOpenOptions {
                load_taxon_mapping: true,
                load_taxon_nodes: true,
                load_taxon_scientific_names: true,
                file_buffer_size: Some(64),
                ..Default::default()
            },
            false,
        )
        .unwrap();
        assert_eq!(db.record_bytes_io(0).unwrap().1, b"ACC raw\xff");
        assert_eq!(db.taxon_scientific_name_bytes(10), b"Name\xff");
        drop(db);
        std::fs::remove_dir_all(root).unwrap();
    }

    use super::*;
    use crate::data::dmnd::{DatabaseFile, DatabaseOpenOptions, MAGIC_NUMBER};
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEMP_PATH_NONCE: AtomicU64 = AtomicU64::new(0);

    fn temp_path(label: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond_rs_{label}_{}_{}",
            std::process::id(),
            TEMP_PATH_NONCE.fetch_add(1, Ordering::Relaxed)
        ))
    }

    #[test]
    fn streamed_builder_hashes_and_writes_full_non_utf8_title_bytes() {
        let input = temp_path("raw_title.faa");
        let output = temp_path("raw_title.dmnd");
        let title = b"q\xff description";
        let mut fasta = b">".to_vec();
        fasta.extend_from_slice(title);
        fasta.extend_from_slice(b"\nARNDCQEGHILKMFPSTWYV\n");
        std::fs::write(&input, &fasta).unwrap();
        build_db_with_block_limit_runtime(
            &[input.to_str().unwrap()],
            output.to_str().unwrap(),
            SequenceType::AminoAcid,
            64,
            &BuildOptions::default(),
            BuildRuntime {
                mask: false,
                ignore_warnings: true,
                ..Default::default()
            },
        )
        .unwrap();
        let mut db =
            DatabaseFile::open_with_options_io(&output, DatabaseOpenOptions::default()).unwrap();
        let (sequence, id) = db.record_bytes_io(0).unwrap();
        assert_eq!(id, title);
        assert!(db.record_io(0).is_err());
        let raw: Vec<u8> = sequence.into_iter().map(|letter| letter as u8).collect();
        let expected = murmurhash3_x64_128(title, &murmurhash3_x64_128(&raw, &[0; 16]));
        assert_eq!(db.header2.hash, expected);
        drop(db);
        std::fs::remove_file(input).unwrap();
        std::fs::remove_file(output).unwrap();
    }

    #[test]
    fn test_build_stats() {
        let records = vec![
            FastaRecord {
                quality: None,
                id: "seq1".to_string(),
                sequence: vec![0, 1, 2, 3],
            },
            FastaRecord {
                quality: None,
                id: "seq2".to_string(),
                sequence: vec![4, 5, 6],
            },
        ];
        let original = records.clone();

        let tmp = std::env::temp_dir().join("test_diamond.dmnd");
        let stats = build_db_from_records(&records, tmp.to_str().unwrap(), SequenceType::AminoAcid)
            .unwrap();

        assert_eq!(stats.sequences, 2);
        assert_eq!(stats.letters, 7);
        for (record, original) in records.iter().zip(&original) {
            assert_eq!(record.id, original.id);
            assert_eq!(record.sequence, original.sequence);
        }

        // Verify it starts with the magic number
        let data = std::fs::read(&tmp).unwrap();
        assert!(data.len() > 8);
        let magic = u64::from_le_bytes(data[0..8].try_into().unwrap());
        assert_eq!(magic, MAGIC_NUMBER);

        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn streamed_blocks_produce_the_same_database_as_owned_records() {
        // The limit deliberately falls inside this small, varied input.  This
        // exercises the production block loop without manufacturing a huge or
        // duplicated fixture.
        let fasta_path = temp_path("streamed_blocks.faa");
        let streamed_path = temp_path("streamed_blocks.dmnd");
        let owned_path = temp_path("owned_records.dmnd");
        std::fs::write(
            &fasta_path,
            b">alpha first\nARND\n>beta\nCQE\n>gamma\nGHILK\n",
        )
        .unwrap();

        let stats = build_db_with_block_limit(
            &[fasta_path.to_str().unwrap()],
            streamed_path.to_str().unwrap(),
            SequenceType::AminoAcid,
            5,
            &BuildOptions::default(),
        )
        .unwrap();
        let records = fasta::read_fasta_file(&fasta_path, SequenceType::AminoAcid).unwrap();
        let owned_stats = build_db_from_records(
            &records,
            owned_path.to_str().unwrap(),
            SequenceType::AminoAcid,
        )
        .unwrap();

        assert_eq!(stats.sequences, 3);
        assert_eq!(stats.letters, 12);
        assert_eq!(owned_stats.sequences, stats.sequences);
        assert_eq!(owned_stats.letters, stats.letters);
        assert_eq!(
            std::fs::read(&streamed_path).unwrap(),
            std::fs::read(&owned_path).unwrap()
        );

        let _ = std::fs::remove_file(fasta_path);
        let _ = std::fs::remove_file(streamed_path);
        let _ = std::fs::remove_file(owned_path);
    }

    #[test]
    fn streaming_failure_removes_incomplete_database() {
        let fasta_path = temp_path("invalid_stream.faa");
        let output_path = temp_path("invalid_stream.dmnd");
        std::fs::write(&fasta_path, b">valid\nARND\n>invalid\nARN?\n").unwrap();

        let error = build_db_with_block_limit(
            &[fasta_path.to_str().unwrap()],
            output_path.to_str().unwrap(),
            SequenceType::AminoAcid,
            2,
            &BuildOptions::default(),
        )
        .unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(!output_path.exists());

        let _ = std::fs::remove_file(fasta_path);
    }

    #[test]
    fn taxonomy_metadata_is_serialized_and_read_back_exactly() {
        let fasta_path = temp_path("taxonomy.faa");
        let mapping_path = temp_path("taxonomy.tsv");
        let nodes_path = temp_path("nodes.dmp");
        let names_path = temp_path("names.dmp");
        let output_path = temp_path("taxonomy.dmnd");
        std::fs::write(
            &fasta_path,
            b">sp|P12345.2|alpha >UniRef90_Q999.1\nMKTAYIAK\n>tr|ACC3.4|beta\nVLSPADKT\n>unmapped.7 gamma\nMPEPTIDE\n",
        )
        .unwrap();
        std::fs::write(
            &mapping_path,
            b"accession.version\ttaxid\nP12345\t10\nQ999\t20\nACC3\t30\nunused\t999\n",
        )
        .unwrap();
        std::fs::write(
            &nodes_path,
            b"1\t|\t1\t|\tno rank\t|\t\n10\t|\t1\t|\tspecies\t|\t\n20\t|\t1\t|\tgenus\t|\t\n30\t|\t20\t|\tstrain\t|\t\n",
        )
        .unwrap();
        std::fs::write(
            &names_path,
            b"1\t|\troot\t|\t\t|\tscientific name\t|\n10\t|\tAlpha species\t|\t\t|\tscientific name\t|\n10\t|\talpha\t|\t\t|\tcommon name\t|\n20\t|\tBeta genus\t|\t\t|\tscientific name\t|\n30\t|\tBeta strain\t|\t\t|\tscientific name\t|\n31\t|\tOpaque\xff taxon\x80\t|\t\t|\tscientific name\t|\n",
        )
        .unwrap();

        let stats = build_db_with_options(
            &[fasta_path.to_str().unwrap()],
            output_path.to_str().unwrap(),
            SequenceType::AminoAcid,
            &BuildOptions {
                taxonmap: Some(mapping_path.to_path_buf()),
                taxonnodes: Some(nodes_path.to_path_buf()),
                taxonnames: Some(names_path.to_path_buf()),
                no_parse_seqids: false,
                file_buffer_size: None,
            },
        )
        .unwrap();
        assert_eq!(stats.sequences, 3);
        assert_eq!(stats.letters, 24);

        let db = DatabaseFile::open_with_options(
            &output_path,
            DatabaseOpenOptions {
                load_taxon_mapping: true,
                load_taxon_nodes: true,
                load_taxon_scientific_names: true,
                ..DatabaseOpenOptions::default()
            },
        )
        .unwrap();
        assert!(db.header2.taxon_array_offset < db.header2.taxon_nodes_offset);
        assert!(db.header2.taxon_nodes_offset < db.header2.taxon_names_offset);
        assert_eq!(db.taxids(0).unwrap(), vec![10, 20]);
        assert_eq!(db.taxids(1).unwrap(), vec![30]);
        assert_eq!(db.taxids(2).unwrap(), Vec::<i32>::new());
        assert_eq!(db.get_parent(30).unwrap(), 20);
        assert_eq!(
            db.rank(10).unwrap(),
            crate::data::taxonomy::Rank::Species as i32
        );
        assert_eq!(db.taxon_scientific_name(10), "Alpha species");
        assert_eq!(db.taxon_scientific_name(11), "11");
        assert_eq!(db.taxon_scientific_name_bytes(31), b"Opaque\xff taxon\x80");
        assert!(db.try_taxon_scientific_name(31).is_err());

        for path in [
            fasta_path,
            mapping_path,
            nodes_path,
            names_path,
            output_path,
        ] {
            let _ = std::fs::remove_file(path);
        }
    }

    #[test]
    fn missing_taxonomy_names_discards_pending_sequence_and_trailer_bytes() {
        let output_path = temp_path("missing_names.dmnd");
        let names_path = temp_path("missing_names.dmp");
        let mut output = DmndBuildWriter::create(
            &output_path,
            SequenceType::AminoAcid,
            BuildOptions {
                taxonnames: Some(names_path.to_path_buf()),
                file_buffer_size: Some(4 << 20),
                ..BuildOptions::default()
            },
        )
        .unwrap();
        assert_eq!(output.current_pos, HEADER_TOTAL);
        // tell commits Serializer bytes to FILE, but ftell need not commit
        // the C buffer to the descriptor. Check the logical constructor bytes
        // after the actual owner closes on the taxonomy error below.
        let mut initial_headers = Vec::new();
        ReferenceHeader::new()
            .write_to(&mut initial_headers)
            .unwrap();
        ReferenceHeader2::new()
            .write_record_to(&mut initial_headers)
            .unwrap();
        output
            .write_block(&mut vec![FastaRecord {
                id: "P1".into(),
                sequence: vec![0, 1, 2, 3],
                quality: None,
            }])
            .unwrap();
        let error = output.finish().unwrap_err();
        assert!(crate::util::io::file_open_error::find(&error).is_some());
        // Source taxonomy.init fails before the mapping tell/seek boundary:
        // only its initial Serializer-committed headers survive fclose.
        assert_eq!(std::fs::read(&output_path).unwrap(), initial_headers);
        let header = ReferenceHeader::read_from(&mut initial_headers.as_slice()).unwrap();
        assert_eq!((header.sequences, header.letters), (0, 0));
        // This field is indeterminate in Source's unfinished Header1; Rust's
        // explicit zero remains the documented safe representation.
        assert_eq!(header.pos_array_offset, 0);
        std::fs::remove_file(output_path).unwrap();
    }

    #[test]
    fn taxonomy_failure_keeps_chunks_already_written_before_finalization() {
        let output_path = temp_path("missing_names_written_chunk.dmnd");
        let names_path = temp_path("missing_names_written_chunk.dmp");
        let mut output = DmndBuildWriter::create(
            &output_path,
            SequenceType::Nucleotide,
            BuildOptions {
                taxonnames: Some(names_path.to_path_buf()),
                file_buffer_size: Some(128),
                ..BuildOptions::default()
            },
        )
        .unwrap();
        output
            .write_block(&mut vec![FastaRecord {
                id: "P1".into(),
                sequence: vec![0; 256],
                quality: None,
            }])
            .unwrap();
        // Query the actual FILE's logical position without fseek/fflush and
        // without flushing the pending Serializer tail.
        assert_eq!(
            output.writer.inner.get_mut().stream_position().unwrap(),
            HEADER_TOTAL + 2 * 128
        );
        output.finish().unwrap_err();
        let committed = std::fs::read(&output_path).unwrap();
        assert_eq!(committed.len(), HEADER_TOTAL as usize + 2 * 128);
        assert_eq!(committed[HEADER_TOTAL as usize], 0xff);
        assert!(committed[HEADER_TOTAL as usize + 1..]
            .iter()
            .all(|&byte| byte == 0));
        std::fs::remove_file(output_path).unwrap();
    }

    #[test]
    fn missing_nodes_retains_committed_mapping_and_file_open_identity() {
        let output_path = temp_path("mapped_missing_nodes.dmnd");
        let mapping_path = temp_path("mapped_missing_nodes.tsv");
        let nodes_path = temp_path("mapped_missing_nodes.dmp");
        std::fs::write(&mapping_path, b"accession.version\ttaxid\nP1\t10\n").unwrap();
        let mut output = DmndBuildWriter::create(
            &output_path,
            SequenceType::AminoAcid,
            BuildOptions {
                taxonmap: Some(mapping_path.to_path_buf()),
                taxonnodes: Some(nodes_path.to_path_buf()),
                file_buffer_size: Some(4 << 20),
                ..BuildOptions::default()
            },
        )
        .unwrap();
        output
            .write_block(&mut vec![FastaRecord {
                id: "P1".into(),
                sequence: vec![0, 1, 2, 3],
                quality: None,
            }])
            .unwrap();
        let error = output.finish().unwrap_err();
        let file_open = crate::util::io::file_open_error::find(&error).unwrap();
        assert_eq!(file_open.path, nodes_path);
        assert!(file_open.raw_os_error.is_some());
        let bytes = std::fs::read(&output_path).unwrap();
        let record_bytes = 4 + "P1".len() + 3;
        let mapping_offset = HEADER_TOTAL as usize + record_bytes + 2 * 16;
        // One mapped taxid: Source TaxonList's varuint count and taxid.
        assert_eq!(&bytes[mapping_offset..], &[3, 21]);
        let header = ReferenceHeader::read_from(&mut bytes.as_slice()).unwrap();
        assert_eq!((header.sequences, header.letters), (0, 0));
        // Source leaves this unfinished field indeterminate, unlike Rust's
        // defined safe zero; all committed payload bytes are asserted above.
        assert_eq!(header.pos_array_offset, 0);
        std::fs::remove_file(mapping_path).unwrap();
        std::fs::remove_file(output_path).unwrap();
    }

    #[test]
    fn taxonomy_failure_retains_incomplete_database_after_sequence_loop() {
        let fasta_path = temp_path("taxonomy_error.faa");
        let mapping_path = temp_path("taxonomy_error.tsv");
        let output_path = temp_path("taxonomy_error.dmnd");
        std::fs::write(&fasta_path, b">P1\nMPEPTIDE\n").unwrap();
        std::fs::write(&mapping_path, b"bad\theader\nP1\t10\n").unwrap();

        let error = build_db_with_options(
            &[fasta_path.to_str().unwrap()],
            output_path.to_str().unwrap(),
            SequenceType::AminoAcid,
            &BuildOptions {
                taxonmap: Some(mapping_path.to_path_buf()),
                ..BuildOptions::default()
            },
        )
        .unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        let bytes = std::fs::read(&output_path).unwrap();
        let header = ReferenceHeader::read_from(&mut bytes.as_slice()).unwrap();
        assert_eq!(header.magic_number, MAGIC_NUMBER);
        assert_eq!((header.sequences, header.letters), (0, 0));
        // Source leaves this field uninitialized until successful finalization;
        // Rust retains its defined zero initialization on the failure path.
        assert_eq!(header.pos_array_offset, 0);
        assert_eq!(bytes.len(), HEADER_TOTAL as usize + 13 + 2 * 16);
        assert_eq!(
            &bytes[bytes.len() - 16..bytes.len() - 8],
            &(HEADER_TOTAL + 13).to_le_bytes()
        );

        let _ = std::fs::remove_file(fasta_path);
        let _ = std::fs::remove_file(mapping_path);
        let _ = std::fs::remove_file(output_path);
    }
}
