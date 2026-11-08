//! Concrete native adapters for the translated intermediate block merger.
use crate::util::io::source_c_path;
use std::collections::{BTreeSet, HashMap};
use std::io;

use super::blastp::{BlastpConfig, NativeSearchCapture, NativeSearchWriter};
use super::search_output::SearchSubjects;
use crate::align::hsp::{Hsp, Match};
use crate::basic::value::{AlignMode, Letter};
use crate::data::block::Block;
use crate::data::fasta::ByteFastaRecord as FastaRecord;
use crate::dp::swipe::HspValues;
use crate::output::intermediate::{IntermediateRecord, OutputFormat};
use crate::output::join_blocks::{
    join_query, JoinBlocksConfig, JoinBlocksDb, JoinQueryCallbacks, JoinStatistics,
    JoinTargetAction,
};
use crate::output::target_culling::{CullingResult, TargetCulling, TargetCullingConfig};
use crate::stats::score_matrix::ScoreMatrix;

#[path = "fasta_reference.rs"]
mod fasta_reference;

#[path = "cluster_search.rs"]
mod cluster_search;
#[path = "native_double_indexed.rs"]
mod native_double_indexed;

#[path = "global_ranking.rs"]
mod global_ranking;
#[path = "multiprocessing_search.rs"]
mod multiprocessing;
pub(crate) use cluster_search::{
    run as run_cluster_edges, run_with_owner as run_cluster_edges_with_owner, ClusterOwnerKind,
    ClusterSearchOwner,
};
#[cfg(feature = "dna")]
#[path = "dna_search.rs"]
mod dna_search;
#[cfg(feature = "dna")]
pub(crate) use dna_search::run as run_nucleotide_files;
#[cfg(feature = "mcl")]
#[path = "mcl_search.rs"]
mod mcl_search;
#[cfg(feature = "mcl")]
pub(crate) use mcl_search::run as run_mcl_search;

struct ReferenceInfo {
    oid: u64,
    title: Vec<u8>,
    len: usize,
    sequence: Option<std::sync::Arc<[Letter]>>,
    self_aln_score: Option<f64>,
    taxids: Vec<crate::basic::value::TaxId>,
}

/// Only actual matched reference metadata survives the reference block owner.
/// Dictionary keys remain separate from source OIDs used for merge precedence.
pub(crate) struct ReferenceDictionary {
    entries: Vec<ReferenceInfo>,
    keep_sequences: bool,
    by_oid: HashMap<u64, usize>,
    sequence_count: u64,
    taxonomy: Option<crate::data::taxonomy::TaxonomyTree>,
    taxids_loaded: usize,
    seen_taxids: BTreeSet<crate::basic::value::TaxId>,
}
impl ReferenceDictionary {
    pub(crate) fn new(sequence_count: u64, keep_sequences: bool) -> Self {
        Self {
            entries: Vec::new(),
            keep_sequences,
            by_oid: HashMap::new(),
            sequence_count,
            taxonomy: None,
            taxids_loaded: 0,
            seen_taxids: BTreeSet::new(),
        }
    }
    fn insert(
        &mut self,
        oid: u64,
        record: &FastaRecord,
        self_aln_score: Option<f64>,
        compact_title: bool,
    ) -> usize {
        if let Some(&id) = self.by_oid.get(&oid) {
            return id;
        }
        let id = self.entries.len();
        self.entries.push(ReferenceInfo {
            oid,
            // Block::dict_id retains only seqid when SAM lacks FULL_TITLES.
            // A single direct reference block keeps its loaded full title.
            title: if compact_title {
                crate::output::format::print_title_bytes(&record.id, false, false, b"\x01", false)
            } else {
                record.id.clone()
            },
            len: record.sequence.len(),
            // Output uses the database sequence, not the masked extension copy.
            sequence: self.keep_sequences.then(|| record.sequence.clone().into()),
            self_aln_score,
            taxids: Vec::new(),
        });
        self.by_oid.insert(oid, id);
        id
    }
    fn species_rank_ids(&self, records: &[IntermediateRecord]) -> BTreeSet<u32> {
        records
            .first()
            .and_then(|record| self.entries.get(record.target_dict_id as usize))
            .and_then(|entry| {
                self.taxonomy.as_ref().map(|tree| {
                    tree.rank_taxids(&entry.taxids, "species")
                        .into_iter()
                        .map(|taxid| taxid as u32)
                        .collect()
                })
            })
            .unwrap_or_default()
    }
    fn remap_direct(
        &mut self,
        capture: &mut NativeSearchCapture,
        references: &[FastaRecord],
        oids: &[u64],
    ) {
        for matches in &mut capture.matches {
            for target in matches {
                let local = capture
                    .target_oids
                    .get(target.target_block_id as usize)
                    .copied()
                    .unwrap_or(target.target_block_id as u64) as usize;
                let original = oids[local];
                target.target_block_id = self.insert(
                    original,
                    &references[local],
                    capture
                        .target_self_aln_scores
                        .get(target.target_block_id as usize)
                        .copied(),
                    false,
                ) as u32;
                target.target_oid = original;
            }
        }
    }
    fn write_block(
        &mut self,
        capture: &mut NativeSearchCapture,
        references: &[FastaRecord],
        oids: &[u64],
        output_format: &OutputFormat,
        writer: &mut impl std::io::Write,
    ) -> io::Result<()> {
        for (query_id, matches) in capture.matches.iter_mut().enumerate() {
            let mut bytes = Vec::new();
            for target in matches.drain(..) {
                let local = capture
                    .target_oids
                    .get(target.target_block_id as usize)
                    .copied()
                    .unwrap_or(target.target_block_id as u64) as usize;
                let original = oids[local];
                let dict = self.insert(
                    original,
                    &references[local],
                    capture
                        .target_self_aln_scores
                        .get(target.target_block_id as usize)
                        .copied(),
                    output_format.code == crate::output::format::FormatCode::Sam,
                );
                for hsp in target.hsps {
                    IntermediateRecord::write(
                        &mut bytes,
                        &hsp,
                        query_id as u32,
                        dict as i64,
                        original,
                        output_format,
                    )?;
                }
            }
            if !bytes.is_empty() {
                writer.write_all(&(query_id as u32).to_ne_bytes())?;
                let size = u32::try_from(bytes.len())
                    .map_err(|_| io::Error::other("intermediate query exceeds uint32"))?;
                writer.write_all(&size.to_ne_bytes())?;
                writer.write_all(&bytes)?;
            }
        }
        IntermediateRecord::finish_file(writer)
    }
}
impl SearchSubjects for ReferenceDictionary {
    fn title_bytes(&self, id: usize) -> io::Result<&[u8]> {
        Ok(&self.entries[id].title)
    }
    fn length(&self, id: usize) -> usize {
        self.entries[id].len
    }
    fn original_oid(&self, id: usize) -> u64 {
        self.entries[id].oid
    }
    fn sequence_count(&self) -> u64 {
        self.sequence_count
    }
    fn self_aln_score(&self, id: usize) -> Option<f64> {
        self.entries[id].self_aln_score
    }
    fn sequence(&self, id: usize) -> Option<&[Letter]> {
        self.entries[id].sequence.as_deref()
    }
    fn taxids(&self, id: usize) -> &[crate::basic::value::TaxId] {
        &self.entries[id].taxids
    }
    fn taxonomy(&self) -> Option<&crate::data::taxonomy::TaxonomyTree> {
        self.taxonomy.as_ref()
    }
}
impl JoinBlocksDb for ReferenceDictionary {
    fn oid(&self, id: i64, _block: i64) -> u64 {
        self.entries[id as usize].oid
    }
}

struct MergedQuery<'a> {
    culling: TargetCulling,
    culling_config: TargetCullingConfig,
    config: &'a BlastpConfig,
    matrix: &'a ScoreMatrix,
    dictionary: &'a ReferenceDictionary,
    frames: &'a [Vec<Letter>],
    source_len: usize,
    format: OutputFormat,
    matches: Vec<Match>,
}
impl JoinQueryCallbacks for MergedQuery<'_> {
    fn cull(&mut self, records: &[IntermediateRecord], _block: i64, _oid: u64) -> JoinTargetAction {
        match self.culling.cull_records(
            records,
            &self.dictionary.species_rank_ids(records),
            &self.culling_config,
            |score| self.matrix.bitscore(score as f64),
        ) {
            CullingResult::Include => JoinTargetAction::Include,
            CullingResult::Next => JoinTargetAction::Next,
            CullingResult::Finished => JoinTargetAction::Finished,
        }
    }
    fn write_match(
        &mut self,
        _out: &mut Vec<u8>,
        record: &IntermediateRecord,
        _block: i64,
        oid: u64,
        target_num: usize,
        _hsp_num: usize,
    ) -> io::Result<()> {
        let translated = self.config.translated_query_layout.is_some();
        let frame = record.frame(
            self.source_len as i32,
            if translated {
                AlignMode::BLASTX
            } else {
                AlignMode::BLASTP
            },
        ) as usize;
        let mut hsp = Hsp::from_intermediate_record(
            record,
            self.source_len as u32,
            self.frames.get(frame).map_or(0, Vec::len) as i32,
            self.dictionary.length(record.target_dict_id as usize) as i32,
            self.format.hsp_values,
            if translated {
                AlignMode::BLASTX
            } else {
                AlignMode::BLASTP
            },
            self.matrix,
        );
        hsp.parse_with_query(
            self.frames,
            self.source_len as i32,
            self.dictionary
                .sequence(record.target_dict_id as usize)
                .unwrap_or(&[]),
            self.format.hsp_values.any(HspValues::TRANSCRIPT),
            false,
            translated,
            self.matrix,
        )
        .map_err(io::Error::other)?;
        if self.matches.len() == target_num {
            let mut target = Match::new(record.target_dict_id as u32, oid);
            target.filter_score = hsp.score;
            target.filter_evalue = hsp.evalue;
            self.matches.push(target);
        }
        self.matches[target_num].hsps.push(hsp);
        Ok(())
    }
    fn add(&mut self, records: &[IntermediateRecord], _block: i64, _oid: u64) {
        self.culling.add_records(
            records,
            &self.dictionary.species_rank_ids(records),
            &self.culling_config,
            |score| self.matrix.bitscore(score as f64),
        );
    }
}

pub(crate) fn merge_query(
    buffers: &[Vec<u8>],
    frames: &[Vec<Letter>],
    source_len: usize,
    dictionary: &ReferenceDictionary,
    matrix: &ScoreMatrix,
    config: &BlastpConfig,
    format: OutputFormat,
) -> io::Result<Vec<Match>> {
    let mut query = MergedQuery {
        culling: TargetCulling::get(config.max_target_seqs, config.range_culling),
        culling_config: TargetCullingConfig {
            query_range_culling: config.range_culling,
            query_range_cover: config.range_cover,
            taxon_k: config.taxon_k.min(u32::MAX as u64) as u32,
            global_ranking_targets: config.global_ranking_targets,
            toppercent: config.toppercent,
        },
        config,
        matrix,
        dictionary,
        frames,
        source_len,
        format,
        matches: Vec::new(),
    };
    let mut statistics = JoinStatistics {
        aligned: 1,
        ..Default::default()
    };
    join_query(
        buffers,
        &mut Vec::new(),
        &mut statistics,
        dictionary,
        &format,
        &JoinBlocksConfig {
            use_score_order: config.toppercent.is_some(),
            ..JoinBlocksConfig::default()
        },
        &mut query,
    )?;
    let mut totals = crate::basic::statistics::Statistics::new();
    use crate::basic::statistics::StatValue;
    totals.inc(StatValue::Pairwise, statistics.pairwise as i64);
    totals.inc(StatValue::Matches, statistics.matches as i64);
    totals.inc(StatValue::Aligned, statistics.aligned as i64);
    crate::basic::statistics::merge_statistics(&totals);
    Ok(query.matches)
}

/// The parser and count/rewind phase retain the very same decoded Source owner.
#[derive(Clone)]
struct PreparedQueryInput(std::sync::Arc<std::sync::Mutex<crate::util::io::InputFile>>);
impl std::io::Read for PreparedQueryInput {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        std::io::Read::read(
            &mut *self
                .0
                .lock()
                .map_err(|_| io::Error::other("query source lock poisoned"))?,
            bytes,
        )
    }
}

/// Forward-only parser with one-record backpressure. The existing native
/// parser owns compression/FASTA/FASTQ rules; block consumers own the budget.
pub(crate) struct RecordStream {
    receiver: Option<std::sync::mpsc::Receiver<io::Result<FastaRecord>>>,
    request: Option<std::sync::mpsc::SyncSender<()>>,
    worker: Option<std::thread::JoinHandle<()>>,
    rewind_file: std::sync::Arc<std::sync::Mutex<Option<PreparedQueryInput>>>,
    path: std::path::PathBuf,
    sequence_type: crate::basic::value::SequenceType,
    keep_quality: bool,
    file_buffer_size: usize,
    ready: bool,
}
impl RecordStream {
    fn open(
        path: std::path::PathBuf,
        sequence_type: crate::basic::value::SequenceType,
        keep_quality: bool,
        file_buffer_size: usize,
    ) -> Self {
        Self::spawn_reader(path, sequence_type, keep_quality, file_buffer_size, None)
    }
    /// Eagerly open the source when Source opens its input before output; the
    /// constructor detects compression and format before opening output. Record
    /// decoding still waits for demand.
    pub(crate) fn open_ready(
        path: std::path::PathBuf,
        sequence_type: crate::basic::value::SequenceType,
        keep_quality: bool,
        file_buffer_size: usize,
    ) -> io::Result<Self> {
        let prepared = PreparedQueryInput(std::sync::Arc::new(std::sync::Mutex::new(
            crate::data::fasta::prepare_source_sequence_reader(&path, file_buffer_size)?,
        )));
        Ok(Self::spawn_reader(
            path,
            sequence_type,
            keep_quality,
            file_buffer_size,
            Some(prepared),
        ))
    }
    fn spawn_reader(
        path: std::path::PathBuf,
        sequence_type: crate::basic::value::SequenceType,
        keep_quality: bool,
        file_buffer_size: usize,
        preopened: Option<PreparedQueryInput>,
    ) -> Self {
        let ready = preopened.is_some();
        let (sender, receiver) = std::sync::mpsc::sync_channel(1);
        let (request, requested) = std::sync::mpsc::sync_channel::<()>(0);
        let rewind_file = std::sync::Arc::new(std::sync::Mutex::new(None));
        let worker_rewind_file = rewind_file.clone();
        let input_path = path.clone();
        let worker = std::thread::spawn(move || {
            // A stream is opened only when the consumer asks for its first
            // record. Recovery and preflight failures can then drop stdin
            // without blocking on a producer which is waiting for input.
            if requested.recv().is_err() {
                return;
            }
            let result = (|| -> io::Result<()> {
                let reader = match preopened {
                    Some(prepared) => prepared,
                    None => PreparedQueryInput(std::sync::Arc::new(std::sync::Mutex::new(
                        crate::data::fasta::prepare_source_sequence_reader(
                            &input_path,
                            file_buffer_size,
                        )?,
                    ))),
                };
                *worker_rewind_file
                    .lock()
                    .map_err(|_| io::Error::other("query source lock poisoned"))? =
                    Some(reader.clone());
                crate::data::fasta::for_each_prepared_sequence_reader_bytes(
                    reader,
                    sequence_type,
                    file_buffer_size,
                    |mut record| {
                        if !keep_quality {
                            record.quality = None;
                        }
                        sender.send(Ok(record)).map_err(|_| {
                            io::Error::new(io::ErrorKind::BrokenPipe, "sequence consumer closed")
                        })?;
                        // Do not read the next record while extension owns the
                        // current query block. Drop closes this request channel.
                        requested.recv().map_err(|_| {
                            io::Error::new(io::ErrorKind::BrokenPipe, "sequence consumer closed")
                        })
                    },
                )
            })();
            if let Err(error) = result {
                let _ = sender.send(Err(error));
            }
        });
        Self {
            receiver: Some(receiver),
            request: Some(request),
            worker: Some(worker),
            rewind_file,
            path,
            sequence_type,
            keep_quality,
            file_buffer_size,
            ready,
        }
    }
    /// Source constructs a FastaFile and detects its format before opening
    /// output. Recovery can retain the unopened request-driven stream.
    fn prepare_for_master(&mut self) -> io::Result<()> {
        if !self.ready {
            let replacement = Self::open_ready(
                self.path.clone(),
                self.sequence_type,
                self.keep_quality,
                self.file_buffer_size,
            )?;
            *self = replacement;
        }
        Ok(())
    }
    fn rewind_after_count(&mut self) -> io::Result<()> {
        if let Some(file) = self
            .rewind_file
            .lock()
            .map_err(|_| io::Error::other("query source lock poisoned"))?
            .as_mut()
        {
            file.0
                .lock()
                .map_err(|_| io::Error::other("query source lock poisoned"))?
                .rewind()
                .map_err(io::Error::other)?;
        }
        Ok(())
    }
    pub(crate) fn next(&mut self) -> io::Result<Option<FastaRecord>> {
        self.ready = true;
        // A failed request means the producer has already reported EOF/error;
        // drain that result and join it through the same completion boundary.
        let _ = self.request.as_ref().expect("live record stream").send(());
        match self.receiver.as_ref().expect("live record stream").recv() {
            Ok(record) => record.map(Some),
            Err(_) => {
                if let Some(worker) = self.worker.take() {
                    worker
                        .join()
                        .map_err(|_| io::Error::other("sequence parser panicked"))?;
                }
                Ok(None)
            }
        }
    }
}
impl Drop for RecordStream {
    fn drop(&mut self) {
        drop(self.receiver.take());
        drop(self.request.take());
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

struct TemporaryDatabase(std::path::PathBuf, bool);
impl Drop for TemporaryDatabase {
    fn drop(&mut self) {
        if !self.1 {
            let _ = std::fs::remove_file(source_c_path(&self.0).as_ref());
        }
    }
}

struct IntermediateFile(std::path::PathBuf, usize, bool);
impl IntermediateFile {
    fn new(config: &BlastpConfig) -> io::Result<(Self, std::io::BufWriter<std::fs::File>)> {
        let directory = crate::util::io::temp_file::source_temp_directory_path(
            &config.tmpdir,
            config.output.as_deref().unwrap_or(std::path::Path::new("")),
        );
        let reserved = crate::util::io::TempFileData::init_in(false, directory, false)
            .map_err(io::Error::other)?;
        let path = std::path::PathBuf::from(&reserved.name);
        drop(reserved);
        let owner = Self(
            path.clone(),
            config.file_buffer_size,
            config.tuning.no_unlink,
        );
        let file = std::fs::OpenOptions::new()
            .write(true)
            .truncate(true)
            .create(cfg!(windows))
            .open(source_c_path(&path).as_ref())
            .map_err(|error| crate::util::io::FileOpenError::into_io(&path, error))?;
        Ok((
            owner,
            std::io::BufWriter::with_capacity(config.file_buffer_size, file),
        ))
    }
    fn reader(&self) -> io::Result<std::io::BufReader<std::fs::File>> {
        Ok(std::io::BufReader::with_capacity(
            self.1,
            std::fs::File::open(source_c_path(&self.0).as_ref())?,
        ))
    }
}
impl Drop for IntermediateFile {
    fn drop(&mut self) {
        if !self.2 {
            let _ = std::fs::remove_file(source_c_path(&self.0).as_ref());
        }
    }
}

struct ReferenceChunk {
    records: Vec<FastaRecord>,
    // Preserve source OIDs even when a PAL filter removes records between
    // matched subjects or a block spans multiple BLAST volumes.
    oids: Vec<u64>,
}

/// Source tracks returned extension targets before any cross-block output
/// culling. Accession rows are ordered by OID and may contain duplicates.
struct UnalignedTargets {
    path: std::path::PathBuf,
    file_buffer_size: usize,
    threads: usize,
    accession_spool: Option<std::sync::Arc<crate::data::accession_spool::AccessionSpool>>,
    aligned: std::sync::Mutex<crate::util::data_structures::BitVector>,
}
impl UnalignedTargets {
    fn new(config: &BlastpConfig, references: &ReferenceSource) -> io::Result<Option<Self>> {
        if config.tuning.unaligned_targets.as_os_str().is_empty() {
            return Ok(None);
        }
        let count = usize::try_from(references.totals().0).map_err(io::Error::other)?;
        Ok(Self::with_count(config, count).map(|mut tracker| {
            tracker.accession_spool = references.accession_spool();
            tracker
        }))
    }
    fn with_count(config: &BlastpConfig, count: usize) -> Option<Self> {
        (!config.tuning.unaligned_targets.as_os_str().is_empty()).then(|| Self {
            path: config.tuning.unaligned_targets.clone(),
            file_buffer_size: config.file_buffer_size,
            threads: if config.threads > 0 {
                config.threads as usize
            } else {
                crate::util::system::hardware_concurrency().max(1) as usize
            },
            accession_spool: None,
            aligned: std::sync::Mutex::new(crate::util::data_structures::BitVector::with_size(
                count,
            )),
        })
    }
    fn from_records(config: &BlastpConfig, records: &[FastaRecord]) -> Option<Self> {
        Self::with_count(config, records.len())
    }
    fn mark_capture(&self, capture: &NativeSearchCapture, oids: &[u64]) -> io::Result<()> {
        let mut aligned = self
            .aligned
            .lock()
            .map_err(|_| io::Error::other("Target tracking lock poisoned."))?;
        for matches in &capture.matches {
            for target in matches {
                let sorted = target.target_block_id as usize;
                let local = capture
                    .target_oids
                    .get(sorted)
                    .copied()
                    .unwrap_or(sorted as u64) as usize;
                let oid = *oids
                    .get(local)
                    .ok_or_else(|| io::Error::other("Target tracking OID outside block."))?;
                if oid >= aligned.size() {
                    return Err(io::Error::other("Target tracking OID outside database."));
                }
                aligned.set(oid as usize);
            }
        }
        Ok(())
    }

    /// Provided original owners already have real OID-to-accession metadata;
    /// stream it to the same disk format without retaining another ID vector.
    fn from_accessions(
        config: &BlastpConfig,
        count: usize,
        mut accession: impl FnMut(u64) -> io::Result<Vec<u8>>,
    ) -> io::Result<Option<Self>> {
        let Some(mut tracker) = Self::with_count(config, count) else {
            return Ok(None);
        };
        let spool = Self::new_spool(config)?;
        let errors = spool.io_errors();
        errors.capture_io((|| {
            for oid in 0..count {
                spool.push_bytes(&accession(oid as u64)?)?;
            }
            Ok(())
        })())?;
        tracker.accession_spool = Some(spool);
        Ok(Some(tracker))
    }
    fn new_spool(
        config: &BlastpConfig,
    ) -> io::Result<std::sync::Arc<crate::data::accession_spool::AccessionSpool>> {
        Self::new_spool_owned(
            config,
            crate::util::io::IoErrorOwner::default().with_source_destructors(),
        )
    }
    fn new_spool_owned(
        config: &BlastpConfig,
        errors: crate::util::io::IoErrorOwner,
    ) -> io::Result<std::sync::Arc<crate::data::accession_spool::AccessionSpool>> {
        let directory = crate::util::io::temp_file::source_temp_directory_path(
            &config.tmpdir,
            config.output.as_deref().unwrap_or(std::path::Path::new("")),
        );
        Ok(std::sync::Arc::new(
            crate::data::accession_spool::AccessionSpool::new_owned(
                &directory,
                config.tuning.no_unlink,
                config.file_buffer_size,
                errors,
            )?,
        ))
    }
    /// Source creates output, then reads/parses the entire real spool before
    /// emitting any row, even if every aligned bit is set.
    fn finish_spooled(self) -> io::Result<()> {
        let retained = self.accession_spool.clone();
        let errors = retained
            .as_ref()
            .map_or_else(Default::default, |spool| spool.io_errors());
        let result = (|| {
            use std::io::Write;
            let aligned = self
                .aligned
                .into_inner()
                .map_err(|_| io::Error::other("Target tracking lock poisoned."))?;
            let mut output = io::BufWriter::with_capacity(
                self.file_buffer_size,
                std::fs::File::create(source_c_path(&(self.path)).as_ref())?,
            );
            let spool = self.accession_spool.ok_or_else(|| io::Error::other(
            "Source provided database accession lookup dereferences a missing OID_TO_ACC_MAPPING spool."))?;
            let table = spool.read_table(self.threads)?;
            for oid in 0..aligned.size() {
                if !aligned.get(oid as usize) {
                    if oid >= table.size() as u64 {
                        return Err(io::Error::other(
                            "Source accession table lookup is outside its actual mapped OID rows.",
                        ));
                    }
                    output.write_all(table.record(oid as i64).string_bytes(0))?;
                    output.write_all(b"\n")?;
                }
            }
            output.flush()
        })();
        errors.capture_io(result)
    }
}

enum ReferenceSource {
    Fasta(fasta_reference::FastaReference),
    Memory {
        records: Vec<FastaRecord>,
        cursor: usize,
    },
    Blast {
        file: Box<crate::data::blastdb::blastdb::BlastDB>,
        load_threads: i32,
        threads: i32,
        minichunk: usize,
        cursor: u64,
        taxonomy_flags: crate::data::sequence_file::SequenceFileFlags,
        filter: Option<crate::data::sequence_file::DbFilter>,
    },
    Database {
        file: Box<crate::data::dmnd::DatabaseFile>,
        cursor: u64,
        _temporary: Option<TemporaryDatabase>,
        taxonomy_flags: crate::data::sequence_file::SequenceFileFlags,
        filter: Option<crate::data::sequence_file::DbFilter>,
    },
}
impl ReferenceSource {
    fn checks_protein_input(&self) -> bool {
        matches!(self, Self::Fasta(_))
    }
    fn accession_spool(
        &self,
    ) -> Option<std::sync::Arc<crate::data::accession_spool::AccessionSpool>> {
        match self {
            Self::Fasta(file) => file.accession_spool(),
            Self::Database { file, .. } => file.accession_spool.clone(),
            Self::Blast { file, .. } => file.accession_spool.clone(),
            Self::Memory { .. } => None,
        }
    }
    /// Resolve one physical OID from its real title owner after output closes.
    /// Filters select search records, not the full unaligned-target domain.
    fn accession(&mut self, oid: u64) -> io::Result<Vec<u8>> {
        let title = match self {
            Self::Fasta(file) => file.title(oid)?,
            Self::Memory { records, .. } => {
                let record = records
                    .get(usize::try_from(oid).map_err(io::Error::other)?)
                    .ok_or_else(|| io::Error::other("Target tracking OID outside database."))?;
                return Ok(crate::util::sequence::seqid_bytes(&record.id));
            }
            Self::Blast { file, .. } => {
                file.source_io(|file| file.title_only_bytes(oid, true, true))?
            }
            Self::Database { file, .. } => file.title_bytes_io(oid)?,
        };
        Ok(crate::util::sequence::seqid_bytes(&title))
    }
    fn record(&mut self, oid: u64) -> io::Result<FastaRecord> {
        match self {
            Self::Fasta(file) => file.record(oid),
            Self::Memory { records, .. } => records
                .get(oid as usize)
                .cloned()
                .ok_or_else(|| io::Error::other("Ranked target OID is outside the database")),
            Self::Blast { file, .. } => {
                let (sequence, id) = file.record_at_bytes_io(oid)?;
                Ok(FastaRecord {
                    id,
                    sequence,
                    quality: None,
                })
            }
            Self::Database { file, .. } => {
                let (sequence, id) = file.record_bytes_io(oid)?;
                Ok(FastaRecord {
                    id,
                    sequence,
                    quality: None,
                })
            }
        }
    }

    fn memory(records: Vec<FastaRecord>) -> Self {
        Self::Memory { records, cursor: 0 }
    }
    fn open(config: &BlastpConfig) -> io::Result<Self> {
        Self::open_with_sequence_type(config, crate::basic::value::SequenceType::AminoAcid)
    }
    fn open_with_sequence_type(
        config: &BlastpConfig,
        sequence_type: crate::basic::value::SequenceType,
    ) -> io::Result<Self> {
        let errors = crate::util::io::IoErrorOwner::default().with_source_destructors();
        let mut retained_spool = None;
        let result = (|| {
            let (kind, path) = crate::data::sequence_file::detect_type_source_with_buffer(
                &[std::path::PathBuf::from(&config.database)],
                crate::data::sequence_file::SequenceFileFlags::NONE,
                crate::data::sequence_file::SequenceFileConfig::default(),
                config.file_buffer_size,
            )?;
            if kind == crate::data::sequence_file::SequenceFileType::Blast {
                if sequence_type != crate::basic::value::SequenceType::AminoAcid {
                    return Err(io::Error::other("Nucleotide search requires a nucleotide database; the native BLAST reader decodes protein databases."));
                }
                use crate::data::sequence_file::SequenceFileFlags as Flags;
                if config.target_indexed {
                    return Err(io::Error::other(
                        "--target-indexed is not compatible with BLAST databases.",
                    ));
                }
                let taxonomy_flags = super::blastp::requested_search_taxonomy_flags(config);
                use crate::output::format::{FormatCode, OutputFlags, OutputFormatSpec};
                let output_flags = match config
                    .outfmt
                    .first()
                    .and_then(|name| FormatCode::parse_bytes(name))
                {
                    Some(FormatCode::Daa) => {
                        OutputFormatSpec::daa_format(
                            config.archive.all_subject_titles,
                            config.archive.all_subject_ids,
                        )
                        .flags
                    }
                    Some(FormatCode::Pairwise) => OutputFormatSpec::pairwise_format().flags,
                    Some(FormatCode::Xml) => OutputFormatSpec::xml_format().flags,
                    Some(FormatCode::Sam) => OutputFormatSpec::sam_format().flags,
                    Some(FormatCode::Tabular | FormatCode::JsonFlat) | None => {
                        super::search_output::tabular_format(config)?.format.flags
                    }
                    _ => OutputFlags::NONE,
                };
                let mut flags = Flags::SEQS | Flags::TITLES | taxonomy_flags;
                if output_flags.any(OutputFlags::ALL_SEQIDS) {
                    flags |= Flags::ALL_SEQIDS;
                }
                if output_flags.any(OutputFlags::FULL_TITLES) || config.no_self_hits {
                    flags |= Flags::FULL_TITLES;
                }
                let accession_spool = if config.tuning.unaligned_targets.as_os_str().is_empty() {
                    None
                } else {
                    Some({
                        let spool = UnalignedTargets::new_spool_owned(config, errors.clone())?;
                        retained_spool = Some(spool.clone());
                        spool
                    })
                };
                let mut file = crate::data::blastdb::blastdb::BlastDB::new_path_with_io_owner(
                    &path, flags, &errors,
                )?;
                // Source BlastDB never calls add_seqid_mapping. Retain the real
                // empty base spool instead of inventing an accession map.
                file.accession_spool = accession_spool;
                let filter = file.search_filter_bytes_path_io(
                    &config.taxonlist,
                    &config.taxon_exclude,
                    &config.seqidlist,
                    config.skip_missing_seqids,
                    config.file_buffer_size,
                )?;
                return Ok(Self::Blast {
                    file: Box::new(file),
                    load_threads: config.tuning.load_threads,
                    threads: if config.threads > 0 {
                        config.threads
                    } else {
                        crate::util::system::hardware_concurrency() as i32
                    },
                    minichunk: config.tuning.minichunk,
                    cursor: 0,
                    taxonomy_flags,
                    filter,
                });
            }
            if kind != crate::data::sequence_file::SequenceFileType::Dmnd {
                use crate::data::sequence_file::SequenceFileFlags as Flags;
                let taxonomy_flags = super::blastp::requested_search_taxonomy_flags(config);
                if taxonomy_flags.contains(
                    Flags::TAXON_MAPPING
                        | Flags::TAXON_NODES
                        | Flags::TAXON_RANKS
                        | Flags::TAXON_SCIENTIFIC_NAMES,
                ) {
                    return Err(io::Error::other(
                        "Fasta database format does not support taxonomic features.",
                    ));
                }
                let source = fasta_reference::FastaReference::open(
                    config,
                    path,
                    sequence_type,
                    errors.clone(),
                )?;
                retained_spool = source.accession_spool();
                if !config.seqidlist.as_os_str().is_empty() {
                    return Err(io::Error::other(
                        "The FASTA database format does not support filtering by accession.",
                    ));
                }
                return Ok(Self::Fasta(source));
            }
            let taxonomy_flags = super::blastp::requested_search_taxonomy_flags(config);
            if !config.taxonlist.is_empty() && !config.taxon_exclude.is_empty() {
                return Err(io::Error::other(
                    "Options --taxonlist and --taxon-exclude are mutually exclusive.",
                ));
            }
            let accession_spool = if config.tuning.unaligned_targets.as_os_str().is_empty() {
                None
            } else {
                Some({
                    let spool = UnalignedTargets::new_spool_owned(config, errors.clone())?;
                    retained_spool = Some(spool.clone());
                    spool
                })
            };
            let mut file = super::blastp::open_search_dmnd_with_buffer(
                &path,
                taxonomy_flags,
                Some(config.file_buffer_size),
            )?;
            if let Some(spool) = &accession_spool {
                file.write_accession_spool(spool)?;
            }
            file.accession_spool = accession_spool;
            if file.header2.sequence_type().map_err(io::Error::other)? != sequence_type {
                return Err(io::Error::other(
                    "Database sequence type does not match the search mode.",
                ));
            }
            let filter = if config.taxonlist.is_empty() && config.taxon_exclude.is_empty() {
                None
            } else {
                struct Access<'a>(std::cell::RefCell<&'a mut crate::data::dmnd::DatabaseFile>);
                impl crate::data::sequence_file::TaxonomyAccess for Access<'_> {
                    fn parent(&self, taxid: i32) -> i32 {
                        self.0.borrow().get_parent(taxid).unwrap_or(0)
                    }
                    fn rank(&self, taxid: i32) -> i32 {
                        self.0.borrow().rank(taxid).unwrap_or(0)
                    }
                    fn max_taxid(&self) -> i32 {
                        self.0.borrow().max_taxid().unwrap_or(0)
                    }
                }
                let count = file.sequence_count();
                let access = Access(std::cell::RefCell::new(&mut file));
                Some(
                    crate::data::sequence_file::filter_by_taxonomy_bytes(
                        &access,
                        if config.taxonlist.is_empty() {
                            &config.taxon_exclude
                        } else {
                            &config.taxonlist
                        },
                        b',',
                        !config.taxon_exclude.is_empty(),
                        count,
                        |oid| access.0.borrow().taxids(oid as usize),
                        |oid| access.0.borrow_mut().seq_length(oid as usize),
                    )
                    .map_err(io::Error::other)?,
                )
            };
            if !config.seqidlist.as_os_str().is_empty() {
                return Err(io::Error::other(if filter.is_some() {
                    "--seqidlist is not compatible with taxonomy filtering."
                } else {
                    "The .dmnd database format does not support filtering by accession."
                }));
            }
            Ok(Self::Database {
                file: Box::new(file),
                cursor: 0,
                _temporary: None,
                taxonomy_flags,
                filter,
            })
        })();
        let result = errors.capture_io(result);
        drop(retained_spool);
        result
    }

    fn source_errors(&self) -> crate::util::io::IoErrorOwner {
        self.accession_spool()
            .map_or_else(Default::default, |spool| spool.io_errors())
    }
    fn species_rank_ids(&mut self, oids: &[u64]) -> io::Result<Vec<Vec<u32>>> {
        let mut output = Vec::with_capacity(oids.len());
        for &oid in oids {
            let taxids = match self {
                Self::Blast { file, .. } => file.taxids(oid as usize),
                Self::Database { file, .. } => {
                    file.taxids(oid as usize).map_err(io::Error::other)?
                }
                Self::Fasta(_) | Self::Memory { .. } => {
                    return Err(io::Error::other(
                        "Taxonomy is not available for in-memory search.",
                    ))
                }
            };
            let mut species = BTreeSet::new();
            for mut taxid in taxids {
                let mut resolved = false;
                for _ in 0..=64 {
                    let rank = match self {
                        Self::Blast { file, .. } => file.rank(taxid),
                        Self::Database { file, .. } => {
                            file.rank(taxid).map_err(io::Error::other)?
                        }
                        Self::Fasta(_) | Self::Memory { .. } => unreachable!(),
                    };
                    if rank == crate::data::taxonomy::Rank::Species as i32 {
                        species.insert(taxid as u32);
                        resolved = true;
                        break;
                    }
                    if taxid <= 1 {
                        species.insert(0);
                        resolved = true;
                        break;
                    }
                    taxid = match self {
                        Self::Blast { file, .. } => file.get_parent(taxid),
                        Self::Database { file, .. } => {
                            file.get_parent(taxid).map_err(io::Error::other)?
                        }
                        Self::Fasta(_) | Self::Memory { .. } => unreachable!(),
                    };
                }
                if !resolved {
                    return Err(io::Error::other("Path in taxonomy too long (rank_taxid)."));
                }
            }
            output.push(species.into_iter().collect());
        }
        Ok(output)
    }

    fn load_matched_taxonomy(&mut self, dictionary: &mut ReferenceDictionary) -> io::Result<()> {
        use crate::data::sequence_file::SequenceFileFlags as Flags;
        if let Self::Blast {
            file,
            taxonomy_flags,
            ..
        } = self
        {
            if !taxonomy_flags.contains(Flags::TAXON_MAPPING) {
                return Ok(());
            }
            let needs_tree = taxonomy_flags.contains(Flags::TAXON_NODES)
                || taxonomy_flags.contains(Flags::TAXON_RANKS)
                || taxonomy_flags.contains(Flags::TAXON_SCIENTIFIC_NAMES);
            if needs_tree && dictionary.taxonomy.is_none() {
                dictionary.taxonomy = Some(crate::data::taxonomy::TaxonomyTree::new());
            }
            for entry in &mut dictionary.entries[dictionary.taxids_loaded..] {
                entry.taxids = file.taxids(entry.oid as usize);
                if let Some(tree) = dictionary.taxonomy.as_mut() {
                    for &assigned in &entry.taxids {
                        let mut taxid = assigned;
                        while taxid > 0 && dictionary.seen_taxids.insert(taxid) {
                            let parent = if taxonomy_flags.contains(Flags::TAXON_NODES) {
                                file.get_parent(taxid)
                            } else {
                                0
                            };
                            let rank = if taxonomy_flags.contains(Flags::TAXON_RANKS) {
                                usize::try_from(file.rank(taxid))
                                    .ok()
                                    .and_then(crate::data::taxonomy::Rank::from_index)
                                    .map_or_else(String::new, |rank| rank.name().to_owned())
                            } else {
                                String::new()
                            };
                            let name = if taxonomy_flags.contains(Flags::TAXON_SCIENTIFIC_NAMES) {
                                file.taxon_scientific_name_bytes(taxid)
                            } else {
                                taxid.to_string().into_bytes()
                            };
                            tree.add_node_bytes(crate::data::taxonomy::TaxonomyNode {
                                taxid,
                                parent,
                                rank,
                                name,
                            });
                            if !taxonomy_flags.contains(Flags::TAXON_NODES)
                                || taxid == 1
                                || parent == taxid
                            {
                                break;
                            }
                            taxid = parent;
                        }
                    }
                }
            }
            dictionary.taxids_loaded = dictionary.entries.len();
            return Ok(());
        }
        let Self::Database {
            file,
            taxonomy_flags,
            ..
        } = self
        else {
            return Ok(());
        };
        if !taxonomy_flags.contains(Flags::TAXON_MAPPING) {
            return Ok(());
        }
        let needs_tree = taxonomy_flags.contains(Flags::TAXON_NODES)
            || taxonomy_flags.contains(Flags::TAXON_RANKS)
            || taxonomy_flags.contains(Flags::TAXON_SCIENTIFIC_NAMES);
        if needs_tree && dictionary.taxonomy.is_none() {
            dictionary.taxonomy = Some(crate::data::taxonomy::TaxonomyTree::new());
        }
        for entry in &mut dictionary.entries[dictionary.taxids_loaded..] {
            entry.taxids = file.taxids(entry.oid as usize).map_err(io::Error::other)?;
            if let Some(tree) = dictionary.taxonomy.as_mut() {
                super::blastp::extend_dmnd_taxonomy(
                    file,
                    *taxonomy_flags,
                    &entry.taxids,
                    tree,
                    &mut dictionary.seen_taxids,
                )?;
            }
        }
        dictionary.taxids_loaded = dictionary.entries.len();
        Ok(())
    }
    fn totals(&self) -> (u64, u64) {
        match self {
            Self::Fasta(file) => file.totals(),
            Self::Memory { records, .. } => (
                records.len() as u64,
                records.iter().map(|r| r.sequence.len() as u64).sum(),
            ),
            Self::Database { file, .. } => (file.sequence_count(), file.letters()),
            Self::Blast { file, .. } => (file.sequence_count(), file.letters()),
        }
    }
    // Source run_query_iteration compares its first block with the active
    // filter's one_count, retaining full database totals for output metadata.
    fn selected_sequence_count(&self) -> u64 {
        match self {
            Self::Blast {
                filter: Some(filter),
                ..
            }
            | Self::Database {
                filter: Some(filter),
                ..
            } => filter
                .oid_filter
                .iter()
                .filter(|&&selected| selected)
                .count() as u64,
            _ => self.totals().0,
        }
    }

    fn scoring_letters(&self) -> u64 {
        match self {
            Self::Blast {
                filter: Some(filter),
                ..
            }
            | Self::Database {
                filter: Some(filter),
                ..
            } if filter.letter_count != 0 => filter.letter_count,
            _ => self.totals().1,
        }
    }
    fn seek_oid(&mut self, oid: u64) {
        match self {
            Self::Fasta(file) => file.cursor = oid,
            Self::Memory { cursor, .. } => *cursor = oid as usize,
            Self::Database { cursor, .. } | Self::Blast { cursor, .. } => *cursor = oid,
        }
    }
    fn query_owner(&self, config: &BlastpConfig) -> io::Result<Self> {
        match self {
            Self::Fasta(file) => Ok(Self::Fasta(file.duplicate()?)),
            Self::Memory { records, .. } => Ok(Self::memory(records.clone())),
            Self::Database {
                file,
                taxonomy_flags,
                filter,
                ..
            } => Ok(Self::Database {
                file: Box::new(super::blastp::open_search_dmnd_with_buffer(
                    file.file_name(),
                    *taxonomy_flags,
                    Some(config.file_buffer_size),
                )?),
                cursor: 0,
                // The reference owner retains any temporary database until the
                // independently positioned query handle has finished reading.
                _temporary: None,
                taxonomy_flags: *taxonomy_flags,
                filter: filter.clone(),
            }),
            Self::Blast { .. } => Self::open(config),
        }
    }
    fn rewind(&mut self) {
        match self {
            Self::Fasta(file) => file.cursor = 0,
            Self::Memory { cursor, .. } => *cursor = 0,
            Self::Database { cursor, .. } | Self::Blast { cursor, .. } => *cursor = 0,
        }
    }
    fn chunk(&mut self, limit: u64) -> io::Result<ReferenceChunk> {
        if let Self::Blast {
            file,
            cursor,
            filter,
            load_threads,
            threads,
            minichunk,
            ..
        } = self
        {
            if *cursor == 0 {
                file.source_io(|file| file.set_seqinfo_ptr(0))?;
            }
            let (decoded, processed) = file.source_io(|file| {
                crate::data::blastdb::parallel_load::load(
                    file,
                    limit,
                    *load_threads,
                    *threads,
                    *minichunk,
                    filter.as_ref(),
                )
            })?;
            *cursor += processed as u64;
            let mut records = Vec::with_capacity(decoded.oids.len());
            for index in 0..decoded.oids.len() {
                records.push(FastaRecord {
                    id: decoded.ids.get(index).to_vec(),
                    sequence: decoded.seqs.get(index).to_vec(),
                    quality: None,
                });
            }
            return Ok(ReferenceChunk {
                records,
                oids: decoded.oids,
            });
        }
        let mut output = ReferenceChunk {
            records: Vec::new(),
            oids: Vec::new(),
        };
        let mut letters = 0;
        while output.records.is_empty() || letters < limit {
            let record = match self {
                Self::Fasta(file) => {
                    if file.cursor >= file.totals().0 {
                        None
                    } else {
                        let oid = file.cursor;
                        let record = file.record(oid)?;
                        file.cursor += 1;
                        output.oids.push(oid);
                        Some(record)
                    }
                }
                Self::Memory { records, cursor } => {
                    let record = records.get(*cursor).cloned();
                    if record.is_some() {
                        output.oids.push(*cursor as u64);
                        *cursor += 1;
                    }
                    record
                }
                Self::Database {
                    file,
                    cursor,
                    filter,
                    ..
                } => {
                    while *cursor < file.sequence_count()
                        && filter.as_ref().is_some_and(|filter| !filter.get(*cursor))
                    {
                        *cursor += 1;
                    }
                    if *cursor == file.sequence_count() {
                        None
                    } else {
                        let (sequence, id) = file.record_bytes_io(*cursor)?;
                        output.oids.push(*cursor);
                        *cursor += 1;
                        Some(FastaRecord {
                            quality: None,
                            id,
                            sequence,
                        })
                    }
                }
                Self::Blast {
                    file,
                    cursor,
                    filter,
                    ..
                } => {
                    while *cursor < file.sequence_count()
                        && filter.as_ref().is_some_and(|filter| !filter.get(*cursor))
                    {
                        *cursor += 1;
                    }
                    if *cursor == file.sequence_count() {
                        None
                    } else {
                        let (sequence, id) = file.record_at_bytes_io(*cursor)?;
                        output.oids.push(*cursor);
                        *cursor += 1;
                        Some(FastaRecord {
                            quality: None,
                            id,
                            sequence,
                        })
                    }
                }
            };
            let Some(record) = record else {
                break;
            };
            letters += record.sequence.len() as u64;
            output.records.push(record);
        }
        Ok(output)
    }
}

enum ProteinQueries {
    Memory(std::vec::IntoIter<FastaRecord>),
    Stream(RecordStream),
    Database(Box<ReferenceSource>),
}
impl ProteinQueries {
    fn chunk(&mut self, limit: u64) -> io::Result<Vec<FastaRecord>> {
        if let Self::Database(database) = self {
            return Ok(database.chunk(limit)?.records);
        }
        let mut records = Vec::new();
        let mut letters = 0;
        while records.is_empty() || letters < limit {
            let record = match self {
                Self::Memory(records) => records.next(),
                Self::Stream(stream) => stream.next()?,
                Self::Database(_) => unreachable!(),
            };
            let Some(record) = record else {
                break;
            };
            if record.sequence.is_empty() {
                continue;
            }
            letters += record.sequence.len() as u64;
            records.push(record);
        }
        Ok(records)
    }
}

struct QueryChunk {
    protein_oids: Vec<u64>,
    self_aln_scores: Vec<f64>,
    records: Vec<FastaRecord>,
    translated: Option<Block>,
    mates: Vec<Option<usize>>,
}
struct QueryView<'a> {
    oid: Option<u64>,
    title: &'a [u8],
    source: &'a [Letter],
    frames: Vec<&'a [Letter]>,
    mate: Option<&'a [Letter]>,
    quality: &'a [u8],
    self_aln_score: Option<f64>,
}
impl QueryChunk {
    fn count(&self) -> usize {
        self.translated
            .as_ref()
            .map_or(self.records.len(), |block| block.source_seqs().len())
    }
    fn query(&self, id: usize) -> io::Result<QueryView<'_>> {
        if let Some(block) = &self.translated {
            Ok(QueryView {
                oid: Some(block.block_id2oid(id as u32)),
                title: block.ids().map_err(io::Error::other)?.get(id),
                source: block.source_seqs().get(id),
                self_aln_score: self.self_aln_scores.get(id).copied(),
                quality: if block.qual().empty() {
                    b""
                } else {
                    block.qual().get(id)
                },
                frames: (id * 6..id * 6 + 6)
                    .map(|frame| block.seqs().get(frame))
                    .collect(),
                mate: self
                    .mates
                    .get(id)
                    .copied()
                    .flatten()
                    .map(|mate| block.source_seqs().get(mate)),
            })
        } else {
            let query = &self.records[id];
            Ok(QueryView {
                oid: self.protein_oids.get(id).copied(),
                title: &query.id,
                source: &query.sequence,
                self_aln_score: self.self_aln_scores.get(id).copied(),
                quality: query.quality.as_deref().unwrap_or(b""),
                frames: vec![&query.sequence],
                mate: None,
            })
        }
    }
}

enum QueryInput<'a> {
    Protein(ProteinQueries),
    ProvidedTranslated {
        block: Block,
        cursor: usize,
        mates: Vec<Option<usize>>,
    },
    Translated {
        streams: Vec<RecordStream>,
        next: usize,
        next_oid: u64,
        config: &'a super::blastx::BlastxConfig,
        frame_mask: i32,
    },
}
impl QueryInput<'_> {
    fn prepare_for_master(&mut self) -> io::Result<()> {
        match self {
            Self::Protein(ProteinQueries::Stream(stream)) => stream.prepare_for_master(),
            Self::Translated { streams, .. } => {
                for stream in streams {
                    stream.prepare_for_master()?;
                }
                Ok(())
            }
            // Provided Blocks and database owners are already open.
            _ => Ok(()),
        }
    }
    fn checks_protein_input(&self) -> bool {
        match self {
            Self::Protein(ProteinQueries::Stream(_)) => true,
            Self::Protein(ProteinQueries::Database(database)) => database.checks_protein_input(),
            _ => false,
        }
    }
    fn rewind_after_count(&mut self) -> io::Result<()> {
        match self {
            Self::Protein(ProteinQueries::Stream(stream)) => stream.rewind_after_count(),
            Self::Translated { streams, .. } => {
                for stream in streams {
                    stream.rewind_after_count()?;
                }
                Ok(())
            }
            Self::ProvidedTranslated { cursor, .. } => {
                *cursor = 0;
                Ok(())
            }
            Self::Protein(ProteinQueries::Memory(_)) => Ok(()),
            Self::Protein(ProteinQueries::Database(database)) => {
                database.rewind();
                Ok(())
            }
        }
    }
    fn chunk(&mut self, budget: u64) -> io::Result<QueryChunk> {
        match self {
            Self::Protein(ProteinQueries::Database(database)) => {
                let chunk = database.chunk(budget)?;
                Ok(QueryChunk {
                    protein_oids: chunk.oids,
                    self_aln_scores: Vec::new(),
                    records: chunk.records,
                    translated: None,
                    mates: Vec::new(),
                })
            }
            Self::Protein(input) => Ok(QueryChunk {
                protein_oids: Vec::new(),
                self_aln_scores: Vec::new(),
                records: input.chunk(budget)?,
                translated: None,
                mates: Vec::new(),
            }),
            Self::ProvidedTranslated {
                block,
                cursor,
                mates,
            } => {
                let begin = *cursor;
                let mut end = begin;
                let mut letters = 0;
                let mut required_end = begin;
                while end < block.source_seqs().len()
                    && (end == begin || letters < budget || end < required_end)
                {
                    letters += block
                        .translated_load_letters(end)
                        .map_err(io::Error::other)?;
                    if let Some(mate) = mates[end] {
                        if mate < begin {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidInput,
                                "Translated mate relation crosses a completed query block.",
                            ));
                        }
                        required_end = required_end.max(mate + 1);
                    }
                    end += 1;
                }
                *cursor = end;
                let output = block
                    .translated_source_range(begin..end)
                    .map_err(io::Error::other)?;
                let local_mates = mates[begin..end]
                    .iter()
                    .map(|mate| mate.map(|id| id - begin))
                    .collect();
                Ok(QueryChunk {
                    protein_oids: Vec::new(),
                    self_aln_scores: Vec::new(),
                    records: Vec::new(),
                    translated: Some(output),
                    mates: local_mates,
                })
            }
            Self::Translated {
                streams,
                next,
                next_oid,
                config,
                frame_mask,
            } => {
                let mut block = Block::new();
                block.set_max_motif_len(config.tuning.max_motif_len);
                let mut letters = 0;
                while block.source_seqs().is_empty()
                    || letters < budget
                    || block.source_seqs().len() % streams.len() != 0
                {
                    let Some(record) = streams[*next].next()? else {
                        // FastaFile stops on the current file's EOF, then checks
                        // that the round-robin cursor is synchronized.
                        if *next != 0 || (streams.len() == 2 && streams[1].next()?.is_some()) {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "Unequal number of sequences in paired read files.",
                            ));
                        }
                        break;
                    };
                    *next = (*next + 1) % streams.len();
                    if record.sequence.is_empty() {
                        continue;
                    }
                    let full_title = config
                        .outfmt
                        .first()
                        .is_some_and(|name| matches!(name.as_slice(), b"0" | b"5" | b"xml"))
                        || config.outfmt.iter().any(|field| field == b"qtitle");
                    let short_title;
                    let title = if full_title {
                        record.id.as_slice()
                    } else {
                        short_title = crate::util::sequence::seqid_bytes(&record.id);
                        &short_title
                    };
                    letters += if let Some(translator) = config.tuning.runtime.translator() {
                        block.push_back_translated_with_translator_bytes(
                            &record.sequence,
                            Some(title),
                            *next_oid,
                            *frame_mask,
                            config.min_orf.unwrap_or(0) as i32,
                            translator,
                        )
                    } else {
                        block.push_back_translated_with_genetic_code_bytes(
                            &record.sequence,
                            Some(title),
                            *next_oid,
                            *frame_mask,
                            config.min_orf.unwrap_or(0) as i32,
                            config.query_gencode,
                        )
                    }
                    .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?
                        as u64;
                    if config
                        .outfmt
                        .iter()
                        .any(|field| matches!(field.as_slice(), b"qqual" | b"full_qqual"))
                        || config.auxiliary_output.needs_quality()
                        || config.tuning.store_query_quality
                    {
                        block.push_source_quality(record.quality.as_deref().unwrap_or(b""));
                    }
                    *next_oid += 1;
                }
                let count = block.source_seqs().len();
                let mates = (0..count)
                    .map(|id| (streams.len() == 2 && (id ^ 1) < count).then_some(id ^ 1))
                    .collect();
                Ok(QueryChunk {
                    protein_oids: Vec::new(),
                    self_aln_scores: Vec::new(),
                    records: Vec::new(),
                    translated: Some(block),
                    mates,
                })
            }
        }
    }
}

pub(crate) fn run_translated_files(
    config: &super::blastx::BlastxConfig,
    block_config: BlastpConfig,
    budget: u64,
) -> io::Result<()> {
    crate::basic::statistics::with_search_statistics_prepared(
        !block_config.tuning.multiprocessing.enabled
            || (!block_config.tuning.multiprocessing.initialize
                && !block_config.tuning.multiprocessing.recover),
        || super::source_config::begin_search_temp(&block_config),
        || {
            reject_translated_self_search(&block_config)?;
            if config.query_files.len() > 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "blastx requires one or two query files",
                ));
            }
            let frame_mask = match config.strand.as_slice() {
                b"both" => 0b111111,
                b"plus" => 0b000111,
                b"minus" => 0b111000,
                value => {
                    let mut message = b"invalid value for --strand: ".to_vec();
                    message.extend_from_slice(value);
                    return Err(crate::util::io::RawDiagnosticError(message).into_io());
                }
            };
            let references = ReferenceSource::open(&block_config)?;
            let query_files = if config.query_files.is_empty() {
                eprintln!(
                    "Query file parameter (--query/-q) is missing. Input will be read from stdin."
                );
                vec![std::path::PathBuf::new()]
            } else {
                config.query_files.clone()
            };
            let streams = query_files
                .iter()
                .map(|path| {
                    RecordStream::open(
                        path.into(),
                        crate::basic::value::SequenceType::Nucleotide,
                        config
                            .outfmt
                            .iter()
                            .any(|field| matches!(field.as_slice(), b"qqual" | b"full_qqual"))
                            || config.auxiliary_output.needs_quality()
                            || config.tuning.store_query_quality,
                        config.file_buffer_size,
                    )
                })
                .collect();
            run_protein(
                &block_config,
                references,
                QueryInput::Translated {
                    streams,
                    next: 0,
                    next_oid: 0,
                    config,
                    frame_mask,
                },
                budget,
            )
        },
    )
}

pub(crate) fn run_translated_block(
    config: &BlastpConfig,
    block: Block,
    budget: u64,
) -> io::Result<()> {
    crate::basic::statistics::with_search_statistics_prepared(
        !config.tuning.multiprocessing.enabled
            || (!config.tuning.multiprocessing.initialize
                && !config.tuning.multiprocessing.recover),
        || super::source_config::begin_search_temp(config),
        || {
            let mates = config
                .translated_query_layout
                .as_ref()
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Translated block requires a translated-query layout.",
                    )
                })?
                .mate_indices
                .clone();
            if mates.len() != block.source_seqs().len()
                || block.seqs().len() != mates.len().saturating_mul(6)
                || mates.iter().flatten().any(|&id| id >= mates.len())
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Invalid translated query block layout.",
                ));
            }
            run_protein(
                config,
                ReferenceSource::open(config)?,
                QueryInput::ProvidedTranslated {
                    block,
                    cursor: 0,
                    mates,
                },
                budget,
            )
        },
    )
}

pub(crate) fn run_protein_records(
    config: &BlastpConfig,
    database: Vec<FastaRecord>,
    queries: Vec<FastaRecord>,
    budget: u64,
) -> io::Result<()> {
    crate::basic::statistics::with_search_statistics_prepared(
        !config.tuning.multiprocessing.enabled
            || (!config.tuning.multiprocessing.initialize
                && !config.tuning.multiprocessing.recover),
        || super::source_config::begin_search_temp(config),
        || {
            run_protein(
                config,
                ReferenceSource::memory(database),
                QueryInput::Protein(ProteinQueries::Memory(queries.into_iter())),
                budget,
            )
        },
    )
}
pub(crate) fn run_protein_queries(
    config: &BlastpConfig,
    queries: Vec<FastaRecord>,
    budget: u64,
) -> io::Result<()> {
    crate::basic::statistics::with_search_statistics_prepared(
        !config.tuning.multiprocessing.enabled
            || (!config.tuning.multiprocessing.initialize
                && !config.tuning.multiprocessing.recover),
        || super::source_config::begin_search_temp(config),
        || {
            run_protein(
                config,
                ReferenceSource::open(config)?,
                QueryInput::Protein(ProteinQueries::Memory(queries.into_iter())),
                budget,
            )
        },
    )
}
pub(crate) fn run_protein_files(config: &BlastpConfig, budget: u64) -> io::Result<()> {
    crate::basic::statistics::with_search_statistics_prepared(
        !config.tuning.multiprocessing.enabled
            || (!config.tuning.multiprocessing.initialize
                && !config.tuning.multiprocessing.recover),
        || super::source_config::begin_search_temp(config),
        || {
            if !config.tuning.self_search && config.query_files.len() > 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Exactly one protein query file is required.",
                ));
            }
            let references = ReferenceSource::open(config)?;
            crate::data::accession_spool::source_scope(references.accession_spool(), || {
                if config.tuning.self_search {
                    let queries = QueryInput::Protein(ProteinQueries::Database(Box::new(
                        references.query_owner(config)?,
                    )));
                    return run_protein(config, references, queries, budget);
                }
                let query_file = config.query_files.first().cloned().unwrap_or_else(|| {
                    eprintln!(
                    "Query file parameter (--query/-q) is missing. Input will be read from stdin."
                );
                    std::path::PathBuf::new()
                });
                let queries = QueryInput::Protein(ProteinQueries::Stream(RecordStream::open(
                    query_file.into(),
                    crate::basic::value::SequenceType::AminoAcid,
                    config
                        .outfmt
                        .iter()
                        .any(|field| matches!(field.as_slice(), b"qqual" | b"full_qqual"))
                        || config.auxiliary_output.needs_quality()
                        || config.tuning.store_query_quality,
                    config.file_buffer_size,
                )));
                run_protein(config, references, queries, budget)
            })
        },
    )
}

/// Source master sorts before printing the first header and before masking.
/// keep_target_id observes the previous query's final reference counter here.
fn prepare_master_query_chunk(
    chunk: &mut QueryChunk,
    config: &BlastpConfig,
    previous_reference: u32,
    min_length_ratio: f64,
) -> io::Result<Vec<u64>> {
    let keep_ids = min_length_ratio != 0.0
        || config.global_ranking_targets != 0
        || (config.tuning.self_search && previous_reference == 0)
        || config.tuning.setup.lin_stage1_combo;
    let sort = min_length_ratio > 0.0
        || (config.tuning.setup.lin_stage1_query && !keep_ids && !config.tuning.kmer_ranking);
    let mut oids = (0..chunk.count() as u64).collect::<Vec<_>>();
    if !sort {
        return Ok(oids);
    }
    if chunk.translated.is_some() {
        // Source Block::length_sorted uses its 6N sequence indices to read
        // N IDs/block2oid, without copying the source DNA arena (block.cpp228).
        return Err(io::Error::other("Translated query length sorting indexes uninitialized source-query metadata in Source's six-context Block path."));
    }
    let mut order = (0..chunk.records.len()).collect::<Vec<_>>();
    order.sort_unstable_by(|&left, &right| {
        chunk.records[right]
            .sequence
            .len()
            .cmp(&chunk.records[left].sequence.len())
            .then_with(|| right.cmp(&left))
    });
    let mut records = std::mem::take(&mut chunk.records)
        .into_iter()
        .map(Some)
        .collect::<Vec<_>>();
    oids = order.iter().map(|&id| id as u64).collect();
    if !chunk.protein_oids.is_empty() {
        let original = std::mem::take(&mut chunk.protein_oids);
        chunk.protein_oids = order.iter().map(|&id| original[id]).collect();
    }
    chunk.records = order
        .into_iter()
        .map(|id| {
            let mut record = records[id].take().unwrap();
            // Source's replacement Block copies seqs/ids/OIDs, not qual_.
            record.quality = None;
            record
        })
        .collect();
    Ok(oids)
}

fn prepare_unsearched_query_chunk(
    config: &BlastpConfig,
    chunk: &mut QueryChunk,
    matrix: &ScoreMatrix,
    masker: Option<&crate::masking::tantan::TantanMasker>,
    prepare: bool,
    database_letters: u64,
) -> io::Result<Option<Vec<u64>>> {
    prepare_query_cache_with_overrides(
        config,
        chunk,
        matrix,
        masker,
        prepare,
        database_letters,
        crate::masking::MaskingAlgo::None,
        config.tuning.setup.freq_masking,
    )
}

fn prepare_query_cache_with_overrides(
    config: &BlastpConfig,
    chunk: &mut QueryChunk,
    matrix: &ScoreMatrix,
    masker: Option<&crate::masking::tantan::TantanMasker>,
    prepare: bool,
    database_letters: u64,
    soft_masking_override: crate::masking::MaskingAlgo,
    frequency_masking: bool,
) -> io::Result<Option<Vec<u64>>> {
    let oids = if prepare {
        let oids = multiprocessing::sort_join_only_queries(chunk, config)?;
        multiprocessing::prepare_join_only_queries(chunk, config, matrix, masker);
        (!config.tuning.query_pre_sorted).then_some(oids)
    } else {
        None
    };
    let mut sequences = if let Some(block) = chunk.translated.as_mut() {
        block.take_seqs()
    } else {
        let mut sequences = crate::data::sequence_set::SequenceSet::new();
        for record in &chunk.records {
            sequences.push(&record.sequence);
        }
        sequences
    };
    let setup = super::blastp::native_search_setup(
        config,
        matrix,
        config.tuning.setup.lin_stage1_query,
        config.round_linear_target || config.tuning.setup.lin_stage1_target,
        frequency_masking,
    )?;
    let result = super::search_algorithm::resolve_query(
        &mut sequences,
        chunk.translated.as_mut(),
        config,
        matrix,
        &setup,
        soft_masking_override,
        frequency_masking,
        database_letters,
    );
    if let Some(block) = chunk.translated.as_mut() {
        block.put_seqs(sequences);
    }
    result?;
    Ok(oids)
}

// Source self loads a protein database block directly into the query owner.
// Such blocks have no source_seqs arena: translated alignment later indexes
// source_seqs[query_id] and six-frame groups that were never initialized.
// Reject this undefined domain before constructing any query stream.
fn reject_translated_self_search(config: &BlastpConfig) -> io::Result<()> {
    if config.tuning.self_search && config.translated_query_layout.is_some() {
        return Err(io::Error::other("Translated --self search has no source-DNA query contexts in the database-owned Source path."));
    }
    Ok(())
}

fn run_protein(
    config: &BlastpConfig,
    references: ReferenceSource,
    input: QueryInput<'_>,
    budget: u64,
) -> io::Result<()> {
    run_protein_with_writer(config, references, input, budget, None)
}

fn run_protein_with_writer(
    config: &BlastpConfig,
    references: ReferenceSource,
    input: QueryInput<'_>,
    budget: u64,
    writer: Option<NativeSearchWriter>,
) -> io::Result<()> {
    run_protein_with_writer_retaining_references(config, references, input, budget, writer)
        .map(drop)
}

/// MCL retains the real database owner through final identifier output.
fn run_protein_with_writer_retaining_references(
    config: &BlastpConfig,
    references: ReferenceSource,
    mut input: QueryInput<'_>,
    budget: u64,
    writer: Option<NativeSearchWriter>,
) -> io::Result<ReferenceSource> {
    crate::data::accession_spool::source_scope(references.accession_spool(), || {
        if rayon::current_thread_index().is_none() {
            return super::blastp::install_search_threads(config.threads, move || {
                run_protein_with_writer_retaining_references(
                    config, references, input, budget, writer,
                )
            });
        }
        let mut normalized = super::blastp::normalize_output_config(config)?.into_owned();
        super::blastp::begin_algorithm_scope(&mut normalized)?;
        let config = &normalized;
        reject_translated_self_search(config)?;
        if config.tuning.self_search
            && !matches!(input, QueryInput::Protein(ProteinQueries::Database(_)))
        {
            input = QueryInput::Protein(ProteinQueries::Database(Box::new(
                references.query_owner(config)?,
            )));
        }
        if config.tuning.multiprocessing.enabled {
            return multiprocessing::run_retaining_references(config, references, input, budget);
        }
        if matches!(references, ReferenceSource::Memory { .. })
            && !config.seqidlist.as_os_str().is_empty()
        {
            return Err(io::Error::other(
                "The .dmnd database format does not support filtering by accession.",
            ));
        }
        if matches!(references, ReferenceSource::Memory { .. })
            && super::blastp::requested_search_taxonomy_flags(config)
                .contains(crate::data::sequence_file::SequenceFileFlags::TAXON_MAPPING)
        {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "taxonomy output is not available for in-memory clustering search",
            ));
        }
        let errors = crate::util::io::IoErrorOwner::default();
        let search_setup = errors.finish(crate::run::config::Config::from_input_with_io_owner(
            &crate::run::config::ConfigInput {
                iterate: config.iterate.clone(),
                self_search: config.tuning.self_search,
                lin_stage1_query: config.tuning.setup.lin_stage1_query,
                lin_stage1_target: config.tuning.setup.lin_stage1_target,
                sensitivity: config.sensitivity,
                target_indexed: config.target_indexed,
                query_translated: config.translated_query_layout.is_some(),
                query_cover: config.query_cover,
                subject_cover: config.subject_cover,
                min_length_ratio: config.tuning.min_length_ratio,
                global_ranking_targets: config.global_ranking_targets,
                frame_shift: config.frame_shift,
                algo: config.tuning.algorithm.source_request(
                    config.query_indexed,
                    config.contiguous_seed_mode,
                    config.swipe_all,
                ),
                ..Default::default()
            },
            &errors,
        ))?;
        if let Some(message) = &search_setup.iteration_message {
            search_message!("{message}");
        }
        native_double_indexed::run(
            config,
            references,
            input,
            budget,
            writer,
            search_setup.sensitivity,
            search_setup.min_length_ratio,
        )
    })
}
