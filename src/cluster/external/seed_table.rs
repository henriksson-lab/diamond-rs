//! Seed-table construction from `diamond/src/cluster/external/seed_table.cpp`.

use crate::util::io::source_c_path;
use std::cell::RefCell;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

use super::external::{ClusterStats, Job};
use super::pair_table::{Bucket, RadixedTable, SeedEntry, RADIX_BITS, RADIX_COUNT};
use super::storage::{BucketWriter, WorkerBucketBuffers};
use crate::basic::reduction::Reduction;
use crate::basic::seed_iterator::SketchIterator;
use crate::basic::shape::Shape;
use crate::basic::value::{Letter, SEED_MASK};
use crate::masking::motifs::mask_motifs;
use crate::masking::{bit_to_hard_mask, mask_sequence, MaskingAlgo, MaskingStat};
use crate::search::seed_complexity::seed_is_complex;
use crate::util::io::CompressedBuffer;
use crate::util::parallel::Atomic;
use crate::util::sequence::seqid_bytes;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Volume {
    pub path: PathBuf,
    pub oid_begin: u64,
    pub oid_end: u64,
    pub record_count: u64,
}

impl Volume {
    pub fn new(path: impl Into<PathBuf>, oid_begin: u64, oid_end: u64, record_count: u64) -> Self {
        Self {
            path: path.into(),
            oid_begin,
            oid_end,
            record_count,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct VolumedFile(pub Vec<Volume>);

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeedSequence<T = String> {
    pub id: T,
    pub sequence: Vec<Letter>,
}

pub type ByteSeedSequence = SeedSequence<Vec<u8>>;
impl From<SeedSequence> for ByteSeedSequence {
    fn from(record: SeedSequence) -> Self {
        Self {
            id: record.id.into_bytes(),
            sequence: record.sequence,
        }
    }
}
impl TryFrom<ByteSeedSequence> for SeedSequence {
    type Error = std::string::FromUtf8Error;
    fn try_from(record: ByteSeedSequence) -> Result<Self, Self::Error> {
        Ok(Self {
            id: String::from_utf8(record.id)?,
            sequence: record.sequence,
        })
    }
}

pub trait SequenceVolumeReader {
    /// Construct an independent reader for an actual Source volume worker.
    /// Legacy stateful library adapters remain explicitly serial when absent;
    /// concrete file readers share only error transport, never input cursors.
    fn fork_reader(&self) -> Option<Box<dyn SequenceVolumeReader + Send>> {
        None
    }
    fn read_volume(&mut self, path: &Path) -> Result<Vec<SeedSequence>, String>;
    fn visit_volume(
        &mut self,
        path: &Path,
        visitor: &mut dyn FnMut(SeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        for record in self.read_volume(path)? {
            visitor(record)?;
        }
        Ok(())
    }
    /// File-backed implementations override this to construct and retain the
    /// input owner before opening dependent output, including empty inputs.
    fn visit_volume_with_open(
        &mut self,
        path: &Path,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(SeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        on_open()?;
        self.visit_volume(path, visitor)
    }
    /// Canonical Source byte title path. Legacy String-only library readers
    /// adapt losslessly; concrete file-backed readers override these methods.
    fn visit_volume_bytes(
        &mut self,
        path: &Path,
        visitor: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        self.visit_volume(path, &mut |record| visitor(record.into()))
    }
    /// Legacy library readers can only adapt their eager read_volume contract.
    /// Actual file-backed readers override this to check before each input read.
    fn visit_volume_bytes_while(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        visitor: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        for record in self.read_volume(path)? {
            if !should_read() {
                break;
            }
            visitor(record.into())?;
        }
        Ok(())
    }
    fn visit_volume_bytes_with_open(
        &mut self,
        path: &Path,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        self.visit_volume_with_open(path, on_open, &mut |record| visitor(record.into()))
    }
    /// Concrete retained readers run the predicate before a record read and
    /// construct the input before on_open. Legacy readers keep their adapter.
    fn visit_volume_bytes_while_with_open(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        on_open()?;
        self.visit_volume_bytes_while(path, should_read, visitor)
    }
    /// Legacy adapters can release dependent owners only after their existing
    /// visitor returns. Concrete retained readers override the unwind boundary.
    fn visit_volume_bytes_while_with_cleanup(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        on_error: &mut dyn FnMut(),
        visitor: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        let result = self.visit_volume_bytes_while_with_open(path, should_read, on_open, visitor);
        if result.is_err() {
            on_error();
        }
        result
    }
}

#[derive(Debug, Default)]
pub struct FastaVolumeReader;

impl SequenceVolumeReader for FastaVolumeReader {
    fn visit_volume_bytes_while_with_cleanup(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        on_error: &mut dyn FnMut(),
        visitor: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        super::reader::visit_volume_bytes_while_owned_with_cleanup(
            path,
            should_read,
            visitor,
            on_open,
            on_error,
            &Default::default(),
        )
    }
    fn visit_volume_bytes_while_with_open(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        super::reader::visit_volume_bytes_while_owned(
            path,
            should_read,
            visitor,
            on_open,
            &Default::default(),
        )
    }
    fn fork_reader(&self) -> Option<Box<dyn SequenceVolumeReader + Send>> {
        Some(Box::new(Self))
    }
    fn visit_volume_bytes(
        &mut self,
        path: &Path,
        visitor: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        super::reader::visit_volume_bytes_owned(path, visitor, &mut || Ok(()), &Default::default())
    }
    fn visit_volume_bytes_while(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        visitor: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        super::reader::visit_volume_bytes_while_owned(
            path,
            should_read,
            visitor,
            &mut || Ok(()),
            &Default::default(),
        )
    }
    fn visit_volume_bytes_with_open(
        &mut self,
        path: &Path,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        super::reader::visit_volume_bytes_owned(path, visitor, on_open, &Default::default())
    }
    fn visit_volume(
        &mut self,
        path: &Path,
        visitor: &mut dyn FnMut(SeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        super::reader::visit_volume(path, visitor)
    }
    fn visit_volume_with_open(
        &mut self,
        path: &Path,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(SeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        super::reader::visit_volume_owned(
            path,
            visitor,
            on_open,
            &crate::util::io::IoErrorOwner::default(),
        )
    }
    fn read_volume(&mut self, path: &Path) -> Result<Vec<SeedSequence>, String> {
        let mut sequences = Vec::new();
        self.visit_volume(path, &mut |record| {
            sequences.push(record);
            Ok(())
        })?;
        Ok(sequences)
    }
}

pub trait SeedMasker {
    fn fork_masker(&self) -> Option<Box<dyn SeedMasker + Send>> {
        None
    }
    fn mask(&mut self, sequence: &mut [Letter], tantan: bool, motif: bool) -> MaskingStat;
}

#[derive(Debug, Default)]
pub struct DefaultSeedMasker;

impl SeedMasker for DefaultSeedMasker {
    fn fork_masker(&self) -> Option<Box<dyn SeedMasker + Send>> {
        Some(Box::new(Self))
    }
    fn mask(&mut self, sequence: &mut [Letter], tantan: bool, motif: bool) -> MaskingStat {
        let mut stats = MaskingStat::new();
        if tantan {
            mask_sequence(sequence, MaskingAlgo::Tantan);
            let masked = sequence
                .iter()
                .filter(|letter| (**letter & SEED_MASK) != 0)
                .count() as u64;
            bit_to_hard_mask(sequence);
            stats.add(MaskingAlgo::Tantan, masked);
        }
        if motif {
            let saved = mask_motifs(sequence);
            stats.add(MaskingAlgo::Motif, saved.len() as u64);
        }
        stats
    }
}

pub struct SeedTableConfig<'a> {
    pub threads: usize,
    pub soft_masking: &'a [u8],
    pub motif_masking: &'a [u8],
    pub sketch_size: i32,
    pub sensitivity_sketch_size: i32,
    pub seed_cut: f64,
    pub sensitivity_seed_cut: f64,
    pub shapes: &'a [Shape],
    pub reduction: &'a Reduction,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SeedTableStats {
    pub volumes_processed: usize,
    pub seeds_considered: u64,
    pub seeds_indexed: u64,
    pub masking_stat: MaskingStat,
}

/// FileArray-compatible seed output, using native packed records inside zlib
/// volumes and one `bucket.tsv` volume manifest per radix.
#[derive(Debug)]
pub struct SeedFileArray {
    base_dir: PathBuf,
    worker_id: i64,
    entries: Vec<Vec<SeedEntry>>,
    bounded: Option<RefCell<BucketWriter>>,
}

impl SeedFileArray {
    pub fn new(base_dir: impl AsRef<Path>, worker_id: i64) -> Result<Self, String> {
        let base_dir = base_dir.as_ref().to_path_buf();
        for radix in 0..RADIX_COUNT {
            fs::create_dir_all(source_c_path(&(base_dir.join(radix.to_string()))).as_ref())
                .map_err(|error| error.to_string())?;
        }
        Ok(Self {
            base_dir,
            worker_id,
            entries: vec![Vec::new(); RADIX_COUNT],
            bounded: None,
        })
    }

    pub fn new_bounded(base_dir: impl AsRef<Path>, worker_id: i64) -> Result<Self, String> {
        Self::new_bounded_owned(
            base_dir,
            worker_id,
            crate::util::io::IoErrorOwner::default(),
        )
    }

    pub fn new_bounded_owned(
        base_dir: impl AsRef<Path>,
        worker_id: i64,
        errors: crate::util::io::IoErrorOwner,
    ) -> Result<Self, String> {
        Ok(Self {
            base_dir: base_dir.as_ref().into(),
            worker_id,
            entries: vec![Vec::new(); RADIX_COUNT],
            bounded: Some(RefCell::new(BucketWriter::new_owned(
                base_dir.as_ref(),
                worker_id,
                RADIX_COUNT,
                errors,
            )?)),
        })
    }

    pub fn write_msb(&mut self, entry: SeedEntry) {
        let radix = (entry.key() >> (64 - RADIX_BITS)) as usize;
        if let Some(writer) = &self.bounded {
            let mut bytes = Vec::with_capacity(SeedEntry::ENCODED_SIZE);
            entry.serialize(&mut bytes);
            writer.borrow_mut().write_deferred(radix, &bytes);
        } else {
            self.entries[radix].push(entry);
        }
    }

    pub fn entries(&self, radix: usize) -> &[SeedEntry] {
        &self.entries[radix]
    }

    pub fn finish(&self) -> Result<RadixedTable, String> {
        if let Some(writer) = &self.bounded {
            return writer.borrow_mut().finish();
        }
        let mut buckets = Vec::with_capacity(RADIX_COUNT);
        for radix in 0..RADIX_COUNT {
            let directory = self.base_dir.join(radix.to_string());
            let manifest = directory.join("bucket.tsv");
            let entries = &self.entries[radix];
            if entries.is_empty() {
                fs::write(source_c_path(&manifest).as_ref(), [])
                    .map_err(|error| error.to_string())?;
            } else {
                let volume = directory.join(format!("worker_{}_volume_0", self.worker_id));
                let mut compressed = CompressedBuffer::new();
                for entry in entries {
                    let mut bytes = Vec::with_capacity(SeedEntry::ENCODED_SIZE);
                    entry.serialize(&mut bytes);
                    compressed
                        .write(&bytes)
                        .map_err(|error| error.to_string())?;
                }
                compressed.finish().map_err(|error| error.to_string())?;
                fs::write(source_c_path(&volume).as_ref(), compressed.data())
                    .map_err(|error| error.to_string())?;
                fs::write(
                    source_c_path(&manifest).as_ref(),
                    super::storage::manifest_line(&volume, &[entries.len() as u64])
                        .map_err(|e| e.to_string())?,
                )
                .map_err(|error| error.to_string())?;
            }
            buckets.push(Bucket::new(manifest, None));
        }
        Ok(RadixedTable(buckets))
    }
}

/// C++ `External::build_seed_table` with sequence-file and masker ownership
/// represented by explicit adapters.
pub fn build_seed_table<R, M>(
    job: &Job,
    volumes: &VolumedFile,
    shape_index: usize,
    config: &SeedTableConfig<'_>,
    reader: &mut R,
    masker: &mut M,
) -> Result<(RadixedTable, SeedTableStats), String>
where
    R: SequenceVolumeReader,
    M: SeedMasker,
{
    let shape = config
        .shapes
        .get(shape_index)
        .ok_or_else(|| format!("Shape index out of bounds: {shape_index}"))?;
    let first_shape = config
        .shapes
        .first()
        .ok_or_else(|| "No seed shapes configured.".to_owned())?;
    let use_tantan = config.soft_masking.is_empty() || config.soft_masking == b"tantan";
    let use_motif = config.motif_masking.is_empty() || config.motif_masking == b"1";
    let mut sketch_size = if config.sketch_size == 0 {
        config.sensitivity_sketch_size
    } else {
        config.sketch_size
    };
    let seed_cut = if config.seed_cut == 0.0 {
        config.sensitivity_seed_cut
    } else {
        config.seed_cut
    };
    let seed_complexity_cut = seed_cut * std::f64::consts::LN_2 * f64::from(first_shape.weight);
    if sketch_size == 0 {
        sketch_size = i32::MAX;
    }

    let base_dir = job.base_dir(None).join(format!("seed_table_{shape_index}"));
    fs::create_dir_all(source_c_path(&base_dir).as_ref()).map_err(|error| error.to_string())?;
    let accessions = job.root_dir().join("accessions");
    let write_accessions = job.round_index() == 0 && shape_index == 0;
    if write_accessions {
        fs::create_dir_all(source_c_path(&accessions).as_ref())
            .map_err(|error| error.to_string())?;
    }
    if config.threads == 0 {
        return Err("Seed-table thread count must be positive.".to_owned());
    }
    let settings = Arc::new(SeedWorkerSettings {
        shape: *shape,
        shape_index,
        shape_count: config.shapes.len(),
        sketch_size,
        complexity_cut: seed_complexity_cut,
        reduction: config.reduction.clone(),
        use_tantan,
        use_motif,
        write_accessions,
        accessions,
    });
    let output = Arc::new(Mutex::new(BucketWriter::new_file_array_owned(
        &base_dir,
        job.worker_id(),
        RADIX_COUNT,
        job.io_errors(),
    )?));
    let failure = Arc::new(Mutex::new(None));
    let queue = Arc::new(Atomic::new(base_dir.join("queue")));
    let totals = Arc::new(Mutex::new(SeedTableStats::default()));
    let factories = reader.fork_reader().zip(masker.fork_masker());
    if let Some((first_reader, first_masker)) = factories {
        let job = Arc::new(job.clone());
        let volumes = Arc::new(volumes.0.clone());
        let mut first = Some((first_reader, first_masker));
        let mut pool = crate::util::parallel::SimpleThreadPool::new();
        // Source starts all configured workers, even when some claim no volume.
        for _ in 0..config.threads {
            let (mut reader, mut masker) = if let Some(first) = first.take() {
                first
            } else {
                (
                    reader
                        .fork_reader()
                        .ok_or("Seed reader factory became unavailable.")?,
                    masker
                        .fork_masker()
                        .ok_or("Seed masker factory became unavailable.")?,
                )
            };
            let job = Arc::clone(&job);
            let volumes = Arc::clone(&volumes);
            let settings = Arc::clone(&settings);
            let queue = Arc::clone(&queue);
            let output = Arc::clone(&output);
            let failure = Arc::clone(&failure);
            let totals = Arc::clone(&totals);
            pool.spawn(move |stop| {
                let mut buffers =
                    WorkerBucketBuffers::new(output, RADIX_COUNT, job.io_errors(), failure);
                let mut encoded = Vec::with_capacity(SeedEntry::ENCODED_SIZE);
                let stats = run_seed_worker(
                    &job,
                    &volumes,
                    &settings,
                    &queue,
                    stop,
                    reader.as_mut(),
                    masker.as_mut(),
                    &mut |entry| {
                        encoded.clear();
                        let radix = (entry.key() >> (64 - RADIX_BITS)) as usize;
                        entry.serialize(&mut encoded);
                        buffers.write(radix, &encoded)
                    },
                )?;
                merge_seed_stats(&mut totals.lock().unwrap(), &stats);
                // BufferArray destruction flushes before the worker is joined.
                Ok(())
            });
        }
        pool.join_all()?;
    } else {
        // Legacy library adapters cannot manufacture independent cursors or
        // mutable masker state. The concrete command always uses the factories.
        let mut buffers = WorkerBucketBuffers::new(
            Arc::clone(&output),
            RADIX_COUNT,
            job.io_errors(),
            Arc::clone(&failure),
        );
        let mut encoded = Vec::with_capacity(SeedEntry::ENCODED_SIZE);
        let stats = run_seed_worker(
            job,
            &volumes.0,
            &settings,
            &queue,
            &AtomicBool::new(false),
            reader,
            masker,
            &mut |entry| {
                encoded.clear();
                let radix = (entry.key() >> (64 - RADIX_BITS)) as usize;
                entry.serialize(&mut encoded);
                buffers.write(radix, &encoded)
            },
        )?;
        merge_seed_stats(&mut totals.lock().unwrap(), &stats);
        drop(buffers);
    }
    if let Some(error) = failure.lock().unwrap().take() {
        return Err(error);
    }
    let stats = *totals.lock().unwrap();
    let buckets = output.lock().unwrap().finish()?;
    drop(output);
    let finished = Atomic::new(base_dir.join("finished"));
    finished.fetch_add(stats.volumes_processed as i64)?;
    finished.await_value(volumes.0.len() as i64)?;
    job.log_stats(&ClusterStats {
        seeds_considered: stats.seeds_considered,
        seeds_indexed: stats.seeds_indexed,
        masking_stat: stats.masking_stat,
        ..ClusterStats::default()
    })?;
    Ok((buckets, stats))
}

struct SeedWorkerSettings {
    shape: Shape,
    shape_index: usize,
    shape_count: usize,
    sketch_size: i32,
    complexity_cut: f64,
    reduction: Reduction,
    use_tantan: bool,
    use_motif: bool,
    write_accessions: bool,
    accessions: PathBuf,
}

fn merge_seed_stats(total: &mut SeedTableStats, local: &SeedTableStats) {
    total.volumes_processed += local.volumes_processed;
    total.seeds_considered += local.seeds_considered;
    total.seeds_indexed += local.seeds_indexed;
    total.masking_stat += local.masking_stat;
}

fn run_seed_worker(
    job: &Job,
    volumes: &[Volume],
    settings: &SeedWorkerSettings,
    queue: &Atomic,
    stop: &AtomicBool,
    reader: &mut dyn SequenceVolumeReader,
    masker: &mut dyn SeedMasker,
    write: &mut dyn FnMut(SeedEntry) -> Result<(), String>,
) -> Result<SeedTableStats, String> {
    let mut stats = SeedTableStats::default();
    // Source's monotonic resource retains reduced-sequence allocations until
    // the worker exits, across all its volumes, rather than per record.
    let mut reduction_arena = Vec::new();
    while !stop.load(Ordering::Relaxed) {
        let volume_index =
            usize::try_from(queue.fetch_add_one()?).map_err(|error| error.to_string())?;
        let Some(volume) = volumes.get(volume_index) else {
            break;
        };
        job.log(&format!(
            "Building seed table. Shape={}/{} Volume={}/{} Records={}",
            settings.shape_index + 1,
            settings.shape_count,
            volume_index + 1,
            volumes.len(),
            crate::util::string::format(volume.record_count),
        ))?;
        let accession_out = RefCell::new(None);
        let mut on_open = || {
            if settings.write_accessions {
                let path = settings.accessions.join(format!("{volume_index}.txt"));
                *accession_out.borrow_mut() = Some(crate::util::io::source_text::TextOutput::open(
                    &path,
                    &job.io_errors(),
                    b"Error opening file ",
                    true,
                )?);
            }
            Ok(())
        };
        let mut oid = volume.oid_begin as i64;
        reader.visit_volume_bytes_while_with_cleanup(
            &volume.path,
            &mut || !stop.load(Ordering::Relaxed),
            &mut on_open,
            &mut || {
                // Source unwinds acc_out before its enclosing input owner.
                accession_out.borrow_mut().take();
            },
            &mut |mut record| {
                if job.round_index() > 0 {
                    oid = c_atoll(&record.id);
                }
                if let Some(out) = accession_out.borrow_mut().as_mut() {
                    out.write(&seqid_bytes(&record.id));
                    out.endl();
                }
                if settings.use_tantan || settings.use_motif {
                    stats.masking_stat += masker.mask(
                        &mut record.sequence,
                        settings.use_tantan,
                        settings.use_motif,
                    );
                }
                reduction_arena.push(
                    record
                        .sequence
                        .iter()
                        .map(|letter| settings.reduction.reduce(*letter) as Letter)
                        .collect::<Vec<_>>(),
                );
                let reduced = reduction_arena.last().unwrap();
                if record.sequence.len() >= settings.shape.length as usize {
                    let mut seeds = SketchIterator::new(
                        reduced,
                        &settings.shape,
                        settings.sketch_size,
                        &settings.reduction,
                    );
                    while seeds.good() {
                        stats.seeds_considered += 1;
                        if seed_is_complex(
                            &record.sequence[seeds.pos() as usize..],
                            &settings.shape,
                            settings.complexity_cut,
                            &settings.reduction,
                        ) {
                            write(SeedEntry::new(
                                seeds.get(),
                                oid,
                                record.sequence.len() as i32,
                            ))?;
                            stats.seeds_indexed += 1;
                        }
                        seeds.increment();
                    }
                }
                oid = oid.wrapping_add(1);
                Ok(())
            },
        )?;
        // Source closes the sequence owner before the per-volume ofstream.
        drop(accession_out);
        stats.volumes_processed += 1;
    }
    Ok(stats)
}

fn c_atoll(value: impl AsRef<[u8]>) -> i64 {
    let mut terminated = value.as_ref().to_vec();
    terminated.push(0);
    unsafe { libc::atoll(terminated.as_ptr().cast()) as i64 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::reduction::Reduction;
    use crate::basic::value::AMINO_ACID_ALPHABET;
    use crate::util::io::InputFile;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

    fn temp_dir(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "diamond_seed_table_{name}_{}_{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    #[derive(Default)]
    struct FixtureReader {
        volumes: Vec<Vec<SeedSequence>>,
        next: usize,
    }

    impl SequenceVolumeReader for FixtureReader {
        fn read_volume(&mut self, _: &Path) -> Result<Vec<SeedSequence>, String> {
            let records = self.volumes[self.next].clone();
            self.next += 1;
            Ok(records)
        }
    }

    #[derive(Default)]
    struct RecordingMasker {
        calls: Vec<(bool, bool)>,
    }

    impl SeedMasker for RecordingMasker {
        fn mask(&mut self, sequence: &mut [Letter], tantan: bool, motif: bool) -> MaskingStat {
            self.calls.push((tantan, motif));
            let mut stats = MaskingStat::new();
            if !sequence.is_empty() {
                sequence[0] = crate::basic::value::MASK_LETTER;
                stats.add(MaskingAlgo::Tantan, 1);
            }
            stats
        }
    }

    fn collect_entries(output: &SeedFileArray) -> Vec<SeedEntry> {
        let mut entries: Vec<_> = output.entries.iter().flatten().copied().collect();
        entries.sort_unstable();
        entries
    }

    #[test]
    fn storage_routes_and_roundtrips_packed_compressed_records() {
        let root = temp_dir("storage");
        let mut output = SeedFileArray::new(&root, 7).unwrap();
        let entry = SeedEntry::new(123, 456, 789);
        output.write_msb(entry);
        let radix = (entry.key() >> 56) as usize;
        assert_eq!(output.entries(radix), &[entry]);
        let table = output.finish().unwrap();
        let manifest = fs::read_to_string(&table[radix].path).unwrap();
        assert!(manifest.ends_with("\t1\n"));
        let path = manifest.split_whitespace().next().unwrap();
        let mut input = InputFile::new(path, 0).unwrap();
        assert_eq!(input.read_u64().unwrap(), entry.seed);
        assert_eq!(input.read_i64().unwrap(), entry.oid);
        assert_eq!(input.read_i32().unwrap(), entry.len);
        input.close().unwrap();
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn seed_workers_keep_independent_real_inputs_and_flush_complete_records() {
        #[derive(Clone)]
        struct ConcurrentReader {
            opened: Arc<std::sync::Barrier>,
            threads: Arc<Mutex<std::collections::HashSet<std::thread::ThreadId>>>,
        }
        impl SequenceVolumeReader for ConcurrentReader {
            fn fork_reader(&self) -> Option<Box<dyn SequenceVolumeReader + Send>> {
                Some(Box::new(self.clone()))
            }
            fn read_volume(&mut self, _: &Path) -> Result<Vec<SeedSequence>, String> {
                Err("The concrete worker must retain its input stream.".into())
            }
            fn visit_volume_bytes_while_with_open(
                &mut self,
                path: &Path,
                should_read: &mut dyn FnMut() -> bool,
                on_open: &mut dyn FnMut() -> Result<(), String>,
                visitor: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
            ) -> Result<(), String> {
                super::super::reader::visit_volume_bytes_while_owned(
                    path,
                    should_read,
                    visitor,
                    &mut || {
                        self.threads
                            .lock()
                            .unwrap()
                            .insert(std::thread::current().id());
                        // Both actual file constructors have returned while
                        // their independent owners are still retained.
                        self.opened.wait();
                        on_open()
                    },
                    &Default::default(),
                )
            }
        }
        let root = temp_dir("parallel-inputs");
        fs::create_dir_all(&root).unwrap();
        let first = root.join("first.faa");
        let second = root.join("second.faa");
        fs::write(&first, b">a\xff description\nARND\n").unwrap();
        fs::write(&second, b">b\x80\nCQEG\n").unwrap();
        let job = Job::new(20, 2, 1 << 20, &root, 0).unwrap();
        let reduction = Reduction::new(
            "A R N D C Q E G H I L K M F P S T W Y V",
            AMINO_ACID_ALPHABET,
        );
        let shapes = [Shape::from_code("1", &reduction)];
        let config = SeedTableConfig {
            threads: 2,
            soft_masking: b"0",
            motif_masking: b"0",
            sketch_size: 64,
            sensitivity_sketch_size: 64,
            seed_cut: -1.0,
            sensitivity_seed_cut: -1.0,
            shapes: &shapes,
            reduction: &reduction,
        };
        let seen = Arc::new(Mutex::new(std::collections::HashSet::new()));
        let mut reader = ConcurrentReader {
            opened: Arc::new(std::sync::Barrier::new(2)),
            threads: Arc::clone(&seen),
        };
        let (table, stats) = build_seed_table(
            &job,
            &VolumedFile(vec![
                Volume::new(first, 10, 11, 1),
                Volume::new(second, 20, 21, 1),
            ]),
            0,
            &config,
            &mut reader,
            &mut DefaultSeedMasker,
        )
        .unwrap();
        assert_eq!(seen.lock().unwrap().len(), 2);
        assert_eq!(stats.volumes_processed, 2);
        assert_eq!(stats.seeds_considered, 8);
        assert_eq!(stats.seeds_indexed, 8);
        let mut oids = Vec::new();
        for bucket in table.iter() {
            for (path, count) in super::super::storage::manifest_records(&bucket.path).unwrap() {
                let mut input = InputFile::new_path(&path, 0).unwrap();
                for _ in 0..count {
                    let _seed = input.read_u64().unwrap();
                    oids.push(input.read_i64().unwrap());
                    assert_eq!(input.read_i32().unwrap(), 4);
                }
                input.close().unwrap();
            }
        }
        oids.sort_unstable();
        assert_eq!(oids, [10, 10, 10, 10, 20, 20, 20, 20]);
        for (index, expected) in [b"a\xff\n".as_slice(), b"b\x80\n".as_slice()]
            .into_iter()
            .enumerate()
        {
            assert_eq!(
                fs::read(root.join(format!("accessions/{index}.txt"))).unwrap(),
                crate::util::io::source_text::expected_text(expected)
            );
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn build_applies_masks_complexity_sketch_and_accession_rules() {
        let root = temp_dir("build");
        let job = Job::new(20, 1, 1 << 20, &root, 3).unwrap();
        let reduction = Reduction::new(
            "A R N D C Q E G H I L K M F P S T W Y V",
            AMINO_ACID_ALPHABET,
        );
        let shapes = [Shape::from_code("11", &reduction)];
        let config = SeedTableConfig {
            threads: 2,
            soft_masking: b"",
            motif_masking: b"",
            sketch_size: 0,
            sensitivity_sketch_size: 100,
            seed_cut: -1.0,
            sensitivity_seed_cut: 1.0,
            shapes: &shapes,
            reduction: &reduction,
        };
        let mut reader = FixtureReader {
            volumes: vec![vec![
                SeedSequence {
                    id: "seq0 description".to_owned(),
                    sequence: vec![0, 1, 2, 3],
                },
                SeedSequence {
                    id: "seq1".to_owned(),
                    sequence: vec![0],
                },
            ]],
            next: 0,
        };
        let mut masker = RecordingMasker::default();
        let volumes = VolumedFile(vec![Volume::new("unused", 10, 12, 2)]);
        let (table, stats) =
            build_seed_table(&job, &volumes, 0, &config, &mut reader, &mut masker).unwrap();

        assert_eq!(table.len(), RADIX_COUNT);
        assert_eq!(masker.calls, [(true, true), (true, true)]);
        assert_eq!(stats.volumes_processed, 1);
        assert_eq!(stats.seeds_considered, 2);
        assert_eq!(stats.seeds_indexed, 2);
        assert_eq!(stats.masking_stat.get(MaskingAlgo::Tantan), 2);
        assert_eq!(
            fs::read(root.join("accessions/0.txt")).unwrap(),
            crate::util::io::source_text::expected_text(b"seq0\nseq1\n")
        );

        let mut entries = Vec::new();
        for bucket in table.iter() {
            let manifest = fs::read_to_string(&bucket.path).unwrap();
            if let Some(path) = manifest.split_whitespace().next() {
                let mut input = InputFile::new(path, 0).unwrap();
                let count = manifest
                    .split_whitespace()
                    .nth(1)
                    .unwrap()
                    .parse::<usize>()
                    .unwrap();
                for _ in 0..count {
                    entries.push(SeedEntry::new(
                        input.read_u64().unwrap(),
                        input.read_i64().unwrap(),
                        input.read_i32().unwrap(),
                    ));
                }
            }
        }
        entries.sort_unstable();
        assert!(entries
            .iter()
            .all(|entry| entry.oid == 10 && entry.len == 4));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn later_round_uses_atoll_ids_and_does_not_write_accessions() {
        let root = temp_dir("round1");
        let mut job = Job::new(100, 1, 1 << 20, &root, 1).unwrap();
        job.next_round().unwrap();
        let reduction = Reduction::default_reduction();
        let shapes = [Shape::from_code("1", &reduction)];
        let config = SeedTableConfig {
            threads: 1,
            soft_masking: b"none",
            motif_masking: b"0",
            sketch_size: 10,
            sensitivity_sketch_size: 0,
            seed_cut: -1.0,
            sensitivity_seed_cut: 0.0,
            shapes: &shapes,
            reduction: &reduction,
        };
        let mut reader = FixtureReader {
            volumes: vec![vec![SeedSequence {
                id: "  -42 suffix".to_owned(),
                sequence: vec![0, 1],
            }]],
            next: 0,
        };
        let mut masker = RecordingMasker::default();
        let volumes = VolumedFile(vec![Volume::new("unused", 8, 9, 1)]);
        let (_, stats) =
            build_seed_table(&job, &volumes, 0, &config, &mut reader, &mut masker).unwrap();
        assert_eq!(stats.seeds_indexed, 2);
        assert!(masker.calls.is_empty());
        assert!(!root.join("accessions").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn effective_defaults_and_shape_errors_are_checked() {
        let reduction = Reduction::default_reduction();
        let root = temp_dir("errors");
        let job = Job::new(0, 0, 1, &root, 0).unwrap();
        let config = SeedTableConfig {
            threads: 1,
            soft_masking: b"none",
            motif_masking: b"0",
            sketch_size: 0,
            sensitivity_sketch_size: 0,
            seed_cut: 0.0,
            sensitivity_seed_cut: 0.5,
            shapes: &[],
            reduction: &reduction,
        };
        let error = build_seed_table(
            &job,
            &VolumedFile::default(),
            0,
            &config,
            &mut FixtureReader::default(),
            &mut RecordingMasker::default(),
        )
        .unwrap_err();
        assert_eq!(error, "Shape index out of bounds: 0");
        assert_eq!(c_atoll(" +123tail"), 123);
        assert_eq!(c_atoll("garbage"), 0);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn collect_helper_observes_native_seed_order() {
        let root = temp_dir("collect");
        let mut output = SeedFileArray::new(&root, 0).unwrap();
        output.write_msb(SeedEntry::new(2, 1, 10));
        output.write_msb(SeedEntry::new(1, 2, 20));
        assert_eq!(collect_entries(&output)[0].seed, 1);
        fs::remove_dir_all(root).unwrap();
    }
}
