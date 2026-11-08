use super::search_output::SearchFormat;
pub use super::search_output::{ArchiveOptions, SamOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use rayon::iter::{IndexedParallelIterator, ParallelIterator};

use crate::align::hsp::{Hsp, Match};
use crate::align::target::{extend as extend_targets, GappedScoreConfig};
use crate::align::ungapped::UngappedStageConfig;
use crate::basic::packed_loc::PackedLoc;
use crate::basic::seed::{seedp_count, seedp_mask};
use crate::basic::shape::Shape;
use crate::basic::shape_config::ShapeConfig;
use crate::basic::statistics::Statistics;
use crate::basic::value::{Letter, SequenceType, TaxId, LETTER_MASK};
use crate::chaining::HammingExtConfig;
use crate::config::Sensitivity;
use crate::data::blastdb::blastdb::BlastDB;
use crate::data::blastdb::volume::SequenceFileFlags;
use crate::data::block::Block;
use crate::data::dmnd::{DatabaseFile, DatabaseOpenOptions};
use crate::data::fasta;
use crate::data::seed_histogram::SeedPartitionRange;
use crate::data::sequence_set::SequenceSet;
use crate::data::taxonomy::{Rank, TaxonomyNode, TaxonomyTree};
use crate::dp::swipe::{Flags, HspValues};
use crate::masking::{MaskingAlgo, MaskingMode};
use crate::output::daa::daa_write::{
    finish_daa_from_sequence_file_bytes, finish_daa_query_record, init_daa,
    write_daa_query_record_bytes, write_daa_record_hsp, DaaRunMetadata, DaaSequenceFile,
};
use crate::output::format::{self, FieldId};
use crate::search::hit::Hit;
use crate::search::hit_buffer::{CompactHit, HitBuffer, HitBufferMode};
use crate::search::left_most::Context as LeftMostContext;
use crate::search::left_most_unclipped::left_most_filter_with_range_backed;
use crate::search::{parallel, sensitivity};
use crate::stats::cbs::CbsMode;
use crate::stats::score_matrix::{CutoffTable2D, ScoreMatrix};
use crate::util::algo::PatternMatcher;
use crate::util::system::trim_freed_heap_pages;

#[inline]
fn stage2_query_bounds(query_len: usize, seed_pos: usize, window: usize) -> (usize, usize) {
    (
        seed_pos.saturating_sub(window),
        seed_pos.saturating_add(window).min(query_len),
    )
}

/// C++ sequence equality for the self-hit predicate ignores the internal
/// masking annotation bits and compares the biological residue codes.
fn same_sequence(left: &[Letter], right: &[Letter]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(&left, &right)| left & LETTER_MASK == right & LETTER_MASK)
}

#[derive(Clone, Copy)]
struct Stage2Filter {
    score: i32,
    window: usize,
    left_most_interval: i32,
    self_search: bool,
}

#[derive(Clone, Copy)]
struct StoredHit {
    subject: u64,
    query_id: u32,
    seed_offset: u32,
    score: u16,
}

struct SeedIndexState {
    query: Option<Block>,
    shapes: ShapeConfig,
    target_seeds: Option<crate::data::seed_set::HashedSeedSet>,
    query_seeds: Option<std::sync::Arc<crate::data::seed_set::HashedSeedSet>>,
    contiguous_seeds: Option<crate::data::seed_set::SeedSet>,
    encoding: crate::data::flags::SeedEncoding,
    seedp_bits: i32,
}

#[inline]
fn packed_locs_from_bytes(bytes: &[u8]) -> &[PackedLoc] {
    let elem_size = std::mem::size_of::<PackedLoc>();
    assert_eq!(elem_size, 5);
    assert_eq!(std::mem::align_of::<PackedLoc>(), 1);
    assert_eq!(bytes.len() % elem_size, 0);
    // SAFETY: `PackedLoc` is `repr(C, packed)`, has alignment one, and is
    // exactly the `high: u8` plus native-endian `low: u32` representation
    // written into DoubleArray by `HashJoinValue::write_ne_bytes`. Every bit
    // pattern is valid for those integer fields and the returned slice cannot
    // outlive the borrowed DoubleArray range.
    unsafe {
        std::slice::from_raw_parts(bytes.as_ptr().cast::<PackedLoc>(), bytes.len() / elem_size)
    }
}

/// Adapt the source-faithful `DoubleArray<PackedLoc>` join result to the
/// native stage-1 worker. Each callback still contains exactly one joined
/// seed partition. Borrow one natural joined-key group in place, run the
/// translated tiled Hamming kernel, and pass its hit indices together with
/// the original target locations to the query-major ungapped worker. This is
/// the same borrowed boundary as upstream `JoinIterator`/`search_query_offset`.
fn map_target_indexed_partition<F>(
    query_hits: &crate::util::data_structures::DoubleArray,
    target_hits: &crate::util::data_structures::DoubleArray,
    absolute_partition: u32,
    query_seqs: &SequenceSet,
    target_seqs: &SequenceSet,
    hamming_filter_id: u32,
    tile_size: usize,
    dispatch: &crate::search::hamming::Stage1DispatchConfig,
    retain_ids: bool,
    ranking: Option<&crate::search::kmer_ranking::KmerRanking>,
    mut map_matches: F,
) -> Result<(Vec<StoredHit>, usize), String>
where
    F: FnMut((u32, u32), &[u32], usize, &[PackedLoc], u32, &mut Vec<StoredHit>),
{
    use crate::search::hamming::{
        run_stage1_packed_loc, run_stage1_packed_loc_id, HitField, WorkSet,
    };
    use crate::search::kmer_ranking::PackedLocId;
    let width = if retain_ids { 9 } else { 5 };
    let mut query = query_hits.begin_with_elem_size(width);
    let mut target = target_hits.begin_with_elem_size(width);
    let mut output = Vec::new();
    let mut raw = 0usize;
    let mut joined_keys = 0u64;
    let mut hamming_passed = 0u64;
    while query.good() && target.good() {
        let stats = if retain_ids {
            let decode = |bytes: &[u8]| -> Vec<PackedLocId> {
                bytes
                    .chunks_exact(9)
                    .map(|value| {
                        PackedLocId::new(
                            packed_locs_from_bytes(&value[..5])[0].as_u64(),
                            u32::from_ne_bytes(value[5..9].try_into().unwrap()),
                        )
                    })
                    .collect()
            };
            let queries = decode(query.range_bytes(width));
            let targets = decode(target.range_bytes(width));
            let target_locs: Vec<_> = targets.iter().map(|hit| PackedLoc::new(hit.loc)).collect();
            let mut work = WorkSet::new(
                query_seqs,
                target_seqs,
                tile_size,
                hamming_filter_id,
                ranking,
                |hits: &mut HitField,
                 query_begin: usize,
                 target_begin: usize,
                 queries: &[PackedLocId],
                 _: &[PackedLocId]| {
                    for local in 0..hits.query_count() {
                        let (id, position) =
                            query_seqs.local_position(queries[query_begin + local].loc as usize);
                        let indices = hits.hits(local);
                        if !indices.is_empty() {
                            hamming_passed = hamming_passed.wrapping_add(indices.len() as u64);
                            map_matches(
                                (id as u32, position as u32),
                                indices,
                                target_begin,
                                &target_locs,
                                absolute_partition,
                                &mut output,
                            );
                        }
                    }
                },
            );
            run_stage1_packed_loc_id(
                [(queries.as_slice(), targets.as_slice())],
                &mut work,
                dispatch,
            )
            .map_err(|error| error.to_string())?;
            work.stats
        } else {
            let queries = packed_locs_from_bytes(query.range_bytes(width));
            let targets = packed_locs_from_bytes(target.range_bytes(width));
            let mut work = WorkSet::new(
                query_seqs,
                target_seqs,
                tile_size,
                hamming_filter_id,
                ranking,
                |hits: &mut HitField,
                 query_begin: usize,
                 target_begin: usize,
                 queries: &[PackedLoc],
                 targets: &[PackedLoc]| {
                    for local in 0..hits.query_count() {
                        let (id, position) = query_seqs
                            .local_position(queries[query_begin + local].as_u64() as usize);
                        let indices = hits.hits(local);
                        if !indices.is_empty() {
                            hamming_passed = hamming_passed.wrapping_add(indices.len() as u64);
                            map_matches(
                                (id as u32, position as u32),
                                indices,
                                target_begin,
                                targets,
                                absolute_partition,
                                &mut output,
                            );
                        }
                    }
                },
            );
            run_stage1_packed_loc([(queries, targets)], &mut work, dispatch)
                .map_err(|error| error.to_string())?;
            work.stats
        };
        raw = raw.saturating_add(stats.seed_hits as usize);
        joined_keys = joined_keys.wrapping_add(stats.seeds_hit);
        query.increment(width);
        target.increment(width);
    }
    if query.good() != target.good() {
        return Err("Mismatched joined seed groups.".into());
    }
    let mut stat = Statistics::new();
    stat.inc(
        crate::basic::statistics::StatValue::SeedsHit,
        joined_keys as i64,
    );
    stat.inc(crate::basic::statistics::StatValue::SeedHits, raw as i64);
    stat.inc(
        crate::basic::statistics::StatValue::TentativeMatches1,
        hamming_passed as i64,
    );
    crate::basic::statistics::merge_statistics(&stat);
    Ok((output, raw))
}

type UngappedKernel = fn(&[Letter], &[&[Letter]], usize, &ScoreMatrix, &mut [i32]);

type SubjectIndexPartitionFilter = for<'a> unsafe fn(
    (u32, u32),
    &[u32],
    &[u64],
    u32,
    &mut Vec<StoredHit>,
    &SequenceSet,
    &SequenceSet,
    &[Stage2Filter],
    &ScoreMatrix,
    &Shape,
    &LeftMostContext<'a>,
    bool,
    bool,
    bool,
    bool,
    usize,
    u32,
    usize,
    UngappedKernel,
) -> usize;

type BorrowedPartitionFilter = for<'a> unsafe fn(
    (u32, u32),
    &[u32],
    usize,
    &[PackedLoc],
    u32,
    &mut Vec<StoredHit>,
    &SequenceSet,
    &SequenceSet,
    &[Stage2Filter],
    &ScoreMatrix,
    &Shape,
    &LeftMostContext<'a>,
    bool,
    bool,
    bool,
    bool,
    usize,
    u32,
    usize,
    UngappedKernel,
) -> usize;

const UNGAPPED_PORTABLE: u8 = 0;
const UNGAPPED_AVX2: u8 = 1;
const UNGAPPED_SSE41: u8 = 2;

#[inline]
fn run_ungapped_kernel<const KERNEL: u8>(
    query: &[Letter],
    subjects: &[&[Letter]],
    window: usize,
    score_matrix: &ScoreMatrix,
    out: &mut [i32],
) {
    if KERNEL == UNGAPPED_AVX2 {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        // SAFETY: this function item is selected only after the one-time AVX2
        // check at the blastp run boundary.
        unsafe {
            crate::dp::simd_ungapped::window_ungapped_best_into_avx2(
                query,
                subjects,
                window,
                score_matrix,
                out,
            );
            return;
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        unreachable!("AVX2 specialization on a non-x86 target");
    }
    if KERNEL == UNGAPPED_SSE41 {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        // SAFETY: this function item is selected only after the one-time
        // SSSE3/SSE4.1 checks at the blastp run boundary.
        unsafe {
            crate::dp::simd_ungapped::window_ungapped_best_into_sse41(
                query,
                subjects,
                window,
                score_matrix,
                out,
            );
            return;
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        unreachable!("SSE4.1 specialization on a non-x86 target");
    }
    crate::dp::simd_ungapped::window_ungapped_best_into(query, subjects, window, score_matrix, out);
}

/// Hit storage that keeps the fast compact representation while RSS is below
/// the requested ceiling, then permanently switches to DIAMOND's compressed
/// temporary-bin format. The ceiling is necessarily soft: loaded input,
/// indices, allocator bookkeeping, and one active disk bin are irreducible.
#[repr(C, packed)]
#[derive(Clone, Copy)]
struct PackedRetainedHit {
    query: u32,
    subject: PackedLoc,
    seed_offset: u32,
    score: u16,
}

struct PackedHitBins {
    bins: Vec<Vec<PackedRetainedHit>>,
    key_partition: Vec<u32>,
}

impl PackedHitBins {
    fn new(key_partition: Vec<u32>) -> Self {
        Self {
            bins: (0..key_partition.len()).map(|_| Vec::new()).collect(),
            key_partition,
        }
    }

    fn source_memory_reserve_hits() -> usize {
        // Source HitBuffer reserves 1 * GIGABYTES / sizeof(Hit) only for
        // trace_pt_membuf. Its pinned packed Hit layout is the same 15 bytes.
        (1usize << 30) / std::mem::size_of::<PackedRetainedHit>()
    }

    fn reserve_source_memory(&mut self) -> io::Result<()> {
        for bin in &mut self.bins {
            bin.try_reserve_exact(Self::source_memory_reserve_hits())
                .map_err(|_| crate::util::io::BadAllocError::default().into_io())?;
        }
        Ok(())
    }

    fn bin(&self, query: u32) -> Option<usize> {
        self.key_partition.iter().position(|&end| query < end)
    }

    fn push(&mut self, hit: StoredHit) {
        if let Some(bin) = self.bin(hit.query_id) {
            self.bins[bin].push(PackedRetainedHit {
                query: hit.query_id,
                subject: PackedLoc::new(hit.subject),
                seed_offset: hit.seed_offset,
                score: hit.score,
            });
        }
    }

    fn total_hits(&self) -> usize {
        self.bins.iter().map(Vec::len).sum()
    }

    fn take_grouped_bin(&mut self, bin: usize) -> (Vec<Vec<CompactHit>>, u32, u32) {
        let begin = if bin == 0 {
            0
        } else {
            self.key_partition[bin - 1]
        };
        let end = self.key_partition[bin];
        let packed = std::mem::take(&mut self.bins[bin]);
        let mut counts = vec![0usize; end.saturating_sub(begin) as usize];
        for hit in &packed {
            let query = hit.query;
            if let Some(count) = counts.get_mut(query.saturating_sub(begin) as usize) {
                *count += 1;
            }
        }
        let mut grouped: Vec<Vec<CompactHit>> =
            counts.into_iter().map(Vec::with_capacity).collect();
        for hit in packed {
            let query = hit.query;
            if let Some(query_hits) = grouped.get_mut(query.saturating_sub(begin) as usize) {
                query_hits.push(CompactHit {
                    subject: hit.subject.as_u64(),
                    seed_offset: hit.seed_offset,
                    score: hit.score,
                });
            }
        }
        (grouped, begin, end)
    }
}

struct AdaptiveHitStore {
    load_threads: i32,
    memory: Option<PackedHitBins>,
    disk: Option<HitBuffer>,
    memory_limit: Option<u64>,
    tmpdir: PathBuf,
    query_count: usize,
    key_partition: Vec<u32>,
    max_subject: u64,
    bytes_since_rss_check: usize,
    rss_checked_once: bool,
    no_unlink: bool,
}

/// Actual retained-hit storage dependency of the shared AlignBackend driver.
/// It borrows the existing owner, so prefetched decoded bins and their recycler
/// stay alive through DP without a flat-hit copy or another buffering layer.
struct NativeGroupedHitSource<'a> {
    store: &'a mut AdaptiveHitStore,
    next_memory_bin: usize,
    pending_memory_bin: Option<usize>,
}

impl crate::align::native::GroupedHitSource for NativeGroupedHitSource<'_> {
    fn prepare(&mut self) -> Result<(), String> {
        if let Some(disk) = self.store.disk.as_mut() {
            disk.try_finish_writing()?;
        }
        Ok(())
    }
    fn load(&mut self, max_bytes: u64) -> bool {
        if let Some(memory) = self.store.memory.as_ref() {
            if self.pending_memory_bin.is_some() || self.next_memory_bin >= memory.bins.len() {
                return false;
            }
            self.pending_memory_bin = Some(self.next_memory_bin);
            self.next_memory_bin += 1;
            true
        } else {
            self.store
                .disk
                .as_mut()
                .expect("retained-hit storage")
                .load_grouped_with_limit(max_bytes)
        }
    }
    fn retrieve(&mut self) -> Result<Option<crate::align::native::GroupedHitBatch>, String> {
        let data = if let Some(memory) = self.store.memory.as_mut() {
            self.pending_memory_bin
                .take()
                .map(|bin| memory.take_grouped_bin(bin))
        } else {
            self.store
                .disk
                .as_mut()
                .ok_or("missing hit storage")?
                .try_retrieve_grouped_owned()?
        };
        Ok(data.map(
            |(hits, query_begin, query_end)| crate::align::native::GroupedHitBatch {
                hits,
                query_begin,
                query_end,
            },
        ))
    }
    fn next_bin_size(&self) -> u64 {
        if let Some(memory) = self.store.memory.as_ref() {
            memory
                .bins
                .get(self.next_memory_bin)
                .map_or(0, |bin| bin.len() as u64)
        } else {
            self.store.disk.as_ref().map_or(0, HitBuffer::next_bin_size)
        }
    }
    fn recycle(&mut self, batch: crate::align::native::GroupedHitBatch) {
        if let Some(disk) = self.store.disk.as_mut() {
            disk.recycle_grouped_hits(batch.hits);
        }
    }
    fn disk_size(&self) -> i64 {
        self.store
            .disk
            .as_ref()
            .map_or(0, HitBuffer::total_disk_size)
    }
    fn finish(&mut self) -> Result<(), String> {
        // All bin vectors have moved through owned batches by this boundary.
        // Release the consumed memory owner before command footer emission.
        self.store.memory.take();
        if let Some(disk) = self.store.disk.as_mut() {
            disk.free_buffer();
            crate::util::log_stream::message_stream()
                .lock()
                .unwrap()
                .write(format!(
                    "Hit buffer: read {:.1} MiB of compressed spill data\n",
                    disk.total_disk_size() as f64 / (1024.0 * 1024.0)
                ))
                .map_err(|error| error.to_string())?;
        }
        Ok(())
    }
}

impl AdaptiveHitStore {
    fn update_ranking_shape(
        &mut self,
        state: &mut NativeRankingState,
        queries: &SequenceSet,
        target: &Block,
        config: &BlastpConfig,
        matrix: &ScoreMatrix,
    ) -> io::Result<()> {
        let memory = self
            .memory
            .as_mut()
            .ok_or_else(|| io::Error::other("global ranking requires its per-shape seed buffer"))?;
        let mut hits = Some(
            memory
                .bins
                .iter_mut()
                .flat_map(std::mem::take)
                .map(|hit| {
                    Hit::with_score(hit.query, hit.subject.as_u64(), hit.seed_offset, hit.score)
                })
                .collect(),
        );
        let sequences: Vec<_> = (0..queries.len()).map(|id| queries.get(id)).collect();
        crate::align::global_ranking::update_table(
            &mut hits,
            &mut state.table,
            &crate::align::global_ranking::UpdateTableConfig {
                query_contexts: if config.translated_query_layout.is_some() {
                    6
                } else {
                    1
                },
                global_ranking_targets: config.global_ranking_targets as usize,
                threads: rayon::current_num_threads(),
                no_reextend: config.tuning.no_reextend,
                xdrop: matrix.rawscore_int(config.ungapped_xdrop_bits),
                query_sequences: &sequences,
                target_block: target,
                score_matrix: matrix,
            },
        )
        .map_err(io::Error::other)?;
        Ok(())
    }

    #[cfg(test)]
    fn new(
        query_count: usize,
        query_bins: usize,
        max_subject: u64,
        query_group_size: usize,
        memory_limit: Option<u64>,
        tmpdir: PathBuf,
    ) -> Self {
        Self::new_partitioned(
            query_count,
            max_subject,
            memory_limit,
            tmpdir,
            Self::key_partition(query_count, query_bins, query_group_size),
        )
    }

    fn new_partitioned(
        query_count: usize,
        max_subject: u64,
        memory_limit: Option<u64>,
        tmpdir: PathBuf,
        key_partition: Vec<u32>,
    ) -> Self {
        Self {
            memory: Some(PackedHitBins::new(key_partition.clone())),
            load_threads: 16,
            disk: None,
            memory_limit,
            tmpdir,
            query_count,
            key_partition,
            max_subject,
            bytes_since_rss_check: 0,
            rss_checked_once: false,
            no_unlink: false,
        }
    }

    fn configure_storage(&mut self, source_memory: bool) -> io::Result<()> {
        // Native's explicit zero ceiling requests disk even when the Source
        // memory flag is supplied. Other soft ceilings describe live RSS;
        // they are never a request to reserve the ceiling before any hit.
        if self.memory_limit == Some(0) {
            self.spill_to_disk()?;
        } else if source_memory {
            self.memory.as_mut().unwrap().reserve_source_memory()?;
            // Source explicitly chooses persistent memory with --hit-membuf.
            // A positive native adaptive ceiling does not silently undo it.
            self.memory_limit = None;
        }
        Ok(())
    }

    fn key_partition(query_count: usize, query_bins: usize, query_group_size: usize) -> Vec<u32> {
        let query_end = query_count.max(1).min(u32::MAX as usize) as u32;
        let query_group_size = query_group_size.max(1);
        let group_count = query_count.div_ceil(query_group_size).max(1);
        let num_bins = group_count.clamp(1, query_bins.max(1));
        let mut key_partition = Vec::with_capacity(num_bins);
        for bin in 1..=num_bins {
            let end_group = (bin * group_count).div_ceil(num_bins);
            let end = end_group
                .saturating_mul(query_group_size)
                .max(1)
                .min(query_end as usize) as u32;
            if key_partition.last().copied() != Some(end) {
                key_partition.push(end);
            }
        }
        if key_partition.last().copied() != Some(query_end) {
            key_partition.push(query_end);
        }
        key_partition
    }

    fn ingest_partitions(&mut self, partitions: Vec<Vec<StoredHit>>) -> io::Result<usize> {
        let incoming = partitions.iter().map(Vec::len).sum::<usize>();
        let incoming_bytes = incoming.saturating_mul(std::mem::size_of::<CompactHit>());
        self.bytes_since_rss_check = self.bytes_since_rss_check.saturating_add(incoming_bytes);
        const RSS_CHECK_GRANULARITY: usize = 8 * 1024 * 1024;
        if self.disk.is_none() && self.memory_limit.is_some() {
            let should_check =
                !self.rss_checked_once || self.bytes_since_rss_check >= RSS_CHECK_GRANULARITY;
            if should_check {
                self.rss_checked_once = true;
                self.bytes_since_rss_check = 0;
                let limit = self.memory_limit.unwrap();
                if (crate::util::system::get_current_rss() as u64)
                    .saturating_add(incoming_bytes as u64)
                    >= limit
                {
                    self.spill_to_disk()?;
                }
            }
        }

        if let Some(memory) = self.memory.as_mut() {
            for partition in partitions {
                for hit in partition {
                    memory.push(hit);
                }
            }
        } else {
            self.write_partitions(partitions)?;
        }
        Ok(incoming)
    }

    fn spill_to_disk(&mut self) -> io::Result<()> {
        if self.disk.is_some() {
            return Ok(());
        }
        // Match upstream's sensitivity-specific query-bin count (16 for most
        // modes, 64 for ultra-sensitive). Each bin is decoded as a unit during
        // alignment, so retaining the default-mode count in ultra-sensitive
        // searches multiplies the live decoded-hit footprint by about four.
        let query_end = self.query_count.max(1).min(u32::MAX as usize) as u32;
        let key_partition = self.key_partition.clone();
        let mut disk = HitBuffer::with_limits(
            key_partition,
            &self.tmpdir,
            self.max_subject > u32::MAX as u64,
            1,
            query_end,
            self.max_subject.max(1),
            HitBufferMode::Disk,
        )
        .map_err(io::Error::other)?;

        disk.set_no_unlink(self.no_unlink);
        disk.set_load_threads(self.load_threads, rayon::current_num_threads().max(1));
        let mut memory = self.memory.take().expect("memory hit store missing");
        let migrated_hits = memory.total_hits();
        for bin in 0..memory.bins.len() {
            let (grouped, begin, _) = memory.take_grouped_bin(bin);
            for (offset, hits) in grouped.into_iter().enumerate() {
                let query = begin as usize + offset;
                for hit in hits {
                    disk.append_disk_hit(Hit::with_score(
                        query as u32,
                        hit.subject,
                        hit.seed_offset,
                        hit.score,
                    ))
                    .map_err(io::Error::other)?;
                }
            }
        }
        disk.take_error().map_err(io::Error::other)?;
        self.disk = Some(disk);
        let mut message = format!(
            "Hit buffer: RSS reached --memory-limit; spilling {migrated_hits} retained hits to "
        )
        .into_bytes();
        let path = if self.tmpdir.as_os_str().is_empty() {
            Path::new(".")
        } else {
            &self.tmpdir
        };
        message.extend_from_slice(
            &crate::util::io::source_os_bytes(path.as_os_str()).map_err(io::Error::other)?,
        );
        message.push(b'\n');
        crate::util::log_stream::message_stream()
            .lock()
            .unwrap()
            .write_bytes(&message)?;
        trim_freed_heap_pages();
        Ok(())
    }

    fn write_partitions(&mut self, partitions: Vec<Vec<StoredHit>>) -> io::Result<()> {
        let disk = self
            .disk
            .as_mut()
            .ok_or_else(|| io::Error::other("disk hit buffer not initialized"))?;
        for partition in partitions {
            for hit in partition {
                disk.append_disk_hit(Hit::with_score(
                    hit.query_id,
                    hit.subject,
                    hit.seed_offset,
                    hit.score,
                ))
                .map_err(io::Error::other)?;
            }
        }
        disk.take_error().map_err(io::Error::other)
    }
}

#[inline(always)]
fn filter_subject_indices_partition_to_hits(
    query: (u32, u32),
    target_indices: &[u32],
    target_subjects: &[u64],
    absolute_partition: u32,
    out: &mut Vec<StoredHit>,
    queries: &SequenceSet,
    targets: &SequenceSet,
    cutoffs: &[Stage2Filter],
    score_matrix: &ScoreMatrix,
    shape: &Shape,
    left_most_context: &LeftMostContext<'_>,
    first_shape: bool,
    chunked: bool,
    use_left_most_range: bool,
    skip_left_most: bool,
    index_chunks: usize,
    min_identities: u32,
    ungapped_lane_count: usize,
    ungapped_kernel: UngappedKernel,
) -> usize {
    let candidates = BorrowedStage1Subjects {
        indices: target_indices,
        subjects: target_subjects,
    };
    macro_rules! run {
        ($first:literal, $chunked:literal, $range:literal, $skip:literal) => {
            filter_query_candidates_impl::<
                BorrowedStage1Subjects<'_>,
                $first,
                $chunked,
                $range,
                $skip,
            >(
                query,
                &candidates,
                absolute_partition,
                out,
                queries,
                targets,
                cutoffs,
                score_matrix,
                shape,
                left_most_context,
                index_chunks,
                min_identities,
                ungapped_lane_count,
                ungapped_kernel,
            )
        };
    }
    if skip_left_most {
        return run!(false, false, false, true);
    }
    match (first_shape, chunked, use_left_most_range) {
        (true, true, true) => run!(true, true, true, false),
        (false, true, true) => run!(false, true, true, false),
        (true, true, false) => run!(true, true, false, false),
        (false, true, false) => run!(false, true, false, false),
        (true, false, false) => run!(true, false, false, false),
        (false, false, false) => run!(false, false, false, false),
        (_, false, true) => unreachable!("left-most range without index chunks"),
    }
}

#[inline]
unsafe fn filter_subject_indices_partition_to_hits_baseline(
    query: (u32, u32),
    target_indices: &[u32],
    target_subjects: &[u64],
    absolute_partition: u32,
    out: &mut Vec<StoredHit>,
    queries: &SequenceSet,
    targets: &SequenceSet,
    cutoffs: &[Stage2Filter],
    score_matrix: &ScoreMatrix,
    shape: &Shape,
    left_most_context: &LeftMostContext<'_>,
    first_shape: bool,
    chunked: bool,
    use_left_most_range: bool,
    skip_left_most: bool,
    index_chunks: usize,
    min_identities: u32,
    ungapped_lane_count: usize,
    ungapped_kernel: UngappedKernel,
) -> usize {
    filter_subject_indices_partition_to_hits(
        query,
        target_indices,
        target_subjects,
        absolute_partition,
        out,
        queries,
        targets,
        cutoffs,
        score_matrix,
        shape,
        left_most_context,
        first_shape,
        chunked,
        use_left_most_range,
        skip_left_most,
        index_chunks,
        min_identities,
        ungapped_lane_count,
        ungapped_kernel,
    )
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn filter_subject_indices_partition_to_hits_avx2(
    query: (u32, u32),
    target_indices: &[u32],
    target_subjects: &[u64],
    absolute_partition: u32,
    out: &mut Vec<StoredHit>,
    queries: &SequenceSet,
    targets: &SequenceSet,
    cutoffs: &[Stage2Filter],
    score_matrix: &ScoreMatrix,
    shape: &Shape,
    left_most_context: &LeftMostContext<'_>,
    first_shape: bool,
    chunked: bool,
    use_left_most_range: bool,
    skip_left_most: bool,
    index_chunks: usize,
    min_identities: u32,
    ungapped_lane_count: usize,
    ungapped_kernel: UngappedKernel,
) -> usize {
    filter_subject_indices_partition_to_hits(
        query,
        target_indices,
        target_subjects,
        absolute_partition,
        out,
        queries,
        targets,
        cutoffs,
        score_matrix,
        shape,
        left_most_context,
        first_shape,
        chunked,
        use_left_most_range,
        skip_left_most,
        index_chunks,
        min_identities,
        ungapped_lane_count,
        ungapped_kernel,
    )
}

#[inline(always)]
fn filter_borrowed_partition_to_hits(
    query: (u32, u32),
    target_indices: &[u32],
    target_begin: usize,
    target_locations: &[PackedLoc],
    absolute_partition: u32,
    out: &mut Vec<StoredHit>,
    queries: &SequenceSet,
    targets: &SequenceSet,
    cutoffs: &[Stage2Filter],
    score_matrix: &ScoreMatrix,
    shape: &Shape,
    left_most_context: &LeftMostContext<'_>,
    first_shape: bool,
    chunked: bool,
    use_left_most_range: bool,
    skip_left_most: bool,
    index_chunks: usize,
    min_identities: u32,
    ungapped_lane_count: usize,
    ungapped_kernel: UngappedKernel,
) -> usize {
    let candidates = BorrowedStage1Targets {
        indices: target_indices,
        target_begin,
        targets: target_locations,
    };
    macro_rules! run {
        ($first:literal, $chunked:literal, $range:literal, $skip:literal) => {
            filter_query_candidates_impl::<
                BorrowedStage1Targets<'_>,
                $first,
                $chunked,
                $range,
                $skip,
            >(
                query,
                &candidates,
                absolute_partition,
                out,
                queries,
                targets,
                cutoffs,
                score_matrix,
                shape,
                left_most_context,
                index_chunks,
                min_identities,
                ungapped_lane_count,
                ungapped_kernel,
            )
        };
    }
    if skip_left_most {
        return run!(false, false, false, true);
    }
    match (first_shape, chunked, use_left_most_range) {
        (true, true, true) => run!(true, true, true, false),
        (false, true, true) => run!(false, true, true, false),
        (true, true, false) => run!(true, true, false, false),
        (false, true, false) => run!(false, true, false, false),
        (true, false, false) => run!(true, false, false, false),
        (false, false, false) => run!(false, false, false, false),
        (_, false, true) => unreachable!("left-most range without index chunks"),
    }
}

#[inline]
unsafe fn filter_borrowed_partition_to_hits_baseline(
    query: (u32, u32),
    target_indices: &[u32],
    target_begin: usize,
    target_locations: &[PackedLoc],
    absolute_partition: u32,
    out: &mut Vec<StoredHit>,
    queries: &SequenceSet,
    targets: &SequenceSet,
    cutoffs: &[Stage2Filter],
    score_matrix: &ScoreMatrix,
    shape: &Shape,
    left_most_context: &LeftMostContext<'_>,
    first_shape: bool,
    chunked: bool,
    use_left_most_range: bool,
    skip_left_most: bool,
    index_chunks: usize,
    min_identities: u32,
    ungapped_lane_count: usize,
    ungapped_kernel: UngappedKernel,
) -> usize {
    filter_borrowed_partition_to_hits(
        query,
        target_indices,
        target_begin,
        target_locations,
        absolute_partition,
        out,
        queries,
        targets,
        cutoffs,
        score_matrix,
        shape,
        left_most_context,
        first_shape,
        chunked,
        use_left_most_range,
        skip_left_most,
        index_chunks,
        min_identities,
        ungapped_lane_count,
        ungapped_kernel,
    )
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn filter_borrowed_partition_to_hits_avx2(
    query: (u32, u32),
    target_indices: &[u32],
    target_begin: usize,
    target_locations: &[PackedLoc],
    absolute_partition: u32,
    out: &mut Vec<StoredHit>,
    queries: &SequenceSet,
    targets: &SequenceSet,
    cutoffs: &[Stage2Filter],
    score_matrix: &ScoreMatrix,
    shape: &Shape,
    left_most_context: &LeftMostContext<'_>,
    first_shape: bool,
    chunked: bool,
    use_left_most_range: bool,
    skip_left_most: bool,
    index_chunks: usize,
    min_identities: u32,
    ungapped_lane_count: usize,
    ungapped_kernel: UngappedKernel,
) -> usize {
    filter_borrowed_partition_to_hits(
        query,
        target_indices,
        target_begin,
        target_locations,
        absolute_partition,
        out,
        queries,
        targets,
        cutoffs,
        score_matrix,
        shape,
        left_most_context,
        first_shape,
        chunked,
        use_left_most_range,
        skip_left_most,
        index_chunks,
        min_identities,
        ungapped_lane_count,
        ungapped_kernel,
    )
}

trait Stage2Subject {
    fn len(&self) -> usize;
    fn absolute_subject(&self, index: usize) -> u64;
}

struct BorrowedStage1Subjects<'a> {
    indices: &'a [u32],
    subjects: &'a [u64],
}

impl Stage2Subject for BorrowedStage1Subjects<'_> {
    #[inline(always)]
    fn len(&self) -> usize {
        self.indices.len()
    }

    #[inline(always)]
    fn absolute_subject(&self, index: usize) -> u64 {
        self.subjects[self.indices[index] as usize]
    }
}

struct BorrowedStage1Targets<'a> {
    indices: &'a [u32],
    target_begin: usize,
    targets: &'a [PackedLoc],
}

impl Stage2Subject for BorrowedStage1Targets<'_> {
    #[inline(always)]
    fn len(&self) -> usize {
        self.indices.len()
    }

    #[inline(always)]
    fn absolute_subject(&self, index: usize) -> u64 {
        self.targets[self.target_begin + self.indices[index] as usize].as_u64()
    }
}

#[inline(always)]
fn filter_query_candidates_impl<
    C: Stage2Subject + ?Sized,
    const FIRST_SHAPE: bool,
    const CHUNKED: bool,
    const USE_LEFT_MOST_RANGE: bool,
    const SKIP_LEFT_MOST: bool,
>(
    query_location: (u32, u32),
    candidates: &C,
    absolute_partition: u32,
    out: &mut Vec<StoredHit>,
    queries: &SequenceSet,
    targets: &SequenceSet,
    cutoffs: &[Stage2Filter],
    score_matrix: &ScoreMatrix,
    shape: &Shape,
    left_most_context: &LeftMostContext<'_>,
    index_chunks: usize,
    min_identities: u32,
    ungapped_lane_count: usize,
    ungapped_kernel: UngappedKernel,
) -> usize {
    let ref_seq_data = targets.data();
    let mut ungapped_count = 0usize;
    let mut left_most_passed = 0usize;
    let current_range = if USE_LEFT_MOST_RANGE {
        let partitions = (left_most_context.seedp_mask + 1) as usize;
        let chunk_size = partitions.div_ceil(index_chunks);
        let begin = (absolute_partition as usize / chunk_size) * chunk_size;
        let end = (begin + chunk_size).min(partitions);
        Some(SeedPartitionRange::with_bounds(begin as u32, end as u32))
    } else {
        None
    };
    let query_id = query_location.0 as usize;
    let q_pos = query_location.1 as usize;
    if query_id >= queries.len() {
        return 0;
    }
    let query = queries.get(query_id);
    let filter = cutoffs[query_id];
    let (q_start, q_end) = stage2_query_bounds(query.len(), q_pos, filter.window);
    let window_left = q_pos - q_start;
    let window_clipped = q_end - q_start;
    let query_window = &query[q_start..q_end];
    let cutoff = filter.score;
    let interval_mod = if filter.left_most_interval > 0 {
        (q_pos % filter.left_most_interval as usize) as i32
    } else {
        window_left as i32
    };
    let interval_overhang = (window_left as i32 - interval_mod).max(0) as usize;
    let left_q_start = q_start + interval_overhang;
    let left_seed_offset = window_left.saturating_sub(interval_overhang);
    if left_q_start >= q_end {
        return 0;
    }
    let left_len = q_end - left_q_start;

    for chunk_begin in (0..candidates.len()).step_by(ungapped_lane_count) {
        let chunk_len = ungapped_lane_count.min(candidates.len() - chunk_begin);
        let mut subject_starts = [0isize; 32];
        for (local_index, start) in subject_starts[..chunk_len].iter_mut().enumerate() {
            *start = candidates.absolute_subject(chunk_begin + local_index) as isize
                - window_left as isize;
        }
        let subject_storage;
        let mut subject_windows = [&[][..]; 32];
        if subject_starts[..chunk_len].iter().all(|&start| start >= 0) {
            for (window, &start) in subject_windows.iter_mut().zip(&subject_starts) {
                *window = &ref_seq_data[start as usize..];
            }
        } else {
            subject_storage = subject_starts[..chunk_len]
                .iter()
                .map(|&start| {
                    (0..window_clipped)
                        .map(|n| {
                            let pos = start + n as isize;
                            if pos < 0 {
                                crate::basic::value::DELIMITER_LETTER
                            } else {
                                ref_seq_data
                                    .get(pos as usize)
                                    .copied()
                                    .unwrap_or(crate::basic::value::DELIMITER_LETTER)
                            }
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            for (window, storage) in subject_windows.iter_mut().zip(&subject_storage) {
                *window = storage;
            }
        }
        let mut scores = [i32::MAX; 32];
        if cutoff != 0 {
            ungapped_kernel(
                query_window,
                &subject_windows[..chunk_len],
                window_clipped,
                score_matrix,
                &mut scores[..chunk_len],
            );
        }
        for (hit_index, &score) in scores[..chunk_len].iter().enumerate() {
            if score <= cutoff {
                continue;
            }
            let subject = (subject_starts[hit_index] + window_left as isize) as usize;
            // Source stage2 rejects the aliased query's own target after the
            // strict ungapped cutoff, before TENTATIVE_MATCHES2/left-most.
            if filter.self_search && targets.local_position(subject).0 as usize == query_id {
                continue;
            }
            ungapped_count += 1;
            let left_subject_start =
                subject as isize - window_left as isize + interval_overhang as isize;
            if !SKIP_LEFT_MOST {
                let keep = left_most_filter_with_range_backed(
                    queries.data(),
                    queries.position(query_id, left_q_start),
                    left_len,
                    ref_seq_data,
                    left_subject_start as usize,
                    left_seed_offset as i32,
                    shape.length,
                    left_most_context,
                    FIRST_SHAPE,
                    shape,
                    cutoff,
                    CHUNKED,
                    min_identities,
                    current_range,
                );
                if !keep {
                    continue;
                }
            }
            left_most_passed += 1;
            out.push(StoredHit {
                subject: subject as u64,
                query_id: query_location.0,
                seed_offset: query_location.1,
                score: if score == i32::MAX {
                    u16::MAX
                } else {
                    score.min(u16::MAX as i32) as u16
                },
            });
        }
    }
    let mut stat = Statistics::new();
    stat.inc(
        crate::basic::statistics::StatValue::TentativeMatches2,
        ungapped_count as i64,
    );
    stat.inc(
        crate::basic::statistics::StatValue::TentativeMatches3,
        left_most_passed as i64,
    );
    crate::basic::statistics::merge_statistics(&stat);
    ungapped_count
}

/// Port of `align.cpp::make_partition` for translated queries. Hit records
/// are already grouped by context, so count all six contexts of a source and
/// extend every threshold crossing through that complete source query.
#[cfg(test)]
fn translated_hit_partitions(
    hits_by_context: &[Vec<CompactHit>],
    source_count: usize,
    query_contexts: usize,
    min_task_trace_points: usize,
) -> Vec<std::ops::Range<usize>> {
    crate::align::native::grouped_partitions(
        hits_by_context,
        source_count,
        query_contexts,
        min_task_trace_points,
    )
}

// Upstream Tantan starts treating 50k residues as an oversized sequence (its
// reusable buffers have a 50k-residue floor).  Queries beyond that size also
// make the final SWIPE pass allocate unusually large temporary profiles and
// traceback matrices.  Keeping 64 waves of those temporaries in the allocator
// before returning to the writer can make freed pages dominate RSS, especially
// for blastx where the six translated contexts are all long.  Bound only this
// unusual case to one wave of queries per worker; ordinary protein workloads
// retain the wider batch used for load balancing.
const LONG_QUERY_RESIDUES: usize = 50_000;

fn output_query_batch_size<T>(records: &[fasta::FastaRecord<T>], threads: usize) -> (usize, bool) {
    let threads = threads.max(1);
    let has_long_query = records
        .iter()
        .any(|record| record.sequence.len() > LONG_QUERY_RESIDUES);
    (
        if has_long_query {
            threads
        } else {
            threads * 64
        },
        has_long_query,
    )
}

pub(crate) fn native_hamming_config(
    config: &BlastpConfig,
    approx_min_id: f64,
    query_or_target_cover: f64,
) -> HammingExtConfig {
    HammingExtConfig {
        hamming_ext: config.tuning.hamming_ext,
        approx_min_id,
        query_or_target_cover,
        query_cover: config.query_cover,
        subject_cover: config.subject_cover,
        max_evalue: config.max_evalue,
        diag_filter_id: config.tuning.diag_filter_id,
        diag_filter_cov: config.tuning.diag_filter_cov,
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MultiprocessingOptions {
    pub enabled: bool,
    pub initialize: bool,
    pub recover: bool,
    pub directory: PathBuf,
    pub query_chunk: i64,
    pub self_compare: bool,
}
impl Default for MultiprocessingOptions {
    fn default() -> Self {
        Self {
            enabled: false,
            initialize: false,
            recover: false,
            directory: PathBuf::new(),
            query_chunk: -1,
            self_compare: false,
        }
    }
}

/// Source formatting controls shared by tabular and contextual output.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SearchOutputOptions {
    pub cluster_threshold: Option<f64>,
    pub clustering_filter: Option<Arc<crate::mcl::clustering_format::ClusteringFilter>>,
    pub qnum_offset: i64,
    pub snum_offset: i64,
    pub xml_blord_format: bool,
    pub no_parse_seqids: bool,
}

/// Logical source block indices for in-process clustered searches.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SearchBlockContext {
    pub query: u32,
    pub reference: u32,
    pub self_search: bool,
}

/// Source search controls consumed by the actual seed/filter/extension stages.
#[derive(Debug, Clone, PartialEq)]
pub struct SearchTuning {
    pub runtime: super::source_config::Runtime,
    pub store_query_quality: bool,
    pub algorithm: super::search_algorithm::AlgorithmControl,
    pub heartbeat: bool,
    pub load_threads: i32,
    pub minichunk: usize,
    pub unaligned_targets: std::path::PathBuf,
    /// Global Source option; zero asks Search::Config to derive it from coverage.
    pub min_length_ratio: f64,
    /// Search::Config output retained across every sensitivity round/block.
    pub search_min_length_ratio: Option<f64>,
    /// The master has already applied Source query sorting before its header.
    pub query_pre_sorted: bool,
    pub rank_ratio: f64,
    pub rank_factor: f64,
    pub output: SearchOutputOptions,
    pub multiprocessing: MultiprocessingOptions,
    pub search_block: SearchBlockContext,
    pub setup: sensitivity::SetupSearchInput,
    pub matrix_file: std::path::PathBuf,
    pub stop_match_score: i32,
    pub matrix_scale: i32,
    pub db_size: u64,
    pub no_ranking: bool,
    pub no_reextend: bool,
    pub self_search: bool,
    pub add_self_aln: bool,
    pub mapany: bool,
    pub target_hard_cap: i64,
    pub ranking_score_drop_factor: f64,
    pub ranking_cutoff_bitscore: f64,
    pub load_balancing: u32,
    pub hash_join: crate::util::algo::HashJoinConfig,
    pub log_query: bool,
    pub log_extend: bool,
    pub target_parallel_verbosity: u32,
    pub threads_align: u32,
    pub min_task_trace_pts: i64,
    pub chaining: crate::chaining::GreedyAlignConfig,
    pub hamming_ext: bool,
    pub diag_filter_id: Option<f64>,
    pub diag_filter_cov: Option<f64>,
    /// Source OptionValue presence disables diagonal filtering, including zero.
    pub mutual_cover: Option<f64>,
    pub trace_pt_fetch_size: u64,
    pub kmer_ranking: bool,
    pub swipe_task_size: i64,
    pub swipe_chunk_size: u32,
    pub transcript_len_estimate: f64,
    pub band_bin: i32,
    pub col_bin: i32,
    pub narrow_band_cov: f64,
    pub narrow_band_factor: f64,
    pub min_band_overlap: f64,
    pub cutoff_score_8bit: i32,
    pub max_swipe_dp: i64,
    pub approx_backtrace: bool,
    pub tantan_min_mask_prob: f64,
    pub cbs_window: i32,
    pub max_motif_len: i32,
    pub tile_size: usize,
    pub padding: Option<i32>,
    pub ignore_warnings: bool,
    pub no_unlink: bool,
    pub no_block_size_limit: bool,
    pub hit_membuf: bool,
    pub short_query_ungapped_bitscore: f64,
    pub short_query_max_len: i32,
    pub ungapped_window: i32,
    pub left_most_interval: i32,
    pub gapped_filter_evalue1: f64,
    pub gapped_filter_window: i32,
    pub matrix_adjust_thresholds: crate::stats::cbs::CbsThresholds,
    pub cbs_err_tolerance: f64,
    pub cbs_it_limit: i32,
    // Source's sole clipping call is commented out; retain parsed settings
    // without enabling a different extension algorithm.
    pub anchor_window: i32,
    pub anchor_score: f64,
}
impl Default for SearchTuning {
    fn default() -> Self {
        Self {
            runtime: Default::default(),
            store_query_quality: false,
            algorithm: Default::default(),
            unaligned_targets: std::path::PathBuf::new(),
            min_length_ratio: 0.0,
            search_min_length_ratio: None,
            query_pre_sorted: false,
            rank_ratio: -1.0,
            rank_factor: -1.0,
            heartbeat: false,
            load_threads: 16,
            minichunk: 262_144,
            output: Default::default(),
            multiprocessing: Default::default(),
            search_block: Default::default(),
            setup: sensitivity::SetupSearchInput {
                gapped_filter_diag_bit_score: 12.0,
                ..Default::default()
            },
            matrix_file: std::path::PathBuf::new(),
            stop_match_score: 1,
            matrix_scale: 1,
            db_size: 0,
            no_ranking: false,
            no_reextend: false,
            self_search: false,
            add_self_aln: false,
            mapany: false,
            target_hard_cap: 0,
            ranking_score_drop_factor: 0.95,
            ranking_cutoff_bitscore: 25.0,
            load_balancing: 0,
            hash_join: Default::default(),
            log_query: false,
            log_extend: false,
            target_parallel_verbosity: u32::MAX,
            threads_align: 0,
            min_task_trace_pts: 1024,
            chaining: crate::chaining::GreedyAlignConfig::default(),
            hamming_ext: false,
            diag_filter_id: None,
            diag_filter_cov: None,
            mutual_cover: None,
            trace_pt_fetch_size: 10_000_000_000,
            kmer_ranking: false,
            swipe_task_size: 100_000_000,
            swipe_chunk_size: 256,
            transcript_len_estimate: 1.0,
            band_bin: 24,
            col_bin: 400,
            narrow_band_cov: 0.0,
            narrow_band_factor: 0.0,
            min_band_overlap: 0.0,
            cutoff_score_8bit: 240,
            max_swipe_dp: 1_000_000,
            approx_backtrace: false,
            tantan_min_mask_prob: 0.9,
            cbs_window: 40,
            max_motif_len: 30,
            tile_size: 1024,
            padding: None,
            ignore_warnings: false,
            no_unlink: false,
            no_block_size_limit: false,
            hit_membuf: false,
            short_query_ungapped_bitscore: 25.0,
            short_query_max_len: 60,
            ungapped_window: 48,
            left_most_interval: 32,
            gapped_filter_evalue1: 2000.0,
            gapped_filter_window: 200,
            matrix_adjust_thresholds: Default::default(),
            cbs_err_tolerance: 1e-8,
            cbs_it_limit: 2000,
            anchor_window: 12,
            anchor_score: 1.0,
        }
    }
}

/// Configuration for a blastp run.
#[derive(Clone)]
pub struct BlastpConfig {
    pub query_files: Vec<std::path::PathBuf>,
    pub database: std::path::PathBuf,
    pub output: Option<std::path::PathBuf>,
    pub matrix: Vec<u8>,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub max_evalue: f64,
    pub min_bit_score: f64,
    pub max_target_seqs: i64,
    /// Maximum HSPs per query/target pair; zero means unlimited.
    pub max_hsps: u32,
    /// Source --unal setting: -1 selects the format default.
    pub report_unaligned: i32,
    /// Include scientific-name lineage in taxonomy format 102.
    pub include_lineage: bool,
    /// Absence differs from a bare --header (the search workflow default).
    pub output_header: Option<Vec<Vec<u8>>>,
    /// Source compression option: empty/0=plain, 1=gzip, zstd=optional build.
    pub compression: Vec<u8>,
    pub file_buffer_size: usize,
    pub sam: SamOptions,
    pub archive: ArchiveOptions,
    pub tuning: SearchTuning,
    /// Command-local scoring snapshot, shared across database chunks.
    pub scoring_matrix: Option<std::sync::Arc<ScoreMatrix>>,
    pub auxiliary_output: super::query_output::AuxiliaryOutputOptions,
    pub taxonlist: Vec<u8>,
    pub taxon_exclude: Vec<u8>,
    pub seqidlist: std::path::PathBuf,
    pub skip_missing_seqids: bool,
    pub range_culling: bool,
    pub range_cover: f64,
    pub culling_overlap: f64,
    pub taxon_k: u64,
    pub query_or_subject_cover: f64,
    pub reference_taxon_rank_ids: Option<std::sync::Arc<[Vec<u32>]>>,
    pub ext_chunk_size: i64,
    pub toppercent: Option<f64>,
    pub global_ranking_targets: i64,
    /// Frameshift penalty; zero selects ordinary independent-frame extension.
    pub frame_shift: i32,
    /// Explicit source extension mode; None selects its workflow/sensitivity default.
    pub extension_mode: Option<crate::search::sensitivity::ExtensionMode>,
    /// Source sensitivity rounds; Some(empty) selects the default schedule.
    pub iterate: Option<Vec<Vec<u8>>>,
    /// Internal current-round target linearization and immutable source-query skip snapshot.
    pub round_linear_target: bool,
    pub round_query_skip: Option<std::sync::Arc<crate::util::data_structures::BitVector>>,
    pub min_id: f64,
    /// Minimum score-derived approximate identity percentage.
    pub approx_min_id: f64,
    /// Explicit seed-only masking, added to sensitivity-default motif masking.
    pub soft_masking: crate::masking::MaskingAlgo,
    pub threads: i32,
    /// Seed-index partitions processed per pass; zero selects the source default.
    pub index_chunks: u32,
    /// Query/reference letter budget in billions; None uses the source default.
    pub block_size: Option<f64>,
    /// Force hashed query-indexed seed filtering (--algo 1).
    pub query_indexed: bool,
    /// Source --algo ctg: dense query seed filter and contiguous encoding.
    pub contiguous_seed_mode: bool,
    /// Exhaustively align every query against every reference without seeds.
    pub swipe_all: bool,
    pub outfmt: Vec<Vec<u8>>,
    pub sensitivity: Sensitivity,
    pub masking: MaskingMode,
    pub motif_masking: Vec<u8>,
    pub min_query_len: i32,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub comp_based_stats: CbsMode,
    pub no_self_hits: bool,
    /// Ungapped extension x-drop threshold in bits (`--xdrop`).
    /// C++ default is 12.3 (`diamond/src/basic/config.cpp:427`) and converts
    /// it to a raw score with the active score matrix at startup.
    pub ungapped_xdrop_bits: f64,
    /// Soft process-RSS ceiling. Once reached, retained seed hits are stored in
    /// compressed temporary bins. `None` keeps the all-in-memory fast path.
    pub memory_limit: Option<u64>,
    /// Directory for spill files; an empty path uses the OS temporary directory.
    pub tmpdir: PathBuf,
    /// Consume the persistent `<database>.seed_idx` built by `makeidx`.
    pub target_indexed: bool,
    /// Internal blastx layout. When present, every source query contributes
    /// six consecutive protein records (one per translated context).
    pub translated_query_layout: Option<TranslatedQueryLayout>,
}

#[derive(Clone, Debug)]
pub struct TranslatedQueryLayout {
    /// Synchronized two-file blastx mate relation. IDs, source DNA, and all
    /// translated contexts live in the authoritative query `Block`.
    pub mate_indices: Vec<Option<usize>>,
}

struct TranslatedExtensionGroup<'a> {
    id: &'a [u8],
    dna_len: i32,
    source_sequence: &'a [Letter],
    quality: &'a [u8],
    mate_sequence: Option<&'a [Letter]>,
    frames: [&'a [Letter]; 6],
}
impl crate::align::native::QueryTitle for TranslatedExtensionGroup<'_> {
    fn query_title_bytes(&self) -> &[u8] {
        self.id
    }
}

/// Owned boundary between native block searching and the source block merger.
/// The search phase moves its final query arena here; completed matches retain
/// their actual HSP values and transcripts rather than reparsed text rows.
#[derive(Default)]
pub(crate) struct NativeSearchCapture {
    /// Freshly computed scores; empty when the caller's prepared query owner
    /// already carries the actual Source self-alignment table.
    pub query_self_aln_scores: Vec<f64>,
    pub target_self_aln_scores: Vec<f64>,
    pub queries: Vec<fasta::ByteFastaRecord>,
    pub query_oids: Vec<u64>,
    pub target_oids: Vec<u64>,
    pub translated: Option<Block>,
    pub matches: Vec<Vec<Match>>,
    pub seeded: Vec<bool>,
    pub matrix: Option<ScoreMatrix>,
    pub ranking: Option<NativeRankingState>,
}

pub(crate) struct NativeRankingState {
    pub table: Vec<crate::align::global_ranking::Hit>,
    pub reference_begin: u64,
    pub reference_oids: Option<std::sync::Arc<[u64]>>,
}

#[derive(Clone, Copy)]
pub(crate) struct ClusterSearchOptions<'a> {
    pub approx_min_id: f64,
    pub linear_query: bool,
    pub linear_target: bool,
    pub self_search: bool,
    pub query_or_target_cover: f64,
    pub mutual_cover: bool,
    pub hamming: HammingExtConfig,
    pub masker: Option<&'a crate::masking::tantan::TantanMasker>,
    pub soft_masking: MaskingAlgo,
    pub freq_masking: bool,
    pub extension: Option<sensitivity::ExtensionMode>,
}

/// Capture an actual clustering block with the same overridden seed/DP
/// settings as run_edges_in_memory, retaining matches for Source block joining.
pub(crate) fn search_cluster_block(
    config: &BlastpConfig,
    database: Vec<fasta::ByteFastaRecord>,
    queries: Option<Vec<fasta::ByteFastaRecord>>,
    database_letters: u64,
    queries_prepared: bool,
    options: ClusterSearchOptions<'_>,
) -> io::Result<NativeSearchCapture> {
    search_cluster_block_ranked(
        config,
        database,
        queries,
        database_letters,
        queries_prepared,
        options,
        None,
    )
}

pub(crate) fn search_cluster_block_ranked(
    config: &BlastpConfig,
    database: Vec<fasta::ByteFastaRecord>,
    queries: Option<Vec<fasta::ByteFastaRecord>>,
    database_letters: u64,
    queries_prepared: bool,
    options: ClusterSearchOptions<'_>,
    ranking: Option<NativeRankingState>,
) -> io::Result<NativeSearchCapture> {
    let capture = std::sync::Mutex::new(NativeSearchCapture {
        ranking,
        ..Default::default()
    });
    // Select the existing edge HSP/alias path; capture takes ownership before
    // formatting, so this sink never retains a second set of output records.
    let edge_sink = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    run_impl(
        config,
        Some(database),
        queries,
        Some(edge_sink),
        Some(options.approx_min_id),
        options.soft_masking,
        options.freq_masking,
        Some(database_letters),
        options.linear_query,
        options.linear_target,
        options.self_search,
        options.query_or_target_cover,
        options.mutual_cover,
        Some(options.hamming),
        options.masker,
        options.extension,
        None,
        Some(&capture),
        queries_prepared,
    )?;
    capture
        .into_inner()
        .map_err(|_| io::Error::other("native search capture poisoned"))
}

pub(crate) fn search_block(
    config: &BlastpConfig,
    database: Vec<fasta::ByteFastaRecord>,
    queries: Vec<fasta::ByteFastaRecord>,
    translated: Option<Block>,
    database_letters: u64,
    queries_prepared: bool,
    masker: Option<&crate::masking::tantan::TantanMasker>,
    ranking: Option<NativeRankingState>,
) -> io::Result<NativeSearchCapture> {
    let capture = std::sync::Mutex::new(NativeSearchCapture {
        ranking,
        ..Default::default()
    });
    run_impl(
        config,
        Some(database),
        if config.tuning.search_block.self_search
            && !config.tuning.multiprocessing.enabled
            && ((!config.tuning.setup.lin_stage1_query
                && config.tuning.search_block.reference == 0)
                || (config.tuning.setup.lin_stage1_query
                    && config.tuning.search_block.reference == config.tuning.search_block.query))
        {
            None
        } else {
            Some(queries)
        },
        None,
        None,
        MaskingAlgo::None,
        false,
        Some(database_letters),
        config.tuning.setup.lin_stage1_query,
        config.round_linear_target,
        config.tuning.search_block.self_search,
        if config.query_cover >= 50.0
            && config.query_cover == config.subject_cover
            && config.translated_query_layout.is_none()
        {
            (config.query_cover / 100.0 - 0.05).max(0.0)
        } else {
            0.0
        },
        false,
        None,
        masker,
        None,
        translated,
        Some(&capture),
        queries_prepared,
    )?;
    capture
        .into_inner()
        .map_err(|_| io::Error::other("native search capture poisoned"))
}

pub(crate) fn validate_added_self_targets(
    matches: &[Match],
    target_count: usize,
) -> io::Result<()> {
    if let Some(target) = matches.iter().find(|target| {
        target.is_synthetic_self_alignment() && target.target_block_id as usize >= target_count
    }) {
        return Err(io::Error::other(format!(
            "--add-self-aln synthetic target index {} is outside the loaded reference block ({} sequences); the source performs an out-of-bounds metadata lookup at this boundary.",
            target.target_block_id, target_count
        )));
    }
    Ok(())
}

// The translated serializer abandons pending bytes on error and is finalized
// explicitly on success, including the gzip trailer.
enum NativeBlockOutput {
    File(crate::util::io::OutputFile),
    Discard(io::Sink),
    Consumer(Box<dyn Write + Send>),
}
impl NativeBlockOutput {
    fn new(config: &BlastpConfig, daa: bool) -> io::Result<Self> {
        if daa && config.output.is_none() {
            return Err(io::Error::other("DAA output requires -a/--daa or -o/--out"));
        }
        let compressor = match config.compression.as_slice() {
            b"" | b"0" => crate::util::io::Compressor::None,
            b"1" => crate::util::io::Compressor::Zlib,
            b"zstd" if cfg!(feature = "zstd") => crate::util::io::Compressor::Zstd,
            b"zstd" => {
                // Source constructs FileSink before discovering this build lacks
                // the optional codec, so preserve the empty file/error ordering.
                let _file = crate::util::io::OutputFile::with_source_profile_path(
                    config.output.as_deref().unwrap_or_else(|| Path::new("")),
                    crate::util::io::Compressor::None,
                    "wb",
                    config.file_buffer_size,
                )
                .map_err(io::Error::other)?;
                return Err(io::Error::other(
                    "Executable was not compiled with ZStd support.",
                ));
            }
            value => {
                let mut message = b"Invalid compression algorithm: ".to_vec();
                message.extend_from_slice(value);
                return Err(crate::util::io::RawDiagnosticError(message).into_io());
            }
        };
        crate::util::io::OutputFile::with_source_profile_path(
            config.output.as_deref().unwrap_or_else(|| Path::new("")),
            compressor,
            "wb",
            config.file_buffer_size,
        )
        .map(Self::File)
        .map_err(io::Error::other)
    }
    fn close(&mut self) -> io::Result<()> {
        match self {
            Self::File(file) => file.close().map_err(io::Error::other),
            Self::Discard(sink) => sink.flush(),
            Self::Consumer(sink) => sink.flush(),
        }
    }
}
impl Write for NativeBlockOutput {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        match self {
            Self::File(file) => {
                file.write_raw(bytes).map_err(io::Error::other)?;
                Ok(bytes.len())
            }
            Self::Discard(sink) => sink.write(bytes),
            Self::Consumer(sink) => sink.write(bytes),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::File(file) => file.flush().map_err(io::Error::other),
            Self::Discard(sink) => sink.flush(),
            Self::Consumer(sink) => sink.flush(),
        }
    }
}
impl std::io::Seek for NativeBlockOutput {
    fn seek(&mut self, position: std::io::SeekFrom) -> io::Result<u64> {
        let Self::File(file) = self else {
            return Err(io::Error::other("output is not seekable"));
        };
        let (offset, origin) = match position {
            std::io::SeekFrom::Start(value) => (
                i64::try_from(value).map_err(io::Error::other)?,
                std::io::SeekFrom::Start(0),
            ),
            std::io::SeekFrom::Current(value) => (value, std::io::SeekFrom::Current(0)),
            std::io::SeekFrom::End(value) => (value, std::io::SeekFrom::End(0)),
        };
        file.seek(offset, origin).map_err(io::Error::other)?;
        file.tell()
            .map(|value| value as u64)
            .map_err(io::Error::other)
    }
}

/// Output lifetime shared by all query/reference blocks. The intermediate
/// merger supplies decoded HSPs and actual compact reference metadata.
/// Source generate_output's final counters; intermediate emission does not
/// count output, while join_query/join_worker count their own final stream.
pub(crate) fn generated_output_statistics(matches: &[Match]) -> Statistics {
    use crate::basic::statistics::StatValue;
    let mut statistics = Statistics::new();
    statistics.inc(StatValue::Pairwise, matches.len() as i64);
    statistics.inc(
        StatValue::Matches,
        matches.iter().map(|target| target.hsps.len() as i64).sum(),
    );
    statistics.inc(
        StatValue::Aligned,
        i64::from(matches.iter().any(|target| !target.hsps.is_empty())),
    );
    statistics
}

pub(crate) struct NativeSearchWriter {
    config: BlastpConfig,
    query_mode: i32,
    output: NativeBlockOutput,
    fields: Vec<FieldId>,
    context_format: Option<SearchFormat>,
    tabular: Option<crate::output::blast_tab_format::ByteBlastTabFormat>,
    daa: bool,
    refs: Vec<(Vec<u8>, u32)>,
    ref_ids: std::collections::HashMap<u64, u32>,
    sequence_count: u64,
    output_database_letters: Option<u64>,
    aligned_queries: u64,
    record_statistics: bool,
}

/// Immutable per-query formatting state shared by blocked-join workers. DAA
/// dictionary entries are reserved at the reference-block boundary before this
/// view is borrowed; workers never mutate or clone the reference dictionary.
pub(crate) struct NativeQueryRenderer<'a> {
    config: &'a BlastpConfig,
    query_mode: i32,
    context_format: Option<SearchFormat>,
    tabular: Option<&'a crate::output::blast_tab_format::ByteBlastTabFormat>,
    daa: bool,
    ref_ids: &'a std::collections::HashMap<u64, u32>,
    database_letters: u64,
}

pub(crate) struct NativeRenderedQuery {
    bytes: Vec<u8>,
    aligned_queries: u64,
    join_task: bool,
}
impl NativeRenderedQuery {
    pub(crate) fn join_task() -> Self {
        Self {
            bytes: Vec::new(),
            aligned_queries: 0,
            join_task: true,
        }
    }
    pub(crate) fn append(&mut self, other: Self) {
        self.bytes.extend_from_slice(&other.bytes);
        self.aligned_queries += other.aligned_queries;
    }
}

/// Sole ordered output owner. Query merging and formatting happen before this
/// commit and do not take its lock or borrow its mutable file stream.
pub(crate) struct NativeJoinedSink<'a> {
    output: &'a mut NativeBlockOutput,
    aligned_queries: &'a mut u64,
    separator: Option<u8>,
    seen_task: bool,
}
impl NativeJoinedSink<'_> {
    pub(crate) fn write(&mut self, rendered: NativeRenderedQuery) -> io::Result<()> {
        if rendered.join_task {
            if self.seen_task {
                if let Some(separator) = self.separator {
                    self.output.write_all(&[separator])?;
                }
            }
            self.seen_task = true;
        }
        self.output.write_all(&rendered.bytes)?;
        *self.aligned_queries += rendered.aligned_queries;
        Ok(())
    }
}

impl NativeQueryRenderer<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn render_query_bytes<S: super::search_output::SearchSubjects + ?Sized>(
        &self,
        query_id: u32,
        oid: u64,
        title: &[u8],
        source: &[Letter],
        frames: &[&[Letter]],
        mate: Option<&[Letter]>,
        quality: &[u8],
        query_self_aln_score: Option<f64>,
        matches: Vec<Match>,
        subjects: &S,
        matrix: &ScoreMatrix,
        aligned: bool,
    ) -> io::Result<NativeRenderedQuery> {
        if matches.iter().any(|target| target.hsps.is_empty()) {
            return Err(io::Error::other("generate_output: target with no hsps."));
        }
        let dna_frames = [source, &[][..], &[][..], &[][..], &[][..], &[][..]];
        let frames = if self.query_mode == 4 {
            &dna_frames[..]
        } else {
            frames
        };
        let mut bytes = Vec::new();
        let mut aligned_queries = 0;
        if self.daa {
            if aligned {
                let begin = write_daa_query_record_bytes(
                    &mut bytes,
                    title,
                    source,
                    if self.query_mode == 2 {
                        SequenceType::AminoAcid
                    } else {
                        SequenceType::Nucleotide
                    },
                );
                for target in matches {
                    let original = subjects.original_oid(target.target_block_id as usize);
                    let dict = self.ref_ids.get(&original).copied().ok_or_else(|| {
                        io::Error::other("Missing reserved DAA reference dictionary entry")
                    })?;
                    for hsp in &target.hsps {
                        write_daa_record_hsp(&mut bytes, hsp, dict);
                    }
                }
                finish_daa_query_record(&mut bytes, begin);
                aligned_queries = 1;
            }
        } else if let Some(format) = self.context_format {
            bytes = format
                .query_with_query_mode_state(
                    query_id,
                    oid,
                    title,
                    source.len() as i32,
                    frames,
                    matches,
                    subjects,
                    matrix,
                    self.query_mode,
                    self.database_letters,
                    self.config.report_unaligned,
                    &self.config.sam,
                    &self.config.tuning.output,
                    query_self_aln_score,
                    aligned,
                )?
                .0;
        } else {
            super::search_output::write_tabular_query_state(
                &mut bytes,
                self.tabular.expect("tabular writer"),
                query_id,
                oid,
                title,
                source,
                frames,
                mate,
                quality,
                query_self_aln_score,
                matches,
                subjects,
                aligned,
            )?;
        }
        Ok(NativeRenderedQuery {
            bytes,
            aligned_queries,
            join_task: false,
        })
    }
}

impl NativeSearchWriter {
    pub(crate) fn joined_parts(
        &mut self,
        matrix: &ScoreMatrix,
    ) -> (NativeQueryRenderer<'_>, NativeJoinedSink<'_>) {
        let database_letters = self.database_letters(matrix);
        let separator = self
            .tabular
            .as_ref()
            .map(|format| format.format.query_separator)
            .filter(|separator| *separator != '\0')
            .map(|separator| separator as u8);
        (
            NativeQueryRenderer {
                config: &self.config,
                query_mode: self.query_mode,
                context_format: self.context_format,
                tabular: self.tabular.as_ref(),
                daa: self.daa,
                ref_ids: &self.ref_ids,
                database_letters,
            },
            NativeJoinedSink {
                output: &mut self.output,
                aligned_queries: &mut self.aligned_queries,
                separator,
                seen_task: false,
            },
        )
    }

    pub(crate) fn set_database_letters(&mut self, full_letters: u64) {
        self.output_database_letters = Some(full_letters);
    }

    pub(crate) fn database_letters(&self, matrix: &ScoreMatrix) -> u64 {
        self.output_database_letters
            .unwrap_or_else(|| matrix.db_letters())
    }

    pub(crate) fn write_ranked_bytes(&mut self, bytes: &[u8]) -> io::Result<()> {
        self.output.write_all(bytes)
    }

    pub(crate) fn new(config: &BlastpConfig, sequence_count: u64) -> io::Result<Self> {
        Self::new_with_query_mode(
            config,
            sequence_count,
            if config.translated_query_layout.is_some() {
                3
            } else {
                2
            },
        )
    }

    pub(crate) fn new_with_query_mode(
        config: &BlastpConfig,
        sequence_count: u64,
        query_mode: i32,
    ) -> io::Result<Self> {
        Self::new_using_output(config, sequence_count, query_mode, None)
    }

    pub(crate) fn new_with_consumer(
        config: &BlastpConfig,
        sequence_count: u64,
        consumer: Box<dyn Write + Send>,
    ) -> io::Result<Self> {
        Self::new_using_output(
            config,
            sequence_count,
            2,
            Some(NativeBlockOutput::Consumer(consumer)),
        )
    }

    fn new_using_output(
        config: &BlastpConfig,
        sequence_count: u64,
        query_mode: i32,
        output: Option<NativeBlockOutput>,
    ) -> io::Result<Self> {
        if !matches!(query_mode, 2 | 3 | 4) {
            return Err(io::Error::other("Invalid native query mode"));
        }
        let context_format = SearchFormat::from_config(config);
        let daa = config
            .outfmt
            .first()
            .is_some_and(|name| name.as_slice() == b"100" || name.as_slice() == b"daa");
        let fields = if daa || context_format.is_some() {
            Vec::new()
        } else if config.outfmt.len() <= 1 {
            format::DEFAULT_TABULAR_FIELDS.to_vec()
        } else {
            config
                .outfmt
                .iter()
                .skip(1)
                .map(|name| {
                    FieldId::from_name_bytes(name).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            crate::util::io::RawDiagnosticError(
                                [b"Invalid output field: ".as_slice(), name.as_slice()].concat(),
                            ),
                        )
                    })
                })
                .collect::<io::Result<Vec<_>>>()?
        };
        let mut tabular = if !daa && context_format.is_none() {
            Some(super::search_output::tabular_format(config)?)
        } else {
            None
        };
        if let Some(format) = &mut tabular {
            format.config.nucleotide_query = query_mode == 4;
        }
        let mut output = match output {
            Some(output) => output,
            None => NativeBlockOutput::new(config, daa)?,
        };
        if daa {
            init_daa(&mut output)?;
        }
        Ok(Self {
            config: config.clone(),
            query_mode,
            output,
            fields,
            context_format,
            tabular,
            daa,
            refs: Vec::new(),
            ref_ids: std::collections::HashMap::new(),
            sequence_count,
            output_database_letters: None,
            aligned_queries: 0,
            record_statistics: true,
        })
    }

    pub(crate) fn needs_reference_sequences(&self) -> bool {
        self.fields
            .iter()
            .any(|field| matches!(field, FieldId::FullSSeq | FieldId::ApproxPIdent))
    }

    pub(crate) fn intermediate_format(&self) -> crate::output::intermediate::OutputFormat {
        let values = if self.config.frame_shift != 0 {
            HspValues::TRANSCRIPT
        } else if let Some(format) = self.context_format {
            format.hsp_values()
        } else if self.daa {
            HspValues::TRANSCRIPT
        } else {
            self.fields
                .iter()
                .fold(HspValues::NONE, |values, field| values | field.hsp_values())
        };
        let values = values
            | self
                .config
                .tuning
                .output
                .clustering_filter
                .as_ref()
                .map_or(HspValues::NONE, |filter| filter.format.hsp_values);
        crate::output::intermediate::OutputFormat::new(
            if self.daa {
                format::FormatCode::Daa
            } else if self.context_format == Some(SearchFormat::Sam) {
                format::FormatCode::Sam
            } else if self.context_format == Some(SearchFormat::Edge) {
                format::FormatCode::Edge
            } else if matches!(
                self.context_format,
                Some(SearchFormat::Bin1 | SearchFormat::Clustering)
            ) {
                format::FormatCode::Bin1
            } else {
                if self
                    .tabular
                    .as_ref()
                    .is_some_and(|format| format.config.json)
                {
                    format::FormatCode::JsonFlat
                } else {
                    format::FormatCode::Tabular
                }
            },
            values,
        )
    }

    pub(crate) fn header_bytes(
        &mut self,
        matrix: &ScoreMatrix,
        title: &[u8],
        len: usize,
    ) -> io::Result<()> {
        if let Some(format) = self.context_format {
            format.header_with_query_mode(
                &mut self.output,
                &self.config,
                matrix,
                title,
                len as i32,
                self.query_mode,
            )?;
        } else if let Some(format) = &self.tabular {
            self.output
                .write_all(&format.header_bytes().map_err(io::Error::other)?)?;
        }
        Ok(())
    }

    /// DAA assigns dictionary IDs while producing each reference block's
    /// intermediate stream, before global culling and query merging.
    pub(crate) fn reserve_archive_subject_bytes(&mut self, oid: u64, title: &[u8], len: usize) {
        if self.daa && !self.ref_ids.contains_key(&oid) {
            let id = self.refs.len() as u32;
            self.refs.push((
                self.config.archive.dictionary_title_bytes(title),
                len as u32,
            ));
            self.ref_ids.insert(oid, id);
        }
    }

    pub(crate) fn write_joined_query_bytes<S: super::search_output::SearchSubjects + ?Sized>(
        &mut self,
        query_id: u32,
        oid: u64,
        title: &[u8],
        source: &[Letter],
        frames: &[&[Letter]],
        mate: Option<&[Letter]>,
        quality: &[u8],
        query_self_aln_score: Option<f64>,
        matches: Vec<Match>,
        subjects: &S,
        matrix: &ScoreMatrix,
    ) -> io::Result<()> {
        let previous = std::mem::replace(&mut self.record_statistics, false);
        let result = self.write_query_bytes(
            query_id,
            oid,
            title,
            source,
            frames,
            mate,
            quality,
            query_self_aln_score,
            matches,
            subjects,
            matrix,
        );
        self.record_statistics = previous;
        result
    }

    pub(crate) fn write_query_bytes<S: super::search_output::SearchSubjects + ?Sized>(
        &mut self,
        query_id: u32,
        oid: u64,
        title: &[u8],
        source: &[Letter],
        frames: &[&[Letter]],
        mate: Option<&[Letter]>,
        quality: &[u8],
        query_self_aln_score: Option<f64>,
        matches: Vec<Match>,
        subjects: &S,
        matrix: &ScoreMatrix,
    ) -> io::Result<()> {
        let output_statistics = generated_output_statistics(&matches);
        let aligned = !matches.is_empty();
        if self.daa && aligned {
            for target in &matches {
                let id = target.target_block_id as usize;
                self.reserve_archive_subject_bytes(
                    subjects.original_oid(id),
                    subjects.title_bytes(id)?,
                    subjects.length(id),
                );
            }
        }
        let (renderer, mut sink) = self.joined_parts(matrix);
        let rendered = renderer.render_query_bytes(
            query_id,
            oid,
            title,
            source,
            frames,
            mate,
            quality,
            query_self_aln_score,
            matches,
            subjects,
            matrix,
            aligned,
        )?;
        sink.write(rendered)?;
        if self.record_statistics {
            crate::basic::statistics::merge_statistics(&output_statistics);
        }
        Ok(())
    }

    pub(crate) fn finish(mut self, matrix: &ScoreMatrix) -> io::Result<()> {
        if self.daa {
            let dictionary = NativeDaaDictionary {
                sequence_count: self.sequence_count,
                refs: self.refs,
            };
            finish_daa_from_sequence_file_bytes(
                &mut self.output,
                &dictionary,
                &DaaRunMetadata {
                    db_letters: matrix.db_letters(),
                    gap_open: matrix.gap_open(),
                    gap_extend: matrix.gap_extend(),
                    reward: 0,
                    penalty: 0,
                    k: matrix.k(),
                    lambda: matrix.lambda(),
                    max_evalue: self.config.max_evalue,
                    matrix: if self.query_mode == 4 {
                        // Source DAA stores config.matrix, not DNA match_reward.
                        self.config.matrix.to_ascii_lowercase()
                    } else {
                        matrix.name().as_bytes().to_vec()
                    },
                    mode: self.query_mode as u32,
                    aligned_queries: self.aligned_queries,
                },
            )?;
        } else if let Some(format) = self.context_format {
            format.footer(&mut self.output)?;
        } else if let Some(format) = &self.tabular {
            self.output.write_all(format.footer().as_bytes())?;
        }
        self.output.close()
    }

    /// MPP query chunks contain only joined query records; Source emits the
    /// command header/footer through its separate master output stream.
    pub(crate) fn finish_multiprocessing_chunk(mut self) -> io::Result<()> {
        self.output.close()
    }
}

enum QueryOutput {
    Skip,
    Tabular { bytes: Vec<u8>, rows: u64 },
    Daa(Vec<Match>),
}
impl crate::align::native::QueryOutput for QueryOutput {
    fn buffered_size(&self) -> usize {
        QueryOutput::buffered_size(self)
    }
}

impl QueryOutput {
    fn buffered_size(&self) -> usize {
        match self {
            Self::Skip => 0,
            Self::Tabular { bytes, .. } => bytes.capacity(),
            Self::Daa(matches) => {
                matches.capacity() * std::mem::size_of::<Match>()
                    + matches
                        .iter()
                        .map(|target| {
                            target.hsps.capacity() * std::mem::size_of::<Hsp>()
                                + target
                                    .hsps
                                    .iter()
                                    .map(|hsp| hsp.transcript.raw_length())
                                    .sum::<usize>()
                        })
                        .sum::<usize>()
            }
        }
    }
}

struct NativeDaaDictionary {
    sequence_count: u64,
    refs: Vec<(Vec<u8>, u32)>,
}

struct NativeDatabaseMetadata {
    subject_taxids: Vec<Vec<TaxId>>,
    taxonomy: Option<TaxonomyTree>,
}

struct NativeBlockSubjects<'a> {
    block: &'a Block,
    metadata: Option<&'a NativeDatabaseMetadata>,
    self_aln_scores: &'a [f64],
    sequences: Option<&'a [std::sync::Arc<[Letter]>]>,
    oids: Option<&'a [u64]>,
}
impl super::search_output::SearchSubjects for NativeBlockSubjects<'_> {
    fn title_bytes(&self, id: usize) -> io::Result<&[u8]> {
        super::search_output::SearchSubjects::title_bytes(self.block, id)
    }
    fn length(&self, id: usize) -> usize {
        self.block.seqs().length(id) as usize
    }
    fn original_oid(&self, id: usize) -> u64 {
        self.oids
            .and_then(|oids| oids.get(id))
            .copied()
            .unwrap_or_else(|| self.block.block_id2oid(id as u32))
    }
    fn sequence_count(&self) -> u64 {
        self.block.seqs().len() as u64
    }
    fn self_aln_score(&self, id: usize) -> Option<f64> {
        self.self_aln_scores.get(id).copied()
    }
    fn sequence(&self, id: usize) -> Option<&[Letter]> {
        Some(
            self.sequences
                .map_or_else(|| self.block.seqs().get(id), |sequences| &sequences[id]),
        )
    }
    fn taxids(&self, id: usize) -> &[TaxId] {
        self.metadata
            .and_then(|metadata| metadata.subject_taxids.get(id))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }
    fn taxonomy(&self) -> Option<&TaxonomyTree> {
        self.metadata
            .and_then(|metadata| metadata.taxonomy.as_ref())
    }
}

fn is_blast_database(path: &Path) -> bool {
    let mut pin = path.as_os_str().to_owned();
    pin.push(".pin");
    let mut pal = path.as_os_str().to_owned();
    pal.push(".pal");
    path.extension().is_some_and(|extension| extension == "pal")
        || crate::util::io::source_c_path(&PathBuf::from(pin)).exists()
        || crate::util::io::source_c_path(&PathBuf::from(pal)).exists()
}

pub(crate) fn requested_search_taxonomy_flags(config: &BlastpConfig) -> SequenceFileFlags {
    let mut flags = requested_taxonomy_flags(&config.outfmt);
    if config.taxon_k != 0 || !config.taxonlist.is_empty() || !config.taxon_exclude.is_empty() {
        flags |= SequenceFileFlags::TAXON_MAPPING | SequenceFileFlags::TAXON_NODES;
    }
    if config.include_lineage
        && config
            .outfmt
            .first()
            .is_some_and(|name| name.as_slice() == b"102")
    {
        flags | SequenceFileFlags::TAXON_SCIENTIFIC_NAMES
    } else {
        flags
    }
}

pub(crate) fn requested_taxonomy_flags(outfmt: &[Vec<u8>]) -> SequenceFileFlags {
    if outfmt.first().is_some_and(|name| name.as_slice() == b"102") {
        return SequenceFileFlags::TAXON_MAPPING | SequenceFileFlags::TAXON_NODES;
    }
    if outfmt
        .first()
        .is_some_and(|name| SearchFormat::from_name_bytes(name).is_some())
    {
        return SequenceFileFlags::NONE;
    }
    let mut flags = SequenceFileFlags::NONE;
    for field in outfmt
        .iter()
        .skip(1)
        .filter_map(|name| FieldId::from_name_bytes(name))
    {
        match field {
            field if field.taxonomy_rank().is_some() => {
                flags |= SequenceFileFlags::TAXON_MAPPING
                    | SequenceFileFlags::TAXON_NODES
                    | SequenceFileFlags::TAXON_RANKS
                    | SequenceFileFlags::TAXON_SCIENTIFIC_NAMES;
            }
            FieldId::STaxIds => flags |= SequenceFileFlags::TAXON_MAPPING,
            FieldId::SSciNames => {
                flags |=
                    SequenceFileFlags::TAXON_MAPPING | SequenceFileFlags::TAXON_SCIENTIFIC_NAMES;
            }
            FieldId::SLineages => {
                flags |= SequenceFileFlags::TAXON_MAPPING
                    | SequenceFileFlags::TAXON_NODES
                    | SequenceFileFlags::TAXON_SCIENTIFIC_NAMES;
            }
            FieldId::SSKingdoms | FieldId::SKingdoms | FieldId::SPhylums => {
                flags |= SequenceFileFlags::TAXON_MAPPING
                    | SequenceFileFlags::TAXON_NODES
                    | SequenceFileFlags::TAXON_RANKS
                    | SequenceFileFlags::TAXON_SCIENTIFIC_NAMES;
            }
            _ => {}
        }
    }
    flags
}

fn load_blast_database(
    path: &Path,
    taxonomy_flags: SequenceFileFlags,
) -> io::Result<(Vec<fasta::ByteFastaRecord>, NativeDatabaseMetadata)> {
    let flags = SequenceFileFlags::SEQS | SequenceFileFlags::TITLES | taxonomy_flags;
    let mut database = BlastDB::new_path(path, flags)?;
    let sequence_count = database.sequence_count() as usize;
    let mut records = Vec::with_capacity(sequence_count);
    let mut subject_taxids = Vec::with_capacity(sequence_count);
    for oid in 0..sequence_count {
        let mut sequence = Vec::new();
        database.source_io(|database| database.seq_data(oid, &mut sequence))?;
        let id = database.source_io(|database| database.seqid_bytes(oid as u64, true, true))?;
        records.push(fasta::ByteFastaRecord {
            quality: None,
            id,
            sequence,
        });
        subject_taxids.push(database.taxids(oid));
    }

    let needs_tree = taxonomy_flags.contains(SequenceFileFlags::TAXON_SCIENTIFIC_NAMES)
        || taxonomy_flags.contains(SequenceFileFlags::TAXON_NODES)
        || taxonomy_flags.contains(SequenceFileFlags::TAXON_RANKS);
    let taxonomy = if needs_tree {
        let mut tree = TaxonomyTree::new();
        let mut seen = std::collections::BTreeSet::new();
        for taxids in &subject_taxids {
            for &assigned_taxid in taxids {
                let mut taxid = assigned_taxid;
                while taxid > 0 && seen.insert(taxid) {
                    let parent = if taxonomy_flags.contains(SequenceFileFlags::TAXON_NODES) {
                        database.get_parent(taxid)
                    } else {
                        0
                    };
                    let rank = if taxonomy_flags.contains(SequenceFileFlags::TAXON_RANKS) {
                        usize::try_from(database.rank(taxid))
                            .ok()
                            .and_then(Rank::from_index)
                            .map_or_else(String::new, |rank| rank.name().to_owned())
                    } else {
                        String::new()
                    };
                    let name = if taxonomy_flags.contains(SequenceFileFlags::TAXON_SCIENTIFIC_NAMES)
                    {
                        database.taxon_scientific_name_bytes(taxid)
                    } else {
                        taxid.to_string().into_bytes()
                    };
                    tree.add_node_bytes(TaxonomyNode {
                        taxid,
                        parent,
                        rank,
                        name,
                    });
                    if !taxonomy_flags.contains(SequenceFileFlags::TAXON_NODES)
                        || taxid == 1
                        || parent == taxid
                    {
                        break;
                    }
                    taxid = parent;
                }
            }
        }
        Some(tree)
    } else {
        None
    };
    Ok((
        records,
        NativeDatabaseMetadata {
            subject_taxids,
            taxonomy,
        },
    ))
}

pub(crate) fn open_search_dmnd(
    path: &Path,
    taxonomy_flags: SequenceFileFlags,
) -> io::Result<DatabaseFile> {
    open_search_dmnd_with_buffer(path, taxonomy_flags, None)
}

pub(crate) fn open_search_dmnd_with_buffer(
    path: &Path,
    taxonomy_flags: SequenceFileFlags,
    file_buffer_size: Option<usize>,
) -> io::Result<DatabaseFile> {
    let needs_mapping = taxonomy_flags.contains(SequenceFileFlags::TAXON_MAPPING);
    let needs_nodes = taxonomy_flags.contains(SequenceFileFlags::TAXON_NODES)
        || taxonomy_flags.contains(SequenceFileFlags::TAXON_RANKS);
    let needs_names = taxonomy_flags.contains(SequenceFileFlags::TAXON_SCIENTIFIC_NAMES);
    let database = DatabaseFile::open_source_with_options_io(
        path,
        DatabaseOpenOptions {
            file_buffer_size,
            load_taxon_mapping: needs_mapping,
            load_taxon_nodes: needs_nodes,
            load_taxon_scientific_names: needs_names,
            ..DatabaseOpenOptions::default()
        },
        taxonomy_flags.contains(SequenceFileFlags::TAXON_RANKS),
    )?;

    Ok(database)
}

/// Add only assigned nodes and their ancestors to the actual renderer tree.
/// The native database still validates/loads the requested components at open.
pub(crate) fn extend_dmnd_taxonomy(
    database: &DatabaseFile,
    taxonomy_flags: SequenceFileFlags,
    taxids: &[TaxId],
    tree: &mut TaxonomyTree,
    seen: &mut std::collections::BTreeSet<TaxId>,
) -> io::Result<()> {
    let needs_nodes = taxonomy_flags.contains(SequenceFileFlags::TAXON_NODES)
        || taxonomy_flags.contains(SequenceFileFlags::TAXON_RANKS);
    let needs_names = taxonomy_flags.contains(SequenceFileFlags::TAXON_SCIENTIFIC_NAMES);
    for &assigned_taxid in taxids {
        let mut taxid = assigned_taxid;
        while taxid > 0 && seen.insert(taxid) {
            let parent = if needs_nodes {
                database.get_parent(taxid).map_err(io::Error::other)?
            } else {
                0
            };
            let rank = if taxonomy_flags.contains(SequenceFileFlags::TAXON_RANKS) {
                usize::try_from(database.rank(taxid).map_err(io::Error::other)?)
                    .ok()
                    .and_then(Rank::from_index)
                    .map_or_else(String::new, |rank| rank.name().to_owned())
            } else {
                String::new()
            };
            let name = if needs_names {
                database.taxon_scientific_name_bytes(taxid)
            } else {
                taxid.to_string().into_bytes()
            };
            tree.add_node_bytes(TaxonomyNode {
                taxid,
                parent,
                rank,
                name,
            });
            if !needs_nodes || taxid == 1 || parent == taxid {
                break;
            }
            taxid = parent;
        }
    }
    Ok(())
}

/// Open a native DIAMOND database with exactly the taxonomy components the
/// selected output fields require.  Upstream `DatabaseFile` validates and
/// loads these components at construction time; doing the same here avoids
/// both silently blank taxonomy columns and loading the usually-large nodes
/// and names arrays for ordinary searches.
fn load_dmnd_database(
    path: &Path,
    taxonomy_flags: SequenceFileFlags,
) -> io::Result<(
    crate::data::dmnd::ReferenceHeader,
    Vec<fasta::ByteFastaRecord>,
    NativeDatabaseMetadata,
)> {
    let needs_mapping = taxonomy_flags.contains(SequenceFileFlags::TAXON_MAPPING);
    let needs_nodes = taxonomy_flags.contains(SequenceFileFlags::TAXON_NODES)
        || taxonomy_flags.contains(SequenceFileFlags::TAXON_RANKS);
    let needs_names = taxonomy_flags.contains(SequenceFileFlags::TAXON_SCIENTIFIC_NAMES);
    let mut database = open_search_dmnd(path, taxonomy_flags)?;

    let header = database.ref_header.clone();
    let sequence_count = database.sequence_count() as usize;
    let mut records = Vec::with_capacity(sequence_count);
    let mut subject_taxids = Vec::with_capacity(sequence_count);
    database.init_seq_access();
    for oid in 0..sequence_count {
        let (sequence, id) = database
            .read_seq_bytes()
            .map_err(io::Error::other)?
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Unexpected end of DIAMOND database",
                )
            })?;
        records.push(fasta::ByteFastaRecord {
            quality: None,
            id,
            sequence,
        });
        subject_taxids.push(if needs_mapping {
            database.taxids(oid).map_err(io::Error::other)?
        } else {
            Vec::new()
        });
    }

    let taxonomy = if needs_nodes || needs_names {
        let mut tree = TaxonomyTree::new();
        let mut seen = std::collections::BTreeSet::new();
        for taxids in &subject_taxids {
            extend_dmnd_taxonomy(&database, taxonomy_flags, taxids, &mut tree, &mut seen)?;
        }
        Some(tree)
    } else {
        None
    };

    Ok((
        header,
        records,
        NativeDatabaseMetadata {
            subject_taxids,
            taxonomy,
        },
    ))
}

impl DaaSequenceFile<Vec<u8>> for NativeDaaDictionary {
    fn sequence_count(&self) -> u64 {
        self.sequence_count
    }

    fn dict_size(&self) -> usize {
        self.refs.len()
    }

    fn dict_title(&self, index: usize) -> Vec<u8> {
        self.refs[index].0.clone()
    }

    fn dict_len(&self, index: usize) -> u32 {
        self.refs[index].1
    }
}

fn encode_daa_query(
    output: &mut Vec<u8>,
    query_name: &[u8],
    query_source: &[Letter],
    input_sequence_type: SequenceType,
    matches: &[Match],
    target_to_dict: &mut [Option<u32>],
    dict_targets: &mut Vec<usize>,
) -> io::Result<u64> {
    let hsp_count = matches
        .iter()
        .map(|target| target.hsps.len() as u64)
        .sum::<u64>();
    if hsp_count == 0 {
        return Ok(0);
    }

    let seek_pos =
        write_daa_query_record_bytes(output, query_name, query_source, input_sequence_type);
    for target in matches {
        let target_id = target.target_block_id as usize;
        let slot = target_to_dict.get_mut(target_id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("DAA target id {target_id} is outside the database block"),
            )
        })?;
        let dict_id = match *slot {
            Some(id) => id,
            None => {
                let id = u32::try_from(dict_targets.len()).map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "DAA dictionary exceeds u32")
                })?;
                *slot = Some(id);
                dict_targets.push(target_id);
                id
            }
        };
        for hsp in &target.hsps {
            write_daa_record_hsp(output, hsp, dict_id);
        }
    }
    finish_daa_query_record(output, seek_pos);
    Ok(hsp_count)
}

#[derive(Default)]
struct ExtensionWorkerScratch {
    /// Flattened six-frame hit list consumed synchronously by extension.
    /// Capacity is retained by the alignment worker across source queries.
    hits: Vec<Hit>,
    query_cbs: [Vec<i8>; 6],
    hauser_values: Vec<f32>,
}

/// Run a simplified blastp search.
///
/// This implements the basic pipeline:
/// 1. Load query and database sequences
/// 2. Extract seeds from both using shapes
/// 3. Find seed matches (hash join)
/// 4. Extend hits with ungapped x-drop
/// 5. Perform gapped Smith-Waterman alignment
/// 6. Filter by e-value and output
pub fn run(config: &BlastpConfig) -> io::Result<()> {
    if let Some(budget) = native_block_budget(config)? {
        return super::search_blocks::run_protein_files(config, budget);
    }
    run_impl(
        config,
        None,
        None,
        None,
        None,
        MaskingAlgo::None,
        false,
        None,
        false,
        false,
        false,
        0.0,
        false,
        None,
        None,
        None,
        None,
        None,
        false,
    )
}

pub(crate) fn native_score_matrix(
    config: &BlastpConfig,
    db_letters: u64,
) -> io::Result<ScoreMatrix> {
    if !config.tuning.matrix_file.as_os_str().is_empty() {
        if config.gap_open == -1 || config.gap_extend == -1 {
            return Err(io::Error::other(
                "Custom scoring matrices require setting the --gapopen and --gapextend options.",
            ));
        }
        if config
            .outfmt
            .first()
            .is_some_and(|format| matches!(format.as_slice(), b"daab" | b"100"))
        {
            return Err(io::Error::other(
                "Custom scoring matrices are not supported for the DAA format.",
            ));
        }
        if config.comp_based_stats as u32 > 1 {
            return Err(io::Error::other("This value for --comp-based-stats is not supported when using a custom scoring matrix."));
        }
    }
    if let Some(matrix) = config.scoring_matrix.as_ref() {
        let mut matrix = matrix.as_ref().clone();
        matrix.set_db_letters(if config.tuning.db_size == 0 {
            db_letters
        } else {
            config.tuning.db_size
        });
        return Ok(matrix);
    }
    native_score_matrix_with_tuning_bytes(
        &config.matrix,
        config.gap_open,
        config.gap_extend,
        config.frame_shift,
        &config.tuning,
        db_letters,
    )
}

pub(crate) fn native_score_matrix_with_tuning(
    matrix: &str,
    gap_open: i32,
    gap_extend: i32,
    frame_shift: i32,
    tuning: &SearchTuning,
    db_letters: u64,
) -> io::Result<ScoreMatrix> {
    native_score_matrix_with_tuning_bytes(
        matrix.as_bytes(),
        gap_open,
        gap_extend,
        frame_shift,
        tuning,
        db_letters,
    )
}

pub(crate) fn native_score_matrix_with_tuning_bytes(
    matrix: &[u8],
    gap_open: i32,
    gap_extend: i32,
    frame_shift: i32,
    tuning: &SearchTuning,
    db_letters: u64,
) -> io::Result<ScoreMatrix> {
    let letters = if tuning.db_size == 0 {
        db_letters
    } else {
        tuning.db_size
    };
    if let Some(matrix) =
        tuning
            .runtime
            .matrix_for(matrix, gap_open, gap_extend, frame_shift, tuning)
    {
        let mut matrix = matrix.as_ref().clone();
        matrix.set_db_letters(letters);
        return Ok(matrix);
    }
    if tuning.matrix_file.as_os_str().is_empty() {
        let matrix = std::str::from_utf8(matrix).map_err(|_| {
            let mut message = b"Unknown scoring matrix: ".to_vec();
            message.extend_from_slice(matrix);
            crate::util::io::RawDiagnosticError(message).into_io()
        })?;
        ScoreMatrix::new_scaled(
            matrix,
            gap_open,
            gap_extend,
            frame_shift,
            tuning.stop_match_score,
            letters,
            tuning.matrix_scale,
        )
        .map_err(io::Error::other)
    } else {
        if gap_open == -1 || gap_extend == -1 {
            return Err(io::Error::other(
                "Custom scoring matrices require setting the --gapopen and --gapextend options.",
            ));
        }
        ScoreMatrix::from_custom_file_path_io(
            &tuning.matrix_file,
            gap_open,
            gap_extend,
            tuning.stop_match_score,
            letters,
        )
    }
}

pub(crate) fn native_tantan_masker(
    config: &BlastpConfig,
    matrix: &ScoreMatrix,
) -> Arc<crate::masking::tantan::TantanMasker> {
    config.tuning.runtime.tantan().unwrap_or_else(|| {
        Arc::new(crate::masking::tantan::TantanMasker::from_score_matrix(
            matrix,
            config.tuning.tantan_min_mask_prob as f32,
        ))
    })
}

pub(crate) fn native_search_setup(
    config: &BlastpConfig,
    matrix: &ScoreMatrix,
    linear_query: bool,
    linear_target: bool,
    freq_masking: bool,
) -> io::Result<sensitivity::SetupSearchResult> {
    let mut input = config.tuning.setup.clone();
    input.threads = if config.threads > 0 {
        config.threads
    } else {
        rayon::current_num_threads().max(1) as i32
    };
    input.contiguous_seed_mode = config.contiguous_seed_mode;
    input.lin_stage1_query = linear_query;
    input.lin_stage1_target = linear_target;
    input.motif_masking = config.motif_masking.clone();
    input.swipe_all = config.swipe_all;
    input.freq_masking |= freq_masking;
    input.frame_shift = config.frame_shift;
    input.global_ranking_targets = config.global_ranking_targets != 0;
    input.approx_min_id = config.approx_min_id;
    input.ext = config
        .extension_mode
        .map_or("", |mode| match mode {
            sensitivity::ExtensionMode::BandedFast => "banded-fast",
            sensitivity::ExtensionMode::BandedSlow => "banded-slow",
            sensitivity::ExtensionMode::Full => "full",
            sensitivity::ExtensionMode::None => "none",
            sensitivity::ExtensionMode::Global => "global",
        })
        .as_bytes()
        .to_vec();
    if input.freq_masking && input.seed_cut != 0.0 {
        return Err(io::Error::other(
            "Incompatible options: --freq-masking, --seed-cut.",
        ));
    }
    if input.freq_sd.is_some_and(|value| value != 0.0) && !input.freq_masking {
        return Err(io::Error::other("--freq-sd requires --freq-masking."));
    }
    if input.minimizer_window.is_some_and(|value| value != 0) && input.contiguous_seed_mode {
        return Err(io::Error::other(
            "Minimizer setting is not compatible with contiguous seed mode.",
        ));
    }
    let errors = crate::util::io::IoErrorOwner::default();
    errors.finish(sensitivity::setup_search_with_io_owner(
        config.sensitivity,
        &input,
        matrix,
        &errors,
    ))
}

pub(crate) fn begin_algorithm_scope(config: &mut BlastpConfig) -> io::Result<()> {
    let requested = config.tuning.algorithm.source_request(
        config.query_indexed,
        config.contiguous_seed_mode,
        config.swipe_all,
    );
    config
        .tuning
        .algorithm
        .begin_source_search(config.sensitivity, requested)
}

/// Source get_output_format settings, applied once before a native search.
pub(crate) fn normalize_output_config(
    config: &BlastpConfig,
) -> io::Result<std::borrow::Cow<'_, BlastpConfig>> {
    normalize_output_config_impl(config, true)
}

#[cfg(feature = "extra")]
pub(crate) fn normalize_output_config_for_dna(
    config: &BlastpConfig,
) -> io::Result<std::borrow::Cow<'_, BlastpConfig>> {
    normalize_output_config_impl(config, false)
}

pub(crate) fn normalize_clustering_filter(
    config: &BlastpConfig,
) -> io::Result<std::borrow::Cow<'_, BlastpConfig>> {
    if config.tuning.output.clustering_filter.is_some() {
        return Ok(std::borrow::Cow::Borrowed(config));
    }
    let clus = cfg!(feature = "mcl")
        && config
            .outfmt
            .first()
            .is_some_and(|name| name.as_slice() == b"clus");
    let filter = if clus {
        if config.outfmt.len() != 2 {
            return Err(io::Error::other(
                "Invalid number of parameters for clustering output format.",
            ));
        }
        Some(
            crate::mcl::clustering_format::ClusteringFilter::new_bytes(
                &config.outfmt[1],
                config.tuning.output.cluster_threshold,
            )
            .map_err(io::Error::other)?,
        )
    } else {
        config
            .tuning
            .output
            .cluster_threshold
            .map(crate::mcl::clustering_format::ClusteringFilter::invalid_output_class)
    };
    if let Some(filter) = filter {
        let mut normalized = config.clone();
        normalized.tuning.output.clustering_filter = Some(Arc::new(filter));
        Ok(std::borrow::Cow::Owned(normalized))
    } else {
        Ok(std::borrow::Cow::Borrowed(config))
    }
}

fn normalize_output_config_impl(
    config: &BlastpConfig,
    protein: bool,
) -> io::Result<std::borrow::Cow<'_, BlastpConfig>> {
    if let std::borrow::Cow::Owned(normalized) = normalize_clustering_filter(config)? {
        return normalize_output_config_impl(&normalized, protein)
            .map(|config| std::borrow::Cow::Owned(config.into_owned()));
    }
    let mp = &config.tuning.multiprocessing;
    if config.global_ranking_targets > 0 && (mp.enabled || mp.initialize || mp.recover) {
        return Err(io::Error::other(
            "Global ranking is not supported in this mode.",
        ));
    }
    if mp.enabled
        && config
            .outfmt
            .first()
            .is_some_and(|format| matches!(format.as_slice(), b"daa" | b"100"))
    {
        return Err(io::Error::other(
            "The DAA format is not supported in multiprocessing mode.",
        ));
    }
    if config.global_ranking_targets > 0 && (config.range_culling || config.taxon_k != 0) {
        return Err(io::Error::other(
            "Global ranking is not supported in this mode.",
        ));
    }
    if config.range_culling && config.frame_shift == 0 {
        return Err(io::Error::other(
            "Query range culling is only supported in frameshift alignment mode (option -F).",
        ));
    }
    if config.range_culling && config.taxon_k != 0 {
        return Err(io::Error::other(
            "--taxon-k is not supported for --range-culling mode.",
        ));
    }
    if protein && config.contiguous_seed_mode && !config.swipe_all {
        if config.target_indexed {
            return Err(io::Error::other("--target-indexed requires --algo 0"));
        }
        if sensitivity::get_traits(config.sensitivity)
            .contiguous_seed
            .is_none()
        {
            return Err(io::Error::other(
                "Contiguous seed mode is not supported for this sensitivity setting.",
            ));
        }
    }
    if protein && config.frame_shift != 0 && config.translated_query_layout.is_none() {
        return Err(io::Error::other(
            "Frameshift alignments are only supported for translated searches.",
        ));
    }
    if config.global_ranking_targets != 0 && config.frame_shift != 0 {
        return Err(io::Error::other(
            "Global ranking mode is not compatible with frameshift alignments.",
        ));
    }
    if protein && config.swipe_all && config.frame_shift > 0 {
        return Err(io::Error::other(
            "Frameshift alignment does not support full matrix extension.",
        ));
    }
    if protein && config.swipe_all {
        sensitivity::soft_masking_algo_bytes(
            &sensitivity::get_traits(config.sensitivity),
            &config.motif_masking,
            true,
            false,
        )
        .map_err(io::Error::other)?;
    }
    if config.global_ranking_targets != 0
        && config
            .outfmt
            .first()
            .is_some_and(|format| matches!(format.as_slice(), b"100" | b"daa"))
    {
        return Err(io::Error::other(
            "The DAA format is not supported in global ranking mode.",
        ));
    }
    preflight_file_buffer_size(config.file_buffer_size)?;
    if !config.compression.is_empty()
        && config
            .outfmt
            .first()
            .is_some_and(|format| matches!(format.as_slice(), b"daa" | b"100"))
    {
        return Err(io::Error::other(
            "Compression is not supported for DAA format.",
        ));
    }
    let mut top = config.toppercent;
    if config
        .outfmt
        .first()
        .is_some_and(|format| format.as_slice() == b"102")
        && top.is_none()
        && config.min_bit_score == 0.0
    {
        top = Some(10.0);
    }
    if top.is_some_and(|value| value < 0.0 || value > 100.0) {
        return Err(io::Error::other(
            "Allowed value range for --top is between 0.0 and 100.0",
        ));
    }
    let max_targets = if top == Some(100.0) {
        top = None;
        i64::MAX
    } else {
        config.max_target_seqs
    };
    if top == config.toppercent
        && max_targets == config.max_target_seqs
        && !(config.swipe_all && (config.query_indexed || config.contiguous_seed_mode))
    {
        Ok(std::borrow::Cow::Borrowed(config))
    } else {
        let mut normalized = config.clone();
        normalized.toppercent = top;
        normalized.max_target_seqs = max_targets;
        if normalized.swipe_all {
            normalized.query_indexed = false;
            normalized.contiguous_seed_mode = false;
        }
        Ok(std::borrow::Cow::Owned(normalized))
    }
}

pub(crate) fn preflight_file_buffer_size(size: usize) -> io::Result<()> {
    Vec::<u8>::new().try_reserve_exact(size).map_err(|_| {
        io::Error::other("Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage.")
    })
}

pub(crate) fn native_block_budget(config: &BlastpConfig) -> io::Result<Option<u64>> {
    // This boundary precedes ReferenceSource/RecordStream allocation in the
    // typed APIs; normalize_output_config runs after their construction.
    preflight_file_buffer_size(config.file_buffer_size)?;
    let mut billions = config.block_size.unwrap_or(0.0);
    if billions == 0.0 {
        billions = sensitivity::get_traits(config.sensitivity).block_size;
    }
    if !billions.is_finite() || billions <= 0.0 || billions * 1e9 < 1.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Invalid block size.",
        ));
    }
    Ok(Some((billions * 1e9) as u64))
}

/// Run the native search with decoded protein query records.
/// Translated callers use `run_with_translated_query_block`, which retains
/// the source DNA, context grouping, IDs and qualities required for extension.
pub fn run_with_query_record_bytes(
    config: &BlastpConfig,
    query_records: Vec<fasta::ByteFastaRecord>,
) -> io::Result<()> {
    if let Some(budget) = native_block_budget(config)? {
        if config.translated_query_layout.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Translated search requires the source-DNA query Block.",
            ));
        }
        return super::search_blocks::run_protein_queries(config, query_records, budget);
    }
    run_impl(
        config,
        None,
        Some(query_records),
        None,
        None,
        MaskingAlgo::None,
        false,
        None,
        false,
        false,
        false,
        0.0,
        false,
        None,
        None,
        None,
        None,
        None,
        false,
    )
}

/// Search decoded protein records using the native seed and extension pipeline.
/// The database name in `config` remains the output-header name; no input file
/// is opened. Each invocation uses its configured worker count, or an existing
/// Rayon pool when the caller explicitly installs one.
pub fn run_with_record_bytes(
    config: &BlastpConfig,
    database: Vec<fasta::ByteFastaRecord>,
    queries: Vec<fasta::ByteFastaRecord>,
) -> io::Result<()> {
    if config.translated_query_layout.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "decoded protein-record search does not accept a translated-query layout",
        ));
    }
    if let Some(budget) = native_block_budget(config)? {
        return super::search_blocks::run_protein_records(config, database, queries, budget);
    }
    run_impl(
        config,
        Some(database),
        Some(queries),
        None,
        None,
        MaskingAlgo::None,
        false,
        None,
        false,
        false,
        false,
        0.0,
        false,
        None,
        None,
        None,
        None,
        None,
        false,
    )
}

/// Run translated search with the upstream-shaped query block as its sole
/// owner of source DNA and six-frame protein contexts.
pub fn run_with_translated_query_block(
    config: &BlastpConfig,
    query_block: Block,
) -> io::Result<()> {
    if let Some(budget) = native_block_budget(config)? {
        return super::search_blocks::run_translated_block(config, query_block, budget);
    }
    run_impl(
        config,
        None,
        None,
        None,
        None,
        MaskingAlgo::None,
        false,
        None,
        false,
        false,
        false,
        0.0,
        false,
        None,
        None,
        None,
        Some(query_block),
        None,
        false,
    )
}

#[derive(Debug, Clone)]
pub struct InMemorySearchEdge {
    pub edge: crate::output::edge::EdgeData,
    pub approx_id: f64,
    pub identity: f64,
}

/// Mirror `Block::remove_soft_masking`: extension consumes the restored query
/// residues (retaining SEED_MASK annotations), not the hard-masked copy used
/// while enumerating seeds.
fn restore_extension_queries<T>(records: &mut [fasta::FastaRecord<T>], restored: &SequenceSet) {
    assert_eq!(records.len(), restored.len());
    for (query_id, record) in records.iter_mut().enumerate() {
        record.sequence.clear();
        record.sequence.extend_from_slice(restored.get(query_id));
    }
}

/// Run the native protein search pipeline on already decoded records and
/// return clustering edge records instead of formatting a file.
pub fn run_edges_in_memory_bytes(
    config: &BlastpConfig,
    database: Vec<fasta::ByteFastaRecord>,
    queries: Option<Vec<fasta::ByteFastaRecord>>,
    approx_min_id: f64,
    linear_stage1_query: bool,
    linear_stage1_target: bool,
    self_search: bool,
    query_or_target_cover: f64,
    mutual_cover: bool,
    hamming_ext: HammingExtConfig,
    soft_tantan_masker: Option<&crate::masking::tantan::TantanMasker>,
    soft_masking_override: MaskingAlgo,
    freq_masking: bool,
    effective_db_letters: Option<u64>,
    extension_mode_override: Option<sensitivity::ExtensionMode>,
) -> io::Result<Vec<InMemorySearchEdge>> {
    let mut source_config = config.clone();
    source_config.tuning.self_search = self_search;
    source_config.tuning.setup.lin_stage1_query |= linear_stage1_query;
    source_config.tuning.setup.lin_stage1_target |= linear_stage1_target;
    let config = &source_config;
    let block_size = config.block_size.or_else(|| {
        (config.global_ranking_targets != 0 || config.iterate.is_some()).then_some(0.0)
    });
    if let Some(mut block_size) = block_size {
        preflight_file_buffer_size(config.file_buffer_size)?;
        if block_size == 0.0 {
            block_size = sensitivity::get_traits(config.sensitivity).block_size;
        }
        let database_letters = effective_db_letters
            .filter(|&letters| letters != 0)
            .unwrap_or_else(|| {
                database
                    .iter()
                    .map(|record| record.sequence.len() as u64)
                    .sum()
            });
        return super::search_blocks::run_cluster_edges(
            config,
            &database,
            queries.as_deref(),
            block_size,
            database_letters,
            ClusterSearchOptions {
                approx_min_id,
                linear_query: linear_stage1_query,
                linear_target: linear_stage1_target,
                self_search,
                query_or_target_cover,
                mutual_cover,
                hamming: hamming_ext,
                masker: soft_tantan_masker,
                soft_masking: soft_masking_override,
                freq_masking,
                extension: extension_mode_override,
            },
        );
    }
    crate::basic::statistics::with_search_statistics_prepared(
        !config.tuning.multiprocessing.enabled
            || (!config.tuning.multiprocessing.initialize
                && !config.tuning.multiprocessing.recover),
        || super::source_config::begin_search_temp(config),
        || {
            let edges = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
            run_impl(
                config,
                Some(database),
                queries,
                Some(edges.clone()),
                Some(approx_min_id),
                soft_masking_override,
                freq_masking,
                effective_db_letters,
                linear_stage1_query,
                linear_stage1_target,
                self_search,
                query_or_target_cover,
                mutual_cover,
                Some(hamming_ext),
                soft_tantan_masker,
                extension_mode_override,
                None,
                None,
                false,
            )?;
            let result = edges.lock().unwrap().clone();
            Ok(result)
        },
    )
}

pub(crate) fn install_search_threads<T: Send>(
    threads: i32,
    work: impl FnOnce() -> io::Result<T> + Send,
) -> io::Result<T> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(0) as usize)
        .build()
        .map_err(io::Error::other)?;
    pool.install(work)
}

fn run_impl(
    config: &BlastpConfig,
    database_records: Option<Vec<fasta::ByteFastaRecord>>,
    query_records_override: Option<Vec<fasta::ByteFastaRecord>>,
    edge_output: Option<std::sync::Arc<std::sync::Mutex<Vec<InMemorySearchEdge>>>>,
    approx_min_id_override: Option<f64>,
    soft_masking_override: MaskingAlgo,
    freq_masking: bool,
    effective_db_letters: Option<u64>,
    linear_stage1_query: bool,
    linear_stage1_target: bool,
    self_search: bool,
    query_or_target_cover: f64,
    mutual_cover: bool,
    hamming_ext: Option<HammingExtConfig>,
    soft_tantan_masker: Option<&crate::masking::tantan::TantanMasker>,
    extension_mode_override: Option<sensitivity::ExtensionMode>,
    mut translated_query_block: Option<Block>,
    capture: Option<&std::sync::Mutex<NativeSearchCapture>>,
    queries_prepared: bool,
) -> io::Result<()> {
    if rayon::current_thread_index().is_none() {
        return install_search_threads(config.threads, move || {
            run_impl(
                config,
                database_records,
                query_records_override,
                edge_output,
                approx_min_id_override,
                soft_masking_override,
                freq_masking,
                effective_db_letters,
                linear_stage1_query,
                linear_stage1_target,
                self_search,
                query_or_target_cover,
                mutual_cover,
                hamming_ext,
                soft_tantan_masker,
                extension_mode_override,
                translated_query_block,
                capture,
                queries_prepared,
            )
        });
    }
    if capture.is_none() && edge_output.is_none() {
        super::source_config::begin_search_temp(config)?;
    }
    let normalized = if capture.is_none() {
        normalize_output_config(config)?
    } else {
        normalize_clustering_filter(config)?
    };
    let config = normalized.as_ref();
    // Source resolves the explicit clustering override before setup validates
    // defaults such as linear full-matrix extension and anchored compatibility.
    let effective_extension_mode = extension_mode_override.or(config.extension_mode);
    let extension_config = if effective_extension_mode != config.extension_mode {
        let mut effective = config.clone();
        effective.extension_mode = effective_extension_mode;
        std::borrow::Cow::Owned(effective)
    } else {
        std::borrow::Cow::Borrowed(config)
    };
    let config = extension_config.as_ref();
    let mut execution_config = config.clone();
    let requested = execution_config.tuning.algorithm.source_request(
        execution_config.query_indexed,
        execution_config.contiguous_seed_mode,
        execution_config.swipe_all,
    );
    execution_config
        .tuning
        .algorithm
        .begin_search(execution_config.sensitivity, requested);
    let selected_algorithm = execution_config.tuning.algorithm.selected()?;
    execution_config.query_indexed = selected_algorithm == crate::run::config::Algo::QueryIndexed;
    execution_config.contiguous_seed_mode = selected_algorithm == crate::run::config::Algo::CtgSeed;
    let config = &execution_config;
    let linear_stage1_query = linear_stage1_query || config.tuning.setup.lin_stage1_query;
    let linear_stage1_target = linear_stage1_target || config.tuning.setup.lin_stage1_target;
    let linear_stage1_combo = config.tuning.setup.lin_stage1_combo;
    if usize::from(linear_stage1_query)
        + usize::from(linear_stage1_target)
        + usize::from(linear_stage1_combo)
        > 1
    {
        return Err(io::Error::other(
            "Multiple linearization options are not allowed to be used together",
        ));
    }
    let mut ranking_state = if config.global_ranking_targets != 0 {
        capture
            .ok_or_else(|| io::Error::other("global ranking requires native block orchestration"))?
            .lock()
            .map_err(|_| io::Error::other("native search capture poisoned"))?
            .ranking
            .take()
    } else {
        None
    };
    let extension_mode_override = extension_mode_override.or(config.extension_mode);
    if config.global_ranking_targets != 0
        && extension_mode_override.is_some_and(|mode| mode != sensitivity::ExtensionMode::Full)
    {
        return Err(io::Error::other(
            "Global ranking only supports full matrix extension.",
        ));
    }
    if config.frame_shift > 0
        && (extension_mode_override == Some(sensitivity::ExtensionMode::Full)
            || (extension_mode_override.is_none()
                && (config.swipe_all
                    || config.global_ranking_targets != 0
                    || ((linear_stage1_query || linear_stage1_target || linear_stage1_combo)
                        && !config.tuning.setup.linclust_banded_ext))))
    {
        return Err(io::Error::other(
            "Frameshift alignment does not support full matrix extension.",
        ));
    }
    if config.swipe_all
        && extension_mode_override.is_some_and(|mode| mode != sensitivity::ExtensionMode::Full)
    {
        return Err(io::Error::other(
            "--swipe only supports full matrix extension.",
        ));
    }
    let approx_min_id = approx_min_id_override.unwrap_or(config.approx_min_id);
    let hamming_ext = hamming_ext.unwrap_or_else(|| {
        native_hamming_config(
            config,
            approx_min_id,
            config.query_or_subject_cover.max(query_or_target_cover),
        )
    });
    let mutual_cover = mutual_cover || config.tuning.mutual_cover.is_some();
    let mut query_seed_indexed = config.query_indexed || config.contiguous_seed_mode;
    let round_skip = config.round_query_skip.as_deref();
    let query_contexts = if config.translated_query_layout.is_some() {
        6
    } else {
        1
    };
    let start = Instant::now();
    if query_seed_indexed && config.target_indexed {
        return Err(io::Error::other("--target-indexed requires --algo 0"));
    }
    let context_format = SearchFormat::from_config(config);
    CbsMode::from_search_code(
        config.comp_based_stats as u32,
        config.translated_query_layout.is_some(),
        config.global_ranking_targets,
    )
    .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?;
    // Upstream linear self-search aliases target to query instead of loading
    // or copying a second block. `run_edges_in_memory` signals that exact mode
    // by supplying the database records without a separate query vector.
    let shared_self_records = database_records.is_some()
        && query_records_override.is_none()
        && (edge_output.is_some() || capture.is_some())
        && config.tuning.self_search;
    // `Search::Config` derives this once from equal query/subject coverage
    // and stage sensitivity. It selects the mutual-coverage stage-1 kernels,
    // not merely length sorting, so carry it through the native seed worker.
    let min_length_ratio = if let Some(ratio) = config.tuning.search_min_length_ratio {
        ratio
    } else {
        crate::run::config::derive_min_length_ratio(
            &crate::run::config::ConfigInput {
                min_length_ratio: config.tuning.min_length_ratio,
                query_cover: config.query_cover,
                subject_cover: config.subject_cover,
                query_translated: config.translated_query_layout.is_some(),
                lin_stage1_query: linear_stage1_query,
                ..Default::default()
            },
            &[crate::run::config::Round::new(
                config.sensitivity,
                linear_stage1_target,
            )],
        )
        .map_err(io::Error::other)?
    };

    // Load database sequences - try DMND first, then FASTA.
    // C++ `auto_append_extension_if_exists` (`config.cpp:770`) APPENDS `.dmnd`
    // only when the file as-given doesn't exist — it never STRIPS an existing
    // extension. So `--db nr.fasta` resolves to `nr.fasta` in C++ (then parsed
    // as FASTA), but `--db nr.fasta` was resolving to `nr.dmnd` in Rust because
    // `with_extension("dmnd")` REPLACES rather than appends. Fix: only consult
    // a `<db>.dmnd` companion when the given path doesn't exist.
    let taxonomy_flags = requested_search_taxonomy_flags(config);
    let taxonomy_requested = taxonomy_flags.contains(SequenceFileFlags::TAXON_MAPPING);
    let mut database_metadata = None;
    let (mut db_records, _db_from_dmnd) = if let Some(records) = database_records {
        if taxonomy_requested && capture.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "taxonomy output is not available for in-memory clustering search",
            ));
        }
        (records, false)
    } else {
        let db_path = Path::new(&config.database);
        let dmnd_path = if db_path.extension().is_some_and(|e| e == "dmnd") {
            db_path.to_path_buf()
        } else if crate::util::io::source_c_path(&db_path).exists() {
            // File-as-given exists (typically a FASTA). Don't try to substitute
            // a `.dmnd` sibling — leave the path alone so the FASTA branch
            // below runs. We set `dmnd_path` to the given path; `exists()` is
            // true but it's not a DMND file, so the DMND read will fail format
            // detection and we fall through to FASTA. To avoid that wasted
            // open, mark with a non-existent suffix instead.
            let mut p = db_path.as_os_str().to_owned();
            p.push(".dmnd");
            std::path::PathBuf::from(p)
        } else {
            // No file as-given — append `.dmnd` (don't replace).
            let mut p = db_path.as_os_str().to_owned();
            p.push(".dmnd");
            std::path::PathBuf::from(p)
        };

        if is_blast_database(&config.database) {
            if config.target_indexed {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--target-indexed is not compatible with BLAST databases.",
                ));
            }
            let (records, metadata) = load_blast_database(&config.database, taxonomy_flags)?;
            search_message!(
                "Database: {} sequences (NCBI BLAST database)",
                records.len()
            );
            database_metadata = Some(metadata);
            (records, false)
        } else if crate::util::io::source_c_path(&dmnd_path).exists() {
            // Keep ordinary search on the lean sequence-only reader.  An
            // empty `Vec<TaxId>` per subject would cost 24 bytes per database
            // sequence even when no taxonomy column was requested.
            let (header, records) = if taxonomy_requested {
                let (header, records, metadata) = load_dmnd_database(&dmnd_path, taxonomy_flags)?;
                database_metadata = Some(metadata);
                (header, records)
            } else {
                crate::data::dmnd_reader::read_dmnd_bytes(&dmnd_path)?
            };
            search_message!(
                "Database: {} sequences, {} letters (DMND)",
                header.sequences,
                header.letters
            );
            (records, true)
        } else {
            if taxonomy_requested {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "taxonomy output requires a database containing taxonomy metadata",
                ));
            }
            let fasta_path = if db_path.extension().is_none() {
                db_path.with_extension("faa")
            } else {
                db_path.to_path_buf()
            };
            let records = fasta::read_fasta_file_bytes(&fasta_path, SequenceType::AminoAcid)?;
            search_message!("Database: {} sequences (FASTA)", records.len());
            (records, false)
        }
    };

    // Compute total database letters for E-value normalization
    let loaded_db_letters: u64 = db_records.iter().map(|r| r.sequence.len() as u64).sum();
    // `Search::run` normalizes E-values with global `config.db_size` when it
    // is set. The low-memory clustering wrapper sets that value to the full
    // input database size before searching any individual super-block.
    let db_letters = effective_db_letters
        .filter(|&letters| letters != 0)
        .unwrap_or(loaded_db_letters);

    let score_matrix = native_score_matrix(config, db_letters)?;
    let setup = native_search_setup(
        config,
        &score_matrix,
        linear_stage1_query,
        linear_stage1_target,
        freq_masking,
    )?;
    let freq_masking = freq_masking || config.tuning.setup.freq_masking;

    // Load query sequences
    let mut query_records = if shared_self_records {
        Vec::new()
    } else if let Some(records) = query_records_override {
        records
    } else {
        let mut records = Vec::new();
        for qf in &config.query_files {
            records.extend(fasta::read_fasta_file_bytes(
                Path::new(qf),
                SequenceType::AminoAcid,
            )?);
        }
        records
    };
    let translated_context_count = translated_query_block
        .as_ref()
        .map_or(0, |block| block.seqs().len());
    search_message!(
        "Queries: {} sequences",
        if let Some(layout) = &config.translated_query_layout {
            layout.mate_indices.len()
        } else if shared_self_records {
            db_records.len()
        } else {
            query_records.len()
        }
    );

    // Edge consumers and structured output address original query OIDs even
    // after length sorting. Retain target permutations for edge consumers.
    let mut edge_query_oids = (edge_output.is_some()
        || context_format.is_some()
        || capture.is_some()
        || config
            .outfmt
            .iter()
            .any(|field| field.as_slice() == b"qnum"))
    .then(|| {
        (0..if shared_self_records {
            db_records.len()
        } else {
            query_records.len()
        })
            .map(|oid| oid as u64)
            .collect::<Vec<_>>()
    });
    let mut edge_target_oids = (edge_output.is_some()
        || context_format == Some(SearchFormat::Edge)
        || capture.is_some()
        || config
            .outfmt
            .iter()
            .any(|field| field.as_slice() == b"snum"))
    .then(|| {
        (0..db_records.len())
            .map(|oid| oid as u64)
            .collect::<Vec<_>>()
    });

    // Upstream derives a minimum length ratio when the two coverage cutoffs
    // are equal and at least 50%, then length-sorts both sides before search
    // (`run/config.cpp` and `run/double_indexed.cpp`).  Besides improving the
    // stage-1 schedule, this determines the observable query order in tabular
    // output.  Retain the original ordinal as the tie breaker: C++ sorts
    // `pair<length, BlockId>` with `greater`, so equal-length records appear
    // in descending input order as well.
    let _ = crate::util::system::log_rss();
    let retain_seed_ids = min_length_ratio != 0.0
        || config.global_ranking_targets != 0
        || (config.tuning.self_search && config.tuning.search_block.reference == 0)
        || linear_stage1_combo;
    let length_sort_queries = min_length_ratio > 0.0
        || (linear_stage1_query && !retain_seed_ids && !config.tuning.kmer_ranking);
    // Iterated rounds linearize targets independently of query coverage.
    // `load_ref_block` sorts that side even when queries retain input order.
    if !shared_self_records && (linear_stage1_target || length_sort_queries) {
        let sort_by_cpp_length_order =
            |records: &mut Vec<fasta::ByteFastaRecord>, oids: &mut Option<Vec<u64>>| {
                let mut indexed: Vec<_> = records.drain(..).enumerate().collect();
                indexed.sort_unstable_by(|(left_id, left), (right_id, right)| {
                    right
                        .sequence
                        .len()
                        .cmp(&left.sequence.len())
                        .then_with(|| right_id.cmp(left_id))
                });
                if let Some(oids) = oids {
                    let previous = std::mem::take(oids);
                    oids.extend(indexed.iter().map(|(oid, _)| previous[*oid]));
                }
                records.extend(indexed.into_iter().map(|(_, mut record)| {
                    record.quality = None;
                    record
                }));
            };
        // The clustering backend has already applied the upstream
        // `Block::length_sort` order before passing an aliased self-search
        // block. Sorting that block again here reverses equal-length groups
        // and invalidates its local-to-database OID mapping.
        // blastx contexts are deliberately consecutive: the spill bins and
        // extension stage consume all six as one source query.
        if length_sort_queries
            && config.translated_query_layout.is_none()
            && !queries_prepared
            && !config.tuning.query_pre_sorted
        {
            sort_by_cpp_length_order(&mut query_records, &mut edge_query_oids);
        }
        if (linear_stage1_target || min_length_ratio > 0.0) && !config.tuning.kmer_ranking {
            if let Some(metadata) = database_metadata.as_mut() {
                let mut indexed: Vec<_> = db_records
                    .drain(..)
                    .zip(metadata.subject_taxids.drain(..))
                    .enumerate()
                    .collect();
                indexed.sort_unstable_by(|(left_id, (left, _)), (right_id, (right, _))| {
                    right
                        .sequence
                        .len()
                        .cmp(&left.sequence.len())
                        .then_with(|| right_id.cmp(left_id))
                });
                let previous_oids = edge_target_oids.as_mut().map(std::mem::take);
                for (oid, (record, taxids)) in indexed {
                    if let Some(oids) = edge_target_oids.as_mut() {
                        oids.push(previous_oids.as_ref().unwrap()[oid]);
                    }
                    db_records.push(record);
                    metadata.subject_taxids.push(taxids);
                }
            } else {
                sort_by_cpp_length_order(&mut db_records, &mut edge_target_oids);
            }
        }
    }

    let output_target_sequences = (capture.is_none()
        && edge_output.is_none()
        && config
            .outfmt
            .iter()
            .any(|field| matches!(field.as_slice(), b"full_sseq" | b"approx_pident")))
    .then(|| {
        db_records
            .iter()
            .map(|record| std::sync::Arc::<[Letter]>::from(record.sequence.clone()))
            .collect::<Vec<_>>()
    });

    use crate::basic::value::{MASK_LETTER, SEED_MASK};
    use rayon::iter::IntoParallelRefMutIterator;
    let owned_tantan_masker = (config.masking == MaskingMode::Tantan
        && soft_tantan_masker.is_none())
    .then(|| native_tantan_masker(config, &score_matrix));
    let tantan_masker = soft_tantan_masker.or(owned_tantan_masker.as_deref());
    // AUTO must inspect the actual hard-masked query seed set before choosing
    // eager versus lazy reference masking. Shared self blocks are already the
    // query owner, and therefore always receive query hard masking here.
    let defer_target_masking =
        selected_algorithm == crate::run::config::Algo::Auto && !shared_self_records;
    match config.masking {
        MaskingMode::None => {
            // Protein databases retain makedb's tantan soft-mask bit.  Search
            // mode 0 disables masking, so expose the stored residues exactly
            // as upstream does instead of letting DP treat them as masked.
            db_records
                .par_iter_mut()
                .for_each(|r| crate::masking::remove_bit_mask(&mut r.sequence));
        }
        MaskingMode::Tantan => {
            // C++ constructs the tantan masker from the active ScoreMatrix and
            // applies it at search time to both targets and queries
            // (`Masking::Masking(const ScoreMatrix&)`, then `mask_seqs` in
            // `run/double_indexed.cpp`). This matters for non-default matrices
            // such as BLOSUM45: reusing makedb's default BLOSUM62 mask leaves
            // low-complexity regions under-masked and inflates self scores.
            // Upstream uses the same configured Masking singleton for hard
            // and soft tantan. Cluster adapters can supply its non-default
            // threshold; ordinary file searches retain the default masker.
            let tantan_masker = tantan_masker.expect("configured tantan masker");
            db_records.par_iter_mut().for_each(|r| {
                crate::masking::remove_bit_mask(&mut r.sequence);
                if shared_self_records
                    || (!defer_target_masking && (!query_seed_indexed || config.frame_shift != 0))
                {
                    tantan_masker.mask(&mut r.sequence);
                }
            });
            let mask_query = |record: &mut fasta::ByteFastaRecord| {
                crate::masking::remove_bit_mask(&mut record.sequence);
                tantan_masker.mask(&mut record.sequence);
            };
            if queries_prepared {
                // Query hard masks persist across reference blocks; seed-only
                // annotations are cleared at the end of each seed-search phase.
            } else if shared_self_records {
                // The database mutation above is the query mutation too.
            } else if let Some(block) = translated_query_block.as_mut() {
                for query_id in 0..block.seqs().len() {
                    crate::masking::remove_bit_mask(block.seqs_mut().get_mut(query_id));
                    tantan_masker.mask(block.seqs_mut().get_mut(query_id));
                }
            } else if query_records
                .iter()
                .any(|record| record.sequence.len() > LONG_QUERY_RESIDUES)
            {
                // Tantan's DP workspace is roughly proportional to sequence
                // length. Six long blastx contexts masked concurrently create
                // one large workspace per Rayon worker, whereas upstream's
                // query-block masking reuses one workspace here.
                query_records.iter_mut().for_each(mask_query);
            } else {
                query_records.par_iter_mut().for_each(mask_query);
            }
        }
        MaskingMode::BlastSeg => {
            // Run::Config maps BLAST SEG hard masking to targets only.
            db_records.par_iter_mut().for_each(|record| {
                crate::masking::remove_bit_mask(&mut record.sequence);
                if shared_self_records
                    || (!defer_target_masking && (!query_seed_indexed || config.frame_shift != 0))
                {
                    crate::masking::mask_sequence(&mut record.sequence, MaskingAlgo::Seg);
                }
            });
            if !shared_self_records && !queries_prepared {
                if let Some(block) = translated_query_block.as_mut() {
                    for query_id in 0..block.seqs().len() {
                        crate::masking::remove_bit_mask(block.seqs_mut().get_mut(query_id));
                    }
                } else {
                    query_records
                        .par_iter_mut()
                        .for_each(|record| crate::masking::remove_bit_mask(&mut record.sequence));
                }
            }
        }
    }

    // Convert tantan soft masks (SEED_MASK bit) to hard masks (MASK_LETTER = X).
    // C++ blastp calls `mask_seqs(..., hard_mask=true)` on both queries and
    // targets (run/double_indexed.cpp:127 and :719), which `Masking::operator()`
    // dispatches to `tantan::mask` in mode 1 — REPLACING masked letters with
    // `value_traits.mask_char` (= MASK_LETTER = 23) rather than OR-ing in the
    // high bit. Downstream the SIMD score profile then scores those positions
    // as `BLOSUM62(X, X) = -1`, not 0.
    //
    // Keeping the soft mask here makes our DP zero those positions instead of
    // applying the -1 BLOSUM penalty — for a self-self of a heavily-masked
    // protein (Q8QZQ8: 120+ masked residues), that's 120+ extra score relative
    // to C++ and shifts which targets fit inside `-k 25`.
    let hard_mask = |seq: &mut [Letter]| {
        for l in seq.iter_mut() {
            if *l & SEED_MASK != 0 {
                *l = MASK_LETTER;
            }
        }
    };
    if config.masking == MaskingMode::Tantan {
        if shared_self_records
            || (!defer_target_masking && (!query_seed_indexed || config.frame_shift != 0))
        {
            db_records
                .par_iter_mut()
                .for_each(|r| hard_mask(&mut r.sequence));
        }
        if !shared_self_records && !queries_prepared {
            if let Some(block) = translated_query_block.as_mut() {
                for query_id in 0..block.seqs().len() {
                    hard_mask(block.seqs_mut().get_mut(query_id));
                }
            } else {
                query_records
                    .par_iter_mut()
                    .for_each(|r| hard_mask(&mut r.sequence));
            }
        }
    }

    // Keep compact, alignment-ready sequence storage before motif masking.
    // Build the final contiguous SequenceSet directly: a Vec<Vec<Letter>>
    // intermediary adds one allocation and 24 bytes of metadata per query
    // context, which is material for six-frame translated searches.
    let query_source = if shared_self_records {
        &db_records
    } else {
        &query_records
    };
    let mut stage2_query_set = if let Some(block) = translated_query_block.as_mut() {
        let sequences = block.take_seqs();
        debug_assert_eq!(sequences.len(), translated_context_count);
        sequences
    } else {
        let mut sequences = SequenceSet::new();
        sequences.reserve_capacity(
            query_source.len(),
            query_source
                .iter()
                .map(|record| record.sequence.len())
                .sum(),
        );
        for record in query_source {
            sequences.push(&record.sequence);
        }
        sequences
    };
    query_seed_indexed = super::search_algorithm::resolve_query(
        &mut stage2_query_set,
        translated_query_block.as_mut(),
        config,
        &score_matrix,
        &setup,
        soft_masking_override,
        freq_masking,
        db_letters,
    )?;
    if defer_target_masking && (!query_seed_indexed || config.frame_shift != 0) {
        db_records
            .par_iter_mut()
            .for_each(|record| match config.masking {
                MaskingMode::None => (),
                MaskingMode::Tantan => {
                    tantan_masker
                        .expect("configured tantan masker")
                        .mask(&mut record.sequence);
                    hard_mask(&mut record.sequence);
                }
                MaskingMode::BlastSeg => {
                    crate::masking::mask_sequence(&mut record.sequence, MaskingAlgo::Seg)
                }
            });
    }
    let source_index_chunks = if config.index_chunks != 0 {
        config.index_chunks as usize
    } else if query_seed_indexed {
        1
    } else {
        sensitivity::get_traits(config.sensitivity).index_chunks as usize
    };
    let mut db_block = Block::new();
    db_block.set_max_motif_len(config.tuning.max_motif_len);
    for (idx, record) in db_records.iter().enumerate() {
        let original_oid = if let Some(state) = ranking_state.as_ref() {
            let local = edge_target_oids
                .as_ref()
                .map_or(idx as u64, |oids| oids[idx]);
            let oid = if let Some(oids) = state.reference_oids.as_ref() {
                *oids.get(local as usize).ok_or_else(|| {
                    io::Error::other("Ranked reference OID mapping has the wrong length.")
                })?
            } else {
                state
                    .reference_begin
                    .checked_add(local)
                    .ok_or_else(|| io::Error::other("Ranked reference OID overflow."))?
            };
            u32::try_from(oid).map_err(|_| {
                io::Error::other("Global ranking subject OID exceeds the uint32 range.")
            })?;
            oid
        } else {
            idx as u64
        };
        db_block
            .push_back_bytes(
                &record.sequence,
                Some(&record.id),
                None,
                original_oid,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .map_err(io::Error::other)?;
    }

    let needs_self_aln = config
        .tuning
        .output
        .clustering_filter
        .as_ref()
        .is_some_and(|filter| {
            filter
                .format
                .flags
                .any(crate::output::format::OutputFlags::SELF_ALN_SCORES)
        })
        || config
            .outfmt
            .iter()
            .any(|field| field.as_slice() == b"normalized_bitscore");
    // At this point eager target masking is complete; query-indexed targets
    // remain unmasked until lazy extension. This is the source table boundary.
    let target_self_aln_scores = if needs_self_aln {
        super::search_output::self_alignment_scores_with_window(
            db_block.seqs(),
            &score_matrix,
            config.comp_based_stats,
            config.tuning.cbs_window,
        )
    } else {
        Vec::new()
    };
    let clustering_target_self_scores = config
        .tuning
        .output
        .clustering_filter
        .as_ref()
        .map(|_| Arc::<[f64]>::from(target_self_aln_scores.clone()));
    let query_self_aln_scores = if needs_self_aln
        && (!queries_prepared || config.tuning.output.clustering_filter.is_some())
    {
        super::search_output::self_alignment_scores_with_window(
            &stage2_query_set,
            &score_matrix,
            config.comp_based_stats,
            config.tuning.cbs_window,
        )
    } else {
        Vec::new()
    };

    let mut hit_store = if config.swipe_all {
        // Source exhaustive mode skips query/reference histograms, soft masks,
        // indices, joins, Hamming/ungapped seed filters and hit ingestion.
        // Global ranking remains accepted: its table stays empty because no
        // seed ranking updates occur, then the normal ranking stage consumes it.
        if shared_self_records {
            query_records = db_records
                .iter_mut()
                .enumerate()
                .map(|(id, record)| fasta::ByteFastaRecord {
                    id: std::mem::take(&mut record.id),
                    sequence: stage2_query_set.get(id).to_vec(),
                    quality: record.quality.take(),
                })
                .collect();
        }
        if let Some(state) = ranking_state {
            let mut output = capture
                .unwrap()
                .lock()
                .map_err(|_| io::Error::other("native search capture poisoned"))?;
            output.query_self_aln_scores = query_self_aln_scores;
            output.queries = query_records;
            output.query_oids = edge_query_oids.unwrap_or_default();
            if let Some(block) = translated_query_block.as_mut() {
                block.put_seqs(stage2_query_set);
            }
            output.translated = translated_query_block;
            output.matrix = Some(score_matrix);
            output.ranking = Some(state);
            return Ok(());
        }
        drop(db_records);
        let contexts = stage2_query_set.len();
        if config.translated_query_layout.is_none() {
            stage2_query_set.clear();
            stage2_query_set.shrink_to_fit();
        }
        // This empty retained-hit container only gives the shared writer its
        // complete query range; no seed hits or fake alignments are fabricated.
        AdaptiveHitStore::new_partitioned(
            contexts,
            db_block.seqs().data().len() as u64,
            None,
            config.tmpdir.clone(),
            vec![contexts.max(1) as u32],
        )
    } else {
        // Ranking x-drop scores use the hard-masked biological residues,
        // before seed-only masks change the enumeration arena.
        let ranking_queries = ranking_state.as_ref().map(|_| stage2_query_set.clone());
        // Explicit search soft masking is a source flag set: retain automatic
        // motif masking and use the same table for every seed-enumeration pass.
        // Alignment/fingerprint arenas keep original residues throughout.
        let automatic_soft_masking = sensitivity::soft_masking_algo_bytes(
            &sensitivity::get_traits(config.sensitivity),
            &config.motif_masking,
            false,
            freq_masking,
        )
        .map_err(io::Error::other)?;
        let configured_soft_masking =
            soft_masking_override != MaskingAlgo::None || config.soft_masking != MaskingAlgo::None;
        let mut configured_query_table = None;
        if configured_soft_masking {
            let bits = automatic_soft_masking as u32
                | config.soft_masking as u32
                | soft_masking_override as u32;
            let remaining_bits = if soft_tantan_masker.is_some() {
                bits & !(MaskingAlgo::Tantan as u32)
            } else {
                bits
            };
            let owned_masking = (config.tuning.runtime.masking().is_none() && remaining_bits != 0)
                .then(|| {
                    crate::masking::Masking::new(
                        &score_matrix,
                        crate::masking::MaskingConfig {
                            tantan_min_mask_prob: config.tuning.tantan_min_mask_prob as f32,
                            max_motif_len: config.tuning.max_motif_len.max(0) as usize,
                            ..Default::default()
                        },
                    )
                });
            let masking = config.tuning.runtime.masking().or(owned_masking.as_ref());
            // A clustering override supplies the actual shared Source kernel.
            // Record its spans in the authoritative seed-block table too; a
            // separate masked FastaRecord view is insufficient for stage0.
            let apply_soft_mask = |sequence: &mut [Letter],
                                   id: usize,
                                   table: &mut crate::masking::MaskingTable|
             -> io::Result<()> {
                if bits & MaskingAlgo::Tantan as u32 != 0 {
                    if let Some(tantan) = soft_tantan_masker {
                        let ranges = tantan.mask_ranges(sequence);
                        let (front, back) = ranges.as_slices();
                        for &(begin, end) in front.iter().chain(back) {
                            table.add(id, begin, end, sequence);
                        }
                    }
                }
                if remaining_bits != 0 {
                    masking
                        .expect("configured masking kernel")
                        .apply_masking(sequence, remaining_bits, id, Some(table))
                        .map_err(io::Error::other)?;
                }
                Ok(())
            };
            let mut target_table = crate::masking::MaskingTable::new();
            for (id, record) in db_records.iter_mut().enumerate() {
                apply_soft_mask(&mut record.sequence, id, &mut target_table)?;
            }
            db_block.set_soft_masking_table(target_table.clone());
            let mut table = crate::masking::MaskingTable::new();
            if shared_self_records {
                table = target_table;
            } else if translated_query_block.is_some() {
                for id in 0..stage2_query_set.len() {
                    apply_soft_mask(stage2_query_set.get_mut(id), id, &mut table)?;
                }
                if query_seed_indexed || config.target_indexed {
                    // Indexed enumeration owns apply/remove on the authoritative
                    // Block. Supply it restored sequences plus this exact table.
                    table.remove(&mut stage2_query_set, 0, false);
                }
            } else {
                for (id, record) in query_records.iter_mut().enumerate() {
                    apply_soft_mask(&mut record.sequence, id, &mut table)?;
                }
            }
            configured_query_table = Some(table);
        }

        // Motif masking — ports C++ `Block::soft_mask(MOTIF)` invoked from
        // `enum_seeds` (enum_seeds.h:202). At default sensitivity DIAMOND
        // hard-masks any 8-letter window matching a curated motif before seed
        // enumeration; the seed iterator then drops seeds overlapping a motif.
        // After enumeration C++ restores the original letters via
        // `Block::remove_soft_masking` so alignment runs against the real
        // sequence — we mirror that with `restore_motifs` further down.
        //
        // Run in parallel across sequences with rayon. C++ also parallelises
        // motif masking via `mask_seqs` (`masking/masking.cpp:172`).
        let soft_masking = if config.soft_masking != MaskingAlgo::None {
            config.soft_masking
        } else if soft_masking_override != MaskingAlgo::None {
            soft_masking_override
        } else {
            automatic_soft_masking
        };
        let use_motif_masking =
            automatic_soft_masking == MaskingAlgo::Motif && !configured_soft_masking;
        let db_motif_saves: Vec<Vec<crate::masking::motifs::MotifMaskEntry>> = if use_motif_masking
        {
            db_records
                .par_iter_mut()
                .map(|r| {
                    crate::masking::motifs::mask_motifs_with_limit(
                        &mut r.sequence,
                        config.tuning.max_motif_len.max(0) as usize,
                    )
                })
                .collect()
        } else {
            vec![Vec::new(); db_records.len()]
        };
        let query_motif_saves: Vec<Vec<crate::masking::motifs::MotifMaskEntry>> =
            if use_motif_masking {
                if shared_self_records {
                    Vec::new()
                } else if translated_query_block.is_some() && !config.target_indexed {
                    (0..stage2_query_set.len())
                        .map(|query_id| {
                            crate::masking::motifs::mask_motifs_with_limit(
                                stage2_query_set.get_mut(query_id),
                                config.tuning.max_motif_len.max(0) as usize,
                            )
                        })
                        .collect()
                } else {
                    query_records
                        .par_iter_mut()
                        .map(|r| {
                            crate::masking::motifs::mask_motifs_with_limit(
                                &mut r.sequence,
                                config.tuning.max_motif_len.max(0) as usize,
                            )
                        })
                        .collect()
                }
            } else {
                vec![Vec::new(); query_records.len()]
            };
        let reduction = setup.reduction;
        let shapes: Vec<Shape> = (0..setup.shapes.count() as usize)
            .map(|id| setup.shapes.get(id).clone())
            .collect();
        search_message!(
            "Sensitivity: {:?}, shapes: {} (weights: {})",
            config.sensitivity,
            shapes.len(),
            shapes
                .iter()
                .map(|s| s.weight.to_string())
                .collect::<Vec<_>>()
                .join(",")
        );

        let source_seedp_bits = sensitivity::seedp_bits(
            setup.shapes.get(0).weight as i32,
            rayon::current_num_threads().max(1) as i32,
            source_index_chunks as i32,
            &reduction,
        );
        // The compact optimized array stores a 32-bit key after exactly ten
        // partition bits, and a 32-bit absolute location. Source PackedLoc
        // retains forty location bits and expands the partition count when
        // necessary to keep a complete seed key in SeedOffset.
        let compact_seed_layout = source_seedp_bits == 10
            && stage2_query_set.data().len() <= u32::MAX as usize
            && db_block.seqs().data().len() <= u32::MAX as usize;

        // Upstream target-indexed mode changes the seed encoding itself: both
        // sides use Murmur-hashed spaced seeds, while the persistent hash tables
        // filter the one-pass query enumeration. Reuse the source-mirrored stage0
        // implementation for that path instead of accepting the flag while still
        // running the ordinary spaced-factor arrays.
        let mut seed_index_state = if !compact_seed_layout
            || config.target_indexed
            || query_seed_indexed
            || linear_stage1_query
            || linear_stage1_combo
            || freq_masking
            || setup.minimizer_window > 0
            || !config.tuning.setup.shape_mask.is_empty()
            || config.tuning.setup.shapes != 0
            || config.tuning.tile_size != 1024
            || config.tuning.kmer_ranking
            || config.tuning.min_length_ratio != 0.0
            || config.tuning.hash_join != crate::util::algo::HashJoinConfig::default()
        {
            // The manual motif pass above prepared the compact seed view.
            // Block enumeration must discover its cached table on original
            // residues, then owns apply/remove for every histogram/chunk.
            if translated_query_block.is_some() && use_motif_masking {
                for (id, saved) in query_motif_saves.iter().enumerate() {
                    crate::masking::motifs::restore_motifs(stage2_query_set.get_mut(id), saved, 0);
                }
            }
            let shape_config = setup.shapes.clone();
            let mut query = if translated_query_block.is_some() {
                None
            } else {
                let mut query = Block::new();
                query.set_max_motif_len(config.tuning.max_motif_len);
                // Source aliases options.target and options.query for the
                // self block. Its query seed owner must still contain those
                // sequences, even though our separate query_records is empty.
                let seed_query_records = if shared_self_records {
                    &db_records
                } else {
                    &query_records
                };
                for (oid, record) in seed_query_records.iter().enumerate() {
                    query
                        .push_back_bytes(
                            stage2_query_set.get(oid),
                            Some(&record.id),
                            None,
                            oid as u64,
                            SequenceType::AminoAcid,
                            0,
                            false,
                        )
                        .map_err(io::Error::other)?;
                }
                Some(query)
            };
            if let Some(block) = translated_query_block.as_mut() {
                block.set_max_motif_len(config.tuning.max_motif_len);
            }
            if let Some(table) = configured_query_table.as_ref() {
                let block = if let Some(block) = translated_query_block.as_mut() {
                    block
                } else {
                    query.as_mut().expect("indexed protein query block")
                };
                block.set_soft_masking_table(table.clone());
            }
            let contexts = if config.translated_query_layout.is_some() {
                6
            } else {
                1
            };
            let encoding = if !query_seed_indexed && !config.target_indexed {
                crate::data::flags::SeedEncoding::SpacedFactor
            } else if config.contiguous_seed_mode {
                crate::data::flags::SeedEncoding::Contiguous
            } else {
                crate::data::flags::SeedEncoding::Hashed
            };
            let seedp_bits = source_seedp_bits;
            let mut histogram_cfg = crate::data::flags::EnumCfg {
                partition: None,
                shape_begin: 0,
                shape_end: shape_config.count(),
                code: encoding,
                skip: None,
                filter_masked_seeds: false,
                mask_seeds: false,
                seed_cut: setup.seed_complexity_cut,
                soft_masking,
                minimizer_window: setup.minimizer_window,
                filter_low_complexity_seeds: false,
                mask_low_complexity_seeds: false,
                sketch_size: setup.sketch_size,
            };
            let mut query_seeds = None;
            let mut contiguous_seeds = None;
            let mut target_seeds = None;
            if config.contiguous_seed_mode {
                if let Some(block) = translated_query_block.as_mut() {
                    block.put_seqs(std::mem::take(&mut stage2_query_set));
                }
                let query_block = if let Some(block) = translated_query_block.as_mut() {
                    block
                } else {
                    query.as_mut().expect("protein contiguous query block")
                };
                histogram_cfg.seed_cut = setup.seed_complexity_cut;
                let seeds = crate::data::seed_set::SeedSet::from_block_with_min_query_len(
                    query_block,
                    1.0,
                    round_skip,
                    histogram_cfg.seed_cut,
                    soft_masking,
                    &shape_config,
                    &reduction,
                    contexts,
                    config.min_query_len,
                )
                .map_err(io::Error::other)?;
                let histogram = crate::data::seed_histogram::SeedHistogram::from_block_with_skip(
                    query_block,
                    false,
                    &crate::data::flags::NO_FILTER,
                    &histogram_cfg,
                    seedp_bits,
                    rayon::current_num_threads().max(1) as u32,
                    &shape_config,
                    &reduction,
                    config.min_query_len,
                    contexts,
                    round_skip,
                )
                .map_err(io::Error::other)?;
                *query_block.hst() = histogram;
                if let Some(block) = translated_query_block.as_mut() {
                    stage2_query_set = block.take_seqs();
                }
                histogram_cfg.soft_masking = MaskingAlgo::None;
                let histogram = crate::data::seed_histogram::SeedHistogram::from_block(
                    &mut db_block,
                    true,
                    &seeds,
                    &histogram_cfg,
                    seedp_bits,
                    rayon::current_num_threads().max(1) as u32,
                    &shape_config,
                    &reduction,
                    config.min_query_len,
                    contexts,
                )
                .map_err(io::Error::other)?;
                *db_block.hst() = histogram;
                contiguous_seeds = Some(seeds);
            } else if query_seed_indexed {
                histogram_cfg.minimizer_window = setup.minimizer_window;
                histogram_cfg.sketch_size = setup.sketch_size;
                if let Some(block) = translated_query_block.as_mut() {
                    block.put_seqs(std::mem::take(&mut stage2_query_set));
                }
                let query_block = if let Some(block) = translated_query_block.as_mut() {
                    block
                } else {
                    query.as_mut().expect("protein query-indexed block")
                };
                let seeds = config
                    .tuning
                    .algorithm
                    .cached_seeds()?
                    .ok_or_else(|| io::Error::other("Query seed cache was not initialized."))?;
                histogram_cfg.seed_cut = setup.seed_complexity_cut;
                histogram_cfg.filter_low_complexity_seeds = true;
                let histogram = crate::data::seed_histogram::SeedHistogram::from_block_with_skip(
                    query_block,
                    false,
                    &crate::data::flags::NO_FILTER,
                    &histogram_cfg,
                    seedp_bits,
                    rayon::current_num_threads().max(1) as u32,
                    &shape_config,
                    &reduction,
                    config.min_query_len as i32,
                    contexts,
                    round_skip,
                )
                .map_err(io::Error::other)?;
                *query_block.hst() = histogram;
                if let Some(block) = translated_query_block.as_mut() {
                    stage2_query_set = block.take_seqs();
                }
                histogram_cfg.filter_low_complexity_seeds = false;
                histogram_cfg.soft_masking = MaskingAlgo::None;
                let histogram = crate::data::seed_histogram::SeedHistogram::from_block(
                    &mut db_block,
                    false,
                    seeds.as_ref(),
                    &histogram_cfg,
                    seedp_bits,
                    rayon::current_num_threads().max(1) as u32,
                    &shape_config,
                    &reduction,
                    config.min_query_len,
                    contexts,
                )
                .map_err(io::Error::other)?;
                *db_block.hst() = histogram;
                query_seeds = Some(seeds);
            } else {
                let histogram = crate::data::seed_histogram::SeedHistogram::from_block(
                    &mut db_block,
                    false,
                    &crate::data::flags::NO_FILTER,
                    &histogram_cfg,
                    seedp_bits,
                    rayon::current_num_threads().max(1) as u32,
                    &shape_config,
                    &reduction,
                    config.min_query_len,
                    contexts,
                )
                .map_err(io::Error::other)?;
                *db_block.hst() = histogram;
                if !config.target_indexed {
                    let query_block = if let Some(block) = translated_query_block.as_mut() {
                        block.put_seqs(std::mem::take(&mut stage2_query_set));
                        block
                    } else {
                        query
                            .as_mut()
                            .expect("protein frequency-filter query block")
                    };
                    let histogram =
                        crate::data::seed_histogram::SeedHistogram::from_block_with_skip(
                            query_block,
                            false,
                            &crate::data::flags::NO_FILTER,
                            &histogram_cfg,
                            seedp_bits,
                            rayon::current_num_threads().max(1) as u32,
                            &shape_config,
                            &reduction,
                            config.min_query_len,
                            contexts,
                            round_skip,
                        )
                        .map_err(io::Error::other)?;
                    *query_block.hst() = histogram;
                    if let Some(block) = translated_query_block.as_mut() {
                        stage2_query_set = block.take_seqs();
                    }
                } else {
                    let mut database_path = PathBuf::from(&config.database);
                    if !crate::util::io::source_c_path(&database_path).exists() {
                        database_path.as_mut_os_string().push(".dmnd");
                    }
                    database_path.as_mut_os_string().push(".seed_idx");
                    target_seeds = Some(
                        crate::data::seed_set::HashedSeedSet::from_index_file(
                            &database_path,
                            &shape_config,
                        )
                        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?,
                    );
                    let mut message = b"Loaded persistent seed index: ".to_vec();
                    message.extend_from_slice(
                        &crate::util::io::source_os_bytes(database_path.as_os_str())
                            .map_err(io::Error::other)?,
                    );
                    message.push(b'\n');
                    crate::util::log_stream::message_stream()
                        .lock()
                        .unwrap()
                        .write_bytes(&message)?;
                }
            }
            Some(SeedIndexState {
                query,
                shapes: shape_config,
                target_seeds,
                query_seeds,
                contiguous_seeds,
                encoding,
                seedp_bits,
            })
        } else {
            None
        };

        // Build partitioned seed arrays and join for each shape (parallel)
        let db_seqs: Vec<&[Letter]> = db_records.iter().map(|r| r.sequence.as_slice()).collect();
        let query_seqs: Vec<&[Letter]> = if shared_self_records {
            db_records
                .iter()
                .map(|record| record.sequence.as_slice())
                .collect()
        } else if translated_query_block.is_some() {
            // The authoritative translated arena is mutated between seed
            // enumeration and Hamming filtering, so it is borrowed only inside
            // each phase below.
            Vec::new()
        } else {
            query_records
                .iter()
                .map(|record| record.sequence.as_slice())
                .collect()
        };
        // Per-sensitivity seed filters matching C++ DIAMOND:
        //   - `seed_cut * ln(2) * shape.weight` is the entropy floor for
        //     `seed_is_complex` (search/setup.cpp).
        //   - Frequent-seed masking is guarded by C++ `config.freq_masking`; the
        //     native default path leaves it disabled.
        let hamming_filter_id = setup.hamming_filter_id;
        // Process shapes in order and let each shape's seed-array build/join use
        // the Rayon pool internally. Running the outer shape loop in parallel
        // competes with the inner reference seed-array work and scales worse on
        // small real query blocks; collecting in shape order preserves the previous
        // output order.
        let ungapped_evalue = setup.ungapped_evalue;
        let short_query_ungapped_cutoff = if ungapped_evalue > 0.0 {
            score_matrix.rawscore_int(config.tuning.short_query_ungapped_bitscore)
        } else {
            0
        };
        let ungapped_cutoffs: Vec<Stage2Filter> = (0..stage2_query_set.len())
            .map(|query_id| {
                let query = stage2_query_set.get(query_id);
                Stage2Filter {
                    score: crate::search::stage2::ungapped_cutoff(
                        query.len() as i32,
                        ungapped_evalue,
                        config.tuning.short_query_max_len,
                        short_query_ungapped_cutoff,
                        config.translated_query_layout.is_some(),
                        |q| setup.cutoff_table_short.call(q),
                        |q| setup.cutoff_table.call(q),
                    ),
                    window: crate::search::stage2::ungapped_window(
                        query.len() as i32,
                        config.translated_query_layout.is_some(),
                        config.tuning.ungapped_window,
                    )
                    .max(0) as usize,
                    left_most_interval: config.tuning.left_most_interval,
                    self_search: config.tuning.self_search
                        && config.tuning.search_block.reference == 0,
                }
            })
            .collect();
        // Left-most filtering consumes SEED_MASK bits. Restore the pre-motif query
        // copy now; each shape will add its low-complexity bits after joining and
        // then search those same prepared arrays, matching upstream stage0.
        let mut low_complexity_query_positions = Vec::new();
        let template_len = shapes.iter().map(|shape| shape.length).max().unwrap_or(0);
        let query_motif_saves_for_restore = if shared_self_records {
            &db_motif_saves
        } else {
            &query_motif_saves
        };
        if translated_query_block.is_none() {
            for (query_id, saved) in query_motif_saves_for_restore.iter().enumerate() {
                crate::masking::motifs::restore_motifs(
                    stage2_query_set.get_mut(query_id),
                    saved,
                    template_len,
                );
            }
            if let Some(table) = configured_query_table.as_ref() {
                table.remove(&mut stage2_query_set, template_len, true);
            }
        }
        // Hamming fingerprints in upstream are loaded after motif masking has
        // been removed, directly from the original contiguous sequences. Use the
        // pre-motif stage-2 copies here as well: seed enumeration still consumes
        // the masked `query_seqs`/`db_seqs`, while fingerprint loading no longer
        // needs a sparse restoration lookup for every candidate window. Query
        // SEED_MASK annotations are harmless because the fingerprint loader
        // strips LETTER_MASK from every byte.
        let effective_index_chunks = source_index_chunks;
        let chunked = effective_index_chunks > 1;
        let use_left_most_range = chunked;
        // Upstream stage2 skips this filter only for minimizer/sketch and linear
        // search modes; sensitivity itself is not a bypass. Keeping all hits for
        // very/ultra-sensitive searches inflated the spill file by two orders of
        // magnitude (and subsequently decoded those unnecessary hits into RAM).
        let skip_left_most =
            setup.minimizer_window > 0 || setup.sketch_size > 0 || linear_stage1_query;
        // Select the stage-2 kernel once. A tiny const-generic trampoline keeps
        // CPU-feature branches out of the batch loop without cloning the much
        // larger partition worker for every ISA.
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        let use_avx2_ungapped = std::arch::is_x86_feature_detected!("avx2");
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        let use_sse41_ungapped = std::arch::is_x86_feature_detected!("ssse3")
            && std::arch::is_x86_feature_detected!("sse4.1");
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        let use_avx2_ungapped = false;
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        let use_sse41_ungapped = false;
        let (ungapped_lane_count, ungapped_kernel): (usize, UngappedKernel) = if use_avx2_ungapped {
            (32, run_ungapped_kernel::<UNGAPPED_AVX2>)
        } else if use_sse41_ungapped {
            (16, run_ungapped_kernel::<UNGAPPED_SSE41>)
        } else {
            (16, run_ungapped_kernel::<UNGAPPED_PORTABLE>)
        };
        // Upstream instantiates the complete stage-1 worker in its AVX2 namespace.
        // Select that worker once per run so the partition hot path contains no
        // CPU-feature branch.
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        let subject_index_partition_filter: SubjectIndexPartitionFilter =
            if std::arch::is_x86_feature_detected!("avx2") {
                filter_subject_indices_partition_to_hits_avx2
            } else {
                filter_subject_indices_partition_to_hits_baseline
            };
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        let subject_index_partition_filter: SubjectIndexPartitionFilter =
            filter_subject_indices_partition_to_hits_baseline;
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        let borrowed_partition_filter: BorrowedPartitionFilter =
            if std::arch::is_x86_feature_detected!("avx2") {
                filter_borrowed_partition_to_hits_avx2
            } else {
                filter_borrowed_partition_to_hits_baseline
            };
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        let borrowed_partition_filter: BorrowedPartitionFilter =
            filter_borrowed_partition_to_hits_baseline;
        let contexts = if config.translated_query_layout.is_some() {
            6
        } else {
            1
        };
        let mut key_partition = if stage2_query_set.letters() > 0 {
            stage2_query_set
                .partition_with_contexts(setup.query_bins, true, true, contexts)
                .into_iter()
                .map(|oid| oid * contexts as u32)
                .collect::<Vec<_>>()
        } else {
            AdaptiveHitStore::key_partition(
                stage2_query_set.len(),
                setup.query_bins as usize,
                contexts,
            )
        };
        // Repeated empty final bins have no observable records and no need for a
        // physical buffer. Source OID boundaries of every populated bin remain.
        key_partition.dedup();
        let mut hit_store = AdaptiveHitStore::new_partitioned(
            stage2_query_set.len(),
            db_block.seqs().data().len() as u64,
            if ranking_state.is_some() {
                None
            } else {
                config.memory_limit
            },
            config.tmpdir.clone(),
            key_partition,
        );
        hit_store.no_unlink = config.tuning.no_unlink;
        hit_store.load_threads = config.tuning.load_threads;
        if ranking_state.is_none() {
            hit_store.configure_storage(config.tuning.hit_membuf)?;
        }
        let mut indexed_hit_targets = vec![
            false;
            if query_seed_indexed {
                db_block.seqs().len()
            } else {
                0
            }
        ];
        let mut retained_hit_count = 0usize;
        let mut raw_seed_matches = 0usize;
        let mut current_left_most_matcher = PatternMatcher::new(&[]);
        let translated_histogram_query_seqs: Vec<&[Letter]> = if translated_query_block.is_some() {
            (0..stage2_query_set.len())
                .map(|query_id| stage2_query_set.get(query_id))
                .collect()
        } else {
            Vec::new()
        };
        let histogram_query_seqs = if translated_query_block.is_some() {
            &translated_histogram_query_seqs
        } else {
            &query_seqs
        };
        // Upstream constructs one worker×shape×partition SeedHistogram for the
        // complete query/reference block, then reuses it for every index chunk.
        // Keep worker sequence ranges in the histogram so parallel fills preserve
        // the same stable record order without atomics.
        let precomputed_worker_histograms = (seed_index_state.is_none()
        && rayon::current_num_threads() > 1
        && setup.sketch_size <= 0)
        .then(|| {
            (
                crate::search::seed_array::SeedArray::count_partitions_by_worker_all_shapes_with_skip(
                    histogram_query_seqs,
                    &shapes,
                    &reduction,
                    10,
                    config.min_query_len,
                    round_skip,
                    query_contexts,
                ),
                crate::search::seed_array::SeedArray::count_partitions_by_worker_all_shapes_with_skip(
                    &db_seqs,
                    &shapes,
                    &reduction,
                    10,
                    config.min_query_len,
                    None,
                    query_contexts,
                ),
            )
        });
        let _ = histogram_query_seqs;
        drop(translated_histogram_query_seqs);
        for (shape_id, shape) in shapes.iter().enumerate() {
            let previous_left_most_matcher = current_left_most_matcher.clone();
            current_left_most_matcher.add_pattern(shape.mask);
            let left_most_context = LeftMostContext {
                previous_matcher: previous_left_most_matcher,
                current_matcher: current_left_most_matcher.clone(),
                short_query_ungapped_cutoff: score_matrix
                    .rawscore_int(config.tuning.short_query_ungapped_bitscore),
                seedp_mask: seedp_mask(source_seedp_bits),
                reduction: &reduction,
            };
            let complexity_cut = setup.seed_complexity_cut;

            if let Some(state) = seed_index_state.as_mut() {
                let mut stage0_config =
                    crate::search::stage0::SearchShapeConfig::new(&state.shapes, &reduction);
                stage0_config.hash_join = config.tuning.hash_join;
                stage0_config.freq_masking = freq_masking;
                stage0_config.freq_sd = setup.freq_sd;
                stage0_config.seedp_bits = state.seedp_bits;
                stage0_config.index_chunks = effective_index_chunks;
                stage0_config.threads = rayon::current_num_threads().max(1);
                stage0_config.seed_encoding = state.encoding;
                stage0_config.seed_complexity_cut = complexity_cut;
                stage0_config.soft_masking = soft_masking;
                stage0_config.minimizer_window = setup.minimizer_window;
                stage0_config.sketch_size = setup.sketch_size;
                stage0_config.min_query_len = config.min_query_len as i32;
                stage0_config.query_contexts = if config.translated_query_layout.is_some() {
                    6
                } else {
                    1
                };
                stage0_config.linear_stage1_query = linear_stage1_query;
                stage0_config.linear_stage1_target = linear_stage1_target;
                stage0_config.keep_target_id = retain_seed_ids;
                stage0_config.kmer_ranking = config.tuning.kmer_ranking;
                stage0_config.target_seeds = state.target_seeds.as_ref();
                if let Some(seeds) = state.contiguous_seeds.as_ref() {
                    stage0_config.reference_filter =
                        crate::search::stage0::ReferenceSeedFilter::BitSet(seeds);
                }
                if let Some(seeds) = state.query_seeds.as_ref() {
                    stage0_config.reference_filter =
                        crate::search::stage0::ReferenceSeedFilter::Hashed(seeds);
                }

                // Translated target-indexed search temporarily reunites the
                // context arena with its source/ID block. `search_shape` performs
                // soft mask -> enumeration -> restore on that same block.
                if let Some(block) = translated_query_block.as_mut() {
                    block.put_seqs(std::mem::take(&mut stage2_query_set));
                }
                let query_stage0_block = if let Some(block) = translated_query_block.as_mut() {
                    block
                } else {
                    state.query.as_mut().expect("protein indexed query block")
                };
                let indexed_report = crate::search::stage0::search_shape_with_skip(
                    shape_id as u32,
                    query_stage0_block,
                    &mut db_block,
                    &stage0_config,
                    round_skip,
                    |partition, query_stage0, target_stage0| {
                        let (hits, raw) = map_target_indexed_partition(
                            partition.query_seed_hits,
                            partition.reference_seed_hits,
                            partition.absolute_partition,
                            query_stage0,
                            target_stage0,
                            hamming_filter_id,
                            config.tuning.tile_size,
                            &crate::search::hamming::Stage1DispatchConfig {
                                lin_stage1_query: linear_stage1_query,
                                lin_stage1_target: linear_stage1_target,
                                lin_stage1_combo: linear_stage1_combo,
                                min_length_ratio,
                                self_search: config.tuning.self_search,
                                current_ref_block: config.tuning.search_block.reference,
                                global_ranking_targets: config.global_ranking_targets != 0,
                                ..Default::default()
                            },
                            retain_seed_ids,
                            partition.kmer_ranking,
                            |query,
                             target_indices,
                             target_begin,
                             target_locations,
                             absolute_partition,
                             out| {
                                // SAFETY: the AVX2 function item is installed only
                                // after the one-time feature check above.
                                unsafe {
                                    borrowed_partition_filter(
                                        query,
                                        target_indices,
                                        target_begin,
                                        target_locations,
                                        absolute_partition,
                                        out,
                                        query_stage0,
                                        target_stage0,
                                        &ungapped_cutoffs,
                                        &score_matrix,
                                        shape,
                                        &left_most_context,
                                        shape_id == 0,
                                        chunked,
                                        use_left_most_range,
                                        skip_left_most,
                                        effective_index_chunks,
                                        hamming_filter_id,
                                        ungapped_lane_count,
                                        ungapped_kernel,
                                    );
                                }
                            },
                        )?;
                        if query_seed_indexed {
                            for hit in &hits {
                                let (target_id, _) =
                                    target_stage0.local_position(hit.subject as usize);
                                indexed_hit_targets[target_id as usize] = true;
                            }
                        }
                        raw_seed_matches = raw_seed_matches.saturating_add(raw);
                        retained_hit_count += hit_store
                            .ingest_partitions(vec![hits])
                            .map_err(|error| error.to_string())?;
                        Ok(())
                    },
                )
                .map_err(io::Error::other)?;
                if let Some(block) = translated_query_block.as_mut() {
                    stage2_query_set = block.take_seqs();
                }
                search_message!(
                    "Indexed seeds (shape {}): query={}, reference={}",
                    shape_id + 1,
                    indexed_report
                        .chunks
                        .iter()
                        .map(|chunk| chunk.query_seed_count)
                        .sum::<usize>(),
                    indexed_report
                        .chunks
                        .iter()
                        .map(|chunk| chunk.reference_seed_count)
                        .sum::<usize>()
                );
                if let Some(state) = ranking_state.as_mut() {
                    hit_store.update_ranking_shape(
                        state,
                        ranking_queries.as_ref().unwrap(),
                        &db_block,
                        config,
                        &score_matrix,
                    )?;
                }
                continue;
            }

            // Keep this borrow inside the ordinary seed-array branch. Persistent
            // index search mutably reuses `db_block` while enumerating hashed
            // seeds, whereas this fingerprint view is only needed by the native
            // spaced-factor join below.
            let db_fingerprint_seqs: Vec<&[Letter]> = (0..db_block.seqs().len())
                .map(|id| db_block.seqs().get(id))
                .collect();

            let partition_count = seedp_count(10) as usize;
            let index_chunks =
                crate::util::algo::Partition::new(partition_count, effective_index_chunks);
            // Upstream computes one histogram for the complete shape, then uses
            // it to allocate/fill every index chunk.  Preserve that ordering in
            // the one-thread non-sketch path; the generic builder otherwise
            // performs a private count pass for every chunk.  Multithread and
            // sketch construction retain their existing worker-local layouts.
            let precomputed_counts = (rayon::current_num_threads() == 1 && setup.sketch_size <= 0)
                .then(|| {
                    let translated_query_seqs = if translated_query_block.is_some() {
                        (0..stage2_query_set.len())
                            .map(|query_id| stage2_query_set.get(query_id))
                            .collect::<Vec<_>>()
                    } else {
                        Vec::new()
                    };
                    let shape_query_seqs = if translated_query_block.is_some() {
                        &translated_query_seqs
                    } else {
                        &query_seqs
                    };
                    (
                        crate::search::seed_array::SeedArray::count_partitions_with_skip(
                            shape_query_seqs,
                            shape,
                            &reduction,
                            10,
                            0.0,
                            config.min_query_len,
                            round_skip,
                            query_contexts,
                        ),
                        crate::search::seed_array::SeedArray::count_partitions_with_skip(
                            &db_seqs,
                            shape,
                            &reduction,
                            10,
                            0.0,
                            config.min_query_len,
                            None,
                            query_contexts,
                        ),
                    )
                });
            let mut shape_raw_seed_matches = 0usize;
            for chunk in 0..index_chunks.parts {
                let translated_query_seqs = if translated_query_block.is_some() {
                    (0..stage2_query_set.len())
                        .map(|query_id| stage2_query_set.get(query_id))
                        .collect::<Vec<_>>()
                } else {
                    Vec::new()
                };
                let shape_query_seqs = if translated_query_block.is_some() {
                    &translated_query_seqs
                } else {
                    &query_seqs
                };
                let prepared = if let Some((query_counts, reference_counts)) =
                    &precomputed_worker_histograms
                {
                    parallel::prepare_seed_join_partition_range_from_worker_counts_min_query_len_with_skip(
                        shape_query_seqs,
                        &db_seqs,
                        shape,
                        &reduction,
                        complexity_cut,
                        config.min_query_len,
                        index_chunks.begin(chunk),
                        index_chunks.end(chunk),
                        &query_counts[shape_id],
                        &reference_counts[shape_id],
                        round_skip,
                        query_contexts,
                    )
                } else if let Some((query_counts, reference_counts)) = &precomputed_counts {
                    parallel::prepare_seed_join_partition_range_from_counts_min_query_len_with_skip(
                        shape_query_seqs,
                        &db_seqs,
                        shape,
                        &reduction,
                        complexity_cut,
                        config.min_query_len,
                        index_chunks.begin(chunk),
                        index_chunks.end(chunk),
                        query_counts,
                        reference_counts,
                        round_skip,
                        query_contexts,
                    )
                } else {
                    parallel::prepare_seed_join_partition_range_min_query_len_with_skip(
                        shape_query_seqs,
                        &db_seqs,
                        shape,
                        &reduction,
                        complexity_cut,
                        config.min_query_len,
                        setup.sketch_size.max(0) as usize,
                        index_chunks.begin(chunk),
                        index_chunks.end(chunk),
                        round_skip,
                        query_contexts,
                    )
                };
                let _ = shape_query_seqs;
                drop(translated_query_seqs);
                if translated_query_block.is_some() {
                    if let Some(table) = configured_query_table.as_ref() {
                        table.remove(&mut stage2_query_set, template_len, true);
                    }
                    for (query_id, saved) in query_motif_saves.iter().enumerate() {
                        crate::masking::motifs::restore_motifs(
                            stage2_query_set.get_mut(query_id),
                            saved,
                            template_len,
                        );
                    }
                }
                low_complexity_query_positions.extend_from_slice(prepared.masked_positions());
                for &(query_id, query_pos) in &low_complexity_query_positions {
                    if (query_id as usize) < stage2_query_set.len() {
                        if let Some(letter) = stage2_query_set
                            .get_mut(query_id as usize)
                            .get_mut(query_pos as usize)
                        {
                            *letter |= SEED_MASK;
                        }
                    }
                }
                let query_fingerprint_seqs: Vec<&[Letter]> = (0..stage2_query_set.len())
                    .map(|id| stage2_query_set.get(id))
                    .collect();
                let map_candidates =
                    |query: (u32, u32),
                     target_indices: &[u32],
                     _target_locations: &[(u32, u32)],
                     target_subjects: &[u64],
                     absolute_partition: u32,
                     out: &mut Vec<StoredHit>| {
                        // SAFETY: the AVX2 function item is installed only after the
                        // one-time feature check above; the baseline item has no extra
                        // instruction-set requirements.
                        unsafe {
                            subject_index_partition_filter(
                                query,
                                target_indices,
                                target_subjects,
                                absolute_partition,
                                out,
                                &stage2_query_set,
                                db_block.seqs(),
                                &ungapped_cutoffs,
                                &score_matrix,
                                shape,
                                &left_most_context,
                                shape_id == 0,
                                chunked,
                                use_left_most_range,
                                skip_left_most,
                                effective_index_chunks,
                                hamming_filter_id,
                                ungapped_lane_count,
                                ungapped_kernel,
                            )
                        }
                    };
                let chunk_raw_seed_matches = if config.memory_limit.is_none() {
                    let mut partitions = Vec::new();
                    let raw =
                        parallel::visit_prepared_seed_hits_streaming_hamming_mode_with_options(
                            &query_fingerprint_seqs,
                            &db_fingerprint_seqs,
                            Some((&stage2_query_set, db_block.seqs())),
                            prepared,
                            hamming_filter_id,
                            usize::MAX,
                            linear_stage1_query,
                            linear_stage1_target,
                            min_length_ratio,
                            config.tuning.self_search && config.tuning.search_block.reference == 0,
                            parallel::PreparedStage1Options {
                                kmer_ranking: config.tuning.kmer_ranking,
                                keep_target_id: retain_seed_ids,
                                tile_size: config.tuning.tile_size,
                            },
                            map_candidates,
                            |batch| partitions.extend(batch),
                        );
                    retained_hit_count += hit_store.ingest_partitions(partitions)?;
                    raw
                } else {
                    let mut spill_error = None;
                    let raw =
                        parallel::visit_prepared_seed_hits_streaming_hamming_mode_with_options(
                            &query_fingerprint_seqs,
                            &db_fingerprint_seqs,
                            Some((&stage2_query_set, db_block.seqs())),
                            prepared,
                            hamming_filter_id,
                            // Keep enough independent partitions in flight that skewed
                            // seed groups do not leave workers idle at every disk-writer
                            // handoff. This remains bounded (and far below a complete
                            // 1024-partition shape) for forced-disk RSS control.
                            rayon::current_num_threads().max(1) * 15,
                            linear_stage1_query,
                            linear_stage1_target,
                            min_length_ratio,
                            config.tuning.self_search && config.tuning.search_block.reference == 0,
                            parallel::PreparedStage1Options {
                                kmer_ranking: config.tuning.kmer_ranking,
                                keep_target_id: retain_seed_ids,
                                tile_size: config.tuning.tile_size,
                            },
                            map_candidates,
                            |partitions| {
                                if spill_error.is_none() {
                                    match hit_store.ingest_partitions(partitions) {
                                        Ok(count) => retained_hit_count += count,
                                        Err(error) => spill_error = Some(error),
                                    }
                                }
                            },
                        );
                    if let Some(error) = spill_error {
                        return Err(error);
                    }
                    raw
                };
                shape_raw_seed_matches += chunk_raw_seed_matches;
                if translated_query_block.is_some()
                    && (chunk + 1 < index_chunks.parts || shape_id + 1 < shapes.len())
                {
                    if let Some(table) = configured_query_table.as_ref() {
                        table.apply(&mut stage2_query_set);
                    }
                    for (query_id, saved) in query_motif_saves.iter().enumerate() {
                        crate::masking::motifs::reapply_motifs(
                            stage2_query_set.get_mut(query_id),
                            saved,
                        );
                    }
                }
            }
            raw_seed_matches += shape_raw_seed_matches;
            if let Some(state) = ranking_state.as_mut() {
                hit_store.update_ranking_shape(
                    state,
                    ranking_queries.as_ref().unwrap(),
                    &db_block,
                    config,
                    &score_matrix,
                )?;
            }
        }
        // C++ stage0/stage2 preserves shape + seed-partition join emission order
        // inside a query. Only group by query so range building is cheap without
        // imposing a Rust-only target/position tie-breaker on ranking boundaries.
        search_message!(
            "Seed matches: {} (raw) -> {} (left-most)",
            raw_seed_matches,
            retained_hit_count,
        );

        let _ = crate::util::system::log_rss();

        // Restore motif-masked positions before alignment so the gapped extension
        // scores against the real letters (C++ `Block::remove_soft_masking`).
        // The earlier `db_seqs`/`query_seqs` slices are dropped here so we can
        // re-borrow the underlying records mutably.
        {
            let _ = &db_seqs;
            let _ = &query_seqs;
        }
        drop(db_seqs);
        drop(query_seqs);
        trim_freed_heap_pages();
        if translated_query_block.is_none() {
            if let Some(block) = seed_index_state
                .as_mut()
                .and_then(|state| state.query.as_mut())
            {
                stage2_query_set = block.take_seqs();
            }
        }
        // C++ `enum_seeds.h:223` calls `seqs.remove_soft_masking(template_len, mask_seeds)`
        // with `mask_seeds=true` only on the query side, propagating SEED_MASK over
        // `[motif_begin - template_len + 1, motif_begin + motif_len)`. The DB-side
        // restore (stage0.cpp:142-144 with `mask_seeds=false`) just rewrites
        // letters. `template_len` is `max(shape.length_)` — derive from active shapes.
        for (record, saved) in db_records.iter_mut().zip(db_motif_saves.iter()) {
            crate::masking::motifs::restore_motifs(&mut record.sequence, saved, 0);
        }
        if shared_self_records {
            // The aliased upstream block becomes the extension-side query only
            // after seed enumeration. Move its IDs and copy the already compact,
            // restored stage-2 sequences; no second FastaRecord set exists during
            // the peak seed-search phase.
            query_records = db_records
                .iter_mut()
                .enumerate()
                .map(|(query_id, record)| fasta::ByteFastaRecord {
                    quality: record.quality.take(),
                    id: std::mem::take(&mut record.id),
                    sequence: stage2_query_set.get(query_id).to_vec(),
                })
                .collect();
            for (query_id, query_pos) in low_complexity_query_positions {
                if let Some(letter) = query_records
                    .get_mut(query_id as usize)
                    .and_then(|record| record.sequence.get_mut(query_pos as usize))
                {
                    *letter |= SEED_MASK;
                }
            }
        } else {
            for (record, saved) in query_records.iter_mut().zip(query_motif_saves.iter()) {
                crate::masking::motifs::restore_motifs(&mut record.sequence, saved, template_len);
            }
            // C++ `Search::mask_seeds` sets the high bit on every query occurrence
            // in a joined seed group rejected by the low-complexity test.
            for (query_id, query_pos) in low_complexity_query_positions {
                if let Some(letter) = query_records
                    .get_mut(query_id as usize)
                    .and_then(|record| record.sequence.get_mut(query_pos as usize))
                {
                    *letter |= SEED_MASK;
                }
            }

            // `Block::remove_soft_masking` restores the same query storage that is
            // later consumed by extension. Move that exact state back before
            // releasing the temporary SequenceSet.
            if config.translated_query_layout.is_none() {
                restore_extension_queries(&mut query_records, &stage2_query_set);
            }
        }

        // Source run_ref_chunk clears seed-only annotations after every shape
        // has finished, before extension and before retaining the query owner
        // for the next reference block. Hard-masked residues remain masked.
        crate::data::frequent_seeds::FrequentSeeds::clear_masking(&mut stage2_query_set);
        for record in &mut query_records {
            crate::masking::remove_bit_mask(&mut record.sequence);
        }
        if shared_self_records {
            // Source aliases the query and reference arenas in this branch.
            crate::data::frequent_seeds::FrequentSeeds::clear_masking(db_block.seqs_mut());
        }
        let _ = crate::util::system::log_rss();

        if let Some(state) = ranking_state {
            let mut output = capture
                .unwrap()
                .lock()
                .map_err(|_| io::Error::other("native search capture poisoned"))?;
            output.query_self_aln_scores = query_self_aln_scores;
            output.queries = query_records;
            output.query_oids = edge_query_oids.unwrap_or_default();
            if let Some(block) = translated_query_block.as_mut() {
                block.put_seqs(stage2_query_set);
            }
            output.translated = translated_query_block;
            output.matrix = Some(score_matrix);
            output.ranking = Some(state);
            return Ok(());
        }

        // Query-indexed seeding reads original targets. Once all retained-hit
        // candidates are known, apply the same pure hard-mask transformation to
        // those targets before any query extension. This ownership boundary keeps
        // extension parallel without concurrently mutating a shared sequence arena.
        if query_seed_indexed && config.frame_shift == 0 {
            for (target_id, selected) in indexed_hit_targets.into_iter().enumerate() {
                if !selected {
                    continue;
                }
                let sequence = db_block.seqs_mut().get_mut(target_id);
                match config.masking {
                    MaskingMode::Tantan => {
                        tantan_masker
                            .expect("configured tantan masker")
                            .mask(sequence);
                        hard_mask(sequence);
                    }
                    MaskingMode::BlastSeg => {
                        crate::masking::mask_sequence(sequence, MaskingAlgo::Seg)
                    }
                    MaskingMode::None => (),
                }
            }
        }

        // Release seed-enumeration record storage and the temporary query copy
        // before extension.
        drop(db_records);
        // Translated extension borrows the same contiguous query arena. Ordinary
        // blastp releases its arena after restoring record-owned extension state.
        if config.translated_query_layout.is_none() {
            stage2_query_set.clear();
            stage2_query_set.shrink_to_fit();
        }
        drop(db_motif_saves);
        drop(query_motif_saves);
        trim_freed_heap_pages();
        hit_store
    };
    let db_ids = db_block.ids().map_err(io::Error::other)?;

    // Once seeding is complete, build lightweight per-source views over the
    // authoritative block and contiguous context arena. Extension borrows all
    // six frames and the source DNA; no sequence bytes are moved or cloned.
    let translated_groups = if let Some(layout) = &config.translated_query_layout {
        let expected = layout.mate_indices.len().saturating_mul(6);
        if stage2_query_set.len() != expected {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "blastx internal layout has {} protein records, expected {expected}",
                    stage2_query_set.len()
                ),
            ));
        }
        let block = translated_query_block
            .as_ref()
            .expect("translated layout requires its source query block");
        let ids = block.ids().map_err(io::Error::other)?;
        let mut groups = Vec::with_capacity(layout.mate_indices.len());
        for (source_index, &mate_index) in layout.mate_indices.iter().enumerate() {
            let context = source_index * 6;
            let frames = std::array::from_fn(|frame| stage2_query_set.get(context + frame));
            groups.push(TranslatedExtensionGroup {
                id: ids.get(source_index),
                dna_len: block.source_seqs().length(source_index) as i32,
                source_sequence: block.source_seqs().get(source_index),
                quality: if block.qual().empty() {
                    b""
                } else {
                    block.qual().get(source_index)
                },
                mate_sequence: mate_index.map(|index| block.source_seqs().get(index)),
                frames,
            });
        }
        Some(groups)
    } else {
        None
    };

    // Parse the two native output paths. DAA uses the same completed `Match`
    // objects as tabular output, but serializes all HSP transcripts and builds
    // the reference dictionary in stable first-use order at the writer.
    let daa_output = config
        .outfmt
        .first()
        .is_some_and(|format| format.as_slice() == b"100" || format.eq_ignore_ascii_case(b"daa"));
    let fields = if daa_output {
        if config.outfmt.len() > 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "DAA output does not accept tabular field names",
            ));
        }
        Vec::new()
    } else if context_format.is_some() {
        Vec::new()
    } else if config.outfmt.is_empty()
        || matches!(
            config.outfmt[0].as_slice(),
            b"6" | b"tab" | b"104" | b"json-flat"
        )
    {
        if config.outfmt.len() > 1 {
            config.outfmt[1..]
                .iter()
                .filter_map(|f| FieldId::from_name_bytes(f))
                .collect()
        } else {
            format::DEFAULT_TABULAR_FIELDS.to_vec()
        }
    } else {
        // Reuse Source's selector/field validation and raw diagnostic, also
        // for direct in-memory library callers without an outer block writer.
        super::search_output::tabular_format(config)?;
        return Err(io::Error::other("Invalid native output format routing"));
    };

    // Match upstream's field-driven DP requirements. In particular positive
    // counts request TRANSCRIPT; coordinates/statistics alone may select a
    // different equal-score alignment or bypass its final traceback.
    let tabular_hsp_values = fields.iter().fold(HspValues::NONE, |values, &field| {
        values | field.hsp_values()
    });

    let tabular_format = if !daa_output && context_format.is_none() && edge_output.is_none() {
        Some(super::search_output::tabular_format(config)?)
    } else {
        None
    };
    let native_subjects = NativeBlockSubjects {
        block: &db_block,
        metadata: database_metadata.as_ref(),
        self_aln_scores: &target_self_aln_scores,
        sequences: output_target_sequences.as_deref(),
        oids: edge_target_oids.as_deref(),
    };

    let mut text_writer = if capture.is_some() || edge_output.is_some() {
        Some(NativeBlockOutput::Discard(io::sink()))
    } else if daa_output {
        None
    } else {
        Some(NativeBlockOutput::new(config, false)?)
    };
    let mut daa_writer = if daa_output && capture.is_none() {
        let mut writer = NativeBlockOutput::new(config, true)?;
        init_daa(&mut writer)?;
        Some(writer)
    } else {
        None
    };
    let mut target_to_dict = vec![None; db_block.seqs().len()];
    let mut dict_targets = Vec::new();
    let mut aligned_queries = 0u64;
    let gapped_filter_evalue = setup.gapped_filter_evalue;
    let cutoff_gapped1 = if gapped_filter_evalue != 0.0 {
        Some(CutoffTable2D::new(
            &score_matrix,
            config.tuning.gapped_filter_evalue1,
        ))
    } else {
        None
    };
    let cutoff_gapped2 = if gapped_filter_evalue != 0.0 {
        Some(CutoffTable2D::new(&score_matrix, gapped_filter_evalue))
    } else {
        None
    };
    let gapped_filter_diag_score = setup.gapped_filter_diag_score;

    if let (Some(format), Some(writer)) = (tabular_format.as_ref(), text_writer.as_mut()) {
        writer.write_all(&format.header_bytes().map_err(io::Error::other)?)?;
    }
    if let (Some(format), Some(writer)) = (context_format, text_writer.as_mut()) {
        let first = if let Some(groups) = translated_groups.as_ref() {
            groups.first().map(|group| (group.id, group.dna_len))
        } else {
            query_records
                .first()
                .map(|record| (record.id.as_slice(), record.sequence.len() as i32))
        };
        if let Some((title, len)) = first {
            format.header(writer, config, &score_matrix, title, len)?;
        }
    }

    if let Some(capture) = capture {
        let count = translated_groups
            .as_ref()
            .map_or(query_records.len(), Vec::len);
        let mut output = capture
            .lock()
            .map_err(|_| io::Error::other("native search capture poisoned"))?;
        output.matches.resize_with(count, Vec::new);
        output.seeded.resize(count, false);
    }

    // `align_worker` skips a source query whose retained-hit iterator is
    // null before invoking the ordinary (non-frameshift) extension pipeline.
    // Only queries with seed hits can emit an unaligned output record.
    // Process each query in input order. Native blastp keeps lazy target
    // masking disabled, so the C++-style extension path can read the shared
    // target block without serializing all queries.
    let swipe_target_parallel = config.swipe_all
        && db_block.seqs().len()
            >= if config.translated_query_layout.is_some() {
                config
                    .translated_query_layout
                    .as_ref()
                    .unwrap()
                    .mate_indices
                    .len()
                    * 6
            } else {
                query_records.len()
            };
    // Source target_parallel serializes the outer query workers. Its SWIPE
    // heuristic additionally enables the DP parallel flag; ordinary searches
    // do not gain that flag merely from the numeric load-balancing option.
    let serial_query_extension = config.tuning.load_balancing == 1 || swipe_target_parallel;
    let extension_flags = if swipe_target_parallel {
        Flags::PARALLEL
    } else {
        Flags::NONE
    };
    let process_query = |(query_idx, query_rec, query_hits): (
        usize,
        &fasta::ByteFastaRecord,
        &[CompactHit],
    )|
     -> io::Result<QueryOutput> {
        if !config.swipe_all && query_hits.is_empty() {
            return Ok(QueryOutput::Skip);
        }
        let query = &query_rec.sequence;

        // CBS (composition-based statistics) correction per query position.
        // C++ default is comp-based-stats=1 (Hauser correction, window=40).
        // `Sequence::operator[]` strips the soft-mask bit while preserving
        // the residue, so compute Hauser from the query as stored rather
        // than converting masked positions to X.
        let query_cbs = if config.comp_based_stats.hauser() {
            crate::stats::cbs::hauser_correction_window(
                query,
                &score_matrix,
                config.tuning.cbs_window as u32 as usize,
            )
        } else {
            Vec::new()
        };
        let query_comp = crate::stats::cbs::compute_composition(query);
        let ungapped_cfg = UngappedStageConfig {
            log_extend: config.tuning.log_extend,
            chaining: config.tuning.chaining,
            cbs_matrix_scale: config.tuning.matrix_scale,
            hamming_ext,
            mutual_cover,
            lin_stage1_query: linear_stage1_query,
            lin_stage1_target: linear_stage1_target,
            comp_based_stats: config.comp_based_stats,
            anchored_swipe: config.tuning.setup.anchored_swipe,
            matrix_adjust_thresholds: config.tuning.matrix_adjust_thresholds,
            matrix_adjust_tolerance: config.tuning.cbs_err_tolerance,
            matrix_adjust_max_iterations: config.tuning.cbs_it_limit,
            xdrop: score_matrix.rawscore_int(config.ungapped_xdrop_bits),
            ..UngappedStageConfig::default()
        };
        let ext_mode = extension_mode_override.unwrap_or(setup.extension_mode);
        let gapped_cfg = GappedScoreConfig {
            rayon_task_pool: true,
            band_bin: config.tuning.band_bin,
            col_bin: config.tuning.col_bin,
            log_query: config.tuning.log_query,
            target_parallel_verbosity: config.tuning.target_parallel_verbosity,
            threads_align: config.tuning.threads_align,
            swipe_error: Some(Default::default()),
            add_self_aln: config.tuning.add_self_aln,
            clustering_filter: config.tuning.output.clustering_filter.clone(),
            clustering_query_self_score: query_self_aln_scores.get(query_idx).copied(),
            clustering_target_self_scores: clustering_target_self_scores.clone(),
            cluster_threshold_present: config
                .tuning
                .output
                .clustering_filter
                .as_ref()
                .is_some_and(|filter| filter.minimum.is_some()),
            cbs_matrix_scale: config.tuning.matrix_scale,
            matrix_adjust_tolerance: config.tuning.cbs_err_tolerance,
            matrix_adjust_max_iterations: config.tuning.cbs_it_limit,
            swipe_task_size: config.tuning.swipe_task_size,
            transcript_len_estimate: config.tuning.transcript_len_estimate,
            narrow_band_cov: config.tuning.narrow_band_cov,
            narrow_band_factor: config.tuning.narrow_band_factor,
            min_band_overlap: config.tuning.min_band_overlap,
            cutoff_score_8bit: config.tuning.cutoff_score_8bit,
            max_swipe_dp: config.tuning.max_swipe_dp,
            approx_backtrace: config.tuning.approx_backtrace,
            tantan_min_mask_prob: config.tuning.tantan_min_mask_prob,
            cbs_window: config.tuning.cbs_window,
            tantan_masker: config.tuning.runtime.tantan(),
            padding: config.tuning.padding.unwrap_or(0),
            no_ranking: config.tuning.no_ranking,
            mapany: config.tuning.mapany,
            target_hard_cap: (config.tuning.target_hard_cap != 0)
                .then_some(config.tuning.target_hard_cap),
            ranking_score_drop_factor: config.tuning.ranking_score_drop_factor,
            ranking_cutoff_bitscore: config.tuning.ranking_cutoff_bitscore,
            anchored_swipe: config.tuning.setup.anchored_swipe,
            inner_culling_overlap: config.culling_overlap,
            comp_based_stats_hauser: config.comp_based_stats.hauser(),
            comp_based_stats_matrix_adjust: config.comp_based_stats.matrix_adjust(),
            query_cover: config.query_cover,
            subject_cover: config.subject_cover,
            query_or_target_cover: config.query_or_subject_cover.max(query_or_target_cover),
            no_self_hits: config.no_self_hits,
            swipe_all: config.swipe_all,
            max_evalue: config.max_evalue,
            min_bit_score: config.min_bit_score,
            min_id: config.min_id,
            approx_min_id: approx_min_id,
            max_target_seqs: config.max_target_seqs,
            max_hsps: config.max_hsps,
            ext_chunk_size: config.ext_chunk_size,
            toppercent: config.toppercent,
            global_ranking_targets: config.global_ranking_targets,
            gapped_filter_evalue,
            sensitivity: config.sensitivity,
            self_: config.tuning.self_search,
            current_ref_block: config.tuning.search_block.reference,
            current_query_block: config.tuning.search_block.query,
            lin_stage1_query: linear_stage1_query,
            ..GappedScoreConfig::default()
        };

        let mut hits = Vec::with_capacity(query_hits.len());
        hits.extend(
            query_hits
                .iter()
                .map(|hit| Hit::with_score(0, hit.subject, hit.seed_offset, hit.score)),
        );
        let mut stat = Statistics::new();
        let output_hsp_values = if edge_output.is_some() {
            // Output::Format::Edge::hsp_values() requests coordinates only.
            HspValues::COORDS
        } else if let Some(format) = context_format {
            format.hsp_values()
        } else if daa_output {
            HspValues::TRANSCRIPT
        } else {
            tabular_hsp_values
        };
        let output_hsp_values = output_hsp_values
            | config
                .tuning
                .output
                .clustering_filter
                .as_ref()
                .map_or(HspValues::NONE, |filter| filter.format.hsp_values);
        let query_cbs_arg: &[Vec<i8>] = if query_cbs.is_empty() {
            &[]
        } else {
            std::slice::from_ref(&query_cbs)
        };
        let query_contexts = [query.as_slice()];
        let mut matches = extend_targets(
            query_idx as u32,
            &mut hits,
            &query_contexts,
            &query_rec.id,
            query.len() as i32,
            query_cbs_arg,
            &query_comp,
            &db_block,
            &mut stat,
            extension_flags,
            ext_mode,
            &gapped_cfg,
            &ungapped_cfg,
            &score_matrix,
            output_hsp_values,
            |query_len, target_len| {
                cutoff_gapped1
                    .as_ref()
                    .map_or(-1, |table| table.call(query_len, target_len))
            },
            |query_len, target_len| {
                cutoff_gapped2
                    .as_ref()
                    .map_or(-1, |table| table.call(query_len, target_len))
            },
            score_matrix.gap_open(),
            score_matrix.gap_extend(),
            gapped_filter_diag_score,
            config.tuning.gapped_filter_window,
            Option::<fn(u32, &crate::align::gapped_filter::SeedHitList) -> Vec<Match>>::None,
        );
        gapped_cfg.check_swipe_error()?;
        crate::basic::statistics::merge_statistics(&stat);
        validate_added_self_targets(&matches, db_block.seqs().len())?;
        if let Some(filter) = &config.tuning.output.clustering_filter {
            filter.check().map_err(io::Error::other)?;
        }
        if config.no_self_hits && ext_mode != sensitivity::ExtensionMode::None {
            // Seed-only mode returns before Source's ordinary HSP filters.
            // Upstream applies this HSP filter after ranking/extension. Keep
            // that placement: a discarded self target is not replaced by a
            // lower-ranked target from beyond the output range.
            matches.retain(|target| {
                let target_id = target.target_block_id as usize;
                target.is_synthetic_self_alignment()
                    || !same_sequence(query, &target.seq)
                    || db_ids.get(target_id) != query_rec.id.as_slice()
            });
        }
        if let Some(capture) = capture {
            let mut output = capture
                .lock()
                .map_err(|_| io::Error::other("native search capture poisoned"))?;
            let rows = matches.iter().map(|target| target.hsps.len() as u64).sum();
            output.matches[query_idx] = matches;
            output.seeded[query_idx] = config.swipe_all || !query_hits.is_empty();
            return Ok(QueryOutput::Tabular {
                bytes: Vec::new(),
                rows,
            });
        }
        let output_statistics = generated_output_statistics(&matches);
        if let Some(edge_output) = &edge_output {
            let mut edges = edge_output.lock().unwrap();
            for target in &matches {
                let Some(hsp) = target.hsps.first() else {
                    continue;
                };
                edges.push(InMemorySearchEdge {
                    edge: crate::output::edge::EdgeData {
                        query: edge_query_oids.as_ref().unwrap()[query_idx],
                        target: edge_target_oids.as_ref().unwrap()[target.target_block_id as usize],
                        qcovhsp: if query.is_empty() {
                            0.0
                        } else {
                            (hsp.query_source_range.length() as f64 * 100.0 / query.len() as f64)
                                as f32
                        },
                        scovhsp: {
                            let target_len =
                                db_block.seqs().length(target.target_block_id as usize);
                            if target_len == 0 {
                                0.0
                            } else {
                                (hsp.subject_range.length() as f64 * 100.0 / target_len as f64)
                                    as f32
                            }
                        },
                        // Output::Format::Edge serializes corrected bit score
                        // into this historical field; GVC consumes it as the
                        // edge weight (output_format.cpp:269).
                        evalue: hsp.corrected_bit_score,
                    },
                    approx_id: hsp.approx_id,
                    identity: if hsp.length == 0 {
                        0.0
                    } else {
                        100.0 * hsp.identities as f64 / hsp.length as f64
                    },
                });
            }
            crate::basic::statistics::merge_statistics(&output_statistics);
            return Ok(QueryOutput::Tabular {
                bytes: Vec::new(),
                rows: matches.len() as u64,
            });
        }
        if daa_output {
            return Ok(QueryOutput::Daa(matches));
        }
        if let Some(format) = context_format {
            let (bytes, rows) = format.query(
                query_idx as u32,
                edge_query_oids.as_ref().unwrap()[query_idx],
                &query_rec.id,
                query.len() as i32,
                &query_contexts,
                matches,
                &native_subjects,
                &score_matrix,
                false,
                loaded_db_letters,
                config.report_unaligned,
                &config.sam,
                &config.tuning.output,
                query_self_aln_scores.get(query_idx).copied(),
            )?;
            crate::basic::statistics::merge_statistics(&output_statistics);
            return Ok(QueryOutput::Tabular { bytes, rows });
        }
        let mut buf = Vec::new();
        let rows = super::search_output::write_tabular_query(
            &mut buf,
            tabular_format.as_ref().expect("tabular format"),
            query_idx as u32,
            edge_query_oids
                .as_ref()
                .map_or(query_idx as u64, |oids| oids[query_idx]),
            &query_rec.id,
            query,
            &query_contexts,
            None,
            query_rec.quality.as_deref().unwrap_or(b""),
            query_self_aln_scores.get(query_idx).copied(),
            matches,
            &native_subjects,
        )?;
        crate::basic::statistics::merge_statistics(&output_statistics);
        Ok(QueryOutput::Tabular { bytes: buf, rows })
    };

    let reference_taxon_rank_ids = config.reference_taxon_rank_ids.as_ref().map(|ranks| {
        (0..db_block.seqs().len())
            .map(|id| {
                let original = edge_target_oids
                    .as_ref()
                    .map_or(id, |oids| oids[id] as usize);
                ranks[original].clone()
            })
            .collect::<Vec<_>>()
    });
    let frameshift_backend = if config.frame_shift != 0 {
        Some(
            crate::align::legacy::native::NativeBackend::new_with_query_sequences(
                translated_query_block.as_ref().ok_or_else(|| {
                    io::Error::other("Frameshift extension requires translated query metadata.")
                })?,
                &stage2_query_set,
                &db_block,
                &score_matrix,
                crate::align::legacy::native::NativeConfig {
                    band_bin: config.tuning.band_bin,
                    col_bin: config.tuning.col_bin,
                    swipe_chunk_size: config.tuning.swipe_chunk_size,
                    transcript_len_estimate: config.tuning.transcript_len_estimate,
                    cbs_matrix_scale: config.tuning.matrix_scale,
                    mapper: crate::align::legacy::query_mapper::QueryMapperConfig {
                        log_query: config.tuning.log_query,
                        cbs_window: config.tuning.cbs_window as u32 as usize,
                        query_contexts: 6,
                        query_translated: true,
                        hauser: config.comp_based_stats.hauser(),
                        taxon_k: config.taxon_k.min(u32::MAX as u64) as u32,
                        culling: crate::output::target_culling::TargetCullingConfig {
                            query_range_culling: config.range_culling,
                            query_range_cover: config.range_cover,
                            taxon_k: config.taxon_k.min(u32::MAX as u64) as u32,
                            toppercent: config.toppercent,
                            ..Default::default()
                        },
                        toppercent: config.toppercent,
                        min_id: config.min_id,
                        query_cover: config.query_cover,
                        subject_cover: config.subject_cover,
                        max_hsps: config.max_hsps,
                        ..Default::default()
                    },
                    pipeline: crate::align::legacy::banded_swipe_pipeline::PipelineConfig {
                        frame_shift: true,
                        rank_ratio: (config.tuning.rank_ratio != -1.0)
                            .then_some(config.tuning.rank_ratio),
                        rank_factor: (config.tuning.rank_factor != -1.0)
                            .then_some(config.tuning.rank_factor),
                        padding: config.tuning.padding,
                        query_range_culling: config.range_culling,
                        threads: if config.threads > 0 {
                            config.threads as usize
                        } else {
                            rayon::current_num_threads().max(1)
                        },
                        ..Default::default()
                    },
                    max_target_seqs: config.max_target_seqs,
                    min_bit_score: config.min_bit_score,
                    max_evalue: config.max_evalue,
                    ungapped_xdrop_bits: config.ungapped_xdrop_bits,
                    discard_output: context_format == Some(SearchFormat::Null),
                    ..Default::default()
                },
                reference_taxon_rank_ids.as_deref(),
            )?,
        )
    } else {
        None
    };

    let process_translated_query =
        |(query_idx, group, query_hits): (usize, &TranslatedExtensionGroup, &[Vec<CompactHit>]),
         scratch: &mut ExtensionWorkerScratch|
         -> io::Result<QueryOutput> {
            if config.frame_shift == 0 && !config.swipe_all && query_hits.iter().all(Vec::is_empty)
            {
                return Ok(QueryOutput::Skip);
            }
            debug_assert_eq!(query_hits.len(), 6);
            let matches = if let Some(backend) = frameshift_backend.as_ref() {
                let mut source_hits = Vec::new();
                for (frame, hits) in query_hits.iter().enumerate() {
                    source_hits.extend(hits.iter().map(|hit| {
                        Hit::with_score(
                            query_idx as u32 * 6 + frame as u32,
                            hit.subject,
                            hit.seed_offset,
                            hit.score,
                        )
                    }));
                }
                let mut statistics = Statistics::new();
                let mut dp_stat = crate::dp::banded_3frame::DpStat::default();
                let result = backend.extend_query_for_capture(
                    query_idx as u32,
                    source_hits,
                    &mut statistics,
                    &mut dp_stat,
                )?;
                crate::basic::statistics::merge_statistics(&statistics);
                crate::dp::banded_3frame::merge_dp_stat(&dp_stat);
                result.matches
            } else {
                let ExtensionWorkerScratch {
                    hits,
                    query_cbs,
                    hauser_values,
                } = scratch;
                let query_cbs: &[Vec<i8>] = if config.comp_based_stats.hauser() {
                    for (&query, output) in group.frames.iter().zip(query_cbs.iter_mut()) {
                        crate::stats::hauser_correction::hauser_correction_into_with_window(
                            query,
                            &score_matrix,
                            config.tuning.cbs_window as u32 as usize,
                            output,
                            hauser_values,
                        );
                    }
                    query_cbs
                } else {
                    &[]
                };
                let query_comp = crate::stats::cbs::compute_composition(group.frames[0]);
                let ungapped_cfg = UngappedStageConfig {
                    log_extend: config.tuning.log_extend,
                    chaining: config.tuning.chaining,
                    cbs_matrix_scale: config.tuning.matrix_scale,
                    query_contexts: 6,
                    query_translated: true,
                    hamming_ext,
                    mutual_cover,
                    lin_stage1_query: linear_stage1_query,
                    lin_stage1_target: linear_stage1_target,
                    comp_based_stats: config.comp_based_stats,
                    anchored_swipe: config.tuning.setup.anchored_swipe,
                    matrix_adjust_thresholds: config.tuning.matrix_adjust_thresholds,
                    matrix_adjust_tolerance: config.tuning.cbs_err_tolerance,
                    matrix_adjust_max_iterations: config.tuning.cbs_it_limit,
                    xdrop: score_matrix.rawscore_int(config.ungapped_xdrop_bits),
                    ..UngappedStageConfig::default()
                };
                let ext_mode = extension_mode_override.unwrap_or(setup.extension_mode);
                let gapped_cfg = GappedScoreConfig {
                    rayon_task_pool: true,
                    band_bin: config.tuning.band_bin,
                    col_bin: config.tuning.col_bin,
                    narrow_band_cov: config.tuning.narrow_band_cov,
                    narrow_band_factor: config.tuning.narrow_band_factor,
                    min_band_overlap: config.tuning.min_band_overlap,
                    log_query: config.tuning.log_query,
                    target_parallel_verbosity: config.tuning.target_parallel_verbosity,
                    threads_align: config.tuning.threads_align,
                    swipe_error: Some(Default::default()),
                    add_self_aln: config.tuning.add_self_aln,
                    clustering_filter: config.tuning.output.clustering_filter.clone(),
                    clustering_query_self_score: query_self_aln_scores.get(query_idx * 6).copied(),
                    clustering_target_self_scores: clustering_target_self_scores.clone(),
                    cluster_threshold_present: config
                        .tuning
                        .output
                        .clustering_filter
                        .as_ref()
                        .is_some_and(|filter| filter.minimum.is_some()),
                    cbs_matrix_scale: config.tuning.matrix_scale,
                    matrix_adjust_tolerance: config.tuning.cbs_err_tolerance,
                    matrix_adjust_max_iterations: config.tuning.cbs_it_limit,
                    query_contexts: 6,
                    query_translated: true,
                    swipe_task_size: config.tuning.swipe_task_size,
                    transcript_len_estimate: config.tuning.transcript_len_estimate,
                    cutoff_score_8bit: config.tuning.cutoff_score_8bit,
                    max_swipe_dp: config.tuning.max_swipe_dp,
                    approx_backtrace: config.tuning.approx_backtrace,
                    tantan_min_mask_prob: config.tuning.tantan_min_mask_prob,
                    cbs_window: config.tuning.cbs_window,
                    tantan_masker: config.tuning.runtime.tantan(),
                    padding: config.tuning.padding.unwrap_or(0),
                    no_ranking: config.tuning.no_ranking,
                    mapany: config.tuning.mapany,
                    target_hard_cap: (config.tuning.target_hard_cap != 0)
                        .then_some(config.tuning.target_hard_cap),
                    ranking_score_drop_factor: config.tuning.ranking_score_drop_factor,
                    ranking_cutoff_bitscore: config.tuning.ranking_cutoff_bitscore,
                    anchored_swipe: config.tuning.setup.anchored_swipe,
                    inner_culling_overlap: config.culling_overlap,
                    comp_based_stats_hauser: config.comp_based_stats.hauser(),
                    comp_based_stats_matrix_adjust: config.comp_based_stats.matrix_adjust(),
                    query_cover: config.query_cover,
                    subject_cover: config.subject_cover,
                    query_or_target_cover: config.query_or_subject_cover.max(query_or_target_cover),
                    no_self_hits: config.no_self_hits,
                    swipe_all: config.swipe_all,
                    max_evalue: config.max_evalue,
                    min_bit_score: config.min_bit_score,
                    min_id: config.min_id,
                    approx_min_id: approx_min_id,
                    max_target_seqs: config.max_target_seqs,
                    max_hsps: config.max_hsps,
                    ext_chunk_size: config.ext_chunk_size,
                    toppercent: config.toppercent,
                    global_ranking_targets: config.global_ranking_targets,
                    gapped_filter_evalue,
                    sensitivity: config.sensitivity,
                    self_: config.tuning.self_search,
                    current_ref_block: config.tuning.search_block.reference,
                    current_query_block: config.tuning.search_block.query,
                    lin_stage1_query: linear_stage1_query,
                    ..GappedScoreConfig::default()
                };

                let total_hits = query_hits.iter().map(Vec::len).sum();
                hits.clear();
                hits.reserve(total_hits);
                for (frame, frame_hits) in query_hits.iter().enumerate() {
                    hits.extend(frame_hits.iter().map(|hit| {
                        Hit::with_score(frame as u32, hit.subject, hit.seed_offset, hit.score)
                    }));
                }
                let mut stat = Statistics::new();
                let output_hsp_values = if let Some(format) = context_format {
                    format.hsp_values()
                } else if daa_output {
                    HspValues::TRANSCRIPT
                } else {
                    tabular_hsp_values
                };
                let output_hsp_values = output_hsp_values
                    | config
                        .tuning
                        .output
                        .clustering_filter
                        .as_ref()
                        .map_or(HspValues::NONE, |filter| filter.format.hsp_values);
                let matches = extend_targets(
                query_idx as u32,
                hits,
                &group.frames,
                group.id,
                group.dna_len,
                query_cbs,
                &query_comp,
                &db_block,
                &mut stat,
                extension_flags,
                ext_mode,
                &gapped_cfg,
                &ungapped_cfg,
                &score_matrix,
                output_hsp_values,
                |query_len, target_len| {
                    cutoff_gapped1
                        .as_ref()
                        .map_or(-1, |table| table.call(query_len, target_len))
                },
                |query_len, target_len| {
                    cutoff_gapped2
                        .as_ref()
                        .map_or(-1, |table| table.call(query_len, target_len))
                },
                score_matrix.gap_open(),
                score_matrix.gap_extend(),
                gapped_filter_diag_score,
                config.tuning.gapped_filter_window,
                Option::<fn(u32, &crate::align::gapped_filter::SeedHitList) -> Vec<Match>>::None,
            );
                gapped_cfg.check_swipe_error()?;
                crate::basic::statistics::merge_statistics(&stat);
                validate_added_self_targets(&matches, db_block.seqs().len())?;
                if let Some(filter) = &config.tuning.output.clustering_filter {
                    filter.check().map_err(io::Error::other)?;
                }
                matches
            };
            if let Some(capture) = capture {
                let mut output = capture
                    .lock()
                    .map_err(|_| io::Error::other("native search capture poisoned"))?;
                output.matches[query_idx] = matches;
                output.seeded[query_idx] = config.frame_shift != 0
                    || config.swipe_all
                    || query_hits.iter().any(|hits| !hits.is_empty());
                return Ok(QueryOutput::Tabular {
                    bytes: Vec::new(),
                    rows: 0,
                });
            }
            let output_statistics = generated_output_statistics(&matches);
            if daa_output {
                return Ok(QueryOutput::Daa(matches));
            }
            if let Some(format) = context_format {
                let (bytes, rows) = format.query(
                    query_idx as u32,
                    query_idx as u64,
                    group.id,
                    group.dna_len,
                    &group.frames,
                    matches,
                    &native_subjects,
                    &score_matrix,
                    true,
                    loaded_db_letters,
                    config.report_unaligned,
                    &config.sam,
                    &config.tuning.output,
                    query_self_aln_scores.get(query_idx).copied(),
                )?;
                crate::basic::statistics::merge_statistics(&output_statistics);
                return Ok(QueryOutput::Tabular { bytes, rows });
            }
            let mut buf = Vec::new();
            let rows = super::search_output::write_tabular_query(
                &mut buf,
                tabular_format.as_ref().expect("tabular format"),
                query_idx as u32,
                edge_query_oids
                    .as_ref()
                    .map_or(query_idx as u64, |oids| oids[query_idx]),
                group.id,
                group.source_sequence,
                &group.frames,
                group.mate_sequence,
                group.quality,
                query_self_aln_scores.get(query_idx).copied(),
                matches,
                &native_subjects,
            )?;
            crate::basic::statistics::merge_statistics(&output_statistics);
            Ok(QueryOutput::Tabular { bytes: buf, rows })
        };
    let heartbeat = config.tuning.heartbeat
        && config.tuning.load_balancing == 0
        && !config.swipe_all
        && crate::util::log_stream::log_stream()
            .lock()
            .unwrap()
            .to_cout();
    if config.tuning.min_task_trace_pts <= 0 {
        return Err(io::Error::other("Source alignment partitioning requires a positive min_task_trace_pts; zero or negative values enter undefined partition arithmetic."));
    }
    let align_pool = if !serial_query_extension && config.tuning.threads_align != 0 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(config.tuning.threads_align as usize)
                .build()
                .map_err(io::Error::other)?,
        )
    } else {
        None
    };
    let resident_query_bytes = translated_query_block.as_ref().map_or_else(
        || {
            query_records
                .iter()
                .map(|record| {
                    (record.sequence.capacity()
                        + record.id.capacity()
                        + record.quality.as_ref().map_or(0, Vec::capacity))
                        as u64
                })
                .sum::<u64>()
        },
        |block| block.mem_size() as u64,
    );
    let hit_resident_bytes = resident_query_bytes.wrapping_add(db_block.mem_size() as u64);
    let alignment_config = crate::align::align::AlignConfig {
        min_task_trace_pts: usize::try_from(config.tuning.min_task_trace_pts)
            .map_err(io::Error::other)?,
        query_contexts: if translated_groups.is_some() { 6 } else { 1 },
        swipe_all: config.swipe_all,
        memory_limit: config.memory_limit.unwrap_or(16 * 1024 * 1024 * 1024),
        trace_pt_fetch_size: config.tuning.trace_pt_fetch_size,
        hit_record_size: 15,
        ..Default::default()
    };
    let mut alignment_state = crate::align::align::AlignState::new(0, 0);
    let mut total_alignments = 0u64;
    if let Some(groups) = translated_groups.as_ref() {
        let process =
            |query: usize, hits: &[Vec<CompactHit>], scratch: &mut ExtensionWorkerScratch| {
                process_translated_query((query, &groups[query], hits), scratch)
            };
        let emit = |query: usize, output: QueryOutput| -> io::Result<()> {
            match output {
                QueryOutput::Skip => (),
                QueryOutput::Tabular { bytes, rows } => {
                    total_alignments += rows;
                    text_writer
                        .as_mut()
                        .expect("tabular writer")
                        .write_all(&bytes)?;
                }
                QueryOutput::Daa(matches) => {
                    let output_statistics = generated_output_statistics(&matches);
                    let group = &groups[query];
                    let mut bytes = Vec::new();
                    let count = encode_daa_query(
                        &mut bytes,
                        group.id,
                        group.source_sequence,
                        SequenceType::Nucleotide,
                        &matches,
                        &mut target_to_dict,
                        &mut dict_targets,
                    )?;
                    crate::basic::statistics::merge_statistics(&output_statistics);
                    if count != 0 {
                        aligned_queries += 1;
                        total_alignments += count;
                        daa_writer.as_mut().expect("DAA writer").write_all(&bytes)?;
                    }
                }
            }
            Ok(())
        };
        let source = NativeGroupedHitSource {
            store: &mut hit_store,
            next_memory_bin: 0,
            pending_memory_bin: None,
        };
        let mut backend = crate::align::native::NativeAlignBackend::new(
            source,
            groups,
            6,
            hit_resident_bytes,
            align_pool.as_ref(),
            serial_query_extension,
            heartbeat,
            process,
            emit,
        );
        crate::align::align::align_queries(&alignment_config, &mut backend, &mut alignment_state)
            .map_err(io::Error::other)?;
    } else {
        let process = |query: usize, hits: &[Vec<CompactHit>], _: &mut ()| {
            process_query((query, &query_records[query], &hits[0]))
        };
        let emit = |query: usize, output: QueryOutput| -> io::Result<()> {
            let record = &query_records[query];
            match output {
                QueryOutput::Skip => (),
                QueryOutput::Tabular { bytes, rows } => {
                    total_alignments += rows;
                    text_writer
                        .as_mut()
                        .expect("tabular writer")
                        .write_all(&bytes)?;
                }
                QueryOutput::Daa(matches) => {
                    let output_statistics = generated_output_statistics(&matches);
                    let mut bytes = Vec::new();
                    let count = encode_daa_query(
                        &mut bytes,
                        &record.id,
                        &record.sequence,
                        SequenceType::AminoAcid,
                        &matches,
                        &mut target_to_dict,
                        &mut dict_targets,
                    )?;
                    crate::basic::statistics::merge_statistics(&output_statistics);
                    if count != 0 {
                        aligned_queries += 1;
                        total_alignments += count;
                        daa_writer.as_mut().expect("DAA writer").write_all(&bytes)?;
                    }
                }
            }
            Ok(())
        };
        let source = NativeGroupedHitSource {
            store: &mut hit_store,
            next_memory_bin: 0,
            pending_memory_bin: None,
        };
        let mut backend = crate::align::native::NativeAlignBackend::new(
            source,
            &query_records,
            1,
            hit_resident_bytes,
            align_pool.as_ref(),
            serial_query_extension,
            heartbeat,
            process,
            emit,
        );
        crate::align::align::align_queries(&alignment_config, &mut backend, &mut alignment_state)
            .map_err(io::Error::other)?;
    }
    // Source disk usage is a maximum across reference blocks, not a sum.
    crate::basic::statistics::max_statistics(
        crate::basic::statistics::StatValue::SearchTempSpace,
        alignment_state
            .statistics
            .get(crate::basic::statistics::StatValue::SearchTempSpace),
    );
    alignment_state.statistics.data
        [crate::basic::statistics::StatValue::SearchTempSpace as usize] = 0;
    crate::basic::statistics::merge_statistics(&alignment_state.statistics);
    if let Some(writer) = text_writer.as_mut() {
        if let Some(format) = context_format {
            format.footer(writer)?;
        } else if let Some(format) = &tabular_format {
            writer.write_all(format.footer().as_bytes())?;
        }
        writer.close()?;
    }
    if let Some(writer) = daa_writer.as_mut() {
        let mut refs = Vec::with_capacity(dict_targets.len());
        for &target_id in &dict_targets {
            let title = db_ids.get(target_id);
            refs.push((
                config.archive.dictionary_title_bytes(&title),
                db_block.seqs().length(target_id) as u32,
            ));
        }
        let dictionary = NativeDaaDictionary {
            sequence_count: db_block.seqs().len() as u64,
            refs,
        };
        let metadata = DaaRunMetadata {
            db_letters: score_matrix.db_letters(),
            gap_open: score_matrix.gap_open(),
            gap_extend: score_matrix.gap_extend(),
            // Protein/translated archives do not use nucleotide reward or
            // penalty. Upstream reads separate uninitialized Config fields
            // here (not match_reward/mismatch_penalty); use deterministic
            // zeroes, matching the pinned comparison executable's output.
            reward: 0,
            penalty: 0,
            k: score_matrix.k(),
            lambda: score_matrix.lambda(),
            max_evalue: config.max_evalue,
            matrix: score_matrix.name().as_bytes().to_vec(),
            mode: if translated_groups.is_some() {
                crate::basic::value::AlignMode::BLASTX as u32
            } else {
                crate::basic::value::AlignMode::BLASTP as u32
            },
            aligned_queries,
        };
        finish_daa_from_sequence_file_bytes(writer, &dictionary, &metadata)?;
        writer.close()?;
    }

    drop(frameshift_backend);
    if let Some(capture) = capture {
        let mut output = capture
            .lock()
            .map_err(|_| io::Error::other("native search capture poisoned"))?;
        output.query_self_aln_scores = query_self_aln_scores;
        output.target_self_aln_scores = target_self_aln_scores;
        output.queries = query_records;
        output.query_oids = edge_query_oids.unwrap_or_default();
        output.target_oids = edge_target_oids.unwrap_or_default();
        if let Some(block) = translated_query_block.as_mut() {
            block.put_seqs(stage2_query_set);
        }
        output.translated = translated_query_block;
        output.matrix = Some(score_matrix);
    }

    let elapsed = start.elapsed();
    search_message!(
        "Reported {} alignments in {:.1}s",
        total_alignments,
        elapsed.as_secs_f64()
    );
    Ok(())
}

// Public UTF-8 record conveniences share the canonical byte search pipeline.
pub fn run_with_query_records(
    config: &BlastpConfig,
    query_records: Vec<fasta::FastaRecord>,
) -> io::Result<()> {
    run_with_query_record_bytes(config, query_records.into_iter().map(Into::into).collect())
}

pub fn run_with_records(
    config: &BlastpConfig,
    database: Vec<fasta::FastaRecord>,
    queries: Vec<fasta::FastaRecord>,
) -> io::Result<()> {
    run_with_record_bytes(
        config,
        database.into_iter().map(Into::into).collect(),
        queries.into_iter().map(Into::into).collect(),
    )
}

pub fn run_edges_in_memory(
    config: &BlastpConfig,
    database: Vec<fasta::FastaRecord>,
    queries: Option<Vec<fasta::FastaRecord>>,
    approx_min_id: f64,
    linear_stage1_query: bool,
    linear_stage1_target: bool,
    self_search: bool,
    query_or_target_cover: f64,
    mutual_cover: bool,
    hamming_ext: HammingExtConfig,
    soft_tantan_masker: Option<&crate::masking::tantan::TantanMasker>,
    soft_masking_override: MaskingAlgo,
    freq_masking: bool,
    effective_db_letters: Option<u64>,
    extension_mode_override: Option<sensitivity::ExtensionMode>,
) -> io::Result<Vec<InMemorySearchEdge>> {
    run_edges_in_memory_bytes(
        config,
        database.into_iter().map(Into::into).collect(),
        queries.map(|records| records.into_iter().map(Into::into).collect()),
        approx_min_id,
        linear_stage1_query,
        linear_stage1_target,
        self_search,
        query_or_target_cover,
        mutual_cover,
        hamming_ext,
        soft_tantan_masker,
        soft_masking_override,
        freq_masking,
        effective_db_letters,
        extension_mode_override,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::{MASK_LETTER, SEED_MASK};

    #[test]
    fn joined_renderer_preserves_fetched_query_after_all_targets_are_culled() {
        let args = ["diamond", "blastp", "-q", "query.faa", "-d", "database"].map(String::from);
        let config =
            crate::commands::dispatch::parse_native_blastp_config(&args, "diamond blastp").unwrap();
        let matrix = ScoreMatrix::new("BLOSUM62", 11, 1, 0, 1, 100).unwrap();
        let ids = std::collections::HashMap::new();
        let renderer = NativeQueryRenderer {
            config: &config,
            query_mode: 2,
            context_format: None,
            tabular: None,
            daa: true,
            ref_ids: &ids,
            database_letters: 100,
        };
        let subjects = Block::new();
        let source = [0, 1, 2];
        let rendered = renderer
            .render_query_bytes(
                0,
                0,
                b"query\xff",
                &source,
                &[&source],
                None,
                &[],
                None,
                Vec::new(),
                &subjects,
                &matrix,
                true,
            )
            .unwrap();
        let mut expected = Vec::new();
        let begin = write_daa_query_record_bytes(
            &mut expected,
            b"query\xff",
            &source,
            SequenceType::AminoAcid,
        );
        finish_daa_query_record(&mut expected, begin);
        assert_eq!(rendered.bytes, expected);
        assert_eq!(rendered.aligned_queries, 1);
        let unaligned = renderer
            .render_query_bytes(
                0,
                0,
                b"query\xff",
                &source,
                &[&source],
                None,
                &[],
                None,
                Vec::new(),
                &subjects,
                &matrix,
                false,
            )
            .unwrap();
        assert!(unaligned.bytes.is_empty());
        assert_eq!(unaligned.aligned_queries, 0);
    }

    fn taxonomy_test_dir() -> PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-native-taxonomy-search-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    #[test]
    fn native_default_block_budget_matches_source_run_defaults() {
        let mut config = super::super::blastx::blastp_config(
            &super::super::blastx::BlastxConfig {
                query_files: Vec::new(),
                database: std::path::PathBuf::new(),
                output: None,
                matrix: "blosum62".into(),
                gap_open: 11,
                gap_extend: 1,
                max_evalue: 0.001,
                min_bit_score: 0.0,
                max_target_seqs: 25,
                max_hsps: 1,
                report_unaligned: -1,
                include_lineage: false,
                output_header: None,
                compression: Vec::new(),
                file_buffer_size: 67_108_864,
                sam: Default::default(),
                archive: Default::default(),
                tuning: Default::default(),
                scoring_matrix: None,
                auxiliary_output: Default::default(),
                taxonlist: Vec::new(),
                taxon_exclude: Vec::new(),
                seqidlist: std::path::PathBuf::new(),
                skip_missing_seqids: false,
                range_culling: false,
                range_cover: 50.0,
                culling_overlap: 50.0,
                taxon_k: 0,
                query_or_subject_cover: 0.0,
                reference_taxon_rank_ids: None,
                ext_chunk_size: 0,
                toppercent: None,
                global_ranking_targets: 0,
                frame_shift: 0,
                extension_mode: None,
                iterate: None,
                min_id: 0.0,
                approx_min_id: 0.0,
                soft_masking: crate::masking::MaskingAlgo::None,
                threads: 1,
                index_chunks: 0,
                block_size: None,
                query_indexed: false,
                contiguous_seed_mode: false,
                min_query_len: 0,
                swipe_all: false,
                outfmt: vec!["6".into()],
                sensitivity: Sensitivity::Default,
                query_gencode: 1,
                strand: "both".into(),
                min_orf: None,
                masking: MaskingMode::None,
                motif_masking: Vec::new(),
                query_cover: 0.0,
                subject_cover: 0.0,
                comp_based_stats: CbsMode::Hauser,
                no_self_hits: false,
                ungapped_xdrop_bits: 12.3,
                memory_limit: None,
                tmpdir: PathBuf::new(),
                target_indexed: false,
            },
            Vec::new(),
        );
        for sensitivity in [
            Sensitivity::Faster,
            Sensitivity::Fast,
            Sensitivity::Default,
            Sensitivity::Linclust40,
            Sensitivity::Linclust20,
            Sensitivity::Shapes6x10,
            Sensitivity::Shapes30x10,
            Sensitivity::MidSensitive,
            Sensitivity::Sensitive,
            Sensitivity::MoreSensitive,
            Sensitivity::VerySensitive,
            Sensitivity::UltraSensitive,
        ] {
            config.sensitivity = sensitivity;
            let expected = if sensitivity >= Sensitivity::VerySensitive {
                400_000_000
            } else {
                2_000_000_000
            };
            assert_eq!(native_block_budget(&config).unwrap(), Some(expected));
            config.block_size = Some(0.0);
            assert_eq!(native_block_budget(&config).unwrap(), Some(expected));
            config.block_size = None;
        }
    }

    #[test]
    fn native_dmnd_taxonomy_metadata_reaches_all_supported_output_fields() {
        let dir = taxonomy_test_dir();
        std::fs::create_dir(&dir).unwrap();
        let fasta = dir.join("subjects.faa");
        let mapping = dir.join("mapping.tsv");
        let nodes = dir.join("nodes.dmp");
        let names = dir.join("names.dmp");
        let dmnd = dir.join("subjects.dmnd");
        std::fs::write(
            &fasta,
            b">sp|P00001.1|one >tr|P00002.1|two\nMKTAYIAK\n>UNMAPPED.1\nVLSPADKT\n",
        )
        .unwrap();
        std::fs::write(
            &mapping,
            b"accession.version\ttaxid\nP00001\t10\nP00002\t11\n",
        )
        .unwrap();
        std::fs::write(
            &nodes,
            b"1\t|\t1\t|\tno rank\t|\t\n2\t|\t1\t|\tsuperkingdom\t|\t\n3\t|\t2\t|\tkingdom\t|\t\n4\t|\t3\t|\tphylum\t|\t\n10\t|\t4\t|\tspecies\t|\t\n11\t|\t4\t|\tspecies\t|\t\n",
        )
        .unwrap();
        std::fs::write(
            &names,
            b"1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\t\t|\tscientific name\t|\n3\t|\tTest kingdom\t|\t\t|\tscientific name\t|\n4\t|\tTest phylum\t|\t\t|\tscientific name\t|\n10\t|\tAlpha species\t|\t\t|\tscientific name\t|\n11\t|\tBeta species\t|\t\t|\tscientific name\t|\n",
        )
        .unwrap();
        crate::data::db_builder::build_db_with_options(
            &[fasta.to_str().unwrap()],
            dmnd.to_str().unwrap(),
            SequenceType::AminoAcid,
            &crate::data::db_builder::BuildOptions {
                taxonmap: Some(mapping.clone()),
                taxonnodes: Some(nodes.clone()),
                taxonnames: Some(names.clone()),
                no_parse_seqids: false,
                file_buffer_size: None,
            },
        )
        .unwrap();

        let flags = requested_taxonomy_flags(&[
            "6".into(),
            "staxids".into(),
            "sscinames".into(),
            "sskingdoms".into(),
            "skingdoms".into(),
            "sphylums".into(),
            "slineages".into(),
        ]);
        let (header, records, metadata) = load_dmnd_database(&dmnd, flags).unwrap();
        assert_eq!(header.sequences, 2);
        assert_eq!(records.len(), 2);
        assert_eq!(metadata.subject_taxids, vec![vec![10, 11], vec![]]);
        let tree = metadata.taxonomy.as_ref().unwrap();
        assert_eq!(
            format::print_staxids(&metadata.subject_taxids[0], false),
            "10;11"
        );
        assert_eq!(
            format::print_taxon_names(metadata.subject_taxids[0].iter().copied(), tree, false),
            "Alpha species;Beta species"
        );
        assert_eq!(
            format::print_rank_taxon_names(
                &metadata.subject_taxids[0],
                tree,
                "superkingdom",
                false
            ),
            "Bacteria"
        );
        assert_eq!(
            format::print_rank_taxon_names(&metadata.subject_taxids[0], tree, "kingdom", false),
            "Test kingdom"
        );
        assert_eq!(
            format::print_rank_taxon_names(&metadata.subject_taxids[0], tree, "phylum", false),
            "Test phylum"
        );
        assert_eq!(
            format::print_lineage(&metadata.subject_taxids[0], tree, false),
            "Bacteria; Test kingdom; Test phylum; Alpha species<>Bacteria; Test kingdom; Test phylum; Beta species"
        );
        assert_eq!(
            format::print_staxids(&metadata.subject_taxids[1], false),
            ""
        );
        assert_eq!(
            format::print_taxon_names(metadata.subject_taxids[1].iter().copied(), tree, false),
            "N/A"
        );
        assert_eq!(
            format::print_lineage(&metadata.subject_taxids[1], tree, false),
            "N/A"
        );

        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn packed_loc_double_array_ranges_are_borrowed_without_copying() {
        let expected = [0x12_3456_789au64, 0x01_0203_0405u64];
        // Start at an intentionally odd byte address. PackedLoc's alignment
        // is one, just like a DoubleArray range following its u32 count.
        let mut storage = vec![0xa5];
        for value in expected {
            storage.push((value >> 32) as u8);
            storage.extend_from_slice(&(value as u32).to_ne_bytes());
        }
        let bytes = &storage[1..];
        let locations = packed_locs_from_bytes(bytes);
        assert_eq!(locations.as_ptr().cast::<u8>(), bytes.as_ptr());
        assert_eq!(locations.len(), expected.len());
        assert_eq!(locations[0].as_u64(), expected[0]);
        assert_eq!(locations[1].as_u64(), expected[1]);
    }

    #[test]
    fn borrowed_stage1_targets_preserve_tile_offsets_and_hit_order() {
        let targets = [
            PackedLoc::new(10),
            PackedLoc::new(20),
            PackedLoc::new(30),
            PackedLoc::new(0x01_0000_0040),
        ];
        let borrowed = BorrowedStage1Targets {
            indices: &[2, 0],
            target_begin: 1,
            targets: &targets,
        };
        assert_eq!(borrowed.len(), 2);
        assert_eq!(borrowed.absolute_subject(0), 0x01_0000_0040);
        assert_eq!(borrowed.absolute_subject(1), 20);
    }

    #[test]
    fn borrowed_stage1_subjects_preserve_hit_field_order() {
        let subjects = [10, 20, 30, 40];
        let borrowed = BorrowedStage1Subjects {
            indices: &[3, 1, 2],
            subjects: &subjects,
        };
        assert_eq!(borrowed.len(), 3);
        assert_eq!(borrowed.absolute_subject(0), 40);
        assert_eq!(borrowed.absolute_subject(1), 20);
        assert_eq!(borrowed.absolute_subject(2), 30);
    }

    #[test]
    #[should_panic]
    fn packed_loc_double_array_ranges_reject_partial_elements() {
        let _ = packed_locs_from_bytes(&[0; 6]);
    }

    #[test]
    fn self_hit_sequence_comparison_ignores_mask_annotations() {
        assert!(same_sequence(&[1, 2, 3], &[1 | SEED_MASK, 2, 3]));
        assert!(!same_sequence(&[1, 2, 3], &[1, 2, 4]));
        assert!(!same_sequence(&[1, 2], &[1, 2, 3]));
    }

    #[test]
    fn extension_queries_use_restored_soft_masking_state() {
        let mut records: Vec<fasta::FastaRecord> = vec![
            fasta::FastaRecord {
                quality: None,
                id: "q0".into(),
                sequence: vec![MASK_LETTER; 4],
            },
            fasta::FastaRecord {
                quality: None,
                id: "q1".into(),
                sequence: vec![MASK_LETTER; 2],
            },
        ];
        let mut restored = SequenceSet::new();
        restored.push(&[1, 2 | SEED_MASK, 3, 4]);
        restored.push(&[5, 6]);

        restore_extension_queries(&mut records, &restored);

        assert_eq!(records[0].sequence, [1, 2 | SEED_MASK, 3, 4]);
        assert_eq!(records[1].sequence, [5, 6]);
    }

    fn context_hits(counts: &[usize], contexts: usize) -> Vec<Vec<CompactHit>> {
        let hit = CompactHit {
            subject: 0,
            seed_offset: 0,
            score: 0,
        };
        counts
            .iter()
            .flat_map(|&count| {
                let mut source = vec![Vec::new(); contexts];
                source[contexts - 1] = vec![hit; count];
                source
            })
            .collect()
    }

    #[test]
    fn translated_partitions_cross_threshold_only_at_source_boundaries() {
        let hits = context_hits(&[0, 700, 400, 1_500, 20], 6);
        assert_eq!(
            translated_hit_partitions(&hits, 5, 6, 1_024),
            vec![0..3, 3..4, 4..5]
        );
    }

    #[test]
    fn translated_partitions_retain_leading_and_trailing_zero_hit_queries() {
        let hits = context_hits(&[0, 0, 1_024, 0, 0], 6);
        assert_eq!(translated_hit_partitions(&hits, 5, 6, 1_024), vec![0..5]);
        let crossed = context_hits(&[1_024, 1, 0], 6);
        assert_eq!(translated_hit_partitions(&crossed, 3, 6, 1_024), vec![0..3]);
        let empty = context_hits(&[0, 0, 0], 6);
        assert_eq!(translated_hit_partitions(&empty, 3, 6, 1_024), vec![0..3]);
    }

    #[test]
    fn translated_partitions_match_upstream_hit_index_algorithm() {
        fn reference(counts: &[usize], threshold: usize) -> Vec<std::ops::Range<usize>> {
            let flat: Vec<usize> = counts
                .iter()
                .enumerate()
                .flat_map(|(source, &count)| std::iter::repeat_n(source, count))
                .collect();
            if flat.is_empty() {
                return (!counts.is_empty())
                    .then_some(0..counts.len())
                    .into_iter()
                    .collect();
            }
            let mut ranges = Vec::new();
            let mut p = 0usize;
            let mut source_begin = 0usize;
            while p < flat.len() {
                let mut q = (p + threshold).min(flat.len() - 1);
                let boundary_source = flat[q];
                while q < flat.len() && flat[q] == boundary_source {
                    q += 1;
                }
                let source_end = if q == flat.len() {
                    counts.len()
                } else {
                    boundary_source + 1
                };
                ranges.push(source_begin..source_end);
                source_begin = source_end;
                p = q;
            }
            ranges
        }

        for a in 0..6 {
            for b in 0..6 {
                for c in 0..6 {
                    for d in 0..6 {
                        let counts = [a, b, c, d];
                        let hits = context_hits(&counts, 6);
                        for threshold in 1..5 {
                            assert_eq!(
                                translated_hit_partitions(&hits, counts.len(), 6, threshold),
                                reference(&counts, threshold),
                                "counts={counts:?}, threshold={threshold}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn long_query_output_batches_are_bounded_by_workers() {
        let record = |len| fasta::FastaRecord {
            quality: None,
            id: "query".to_owned(),
            sequence: vec![0; len],
        };
        let ordinary = vec![record(LONG_QUERY_RESIDUES); 5];
        let translated_contigs = vec![record(LONG_QUERY_RESIDUES + 1); 6];

        assert_eq!(output_query_batch_size(&ordinary, 4), (256, false));
        assert_eq!(output_query_batch_size(&translated_contigs, 4), (4, true));
        assert_eq!(output_query_batch_size(&translated_contigs, 1), (1, true));
        assert_eq!(output_query_batch_size(&translated_contigs, 0), (1, true));
    }

    #[test]
    fn automatic_hit_store_ceiling_is_lazy() {
        let mut store = AdaptiveHitStore::new(32, 16, 100, 1, Some(16_000_000_000), PathBuf::new());
        store.configure_storage(false).unwrap();
        assert_eq!(store.memory_limit, Some(16_000_000_000));
        assert!(store.disk.is_none());
        assert!(store
            .memory
            .as_ref()
            .unwrap()
            .bins
            .iter()
            .all(|bin| bin.capacity() == 0));
        assert_eq!(PackedHitBins::source_memory_reserve_hits(), (1 << 30) / 15);
    }

    #[test]
    fn zero_hit_store_ceiling_overrides_source_memory_flag() {
        let tmpdir = std::env::temp_dir().join(format!(
            "diamond-hit-store-zero-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ));
        std::fs::create_dir(&tmpdir).unwrap();
        let mut store = AdaptiveHitStore::new(2, 2, 100, 1, Some(0), tmpdir.clone());
        store.configure_storage(true).unwrap();
        assert!(store.memory.is_none());
        assert!(store.disk.is_some());
        drop(store);
        std::fs::remove_dir(tmpdir).unwrap();
    }

    #[test]
    fn adaptive_hit_store_migrates_existing_hits_and_keeps_writing_to_disk() {
        let tmpdir = std::env::temp_dir().join(format!(
            "diamond-adaptive-spill-test-{}",
            std::process::id()
        ));
        // Source requires an existing temporary directory; migration itself
        // must not create the directory supplied by the caller.
        std::fs::create_dir_all(&tmpdir).unwrap();
        let mut store = AdaptiveHitStore::new(2, 16, 100, 1, None, tmpdir.clone());
        store
            .ingest_partitions(vec![vec![StoredHit {
                query_id: 0,
                subject: 10,
                seed_offset: 2,
                score: 30,
            }]])
            .unwrap();
        assert_eq!(store.memory.as_ref().unwrap().total_hits(), 1);
        assert!(store.disk.is_none());

        // Cross the ceiling only after a hit already resides in memory. The
        // first hit must migrate and the incoming hit must use the persistent
        // disk writer rather than recreating an in-memory buffer.
        store.memory_limit = Some(0);
        store
            .ingest_partitions(vec![vec![StoredHit {
                query_id: 1,
                subject: 20,
                seed_offset: 3,
                score: 40,
            }]])
            .unwrap();
        assert!(store.memory.is_none());
        assert_eq!(store.disk.as_ref().unwrap().total_hits(), 2);

        store
            .ingest_partitions(vec![vec![StoredHit {
                query_id: 0,
                subject: 30,
                seed_offset: 4,
                score: 50,
            }]])
            .unwrap();
        assert!(store.memory.is_none());
        assert_eq!(store.disk.as_ref().unwrap().total_hits(), 3);
        drop(store);
        let _ = std::fs::remove_dir(tmpdir);
    }

    #[test]
    fn packed_memory_bins_preserve_query_hit_order() {
        assert_eq!(std::mem::size_of::<PackedRetainedHit>(), 15);
        let mut store = AdaptiveHitStore::new(4, 2, 100, 1, None, PathBuf::new());
        store
            .ingest_partitions(vec![vec![
                StoredHit {
                    query_id: 3,
                    subject: 30,
                    seed_offset: 3,
                    score: 13,
                },
                StoredHit {
                    query_id: 0,
                    subject: 10,
                    seed_offset: 1,
                    score: 11,
                },
                StoredHit {
                    query_id: 3,
                    subject: 31,
                    seed_offset: 4,
                    score: 14,
                },
                StoredHit {
                    query_id: 1,
                    subject: 20,
                    seed_offset: 2,
                    score: 12,
                },
            ]])
            .unwrap();

        let mut memory = store.memory.take().unwrap();
        let (first, begin, end) = memory.take_grouped_bin(0);
        assert_eq!((begin, end), (0, 2));
        assert_eq!(first[0][0].subject, 10);
        assert_eq!(first[1][0].subject, 20);
        let (second, begin, end) = memory.take_grouped_bin(1);
        assert_eq!((begin, end), (2, 4));
        assert!(second[0].is_empty());
        assert_eq!(
            second[1].iter().map(|hit| hit.subject).collect::<Vec<_>>(),
            vec![30, 31]
        );
    }

    #[test]
    fn stage2_window_keeps_anchor_relative_right_edge_when_left_clipped() {
        assert_eq!(stage2_query_bounds(232, 3, 48), (0, 51));
        assert_eq!(stage2_query_bounds(232, 100, 48), (52, 148));
        assert_eq!(stage2_query_bounds(120, 100, 48), (52, 120));
    }

    #[test]
    fn masking_disabled_clears_makedb_soft_masks_before_alignment() {
        const QUERY: &str = ">BAF52360.1\nMINSINSFFSSIPRSISSVTRNSSFTASQHKSTPNTVKTSSPLSPSNSPASATTIFKVKNSYTESGLQRSTSYTQSSIEKNALHRPLPDVAQRLVQHLAEHGIQPARNMAEHIPPAPNWPAPTPPVQNEQSRPLPDVAQRLVQHLAEHGIQPARNMAEHIPPAPNWPAPTPPVQNEQSRPLPDVAQRLVQHLAEHGIQPARNMAEHIPPAPNWPAPTPPVQNEQSRPLPDVAQRLMQHLAEHGIQPARNMAEHIPPAPNWPAPTPPVQNEQSRPLPDVAQRLMQHLAEHGINTSKRS*\n";
        const TARGET: &str = ">BAF52367.1\nMINSINSFFSSIPRSISSVMRNSSFTASQHKSTPNTVKTSSPLSPSNSPASATTIFKVKNSYTESGLQRSTSYTQSSIEKNALHRPLPDVAKRLVQHLAEHGIQPARNMAEHIPPAPNWPAPPPPVQNEQSRPLPDVAQRLMQHLAEHGIQPARNMAEHIPPAPNWPAPPPPVQNEQSRPLPDVAQRLMQHLAEHGIQPARNMAEHIPPAPNWPAPPPPVQNEQSRPLPDVAQRLVQHLAEHGIQPARNMAEHIPPAPNWPAPPPPVQNEQSRPLPDVAQRLMQHLAEHGIQPARNMAEHIPPAPNWPAPTPPVQNEQSRPLPDVAQRLMQHLAEHGIQPARNMAEHIPPAPNWPAPTPPVQNEQSRPLPDVAQRLMQHLAEHGINTSKRS*\n";

        let tmpdir = std::env::temp_dir().join(format!(
            "diamond-mask0-soft-mask-regression-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&tmpdir);
        std::fs::create_dir(&tmpdir).unwrap();
        let query_path = tmpdir.join("query.faa");
        let target_path = tmpdir.join("target.faa");
        let db_path = tmpdir.join("target.dmnd");
        let output_path = tmpdir.join("output.tsv");
        std::fs::write(&query_path, QUERY).unwrap();
        std::fs::write(&target_path, TARGET).unwrap();
        crate::data::db_builder::build_db(
            &[target_path.to_str().unwrap()],
            db_path.to_str().unwrap(),
            crate::basic::value::SequenceType::AminoAcid,
        )
        .unwrap();

        run(&BlastpConfig {
            query_files: vec![query_path.clone()],
            database: db_path.clone(),
            output: Some(output_path.clone()),
            matrix: b"blosum62".to_vec(),
            gap_open: 11,
            gap_extend: 1,
            max_evalue: 0.001,
            min_bit_score: 0.0,
            max_target_seqs: 25,
            max_hsps: 1,
            report_unaligned: -1,
            include_lineage: false,
            output_header: None,
            compression: Vec::new(),
            file_buffer_size: 67_108_864,
            sam: Default::default(),
            archive: Default::default(),
            tuning: Default::default(),
            scoring_matrix: None,
            auxiliary_output: Default::default(),
            taxonlist: Vec::new(),
            taxon_exclude: Vec::new(),
            seqidlist: std::path::PathBuf::new(),
            skip_missing_seqids: false,
            range_culling: false,
            range_cover: 50.0,
            culling_overlap: 50.0,
            taxon_k: 0,
            query_or_subject_cover: 0.0,
            reference_taxon_rank_ids: None,
            ext_chunk_size: 0,
            toppercent: None,
            global_ranking_targets: 0,
            frame_shift: 0,
            extension_mode: None,
            iterate: None,
            round_linear_target: false,
            round_query_skip: None,
            min_id: 0.0,
            approx_min_id: 0.0,
            soft_masking: crate::masking::MaskingAlgo::None,
            threads: 1,
            index_chunks: 0,
            block_size: None,
            query_indexed: false,
            contiguous_seed_mode: false,
            swipe_all: false,
            outfmt: vec![
                b"6".to_vec(),
                b"qseqid".to_vec(),
                b"sseqid".to_vec(),
                b"pident".to_vec(),
                b"length".to_vec(),
                b"mismatch".to_vec(),
                b"gapopen".to_vec(),
                b"qstart".to_vec(),
                b"qend".to_vec(),
                b"sstart".to_vec(),
                b"send".to_vec(),
                b"evalue".to_vec(),
                b"bitscore".to_vec(),
            ],
            sensitivity: Sensitivity::Default,
            masking: MaskingMode::None,
            motif_masking: Vec::new(),
            min_query_len: 0,
            query_cover: 0.0,
            subject_cover: 0.0,
            comp_based_stats: CbsMode::Disabled,
            no_self_hits: false,
            ungapped_xdrop_bits: 12.3,
            memory_limit: None,
            tmpdir: tmpdir.clone(),
            target_indexed: false,
            translated_query_layout: None,
        })
        .unwrap();

        // Exact upstream 2.1.24 result for this real TccP2 pair. Before the
        // fix, stored tantan bits truncated the alignment to 131 residues.
        assert_eq!(
            std::fs::read_to_string(&output_path).unwrap(),
            "BAF52360.1\tBAF52367.1\t95.6\t295\t13\t0\t1\t295\t1\t295\t3.24e-218\t587\n"
        );
        let _ = std::fs::remove_dir_all(tmpdir);
    }
}
