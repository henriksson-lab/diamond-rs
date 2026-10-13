use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::align::hsp::Match;
use crate::align::target::{extend as extend_targets, GappedScoreConfig};
use crate::align::ungapped::UngappedStageConfig;
use crate::basic::packed_loc::PackedLoc;
use crate::basic::reduction::Reduction;
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
    finish_daa_from_sequence_file, finish_daa_query_record, init_daa, write_daa_query_record,
    write_daa_record_hsp, DaaRunMetadata, DaaSequenceFile,
};
use crate::output::format::{self, FieldId, Hsp as OutputHsp};
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
fn stage2_query_bounds(query_len: usize, seed_pos: usize) -> (usize, usize) {
    let window = crate::dp::ungapped_window::UNGAPPED_WINDOW;
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
struct StoredHit {
    subject: u64,
    query_id: u32,
    seed_offset: u32,
    score: u16,
}

struct TargetIndexedState {
    query: Block,
    shapes: ShapeConfig,
    seeds: crate::data::seed_set::HashedSeedSet,
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
    mut map_matches: F,
) -> (Vec<StoredHit>, usize)
where
    F: FnMut((u32, u32), &[u32], usize, &[PackedLoc], u32, &mut Vec<StoredHit>),
{
    const PACKED_LOC_SIZE: usize = 5;
    let mut query = query_hits.begin_with_elem_size(PACKED_LOC_SIZE);
    let mut target = target_hits.begin_with_elem_size(PACKED_LOC_SIZE);
    let mut output = Vec::new();
    let mut raw_count = 0usize;

    while query.good() && target.good() {
        let query_bytes = query.range_bytes(PACKED_LOC_SIZE);
        let target_bytes = target.range_bytes(PACKED_LOC_SIZE);
        raw_count = raw_count
            .saturating_add((query.count() as usize).saturating_mul(target.count() as usize));
        let query_locs = packed_locs_from_bytes(query_bytes);
        let target_locs = packed_locs_from_bytes(target_bytes);
        crate::search::hamming::stage1(
            query_locs,
            target_locs,
            query_seqs,
            target_seqs,
            1024,
            hamming_filter_id,
            |hits, query_begin, target_begin, queries, targets| {
                for local_query in 0..hits.query_count() {
                    let query_position = queries[query_begin + local_query].as_u64() as usize;
                    let (query_id, query_pos) = query_seqs.local_position(query_position);
                    let target_indices = hits.hits(local_query);
                    if !target_indices.is_empty() {
                        map_matches(
                            (query_id as u32, query_pos as u32),
                            target_indices,
                            target_begin,
                            targets,
                            absolute_partition,
                            &mut output,
                        );
                    }
                }
            },
        );
        query.increment(PACKED_LOC_SIZE);
        target.increment(PACKED_LOC_SIZE);
    }
    debug_assert_eq!(query.good(), target.good());
    (output, raw_count)
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
    &[i32],
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
    &[i32],
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
    fn new(key_partition: Vec<u32>, memory_limit: Option<u64>) -> Self {
        const UPSTREAM_BIN_RESERVE_BYTES: u64 = 1 << 30;
        let bin_count = key_partition.len().max(1) as u64;
        let reserve_bytes = memory_limit
            .map(|limit| (limit / bin_count).min(UPSTREAM_BIN_RESERVE_BYTES))
            .unwrap_or(0);
        let reserve_hits = (reserve_bytes / std::mem::size_of::<PackedRetainedHit>() as u64)
            .min(usize::MAX as u64) as usize;
        Self {
            bins: (0..key_partition.len())
                .map(|_| Vec::with_capacity(reserve_hits))
                .collect(),
            key_partition,
        }
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
    memory: Option<PackedHitBins>,
    disk: Option<HitBuffer>,
    memory_limit: Option<u64>,
    tmpdir: PathBuf,
    query_count: usize,
    query_bins: usize,
    max_subject: u64,
    query_group_size: usize,
    bytes_since_rss_check: usize,
    rss_checked_once: bool,
}

impl AdaptiveHitStore {
    fn new(
        query_count: usize,
        query_bins: usize,
        max_subject: u64,
        query_group_size: usize,
        memory_limit: Option<u64>,
        tmpdir: PathBuf,
    ) -> Self {
        let key_partition = Self::key_partition(query_count, query_bins, query_group_size);
        Self {
            memory: Some(PackedHitBins::new(key_partition, memory_limit)),
            disk: None,
            memory_limit,
            tmpdir,
            query_count,
            query_bins,
            max_subject,
            query_group_size: query_group_size.max(1),
            bytes_since_rss_check: 0,
            rss_checked_once: false,
        }
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
        let key_partition =
            Self::key_partition(self.query_count, self.query_bins, self.query_group_size);
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
        eprintln!(
            "Hit buffer: RSS reached --memory-limit; spilling {migrated_hits} retained hits to {}",
            if self.tmpdir.as_os_str().is_empty() {
                std::env::temp_dir().display().to_string()
            } else {
                self.tmpdir.display().to_string()
            }
        );
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
    cutoffs: &[i32],
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
    cutoffs: &[i32],
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
    cutoffs: &[i32],
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
    cutoffs: &[i32],
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
    cutoffs: &[i32],
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
    cutoffs: &[i32],
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
    cutoffs: &[i32],
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
    let current_range = if USE_LEFT_MOST_RANGE {
        let partitions = seedp_count(10) as usize;
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
    let (q_start, q_end) = stage2_query_bounds(query.len(), q_pos);
    let window_left = q_pos - q_start;
    let window_clipped = q_end - q_start;
    let query_window = &query[q_start..q_end];
    let cutoff = cutoffs[query_id];
    let interval_mod = (q_pos % 32) as i32;
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
            ungapped_count += 1;
            let subject = (subject_starts[hit_index] + window_left as isize) as usize;
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
    ungapped_count
}

/// Port of `align.cpp::make_partition` for translated queries. Hit records
/// are already grouped by context, so count all six contexts of a source and
/// extend every threshold crossing through that complete source query.
fn translated_hit_partitions(
    hits_by_context: &[Vec<CompactHit>],
    source_count: usize,
    query_contexts: usize,
    min_task_trace_points: usize,
) -> Vec<std::ops::Range<usize>> {
    assert!(query_contexts > 0 && min_task_trace_points > 0);
    debug_assert!(hits_by_context.len() >= source_count.saturating_mul(query_contexts));
    if source_count == 0 {
        return Vec::new();
    }
    let mut ranges = Vec::new();
    let mut task_begin = 0usize;
    let mut trace_points = 0usize;
    for source in 0..source_count {
        let context_begin = source * query_contexts;
        trace_points = trace_points.saturating_add(
            hits_by_context[context_begin..context_begin + query_contexts]
                .iter()
                .map(Vec::len)
                .sum::<usize>(),
        );
        // Upstream probes `p + min_task_trace_pts` (the 1025th record for
        // the default 1024), then extends through that record's whole query.
        if trace_points > min_task_trace_points {
            ranges.push(task_begin..source + 1);
            task_begin = source + 1;
            trace_points = 0;
        }
    }
    if task_begin < source_count {
        if trace_points == 0 {
            if let Some(last) = ranges.last_mut() {
                // HitIterator appends trailing no-hit queries while fetching
                // the final hit partition; they are not a separate task.
                last.end = source_count;
            } else {
                ranges.push(0..source_count);
            }
        } else {
            ranges.push(task_begin..source_count);
        }
    }
    ranges
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

fn output_query_batch_size(records: &[fasta::FastaRecord], threads: usize) -> (usize, bool) {
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

/// Configuration for a blastp run.
#[derive(Clone)]
pub struct BlastpConfig {
    pub query_files: Vec<String>,
    pub database: String,
    pub output: Option<String>,
    pub matrix: String,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub max_evalue: f64,
    pub max_target_seqs: i64,
    pub ext_chunk_size: i64,
    pub toppercent: Option<f64>,
    pub global_ranking_targets: i64,
    pub min_id: f64,
    pub threads: i32,
    pub outfmt: Vec<String>,
    pub sensitivity: Sensitivity,
    pub masking: MaskingMode,
    pub motif_masking: String,
    pub min_query_len: usize,
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
pub struct TranslatedQuerySource {
    pub id: String,
    pub dna_len: i32,
    /// Original nucleotide query, retained for native DAA records. The six
    /// translated protein contexts are owned separately by the extension
    /// groups, so this is the only stored source-DNA copy.
    pub sequence: Vec<Letter>,
    /// Adjacent source-query ordinal for synchronized two-file blastx input.
    /// Storing an index instead of a second sequence keeps the mate relation
    /// zero-copy, matching upstream's shared source-sequence block.
    pub mate_index: Option<usize>,
}

#[derive(Clone, Debug)]
pub struct TranslatedQueryLayout {
    pub sources: Vec<TranslatedQuerySource>,
}

struct TranslatedExtensionGroup<'a> {
    id: String,
    dna_len: i32,
    source_sequence: &'a [Letter],
    mate_sequence: Option<&'a [Letter]>,
    frames: Vec<Vec<Letter>>,
}

enum QueryOutput {
    Tabular { bytes: Vec<u8>, rows: u64 },
    Daa(Vec<Match>),
}

struct NativeDaaDictionary {
    sequence_count: u64,
    refs: Vec<(String, u32)>,
}

struct NativeDatabaseMetadata {
    subject_taxids: Vec<Vec<TaxId>>,
    taxonomy: Option<TaxonomyTree>,
}

fn is_blast_database(path: &str) -> bool {
    let path = Path::new(path);
    let mut pin = path.as_os_str().to_owned();
    pin.push(".pin");
    let mut pal = path.as_os_str().to_owned();
    pal.push(".pal");
    path.extension().is_some_and(|extension| extension == "pal")
        || PathBuf::from(pin).exists()
        || PathBuf::from(pal).exists()
}

fn requested_taxonomy_flags(outfmt: &[String]) -> SequenceFileFlags {
    let mut flags = SequenceFileFlags::NONE;
    for field in outfmt
        .iter()
        .skip(1)
        .filter_map(|name| FieldId::from_name(name))
    {
        match field {
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
    path: &str,
    taxonomy_flags: SequenceFileFlags,
) -> io::Result<(Vec<fasta::FastaRecord>, NativeDatabaseMetadata)> {
    let flags = SequenceFileFlags::SEQS | SequenceFileFlags::TITLES | taxonomy_flags;
    let mut database = BlastDB::new(path, flags).map_err(io::Error::other)?;
    let sequence_count = database.sequence_count() as usize;
    let mut records = Vec::with_capacity(sequence_count);
    let mut subject_taxids = Vec::with_capacity(sequence_count);
    for oid in 0..sequence_count {
        let mut sequence = Vec::new();
        database
            .seq_data(oid, &mut sequence)
            .map_err(io::Error::other)?;
        let id = database
            .seqid(oid as u64, true, true)
            .map_err(io::Error::other)?;
        records.push(fasta::FastaRecord { id, sequence });
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
                        database.taxon_scientific_name(taxid)
                    } else {
                        taxid.to_string()
                    };
                    tree.add_node(TaxonomyNode {
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
    Vec<fasta::FastaRecord>,
    NativeDatabaseMetadata,
)> {
    let needs_mapping = taxonomy_flags.contains(SequenceFileFlags::TAXON_MAPPING);
    let needs_nodes = taxonomy_flags.contains(SequenceFileFlags::TAXON_NODES)
        || taxonomy_flags.contains(SequenceFileFlags::TAXON_RANKS);
    let needs_names = taxonomy_flags.contains(SequenceFileFlags::TAXON_SCIENTIFIC_NAMES);
    let mut database = DatabaseFile::open_with_options(
        path,
        DatabaseOpenOptions {
            load_taxon_mapping: needs_mapping,
            load_taxon_nodes: needs_nodes,
            load_taxon_scientific_names: needs_names,
            ..DatabaseOpenOptions::default()
        },
    )
    .map_err(io::Error::other)?;

    // Ranks were added to the native node payload at build 131.  C++ rejects
    // an older database as soon as an output field requests rank lookup.
    if taxonomy_flags.contains(SequenceFileFlags::TAXON_RANKS) && database.ref_header.build < 131 {
        return Err(io::Error::other(
            "Options require taxonomy information included in the database. Please use the respective options to build this information into the database when running diamond makedb: taxonomy ranks information (database needs to be built with diamond version >= 0.9.30",
        ));
    }

    let header = database.ref_header.clone();
    let sequence_count = database.sequence_count() as usize;
    let mut records = Vec::with_capacity(sequence_count);
    let mut subject_taxids = Vec::with_capacity(sequence_count);
    database.init_seq_access();
    for oid in 0..sequence_count {
        let (sequence, id) = database
            .read_seq()
            .map_err(io::Error::other)?
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Unexpected end of DIAMOND database",
                )
            })?;
        records.push(fasta::FastaRecord { id, sequence });
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
                        database.taxon_scientific_name(taxid)
                    } else {
                        taxid.to_string()
                    };
                    tree.add_node(TaxonomyNode {
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

impl DaaSequenceFile for NativeDaaDictionary {
    fn sequence_count(&self) -> u64 {
        self.sequence_count
    }

    fn dict_size(&self) -> usize {
        self.refs.len()
    }

    fn dict_title(&self, index: usize) -> String {
        self.refs[index].0.clone()
    }

    fn dict_len(&self, index: usize) -> u32 {
        self.refs[index].1
    }
}

fn encode_daa_query(
    output: &mut Vec<u8>,
    query_name: &str,
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

    let seek_pos = write_daa_query_record(output, query_name, query_source, input_sequence_type);
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
    run_impl(
        config,
        None,
        None,
        None,
        None,
        MaskingAlgo::None,
        false,
        false,
        false,
        0.0,
        false,
        None,
        None,
        None,
    )
}

/// Run the native search with already translated query records.
///
/// `blastx` owns these records and transfers them directly, matching the
/// upstream in-process query-block handoff without serializing and reparsing
/// an intermediate FASTA file.
pub fn run_with_query_records(
    config: &BlastpConfig,
    query_records: Vec<fasta::FastaRecord>,
) -> io::Result<()> {
    run_impl(
        config,
        None,
        Some(query_records),
        None,
        None,
        MaskingAlgo::None,
        false,
        false,
        false,
        0.0,
        false,
        None,
        None,
        None,
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
fn restore_extension_queries(records: &mut [fasta::FastaRecord], restored: &SequenceSet) {
    assert_eq!(records.len(), restored.len());
    for (query_id, record) in records.iter_mut().enumerate() {
        record.sequence.clear();
        record.sequence.extend_from_slice(restored.get(query_id));
    }
}

/// Run the native protein search pipeline on already decoded records and
/// return clustering edge records instead of formatting a file.
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
    extension_mode_override: Option<sensitivity::ExtensionMode>,
) -> io::Result<Vec<InMemorySearchEdge>> {
    let edges = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    run_impl(
        config,
        Some(database),
        queries,
        Some(edges.clone()),
        Some(approx_min_id),
        soft_masking_override,
        linear_stage1_query,
        linear_stage1_target,
        self_search,
        query_or_target_cover,
        mutual_cover,
        Some(hamming_ext),
        soft_tantan_masker,
        extension_mode_override,
    )?;
    let result = edges.lock().unwrap().clone();
    Ok(result)
}

fn run_impl(
    config: &BlastpConfig,
    database_records: Option<Vec<fasta::FastaRecord>>,
    query_records_override: Option<Vec<fasta::FastaRecord>>,
    edge_output: Option<std::sync::Arc<std::sync::Mutex<Vec<InMemorySearchEdge>>>>,
    approx_min_id_override: Option<f64>,
    soft_masking_override: MaskingAlgo,
    linear_stage1_query: bool,
    linear_stage1_target: bool,
    self_search: bool,
    query_or_target_cover: f64,
    mutual_cover: bool,
    hamming_ext: Option<HammingExtConfig>,
    soft_tantan_masker: Option<&crate::masking::tantan::TantanMasker>,
    extension_mode_override: Option<sensitivity::ExtensionMode>,
) -> io::Result<()> {
    let start = Instant::now();
    // Upstream linear self-search aliases target to query instead of loading
    // or copying a second block. `run_edges_in_memory` signals that exact mode
    // by supplying the database records without a separate query vector.
    let shared_self_records = database_records.is_some()
        && query_records_override.is_none()
        && edge_output.is_some()
        && self_search;
    // `Search::Config` derives this once from equal query/subject coverage
    // and stage sensitivity. It selects the mutual-coverage stage-1 kernels,
    // not merely length sorting, so carry it through the native seed worker.
    let min_length_ratio =
        if config.query_cover >= 50.0 && config.query_cover == config.subject_cover {
            if linear_stage1_query && config.sensitivity < Sensitivity::Linclust40 {
                (config.query_cover / 100.0 + 0.05).min(0.92)
            } else {
                (config.query_cover / 100.0 - 0.05).max(0.0)
            }
        } else {
            0.0
        };

    // Honor `--threads N`. Previously the value was parsed but ignored, so
    // rayon fell back to its global pool (= `RAYON_NUM_THREADS` env or core
    // count). `try_build_global` is a no-op once the global pool exists, so
    // first call wins — subsequent invocations (e.g. multiple blastp runs in
    // the same process) silently keep the original count. For a one-shot CLI
    // run this is correct; for library use callers should install their own
    // pool. `config.threads <= 0` keeps the default.
    if config.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads as usize)
            .build_global();
    }

    // Load database sequences - try DMND first, then FASTA.
    // C++ `auto_append_extension_if_exists` (`config.cpp:770`) APPENDS `.dmnd`
    // only when the file as-given doesn't exist — it never STRIPS an existing
    // extension. So `--db nr.fasta` resolves to `nr.fasta` in C++ (then parsed
    // as FASTA), but `--db nr.fasta` was resolving to `nr.dmnd` in Rust because
    // `with_extension("dmnd")` REPLACES rather than appends. Fix: only consult
    // a `<db>.dmnd` companion when the given path doesn't exist.
    let taxonomy_flags = requested_taxonomy_flags(&config.outfmt);
    let taxonomy_requested = taxonomy_flags.contains(SequenceFileFlags::TAXON_MAPPING);
    let mut database_metadata = None;
    let (mut db_records, _db_from_dmnd) = if let Some(records) = database_records {
        if taxonomy_requested {
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
        } else if db_path.exists() {
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
            eprintln!(
                "Database: {} sequences (NCBI BLAST database)",
                records.len()
            );
            database_metadata = Some(metadata);
            (records, false)
        } else if dmnd_path.exists() {
            // Keep ordinary search on the lean sequence-only reader.  An
            // empty `Vec<TaxId>` per subject would cost 24 bytes per database
            // sequence even when no taxonomy column was requested.
            let (header, records) = if taxonomy_requested {
                let (header, records, metadata) = load_dmnd_database(&dmnd_path, taxonomy_flags)?;
                database_metadata = Some(metadata);
                (header, records)
            } else {
                crate::data::dmnd_reader::read_dmnd(&dmnd_path)?
            };
            eprintln!(
                "Database: {} sequences, {} letters (DMND)",
                header.sequences, header.letters
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
            let records = fasta::read_fasta_file(&fasta_path, SequenceType::AminoAcid)?;
            eprintln!("Database: {} sequences (FASTA)", records.len());
            (records, false)
        }
    };

    // Compute total database letters for E-value normalization
    let db_letters: u64 = db_records.iter().map(|r| r.sequence.len() as u64).sum();

    // Load scoring matrix with database size
    let score_matrix = ScoreMatrix::new(
        &config.matrix,
        config.gap_open,
        config.gap_extend,
        0,
        1,
        db_letters,
    )
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    // Load query sequences
    let mut query_records = if shared_self_records {
        Vec::new()
    } else if let Some(records) = query_records_override {
        records
    } else {
        let mut records = Vec::new();
        for qf in &config.query_files {
            records.extend(fasta::read_fasta_file(
                Path::new(qf),
                SequenceType::AminoAcid,
            )?);
        }
        records
    };
    eprintln!(
        "Queries: {} sequences",
        if shared_self_records {
            db_records.len()
        } else {
            query_records.len()
        }
    );

    // Upstream derives a minimum length ratio when the two coverage cutoffs
    // are equal and at least 50%, then length-sorts both sides before search
    // (`run/config.cpp` and `run/double_indexed.cpp`).  Besides improving the
    // stage-1 schedule, this determines the observable query order in tabular
    // output.  Retain the original ordinal as the tie breaker: C++ sorts
    // `pair<length, BlockId>` with `greater`, so equal-length records appear
    // in descending input order as well.
    if !shared_self_records
        && config.query_cover >= 50.0
        && config.query_cover == config.subject_cover
    {
        let sort_by_cpp_length_order = |records: &mut Vec<fasta::FastaRecord>| {
            let mut indexed: Vec<_> = records.drain(..).enumerate().collect();
            indexed.sort_unstable_by(|(left_id, left), (right_id, right)| {
                right
                    .sequence
                    .len()
                    .cmp(&left.sequence.len())
                    .then_with(|| right_id.cmp(left_id))
            });
            records.extend(indexed.into_iter().map(|(_, record)| record));
        };
        // The clustering backend has already applied the upstream
        // `Block::length_sort` order before passing an aliased self-search
        // block. Sorting that block again here reverses equal-length groups
        // and invalidates its local-to-database OID mapping.
        // blastx contexts are deliberately consecutive: the spill bins and
        // extension stage consume all six as one source query.
        if !shared_self_records && config.translated_query_layout.is_none() {
            sort_by_cpp_length_order(&mut query_records);
        }
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
            for (_, (record, taxids)) in indexed {
                db_records.push(record);
                metadata.subject_taxids.push(taxids);
            }
        } else {
            sort_by_cpp_length_order(&mut db_records);
        }
    }

    use crate::basic::value::{MASK_LETTER, SEED_MASK};
    use rayon::iter::IntoParallelRefMutIterator;
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
            let tantan_masker =
                crate::masking::tantan::TantanMasker::from_score_matrix(&score_matrix, 0.9);
            db_records.par_iter_mut().for_each(|r| {
                crate::masking::remove_bit_mask(&mut r.sequence);
                tantan_masker.mask(&mut r.sequence);
            });
            let mask_query = |record: &mut fasta::FastaRecord| {
                crate::masking::remove_bit_mask(&mut record.sequence);
                tantan_masker.mask(&mut record.sequence);
            };
            if shared_self_records {
                // The database mutation above is the query mutation too.
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
                crate::masking::mask_sequence(&mut record.sequence, MaskingAlgo::Seg);
            });
            if !shared_self_records {
                query_records
                    .par_iter_mut()
                    .for_each(|record| crate::masking::remove_bit_mask(&mut record.sequence));
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
        db_records
            .par_iter_mut()
            .for_each(|r| hard_mask(&mut r.sequence));
        if !shared_self_records {
            query_records
                .par_iter_mut()
                .for_each(|r| hard_mask(&mut r.sequence));
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
    let mut stage2_query_set = SequenceSet::new();
    stage2_query_set.reserve_capacity(
        query_source.len(),
        query_source
            .iter()
            .map(|record| record.sequence.len())
            .sum(),
    );
    for record in query_source {
        stage2_query_set.push(&record.sequence);
    }
    let mut db_block = Block::new();
    for (idx, record) in db_records.iter().enumerate() {
        db_block
            .push_back(
                &record.sequence,
                Some(&record.id),
                None,
                idx as u64,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .map_err(io::Error::other)?;
    }

    // Clustering config sets hard query/target masking to 0 but adds TANTAN
    // to Search::Config::soft_masking. Apply it only to seed-enumeration
    // records: `db_block` and `stage2_query_set` above deliberately retain the
    // original residues for ungapped/gapped extension.
    let mut soft_query_spans = if soft_masking_override != MaskingAlgo::None {
        vec![Vec::<(usize, usize)>::new(); stage2_query_set.len()]
    } else {
        Vec::new()
    };
    if soft_masking_override == MaskingAlgo::Tantan {
        // Upstream constructs its score-matrix-derived `Masking` singleton
        // once before entering the cascaded rounds. Ordinary searches retain
        // the local construction, while clustering supplies that command-wide
        // state so PAM matrices do not repeat the expensive lambda solve in
        // every round.
        let owned_tantan_masker = soft_tantan_masker
            .is_none()
            .then(|| crate::masking::tantan::TantanMasker::from_score_matrix(&score_matrix, 0.9));
        let tantan_masker = soft_tantan_masker
            .or(owned_tantan_masker.as_ref())
            .expect("soft masking requires a TANTAN masker");
        if shared_self_records {
            db_records
                .par_iter_mut()
                .zip(soft_query_spans.par_iter_mut())
                .for_each(|(record, spans)| {
                    let ranges = tantan_masker.mask_ranges(&mut record.sequence);
                    let (front, back) = ranges.as_slices();
                    for &(begin, end) in front.iter().chain(back) {
                        let begin = begin as usize;
                        let end = end as usize;
                        spans.push((begin, end));
                        record.sequence[begin..end].fill(crate::basic::value::MASK_LETTER);
                    }
                });
        } else {
            db_records.par_iter_mut().for_each(|record| {
                let ranges = tantan_masker.mask_ranges(&mut record.sequence);
                let (front, back) = ranges.as_slices();
                for &(begin, end) in front.iter().chain(back) {
                    record.sequence[begin as usize..end as usize]
                        .fill(crate::basic::value::MASK_LETTER);
                }
            });
            query_records
                .par_iter_mut()
                .zip(soft_query_spans.par_iter_mut())
                .for_each(|(record, spans)| {
                    let ranges = tantan_masker.mask_ranges(&mut record.sequence);
                    let (front, back) = ranges.as_slices();
                    for &(begin, end) in front.iter().chain(back) {
                        let begin = begin as usize;
                        let end = end as usize;
                        spans.push((begin, end));
                        record.sequence[begin..end].fill(crate::basic::value::MASK_LETTER);
                    }
                });
        }
    } else if soft_masking_override == MaskingAlgo::Seg {
        let mask_record = |record: &mut fasta::FastaRecord, spans: &mut Vec<(usize, usize)>| {
            let encoded = record
                .sequence
                .iter()
                .map(|&letter| (letter & LETTER_MASK) as u8)
                .collect::<Vec<_>>();
            for (begin, end) in crate::masking::blast_seg::seq_buffer_seg(&encoded, 0, None) {
                let begin = begin as usize;
                let end = end as usize + 1;
                spans.push((begin, end));
                record.sequence[begin..end].fill(MASK_LETTER);
            }
        };
        if shared_self_records {
            db_records
                .par_iter_mut()
                .zip(soft_query_spans.par_iter_mut())
                .for_each(|(record, spans)| mask_record(record, spans));
        } else {
            db_records.par_iter_mut().for_each(|record| {
                let mut ignored = Vec::new();
                mask_record(record, &mut ignored);
            });
            query_records
                .par_iter_mut()
                .zip(soft_query_spans.par_iter_mut())
                .for_each(|(record, spans)| mask_record(record, spans));
        }
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
    let soft_masking = sensitivity::soft_masking_algo(
        &sensitivity::get_traits(config.sensitivity),
        &config.motif_masking,
        false,
        false,
    )
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let use_motif_masking = soft_masking == MaskingAlgo::Motif;
    let db_motif_saves: Vec<Vec<crate::masking::motifs::MotifMaskEntry>> = if use_motif_masking {
        db_records
            .par_iter_mut()
            .map(|r| crate::masking::motifs::mask_motifs(&mut r.sequence))
            .collect()
    } else {
        vec![Vec::new(); db_records.len()]
    };
    let query_motif_saves: Vec<Vec<crate::masking::motifs::MotifMaskEntry>> = if use_motif_masking {
        if shared_self_records {
            Vec::new()
        } else {
            query_records
                .par_iter_mut()
                .map(|r| crate::masking::motifs::mask_motifs(&mut r.sequence))
                .collect()
        }
    } else {
        vec![Vec::new(); query_records.len()]
    };
    // Set up seed extraction using sensitivity-appropriate shapes
    let reduction = Reduction::default_reduction();
    let shape_codes = sensitivity::get_shape_codes(config.sensitivity);
    let shapes: Vec<Shape> = shape_codes
        .iter()
        .map(|code| Shape::from_code(code, &reduction))
        .collect();
    eprintln!(
        "Sensitivity: {:?}, shapes: {} (weights: {})",
        config.sensitivity,
        shapes.len(),
        shapes
            .iter()
            .map(|s| s.weight.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );

    // Upstream target-indexed mode changes the seed encoding itself: both
    // sides use Murmur-hashed spaced seeds, while the persistent hash tables
    // filter the one-pass query enumeration. Reuse the source-mirrored stage0
    // implementation for that path instead of accepting the flag while still
    // running the ordinary spaced-factor arrays.
    let mut target_indexed_state = if config.target_indexed {
        let shape_codes = shape_codes
            .iter()
            .map(|code| (*code).to_owned())
            .collect::<Vec<_>>();
        let shape_config = ShapeConfig::from_codes(&shape_codes, 0, &reduction)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?;
        let mut query = Block::new();
        for (oid, record) in query_records.iter().enumerate() {
            query
                .push_back(
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
        let histogram_cfg = crate::data::flags::EnumCfg {
            partition: None,
            shape_begin: 0,
            shape_end: shape_config.count(),
            code: crate::data::flags::SeedEncoding::Hashed,
            skip: None,
            filter_masked_seeds: false,
            mask_seeds: false,
            seed_cut: 0.0,
            soft_masking,
            minimizer_window: 0,
            filter_low_complexity_seeds: false,
            mask_low_complexity_seeds: false,
            sketch_size: 0,
        };
        let histogram = crate::data::seed_histogram::SeedHistogram::from_block(
            &mut db_block,
            false,
            &crate::data::flags::NO_FILTER,
            &histogram_cfg,
            10,
            rayon::current_num_threads().max(1) as u32,
            &shape_config,
            &reduction,
            0,
            1,
        )
        .map_err(io::Error::other)?;
        *db_block.hst() = histogram;

        let database_path = {
            let given = Path::new(&config.database);
            if given.exists() {
                given.to_path_buf()
            } else {
                let mut path = given.as_os_str().to_owned();
                path.push(".dmnd");
                PathBuf::from(path)
            }
        };
        let index_path = {
            let mut path = database_path.as_os_str().to_owned();
            path.push(".seed_idx");
            PathBuf::from(path)
        };
        let seeds =
            crate::data::seed_set::HashedSeedSet::from_index_file(&index_path, &shape_config)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
        eprintln!("Loaded persistent seed index: {}", index_path.display());
        Some(TargetIndexedState {
            query,
            shapes: shape_config,
            seeds,
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
    let traits = sensitivity::get_traits(config.sensitivity);
    let hamming_filter_id = traits.min_identities.max(sensitivity::hamming_id_cutoff(
        approx_min_id_override.unwrap_or(0.0),
    ));
    // Process shapes in order and let each shape's seed-array build/join use
    // the Rayon pool internally. Running the outer shape loop in parallel
    // competes with the inner reference seed-array work and scales worse on
    // small real query blocks; collecting in shape order preserves the previous
    // output order.
    let ungapped_evalue = traits.ungapped_evalue as f64;
    let short_query_ungapped_cutoff = if ungapped_evalue > 0.0 {
        score_matrix.rawscore_int(25.0)
    } else {
        0
    };
    let ungapped_cutoffs: Vec<i32> = (0..stage2_query_set.len())
        .map(|query_id| {
            let query = stage2_query_set.get(query_id);
            crate::search::stage2::ungapped_cutoff(
                query.len() as i32,
                ungapped_evalue,
                60,
                short_query_ungapped_cutoff,
                false,
                |q| score_matrix.ungapped_cutoff(q as usize, ungapped_evalue),
                |q| score_matrix.ungapped_cutoff(q as usize, ungapped_evalue),
            )
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
    for (query_id, saved) in query_motif_saves_for_restore.iter().enumerate() {
        crate::masking::motifs::restore_motifs(
            stage2_query_set.get_mut(query_id),
            saved,
            template_len,
        );
    }
    // C++ MaskingTable::remove restores TANTAN-masked residues for extension,
    // then annotates every seed start whose template overlaps the restored
    // range. Preserve the residues already stored in stage2_query_set and apply
    // precisely that left-extended SEED_MASK interval.
    for (query_id, spans) in soft_query_spans.iter().enumerate() {
        let query = stage2_query_set.get_mut(query_id);
        for &(begin, end) in spans {
            let mask_begin = begin.saturating_sub((template_len.max(1) - 1) as usize);
            for letter in &mut query[mask_begin..end] {
                *letter |= SEED_MASK;
            }
        }
    }
    // Hamming fingerprints in upstream are loaded after motif masking has
    // been removed, directly from the original contiguous sequences. Use the
    // pre-motif stage-2 copies here as well: seed enumeration still consumes
    // the masked `query_seqs`/`db_seqs`, while fingerprint loading no longer
    // needs a sparse restoration lookup for every candidate window. Query
    // SEED_MASK annotations are harmless because the fingerprint loader
    // strips LETTER_MASK from every byte.
    let effective_index_chunks = if config.target_indexed {
        1
    } else {
        traits.index_chunks as usize
    };
    let chunked = effective_index_chunks > 1;
    let use_left_most_range = chunked;
    // Upstream stage2 skips this filter only for minimizer/sketch and linear
    // search modes; sensitivity itself is not a bypass. Keeping all hits for
    // very/ultra-sensitive searches inflated the spill file by two orders of
    // magnitude (and subsequently decoded those unnecessary hits into RAM).
    let skip_left_most =
        traits.minimizer_window > 0 || traits.sketch_size > 0 || linear_stage1_query;
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
    let mut hit_store = AdaptiveHitStore::new(
        stage2_query_set.len(),
        traits.query_bins as usize,
        db_block.seqs().data().len() as u64,
        if config.translated_query_layout.is_some() {
            6
        } else {
            1
        },
        config.memory_limit,
        config.tmpdir.clone(),
    );
    // A zero ceiling explicitly requests upstream-compatible disk buffering;
    // enter that mode before stage 1 rather than waiting for the first batch.
    if config.memory_limit == Some(0) {
        hit_store.spill_to_disk()?;
    }
    let mut retained_hit_count = 0usize;
    let mut raw_seed_matches = 0usize;
    let mut current_left_most_matcher = PatternMatcher::new(&[]);
    // Upstream constructs one worker×shape×partition SeedHistogram for the
    // complete query/reference block, then reuses it for every index chunk.
    // Keep worker sequence ranges in the histogram so parallel fills preserve
    // the same stable record order without atomics.
    let precomputed_worker_histograms = (target_indexed_state.is_none()
        && rayon::current_num_threads() > 1
        && traits.sketch_size <= 0)
        .then(|| {
            (
                crate::search::seed_array::SeedArray::count_partitions_by_worker_all_shapes(
                    &query_seqs,
                    &shapes,
                    &reduction,
                    10,
                    config.min_query_len,
                ),
                crate::search::seed_array::SeedArray::count_partitions_by_worker_all_shapes(
                    &db_seqs,
                    &shapes,
                    &reduction,
                    10,
                    config.min_query_len,
                ),
            )
        });
    for (shape_id, shape) in shapes.iter().enumerate() {
        let previous_left_most_matcher = current_left_most_matcher.clone();
        current_left_most_matcher.add_pattern(shape.mask);
        let left_most_context = LeftMostContext {
            previous_matcher: previous_left_most_matcher,
            current_matcher: current_left_most_matcher.clone(),
            short_query_ungapped_cutoff: score_matrix.rawscore_int(25.0),
            seedp_mask: seedp_mask(10),
            reduction: &reduction,
        };
        let complexity_cut = traits.seed_cut * std::f64::consts::LN_2 * shape.weight as f64;

        if let Some(state) = target_indexed_state.as_mut() {
            let mut stage0_config =
                crate::search::stage0::SearchShapeConfig::new(&state.shapes, &reduction);
            stage0_config.seedp_bits = 10;
            stage0_config.index_chunks = 1;
            stage0_config.threads = rayon::current_num_threads().max(1);
            stage0_config.seed_encoding = crate::data::flags::SeedEncoding::Hashed;
            stage0_config.seed_complexity_cut = complexity_cut;
            stage0_config.soft_masking = soft_masking;
            stage0_config.minimizer_window = traits.minimizer_window;
            stage0_config.sketch_size = traits.sketch_size;
            stage0_config.min_query_len = config.min_query_len as i32;
            stage0_config.query_contexts = if config.translated_query_layout.is_some() {
                6
            } else {
                1
            };
            stage0_config.target_seeds = Some(&state.seeds);

            let indexed_report = crate::search::stage0::search_shape(
                shape_id as u32,
                &mut state.query,
                &mut db_block,
                &stage0_config,
                |partition, query_stage0, target_stage0| {
                    let (hits, raw) = map_target_indexed_partition(
                        partition.query_seed_hits,
                        partition.reference_seed_hits,
                        partition.absolute_partition,
                        query_stage0,
                        target_stage0,
                        hamming_filter_id,
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
                                    false,
                                    false,
                                    skip_left_most,
                                    1,
                                    hamming_filter_id,
                                    ungapped_lane_count,
                                    ungapped_kernel,
                                );
                            }
                        },
                    );
                    raw_seed_matches = raw_seed_matches.saturating_add(raw);
                    retained_hit_count += hit_store
                        .ingest_partitions(vec![hits])
                        .map_err(|error| error.to_string())?;
                    Ok(())
                },
            )
            .map_err(io::Error::other)?;
            eprintln!(
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
        let precomputed_counts = (rayon::current_num_threads() == 1 && traits.sketch_size <= 0)
            .then(|| {
                (
                    crate::search::seed_array::SeedArray::count_partitions(
                        &query_seqs,
                        shape,
                        &reduction,
                        10,
                        0.0,
                        config.min_query_len,
                    ),
                    crate::search::seed_array::SeedArray::count_partitions(
                        &db_seqs,
                        shape,
                        &reduction,
                        10,
                        0.0,
                        config.min_query_len,
                    ),
                )
            });
        let mut shape_raw_seed_matches = 0usize;
        for chunk in 0..index_chunks.parts {
            let prepared =
                if let Some((query_counts, reference_counts)) = &precomputed_worker_histograms {
                    parallel::prepare_seed_join_partition_range_from_worker_counts_min_query_len(
                        &query_seqs,
                        &db_seqs,
                        shape,
                        &reduction,
                        complexity_cut,
                        config.min_query_len,
                        index_chunks.begin(chunk),
                        index_chunks.end(chunk),
                        &query_counts[shape_id],
                        &reference_counts[shape_id],
                    )
                } else if let Some((query_counts, reference_counts)) = &precomputed_counts {
                    parallel::prepare_seed_join_partition_range_from_counts_min_query_len(
                        &query_seqs,
                        &db_seqs,
                        shape,
                        &reduction,
                        complexity_cut,
                        config.min_query_len,
                        index_chunks.begin(chunk),
                        index_chunks.end(chunk),
                        query_counts,
                        reference_counts,
                    )
                } else {
                    parallel::prepare_seed_join_partition_range_min_query_len(
                        &query_seqs,
                        &db_seqs,
                        shape,
                        &reduction,
                        complexity_cut,
                        config.min_query_len,
                        traits.sketch_size.max(0) as usize,
                        index_chunks.begin(chunk),
                        index_chunks.end(chunk),
                    )
                };
            for &(query_id, query_pos) in prepared.masked_positions() {
                if query_id as usize >= stage2_query_set.len() {
                    continue;
                }
                if let Some(letter) = stage2_query_set
                    .get_mut(query_id as usize)
                    .get_mut(query_pos as usize)
                {
                    *letter |= SEED_MASK;
                }
            }
            low_complexity_query_positions.extend_from_slice(prepared.masked_positions());
            let query_fingerprint_seqs: Vec<&[Letter]> = (0..stage2_query_set.len())
                .map(|id| stage2_query_set.get(id))
                .collect();
            let map_candidates = |query: (u32, u32),
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
                let raw = parallel::visit_prepared_seed_hits_streaming_hamming_mode(
                    &query_fingerprint_seqs,
                    &db_fingerprint_seqs,
                    Some((&stage2_query_set, db_block.seqs())),
                    prepared,
                    hamming_filter_id,
                    usize::MAX,
                    linear_stage1_query,
                    min_length_ratio,
                    self_search,
                    map_candidates,
                    |batch| partitions.extend(batch),
                );
                retained_hit_count += hit_store.ingest_partitions(partitions)?;
                raw
            } else {
                let mut spill_error = None;
                let raw = parallel::visit_prepared_seed_hits_streaming_hamming_mode(
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
                    min_length_ratio,
                    self_search,
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
        }
        raw_seed_matches += shape_raw_seed_matches;
    }
    // C++ stage0/stage2 preserves shape + seed-partition join emission order
    // inside a query. Only group by query so range building is cheap without
    // imposing a Rust-only target/position tie-breaker on ranking boundaries.
    eprintln!(
        "Seed matches: {} (raw) -> {} (left-most)",
        raw_seed_matches, retained_hit_count,
    );

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
            .map(|(query_id, record)| fasta::FastaRecord {
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
        restore_extension_queries(&mut query_records, &stage2_query_set);
    }

    // Release seed-enumeration record storage and the temporary query copy
    // before extension.
    drop(db_records);
    drop(stage2_query_set);
    drop(db_motif_saves);
    drop(query_motif_saves);
    trim_freed_heap_pages();
    let db_ids = db_block.ids().map_err(io::Error::other)?;

    // Once seeding is complete, move (do not clone) the six translated
    // protein sequences into one source-query object. This is the shape the
    // extension code expects for translated search and avoids retaining a
    // second copy of long contigs.
    let translated_groups = if let Some(layout) = &config.translated_query_layout {
        let expected = layout.sources.len().saturating_mul(6);
        if query_records.len() != expected {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "blastx internal layout has {} protein records, expected {expected}",
                    query_records.len()
                ),
            ));
        }
        let mut records = std::mem::take(&mut query_records).into_iter();
        let mut groups = Vec::with_capacity(layout.sources.len());
        for source in &layout.sources {
            let frames = records
                .by_ref()
                .take(6)
                .map(|record| record.sequence)
                .collect();
            groups.push(TranslatedExtensionGroup {
                id: source.id.clone(),
                dna_len: source.dna_len,
                source_sequence: &source.sequence,
                mate_sequence: source
                    .mate_index
                    .and_then(|index| layout.sources.get(index))
                    .map(|mate| mate.sequence.as_slice()),
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
        .is_some_and(|format| format == "100" || format.eq_ignore_ascii_case("daa"));
    let fields = if daa_output {
        if config.outfmt.len() > 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "DAA output does not accept tabular field names",
            ));
        }
        Vec::new()
    } else if config.outfmt.is_empty() || config.outfmt[0] == "6" || config.outfmt[0] == "tab" {
        if config.outfmt.len() > 1 {
            config.outfmt[1..]
                .iter()
                .filter_map(|f| FieldId::from_name(f))
                .collect()
        } else {
            format::DEFAULT_TABULAR_FIELDS.to_vec()
        }
    } else {
        return Err(io::Error::other(format!(
            "Output format '{}' is not implemented in the native blastp pipeline. \
             Supported: 6/tab and 100/daa. Use --legacy to route to the C++ engine for other formats.",
            config.outfmt[0]
        )));
    };

    let mut text_writer: Option<Box<dyn Write>> = if daa_output {
        None
    } else if edge_output.is_some() {
        Some(Box::new(io::sink()))
    } else {
        Some(match &config.output {
            Some(path) => Box::new(BufWriter::new(std::fs::File::create(path)?)),
            None => Box::new(BufWriter::new(io::stdout())),
        })
    };
    let mut daa_writer = if daa_output {
        let path = config.output.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "DAA output requires -a/--daa or -o/--out",
            )
        })?;
        let mut writer = BufWriter::new(std::fs::File::create(path)?);
        init_daa(&mut writer)?;
        Some(writer)
    } else {
        None
    };
    let mut target_to_dict = vec![None; db_block.seqs().len()];
    let mut dict_targets = Vec::new();
    let mut aligned_queries = 0u64;
    const GAPPED_FILTER_EVALUE1: f64 = 2000.0;
    const GAPPED_FILTER_DIAG_BITS: f64 = 12.0;
    const GAPPED_FILTER_WINDOW: i32 = 200;

    let gapped_filter_evalue = traits.gapped_filter_evalue as f64;
    let cutoff_gapped1 = if gapped_filter_evalue != 0.0 {
        Some(CutoffTable2D::new(&score_matrix, GAPPED_FILTER_EVALUE1))
    } else {
        None
    };
    let cutoff_gapped2 = if gapped_filter_evalue != 0.0 {
        Some(CutoffTable2D::new(&score_matrix, gapped_filter_evalue))
    } else {
        None
    };
    let gapped_filter_diag_score = score_matrix.rawscore_int(GAPPED_FILTER_DIAG_BITS);

    // Process each query in input order. Native blastp keeps lazy target
    // masking disabled, so the C++-style extension path can read the shared
    // target block without serializing all queries.
    let process_query = |(query_idx, query_rec, query_hits): (
        usize,
        &fasta::FastaRecord,
        &[CompactHit],
    )|
     -> io::Result<QueryOutput> {
        let query = &query_rec.sequence;

        // CBS (composition-based statistics) correction per query position.
        // C++ default is comp-based-stats=1 (Hauser correction, window=40).
        // `Sequence::operator[]` strips the soft-mask bit while preserving
        // the residue, so compute Hauser from the query as stored rather
        // than converting masked positions to X.
        let query_cbs = if config.comp_based_stats.hauser() {
            crate::stats::cbs::hauser_correction(query, &score_matrix)
        } else {
            Vec::new()
        };
        let query_comp = crate::stats::cbs::compute_composition(query);
        let ungapped_cfg = UngappedStageConfig {
            hamming_ext: hamming_ext.unwrap_or_default(),
            mutual_cover,
            lin_stage1_query: linear_stage1_query,
            lin_stage1_target: linear_stage1_target,
            comp_based_stats: config.comp_based_stats,
            xdrop: score_matrix.rawscore_int(config.ungapped_xdrop_bits),
            ..UngappedStageConfig::default()
        };
        let ext_mode = extension_mode_override.unwrap_or_else(|| {
            if linear_stage1_query || linear_stage1_target {
                sensitivity::ExtensionMode::Full
            } else {
                sensitivity::default_ext_mode(config.sensitivity)
            }
        });
        let gapped_cfg = GappedScoreConfig {
            comp_based_stats_hauser: config.comp_based_stats.hauser(),
            comp_based_stats_matrix_adjust: config.comp_based_stats.matrix_adjust(),
            query_cover: config.query_cover,
            subject_cover: config.subject_cover,
            query_or_target_cover,
            no_self_hits: config.no_self_hits,
            max_evalue: config.max_evalue,
            min_id: config.min_id,
            approx_min_id: approx_min_id_override.unwrap_or(0.0),
            max_target_seqs: config.max_target_seqs,
            ext_chunk_size: config.ext_chunk_size,
            toppercent: config.toppercent,
            global_ranking_targets: config.global_ranking_targets,
            gapped_filter_evalue,
            sensitivity: config.sensitivity,
            self_: self_search,
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
        } else if daa_output {
            HspValues::TRANSCRIPT
        } else {
            HspValues::COORDS
                | HspValues::IDENT
                | HspValues::LENGTH
                | HspValues::MISMATCHES
                | HspValues::GAP_OPENINGS
        };
        let query_cbs_arg: &[Vec<i8>] = if query_cbs.is_empty() {
            &[]
        } else {
            std::slice::from_ref(&query_cbs)
        };
        let mut matches = extend_targets(
            query_idx as u32,
            &mut hits,
            std::slice::from_ref(query),
            &query_rec.id,
            query.len() as i32,
            query_cbs_arg,
            &query_comp,
            &db_block,
            &mut stat,
            Flags::NONE,
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
            GAPPED_FILTER_WINDOW,
            Option::<fn(u32, &crate::align::gapped_filter::SeedHitList) -> Vec<Match>>::None,
        );
        if config.no_self_hits {
            // Upstream applies this HSP filter after ranking/extension. Keep
            // that placement: a discarded self target is not replaced by a
            // lower-ranked target from beyond the output range.
            matches.retain(|target| {
                let target_id = target.target_block_id as usize;
                !same_sequence(query, db_block.seqs().get(target_id))
                    || db_ids.get(target_id) != query_rec.id.as_bytes()
            });
        }
        if let Some(edge_output) = &edge_output {
            let mut edges = edge_output.lock().unwrap();
            for target in &matches {
                let Some(hsp) = target.hsps.first() else {
                    continue;
                };
                edges.push(InMemorySearchEdge {
                    edge: crate::output::edge::EdgeData {
                        query: query_idx as u64,
                        target: target.target_block_id as u64,
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
            return Ok(QueryOutput::Tabular {
                bytes: Vec::new(),
                rows: matches.len() as u64,
            });
        }
        if daa_output {
            return Ok(QueryOutput::Daa(matches));
        }
        // Format each completed match directly into its query buffer. The
        // match list is already in output order, so retaining a second
        // vector of copied HSP summaries only adds allocation and traffic.
        let mut buf: Vec<u8> = Vec::new();
        let mut rows = 0u64;
        for m in matches {
            let Some(best_hsp) = m.hsps.first() else {
                continue;
            };
            let target_id = m.target_block_id as usize;

            let hsp = OutputHsp {
                score: best_hsp.score,
                evalue: best_hsp.evalue,
                bit_score: best_hsp.bit_score,
                query_range: (best_hsp.query_range.begin, best_hsp.query_range.end),
                subject_range: (best_hsp.subject_range.begin, best_hsp.subject_range.end),
                query_source_range: (best_hsp.query_range.begin, best_hsp.query_range.end),
                subject_source_range: (best_hsp.subject_range.begin, best_hsp.subject_range.end),
                frame: best_hsp.frame,
                length: best_hsp.length,
                identities: best_hsp.identities,
                mismatches: best_hsp.mismatches,
                positives: best_hsp.positives,
                gap_openings: best_hsp.gap_openings,
                gaps: best_hsp.gaps,
            };
            let target_title = std::str::from_utf8(db_ids.get(target_id)).map_err(|error| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid UTF-8 database title: {error}"),
                )
            })?;
            let subject_taxids = database_metadata
                .as_ref()
                .and_then(|metadata| metadata.subject_taxids.get(target_id))
                .map(Vec::as_slice)
                .unwrap_or(&[]);
            let taxonomy = database_metadata
                .as_ref()
                .and_then(|metadata| metadata.taxonomy.as_ref());
            format::write_tabular_row_with_query_info(
                &mut buf,
                &query_rec.id,
                target_title,
                &hsp,
                &fields,
                query.len() as i32,
                db_block.seqs().length(target_id) as i32,
                None,
                subject_taxids,
                taxonomy,
            )?;
            rows += 1;
        }
        Ok(QueryOutput::Tabular { bytes: buf, rows })
    };

    let process_translated_query =
        |(query_idx, group, query_hits): (usize, &TranslatedExtensionGroup, &[Vec<CompactHit>]),
         scratch: &mut ExtensionWorkerScratch|
         -> io::Result<QueryOutput> {
            debug_assert_eq!(group.frames.len(), 6);
            debug_assert_eq!(query_hits.len(), 6);
            let ExtensionWorkerScratch {
                hits,
                query_cbs,
                hauser_values,
            } = scratch;
            let query_cbs: &[Vec<i8>] = if config.comp_based_stats.hauser() {
                for (query, output) in group.frames.iter().zip(query_cbs.iter_mut()) {
                    crate::stats::hauser_correction::hauser_correction_into(
                        query,
                        &score_matrix,
                        output,
                        hauser_values,
                    );
                }
                query_cbs
            } else {
                &[]
            };
            let query_comp = crate::stats::cbs::compute_composition(&group.frames[0]);
            let ungapped_cfg = UngappedStageConfig {
                query_contexts: 6,
                query_translated: true,
                hamming_ext: hamming_ext.unwrap_or_default(),
                mutual_cover,
                lin_stage1_query: linear_stage1_query,
                lin_stage1_target: linear_stage1_target,
                comp_based_stats: config.comp_based_stats,
                xdrop: score_matrix.rawscore_int(config.ungapped_xdrop_bits),
                ..UngappedStageConfig::default()
            };
            let ext_mode = extension_mode_override.unwrap_or_else(|| {
                if linear_stage1_query || linear_stage1_target {
                    sensitivity::ExtensionMode::Full
                } else {
                    sensitivity::default_ext_mode(config.sensitivity)
                }
            });
            let gapped_cfg = GappedScoreConfig {
                query_contexts: 6,
                query_translated: true,
                comp_based_stats_hauser: config.comp_based_stats.hauser(),
                comp_based_stats_matrix_adjust: config.comp_based_stats.matrix_adjust(),
                query_cover: config.query_cover,
                subject_cover: config.subject_cover,
                query_or_target_cover,
                no_self_hits: config.no_self_hits,
                max_evalue: config.max_evalue,
                min_id: config.min_id,
                approx_min_id: approx_min_id_override.unwrap_or(0.0),
                max_target_seqs: config.max_target_seqs,
                ext_chunk_size: config.ext_chunk_size,
                toppercent: config.toppercent,
                global_ranking_targets: config.global_ranking_targets,
                gapped_filter_evalue,
                sensitivity: config.sensitivity,
                self_: self_search,
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
            let output_hsp_values = if daa_output {
                HspValues::TRANSCRIPT
            } else {
                HspValues::COORDS
                    | HspValues::IDENT
                    | HspValues::LENGTH
                    | HspValues::MISMATCHES
                    | HspValues::GAP_OPENINGS
            };
            let matches = extend_targets(
                query_idx as u32,
                hits,
                &group.frames,
                &group.id,
                group.dna_len,
                query_cbs,
                &query_comp,
                &db_block,
                &mut stat,
                Flags::NONE,
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
                GAPPED_FILTER_WINDOW,
                Option::<fn(u32, &crate::align::gapped_filter::SeedHitList) -> Vec<Match>>::None,
            );
            if daa_output {
                return Ok(QueryOutput::Daa(matches));
            }
            let mut buf = Vec::new();
            let mut rows = 0u64;
            for m in matches {
                let Some(best_hsp) = m.hsps.first() else {
                    continue;
                };
                let target_id = m.target_block_id as usize;
                let source_range = if best_hsp.frame >= 3 {
                    // `absolute_interval` is normalized to an ascending half-open
                    // interval. BLAST tabular coordinates are strand-oriented,
                    // so reverse contexts print the high coordinate first.
                    (
                        best_hsp.query_source_range.end - 1,
                        best_hsp.query_source_range.begin + 1,
                    )
                } else {
                    (
                        best_hsp.query_source_range.begin,
                        best_hsp.query_source_range.end,
                    )
                };
                let hsp = OutputHsp {
                    score: best_hsp.score,
                    evalue: best_hsp.evalue,
                    bit_score: best_hsp.bit_score,
                    // Tabular qcovhsp is defined in source-query coordinates for
                    // translated search, not in amino-acid frame coordinates.
                    query_range: (
                        best_hsp.query_source_range.begin,
                        best_hsp.query_source_range.end,
                    ),
                    subject_range: (best_hsp.subject_range.begin, best_hsp.subject_range.end),
                    query_source_range: source_range,
                    subject_source_range: (
                        best_hsp.subject_range.begin,
                        best_hsp.subject_range.end,
                    ),
                    frame: crate::basic::translate::Frame::from_index(best_hsp.frame)
                        .signed_frame(),
                    length: best_hsp.length,
                    identities: best_hsp.identities,
                    mismatches: best_hsp.mismatches,
                    positives: best_hsp.positives,
                    gap_openings: best_hsp.gap_openings,
                    gaps: best_hsp.gaps,
                };
                let target_title = std::str::from_utf8(db_ids.get(target_id)).map_err(|error| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid UTF-8 database title: {error}"),
                    )
                })?;
                format::write_tabular_row_with_query_info(
                    &mut buf,
                    &group.id,
                    target_title,
                    &hsp,
                    &fields,
                    group.dna_len,
                    db_block.seqs().length(target_id) as i32,
                    group.mate_sequence,
                    database_metadata
                        .as_ref()
                        .and_then(|metadata| metadata.subject_taxids.get(target_id))
                        .map(Vec::as_slice)
                        .unwrap_or(&[]),
                    database_metadata
                        .as_ref()
                        .and_then(|metadata| metadata.taxonomy.as_ref()),
                )?;
                rows += 1;
            }
            Ok(QueryOutput::Tabular { bytes: buf, rows })
        };
    let mut total_alignments = 0u64;
    // Preserve input order without retaining every query's formatted output.
    // A moderately sized batch gives Rayon enough work to balance variable
    // protein lengths while keeping buffered output proportional to threads.
    let mut write_query_range = |range_begin: usize,
                                 records: &[fasta::FastaRecord],
                                 hits_by_query: &[Vec<CompactHit>]|
     -> io::Result<()> {
        let (output_batch_size, trim_between_batches) =
            output_query_batch_size(records, rayon::current_num_threads());
        for (batch_idx, query_batch) in records.chunks(output_batch_size).enumerate() {
            let local_begin = batch_idx * output_batch_size;
            let query_begin = range_begin + local_begin;
            let hit_batch = &hits_by_query[local_begin..local_begin + query_batch.len()];
            let batch_output: Vec<QueryOutput> = if rayon::current_num_threads() == 1 {
                query_batch
                    .iter()
                    .zip(hit_batch)
                    .enumerate()
                    .map(|(offset, (record, hits))| {
                        process_query((query_begin + offset, record, hits))
                    })
                    .collect::<io::Result<Vec<_>>>()?
            } else {
                query_batch
                    .par_iter()
                    .zip(hit_batch.par_iter())
                    .enumerate()
                    .map(|(offset, (record, hits))| {
                        process_query((query_begin + offset, record, hits))
                    })
                    .collect::<io::Result<Vec<_>>>()?
            };
            for (record, output) in query_batch.iter().zip(batch_output) {
                match output {
                    QueryOutput::Tabular { bytes, rows } => {
                        total_alignments += rows;
                        text_writer
                            .as_mut()
                            .expect("tabular writer")
                            .write_all(&bytes)?;
                    }
                    QueryOutput::Daa(matches) => {
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
                        if count != 0 {
                            aligned_queries += 1;
                            total_alignments += count;
                            daa_writer.as_mut().expect("DAA writer").write_all(&bytes)?;
                        }
                    }
                }
            }
            if trim_between_batches && local_begin.saturating_add(query_batch.len()) < records.len()
            {
                trim_freed_heap_pages();
            }
        }
        Ok(())
    };

    if let Some(groups) = translated_groups.as_ref() {
        // Release the normal-query writer closure's borrow before selecting
        // the six-context path.
        drop(write_query_range);
        let mut write_translated_range = |source_begin: usize,
                                          source_groups: &[TranslatedExtensionGroup],
                                          hits_by_context: &[Vec<CompactHit>]|
         -> io::Result<()> {
            const MIN_TASK_TRACE_POINTS: usize = 1024;
            let partitions = translated_hit_partitions(
                hits_by_context,
                source_groups.len(),
                6,
                MIN_TASK_TRACE_POINTS,
            );
            let process_partition = |partition: usize, scratch: &mut ExtensionWorkerScratch| {
                partitions[partition]
                    .clone()
                    .map(|offset| {
                        let begin = offset * 6;
                        process_translated_query(
                            (
                                source_begin + offset,
                                &source_groups[offset],
                                &hits_by_context[begin..begin + 6],
                            ),
                            scratch,
                        )
                    })
                    .collect::<io::Result<Vec<QueryOutput>>>()
            };
            let partition_output: Vec<io::Result<Vec<QueryOutput>>> =
                if rayon::current_num_threads() == 1 || partitions.len() <= 1 {
                    let mut scratch = ExtensionWorkerScratch::default();
                    (0..partitions.len())
                        .map(|partition| process_partition(partition, &mut scratch))
                        .collect()
                } else {
                    use std::sync::atomic::{AtomicUsize, Ordering};
                    use std::sync::OnceLock;

                    let next = AtomicUsize::new(0);
                    let output: Vec<OnceLock<io::Result<Vec<QueryOutput>>>> =
                        (0..partitions.len()).map(|_| OnceLock::new()).collect();
                    let worker_count = rayon::current_num_threads().min(partitions.len());
                    rayon::scope(|scope| {
                        for _ in 0..worker_count {
                            let next = &next;
                            let output = &output;
                            let process_partition = &process_partition;
                            scope.spawn(move |_| {
                                let mut scratch = ExtensionWorkerScratch::default();
                                loop {
                                    let partition = next.fetch_add(1, Ordering::Relaxed);
                                    if partition >= output.len() {
                                        break;
                                    }
                                    let result = process_partition(partition, &mut scratch);
                                    let was_empty = output[partition].set(result).is_ok();
                                    debug_assert!(was_empty);
                                }
                            });
                        }
                    });
                    output
                        .into_iter()
                        .map(|slot| {
                            slot.into_inner()
                                .expect("translated partition worker left an empty output slot")
                        })
                        .collect()
                };
            for (partition_index, partition) in partition_output.into_iter().enumerate() {
                let outputs = partition?;
                for (offset, output) in partitions[partition_index].clone().zip(outputs) {
                    match output {
                        QueryOutput::Tabular { bytes, rows } => {
                            total_alignments += rows;
                            text_writer
                                .as_mut()
                                .expect("tabular writer")
                                .write_all(&bytes)?;
                        }
                        QueryOutput::Daa(matches) => {
                            let group = &source_groups[offset];
                            let mut bytes = Vec::new();
                            let count = encode_daa_query(
                                &mut bytes,
                                &group.id,
                                group.source_sequence,
                                SequenceType::Nucleotide,
                                &matches,
                                &mut target_to_dict,
                                &mut dict_targets,
                            )?;
                            if count != 0 {
                                aligned_queries += 1;
                                total_alignments += count;
                                daa_writer.as_mut().expect("DAA writer").write_all(&bytes)?;
                            }
                        }
                    }
                }
            }
            Ok(())
        };

        if let Some(mut memory) = hit_store.memory.take() {
            for bin in 0..memory.bins.len() {
                let (hits_by_context, begin, end) = memory.take_grouped_bin(bin);
                let begin_idx = begin as usize;
                let end_idx = (end as usize).min(groups.len() * 6);
                debug_assert_eq!(begin_idx % 6, 0);
                debug_assert_eq!(end_idx % 6, 0);
                write_translated_range(
                    begin_idx / 6,
                    &groups[begin_idx / 6..end_idx / 6],
                    &hits_by_context[..end_idx - begin_idx],
                )?;
            }
        } else {
            let disk = hit_store
                .disk
                .as_mut()
                .ok_or_else(|| io::Error::other("missing hit storage"))?;
            disk.try_finish_writing().map_err(io::Error::other)?;
            let mut load_pending = disk.load_grouped();
            while load_pending {
                let Some((mut bin_hits, begin, mut end)) = disk
                    .try_retrieve_grouped_owned()
                    .map_err(io::Error::other)?
                else {
                    break;
                };
                load_pending = disk.load_grouped();
                // Upstream intended `HitBuffer::load(max_size)` to combine
                // adjacent bins, but its current loop condition makes that
                // read exactly one bin. A translated bin often contains only
                // two or three source queries, which strands alignment
                // workers behind a barrier. Restore the intended read-ahead,
                // bounded by live decoded-hit storage. Inner vectors are
                // moved into the window; individual hits are never copied.
                const MAX_DECODED_HIT_WINDOW: usize = 64 * 1024 * 1024;
                let mut decoded_bytes = bin_hits
                    .iter()
                    .map(|hits| hits.capacity() * std::mem::size_of::<CompactHit>())
                    .sum::<usize>();
                while load_pending && decoded_bytes < MAX_DECODED_HIT_WINDOW {
                    let Some((mut next_hits, next_begin, next_end)) = disk
                        .try_retrieve_grouped_owned()
                        .map_err(io::Error::other)?
                    else {
                        break;
                    };
                    debug_assert_eq!(next_begin, end);
                    decoded_bytes = decoded_bytes.saturating_add(
                        next_hits
                            .iter()
                            .map(|hits| hits.capacity() * std::mem::size_of::<CompactHit>())
                            .sum::<usize>(),
                    );
                    bin_hits.append(&mut next_hits);
                    end = next_end;
                    load_pending = disk.load_grouped();
                }
                let begin_idx = begin as usize;
                let end_idx = (end as usize).min(groups.len() * 6);
                debug_assert_eq!(begin_idx % 6, 0);
                debug_assert_eq!(end_idx % 6, 0);
                write_translated_range(
                    begin_idx / 6,
                    &groups[begin_idx / 6..end_idx / 6],
                    &bin_hits[..end_idx - begin_idx],
                )?;
                disk.recycle_grouped_hits(bin_hits);
            }
            disk.free_buffer();
            eprintln!(
                "Hit buffer: read {:.1} MiB of compressed spill data",
                disk.total_disk_size() as f64 / (1024.0 * 1024.0)
            );
        }
    } else if let Some(mut memory) = hit_store.memory.take() {
        for bin in 0..memory.bins.len() {
            let (hits_by_query, begin, end) = memory.take_grouped_bin(bin);
            let begin_idx = begin as usize;
            let end_idx = (end as usize).min(query_records.len());
            write_query_range(
                begin_idx,
                &query_records[begin_idx..end_idx],
                &hits_by_query[..end_idx - begin_idx],
            )?;
        }
    } else {
        let disk = hit_store
            .disk
            .as_mut()
            .ok_or_else(|| io::Error::other("missing hit storage"))?;
        disk.try_finish_writing().map_err(io::Error::other)?;
        let mut load_pending = disk.load_grouped();
        while load_pending {
            let Some((bin_hits, begin, end)) = disk
                .try_retrieve_grouped_owned()
                .map_err(io::Error::other)?
            else {
                break;
            };
            // Start decoding the next spill bin before alignment of this bin.
            load_pending = disk.load_grouped();
            let begin_idx = begin as usize;
            let end_idx = (end as usize).min(query_records.len());
            write_query_range(
                begin_idx,
                &query_records[begin_idx..end_idx],
                &bin_hits[..end_idx - begin_idx],
            )?;
            disk.recycle_grouped_hits(bin_hits);
        }
        disk.free_buffer();
        eprintln!(
            "Hit buffer: read {:.1} MiB of compressed spill data",
            disk.total_disk_size() as f64 / (1024.0 * 1024.0)
        );
    }
    if let Some(writer) = text_writer.as_mut() {
        writer.flush()?;
    }
    if let Some(writer) = daa_writer.as_mut() {
        let mut refs = Vec::with_capacity(dict_targets.len());
        for &target_id in &dict_targets {
            let title = std::str::from_utf8(db_ids.get(target_id))
                .map_err(|error| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid UTF-8 database title: {error}"),
                    )
                })?
                .to_owned();
            refs.push((title, db_block.seqs().length(target_id) as u32));
        }
        let dictionary = NativeDaaDictionary {
            sequence_count: db_block.seqs().len() as u64,
            refs,
        };
        let metadata = DaaRunMetadata {
            db_letters: score_matrix.db_letters(),
            gap_open: score_matrix.gap_open(),
            gap_extend: score_matrix.gap_extend(),
            // These configuration values are meaningful only for blastn but
            // upstream writes their defaults into every DAA header.
            reward: 2,
            penalty: -3,
            k: score_matrix.k(),
            lambda: score_matrix.lambda(),
            max_evalue: config.max_evalue,
            matrix: score_matrix.name().to_owned(),
            mode: if translated_groups.is_some() {
                crate::basic::value::AlignMode::BLASTX as u32
            } else {
                crate::basic::value::AlignMode::BLASTP as u32
            },
            aligned_queries,
        };
        finish_daa_from_sequence_file(writer, &dictionary, &metadata)?;
        writer.flush()?;
    }

    let elapsed = start.elapsed();
    eprintln!(
        "Reported {} alignments in {:.1}s",
        total_alignments,
        elapsed.as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::{MASK_LETTER, SEED_MASK};

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
                taxonmap: Some(mapping.to_string_lossy().into_owned()),
                taxonnodes: Some(nodes.to_string_lossy().into_owned()),
                taxonnames: Some(names.to_string_lossy().into_owned()),
                no_parse_seqids: false,
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
        let mut records = vec![
            fasta::FastaRecord {
                id: "q0".into(),
                sequence: vec![MASK_LETTER; 4],
            },
            fasta::FastaRecord {
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
    fn adaptive_hit_store_migrates_existing_hits_and_keeps_writing_to_disk() {
        let tmpdir = std::env::temp_dir().join(format!(
            "diamond-adaptive-spill-test-{}",
            std::process::id()
        ));
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
        assert_eq!(stage2_query_bounds(232, 3), (0, 51));
        assert_eq!(stage2_query_bounds(232, 100), (52, 148));
        assert_eq!(stage2_query_bounds(120, 100), (52, 120));
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
            query_files: vec![query_path.to_string_lossy().into_owned()],
            database: db_path.to_string_lossy().into_owned(),
            output: Some(output_path.to_string_lossy().into_owned()),
            matrix: "blosum62".to_owned(),
            gap_open: 11,
            gap_extend: 1,
            max_evalue: 0.001,
            max_target_seqs: 25,
            ext_chunk_size: 0,
            toppercent: None,
            global_ranking_targets: 0,
            min_id: 0.0,
            threads: 1,
            outfmt: vec![
                "6".to_owned(),
                "qseqid".to_owned(),
                "sseqid".to_owned(),
                "pident".to_owned(),
                "length".to_owned(),
                "mismatch".to_owned(),
                "gapopen".to_owned(),
                "qstart".to_owned(),
                "qend".to_owned(),
                "sstart".to_owned(),
                "send".to_owned(),
                "evalue".to_owned(),
                "bitscore".to_owned(),
            ],
            sensitivity: Sensitivity::Default,
            masking: MaskingMode::None,
            motif_masking: String::new(),
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
