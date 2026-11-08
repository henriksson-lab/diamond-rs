//! Concrete native nucleotide reference-block search.
//!
//! Reference storage, minimizer index, Karlin model and both extension engines
//! are owned for the full block lifetime. Command orchestration can retain
//! source query records independently for reporting and iterative passes.
use super::build_score::Blastn_Score;
use super::chain::ChainingParameters;
use super::dna_index::Index;
use super::extension::{DnaExtensionAlgo, ExtensionConfig};
use super::extension_chain::{chaining_and_extension, compute_chains, ChainExtensionConfig};
use super::seed_set_dna::seed_lookup;
use super::setup::{self, DnaSearchConfig};
use crate::align::hsp::Match;
use crate::basic::packed_loc::PackedLoc;
use crate::basic::reduction::Reduction;
use crate::basic::seed::{seed_partition, seed_partition_offset, seedp_mask};
use crate::basic::seed_iterator::MinimizerIterator;
use crate::basic::shape::Shape;
use crate::basic::value::{Letter, NUCLEOTIDE_ALPHABET};
use crate::config::Sensitivity;
use crate::data::seed_array::SeedArrayEntry;
use crate::data::sequence_set::LetterStringSet;

/// Block::length_sorted orders by (length, previous block ID), descending.
/// Move the sequence/title owner and its Source OID together. Source's
/// replacement Block omits quality storage; qqual/full_qqual then print '*'.
pub(crate) fn sort_records_by_length<T>(
    records: &mut Vec<crate::data::fasta::FastaRecord<T>>,
    oids: &mut Vec<u64>,
) -> Result<(), String> {
    if records.len() != oids.len() {
        return Err("Nucleotide record/OID permutation size mismatch".into());
    }
    let mut ordered = records
        .drain(..)
        .zip(oids.drain(..))
        .enumerate()
        .collect::<Vec<_>>();
    ordered.sort_unstable_by(|(left_id, (left, _)), (right_id, (right, _))| {
        right
            .sequence
            .len()
            .cmp(&left.sequence.len())
            .then_with(|| right_id.cmp(left_id))
    });
    for (_, (mut record, oid)) in ordered {
        record.quality = None;
        records.push(record);
        oids.push(oid);
    }
    Ok(())
}

thread_local! {
    static ACTIVE_TIMERS: std::cell::RefCell<Vec<std::sync::Arc<std::sync::Mutex<super::timer::ExtensionTimer>>>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// One Source TotalTime owner per command, with scoped worker access. The
/// destructor prints even when a search returns an error, as Source does.
pub(crate) struct NativeDnaTiming(std::sync::Arc<std::sync::Mutex<super::timer::ExtensionTimer>>);
pub(crate) struct DnaTimingScope;
impl NativeDnaTiming {
    pub(crate) fn new() -> Self {
        Self(std::sync::Arc::new(std::sync::Mutex::new(
            super::timer::ExtensionTimer::new(),
        )))
    }
    pub(crate) fn enter(&self) -> DnaTimingScope {
        ACTIVE_TIMERS.with(|timers| timers.borrow_mut().push(self.0.clone()));
        DnaTimingScope
    }
}
impl Drop for DnaTimingScope {
    fn drop(&mut self) {
        ACTIVE_TIMERS.with(|timers| {
            timers.borrow_mut().pop();
        });
    }
}
impl Drop for NativeDnaTiming {
    fn drop(&mut self) {
        let timer = self
            .0
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .clone();
        super::timer::TotalTime { timer }.print();
    }
}
pub(crate) fn record_time(operation: usize, duration: std::time::Duration) {
    ACTIVE_TIMERS.with(|timers| {
        if let Some(timer) = timers.borrow().last() {
            timer
                .lock()
                .unwrap_or_else(|error| error.into_inner())
                .update(operation as i32, duration);
        }
    });
}

/// WITH_DNA option values. Gap penalties remain supplied by the shared command
/// configuration: Source run/config.cpp replaces each -1 sentinel by 5/2
/// before constructing Blastn_Score. This engine receives the resulting model.
#[derive(Debug, Clone)]
pub struct DnaOptions {
    pub zdrop: i32,
    pub repetitive_cutoff: f64,
    pub algorithm: DnaExtensionAlgo,
    pub chaining_out: bool,
    pub align_long_reads: bool,
    pub best_hsp_only: bool,
    pub chain_pen_gap_scale: f64,
    pub chain_pen_skip_scale: f64,
    pub penalty: i32,
    pub reward: i32,
    pub chain_align_cutoff: Option<f64>,
    pub min_chain_score: Option<i32>,
    pub max_overlap_extension: Option<f64>,
    pub zdrop_extension: i32,
    pub zdrop_global: i32,
    pub band_extension: i32,
    pub band_global: i32,
}
impl Default for DnaOptions {
    fn default() -> Self {
        Self {
            zdrop: 40,
            repetitive_cutoff: 0.0002,
            algorithm: DnaExtensionAlgo::Ksw,
            chaining_out: false,
            align_long_reads: false,
            best_hsp_only: false,
            chain_pen_gap_scale: 0.8,
            chain_pen_skip_scale: 0.0,
            penalty: -3,
            reward: 2,
            chain_align_cutoff: None,
            min_chain_score: None,
            max_overlap_extension: None,
            zdrop_extension: 80,
            zdrop_global: 150,
            band_extension: 40,
            band_global: 40,
        }
    }
}

pub struct NativeDnaBlock {
    targets: LetterStringSet,
    chain_targets: Option<Vec<Vec<Letter>>>,
    index: Index,
    shape: Shape,
    reduction: Reduction,
    statistics: Blastn_Score,
    options: DnaOptions,
    setup: DnaSearchConfig,
    minimizer_window: i32,
    max_evalue: f64,
}
impl NativeDnaBlock {
    /// Build the actual Source minimizer index. `min_query_len` applies to
    /// reference enumeration with global context=1. DNA query lookup bypasses
    /// query histograms (HitIterator::single_query), so it has no such filter.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        targets: LetterStringSet,
        statistics: Blastn_Score,
        sensitivity: Sensitivity,
        options: DnaOptions,
        minimizer_window: Option<i32>,
        min_query_len: i32,
        threads: i32,
        index_chunks: i32,
        max_evalue: f64,
    ) -> Result<Self, String> {
        Self::new_with_seed_shapes(
            targets,
            statistics,
            sensitivity,
            options,
            minimizer_window,
            min_query_len,
            threads,
            index_chunks,
            max_evalue,
            &[],
            0,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new_with_seed_shapes(
        targets: LetterStringSet,
        statistics: Blastn_Score,
        sensitivity: Sensitivity,
        options: DnaOptions,
        minimizer_window: Option<i32>,
        min_query_len: i32,
        threads: i32,
        index_chunks: i32,
        max_evalue: f64,
        shape_mask: &[String],
        shape_count: u32,
    ) -> Result<Self, String> {
        let shape_mask = shape_mask
            .iter()
            .map(|code| code.as_bytes().to_vec())
            .collect::<Vec<_>>();
        Self::new_with_seed_shapes_bytes(
            targets,
            statistics,
            sensitivity,
            options,
            minimizer_window,
            min_query_len,
            threads,
            index_chunks,
            max_evalue,
            &shape_mask,
            shape_count,
        )
    }

    /// Canonical Source shape-byte constructor. String convenience callers
    /// share the same index and storage validation through the adapter above.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_seed_shapes_bytes(
        targets: LetterStringSet,
        statistics: Blastn_Score,
        sensitivity: Sensitivity,
        options: DnaOptions,
        minimizer_window: Option<i32>,
        min_query_len: i32,
        threads: i32,
        index_chunks: i32,
        max_evalue: f64,
        shape_mask: &[Vec<u8>],
        shape_count: u32,
    ) -> Result<Self, String> {
        if threads < 1 || index_chunks < 1 {
            return Err("Invalid nucleotide index thread/chunk count".into());
        }
        if !options.repetitive_cutoff.is_finite()
            || !(0.0..=1.0).contains(&options.repetitive_cutoff)
        {
            return Err("Invalid repetitive minimizer fraction".into());
        }
        let traits =
            setup::sensitivity_traits(sensitivity).ok_or("Unsupported nucleotide sensitivity")?;
        let reduction = Reduction::new("A C G T", NUCLEOTIDE_ALPHABET);
        let codes = if shape_mask.is_empty() {
            setup::shape_codes(sensitivity)
                .ok_or("Unsupported nucleotide seed shapes")?
                .iter()
                .map(|code| code.as_bytes().to_vec())
                .collect::<Vec<_>>()
        } else {
            shape_mask.to_vec()
        };
        let selected = if shape_count == 0 {
            codes.len()
        } else {
            (shape_count as usize).min(codes.len())
        };
        if selected == 0 || selected > crate::basic::consts::MAX_SHAPES {
            return Err("Nucleotide seed shape count exceeds supported storage".into());
        }
        for code in codes.iter().take(selected) {
            let code = crate::util::sequence::c_string_bytes(code);
            let weight = code.iter().filter(|&&c| c == b'1').count();
            if code.len() > 32 || !(2..=31).contains(&weight) {
                // Source Shape asserts these storage/coordinate bounds; its
                // midpoint index underflows for weights0/1 and 4^32 overflows.
                return Err("Nucleotide seed shape exceeds supported length/weight".into());
            }
        }
        let shapes = crate::basic::shape_config::ShapeConfig::from_codes_bytes(
            &codes,
            shape_count,
            &reduction,
        )?;
        // DNA Index and seed_lookup explicitly use shape0, after validating
        // all selected shapes for the uniform-weight Source contract.
        let shape = shapes[0];
        let minimizer_window = minimizer_window
            .filter(|&n| n != 0)
            .unwrap_or(traits.minimizer_window);
        if minimizer_window < 1 {
            return Err("Nucleotide minimizer window must be positive".into());
        }
        let mut search_setup = DnaSearchConfig {
            chain_fraction_align_option: options.chain_align_cutoff.filter(|&n| n != 0.0),
            min_chain_score_option: options.min_chain_score.filter(|&n| n != 0),
            max_overlap_extension_option: options.max_overlap_extension.filter(|&n| n != 0.0),
            chain_pen_gap_scale: options.chain_pen_gap_scale,
            chain_pen_skip_scale: options.chain_pen_skip_scale,
            ..DnaSearchConfig::default()
        };
        if !setup::setup_search(sensitivity, &mut search_setup, shape.length) {
            return Err("Unsupported nucleotide sensitivity".into());
        }
        let bits =
            crate::search::setup::seedp_bits(shape.weight, threads, index_chunks, &reduction);
        let count = 1usize
            .checked_shl(bits as u32)
            .ok_or("Nucleotide seed partition count overflow")?;
        let mut partitions = Vec::new();
        partitions
            .try_reserve_exact(count)
            .map_err(|_| "Unable to allocate nucleotide index partitions")?;
        partitions.resize_with(count, Vec::new);
        let mask = seedp_mask(bits);
        for id in 0..targets.size() {
            let sequence = targets.get(id as usize);
            if sequence.len() > i32::MAX as usize {
                return Err("Nucleotide reference exceeds Loc coordinate range".into());
            }
            if sequence.len() < shape.length as usize
                || (sequence.len() as i64) < i64::from(min_query_len)
            {
                continue;
            }
            // Reference enumeration reduces first; DNA query lookup keeps raw
            // codes. In particular Source's otherwise unmapped N4 becomes0
            // only on this reference side. Preserve that asymmetric mapping.
            let reduced = sequence
                .iter()
                .map(|&letter| {
                    reduction.reduce(letter & crate::basic::value::LETTER_MASK) as Letter
                })
                .collect::<Vec<_>>();
            let mut it = MinimizerIterator::new(&reduced, &shape, minimizer_window, &reduction);
            while it.good() {
                let seed = it.get();
                let position = targets.position(id, it.pos());
                if position < 0 || position as u64 > 0xffffffffff {
                    return Err("Nucleotide seed position exceeds PackedLoc".into());
                }
                partitions[seed_partition(seed, mask) as usize].push(SeedArrayEntry::new(
                    seed_partition_offset(seed, bits as u64),
                    PackedLoc::new(position as u64),
                ));
                it.increment();
            }
        }
        let index = Index::new(partitions, bits, options.repetitive_cutoff);
        let chain_targets = (options.chaining_out || options.align_long_reads).then(|| {
            (0..targets.size())
                .map(|id| targets.get(id as usize).to_vec())
                .collect()
        });
        Ok(Self {
            targets,
            chain_targets,
            index,
            shape,
            reduction,
            statistics,
            options,
            setup: search_setup,
            minimizer_window,
            max_evalue,
        })
    }
    pub fn statistics(&self) -> &Blastn_Score {
        &self.statistics
    }
    pub fn targets(&self) -> &LetterStringSet {
        &self.targets
    }

    /// Return original Blastn_Score HSP statistics. No protein score recomputation
    /// occurs; errors from KSW/WFA/chaining propagate to the command caller.
    pub fn search(&self, query: &[Letter]) -> Result<Vec<Match>, String> {
        if query.len() > i32::MAX as usize {
            return Err("Nucleotide query exceeds Loc coordinate range".into());
        }
        if query
            .iter()
            .any(|&l| !(0..=4).contains(&(l & crate::basic::value::LETTER_MASK)))
        {
            return Err("Invalid encoded nucleotide query".into());
        }
        if let Some(targets) = &self.chain_targets {
            let reverse = crate::basic::sequence_utils::reverse_complement(query);
            let mut params = ChainingParameters::new(
                self.setup.chain_pen_gap as f32,
                self.setup.chain_pen_skip as f32,
                self.setup.min_chain_score,
                self.setup.max_overlap_extension as f32,
            );
            params.best_hsp_only = self.options.best_hsp_only;
            let seed_started = std::time::Instant::now();
            let forward_hits = seed_lookup(
                query,
                &self.targets,
                &self.index,
                self.minimizer_window,
                &self.shape,
                &self.reduction,
            );
            record_time(5, seed_started.elapsed());
            let mut chains = compute_chains(
                forward_hits,
                query,
                targets,
                self.shape.length,
                false,
                &params,
            );
            let seed_started = std::time::Instant::now();
            let reverse_hits = seed_lookup(
                &reverse,
                &self.targets,
                &self.index,
                self.minimizer_window,
                &self.shape,
                &self.reduction,
            );
            record_time(5, seed_started.elapsed());
            chains.extend(compute_chains(
                reverse_hits,
                &reverse,
                targets,
                self.shape.length,
                true,
                &params,
            ));
            chaining_and_extension(
                chains,
                targets,
                self.setup.chain_fraction_align as f32,
                self.options.chaining_out,
                self.shape.length,
                query,
                &reverse,
                &params,
                ChainExtensionConfig {
                    score_builder: &self.statistics,
                    algorithm: self.options.algorithm,
                    zdrop_extension: self.options.zdrop_extension,
                    zdrop_global: self.options.zdrop_global,
                    band_extension: self.options.band_extension,
                    band_global: self.options.band_global,
                    max_evalue: self.max_evalue,
                    best_hsp_only: self.options.best_hsp_only,
                },
            )
        } else {
            // The infallible legacy callback is unreachable for the short-read
            // branch; chained execution above has its own fallible call graph.
            let unused = |_: &[Letter], _: &[Letter]| Vec::new();
            let cfg = ExtensionConfig {
                targets: &self.targets,
                index: &self.index,
                seed_shape: &self.shape,
                reduction: &self.reduction,
                score_builder: &self.statistics,
                minimizer_window: self.minimizer_window,
                kmer_size: self.shape.length,
                max_evalue: self.max_evalue,
                dna_extension: self.options.algorithm,
                zdrop: self.options.zdrop,
                chaining_out: false,
                align_long_reads: false,
                chaining_extension: &unused,
                timing: Some(&record_time),
            };
            super::extension::extend(&cfg, query).map(|(matches, _)| matches)
        }
    }
}
