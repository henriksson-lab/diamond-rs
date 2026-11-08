//! Source's command-wide AUTO selection and per-query-iteration seed cache.
use crate::config::Sensitivity;
use crate::data::seed_set::HashedSeedSet;
use crate::run::config::Algo;
use std::io;
use std::sync::{Arc, Mutex};

#[derive(Debug)]
struct Runtime {
    selected: Algo,
    final_sensitivity: Sensitivity,
    query_seeds: Option<Arc<HashedSeedSet>>,
}
#[derive(Debug, Clone)]
pub struct AlgorithmControl {
    pub requested: Algo,
    runtime: Option<Arc<Mutex<Runtime>>>,
}
impl Default for AlgorithmControl {
    fn default() -> Self {
        Self {
            requested: Algo::Auto,
            runtime: None,
        }
    }
}
// Runtime state/caches are not user configuration and do not affect config
// equality. Clones within one search intentionally share that runtime owner.
impl PartialEq for AlgorithmControl {
    fn eq(&self, other: &Self) -> bool {
        self.requested == other.requested
    }
}
impl AlgorithmControl {
    /// Clustering writes Source's global algorithm before each search phase.
    /// Update the resolved state as well as the requested scalar; an earlier
    /// AUTO/cache decision must not survive that explicit assignment.
    pub(crate) fn force_source_algorithm(
        &mut self,
        sensitivity: Sensitivity,
        requested: Algo,
    ) -> io::Result<()> {
        self.requested = requested;
        self.begin_search(sensitivity, requested);
        let mut state = self
            .runtime
            .as_ref()
            .unwrap()
            .lock()
            .map_err(|_| io::Error::other("Algorithm state lock poisoned."))?;
        state.selected = requested;
        state.final_sensitivity = sensitivity;
        state.query_seeds = None;
        Ok(())
    }
    pub(crate) fn begin_source_search(
        &mut self,
        sensitivity: Sensitivity,
        requested: Algo,
    ) -> io::Result<()> {
        self.begin_search(sensitivity, requested);
        let mut state = self
            .runtime
            .as_ref()
            .unwrap()
            .lock()
            .map_err(|_| io::Error::other("Algorithm state lock poisoned."))?;
        if state.selected == Algo::Auto {
            state.final_sensitivity = sensitivity;
        }
        state.query_seeds = None;
        Ok(())
    }
    pub(crate) fn begin_search(&mut self, sensitivity: Sensitivity, requested: Algo) {
        if self.runtime.is_none() {
            self.runtime = Some(Arc::new(Mutex::new(Runtime {
                selected: requested,
                final_sensitivity: sensitivity,
                query_seeds: None,
            })));
        }
    }
    pub(crate) fn source_request(
        &self,
        query_indexed: bool,
        contiguous: bool,
        swipe: bool,
    ) -> Algo {
        if swipe {
            Algo::DoubleIndexed
        } else if contiguous {
            Algo::CtgSeed
        } else if query_indexed {
            Algo::QueryIndexed
        } else {
            self.requested
        }
    }
    pub(crate) fn begin_round(&self) -> io::Result<()> {
        if let Some(runtime) = &self.runtime {
            runtime
                .lock()
                .map_err(|_| io::Error::other("Algorithm state lock poisoned."))?
                .query_seeds = None;
        }
        Ok(())
    }
    pub(crate) fn selected(&self) -> io::Result<Algo> {
        match &self.runtime {
            Some(runtime) => Ok(runtime
                .lock()
                .map_err(|_| io::Error::other("Algorithm state lock poisoned."))?
                .selected),
            None => Ok(self.requested),
        }
    }
    pub(crate) fn cached_seeds(&self) -> io::Result<Option<Arc<HashedSeedSet>>> {
        match &self.runtime {
            Some(runtime) => Ok(runtime
                .lock()
                .map_err(|_| io::Error::other("Algorithm state lock poisoned."))?
                .query_seeds
                .clone()),
            None => Ok(None),
        }
    }
    pub(crate) fn auto_candidate(
        &self,
        query_letters: u64,
        database_letters: u64,
        target_indexed: bool,
        swipe: bool,
        minimizer: i32,
        sketch: i32,
    ) -> io::Result<bool> {
        let runtime = self
            .runtime
            .as_ref()
            .ok_or_else(|| io::Error::other("Algorithm state was not initialized."))?;
        let state = runtime
            .lock()
            .map_err(|_| io::Error::other("Algorithm state lock poisoned."))?;
        Ok(
            crate::search::sensitivity::get_traits(state.final_sensitivity).query_indexed
                && query_letters <= crate::run::double_indexed::MAX_INDEX_QUERY_SIZE
                && database_letters >= crate::run::double_indexed::MIN_QUERY_INDEXED_DB_SIZE
                && !target_indexed
                && !swipe
                && minimizer == 0
                && sketch == 0,
        )
    }
    pub(crate) fn resolve(&self, mode: Algo, seeds: Option<Arc<HashedSeedSet>>) -> io::Result<()> {
        let runtime = self
            .runtime
            .as_ref()
            .ok_or_else(|| io::Error::other("Algorithm state was not initialized."))?;
        let mut state = runtime
            .lock()
            .map_err(|_| io::Error::other("Algorithm state lock poisoned."))?;
        state.selected = mode;
        state.query_seeds = seeds;
        Ok(())
    }
}

/// Enumerate the real hard-masked query arena; no translated contexts are
/// copied or translated again to make an AUTO estimate.
pub(crate) fn build_query_seeds(
    sequences: &mut crate::data::sequence_set::SequenceSet,
    translated: Option<&mut crate::data::block::Block>,
    config: &super::blastp::BlastpConfig,
    matrix: &crate::stats::score_matrix::ScoreMatrix,
    setup: &crate::search::sensitivity::SetupSearchResult,
    soft_override: crate::masking::MaskingAlgo,
    frequency_masking: bool,
) -> io::Result<Arc<HashedSeedSet>> {
    use crate::masking::{Masking, MaskingAlgo, MaskingConfig, MaskingTable};
    let mut temporary = crate::data::block::Block::new();
    let block = translated.unwrap_or(&mut temporary);
    block.set_max_motif_len(config.tuning.max_motif_len);
    block.put_seqs(std::mem::take(sequences));
    let result = (|| {
        let automatic = crate::search::sensitivity::soft_masking_algo_bytes(
            &crate::search::sensitivity::get_traits(config.sensitivity),
            &config.motif_masking,
            false,
            frequency_masking,
        )
        .map_err(io::Error::other)?;
        let bits = automatic as u32 | soft_override as u32 | config.soft_masking as u32;
        // The installed table contains the full union; a non-NONE selector
        // makes Block enumeration apply/remove that exact table on BOTH hash
        // passes, including iteration's present (possibly all-false) skip mask.
        let soft_masking = if config.soft_masking != MaskingAlgo::None {
            config.soft_masking
        } else if soft_override != MaskingAlgo::None {
            soft_override
        } else {
            automatic
        };
        let owned_masking = config.tuning.runtime.masking().is_none().then(|| {
            Masking::new(
                matrix,
                MaskingConfig {
                    tantan_min_mask_prob: config.tuning.tantan_min_mask_prob as f32,
                    max_motif_len: config.tuning.max_motif_len.max(0) as usize,
                    ..Default::default()
                },
            )
        });
        let masking = config
            .tuning
            .runtime
            .masking()
            .or(owned_masking.as_ref())
            .unwrap();
        let mut table = MaskingTable::new();
        for id in 0..block.seqs().len() {
            masking
                .apply_masking(block.seqs_mut().get_mut(id), bits, id, Some(&mut table))
                .map_err(io::Error::other)?;
        }
        table.remove(block.seqs_mut(), 0, false);
        block.set_soft_masking_table(table);
        let seeds = HashedSeedSet::from_block_with_min_query_len(
            block,
            config.round_query_skip.as_deref(),
            setup.seed_complexity_cut,
            soft_masking,
            &setup.shapes,
            &setup.reduction,
            if config.translated_query_layout.is_some() {
                6
            } else {
                1
            },
            config.min_query_len,
        )
        .map_err(io::Error::other)?;
        Ok(Arc::new(seeds))
    })();
    block.remove_soft_masking(0, false);
    *sequences = block.take_seqs();
    result
}

pub(crate) fn resolve_query(
    sequences: &mut crate::data::sequence_set::SequenceSet,
    mut translated: Option<&mut crate::data::block::Block>,
    config: &super::blastp::BlastpConfig,
    matrix: &crate::stats::score_matrix::ScoreMatrix,
    setup: &crate::search::sensitivity::SetupSearchResult,
    soft_override: crate::masking::MaskingAlgo,
    frequency_masking: bool,
    database_letters: u64,
) -> io::Result<bool> {
    let mut mode = config.tuning.algorithm.selected()?;
    if mode == Algo::Auto {
        if config.tuning.algorithm.auto_candidate(
            sequences.letters() as u64,
            database_letters,
            config.target_indexed,
            config.swipe_all,
            setup.minimizer_window,
            setup.sketch_size,
        )? {
            let seeds = build_query_seeds(
                sequences,
                translated.as_deref_mut(),
                config,
                matrix,
                setup,
                soft_override,
                frequency_masking,
            )?;
            if crate::run::double_indexed::use_query_index(
                seeds.max_table_size(),
                crate::util::system::l3_cache_size(),
            ) {
                mode = Algo::QueryIndexed;
                config.tuning.algorithm.resolve(mode, Some(seeds))?;
            } else {
                mode = Algo::DoubleIndexed;
                config.tuning.algorithm.resolve(mode, None)?;
            }
        } else {
            mode = Algo::DoubleIndexed;
            config.tuning.algorithm.resolve(mode, None)?;
        }
    }
    if mode == Algo::QueryIndexed && config.tuning.algorithm.cached_seeds()?.is_none() {
        let seeds = build_query_seeds(
            sequences,
            translated.as_deref_mut(),
            config,
            matrix,
            setup,
            soft_override,
            frequency_masking,
        )?;
        config.tuning.algorithm.resolve(mode, Some(seeds))?;
    }
    Ok(matches!(mode, Algo::QueryIndexed | Algo::CtgSeed))
}
