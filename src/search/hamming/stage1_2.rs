//! Translation of `diamond/src/search/hamming/stage1_2.cpp`.
//!
//! The C++ overloads consume a `JoinIterator` and a large search `WorkSet`.
//! Rust represents each joined key as a pair of borrowed slices and keeps the
//! stage-2 tile callback explicit, avoiding a dependency on global config or
//! the monolithic C++ ownership graph.

use std::fmt;
use std::marker::PhantomData;

use crate::basic::packed_loc::PackedLoc;
use crate::data::sequence_set::SequenceSet;
use crate::search::kmer_ranking::{KmerRanking, PackedLocId};

use super::{
    stage1, stage1_longest_combo_lin, stage1_mutual_cov, stage1_mutual_cov_query_lin,
    stage1_mutual_cov_target_lin, stage1_query_lin, stage1_query_lin_ranked, stage1_self,
    stage1_self_mutual_cov, stage1_target_lin, HitField,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage1KernelKind {
    Stage1,
    Stage1Self,
    Stage1QueryLin,
    Stage1QueryLinRanked,
    Stage1TargetLin,
    Stage1LongestComboLin,
    Stage1MutualCov,
    Stage1SelfMutualCov,
    Stage1MutualCovQueryLin,
    Stage1MutualCovTargetLin,
}

/// Explicit replacement for the fields read from C++ `config` and
/// `Search::Config` by this translation unit.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Stage1DispatchConfig {
    pub lin_stage1_combo: bool,
    pub lin_stage1_query: bool,
    pub lin_stage1_target: bool,
    pub min_length_ratio: f64,
    pub self_search: bool,
    pub current_ref_block: u32,
    pub global_ranking_targets: bool,
    /// Models builds that define C++ `HIT_KEEP_TARGET_ID`.
    pub hit_keep_target_id: bool,
}

/// C++ `stage1_dispatch(const Search::Config*, PackedLocId)`.
pub fn stage1_dispatch_packed_loc_id(cfg: &Stage1DispatchConfig) -> Stage1KernelKind {
    if cfg.lin_stage1_combo {
        return Stage1KernelKind::Stage1LongestComboLin;
    }
    if cfg.lin_stage1_query {
        return if cfg.min_length_ratio > 0.0 {
            Stage1KernelKind::Stage1MutualCovQueryLin
        } else {
            Stage1KernelKind::Stage1QueryLinRanked
        };
    }
    if cfg.lin_stage1_target {
        return if cfg.min_length_ratio > 0.0 {
            Stage1KernelKind::Stage1MutualCovTargetLin
        } else {
            Stage1KernelKind::Stage1TargetLin
        };
    }
    if cfg.min_length_ratio > 0.0 {
        return if cfg.self_search && cfg.current_ref_block == 0 {
            Stage1KernelKind::Stage1SelfMutualCov
        } else {
            Stage1KernelKind::Stage1MutualCov
        };
    }
    if cfg.self_search && cfg.current_ref_block == 0 {
        Stage1KernelKind::Stage1Self
    } else {
        Stage1KernelKind::Stage1
    }
}

/// C++ `stage1_dispatch(const Search::Config*, PackedLoc)`.
pub fn stage1_dispatch_packed_loc(cfg: &Stage1DispatchConfig) -> Stage1KernelKind {
    if cfg.lin_stage1_query {
        Stage1KernelKind::Stage1QueryLin
    } else if cfg.lin_stage1_target {
        Stage1KernelKind::Stage1TargetLin
    } else if cfg.self_search && cfg.current_ref_block == 0 {
        Stage1KernelKind::Stage1Self
    } else {
        Stage1KernelKind::Stage1
    }
}

/// C++ `keep_target_id(const Search::Config&)`, including the compile-time
/// `HIT_KEEP_TARGET_ID` branch as explicit configuration.
pub fn keep_target_id(cfg: &Stage1DispatchConfig) -> bool {
    cfg.hit_keep_target_id
        || cfg.min_length_ratio != 0.0
        || cfg.global_ranking_targets
        || (cfg.self_search && cfg.current_ref_block == 0)
        || cfg.lin_stage1_combo
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Stage1Statistics {
    /// C++ `Statistics::SEEDS_HIT`: joined seed keys processed.
    pub seeds_hit: u64,
    /// C++ `Statistics::SEED_HITS`: candidate location pairs considered by
    /// the selected kernel.
    pub seed_hits: u64,
}

pub struct WorkSet<'a, L, F> {
    pub query_sequences: &'a SequenceSet,
    pub target_sequences: &'a SequenceSet,
    pub tile_size: usize,
    pub hamming_filter_id: u32,
    pub kmer_ranking: Option<&'a KmerRanking>,
    pub stats: Stage1Statistics,
    pub search_tile: F,
    marker: PhantomData<fn(L)>,
}

impl<'a, L, F> WorkSet<'a, L, F> {
    pub fn new(
        query_sequences: &'a SequenceSet,
        target_sequences: &'a SequenceSet,
        tile_size: usize,
        hamming_filter_id: u32,
        kmer_ranking: Option<&'a KmerRanking>,
        search_tile: F,
    ) -> Self {
        Self {
            query_sequences,
            target_sequences,
            tile_size,
            hamming_filter_id,
            kmer_ranking,
            stats: Stage1Statistics::default(),
            search_tile,
            marker: PhantomData,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage1RunError {
    MissingKmerRanking,
    UnsupportedKernel(Stage1KernelKind),
}

impl fmt::Display for Stage1RunError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingKmerRanking => formatter.write_str("stage1 kernel requires k-mer ranking"),
            Self::UnsupportedKernel(kernel) => {
                write!(
                    formatter,
                    "stage1 kernel {kernel:?} is invalid for location type"
                )
            }
        }
    }
}

impl std::error::Error for Stage1RunError {}

/// C++ `run_stage1(JoinIterator<PackedLoc>&, ...)`.
pub fn run_stage1_packed_loc<'group, I, F>(
    groups: I,
    work_set: &mut WorkSet<'_, PackedLoc, F>,
    cfg: &Stage1DispatchConfig,
) -> Result<(), Stage1RunError>
where
    I: IntoIterator<Item = (&'group [PackedLoc], &'group [PackedLoc])>,
    F: FnMut(&mut HitField, usize, usize, &[PackedLoc], &[PackedLoc]),
{
    let kernel = stage1_dispatch_packed_loc(cfg);
    for (query, target) in groups {
        work_set.stats.seeds_hit = work_set.stats.seeds_hit.wrapping_add(1);
        let seed_hits = execute_packed_loc(kernel, query, target, work_set)?;
        work_set.stats.seed_hits = work_set.stats.seed_hits.wrapping_add(seed_hits);
    }
    Ok(())
}

fn execute_packed_loc<F>(
    kernel: Stage1KernelKind,
    query: &[PackedLoc],
    target: &[PackedLoc],
    work_set: &mut WorkSet<'_, PackedLoc, F>,
) -> Result<u64, Stage1RunError>
where
    F: FnMut(&mut HitField, usize, usize, &[PackedLoc], &[PackedLoc]),
{
    let WorkSet {
        query_sequences,
        target_sequences,
        tile_size,
        hamming_filter_id,
        search_tile,
        ..
    } = work_set;
    let count = match kernel {
        Stage1KernelKind::Stage1 => stage1(
            query,
            target,
            query_sequences,
            target_sequences,
            *tile_size,
            *hamming_filter_id,
            search_tile,
        ),
        Stage1KernelKind::Stage1Self => stage1_self(
            target,
            target_sequences,
            *tile_size,
            *hamming_filter_id,
            search_tile,
        ),
        Stage1KernelKind::Stage1QueryLin => stage1_query_lin(
            query,
            target,
            query_sequences,
            target_sequences,
            *tile_size,
            *hamming_filter_id,
            search_tile,
        ),
        Stage1KernelKind::Stage1TargetLin => stage1_target_lin(
            query,
            target,
            query_sequences,
            target_sequences,
            *tile_size,
            *hamming_filter_id,
            search_tile,
        ),
        invalid => return Err(Stage1RunError::UnsupportedKernel(invalid)),
    };
    Ok(count)
}

/// C++ `run_stage1(JoinIterator<PackedLocId>&, ...)`.
pub fn run_stage1_packed_loc_id<'group, I, F>(
    groups: I,
    work_set: &mut WorkSet<'_, PackedLocId, F>,
    cfg: &Stage1DispatchConfig,
) -> Result<(), Stage1RunError>
where
    I: IntoIterator<Item = (&'group [PackedLocId], &'group [PackedLocId])>,
    F: FnMut(&mut HitField, usize, usize, &[PackedLocId], &[PackedLocId]),
{
    let kernel = stage1_dispatch_packed_loc_id(cfg);
    for (query, target) in groups {
        work_set.stats.seeds_hit = work_set.stats.seeds_hit.wrapping_add(1);
        let seed_hits = execute_packed_loc_id(kernel, query, target, work_set, cfg)?;
        work_set.stats.seed_hits = work_set.stats.seed_hits.wrapping_add(seed_hits);
    }
    Ok(())
}

fn execute_packed_loc_id<F>(
    kernel: Stage1KernelKind,
    query: &[PackedLocId],
    target: &[PackedLocId],
    work_set: &mut WorkSet<'_, PackedLocId, F>,
    cfg: &Stage1DispatchConfig,
) -> Result<u64, Stage1RunError>
where
    F: FnMut(&mut HitField, usize, usize, &[PackedLocId], &[PackedLocId]),
{
    let WorkSet {
        query_sequences,
        target_sequences,
        tile_size,
        hamming_filter_id,
        kmer_ranking,
        search_tile,
        ..
    } = work_set;
    let count = match kernel {
        Stage1KernelKind::Stage1 => stage1(
            query,
            target,
            query_sequences,
            target_sequences,
            *tile_size,
            *hamming_filter_id,
            search_tile,
        ),
        Stage1KernelKind::Stage1Self => stage1_self(
            target,
            target_sequences,
            *tile_size,
            *hamming_filter_id,
            search_tile,
        ),
        Stage1KernelKind::Stage1QueryLinRanked => stage1_query_lin_ranked(
            query,
            target,
            query_sequences,
            target_sequences,
            (*kmer_ranking).ok_or(Stage1RunError::MissingKmerRanking)?,
            *tile_size,
            *hamming_filter_id,
            search_tile,
        ),
        Stage1KernelKind::Stage1TargetLin => stage1_target_lin(
            query,
            target,
            query_sequences,
            target_sequences,
            *tile_size,
            *hamming_filter_id,
            search_tile,
        ),
        Stage1KernelKind::Stage1LongestComboLin => stage1_longest_combo_lin(
            query,
            target,
            query_sequences,
            target_sequences,
            *tile_size,
            *hamming_filter_id,
            search_tile,
        ),
        Stage1KernelKind::Stage1MutualCov => stage1_mutual_cov(
            query,
            target,
            query_sequences,
            target_sequences,
            *tile_size,
            *hamming_filter_id,
            cfg.min_length_ratio,
            search_tile,
        ),
        Stage1KernelKind::Stage1SelfMutualCov => stage1_self_mutual_cov(
            target,
            target_sequences,
            *hamming_filter_id,
            cfg.min_length_ratio,
            search_tile,
        ),
        Stage1KernelKind::Stage1MutualCovQueryLin => stage1_mutual_cov_query_lin(
            query,
            target,
            query_sequences,
            target_sequences,
            *hamming_filter_id,
            cfg.min_length_ratio,
            cfg.self_search,
            search_tile,
        ),
        Stage1KernelKind::Stage1MutualCovTargetLin => stage1_mutual_cov_target_lin(
            query,
            target,
            query_sequences,
            target_sequences,
            *hamming_filter_id,
            cfg.min_length_ratio,
            search_tile,
        ),
        invalid => return Err(Stage1RunError::UnsupportedKernel(invalid)),
    };
    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::Letter;

    fn sequences(lengths: &[usize]) -> SequenceSet {
        let mut sequences = SequenceSet::new();
        for &length in lengths {
            let sequence: Vec<Letter> = (0..length).map(|index| (index % 20) as Letter).collect();
            sequences.push(&sequence);
        }
        sequences
    }

    #[test]
    fn packed_loc_runner_processes_each_join_group_and_accumulates_stats() {
        let query_sequences = sequences(&[96]);
        let target_sequences = sequences(&[96]);
        let q1 = [
            PackedLoc::new(query_sequences.position(0, 20) as u64),
            PackedLoc::new(query_sequences.position(0, 30) as u64),
        ];
        let s1 = [PackedLoc::new(target_sequences.position(0, 20) as u64)];
        let q2 = [PackedLoc::new(query_sequences.position(0, 40) as u64)];
        let s2 = [
            PackedLoc::new(target_sequences.position(0, 30) as u64),
            PackedLoc::new(target_sequences.position(0, 40) as u64),
        ];
        let mut tile_calls = Vec::new();
        let mut work_set = WorkSet::new(
            &query_sequences,
            &target_sequences,
            2,
            48,
            None,
            |hits: &mut HitField, i, j, _: &[PackedLoc], _: &[PackedLoc]| {
                tile_calls.push((i, j, hits.query_count()));
            },
        );

        run_stage1_packed_loc(
            [(&q1[..], &s1[..]), (&q2[..], &s2[..])],
            &mut work_set,
            &Stage1DispatchConfig::default(),
        )
        .unwrap();

        assert_eq!(work_set.stats.seeds_hit, 2);
        assert_eq!(work_set.stats.seed_hits, 4);
        assert_eq!(tile_calls, vec![(0, 0, 2), (0, 0, 1)]);
    }

    #[test]
    fn self_dispatch_uses_target_side_and_triangular_candidate_count() {
        let query_sequences = sequences(&[96]);
        let target_sequences = sequences(&[96]);
        let query = [PackedLoc::new(query_sequences.position(0, 20) as u64)];
        let target = [
            PackedLoc::new(target_sequences.position(0, 20) as u64),
            PackedLoc::new(target_sequences.position(0, 30) as u64),
            PackedLoc::new(target_sequences.position(0, 40) as u64),
        ];
        let mut work_set = WorkSet::new(
            &query_sequences,
            &target_sequences,
            2,
            48,
            None,
            |_: &mut HitField, _, _, _: &[PackedLoc], _: &[PackedLoc]| {},
        );
        run_stage1_packed_loc(
            [(&query[..], &target[..])],
            &mut work_set,
            &Stage1DispatchConfig {
                self_search: true,
                ..Stage1DispatchConfig::default()
            },
        )
        .unwrap();
        assert_eq!(
            work_set.stats,
            Stage1Statistics {
                seeds_hit: 1,
                seed_hits: 3
            }
        );
    }

    #[test]
    fn packed_loc_id_runner_executes_ranked_dispatch_and_requires_ranking() {
        let query_sequences = sequences(&[80, 100]);
        let target_sequences = sequences(&[90]);
        let query = [
            PackedLocId::new(query_sequences.position(0, 20) as u64, 0),
            PackedLocId::new(query_sequences.position(1, 20) as u64, 1),
        ];
        let target = [PackedLocId::new(target_sequences.position(0, 20) as u64, 0)];
        let cfg = Stage1DispatchConfig {
            lin_stage1_query: true,
            ..Stage1DispatchConfig::default()
        };

        let mut missing = WorkSet::new(
            &query_sequences,
            &target_sequences,
            4,
            48,
            None,
            |_: &mut HitField, _, _, _: &[PackedLocId], _: &[PackedLocId]| {},
        );
        assert_eq!(
            run_stage1_packed_loc_id([(&query[..], &target[..])], &mut missing, &cfg),
            Err(Stage1RunError::MissingKmerRanking)
        );
        assert_eq!(missing.stats.seeds_hit, 1);

        let ranking = KmerRanking::from_queries(&query_sequences);
        let mut calls = Vec::new();
        let mut ranked = WorkSet::new(
            &query_sequences,
            &target_sequences,
            4,
            48,
            Some(&ranking),
            |_: &mut HitField, i, j, _: &[PackedLocId], _: &[PackedLocId]| calls.push((i, j)),
        );
        run_stage1_packed_loc_id([(&query[..], &target[..])], &mut ranked, &cfg).unwrap();
        assert_eq!(
            ranked.stats,
            Stage1Statistics {
                seeds_hit: 1,
                seed_hits: 1
            }
        );
        assert_eq!(calls, vec![(1, 0)]);
    }

    #[test]
    fn dispatch_precedence_and_target_id_conditions_match_cpp() {
        let cfg = Stage1DispatchConfig {
            lin_stage1_combo: true,
            lin_stage1_query: true,
            lin_stage1_target: true,
            min_length_ratio: 0.8,
            self_search: true,
            current_ref_block: 0,
            global_ranking_targets: true,
            hit_keep_target_id: false,
        };
        assert_eq!(
            stage1_dispatch_packed_loc_id(&cfg),
            Stage1KernelKind::Stage1LongestComboLin
        );
        assert_eq!(
            stage1_dispatch_packed_loc(&cfg),
            Stage1KernelKind::Stage1QueryLin
        );
        assert!(keep_target_id(&cfg));
        assert!(!keep_target_id(&Stage1DispatchConfig::default()));
    }
}
