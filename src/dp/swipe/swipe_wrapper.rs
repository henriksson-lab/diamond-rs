//! Dispatch and scheduling facade for `dp/swipe/swipe_wrapper.cpp`.
//!
//! The score kernels are architecture-neutral implementations in the parent
//! module. This file mirrors the C++ wrapper's overload selection, worker/task
//! partitioning, and public entry points while retaining those scalar kernels.

use crate::align::hsp::Hsp;
use crate::basic::statistics::StatValue;
use crate::data::sequence_set::SequenceSet;
use crate::{stats, util::geo};

use super::{
    targets, Anchor, CarryOver, DpTarget, Flags, HspValues, Params, TargetVec, Targets, BINS,
    SCORE_BINS,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowCounterKind {
    Dummy,
    Vector,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DispatchCell {
    Score,
    Forward,
    Backward,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IdMaskKind {
    Dummy,
    Vector,
}

/// Runtime representation of the C++ `SwipeConfig<tb, RC, C, IdM>` template.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DispatchConfig {
    pub traceback: bool,
    pub row_counter: RowCounterKind,
    pub cell: DispatchCell,
    pub id_mask: IdMaskKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SwipeRuntimeConfig {
    pub threads: usize,
    pub threads_align: usize,
    pub thread_pool: bool,
    pub swipe_task_size: i64,
    /// Architecture-neutral stand-in for `ScoreTraits<Sv>::CHANNELS`.
    pub channels: usize,
}

impl Default for SwipeRuntimeConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            threads_align: 0,
            thread_pool: false,
            swipe_task_size: i64::MAX,
            channels: 1,
        }
    }
}

/// C++ target-vector `sort` overload.
pub fn sort_targets(targets: &mut [DpTarget], band_bin: i32, col_bin: i32) {
    super::sort(targets, band_bin, col_bin);
}

/// C++ `SequenceSet::ConstIterator` `sort` overload (intentionally a no-op).
pub fn sort_sequence_set(_subjects: &SequenceSet) {}

/// C++ scalar `bin(int)` overload.
pub fn score_bin(score: i32) -> usize {
    super::bin_score(score)
}

/// C++ seven-argument `bin` dispatch entry point.
#[allow(clippy::too_many_arguments)]
pub fn bin(
    values: HspValues,
    query_len: i32,
    score: i32,
    ungapped_score: i32,
    dp_size: i64,
    score_width: usize,
    mismatch_estimate: i32,
    params: &Params<'_>,
) -> usize {
    super::bin(
        values,
        query_len,
        score,
        ungapped_score,
        dp_size,
        score_width,
        mismatch_estimate,
        params.cutoff_score_8bit,
        params.max_swipe_dp,
        params.approx_backtrace,
    )
}

/// C++ `matrix_size<TargetVec>` overload.
pub fn matrix_size_targets(
    query_len: i32,
    targets: &[DpTarget],
    flags: Flags,
    channels: i64,
) -> i64 {
    super::matrix_size(query_len, targets, flags, channels)
}

/// C++ `matrix_size<SequenceSet>` overload always returns zero.
pub fn matrix_size_sequence_set(
    _query_len: i32,
    _subjects: &SequenceSet,
    _flags: Flags,
    _channels: i64,
) -> usize {
    0
}

pub fn reversed(values: HspValues) -> bool {
    super::reversed(values)
}

/// Selects the concrete template configuration used by the C++ overload chain.
pub fn dispatch_config(
    values: HspValues,
    round: i32,
    bin: usize,
) -> Result<DispatchConfig, String> {
    if values == HspValues::NONE {
        return Ok(DispatchConfig {
            traceback: false,
            row_counter: RowCounterKind::Dummy,
            cell: DispatchCell::Score,
            id_mask: IdMaskKind::Dummy,
        });
    }
    if bin < SCORE_BINS {
        return Ok(DispatchConfig {
            traceback: true,
            row_counter: RowCounterKind::Vector,
            cell: DispatchCell::Score,
            id_mask: IdMaskKind::Dummy,
        });
    }
    match round {
        0 if values.any(HspValues::IDENT | HspValues::LENGTH) => Ok(DispatchConfig {
            traceback: false,
            row_counter: RowCounterKind::Vector,
            cell: DispatchCell::Forward,
            id_mask: IdMaskKind::Vector,
        }),
        0 => Ok(DispatchConfig {
            traceback: false,
            row_counter: RowCounterKind::Vector,
            cell: DispatchCell::Score,
            id_mask: IdMaskKind::Dummy,
        }),
        1 if values.any(HspValues::MISMATCHES | HspValues::GAP_OPENINGS) => Ok(DispatchConfig {
            traceback: false,
            row_counter: RowCounterKind::Vector,
            cell: DispatchCell::Backward,
            id_mask: IdMaskKind::Vector,
        }),
        1 => Ok(DispatchConfig {
            traceback: false,
            row_counter: RowCounterKind::Vector,
            cell: DispatchCell::Score,
            id_mask: IdMaskKind::Dummy,
        }),
        _ => Err("Unreachable".to_owned()),
    }
}

/// Architecture-neutral counterpart of the C++ `dispatch_swipe` overload set.
pub fn dispatch_swipe(
    targets: &[DpTarget],
    overflow: &mut TargetVec,
    round: i32,
    bin: usize,
    params: &Params<'_>,
) -> Result<Vec<Hsp>, String> {
    let _configuration = dispatch_config(params.v, round, bin)?;
    Ok(super::dispatch_swipe(targets, overflow, params))
}

/// C++ `swipe_worker`; `worker` and `worker_count` model its shared atomic
/// channel allocator deterministically. Every target is evaluated once.
#[allow(clippy::too_many_arguments)]
pub fn swipe_worker(
    targets: &[DpTarget],
    worker: usize,
    worker_count: usize,
    channels: usize,
    round: i32,
    bin: usize,
    params: &Params<'_>,
) -> Result<(Vec<Hsp>, TargetVec), String> {
    if worker_count == 0 || channels == 0 {
        return Err("SWIPE worker count and channel count must be positive".to_owned());
    }
    let mut output = Vec::new();
    let mut overflow = TargetVec::default();
    let stride = worker_count
        .checked_mul(channels)
        .ok_or_else(|| "SWIPE worker stride overflow".to_owned())?;
    let mut begin = worker
        .checked_mul(channels)
        .ok_or_else(|| "SWIPE worker offset overflow".to_owned())?;
    while begin < targets.len() {
        let end = (begin + channels).min(targets.len());
        output.extend(dispatch_swipe(
            &targets[begin..end],
            &mut overflow,
            round,
            bin,
            params,
        )?);
        begin += stride;
    }
    Ok((output, overflow))
}

/// C++ `swipe_task` for one already-sized task range.
pub fn swipe_task(
    targets: &[DpTarget],
    round: i32,
    bin: usize,
    params: &Params<'_>,
) -> Result<(Vec<Hsp>, TargetVec), String> {
    let mut overflow = TargetVec::default();
    let output = dispatch_swipe(targets, &mut overflow, round, bin, params)?;
    Ok((output, overflow))
}

fn append_result(
    output: &mut Vec<Hsp>,
    overflow: &mut TargetVec,
    mut result: (Vec<Hsp>, TargetVec),
) {
    output.append(&mut result.0);
    overflow.push_back_vec(&result.1);
}

/// C++ `swipe_threads` with explicit replacements for global thread settings
/// and the optional thread-pool pointer. Work is deterministic but preserves
/// worker merge order, task boundaries, counters, and overflow aggregation.
pub fn swipe_threads(
    targets: &[DpTarget],
    overflow: &mut TargetVec,
    round: i32,
    bin: usize,
    params: &Params<'_>,
    runtime: SwipeRuntimeConfig,
) -> Result<Vec<Hsp>, String> {
    if targets.is_empty() {
        return Ok(Vec::new());
    }
    let channels = runtime.channels.max(1);
    if params.flags.any(Flags::PARALLEL) {
        let workers = if runtime.threads_align != 0 {
            runtime.threads_align
        } else {
            runtime.threads
        };
        if workers == 0 {
            return Ok(Vec::new());
        }
        let mut output = Vec::new();
        for worker in 0..workers {
            append_result(
                &mut output,
                overflow,
                swipe_worker(targets, worker, workers, channels, round, bin, params)?,
            );
        }
        return Ok(output);
    }
    if !runtime.thread_pool {
        return dispatch_swipe(targets, overflow, round, bin, params);
    }

    let mut output = Vec::new();
    let mut task_begin = 0usize;
    let mut task_cells = 0i64;
    let task_limit = runtime.swipe_task_size.max(1);
    let mut i = 0usize;
    while i < targets.len() {
        let end = (i + channels).min(targets.len());
        task_cells = task_cells.saturating_add(
            targets[i..end]
                .iter()
                .map(|target| target.cells(params.flags, params.query.len() as i32))
                .sum::<i64>(),
        );
        i = end;
        if task_cells >= task_limit {
            params.inc_stat(StatValue::SwipeTasksTotal, 1);
            params.inc_stat(StatValue::SwipeTasksAsync, 1);
            append_result(
                &mut output,
                overflow,
                swipe_task(&targets[task_begin..i], round, bin, params)?,
            );
            task_begin = i;
            task_cells = 0;
        }
    }
    if task_begin == 0 {
        params.inc_stat(StatValue::SwipeTasksTotal, 1);
        return dispatch_swipe(targets, overflow, round, bin, params);
    }
    if task_begin < targets.len() {
        params.inc_stat(StatValue::SwipeTasksTotal, 1);
        params.inc_stat(StatValue::SwipeTasksAsync, 1);
        append_result(
            &mut output,
            overflow,
            swipe_task(&targets[task_begin..], round, bin, params)?,
        );
    }
    Ok(output)
}

pub fn swipe_bin(
    bin: usize,
    targets: &mut [DpTarget],
    round: i32,
    params: &mut Params<'_>,
) -> (Vec<Hsp>, TargetVec) {
    super::swipe_bin(bin, targets, round, params)
}

pub fn mismatch_est(query_len: i32, target_len: i32, alignment_len: i32, values: HspValues) -> i32 {
    super::mismatch_est(query_len, target_len, alignment_len, values)
}

pub fn recompute_reversed(hsps: &mut [Hsp], params: &mut Params<'_>) -> Vec<Hsp> {
    let mut dp_targets = targets();
    let query_len = params.query.len() as i32;
    for hsp in hsps.iter() {
        let query_cover = hsp.query_cover(params.query_source_len);
        let subject_cover = hsp.subject_cover(hsp.target_seq.len() as i32);
        let query_cutoff = if params.query_or_target_cover > 0.0 {
            params.query_or_target_cover
        } else {
            params.query_cover
        };
        let subject_cutoff = if params.query_or_target_cover > 0.0 {
            params.query_or_target_cover
        } else {
            params.subject_cover
        };
        let (query_min_len, subject_min_len) = hsp.min_range_len(
            query_cutoff,
            subject_cutoff,
            query_len,
            hsp.target_seq.len() as i32,
        );
        let query_approx_id = stats::approx_id(hsp.score, query_min_len, 0);
        let subject_approx_id = stats::approx_id(hsp.score, subject_min_len, 0);
        if query_cover < params.query_cover
            || subject_cover < params.subject_cover
            || query_cover.max(subject_cover) < params.query_or_target_cover
            || (params.query_or_target_cover == 0.0
                && query_approx_id.min(subject_approx_id) < params.approx_min_id)
            || (params.query_or_target_cover > 0.0
                && query_approx_id.max(subject_approx_id) < params.approx_min_id)
        {
            continue;
        }
        let target_len = hsp.subject_range.end;
        if target_len <= 0 || hsp.swipe_bin < 0 {
            continue;
        }
        let reversed_target = hsp.target_seq[..target_len as usize].to_vec();
        let band = if params.flags.any(Flags::FULL_MATRIX) {
            query_len
        } else {
            hsp.d_end - hsp.d_begin
        };
        let selected_bin = bin(
            params.v,
            band,
            hsp.score,
            0,
            i64::MAX,
            0,
            mismatch_est(hsp.query_range.end, target_len, hsp.length, params.v),
            params,
        )
        .max(hsp.swipe_bin as usize);
        debug_assert!(selected_bin >= SCORE_BINS);
        if selected_bin >= BINS {
            continue;
        }
        let carry_over = CarryOver::new(
            hsp.query_range.end,
            hsp.subject_range.end,
            hsp.identities,
            hsp.length,
        );
        let mut target = DpTarget::new(
            reversed_target,
            hsp.target_seq.len() as i32,
            geo::rev_diag(hsp.d_end - 1, query_len, target_len),
            geo::rev_diag(hsp.d_begin, query_len, target_len) + 1,
            hsp.swipe_target as i64,
            query_len,
            carry_over,
            Anchor::default(),
        );
        if let Some(matrix) = &hsp.matrix {
            target = target.with_matrix(matrix.clone(), params.cbs_matrix_scale);
        }
        dp_targets[selected_bin].push_back(target);
    }

    let mut reversed_query = params.query.to_vec();
    reversed_query.reverse();
    let reversed_bias;
    let composition_bias = match params.composition_bias {
        Some(bias) => {
            reversed_bias = reverse_composition_bias(bias, params.query.len());
            Some(reversed_bias.as_slice())
        }
        None => None,
    };
    let mut reversed_params = Params {
        query: &reversed_query,
        composition_bias,
        reverse_targets: true,
        target_max_len: 0,
        ..params.clone()
    };
    let mut output = Vec::new();
    let mut overflow_targets: Option<Targets> = None;
    for selected_bin in SCORE_BINS..BINS {
        reversed_params.target_max_len = dp_targets[selected_bin].max_len();
        reversed_params.swipe_bin = selected_bin as i32;
        let (mut bin_output, overflow) = swipe_bin(
            selected_bin,
            dp_targets[selected_bin].as_mut_slice(),
            1,
            &mut reversed_params,
        );
        if !overflow.empty() && selected_bin + 1 < BINS {
            let overflow_bins = overflow_targets.get_or_insert_with(targets);
            for target in overflow.as_slice() {
                for hsp in hsps.iter() {
                    if hsp.swipe_target as i64 == target.target_idx {
                        let mut overflow_target = DpTarget::new(
                            hsp.target_seq.clone(),
                            hsp.target_seq.len() as i32,
                            hsp.d_begin,
                            hsp.d_end,
                            hsp.swipe_target as i64,
                            params.query.len() as i32,
                            CarryOver::default(),
                            Anchor::default(),
                        );
                        if let Some(matrix) = &hsp.matrix {
                            overflow_target = overflow_target
                                .with_matrix(matrix.clone(), params.cbs_matrix_scale);
                        }
                        overflow_bins[selected_bin + 1].push_back(overflow_target);
                    }
                }
            }
        }
        output.append(&mut bin_output);
    }
    if let Some(overflow_targets) = overflow_targets {
        output.extend(swipe(&overflow_targets, params));
    }
    output
}

pub fn swipe(targets: &Targets, params: &mut Params<'_>) -> Vec<Hsp> {
    let mut previous = (Vec::new(), TargetVec::default());
    let mut forward = Vec::new();
    let mut reverse = Vec::new();
    for algorithm_bin in 0..super::ALGO_BINS {
        for score_bin in 0..SCORE_BINS {
            let selected_bin = algorithm_bin * SCORE_BINS + score_bin;
            let mut round_targets = TargetVec::default();
            round_targets.reserve(targets[selected_bin].size() + previous.1.size());
            round_targets.push_back_vec(&targets[selected_bin]);
            round_targets.push_back_vec(&previous.1);
            params.target_max_len = round_targets.max_len();
            params.swipe_bin = selected_bin as i32;
            previous = swipe_bin(selected_bin, round_targets.as_mut_slice(), 0, params);
            if algorithm_bin == 0 {
                forward.append(&mut previous.0);
            } else {
                reverse.append(&mut previous.0);
            }
        }
        debug_assert!(previous.1.empty());
    }
    if !reverse.is_empty() {
        forward.append(&mut recompute_reversed(&mut reverse, params));
    }
    forward
}

pub fn swipe_set(subjects: &SequenceSet, params: &mut Params<'_>) -> Vec<Hsp> {
    let selected_bin = bin(params.v, 0, 0, 0, 0, 0, 0, params);
    let mut round_targets = TargetVec::default();
    round_targets.reserve(subjects.len() as i64);
    for target_id in 0..subjects.len() {
        let sequence = subjects.get(target_id).to_vec();
        let query_len = params.query.len() as i32;
        let target_len = sequence.len() as i32;
        round_targets.push_back(DpTarget::new(
            sequence,
            target_len,
            -(target_len - 1),
            query_len,
            target_id as i64,
            query_len,
            CarryOver::default(),
            Anchor::default(),
        ));
    }
    let (mut output, overflow) = swipe_bin(
        selected_bin.min(BINS - 1),
        round_targets.as_mut_slice(),
        0,
        params,
    );
    if reversed(params.v) {
        output = recompute_reversed(&mut output, params);
    }
    if selected_bin < BINS - 1 && !overflow.empty() {
        let mut overflow_bins = targets();
        overflow_bins[selected_bin + 1] = overflow;
        output.extend(swipe(&overflow_bins, params));
    }
    output
}

fn reverse_composition_bias(bias: &[i8], query_len: usize) -> Vec<i8> {
    let mut output: Vec<i8> = bias[..query_len.min(bias.len())].to_vec();
    output.reverse();
    output.extend(std::iter::repeat_n(0, 32));
    output
}

#[cfg(test)]
mod tests {
    use std::sync::{Arc, Mutex};

    use crate::basic::statistics::Statistics;
    use crate::stats::score_matrix::ScoreMatrix;

    use super::*;
    use crate::dp::swipe::{targets, CarryOver};

    fn matrix() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    fn target(id: i64, len: usize) -> DpTarget {
        DpTarget::full(vec![0; len], len as i32, id, CarryOver::default())
    }

    #[test]
    fn dispatch_config_matches_all_cpp_template_branches() {
        assert_eq!(
            dispatch_config(HspValues::NONE, 0, 0).unwrap().row_counter,
            RowCounterKind::Dummy
        );
        assert!(dispatch_config(HspValues::COORDS, 0, 0).unwrap().traceback);
        assert_eq!(
            dispatch_config(HspValues::IDENT, 0, 3).unwrap().cell,
            DispatchCell::Forward
        );
        assert_eq!(
            dispatch_config(HspValues::MISMATCHES, 1, 3).unwrap().cell,
            DispatchCell::Backward
        );
        assert_eq!(
            dispatch_config(HspValues::QUERY_END, 1, 3).unwrap().cell,
            DispatchCell::Score
        );
        assert_eq!(
            dispatch_config(HspValues::COORDS, 2, 3).unwrap_err(),
            "Unreachable"
        );
    }

    #[test]
    fn overloads_preserve_sort_and_matrix_size_behavior() {
        let mut values = vec![
            DpTarget::new(
                vec![0; 8],
                8,
                -2,
                4,
                0,
                6,
                CarryOver::default(),
                Default::default(),
            ),
            DpTarget::new(
                vec![0; 3],
                3,
                0,
                2,
                1,
                6,
                CarryOver::default(),
                Default::default(),
            ),
        ];
        sort_targets(&mut values, 2, 2);
        assert_eq!(values[0].target_idx, 1);
        assert_eq!(matrix_size_targets(6, &values, Flags::FULL_MATRIX, 8), 192);
        assert_eq!(
            matrix_size_sequence_set(6, &SequenceSet::new(), Flags::NONE, 8),
            0
        );
        assert_eq!(score_bin(254), 0);
        assert_eq!(score_bin(255), 1);
        assert_eq!(score_bin(65535), 2);
    }

    #[test]
    fn parallel_workers_cover_targets_once_and_merge_by_worker() {
        let sm = matrix();
        let query = vec![0; 4];
        let mut params = Params::new(&query, &sm);
        params.flags = Flags::FULL_MATRIX | Flags::PARALLEL;
        params.v = HspValues::NONE;
        let input = vec![target(0, 4), target(1, 4), target(2, 4), target(3, 4)];
        let mut overflow = TargetVec::default();
        let output = swipe_threads(
            &input,
            &mut overflow,
            0,
            0,
            &params,
            SwipeRuntimeConfig {
                threads: 2,
                channels: 1,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            output.iter().map(|h| h.swipe_target).collect::<Vec<_>>(),
            vec![0, 2, 1, 3]
        );
        assert!(overflow.empty());
    }

    #[test]
    fn task_pool_batches_by_cells_and_updates_cpp_counters() {
        let sm = matrix();
        let query = vec![0; 4];
        let statistics = Arc::new(Mutex::new(Statistics::new()));
        let mut params = Params::new(&query, &sm);
        params.flags = Flags::FULL_MATRIX;
        params.v = HspValues::NONE;
        params.statistics = Some(statistics.clone());
        let input = vec![target(0, 4), target(1, 4), target(2, 4)];
        let mut overflow = TargetVec::default();
        let output = swipe_threads(
            &input,
            &mut overflow,
            0,
            0,
            &params,
            SwipeRuntimeConfig {
                thread_pool: true,
                swipe_task_size: 20,
                channels: 1,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(output.len(), 3);
        let statistics = statistics.lock().unwrap();
        assert_eq!(statistics.get(StatValue::SwipeTasksTotal), 2);
        assert_eq!(statistics.get(StatValue::SwipeTasksAsync), 2);
    }

    #[test]
    fn public_swipe_facades_execute_parent_kernels() {
        let sm = matrix();
        let query = vec![0; 4];
        let mut params = Params::new(&query, &sm);
        params.flags = Flags::FULL_MATRIX;
        params.v = HspValues::NONE;
        let mut bins = targets();
        bins[0].push_back(target(7, 4));
        let result = swipe(&bins, &mut params);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].swipe_target, 7);
    }
}
