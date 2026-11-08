//! Cluster-member realignment from `diamond/src/cluster/output.cpp`.

use std::io::Write;
use std::sync::{Arc, Mutex};
#[path = "context_stream.rs"]
mod context_stream;
pub use context_stream::ClusterOutputIo;

use rayon::prelude::*;

use crate::align::hsp::{ByteHspContext, HspContext};
use crate::basic::statistics::Statistics;
use crate::basic::value::{CentroidId, Letter, OId};
use crate::data::block::Block;
use crate::dp::swipe::{
    bin as swipe_bin, swipe, targets as make_dp_targets, CarryOver, DpTarget, Flags, HspValues,
    Params, Targets,
};
use crate::stats::cbs::hauser_correction_window;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::data_structures::{make_flat_array, FlatArray};

pub trait ClusterRealignDatabase: Send {
    fn sequence_count(&self) -> OId;
    fn letters(&self) -> u64;
    fn titles_lazy(&self) -> bool;
    fn set_load_flags(&mut self, sequences: bool, titles: bool);
    fn set_seqinfo_ptr(&mut self, oid: OId) -> Result<(), String>;
    fn tell_seq(&self) -> OId;
    fn load_seqs(&mut self, block_size: u64) -> Result<Arc<Block>, String>;
    fn seqid(&mut self, oid: OId) -> Result<String, String>;
    fn seqid_bytes(&mut self, oid: OId) -> Result<Vec<u8>, String> {
        self.seqid(oid).map(String::into_bytes)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClusterOutputConfig {
    pub block_size: u64,
    pub db_size: Option<u64>,
    pub comp_based_stats_hauser: bool,
    pub cbs_window: i32,
    pub cutoff_score_8bit: i32,
    pub max_swipe_dp: i64,
    pub band_bin: i32,
    pub col_bin: i32,
    pub swipe_task_size: i64,
    pub transcript_len_estimate: f64,
    pub threads_align: u32,
    pub target_parallel_verbosity: u32,
    pub approx_backtrace: bool,
    pub max_evalue: f64,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_or_target_cover: f64,
    pub approx_min_id: f64,
    pub cbs_matrix_scale: i32,
}

impl Default for ClusterOutputConfig {
    fn default() -> Self {
        Self {
            block_size: 1 << 29,
            db_size: None,
            comp_based_stats_hauser: false,
            cbs_window: 40,
            cutoff_score_8bit: 240,
            max_swipe_dp: 1_000_000,
            band_bin: 24,
            col_bin: 400,
            swipe_task_size: 100_000_000,
            transcript_len_estimate: 1.0,
            threads_align: 0,
            target_parallel_verbosity: u32::MAX,
            approx_backtrace: false,
            max_evalue: f64::MAX,
            query_cover: 0.0,
            subject_cover: 0.0,
            query_or_target_cover: 0.0,
            approx_min_id: 0.0,
            cbs_matrix_scale: 1,
        }
    }
}

struct Cfg<'a, D> {
    hsp_values: HspValues,
    lazy_titles: bool,
    clusters: &'a FlatArray<OId>,
    centroids: &'a [OId],
    database: Mutex<&'a mut D>,
    centroid_block: Arc<Block>,
    member_block: Arc<Block>,
}

impl<'a, D> Cfg<'a, D> {
    fn new(
        hsp_values: HspValues,
        lazy_titles: bool,
        clusters: &'a FlatArray<OId>,
        centroids: &'a [OId],
        database: &'a mut D,
        centroid_block: Arc<Block>,
        member_block: Arc<Block>,
    ) -> Self {
        Self {
            hsp_values,
            lazy_titles,
            clusters,
            centroids,
            database: Mutex::new(database),
            centroid_block,
            member_block,
        }
    }
}

/// C++ file-local `align_centroid`.
fn align_centroid<D: ClusterRealignDatabase>(
    centroid: CentroidId,
    statistics: &mut Statistics,
    config: &ClusterOutputConfig,
    score_matrix: &ScoreMatrix,
    cfg: &Cfg<'_, D>,
) -> Result<Vec<ByteHspContext>, String> {
    let centroid_oid = cfg.centroids[centroid as usize];
    let centroid_id = cfg.centroid_block.oid2block_id(centroid_oid)?;
    let centroid_seq = cfg.centroid_block.seqs().get(centroid_id as usize);
    let members = cfg.clusters.range(centroid as u64);
    let member_begin = members.partition_point(|oid| *oid < cfg.member_block.oid_begin());
    let member_end = members.partition_point(|oid| *oid < cfg.member_block.oid_end());
    let mut dp_targets: Targets = make_dp_targets();
    for member_oid in &members[member_begin..member_end] {
        let block_id = cfg.member_block.oid2block_id(*member_oid)?;
        let sequence = cfg.member_block.seqs().get(block_id as usize);
        let bin = swipe_bin(
            cfg.hsp_values,
            centroid_seq.len() as i32,
            0,
            0,
            sequence.len() as i64 * centroid_seq.len() as i64,
            0,
            0,
            config.cutoff_score_8bit,
            config.max_swipe_dp,
            config.approx_backtrace,
        );
        dp_targets[bin].push_back(DpTarget::full(
            sequence.to_vec(),
            sequence.len() as i32,
            block_id as i64,
            CarryOver::default(),
        ));
    }

    let centroid_title = if cfg.lazy_titles {
        cfg.database.lock().unwrap().seqid_bytes(centroid_oid)?
    } else {
        cfg.centroid_block.ids()?.get(centroid_id as usize).to_vec()
    };
    let correction = if config.comp_based_stats_hauser {
        Some(hauser_correction_window(
            centroid_seq,
            score_matrix,
            config.cbs_window as u32 as usize,
        ))
    } else {
        None
    };
    let mut params = Params::new(centroid_seq, score_matrix);
    params.query_id_bytes = Some(&centroid_title);
    params.frame = 0;
    params.query_source_len = centroid_seq.len() as i32;
    params.composition_bias = correction.as_deref();
    params.flags = Flags::FULL_MATRIX;
    params.v = cfg.hsp_values;
    params.cutoff_score_8bit = config.cutoff_score_8bit;
    params.max_swipe_dp = config.max_swipe_dp;
    params.band_bin = config.band_bin;
    params.col_bin = config.col_bin;
    params.swipe_task_size = config.swipe_task_size;
    params.transcript_len_estimate = config.transcript_len_estimate;
    params.threads_align = config.threads_align;
    params.target_parallel_verbosity = config.target_parallel_verbosity;
    let swipe_error = crate::dp::swipe::SwipeError::default();
    params.error_sink = Some(swipe_error.clone());
    params.approx_backtrace = config.approx_backtrace;
    params.max_evalue = config.max_evalue;
    params.query_cover = config.query_cover;
    params.subject_cover = config.subject_cover;
    params.query_or_target_cover = config.query_or_target_cover;
    params.approx_min_id = config.approx_min_id;
    params.cbs_matrix_scale = config.cbs_matrix_scale;
    let swipe_statistics = Arc::new(Mutex::new(Statistics::new()));
    params.statistics = Some(swipe_statistics.clone());

    let mut output = Vec::new();
    let hsps = swipe(&dp_targets, &mut params);
    if let Some(error) = swipe_error.take() {
        return Err(error);
    }
    *statistics += &swipe_statistics.lock().unwrap();
    for mut hsp in hsps {
        let subject_block_id = hsp.swipe_target as usize;
        let subject_oid = cfg.member_block.block_id2oid(subject_block_id as u32);
        let subject_title = if cfg.lazy_titles {
            cfg.database.lock().unwrap().seqid_bytes(subject_oid)?
        } else {
            cfg.member_block.ids()?.get(subject_block_id).to_vec()
        };
        // C++ serializes each context to its intermediate stream and invokes
        // the callback on a fresh deserialized context. Query sequences and
        // corrected scores are absent from that record. In particular the
        // normalized-identity callback observes an empty query sequence.
        hsp.corrected_bit_score = 0.0;
        let mut context = ByteHspContext::new_bytes(
            hsp,
            centroid_id,
            centroid_oid,
            vec![Vec::new()],
            0,
            centroid_title.clone(),
            subject_oid,
            cfg.member_block.seqs().length(subject_block_id),
            subject_title,
            0,
            0,
            Vec::<Letter>::new(),
            0.0,
            0.0,
        );
        context.query_len = centroid_seq.len() as i32;
        output.push(context);
    }
    Ok(output)
}

/// C++ file-local `run_block_pair`; deterministic ownership replaces its
/// reorder queue and temporary serialized file.
fn run_block_pair<D: ClusterRealignDatabase>(
    begin: CentroidId,
    end: CentroidId,
    config: &ClusterOutputConfig,
    score_matrix: &ScoreMatrix,
    cfg: &Cfg<'_, D>,
    io_config: &ClusterOutputIo,
    errors: Option<&crate::util::io::IoErrorOwner>,
) -> Result<(context_stream::ContextFile, Statistics), String> {
    let (file, mut writer) =
        context_stream::capture_io(context_stream::ContextFile::create(io_config), errors)?;
    let statistics = Mutex::new(Statistics::new());
    let sink = crate::output::live_sink::LiveSink::new(
        begin as usize,
        false,
        |contexts: Vec<ByteHspContext>| {
            for context in contexts {
                context_stream::write_context(&mut writer, &context)?;
            }
            Ok(())
        },
    );
    let titles = vec![""; (end - begin) as usize];
    context_stream::capture_io(
        sink.execute(begin as usize, &titles, || {
            (begin..end).into_par_iter().try_for_each(|centroid| {
                let mut stats = Statistics::new();
                let contexts = align_centroid(centroid, &mut stats, config, score_matrix, cfg)
                    .map_err(std::io::Error::other)?;
                *statistics.lock().unwrap() += &stats;
                let bytes = contexts.len() * std::mem::size_of::<ByteHspContext>();
                sink.push(centroid as usize, contexts, bytes)
            })
        }),
        errors,
    )?;
    context_stream::capture_io(sink.finish(end as usize), errors)?;
    context_stream::capture_io(writer.flush(), errors)?;
    drop(writer);
    Ok((file, statistics.into_inner().unwrap()))
}

/// C++ structured-cluster `realign` overload.
pub fn realign_clusters<D, F>(
    clusters: &FlatArray<OId>,
    centroids: &[OId],
    database: &mut D,
    callback: F,
    hsp_values: HspValues,
    config: &ClusterOutputConfig,
    score_matrix: &mut ScoreMatrix,
) -> Result<Statistics, String>
where
    D: ClusterRealignDatabase,
    F: FnMut(&HspContext) -> Result<(), String>,
{
    realign_clusters_with_io(
        clusters,
        centroids,
        database,
        callback,
        hsp_values,
        config,
        score_matrix,
        &ClusterOutputIo::default(),
    )
}

/// Native command adapter preserving Source's bounded temporary stream lifecycle.
pub fn realign_clusters_with_io<D, F>(
    clusters: &FlatArray<OId>,
    centroids: &[OId],
    database: &mut D,
    mut callback: F,
    hsp_values: HspValues,
    config: &ClusterOutputConfig,
    score_matrix: &mut ScoreMatrix,
    io_config: &ClusterOutputIo,
) -> Result<Statistics, String>
where
    D: ClusterRealignDatabase,
    F: FnMut(&HspContext) -> Result<(), String>,
{
    realign_clusters_with_io_owner(
        clusters,
        centroids,
        database,
        |context| {
            callback(
                &context
                    .try_to_string_context()
                    .map_err(|error| error.to_string())?,
            )
        },
        hsp_values,
        config,
        score_matrix,
        io_config,
        None,
    )
}

pub fn realign_clusters_with_io_errors<D, F>(
    clusters: &FlatArray<OId>,
    centroids: &[OId],
    database: &mut D,
    mut callback: F,
    hsp_values: HspValues,
    config: &ClusterOutputConfig,
    score_matrix: &mut ScoreMatrix,
    io_config: &ClusterOutputIo,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<Statistics, String>
where
    D: ClusterRealignDatabase,
    F: FnMut(&HspContext) -> Result<(), String>,
{
    realign_clusters_with_io_owner(
        clusters,
        centroids,
        database,
        |context| {
            callback(
                &context
                    .try_to_string_context()
                    .map_err(|error| error.to_string())?,
            )
        },
        hsp_values,
        config,
        score_matrix,
        io_config,
        Some(errors),
    )
}

/// Byte-canonical command entry into the same Source realignment engine.
pub fn realign_clusters_with_io_errors_bytes<D, F>(
    clusters: &FlatArray<OId>,
    centroids: &[OId],
    database: &mut D,
    callback: F,
    hsp_values: HspValues,
    config: &ClusterOutputConfig,
    score_matrix: &mut ScoreMatrix,
    io_config: &ClusterOutputIo,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<Statistics, String>
where
    D: ClusterRealignDatabase,
    F: FnMut(&ByteHspContext) -> Result<(), String>,
{
    realign_clusters_with_io_owner(
        clusters,
        centroids,
        database,
        callback,
        hsp_values,
        config,
        score_matrix,
        io_config,
        Some(errors),
    )
}

fn realign_clusters_with_io_owner<D, F>(
    clusters: &FlatArray<OId>,
    centroids: &[OId],
    database: &mut D,
    mut callback: F,
    hsp_values: HspValues,
    config: &ClusterOutputConfig,
    score_matrix: &mut ScoreMatrix,
    io_config: &ClusterOutputIo,
    errors: Option<&crate::util::io::IoErrorOwner>,
) -> Result<Statistics, String>
where
    D: ClusterRealignDatabase,
    F: FnMut(&ByteHspContext) -> Result<(), String>,
{
    let mut statistics = Statistics::new();
    database.set_seqinfo_ptr(0)?;
    score_matrix.set_db_letters(config.db_size.unwrap_or_else(|| database.letters()));
    let lazy_titles = database.titles_lazy();
    database.set_load_flags(true, !lazy_titles);
    let mut centroid_offset = 0;
    while centroid_offset < database.sequence_count() {
        database.set_seqinfo_ptr(centroid_offset)?;
        let centroid_block = database.load_seqs(config.block_size)?;
        centroid_offset = database.tell_seq();
        database.set_seqinfo_ptr(0)?;
        let begin = centroids.partition_point(|oid| *oid < centroid_block.oid_begin()) as OId;
        let end = centroids.partition_point(|oid| *oid < centroid_block.oid_end()) as OId;
        let mut temporary_outputs = Vec::new();
        loop {
            let member_block = if centroid_block.seqs().len() as OId == database.sequence_count() {
                centroid_block.clone()
            } else {
                database.load_seqs(config.block_size)?
            };
            if member_block.empty() {
                break;
            }
            let cfg = Cfg::new(
                hsp_values,
                lazy_titles,
                clusters,
                centroids,
                database,
                centroid_block.clone(),
                member_block,
            );
            let (block_output, block_statistics) =
                run_block_pair(begin, end, config, score_matrix, &cfg, io_config, errors)?;
            temporary_outputs.push(block_output);
            statistics += &block_statistics;
            if centroid_block.seqs().len() as OId == cfg.database.lock().unwrap().sequence_count() {
                break;
            }
        }
        context_stream::merge_files(temporary_outputs, &mut callback, errors)?;
    }
    Ok(statistics)
}

/// Deterministic equivalent of C++ `cluster_sorted`, used by the mapping
/// overload declared beside this source in `cluster.h`.
pub fn cluster_sorted(clustering: &[OId]) -> (FlatArray<OId>, Vec<OId>) {
    let mut pairs: Vec<(OId, OId)> = clustering
        .iter()
        .copied()
        .enumerate()
        .map(|(member, centroid)| (centroid, member as OId))
        .collect();
    make_flat_array(&mut pairs)
}

/// C++ mapping-vector `realign` overload.
pub fn realign_mapping<D, F>(
    clustering: &[OId],
    database: &mut D,
    callback: F,
    hsp_values: HspValues,
    config: &ClusterOutputConfig,
    score_matrix: &mut ScoreMatrix,
) -> Result<Statistics, String>
where
    D: ClusterRealignDatabase,
    F: FnMut(&HspContext) -> Result<(), String>,
{
    let (clusters, centroids) = cluster_sorted(clustering);
    realign_clusters(
        &clusters,
        &centroids,
        database,
        callback,
        hsp_values,
        config,
        score_matrix,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::SequenceType;

    struct MemoryDatabase {
        block: Arc<Block>,
        cursor: OId,
        loaded: bool,
        lazy: bool,
        flags: Vec<(bool, bool)>,
        titles: Vec<String>,
    }

    impl MemoryDatabase {
        fn new(lazy: bool) -> Self {
            let mut block = Block::new();
            for (oid, (id, sequence)) in
                [("centroid", vec![0, 1, 2, 3]), ("member", vec![0, 1, 2, 3])]
                    .into_iter()
                    .enumerate()
            {
                block
                    .push_back(
                        &sequence,
                        Some(id),
                        None,
                        oid as OId,
                        SequenceType::AminoAcid,
                        1,
                        false,
                    )
                    .unwrap();
            }
            Self {
                block: Arc::new(block),
                cursor: 0,
                loaded: false,
                lazy,
                flags: Vec::new(),
                titles: vec!["lazy-centroid".to_owned(), "lazy-member".to_owned()],
            }
        }
    }

    impl ClusterRealignDatabase for MemoryDatabase {
        fn sequence_count(&self) -> OId {
            2
        }
        fn letters(&self) -> u64 {
            8
        }
        fn titles_lazy(&self) -> bool {
            self.lazy
        }
        fn set_load_flags(&mut self, sequences: bool, titles: bool) {
            self.flags.push((sequences, titles));
        }
        fn set_seqinfo_ptr(&mut self, oid: OId) -> Result<(), String> {
            self.cursor = oid;
            self.loaded = false;
            Ok(())
        }
        fn tell_seq(&self) -> OId {
            self.cursor
        }
        fn load_seqs(&mut self, _: u64) -> Result<Arc<Block>, String> {
            if self.loaded || self.cursor >= 2 {
                return Ok(Arc::new(Block::new()));
            }
            self.loaded = true;
            self.cursor = 2;
            Ok(self.block.clone())
        }
        fn seqid(&mut self, oid: OId) -> Result<String, String> {
            Ok(self.titles[oid as usize].clone())
        }
    }

    #[test]
    fn cluster_sorted_groups_members_and_orders_centroids() {
        let (clusters, centroids) = cluster_sorted(&[5, 2, 5, 2]);
        assert_eq!(centroids, [2, 5]);
        assert_eq!(clusters.range(0), &[1, 3]);
        assert_eq!(clusters.range(1), &[0, 2]);
    }

    #[test]
    fn whole_database_realign_emits_ordered_contexts_and_sets_flags() {
        let mut database = MemoryDatabase::new(false);
        let mut matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let (clusters, centroids) = cluster_sorted(&[0, 0]);
        let mut output = Vec::new();
        let statistics = realign_clusters(
            &clusters,
            &centroids,
            &mut database,
            |context| {
                output.push((
                    context.query_oid,
                    context.subject_oid,
                    context.query_title.clone(),
                    context.target_title.clone(),
                ));
                Ok(())
            },
            HspValues::COORDS,
            &ClusterOutputConfig::default(),
            &mut matrix,
        )
        .unwrap();
        assert_eq!(matrix.db_letters(), 8);
        assert_eq!(database.flags, [(true, true)]);
        assert_eq!(output.len(), 2);
        assert_eq!(output[0].0, 0);
        assert_eq!(output[0].2, "centroid");
        assert_eq!(output[1].3, "member");
        assert!(statistics.get(crate::basic::statistics::StatValue::SwipeTasksTotal) > 0);
    }

    #[test]
    fn lazy_titles_and_mapping_overload_are_preserved() {
        let mut database = MemoryDatabase::new(true);
        let mut matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let mut titles = Vec::new();
        realign_mapping(
            &[0, 0],
            &mut database,
            |context| {
                titles.push((context.query_title.clone(), context.target_title.clone()));
                Ok(())
            },
            HspValues::COORDS,
            &ClusterOutputConfig {
                db_size: Some(99),
                ..ClusterOutputConfig::default()
            },
            &mut matrix,
        )
        .unwrap();
        assert_eq!(database.flags, [(true, false)]);
        assert_eq!(matrix.db_letters(), 99);
        assert!(titles.iter().all(|(query, _)| query == "lazy-centroid"));
        assert!(titles.iter().any(|(_, target)| target == "lazy-member"));
    }
}
