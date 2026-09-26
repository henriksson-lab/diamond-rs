//! Cluster-member realignment from `diamond/src/cluster/output.cpp`.

use std::sync::{Arc, Mutex};

use crate::align::hsp::HspContext;
use crate::basic::statistics::Statistics;
use crate::basic::value::{CentroidId, Letter, OId};
use crate::data::block::Block;
use crate::dp::swipe::{
    bin as swipe_bin, swipe, targets as make_dp_targets, CarryOver, DpTarget, Flags, HspValues,
    Params, Targets,
};
use crate::stats::cbs::hauser_correction;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::data_structures::{make_flat_array, FlatArray};

pub trait ClusterRealignDatabase {
    fn sequence_count(&self) -> OId;
    fn letters(&self) -> u64;
    fn titles_lazy(&self) -> bool;
    fn set_load_flags(&mut self, sequences: bool, titles: bool);
    fn set_seqinfo_ptr(&mut self, oid: OId) -> Result<(), String>;
    fn tell_seq(&self) -> OId;
    fn load_seqs(&mut self, block_size: u64) -> Result<Arc<Block>, String>;
    fn seqid(&mut self, oid: OId) -> Result<String, String>;
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClusterOutputConfig {
    pub block_size: u64,
    pub db_size: Option<u64>,
    pub comp_based_stats_hauser: bool,
    pub cutoff_score_8bit: i32,
    pub max_swipe_dp: i64,
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
            cutoff_score_8bit: 240,
            max_swipe_dp: 1_000_000,
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
    database: &'a mut D,
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
            database,
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
    cfg: &mut Cfg<'_, D>,
) -> Result<Vec<HspContext>, String> {
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
        cfg.database.seqid(centroid_oid)?
    } else {
        String::from_utf8_lossy(cfg.centroid_block.ids()?.get(centroid_id as usize)).into_owned()
    };
    let correction = if config.comp_based_stats_hauser {
        Some(hauser_correction(centroid_seq, score_matrix))
    } else {
        None
    };
    let mut params = Params::new(centroid_seq, score_matrix);
    params.query_id = Some(&centroid_title);
    params.frame = 0;
    params.query_source_len = centroid_seq.len() as i32;
    params.composition_bias = correction.as_deref();
    params.flags = Flags::FULL_MATRIX;
    params.v = cfg.hsp_values;
    params.cutoff_score_8bit = config.cutoff_score_8bit;
    params.max_swipe_dp = config.max_swipe_dp;
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
    *statistics += &swipe_statistics.lock().unwrap();
    for hsp in hsps {
        let subject_block_id = hsp.swipe_target as usize;
        let subject_oid = cfg.member_block.block_id2oid(subject_block_id as u32);
        let subject_title = if cfg.lazy_titles {
            cfg.database.seqid(subject_oid)?
        } else {
            String::from_utf8_lossy(cfg.member_block.ids()?.get(subject_block_id)).into_owned()
        };
        output.push(HspContext::new(
            hsp,
            centroid_id,
            centroid_oid,
            vec![centroid_seq.to_vec()],
            centroid_seq.len() as i32,
            centroid_title.clone(),
            subject_oid,
            cfg.member_block.seqs().length(subject_block_id),
            subject_title,
            0,
            0,
            Vec::<Letter>::new(),
            0.0,
            0.0,
        ));
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
    cfg: &mut Cfg<'_, D>,
) -> Result<(Vec<HspContext>, Statistics), String> {
    let mut output = Vec::new();
    let mut total_statistics = Statistics::new();
    for centroid in begin..end {
        let mut statistics = Statistics::new();
        output.extend(align_centroid(
            centroid,
            &mut statistics,
            config,
            score_matrix,
            cfg,
        )?);
        total_statistics += &statistics;
    }
    Ok((output, total_statistics))
}

/// C++ structured-cluster `realign` overload.
pub fn realign_clusters<D, F>(
    clusters: &FlatArray<OId>,
    centroids: &[OId],
    database: &mut D,
    mut callback: F,
    hsp_values: HspValues,
    config: &ClusterOutputConfig,
    score_matrix: &mut ScoreMatrix,
) -> Result<Statistics, String>
where
    D: ClusterRealignDatabase,
    F: FnMut(&HspContext) -> Result<(), String>,
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
            let mut cfg = Cfg::new(
                hsp_values,
                lazy_titles,
                clusters,
                centroids,
                database,
                centroid_block.clone(),
                member_block,
            );
            let (block_output, block_statistics) =
                run_block_pair(begin, end, config, score_matrix, &mut cfg)?;
            temporary_outputs.push(block_output);
            statistics += &block_statistics;
            if centroid_block.seqs().len() as OId == cfg.database.sequence_count() {
                break;
            }
        }
        let mut merged: Vec<HspContext> = temporary_outputs.into_iter().flatten().collect();
        merged.sort_by_key(|context| context.query_oid);
        for context in &merged {
            callback(context)?;
        }
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
