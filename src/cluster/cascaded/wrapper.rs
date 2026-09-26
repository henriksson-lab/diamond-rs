//! Translation of `diamond/src/cluster/cascaded/wrapper.cpp`.

use crate::basic::value::{Loc, OId, SuperBlockId};
use crate::cluster::helpers::{init_thresholds, ClusterCommand, ThresholdConfig};
use crate::config::{block_size, Sensitivity};
use crate::output::edge::EdgeData;
use crate::search::sensitivity::get_traits;

use super::cascaded::CascadedConfig;
use super::helpers::{cluster_steps_with_config, Cascaded};

pub fn seq_mem_use(len: Loc, id_len: Loc, c: i32, min: i32, sketch_size: i32) -> i64 {
    assert!(min > 1 || sketch_size > 0);
    let seed_count = if min > 1 {
        len / (min / 2)
    } else {
        sketch_size
    };
    let mut extend_stage = seed_count * (15 + 16) + 12 + 2 * len;
    extend_stage /= 2;
    i64::from((len + 8 + id_len + 8 + seed_count * 9 / c + 8 + 8 + 4).max(extend_stage))
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BestCentroid {
    pub targets: Vec<OId>,
}

impl BestCentroid {
    pub fn new(block_size: OId) -> Result<Self, String> {
        Ok(Self {
            targets: vec![
                OId::MAX;
                usize::try_from(block_size)
                    .map_err(|_| "Super block size exceeds address space".to_owned())?
            ],
        })
    }

    pub fn consume(&mut self, bytes: &[u8]) -> Result<(), String> {
        if bytes.len() % EdgeData::SIZE != 0 {
            return Err("Invalid edge buffer size".to_owned());
        }
        for record in bytes.chunks_exact(EdgeData::SIZE) {
            let query = OId::from_ne_bytes(record[0..8].try_into().unwrap());
            let target = OId::from_ne_bytes(record[8..16].try_into().unwrap());
            let slot = self
                .targets
                .get_mut(query as usize)
                .ok_or_else(|| format!("Query OID out of super block range: {query}"))?;
            *slot = target;
        }
        Ok(())
    }

    pub fn finalize(&mut self) {}
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WrapperCommand {
    Cascaded,
    LinClust,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CascadedWrapperConfig {
    pub database: Option<String>,
    pub parallel_tmpdir: String,
    pub command: WrapperCommand,
    pub thresholds: ThresholdConfig,
    pub core: CascadedConfig,
    pub memory_limit: u64,
    pub threads: i32,
    pub hamming_ext: bool,
    pub freq_masking: bool,
    pub db_size: Option<u64>,
    pub oid_output: bool,
    pub simple_header: bool,
}

impl Default for CascadedWrapperConfig {
    fn default() -> Self {
        Self {
            database: None,
            parallel_tmpdir: String::new(),
            command: WrapperCommand::Cascaded,
            thresholds: ThresholdConfig {
                member_cover: None,
                mutual_cover: None,
                approx_min_id: None,
                soft_masking: None,
                masking: None,
                diag_filter_id: None,
                diag_filter_cov: None,
                command: ClusterCommand::Other,
            },
            core: CascadedConfig::default(),
            memory_limit: 16_000_000_000,
            threads: 1,
            hamming_ext: false,
            freq_masking: false,
            db_size: None,
            oid_output: false,
            simple_header: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CentroidSearchConfig {
    pub output_format: &'static str,
    pub self_search: bool,
    pub max_target_seqs: i64,
    pub top_percent: Option<f64>,
    pub sensitivity: Sensitivity,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_or_target_cover: f64,
    pub lin_stage1_target: bool,
    pub iterate: Vec<String>,
    pub lin_stage1_query: bool,
    pub chunk_size: f64,
    pub lowmem: i32,
}

pub struct LengthSortedBlock<B> {
    pub block: B,
    pub super_block_id_to_oid: Vec<OId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WrapperOutput {
    Dense(Vec<SuperBlockId>),
    Pairs(Vec<(OId, OId)>),
}

pub trait CascadedWrapperBackend {
    type Database;
    type Block;

    fn external(&mut self, config: &CascadedWrapperConfig) -> Result<(), String>;
    fn open_database(&mut self, path: &str) -> Result<Self::Database, String>;
    fn database_name(&self, database: &Self::Database) -> String;
    fn database_is_blast(&self, database: &Self::Database) -> bool;
    fn database_sequence_count(&self, database: &Self::Database) -> OId;
    fn database_letters(&self, database: &Self::Database) -> u64;
    fn length_sort(
        &mut self,
        database: &mut Self::Database,
        memory_limit: i64,
        minimizer_window: i32,
        sketch_size: i32,
    ) -> Result<Vec<LengthSortedBlock<Self::Block>>, String>;
    fn block_sequence_count(&self, block: &Self::Block) -> OId;
    fn block_letters(&self, block: &Self::Block) -> u64;
    fn centroid_sequence_count(&self) -> OId;
    fn centroid_letters(&self) -> u64;
    fn run_centroid_search(
        &mut self,
        block: &mut Self::Block,
        config: &CentroidSearchConfig,
        consumer: &mut BestCentroid,
    ) -> Result<(), String>;
    fn sub_database(
        &mut self,
        block: &Self::Block,
        ids: &[SuperBlockId],
    ) -> Result<Self::Block, String>;
    fn cascaded_database(
        &mut self,
        database: &mut Self::Database,
        linear: bool,
        config: &mut CascadedConfig,
    ) -> Result<Vec<SuperBlockId>, String>;
    fn cascaded_block(
        &mut self,
        block: &mut Self::Block,
        linear: bool,
        config: &mut CascadedConfig,
    ) -> Result<Vec<SuperBlockId>, String>;
    fn append_centroids(&mut self, block: &Self::Block, ids: &[SuperBlockId])
        -> Result<(), String>;
    fn close_block(&mut self, block: Self::Block) -> Result<(), String>;
    fn write_output(
        &mut self,
        database: &Self::Database,
        output: WrapperOutput,
        oid_output: bool,
        simple_header: bool,
    ) -> Result<(), String>;
    fn close_database(&mut self, database: Self::Database) -> Result<(), String>;
}

#[derive(Debug, Clone, PartialEq)]
pub struct WrapperState {
    pub linclust: bool,
    pub sensitivity: Sensitivity,
    pub centroid2oid: Vec<OId>,
    pub oid_to_centroid_oid: Vec<(OId, OId)>,
}

impl WrapperState {
    pub fn new(config: &CascadedWrapperConfig) -> Result<Self, String> {
        let linclust = config.command == WrapperCommand::LinClust;
        let steps = cluster_steps_with_config(
            config
                .thresholds
                .approx_min_id
                .unwrap_or(config.core.approx_min_id),
            linclust,
            &config.core.helpers,
        );
        let last = steps
            .last()
            .ok_or_else(|| "No cluster steps configured.".to_owned())?;
        let sensitivity = sensitivity_from_name(last.strip_suffix("_lin").unwrap_or(last))?;
        Ok(Self {
            linclust,
            sensitivity,
            centroid2oid: Vec::new(),
            oid_to_centroid_oid: Vec::new(),
        })
    }
}

fn sensitivity_from_name(name: &str) -> Result<Sensitivity, String> {
    match name {
        "faster" => Ok(Sensitivity::Faster),
        "fast" => Ok(Sensitivity::Fast),
        "default" => Ok(Sensitivity::Default),
        "linclust-40" => Ok(Sensitivity::Linclust40),
        "linclust-20" => Ok(Sensitivity::Linclust20),
        "shapes-6x10" => Ok(Sensitivity::Shapes6x10),
        "shapes-30x10" => Ok(Sensitivity::Shapes30x10),
        "mid-sensitive" => Ok(Sensitivity::MidSensitive),
        "sensitive" => Ok(Sensitivity::Sensitive),
        "more-sensitive" => Ok(Sensitivity::MoreSensitive),
        "very-sensitive" => Ok(Sensitivity::VerySensitive),
        "ultra-sensitive" => Ok(Sensitivity::UltraSensitive),
        _ => Err(format!("Invalid sensitivity level: {name}")),
    }
}

fn search_vs_centroids<B: CascadedWrapperBackend>(
    backend: &mut B,
    super_block: &mut B::Block,
    super_block_id_to_oid: &[OId],
    state: &mut WrapperState,
    config: &CascadedWrapperConfig,
) -> Result<Vec<SuperBlockId>, String> {
    let sequence_count = backend.block_sequence_count(super_block);
    if super_block_id_to_oid.len() != sequence_count as usize {
        return Err("Super block OID mapping has the wrong size.".to_owned());
    }
    let letters = backend
        .block_letters(super_block)
        .max(backend.centroid_letters());
    let (chunk_size, lowmem) = block_size(
        config.memory_limit.min(i64::MAX as u64) as i64,
        letters.min(i64::MAX as u64) as i64,
        state.sensitivity,
        state.linclust,
        config.threads,
    );
    let (query_cover, subject_cover) = if let Some(mutual) = config.thresholds.mutual_cover {
        (mutual, mutual)
    } else {
        (config.thresholds.member_cover.unwrap(), 0.0)
    };
    let search_config = CentroidSearchConfig {
        output_format: "edge",
        self_search: false,
        max_target_seqs: 1,
        top_percent: None,
        sensitivity: state.sensitivity,
        query_cover,
        subject_cover,
        query_or_target_cover: 0.0,
        lin_stage1_target: state.linclust,
        iterate: Vec::new(),
        lin_stage1_query: false,
        chunk_size,
        lowmem,
    };
    let mut best = BestCentroid::new(sequence_count)?;
    backend.run_centroid_search(super_block, &search_config, &mut best)?;
    best.finalize();

    let mut unaligned = Vec::new();
    for (i, (&oid, &target)) in super_block_id_to_oid.iter().zip(&best.targets).enumerate() {
        if target == OId::MAX {
            unaligned.push(i as SuperBlockId);
        } else {
            let centroid_oid = *state
                .centroid2oid
                .get(target as usize)
                .ok_or_else(|| format!("Centroid index out of range: {target}"))?;
            state.oid_to_centroid_oid.push((centroid_oid, oid));
        }
    }
    Ok(unaligned)
}

pub fn run<B: CascadedWrapperBackend>(
    backend: &mut B,
    config: &mut CascadedWrapperConfig,
) -> Result<(), String> {
    let database_path = config
        .database
        .clone()
        .ok_or_else(|| "Database is required.".to_owned())?;
    config.thresholds.command = match config.command {
        WrapperCommand::LinClust => ClusterCommand::LinClust,
        WrapperCommand::Cascaded => ClusterCommand::Other,
    };
    init_thresholds(&mut config.thresholds)?;
    config.core.member_cover = config.thresholds.member_cover.unwrap_or(80.0);
    config.core.mutual_cover = config.thresholds.mutual_cover;
    config.core.approx_min_id = config.thresholds.approx_min_id.unwrap();

    if !config.parallel_tmpdir.is_empty() {
        if config.command == WrapperCommand::LinClust {
            return backend.external(config);
        }
        return Err("Option is not permitted for this workflow: --parallel-tmpdir".to_owned());
    }
    config.hamming_ext = config.core.approx_min_id >= 50.0;

    let mut database = backend.open_database(&database_path)?;
    if backend.database_is_blast(&database) {
        return Err("Clustering is not supported for BLAST databases.".to_owned());
    }
    let letters = backend.database_letters(&database);
    let sequence_count = backend.database_sequence_count(&database);
    let memory_limit = config.memory_limit.min(i64::MAX as u64) as i64;
    let block_letters = (block_size(
        memory_limit,
        letters.min(i64::MAX as u64) as i64,
        Sensitivity::Faster,
        true,
        config.threads,
    )
    .0 * 1e9) as i64;

    let linear = config.command == WrapperCommand::LinClust;
    if block_letters as u64 >= letters && sequence_count < SuperBlockId::MAX as OId {
        let centroids = backend.cascaded_database(&mut database, linear, &mut config.core)?;
        backend.write_output(
            &database,
            WrapperOutput::Dense(centroids),
            config.oid_output,
            config.simple_header,
        )?;
    } else {
        let mut state = WrapperState::new(config)?;
        config.db_size = Some(letters);
        config.core.db_size = letters;
        let traits = get_traits(Sensitivity::Faster);
        let super_blocks = backend.length_sort(
            &mut database,
            memory_limit / 2,
            traits.minimizer_window,
            traits.sketch_size,
        )?;
        config.freq_masking = true;

        for (index, mut item) in super_blocks.into_iter().enumerate() {
            let block_count = backend.block_sequence_count(&item.block);
            if item.super_block_id_to_oid.len() != block_count as usize {
                return Err("Super block OID mapping has the wrong size.".to_owned());
            }
            let (mut unaligned_db, unaligned) = if index == 0 {
                let block_count = SuperBlockId::try_from(block_count)
                    .map_err(|_| "Super block exceeds SuperBlockId range.".to_owned())?;
                let ids = (0..block_count).collect::<Vec<_>>();
                (item.block, ids)
            } else {
                let unaligned = search_vs_centroids(
                    backend,
                    &mut item.block,
                    &item.super_block_id_to_oid,
                    &mut state,
                    config,
                )?;
                let sub = backend.sub_database(&item.block, &unaligned)?;
                backend.close_block(item.block)?;
                (sub, unaligned)
            };
            let clustering =
                backend.cascaded_block(&mut unaligned_db, state.linclust, &mut config.core)?;
            if clustering.len() != unaligned.len() {
                return Err("Cascaded clustering has the wrong size.".to_owned());
            }
            let mut centroids = Vec::new();
            for i in 0..unaligned.len() {
                let member_index = unaligned[i] as usize;
                let centroid_unaligned = clustering[i] as usize;
                let centroid_index = *unaligned.get(centroid_unaligned).ok_or_else(|| {
                    format!("Cascaded centroid index out of range: {}", clustering[i])
                })? as usize;
                let member_oid =
                    *item
                        .super_block_id_to_oid
                        .get(member_index)
                        .ok_or_else(|| {
                            format!("Super block member index out of range: {member_index}")
                        })?;
                let centroid_oid =
                    *item
                        .super_block_id_to_oid
                        .get(centroid_index)
                        .ok_or_else(|| {
                            format!("Super block centroid index out of range: {centroid_index}")
                        })?;
                state.oid_to_centroid_oid.push((centroid_oid, member_oid));
                if member_oid == centroid_oid {
                    state.centroid2oid.push(centroid_oid);
                    centroids.push(i as SuperBlockId);
                }
            }
            backend.append_centroids(&unaligned_db, &centroids)?;
            backend.close_block(unaligned_db)?;
        }
        backend.write_output(
            &database,
            WrapperOutput::Pairs(state.oid_to_centroid_oid),
            config.oid_output,
            config.simple_header,
        )?;
    }
    backend.close_database(database)
}

impl Cascaded {
    /// Explicit-backend counterpart of C++ `Cascaded::run`.
    pub fn run<B: CascadedWrapperBackend>(
        backend: &mut B,
        config: &mut CascadedWrapperConfig,
    ) -> Result<(), String> {
        run(backend, config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone)]
    struct Block {
        ids: Vec<OId>,
        letters: u64,
    }
    struct Db {
        count: OId,
        letters: u64,
        blast: bool,
    }

    #[derive(Default)]
    struct Backend {
        db_count: OId,
        db_letters: u64,
        blast: bool,
        blocks: Vec<LengthSortedBlock<Block>>,
        block_clusterings: Vec<Vec<SuperBlockId>>,
        dense: Vec<SuperBlockId>,
        centroid_ids: Vec<OId>,
        searches: Vec<CentroidSearchConfig>,
        output: Option<WrapperOutput>,
        external_calls: usize,
        closed_db: usize,
    }

    impl CascadedWrapperBackend for Backend {
        type Database = Db;
        type Block = Block;
        fn external(&mut self, _: &CascadedWrapperConfig) -> Result<(), String> {
            self.external_calls += 1;
            Ok(())
        }
        fn open_database(&mut self, _: &str) -> Result<Db, String> {
            Ok(Db {
                count: self.db_count,
                letters: self.db_letters,
                blast: self.blast,
            })
        }
        fn database_name(&self, _: &Db) -> String {
            "db".into()
        }
        fn database_is_blast(&self, db: &Db) -> bool {
            db.blast
        }
        fn database_sequence_count(&self, db: &Db) -> OId {
            db.count
        }
        fn database_letters(&self, db: &Db) -> u64 {
            db.letters
        }
        fn length_sort(
            &mut self,
            _: &mut Db,
            _: i64,
            min: i32,
            sketch: i32,
        ) -> Result<Vec<LengthSortedBlock<Block>>, String> {
            assert_eq!(min, 0);
            assert_eq!(sketch, 21);
            Ok(std::mem::take(&mut self.blocks))
        }
        fn block_sequence_count(&self, b: &Block) -> OId {
            b.ids.len() as OId
        }
        fn block_letters(&self, b: &Block) -> u64 {
            b.letters
        }
        fn centroid_sequence_count(&self) -> OId {
            self.centroid_ids.len() as OId
        }
        fn centroid_letters(&self) -> u64 {
            self.centroid_ids.len() as u64 * 10
        }
        fn run_centroid_search(
            &mut self,
            _: &mut Block,
            cfg: &CentroidSearchConfig,
            consumer: &mut BestCentroid,
        ) -> Result<(), String> {
            self.searches.push(cfg.clone());
            let edge = EdgeData {
                query: 0,
                target: 0,
                qcovhsp: 100.0,
                scovhsp: 100.0,
                evalue: 0.0,
            };
            let mut bytes = Vec::new();
            edge.write(&mut bytes).unwrap();
            consumer.consume(&bytes)
        }
        fn sub_database(&mut self, b: &Block, ids: &[SuperBlockId]) -> Result<Block, String> {
            Ok(Block {
                ids: ids.iter().map(|&i| b.ids[i as usize]).collect(),
                letters: ids.len() as u64 * 10,
            })
        }
        fn cascaded_database(
            &mut self,
            _: &mut Db,
            _: bool,
            _: &mut CascadedConfig,
        ) -> Result<Vec<SuperBlockId>, String> {
            Ok(self.dense.clone())
        }
        fn cascaded_block(
            &mut self,
            _: &mut Block,
            _: bool,
            _: &mut CascadedConfig,
        ) -> Result<Vec<SuperBlockId>, String> {
            Ok(self.block_clusterings.remove(0))
        }
        fn append_centroids(&mut self, b: &Block, ids: &[SuperBlockId]) -> Result<(), String> {
            self.centroid_ids
                .extend(ids.iter().map(|&i| b.ids[i as usize]));
            Ok(())
        }
        fn close_block(&mut self, _: Block) -> Result<(), String> {
            Ok(())
        }
        fn write_output(
            &mut self,
            _: &Db,
            output: WrapperOutput,
            _: bool,
            _: bool,
        ) -> Result<(), String> {
            self.output = Some(output);
            Ok(())
        }
        fn close_database(&mut self, _: Db) -> Result<(), String> {
            self.closed_db += 1;
            Ok(())
        }
    }

    fn config() -> CascadedWrapperConfig {
        CascadedWrapperConfig {
            database: Some("db".into()),
            memory_limit: 1_000_000,
            ..Default::default()
        }
    }

    #[test]
    fn memory_estimator_and_best_centroid_match_cpp() {
        assert_eq!(seq_mem_use(100, 20, 1, 0, 21), 431);
        let mut best = BestCentroid::new(2).unwrap();
        let edge = EdgeData {
            query: 1,
            target: 7,
            qcovhsp: 0.0,
            scovhsp: 0.0,
            evalue: 0.0,
        };
        let mut bytes = Vec::new();
        edge.write(&mut bytes).unwrap();
        best.consume(&bytes).unwrap();
        assert_eq!(best.targets, vec![OId::MAX, 7]);
    }

    #[test]
    fn external_and_in_memory_paths() {
        let mut external_cfg = config();
        external_cfg.command = WrapperCommand::LinClust;
        external_cfg.parallel_tmpdir = "tmp".into();
        let mut backend = Backend::default();
        run(&mut backend, &mut external_cfg).unwrap();
        assert_eq!(backend.external_calls, 1);
        let mut bad = config();
        bad.parallel_tmpdir = "tmp".into();
        assert_eq!(
            run(&mut backend, &mut bad).unwrap_err(),
            "Option is not permitted for this workflow: --parallel-tmpdir"
        );

        let mut backend = Backend {
            db_count: 3,
            db_letters: 30,
            dense: vec![0, 0, 2],
            ..Default::default()
        };
        let mut cfg = config();
        cfg.memory_limit = 16_000_000_000;
        run(&mut backend, &mut cfg).unwrap();
        assert_eq!(backend.output, Some(WrapperOutput::Dense(vec![0, 0, 2])));
        assert_eq!(backend.closed_db, 1);
    }

    #[test]
    fn length_sorted_path_composes_old_and_new_centroids() {
        let mut backend = Backend {
            db_count: 4,
            db_letters: 10_000_000,
            blocks: vec![
                LengthSortedBlock {
                    block: Block {
                        ids: vec![0, 1],
                        letters: 20,
                    },
                    super_block_id_to_oid: vec![0, 1],
                },
                LengthSortedBlock {
                    block: Block {
                        ids: vec![2, 3],
                        letters: 20,
                    },
                    super_block_id_to_oid: vec![2, 3],
                },
            ],
            block_clusterings: vec![vec![0, 0], vec![0]],
            ..Default::default()
        };
        let mut cfg = config();
        run(&mut backend, &mut cfg).unwrap();
        assert_eq!(backend.centroid_ids, vec![0, 3]);
        assert_eq!(backend.searches.len(), 1);
        assert!(!backend.searches[0].self_search);
        assert_eq!(backend.searches[0].max_target_seqs, 1);
        assert_eq!(
            backend.output,
            Some(WrapperOutput::Pairs(vec![(0, 0), (0, 1), (0, 2), (3, 3)]))
        );
        assert!(cfg.freq_masking);
        assert_eq!(cfg.db_size, Some(10_000_000));
    }
}
