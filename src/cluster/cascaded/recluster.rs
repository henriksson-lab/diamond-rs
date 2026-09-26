//! Recursive cluster refinement from `diamond/src/cluster/cascaded/recluster.cpp`.
//!
//! Database/search ownership and the mutable C++ configuration are explicit;
//! the recursive mapping algorithm and all search-policy mutations remain in
//! this source-mirrored module.

use crate::align::hsp::HspContext;
use crate::basic::value::{OId, SuperBlockId};
use crate::cluster::cascaded::helpers::cluster_steps;
use crate::cluster::realign::DatabaseMetadata;
use crate::cluster::reassign::{update_clustering, Mapback, DEFAULT_MEMBER_COVER};
use crate::config::Sensitivity;
use crate::data::sequence_file::SequenceFileFlags;
use crate::dp::swipe::HspValues;
use crate::util::data_structures::{make_flat_array, FlatArray};

#[derive(Debug, Clone, PartialEq)]
pub struct ReclusterConfig {
    pub database: String,
    pub clustering: String,
    pub member_cover: Option<f64>,
    pub mutual_cover: Option<f64>,
    pub approx_min_id: Option<f64>,
    pub soft_masking: Option<String>,
    pub masking: Option<String>,
    pub diag_filter_id: Option<f64>,
    pub diag_filter_cov: Option<f64>,
    pub cluster_steps: Vec<String>,
    pub db_size: Option<u64>,
}

impl Default for ReclusterConfig {
    fn default() -> Self {
        Self {
            database: String::new(),
            clustering: String::new(),
            member_cover: None,
            mutual_cover: None,
            approx_min_id: None,
            soft_masking: None,
            masking: None,
            diag_filter_id: None,
            diag_filter_cov: None,
            cluster_steps: Vec::new(),
            db_size: None,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReclusterSearchConfig {
    pub command: &'static str,
    pub max_target_seqs: usize,
    pub iterate: Vec<String>,
    pub output_format: &'static str,
    pub self_search: bool,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_or_target_cover: f64,
    pub sensitivity: Sensitivity,
    pub lowmem: usize,
    pub chunk_size: f64,
    pub lin_stage1_query: bool,
    pub lin_stage1_target: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReclusterSummary {
    pub database: DatabaseMetadata,
    pub coverage_cutoff: f64,
    pub clustering: Vec<OId>,
}

pub trait ReclusterBackend {
    type Database;

    fn open_database(
        &mut self,
        path: &str,
        flags: SequenceFileFlags,
    ) -> Result<(Self::Database, DatabaseMetadata), String>;
    fn close_database(&mut self, database: Self::Database) -> Result<(), String>;
    fn open_output(&mut self) -> Result<(), String>;
    fn read_clustering(
        &mut self,
        path: &str,
        database: &Self::Database,
    ) -> Result<Vec<OId>, String>;
    fn sequence_count(&self, database: &Self::Database) -> usize;
    fn create_subdatabase(
        &mut self,
        database: &Self::Database,
        ids: &[OId],
    ) -> Result<Self::Database, String>;
    fn realign(
        &mut self,
        database: &Self::Database,
        clusters: &FlatArray<OId>,
        centroids: &[OId],
        hsp_values: HspValues,
        callback: &mut dyn FnMut(&HspContext),
    ) -> Result<(), String>;
    fn reset_statistics(&mut self);
    fn run_search(
        &mut self,
        centroids: &Self::Database,
        queries: &Self::Database,
        config: &ReclusterSearchConfig,
        consume: &mut dyn FnMut(&[u8]) -> Result<(), String>,
    ) -> Result<(), String>;
    fn cascaded(
        &mut self,
        database: &Self::Database,
        linear: bool,
    ) -> Result<Vec<SuperBlockId>, String>;
    fn init_random_access(&mut self, database: &mut Self::Database) -> Result<(), String>;
    fn write_clustering(
        &mut self,
        database: &Self::Database,
        clustering: &[OId],
    ) -> Result<(), String>;
}

fn cluster_sorted(clustering: &[OId]) -> (FlatArray<OId>, Vec<OId>) {
    let mut pairs: Vec<(OId, OId)> = clustering
        .iter()
        .enumerate()
        .map(|(member, &centroid)| (centroid, member as OId))
        .collect();
    make_flat_array(&mut pairs)
}

fn initialize_thresholds(config: &mut ReclusterConfig) -> Result<(), String> {
    if config.member_cover.is_some() && config.mutual_cover.is_some() {
        return Err("--member-cover and --mutual-cover are mutually exclusive.".to_owned());
    }
    if config.mutual_cover.is_none() && config.member_cover.is_none() {
        config.member_cover = Some(DEFAULT_MEMBER_COVER);
    }
    let approx_min_id = *config.approx_min_id.get_or_insert(50.0);
    config
        .soft_masking
        .get_or_insert_with(|| "tantan".to_owned());
    config.masking.get_or_insert_with(|| "0".to_owned());
    if approx_min_id >= 90.0 && config.mutual_cover.is_none() {
        config.diag_filter_id.get_or_insert(approx_min_id - 10.0);
        let member_cover = config.member_cover.unwrap_or(DEFAULT_MEMBER_COVER);
        config
            .diag_filter_cov
            .get_or_insert(if member_cover > 50.0 {
                member_cover - 10.0
            } else {
                0.0
            });
    }
    Ok(())
}

fn parse_sensitivity(name: &str) -> Result<Sensitivity, String> {
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

fn search_config(config: &ReclusterConfig) -> Result<ReclusterSearchConfig, String> {
    let coverage = config.mutual_cover.or(config.member_cover).unwrap();
    let steps = if config.cluster_steps.is_empty() {
        cluster_steps(config.approx_min_id.unwrap(), false)
    } else {
        config.cluster_steps.clone()
    };
    Ok(ReclusterSearchConfig {
        command: "blastp",
        max_target_seqs: 1,
        iterate: Vec::new(),
        output_format: "edge",
        self_search: false,
        query_cover: coverage,
        subject_cover: config.mutual_cover.map_or(0.0, |_| coverage),
        query_or_target_cover: 0.0,
        sensitivity: parse_sensitivity(
            steps
                .last()
                .ok_or_else(|| "No clustering sensitivity step configured.".to_owned())?,
        )?,
        lowmem: 1,
        chunk_size: 4.0,
        lin_stage1_query: false,
        lin_stage1_target: false,
    })
}

/// Rust namespace for the private C++ `recluster(db, clustering, iteration)`
/// overload, since Rust does not support free-function overloading.
struct RecursiveRecluster;

impl RecursiveRecluster {
    fn recluster<B: ReclusterBackend>(
        backend: &mut B,
        database: &B::Database,
        clustering: &[OId],
        config: &ReclusterConfig,
        iteration: usize,
    ) -> Result<Vec<OId>, String> {
        let sequence_count = backend.sequence_count(database);
        if clustering.len() != sequence_count {
            return Err(format!(
            "Invalid clustering size in recluster iteration {}: expected {sequence_count}, got {}.",
            iteration + 1,
            clustering.len()
        ));
        }
        let (clusters, centroids) = cluster_sorted(clustering);
        let mut centroid_aligned = vec![false; sequence_count];
        for &centroid in &centroids {
            let slot = centroid_aligned
                .get_mut(centroid as usize)
                .ok_or_else(|| format!("Centroid OID out of range: {centroid}"))?;
            *slot = true;
        }
        let mut hsp_values = HspValues::TARGET_COORDS;
        if config.mutual_cover.is_some() {
            hsp_values = hsp_values | HspValues::QUERY_COORDS;
        }
        if config.approx_min_id.unwrap() > 0.0 {
            hsp_values =
                hsp_values | HspValues::QUERY_COORDS | HspValues::IDENT | HspValues::LENGTH;
        }
        let mutual_cover = config.mutual_cover;
        let member_cover = config.member_cover.unwrap_or(DEFAULT_MEMBER_COVER);
        let approx_min_id = config.approx_min_id.unwrap();
        let mut callback_error = None;
        let mut callback = |hsp: &HspContext| {
            let coverage = mutual_cover
                .is_some_and(|cutoff| hsp.qcovhsp() >= cutoff && hsp.scovhsp() >= cutoff)
                || (mutual_cover.is_none() && hsp.scovhsp() >= member_cover);
            if coverage && (hsp.approx_id() >= approx_min_id || hsp.id_percent() >= approx_min_id) {
                match centroid_aligned.get_mut(hsp.subject_oid as usize) {
                    Some(slot) => *slot = true,
                    None => {
                        callback_error = Some(format!(
                            "Realignment subject OID out of range: {}",
                            hsp.subject_oid
                        ));
                    }
                }
            }
        };
        backend.realign(database, &clusters, &centroids, hsp_values, &mut callback)?;
        if let Some(error) = callback_error {
            return Err(error);
        }

        let unaligned_members: Vec<OId> = centroid_aligned
            .iter()
            .enumerate()
            .filter_map(|(oid, &aligned)| (!aligned).then_some(oid as OId))
            .collect();
        if unaligned_members.is_empty() {
            return Ok(clustering.to_vec());
        }
        let unaligned = backend.create_subdatabase(database, &unaligned_members)?;
        let centroid_database = backend.create_subdatabase(database, &centroids)?;
        backend.reset_statistics();
        let search = search_config(config)?;
        let mut mapback = Mapback::new(unaligned_members.len());
        backend.run_search(&centroid_database, &unaligned, &search, &mut |bytes| {
            mapback.consume(bytes, config.mutual_cover, member_cover)
        })?;
        let mut output = clustering.to_vec();
        update_clustering(
            &mut output,
            &mapback.centroid_id,
            &unaligned_members,
            &centroids,
        )?;
        let unmapped_members = mapback.unmapped();
        if unmapped_members.is_empty() {
            return Ok(output);
        }

        let unmapped = backend.create_subdatabase(&unaligned, &unmapped_members)?;
        let initial: Vec<OId> = backend
            .cascaded(&unmapped, false)?
            .into_iter()
            .map(OId::from)
            .collect();
        let recursive = Self::recluster(backend, &unmapped, &initial, config, iteration + 1)?;
        if recursive.len() != unmapped_members.len() {
            return Err("Recursive reclustering mapping size mismatch.".to_owned());
        }
        for (i, &recursive_centroid) in recursive.iter().enumerate() {
            let member_index = *unmapped_members
                .get(i)
                .ok_or_else(|| "Recursive reclustering member index mismatch.".to_owned())?;
            let centroid_index = *unmapped_members
                .get(recursive_centroid as usize)
                .ok_or_else(|| {
                    format!("Recursive centroid index out of range: {recursive_centroid}")
                })?;
            let member = unaligned_members[member_index as usize];
            let centroid = unaligned_members[centroid_index as usize];
            output[member as usize] = centroid;
        }
        Ok(output)
    }
}

/// Execute public C++ `Cluster::recluster()`.
pub fn recluster<B: ReclusterBackend>(
    config: &mut ReclusterConfig,
    backend: &mut B,
) -> Result<ReclusterSummary, String> {
    if config.database.is_empty() {
        return Err("Database file is required for reclustering.".to_owned());
    }
    if config.clustering.is_empty() {
        return Err("Clustering file is required for reclustering.".to_owned());
    }
    initialize_thresholds(config)?;
    let coverage_cutoff = config.mutual_cover.or(config.member_cover).unwrap();
    let flags = SequenceFileFlags::NEED_LETTER_COUNT
        | SequenceFileFlags::ACC_TO_OID_MAPPING
        | SequenceFileFlags::OID_TO_ACC_MAPPING;
    let (mut database, metadata) = backend.open_database(&config.database, flags)?;
    config.db_size = Some(metadata.letters);
    let execution = (|| {
        backend.open_output()?;
        let clustering = backend.read_clustering(&config.clustering, &database)?;
        let clustering = RecursiveRecluster::recluster(backend, &database, &clustering, config, 0)?;
        if metadata.titles_lazy {
            backend.init_random_access(&mut database)?;
        }
        backend.write_clustering(&database, &clustering)?;
        Ok(ReclusterSummary {
            database: metadata,
            coverage_cutoff,
            clustering,
        })
    })();
    let close = backend.close_database(database);
    match (execution, close) {
        (Err(error), _) => Err(error),
        (Ok(_), Err(error)) => Err(error),
        (Ok(summary), Ok(())) => Ok(summary),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;
    use crate::output::edge::EdgeData;
    use crate::util::interval::Interval;

    #[derive(Debug)]
    struct TestDatabase {
        id: usize,
        count: usize,
    }

    #[derive(Default)]
    struct FakeBackend {
        calls: Vec<String>,
        input: Vec<OId>,
        next_database: usize,
        search_edges: Vec<Vec<u8>>,
        search_configs: Vec<ReclusterSearchConfig>,
        realign_values: Vec<HspValues>,
        realign_hits: Vec<HspContext>,
        cascaded_mapping: Vec<SuperBlockId>,
        written: Vec<OId>,
        flags: Option<SequenceFileFlags>,
        titles_lazy: bool,
    }

    fn edge(query: OId, target: OId, qcov: f32, scov: f32) -> Vec<u8> {
        let mut bytes = Vec::new();
        EdgeData {
            query,
            target,
            qcovhsp: qcov,
            scovhsp: scov,
            evalue: 0.0,
        }
        .write(&mut bytes)
        .unwrap();
        bytes
    }

    fn hit(subject: OId, qcov: i32, scov: i32, approx_id: f64, identity: i32) -> HspContext {
        let mut hsp = Hsp::new();
        hsp.query_source_range = Interval::new(0, qcov);
        hsp.subject_range = Interval::new(0, scov);
        hsp.approx_id = approx_id;
        hsp.identities = identity;
        hsp.length = 100;
        HspContext {
            subject_oid: subject,
            query_len: 100,
            subject_len: 100,
            hsp,
            ..HspContext::default()
        }
    }

    impl ReclusterBackend for FakeBackend {
        type Database = TestDatabase;

        fn open_database(
            &mut self,
            path: &str,
            flags: SequenceFileFlags,
        ) -> Result<(Self::Database, DatabaseMetadata), String> {
            self.calls.push(format!("open:{path}"));
            self.flags = Some(flags);
            Ok((
                TestDatabase {
                    id: 0,
                    count: self.input.len(),
                },
                DatabaseMetadata {
                    sequence_count: self.input.len() as u64,
                    letters: 999,
                    titles_lazy: self.titles_lazy,
                },
            ))
        }

        fn close_database(&mut self, database: Self::Database) -> Result<(), String> {
            self.calls.push(format!("close:{}", database.id));
            Ok(())
        }

        fn open_output(&mut self) -> Result<(), String> {
            self.calls.push("open_output".to_owned());
            Ok(())
        }

        fn read_clustering(&mut self, path: &str, _: &Self::Database) -> Result<Vec<OId>, String> {
            self.calls.push(format!("read:{path}"));
            Ok(self.input.clone())
        }

        fn sequence_count(&self, database: &Self::Database) -> usize {
            database.count
        }

        fn create_subdatabase(
            &mut self,
            database: &Self::Database,
            ids: &[OId],
        ) -> Result<Self::Database, String> {
            self.next_database += 1;
            self.calls.push(format!("sub:{}:{ids:?}", database.id));
            Ok(TestDatabase {
                id: self.next_database,
                count: ids.len(),
            })
        }

        fn realign(
            &mut self,
            database: &Self::Database,
            _: &FlatArray<OId>,
            _: &[OId],
            hsp_values: HspValues,
            callback: &mut dyn FnMut(&HspContext),
        ) -> Result<(), String> {
            self.calls.push(format!("realign:{}", database.id));
            self.realign_values.push(hsp_values);
            if database.id == 0 {
                for hit in &self.realign_hits {
                    callback(hit);
                }
            }
            Ok(())
        }

        fn reset_statistics(&mut self) {
            self.calls.push("reset_statistics".to_owned());
        }

        fn run_search(
            &mut self,
            centroids: &Self::Database,
            queries: &Self::Database,
            config: &ReclusterSearchConfig,
            consume: &mut dyn FnMut(&[u8]) -> Result<(), String>,
        ) -> Result<(), String> {
            self.calls
                .push(format!("search:{}:{}", centroids.id, queries.id));
            self.search_configs.push(config.clone());
            for bytes in self.search_edges.clone() {
                consume(&bytes)?;
            }
            Ok(())
        }

        fn cascaded(
            &mut self,
            database: &Self::Database,
            linear: bool,
        ) -> Result<Vec<SuperBlockId>, String> {
            self.calls
                .push(format!("cascaded:{}:{linear}", database.id));
            Ok(self.cascaded_mapping.clone())
        }

        fn init_random_access(&mut self, database: &mut Self::Database) -> Result<(), String> {
            self.calls.push(format!("random_access:{}", database.id));
            Ok(())
        }

        fn write_clustering(
            &mut self,
            database: &Self::Database,
            clustering: &[OId],
        ) -> Result<(), String> {
            self.calls.push(format!("write:{}", database.id));
            self.written = clustering.to_vec();
            Ok(())
        }
    }

    fn config() -> ReclusterConfig {
        ReclusterConfig {
            database: "db.dmnd".to_owned(),
            clustering: "clusters.tsv".to_owned(),
            ..ReclusterConfig::default()
        }
    }

    #[test]
    fn recursively_reclusters_unmapped_members_and_merges_local_oids() {
        let mut backend = FakeBackend {
            input: vec![0, 0, 2, 2, 4],
            search_edges: vec![edge(0, 2, 90.0, 1.0)],
            cascaded_mapping: vec![0],
            titles_lazy: true,
            ..FakeBackend::default()
        };
        let mut cfg = config();
        cfg.approx_min_id = Some(95.0);
        let summary = recluster(&mut cfg, &mut backend).unwrap();

        assert_eq!(summary.clustering, vec![0, 4, 2, 3, 4]);
        assert_eq!(backend.written, summary.clustering);
        assert_eq!(cfg.db_size, Some(999));
        let search = &backend.search_configs[0];
        assert_eq!(search.command, "blastp");
        assert_eq!(search.max_target_seqs, 1);
        assert_eq!(search.output_format, "edge");
        assert!(!search.self_search && search.iterate.is_empty());
        assert_eq!((search.query_cover, search.subject_cover), (80.0, 0.0));
        assert_eq!(search.query_or_target_cover, 0.0);
        assert_eq!(search.sensitivity, Sensitivity::Fast);
        assert_eq!((search.lowmem, search.chunk_size), (1, 4.0));
        assert!(!search.lin_stage1_query && !search.lin_stage1_target);
        assert!(backend.realign_values[0].all(
            HspValues::TARGET_COORDS
                | HspValues::QUERY_COORDS
                | HspValues::IDENT
                | HspValues::LENGTH
        ));
        assert!(backend.calls.contains(&"sub:1:[1]".to_owned()));
        assert!(backend.calls.contains(&"cascaded:3:false".to_owned()));
        assert_eq!(
            backend.calls[backend.calls.len() - 3..],
            ["random_access:0", "write:0", "close:0"]
        );
        let flags = backend.flags.unwrap();
        assert!(flags.contains(SequenceFileFlags::NEED_LETTER_COUNT));
        assert!(flags.contains(SequenceFileFlags::ACC_TO_OID_MAPPING));
        assert!(flags.contains(SequenceFileFlags::OID_TO_ACC_MAPPING));
    }

    #[test]
    fn realign_callback_accepts_identity_fallback_and_mutual_coverage() {
        let mut backend = FakeBackend {
            input: vec![0, 0, 2],
            realign_hits: vec![hit(1, 90, 90, 10.0, 95)],
            ..FakeBackend::default()
        };
        let mut cfg = config();
        cfg.mutual_cover = Some(80.0);
        cfg.approx_min_id = Some(90.0);
        let summary = recluster(&mut cfg, &mut backend).unwrap();
        assert_eq!(summary.clustering, vec![0, 0, 2]);
        assert!(backend.search_configs.is_empty());
        assert!(backend.realign_values[0].all(
            HspValues::TARGET_COORDS
                | HspValues::QUERY_COORDS
                | HspValues::IDENT
                | HspValues::LENGTH
        ));
        assert_eq!(cfg.diag_filter_id, None);
    }

    #[test]
    fn coverage_and_identity_failures_are_sent_to_recursive_clustering() {
        let mut backend = FakeBackend {
            input: vec![0, 0],
            realign_hits: vec![hit(1, 100, 79, 99.0, 99)],
            cascaded_mapping: vec![0],
            ..FakeBackend::default()
        };
        let mut cfg = config();
        cfg.member_cover = Some(80.0);
        cfg.approx_min_id = Some(90.0);
        let summary = recluster(&mut cfg, &mut backend).unwrap();
        assert_eq!(summary.clustering, vec![0, 1]);
        assert_eq!(backend.search_configs.len(), 1);
    }

    #[test]
    fn rejects_configuration_and_invalid_recursive_mapping() {
        let mut missing = ReclusterConfig::default();
        assert!(recluster(&mut missing, &mut FakeBackend::default())
            .unwrap_err()
            .contains("Database"));
        let mut conflict = config();
        conflict.member_cover = Some(80.0);
        conflict.mutual_cover = Some(80.0);
        assert!(recluster(&mut conflict, &mut FakeBackend::default())
            .unwrap_err()
            .contains("mutually exclusive"));

        let mut backend = FakeBackend {
            input: vec![0, 0],
            cascaded_mapping: vec![1],
            ..FakeBackend::default()
        };
        assert!(recluster(&mut config(), &mut backend)
            .unwrap_err()
            .contains("Centroid OID out of range"));
        assert_eq!(backend.calls.last().map(String::as_str), Some("close:0"));
    }
}
