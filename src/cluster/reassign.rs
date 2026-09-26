//! Cluster reassignment workflow from `diamond/src/cluster/reassign.cpp`.
//!
//! Upstream obtains configuration and the search/file implementations through
//! globals. Rust keeps those inputs explicit and uses [`ReassignBackend`] for
//! facilities implemented in other source files.

use crate::basic::value::OId;
use crate::cluster::cascaded::helpers::cluster_steps;
use crate::cluster::realign::DatabaseMetadata;
use crate::config::Sensitivity;
use crate::data::sequence_file::SequenceFileFlags;
use crate::output::edge::EdgeData;

pub const DEFAULT_MEMBER_COVER: f64 = 80.0;
pub const MAPBACK_NIL: OId = OId::MAX;

#[derive(Debug, Clone, PartialEq)]
pub struct ReassignConfig {
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

impl Default for ReassignConfig {
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
pub struct ReassignSearchConfig {
    pub command: &'static str,
    pub max_target_seqs: usize,
    pub output_format: &'static str,
    pub self_search: bool,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub sensitivity: Sensitivity,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubDatabase {
    Members,
    Centroids,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReassignSummary {
    pub database: DatabaseMetadata,
    pub coverage_cutoff: f64,
    pub centroid_count: usize,
    pub member_count: usize,
    pub reassigned_count: usize,
    pub search: ReassignSearchConfig,
}

pub trait ReassignBackend {
    fn open_database(
        &mut self,
        path: &str,
        flags: SequenceFileFlags,
    ) -> Result<DatabaseMetadata, String>;
    fn open_output(&mut self) -> Result<(), String>;
    fn read_clustering(&mut self, path: &str) -> Result<Vec<OId>, String>;
    /// Mirrors `sub_db(ids)` followed by `set_seqinfo_ptr(0)`.
    fn create_subdatabase(&mut self, kind: SubDatabase, ids: &[OId]) -> Result<(), String>;
    fn reset_statistics(&mut self);
    fn run_search(
        &mut self,
        config: &ReassignSearchConfig,
        consume: &mut dyn FnMut(&[u8]) -> Result<(), String>,
    ) -> Result<(), String>;
    fn init_random_access(&mut self) -> Result<(), String>;
    fn write_clustering(&mut self, clustering: &[OId]) -> Result<(), String>;
    fn close_database(&mut self) -> Result<(), String>;
}

/// C++ `split`: retain database order in both result vectors.
pub fn split(clustering: &[OId]) -> (Vec<OId>, Vec<OId>) {
    let mut centroids = Vec::new();
    let mut members = Vec::with_capacity(clustering.len());
    for (oid, &centroid) in clustering.iter().enumerate() {
        let oid = oid as OId;
        if centroid == oid {
            centroids.push(oid);
        } else {
            members.push(oid);
        }
    }
    (centroids, members)
}

/// Safe counterpart of the header-inline C++ `Mapback` consumer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Mapback {
    pub centroid_id: Vec<OId>,
}

impl Mapback {
    pub fn new(count: usize) -> Self {
        Self {
            centroid_id: vec![MAPBACK_NIL; count],
        }
    }

    pub fn consume(
        &mut self,
        bytes: &[u8],
        mutual_cover: Option<f64>,
        member_cover: f64,
    ) -> Result<(), String> {
        if bytes.len() % EdgeData::SIZE != 0 {
            return Err("Truncated edge record in reassignment output.".to_owned());
        }
        let cutoff = mutual_cover.unwrap_or(member_cover);
        let mut query = None;
        for record in bytes.chunks_exact(EdgeData::SIZE) {
            let edge = decode_edge(record);
            if query.is_some_and(|previous| previous != edge.query) {
                return Err("Reassignment edge batch contains multiple queries.".to_owned());
            }
            query = Some(edge.query);
            let slot = self
                .centroid_id
                .get_mut(edge.query as usize)
                .ok_or_else(|| format!("Reassignment query OID out of range: {}", edge.query))?;
            if f64::from(edge.qcovhsp) >= cutoff
                && mutual_cover.is_none_or(|_| f64::from(edge.scovhsp) >= cutoff)
            {
                *slot = edge.target;
            }
        }
        Ok(())
    }

    pub fn unmapped(&self) -> Vec<OId> {
        self.centroid_id
            .iter()
            .enumerate()
            .filter_map(|(i, &centroid)| (centroid == MAPBACK_NIL).then_some(i as OId))
            .collect()
    }
}

fn decode_edge(record: &[u8]) -> EdgeData {
    EdgeData {
        query: OId::from_ne_bytes(record[0..8].try_into().unwrap()),
        target: OId::from_ne_bytes(record[8..16].try_into().unwrap()),
        qcovhsp: f32::from_ne_bytes(record[16..20].try_into().unwrap()),
        scovhsp: f32::from_ne_bytes(record[20..24].try_into().unwrap()),
        evalue: f64::from_ne_bytes(record[24..32].try_into().unwrap()),
    }
}

/// Header-inline C++ `update_clustering`, with checked indices and intended
/// `Mapback::NIL_VALUE` handling made explicit.
pub fn update_clustering(
    clustering: &mut [OId],
    mapping: &[OId],
    members: &[OId],
    centroids: &[OId],
) -> Result<usize, String> {
    if mapping.len() != members.len() {
        return Err("Reassignment mapping/member count mismatch.".to_owned());
    }
    let mut changed = 0;
    for (&member, &centroid_index) in members.iter().zip(mapping) {
        if centroid_index == MAPBACK_NIL {
            continue;
        }
        let &centroid = centroids
            .get(centroid_index as usize)
            .ok_or_else(|| format!("Reassignment centroid index out of range: {centroid_index}"))?;
        let slot = clustering
            .get_mut(member as usize)
            .ok_or_else(|| format!("Reassignment member OID out of range: {member}"))?;
        if *slot != centroid {
            *slot = centroid;
            changed += 1;
        }
    }
    Ok(changed)
}

fn initialize_thresholds(config: &mut ReassignConfig) -> Result<(), String> {
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

/// Execute C++ `Cluster::reassign()` with global state made explicit.
pub fn reassign<B: ReassignBackend>(
    config: &mut ReassignConfig,
    backend: &mut B,
) -> Result<ReassignSummary, String> {
    if config.database.is_empty() {
        return Err("Database file is required for reassignment.".to_owned());
    }
    if config.clustering.is_empty() {
        return Err("Clustering file is required for reassignment.".to_owned());
    }
    initialize_thresholds(config)?;
    let coverage_cutoff = config
        .mutual_cover
        .or(config.member_cover)
        .expect("threshold initialized");
    let flags = SequenceFileFlags::NEED_LETTER_COUNT
        | SequenceFileFlags::ACC_TO_OID_MAPPING
        | SequenceFileFlags::OID_TO_ACC_MAPPING;
    let database = backend.open_database(&config.database, flags)?;
    config.db_size = Some(database.letters);
    let execution = (|| {
        backend.open_output()?;
        let mut clustering = backend.read_clustering(&config.clustering)?;
        if clustering.len() != database.sequence_count as usize {
            return Err("Invalid/incomplete clustering.".to_owned());
        }
        let (centroids, members) = split(&clustering);
        backend.create_subdatabase(SubDatabase::Members, &members)?;
        backend.create_subdatabase(SubDatabase::Centroids, &centroids)?;
        backend.reset_statistics();
        let steps = if config.cluster_steps.is_empty() {
            cluster_steps(config.approx_min_id.unwrap(), false)
        } else {
            config.cluster_steps.clone()
        };
        let sensitivity = parse_sensitivity(
            steps
                .last()
                .ok_or_else(|| "No clustering sensitivity step configured.".to_owned())?,
        )?;
        let search = ReassignSearchConfig {
            command: "blastp",
            max_target_seqs: 1,
            output_format: "edge",
            self_search: false,
            query_cover: coverage_cutoff,
            subject_cover: config.mutual_cover.map_or(0.0, |_| coverage_cutoff),
            sensitivity,
        };
        let mut mapback = Mapback::new(members.len());
        backend.run_search(&search, &mut |bytes| {
            mapback.consume(
                bytes,
                config.mutual_cover,
                config.member_cover.unwrap_or(DEFAULT_MEMBER_COVER),
            )
        })?;
        let reassigned_count =
            update_clustering(&mut clustering, &mapback.centroid_id, &members, &centroids)?;
        if database.titles_lazy {
            backend.init_random_access()?;
        }
        backend.write_clustering(&clustering)?;
        Ok(ReassignSummary {
            database,
            coverage_cutoff,
            centroid_count: centroids.len(),
            member_count: members.len(),
            reassigned_count,
            search,
        })
    })();
    let close = backend.close_database();
    match (execution, close) {
        (Err(error), _) => Err(error),
        (Ok(_), Err(error)) => Err(error),
        (Ok(summary), Ok(())) => Ok(summary),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Default)]
    struct FakeBackend {
        calls: Vec<String>,
        flags: Option<SequenceFileFlags>,
        clustering: Vec<OId>,
        edges: Vec<Vec<u8>>,
        written: Vec<OId>,
        search: Option<ReassignSearchConfig>,
        titles_lazy: bool,
        fail_search: bool,
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
    impl ReassignBackend for FakeBackend {
        fn open_database(
            &mut self,
            path: &str,
            flags: SequenceFileFlags,
        ) -> Result<DatabaseMetadata, String> {
            self.calls.push(format!("open:{path}"));
            self.flags = Some(flags);
            Ok(DatabaseMetadata {
                sequence_count: self.clustering.len() as u64,
                letters: 456,
                titles_lazy: self.titles_lazy,
            })
        }
        fn open_output(&mut self) -> Result<(), String> {
            self.calls.push("open_output".into());
            Ok(())
        }
        fn read_clustering(&mut self, path: &str) -> Result<Vec<OId>, String> {
            self.calls.push(format!("read:{path}"));
            Ok(self.clustering.clone())
        }
        fn create_subdatabase(&mut self, kind: SubDatabase, ids: &[OId]) -> Result<(), String> {
            self.calls.push(format!("sub:{kind:?}:{ids:?}"));
            Ok(())
        }
        fn reset_statistics(&mut self) {
            self.calls.push("reset_statistics".into());
        }
        fn run_search(
            &mut self,
            config: &ReassignSearchConfig,
            consume: &mut dyn FnMut(&[u8]) -> Result<(), String>,
        ) -> Result<(), String> {
            self.calls.push("search".into());
            self.search = Some(config.clone());
            if self.fail_search {
                return Err("search failed".into());
            }
            for bytes in self.edges.clone() {
                consume(&bytes)?;
            }
            Ok(())
        }
        fn init_random_access(&mut self) -> Result<(), String> {
            self.calls.push("random_access".into());
            Ok(())
        }
        fn write_clustering(&mut self, clustering: &[OId]) -> Result<(), String> {
            self.calls.push("write".into());
            self.written = clustering.to_vec();
            Ok(())
        }
        fn close_database(&mut self) -> Result<(), String> {
            self.calls.push("close".into());
            Ok(())
        }
    }
    fn config() -> ReassignConfig {
        ReassignConfig {
            database: "db.dmnd".into(),
            clustering: "clusters.tsv".into(),
            ..ReassignConfig::default()
        }
    }

    #[test]
    fn reassign_mirrors_partition_search_update_and_lazy_output_order() {
        let mut backend = FakeBackend {
            clustering: vec![0, 0, 2, 2, 4],
            edges: vec![edge(0, 2, 90.0, 5.0), edge(1, 0, 70.0, 99.0)],
            titles_lazy: true,
            ..FakeBackend::default()
        };
        let mut cfg = config();
        cfg.approx_min_id = Some(95.0);
        let summary = reassign(&mut cfg, &mut backend).unwrap();
        assert_eq!(backend.written, vec![0, 4, 2, 2, 4]);
        assert_eq!(
            (
                summary.reassigned_count,
                summary.centroid_count,
                summary.member_count
            ),
            (1, 3, 2)
        );
        assert_eq!(summary.search.sensitivity, Sensitivity::Fast);
        assert_eq!(summary.search.command, "blastp");
        assert_eq!(
            (summary.search.query_cover, summary.search.subject_cover),
            (80.0, 0.0)
        );
        assert_eq!(
            (cfg.db_size, cfg.diag_filter_id, cfg.diag_filter_cov),
            (Some(456), Some(85.0), Some(70.0))
        );
        assert_eq!(
            (cfg.soft_masking.as_deref(), cfg.masking.as_deref()),
            (Some("tantan"), Some("0"))
        );
        let flags = backend.flags.unwrap();
        assert!(
            flags.contains(SequenceFileFlags::NEED_LETTER_COUNT)
                && flags.contains(SequenceFileFlags::ACC_TO_OID_MAPPING)
                && flags.contains(SequenceFileFlags::OID_TO_ACC_MAPPING)
        );
        assert_eq!(
            backend.calls,
            [
                "open:db.dmnd",
                "open_output",
                "read:clusters.tsv",
                "sub:Members:[1, 3]",
                "sub:Centroids:[0, 2, 4]",
                "reset_statistics",
                "search",
                "random_access",
                "write",
                "close"
            ]
        );
    }

    #[test]
    fn mutual_coverage_requires_both_sides_and_suppresses_diag_defaults() {
        let mut backend = FakeBackend {
            clustering: vec![0, 0, 2],
            edges: vec![edge(0, 1, 95.0, 89.0)],
            ..FakeBackend::default()
        };
        let mut cfg = config();
        cfg.mutual_cover = Some(90.0);
        cfg.approx_min_id = Some(95.0);
        let summary = reassign(&mut cfg, &mut backend).unwrap();
        assert_eq!(backend.written, vec![0, 0, 2]);
        assert_eq!(summary.reassigned_count, 0);
        assert_eq!(
            (summary.search.query_cover, summary.search.subject_cover),
            (90.0, 90.0)
        );
        assert_eq!(
            (cfg.member_cover, cfg.diag_filter_id, cfg.diag_filter_cov),
            (None, None, None)
        );
    }

    #[test]
    fn mapback_validates_layout_grouping_and_bounds() {
        let mut mapback = Mapback::new(2);
        assert!(mapback
            .consume(&[0; 31], None, 80.0)
            .unwrap_err()
            .contains("Truncated"));
        let mut mixed = edge(0, 0, 90.0, 90.0);
        mixed.extend(edge(1, 0, 90.0, 90.0));
        assert!(mapback
            .consume(&mixed, None, 80.0)
            .unwrap_err()
            .contains("multiple"));
        assert!(Mapback::new(1)
            .consume(&edge(1, 0, 90.0, 90.0), None, 80.0)
            .unwrap_err()
            .contains("out of range"));
        assert_eq!(Mapback::new(3).unmapped(), vec![0, 1, 2]);
    }

    #[test]
    fn update_skips_nil_and_rejects_invalid_indices() {
        let mut clustering = vec![0, 0, 2];
        assert_eq!(
            update_clustering(&mut clustering, &[MAPBACK_NIL], &[1], &[0, 2]).unwrap(),
            0
        );
        assert!(update_clustering(&mut clustering, &[2], &[1], &[0, 2])
            .unwrap_err()
            .contains("centroid index"));
    }

    #[test]
    fn validates_config_and_closes_after_search_error() {
        let mut missing = ReassignConfig::default();
        assert!(reassign(&mut missing, &mut FakeBackend::default())
            .unwrap_err()
            .contains("Database"));
        let mut conflicting = config();
        conflicting.member_cover = Some(80.0);
        conflicting.mutual_cover = Some(80.0);
        assert!(reassign(&mut conflicting, &mut FakeBackend::default())
            .unwrap_err()
            .contains("mutually exclusive"));
        let mut backend = FakeBackend {
            clustering: vec![0, 0],
            fail_search: true,
            ..FakeBackend::default()
        };
        assert_eq!(
            reassign(&mut config(), &mut backend),
            Err("search failed".into())
        );
        assert_eq!(backend.calls.last().map(String::as_str), Some("close"));
    }
}
