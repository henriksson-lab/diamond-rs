//! Multi-node clustering orchestration from
//! `diamond/src/cluster/multinode/multinode.cpp`.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use super::data::{get_reps, GetRepsConfig, RepresentativeJob, SequenceRecord, SequenceSource};
use super::output::{merge, MergeConfig, MultinodeOutputJob, Volume};
use crate::basic::value::{OId, AMINO_ACID_ALPHABET, LETTER_MASK};
use crate::cluster::cascaded::helpers::{
    cluster_steps, default_round_approx_id, default_round_cov,
};
use crate::cluster::external::external::ClusterStats;
use crate::cluster::helpers::round_value;
use crate::config::Sensitivity;
use crate::masking::MaskingAlgo;
use crate::util::parallel::{Atomic, FileStack};
use crate::util::sequence::seqid;

pub const CASCADED_ROUND_MAX_EVALUE: f64 = 0.001;

#[derive(Debug)]
pub struct Job {
    pub max_oid: OId,
    pub mem_limit: u64,
    root_dir: PathBuf,
    worker_id: i64,
    round: i32,
    round_count: i32,
    input_count: Vec<u64>,
    start: Instant,
    log_file: FileStack,
}

impl Job {
    pub fn new(max_oid: OId, mem_limit: u64, root_dir: impl AsRef<Path>) -> Result<Self, String> {
        let root_dir = root_dir.as_ref().to_path_buf();
        fs::create_dir_all(root_dir.join("round0")).map_err(|error| error.to_string())?;
        let log_file = FileStack::new(root_dir.join("diamond_job.log"));
        let worker_id = Atomic::new(root_dir.join("worker_id")).fetch_add_one()?;
        Ok(Self {
            max_oid,
            mem_limit,
            root_dir,
            worker_id,
            round: 0,
            round_count: 0,
            input_count: Vec::new(),
            start: Instant::now(),
            log_file,
        })
    }

    pub const fn worker_id(&self) -> i64 {
        self.worker_id
    }

    pub fn root_dir(&self) -> &Path {
        &self.root_dir
    }

    pub fn base_dir(&self) -> PathBuf {
        self.base_dir_for_round(self.round)
    }

    pub fn base_dir_for_round(&self, round: i32) -> PathBuf {
        self.root_dir.join(format!("round{round}"))
    }

    /// C++ `Job::log(const char*, ...)` after formatting its varargs.
    pub fn log(&self, message: &str) -> Result<(), String> {
        let line = format!(
            "[{}, {}] {}\n",
            self.worker_id,
            self.start.elapsed().as_secs(),
            message
        );
        self.log_file.push_string(&line).map(|_| ())
    }

    /// C++ `Job::log(const ClusterStats&)` overload.
    pub fn log_stats(&self, stats: &ClusterStats) -> Result<(), String> {
        self.log(&format!(
            "Masked letters:   tantan: {}  seg: {}  motif: {}",
            stats.masking_stat.get(MaskingAlgo::Tantan),
            stats.masking_stat.get(MaskingAlgo::Seg),
            stats.masking_stat.get(MaskingAlgo::Motif)
        ))?;
        self.log(&format!("Seeds considered: {}", stats.seeds_considered))?;
        self.log(&format!("Seeds indexed: {}", stats.seeds_indexed))?;
        self.log(&format!(
            "Extensions computed: {}",
            stats.extensions_computed
        ))?;
        self.log(&format!(
            "Alignments passing e-value filter: {}",
            stats.hits_evalue_filtered
        ))?;
        self.log(&format!(
            "Alignments passing all filters: {}",
            stats.hits_filtered
        ))
    }

    pub fn next_round(&mut self) -> Result<(), String> {
        self.round += 1;
        fs::create_dir_all(self.base_dir()).map_err(|error| error.to_string())
    }

    pub const fn round(&self) -> i32 {
        self.round
    }

    pub fn set_round(&mut self, input_count: u64) {
        self.input_count.push(input_count);
    }

    pub fn sparse_input_count(&self, round: usize) -> u64 {
        self.input_count[round]
    }

    pub fn set_round_count(&mut self, count: i32) {
        self.round_count = count;
    }

    pub const fn round_count(&self) -> i32 {
        self.round_count
    }

    pub fn last_round(&self) -> bool {
        self.round == self.round_count - 1
    }
}

impl RepresentativeJob for Job {
    fn last_round(&self) -> bool {
        self.last_round()
    }

    fn base_dir(&self) -> PathBuf {
        self.base_dir()
    }

    fn round(&self) -> i32 {
        self.round()
    }

    fn log(&self, message: &str) -> Result<(), String> {
        self.log(message)
    }
}

impl MultinodeOutputJob for Job {
    fn max_oid(&self) -> OId {
        self.max_oid
    }

    fn root_dir(&self) -> &Path {
        self.root_dir()
    }

    fn base_dir(&self, round: i32) -> PathBuf {
        self.base_dir_for_round(round)
    }

    fn round(&self) -> i32 {
        self.round()
    }

    fn log(&mut self, message: &str) -> Result<(), String> {
        Job::log(self, message)
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct VolumedFile(pub Vec<Volume>);

impl VolumedFile {
    pub fn from_manifest(path: &Path) -> Result<Self, String> {
        let text = fs::read_to_string(path)
            .map_err(|_| format!("Error opening file {}", path.display()))?;
        let mut volumes = Vec::new();
        let mut next_oid = 0;
        for line in text.lines().filter(|line| !line.trim().is_empty()) {
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.len() < 2 {
                return Err("Format error in VolumedFile".to_owned());
            }
            let record_count = fields[1]
                .parse::<OId>()
                .map_err(|_| "Format error in VolumedFile".to_owned())?;
            let (oid_begin, oid_end) = if fields.len() >= 4 {
                (
                    fields[2]
                        .parse::<OId>()
                        .map_err(|_| "Format error in VolumedFile".to_owned())?,
                    fields[3]
                        .parse::<OId>()
                        .map_err(|_| "Format error in VolumedFile".to_owned())?,
                )
            } else {
                (next_oid, next_oid + record_count)
            };
            next_oid += record_count;
            volumes.push(Volume::new(fields[0], oid_begin, oid_end, record_count));
        }
        volumes.sort_by_key(|volume| volume.oid_begin);
        Ok(Self(volumes))
    }

    pub fn sparse_records(&self) -> OId {
        self.0.iter().map(|volume| volume.record_count).sum()
    }

    pub fn max_oid(&self) -> OId {
        self.0
            .iter()
            .filter_map(|volume| volume.oid_end.checked_sub(1))
            .max()
            .unwrap_or(0)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MultinodeConfig {
    pub output_file: PathBuf,
    pub oid_output: bool,
    pub root_dir: PathBuf,
    pub memory_limit: u64,
    pub threads: usize,
    pub mutual_cover: Option<f64>,
    pub member_cover: f64,
    pub approx_min_id: f64,
    pub sensitivity: Sensitivity,
    pub round_coverage: Vec<String>,
    pub round_approx_id: Vec<String>,
    pub max_evalue: f64,
    pub min_length_ratio: f64,
    pub file_buffer_size: usize,
}

impl Default for MultinodeConfig {
    fn default() -> Self {
        Self {
            output_file: PathBuf::new(),
            oid_output: false,
            root_dir: PathBuf::from("diamond-tmp"),
            memory_limit: 16_000_000_000,
            threads: 1,
            mutual_cover: None,
            member_cover: 80.0,
            approx_min_id: 50.0,
            sensitivity: Sensitivity::Default,
            round_coverage: Vec::new(),
            round_approx_id: Vec::new(),
            max_evalue: 0.001,
            min_length_ratio: 0.0,
            file_buffer_size: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BlockSearchRequest {
    pub command: &'static str,
    pub query_volume: usize,
    pub target_volume: usize,
    pub query_file: Option<PathBuf>,
    pub database: PathBuf,
    pub output_file: PathBuf,
    pub query_offset: OId,
    pub subject_offset: OId,
    pub numeric_ids: bool,
    pub output_fields: [&'static str; 6],
    pub query_cover: f64,
    pub subject_cover: f64,
    pub query_or_target_cover: f64,
    pub self_search: bool,
    pub lin_stage1_query: bool,
    pub lin_stage1_combo: bool,
    pub max_target_seqs: i64,
    pub top_percent: Option<f64>,
    pub double_indexed: bool,
    pub iterate_empty: bool,
    pub map_any: bool,
    pub lin_stage1_target: bool,
    pub lowmem: i32,
    pub chunk_size: u64,
    pub db_size: u64,
    pub comp_based_stats: i32,
    pub min_length_ratio: f64,
    pub sensitivity: Sensitivity,
    pub approx_min_id: f64,
    pub max_evalue: f64,
}

pub trait MultinodeBackend: SequenceSource {
    fn run_search(&mut self, request: &BlockSearchRequest) -> Result<(), String>;
    fn greedy_vertex_cover(
        &mut self,
        database: &Path,
        edges: &Path,
        output: &Path,
    ) -> Result<(), String>;
}

/// Number of upper-triangular block combinations, including the diagonal.
pub const fn combos(n: i64) -> i64 {
    n * (n + 1) / 2
}

pub const fn combo_to_rank(i: i64, j: i64, n: i64) -> i64 {
    (2 * n - i + 1) * i / 2 + j - i
}

pub fn rank_to_combo(rank: i64, n: i64) -> (i64, i64) {
    let a = 2 * n + 1;
    let i = (((2 * n + 1) as f64 - ((a * a - 8 * rank) as f64).sqrt()) / 2.0).floor() as i64;
    let j = i + rank - (2 * n - i + 1) * i / 2;
    (i, j)
}

pub fn run_block_combo<B: MultinodeBackend>(
    job: &Job,
    volumes: &VolumedFile,
    query_volume: usize,
    target_volume: usize,
    base_dir: &Path,
    config: &MultinodeConfig,
    backend: &mut B,
) -> Result<(), String> {
    let coverage_values = if config.round_coverage.is_empty() {
        default_round_cov(job.round_count())
    } else {
        config.round_coverage.clone()
    };
    let coverage = config.mutual_cover.unwrap_or(config.member_cover);
    // The source computes this value and then deliberately overwrites the
    // search coverage below; retain validation and arithmetic side effects.
    let _round_coverage = coverage.max(round_value(
        &coverage_values,
        "--round-coverage",
        job.round() as usize,
        job.round_count() as usize,
    )?);
    let mutual = config.mutual_cover;
    let same_block = query_volume == target_volume;
    let round_zero = job.round() == 0;
    let query = &volumes.0[query_volume];
    let target = &volumes.0[target_volume];
    let request = BlockSearchRequest {
        command: "blastp",
        query_volume,
        target_volume,
        query_file: (!same_block).then(|| query.path.clone()),
        database: target.path.clone(),
        output_file: base_dir.join(format!("{query_volume}_{target_volume}.tsv")),
        query_offset: if round_zero { query.oid_begin } else { 0 },
        subject_offset: if round_zero { target.oid_begin } else { 0 },
        numeric_ids: round_zero,
        output_fields: if round_zero {
            [
                "tab",
                "qnum",
                "snum",
                "qcovhsp",
                "scovhsp",
                "corrected_bitscore",
            ]
        } else {
            [
                "tab",
                "qseqid",
                "sseqid",
                "qcovhsp",
                "scovhsp",
                "corrected_bitscore",
            ]
        },
        query_cover: mutual.unwrap_or(config.member_cover),
        subject_cover: mutual.unwrap_or(0.0),
        query_or_target_cover: 0.0,
        self_search: same_block,
        lin_stage1_query: same_block,
        lin_stage1_combo: !same_block,
        max_target_seqs: i64::MAX,
        top_percent: None,
        double_indexed: true,
        iterate_empty: true,
        map_any: false,
        lin_stage1_target: false,
        lowmem: 1,
        chunk_size: 1024,
        db_size: 1_000_000_000,
        comp_based_stats: 0,
        min_length_ratio: config.min_length_ratio,
        sensitivity: config.sensitivity,
        approx_min_id: config.approx_min_id,
        max_evalue: config.max_evalue,
    };
    backend.run_search(&request)
}

pub fn run_block_combos<B: MultinodeBackend>(
    job: &Job,
    volumes: &VolumedFile,
    config: &MultinodeConfig,
    backend: &mut B,
) -> Result<String, String> {
    let base_dir = job.base_dir().join("alignments");
    fs::create_dir_all(&base_dir).map_err(|error| error.to_string())?;
    if job.round() == 0 {
        fs::create_dir_all(job.root_dir().join("accessions")).map_err(|error| error.to_string())?;
    }
    let n = volumes.0.len() as i64;
    let queue = Atomic::new(base_dir.join("queue"));
    let mut processed = 0_i64;
    loop {
        let rank = queue.fetch_add_one()?;
        if rank >= combos(n) {
            break;
        }
        let (query_volume, target_volume) = rank_to_combo(rank, n);
        job.log(&format!(
            "Searching blocks. Rank={}/{} Blocks={},{}",
            rank + 1,
            combos(n),
            query_volume,
            target_volume
        ))?;
        if job.round() == 0 && query_volume == target_volume {
            write_accessions(
                job,
                &volumes.0[query_volume as usize],
                query_volume,
                backend,
            )?;
        }
        run_block_combo(
            job,
            volumes,
            query_volume as usize,
            target_volume as usize,
            &base_dir,
            config,
            backend,
        )?;
        processed += 1;
    }
    let finished = Atomic::new(base_dir.join("finished"));
    finished.fetch_add(processed)?;
    finished.await_value(combos(n))?;
    let concat_lock = Atomic::new(base_dir.join("concat_lock"));
    let concat_done = Atomic::new(base_dir.join("concat_done"));
    if concat_lock.fetch_add_one()? == 0 {
        finalize_block_combos(job, volumes, backend, &base_dir, n)?;
        concat_done.fetch_add_one()?;
    } else {
        concat_done.await_value(1)?;
    }
    get_reps(
        job,
        &volumes.0,
        GetRepsConfig {
            threads: config.threads,
        },
        backend,
    )
}

fn finalize_block_combos<B: MultinodeBackend>(
    job: &Job,
    volumes: &VolumedFile,
    backend: &mut B,
    base_dir: &Path,
    n: i64,
) -> Result<(), String> {
    let alignment_path = job.base_dir().join("alignments.tsv");
    job.log(&format!(
        "Concatenating alignment files to {}",
        alignment_path.display()
    ))?;
    let mut alignment = BufWriter::new(
        File::create(&alignment_path)
            .map_err(|_| format!("Error opening file {}", alignment_path.display()))?,
    );
    for i in 0..n {
        for j in i..n {
            let path = base_dir.join(format!("{i}_{j}.tsv"));
            let bytes =
                fs::read(&path).map_err(|_| format!("Error opening file {}", path.display()))?;
            alignment
                .write_all(&bytes)
                .map_err(|_| format!("Error writing {}", alignment_path.display()))?;
        }
    }
    alignment.flush().map_err(|error| error.to_string())?;

    let oid_path = job.root_dir().join("oids.txt");
    if job.round() == 0 {
        job.log("Writing oid file")?;
        let mut text = String::new();
        for oid in 0..=volumes.max_oid() {
            text.push_str(&format!("{oid}\n"));
        }
        fs::write(&oid_path, text).map_err(|error| error.to_string())?;
    }
    job.log("Running greedy vertex cover")?;
    let database = if job.round() == 0 {
        oid_path
    } else {
        job.base_dir_for_round(job.round() - 1).join("rep_ids")
    };
    let clusters = job.base_dir().join("clusters.tsv");
    backend.greedy_vertex_cover(&database, &alignment_path, &clusters)?;
    if !job.last_round() {
        write_representative_ids(job, volumes, &clusters)?;
    }
    Ok(())
}

pub fn round<B: MultinodeBackend>(
    job: &mut Job,
    volumes: &VolumedFile,
    config: &mut MultinodeConfig,
    backend: &mut B,
) -> Result<String, String> {
    if let Some(mutual_cover) = config.mutual_cover {
        config.min_length_ratio = if config.sensitivity < Sensitivity::Linclust40 {
            (mutual_cover / 100.0 + 0.05).min(1.0)
        } else {
            mutual_cover / 100.0 - 0.05
        };
    }
    job.log(&format!(
        "Starting round {} sensitivity {:?}",
        job.round(),
        config.sensitivity
    ))?;
    job.set_round(volumes.sparse_records());
    run_block_combos(job, volumes, config, backend)
}

pub fn multinode<B: MultinodeBackend>(
    config: &mut MultinodeConfig,
    input: &VolumedFile,
    backend: &mut B,
) -> Result<OId, String> {
    if config.output_file.as_os_str().is_empty() {
        return Err("Option missing: output file (--out/-o)".to_owned());
    }
    if config.threads == 0 {
        return Err("Thread count must be positive.".to_owned());
    }
    config.file_buffer_size = 64 * 1024;
    let mut job = Job::new(input.max_oid(), config.memory_limit, &config.root_dir)?;
    if job.worker_id() == 0 {
        job.log(&format!(
            "{} coverage = {}",
            if config.mutual_cover.is_some() {
                "Bi-directional"
            } else {
                "Uni-directional"
            },
            config.mutual_cover.unwrap_or(config.member_cover)
        ))?;
        job.log(&format!("Approx. id = {}", config.approx_min_id))?;
        job.log(&format!("#Volumes = {}", input.0.len()))?;
        job.log(&format!("#Sequences = {}", input.sparse_records()))?;
    }
    let steps = cluster_steps(config.approx_min_id, true);
    job.set_round_count(steps.len() as i32);
    let evalue_cutoff = config.max_evalue;
    let target_approx_id = config.approx_min_id;
    let startup_lock = Atomic::new(job.root_dir().join("startup_lock"));
    let startup_done = Atomic::new(job.root_dir().join("startup_done"));
    if startup_lock.fetch_add_one()? == 0 {
        prepare_inputs(&job, input, backend)?;
        startup_done.fetch_add_one()?;
    } else {
        startup_done.await_value(1)?;
    }
    let input_volumes = VolumedFile::from_manifest(&job.root_dir().join("input.tsv"))?;
    let mut representatives = String::new();
    for (index, step) in steps.iter().enumerate() {
        config.sensitivity = sensitivity_from_step(step)?;
        let round_values = if config.round_approx_id.is_empty() {
            default_round_approx_id(job.round_count())
        } else {
            config.round_approx_id.clone()
        };
        config.approx_min_id = target_approx_id.max(round_value(
            &round_values,
            "--round-approx-id",
            job.round() as usize,
            job.round_count() as usize,
        )?);
        config.max_evalue = if index + 1 == steps.len() {
            evalue_cutoff
        } else {
            evalue_cutoff.min(CASCADED_ROUND_MAX_EVALUE)
        };
        let volumes = if index == 0 {
            input_volumes.clone()
        } else {
            VolumedFile::from_manifest(Path::new(&representatives))?
        };
        representatives = round(&mut job, &volumes, config, backend)?;
        if index + 1 < steps.len() {
            job.next_round()?;
        }
    }
    let output_lock = Atomic::new(job.root_dir().join("output_lock"));
    if output_lock.fetch_add_one()? == 0 {
        merge(
            &mut job,
            &input_volumes.0,
            &MergeConfig {
                output_file: config.output_file.clone(),
                oid_output: config.oid_output,
            },
        )
    } else {
        Ok(0)
    }
}

fn write_accessions<B: SequenceSource>(
    job: &Job,
    volume: &Volume,
    index: i64,
    backend: &B,
) -> Result<(), String> {
    let records = backend.read_records(&volume.path)?;
    let mut text = records
        .iter()
        .map(|record| seqid(&record.id))
        .collect::<Vec<_>>()
        .join("\n");
    if !records.is_empty() {
        text.push('\n');
    }
    fs::write(
        job.root_dir()
            .join("accessions")
            .join(format!("{index}.txt")),
        text,
    )
    .map_err(|error| error.to_string())
}

fn write_representative_ids(job: &Job, volumes: &VolumedFile, path: &Path) -> Result<(), String> {
    let text =
        fs::read_to_string(path).map_err(|_| format!("Error opening file {}", path.display()))?;
    let mut per_volume = vec![String::new(); volumes.0.len()];
    let mut all = String::new();
    let mut fields = text.split_whitespace();
    while let (Some(rep), Some(member)) = (fields.next(), fields.next()) {
        let rep = rep
            .parse::<OId>()
            .map_err(|_| "Invalid cluster representative.".to_owned())?;
        let member = member
            .parse::<OId>()
            .map_err(|_| "Invalid cluster member.".to_owned())?;
        if rep != member {
            continue;
        }
        let volume = volumes
            .0
            .iter()
            .position(|volume| member >= volume.oid_begin && member < volume.oid_end)
            .ok_or_else(|| format!("Cluster member OID out of volume range: {member}"))?;
        per_volume[volume].push_str(&format!("{rep}\n"));
        all.push_str(&format!("{rep}\n"));
    }
    for (index, text) in per_volume.into_iter().enumerate() {
        fs::write(job.base_dir().join(format!("rep_ids{index}")), text)
            .map_err(|error| error.to_string())?;
    }
    fs::write(job.base_dir().join("rep_ids"), all).map_err(|error| error.to_string())
}

fn prepare_inputs<B: SequenceSource>(
    job: &Job,
    volumes: &VolumedFile,
    backend: &B,
) -> Result<VolumedFile, String> {
    job.log(&format!("Memory limit = {}", job.mem_limit))?;
    let block_size = (job.mem_limit / 20).max(1);
    job.log(&format!("Block size = {block_size}"))?;
    let mut records = Vec::<SequenceRecord>::new();
    let mut letters = 0_u64;
    for (index, volume) in volumes.0.iter().enumerate() {
        job.log(&format!("Indexing volume {index}/{}", volumes.0.len()))?;
        let volume_records = backend.read_records(&volume.path)?;
        letters = letters.saturating_add(
            volume_records
                .iter()
                .map(|record| record.sequence.len() as u64)
                .sum::<u64>(),
        );
        records.extend(volume_records);
    }
    let mut blocks: Vec<Vec<SequenceRecord>> = vec![Vec::new()];
    let mut sizes = vec![0_u64];
    for record in records {
        let length = record.sequence.len() as u64;
        let last = sizes.len() - 1;
        if !blocks[last].is_empty() && sizes[last].saturating_add(length) > block_size {
            blocks.push(Vec::new());
            sizes.push(0);
        }
        let last = sizes.len() - 1;
        sizes[last] = sizes[last].saturating_add(length);
        blocks[last].push(record);
    }
    if blocks.len() == 1 && blocks[0].is_empty() {
        blocks.clear();
        sizes.clear();
    }
    let mut manifest = String::new();
    let mut output_volumes = Vec::new();
    let mut oid = 0_u64;
    let mut boundaries = vec![0_u64];
    for (index, block) in blocks.iter().enumerate() {
        let path = job.root_dir().join(format!("input{index}.faa"));
        let begin = oid;
        let mut output = Vec::new();
        for record in block {
            output.extend_from_slice(format!(">{oid}\n").as_bytes());
            for &letter in &record.sequence {
                output.push(
                    *AMINO_ACID_ALPHABET
                        .get((letter & LETTER_MASK) as usize)
                        .unwrap_or(&b'X'),
                );
            }
            output.push(b'\n');
            oid += 1;
        }
        fs::write(&path, output).map_err(|error| error.to_string())?;
        manifest.push_str(&format!("{}\t{}\n", path.display(), block.len()));
        output_volumes.push(Volume::new(path, begin, oid, block.len() as OId));
        boundaries.push(oid);
    }
    fs::write(job.root_dir().join("input.tsv"), manifest).map_err(|error| error.to_string())?;
    let mut boundary_bytes = Vec::with_capacity(boundaries.len() * 8);
    for boundary in &boundaries {
        boundary_bytes.extend_from_slice(&boundary.to_ne_bytes());
    }
    fs::write(job.base_dir().join("blocks"), boundary_bytes).map_err(|error| error.to_string())?;
    job.log(&format!(
        "Sequences in database = {}\nLetters in database = {}\nDatabase blocks:\n{}",
        volumes.max_oid() + 1,
        letters,
        output_volumes
            .iter()
            .zip(sizes)
            .map(|(volume, size)| format!("{}\t{}", volume.oid_begin, size))
            .collect::<Vec<_>>()
            .join("\n")
    ))?;
    Ok(VolumedFile(output_volumes))
}

fn sensitivity_from_step(step: &str) -> Result<Sensitivity, String> {
    match step.strip_suffix("_lin").unwrap_or(step) {
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
        value => Err(format!("Invalid sensitivity level: {value}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::Letter;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::sync::Mutex;

    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

    fn temp_dir(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "diamond_multinode_{name}_{}_{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    fn letters(text: &str) -> Vec<Letter> {
        text.bytes()
            .map(|byte| {
                AMINO_ACID_ALPHABET
                    .iter()
                    .position(|&aa| aa == byte)
                    .unwrap() as Letter
            })
            .collect()
    }

    #[derive(Default)]
    struct Backend {
        records: Mutex<std::collections::HashMap<PathBuf, Vec<SequenceRecord>>>,
        requests: Vec<BlockSearchRequest>,
    }

    impl SequenceSource for Backend {
        fn read_records(&self, path: &Path) -> Result<Vec<SequenceRecord>, String> {
            self.records
                .lock()
                .unwrap()
                .get(path)
                .cloned()
                .ok_or_else(|| format!("Missing fixture {}", path.display()))
        }
    }

    impl MultinodeBackend for Backend {
        fn run_search(&mut self, request: &BlockSearchRequest) -> Result<(), String> {
            self.requests.push(request.clone());
            fs::write(
                &request.output_file,
                format!("{}\t{}\n", request.query_volume, request.target_volume),
            )
            .map_err(|error| error.to_string())
        }

        fn greedy_vertex_cover(
            &mut self,
            database: &Path,
            _: &Path,
            output: &Path,
        ) -> Result<(), String> {
            let ids = fs::read_to_string(database).map_err(|error| error.to_string())?;
            let mut clustering = String::new();
            for id in ids.lines() {
                clustering.push_str(&format!("{id}\t{id}\n"));
            }
            fs::write(output, clustering).map_err(|error| error.to_string())
        }
    }

    #[test]
    fn combination_rank_is_bijective_in_source_order() {
        let mut seen = Vec::new();
        for rank in 0..combos(4) {
            let (i, j) = rank_to_combo(rank, 4);
            assert!(i <= j && j < 4);
            assert_eq!(combo_to_rank(i, j, 4), rank);
            seen.push((i, j));
        }
        assert_eq!(
            seen,
            [
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 3),
                (1, 1),
                (1, 2),
                (1, 3),
                (2, 2),
                (2, 3),
                (3, 3)
            ]
        );
    }

    #[test]
    fn job_methods_track_rounds_counts_and_both_log_overloads() {
        let root = temp_dir("job");
        let mut job = Job::new(9, 4096, &root).unwrap();
        job.set_round_count(2);
        job.set_round(10);
        job.log("hello").unwrap();
        job.log_stats(&ClusterStats {
            seeds_considered: 3,
            ..ClusterStats::default()
        })
        .unwrap();
        assert_eq!(job.worker_id(), 0);
        assert_eq!(job.sparse_input_count(0), 10);
        assert!(!job.last_round());
        job.next_round().unwrap();
        assert!(job.last_round());
        assert!(job.base_dir().is_dir());
        let log = fs::read_to_string(root.join("diamond_job.log")).unwrap();
        assert!(log.contains("hello"));
        assert!(log.contains("Seeds considered: 3"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn block_combo_preserves_round_zero_offsets_and_search_modes() {
        let root = temp_dir("combo");
        let mut job = Job::new(7, 4096, &root).unwrap();
        job.set_round_count(2);
        let volumes = VolumedFile(vec![
            Volume::new("a.faa", 0, 4, 4),
            Volume::new("b.faa", 4, 8, 4),
        ]);
        let config = MultinodeConfig {
            root_dir: root.clone(),
            member_cover: 73.0,
            ..MultinodeConfig::default()
        };
        let mut backend = Backend::default();
        let out = job.base_dir().join("alignments");
        fs::create_dir_all(&out).unwrap();
        run_block_combo(&job, &volumes, 0, 0, &out, &config, &mut backend).unwrap();
        run_block_combo(&job, &volumes, 0, 1, &out, &config, &mut backend).unwrap();
        let diagonal = &backend.requests[0];
        assert!(diagonal.self_search && diagonal.lin_stage1_query);
        assert_eq!(diagonal.query_file, None);
        assert_eq!(diagonal.command, "blastp");
        assert_eq!(diagonal.output_fields[1..3], ["qnum", "snum"]);
        assert!(diagonal.double_indexed && diagonal.iterate_empty);
        assert!(!diagonal.map_any && !diagonal.lin_stage1_target);
        assert_eq!((diagonal.lowmem, diagonal.chunk_size), (1, 1024));
        let off_diagonal = &backend.requests[1];
        assert!(!off_diagonal.self_search && off_diagonal.lin_stage1_combo);
        assert_eq!(off_diagonal.query_offset, 0);
        assert_eq!(off_diagonal.subject_offset, 4);
        assert_eq!(off_diagonal.query_cover, 73.0);
        assert_eq!(off_diagonal.query_or_target_cover, 0.0);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn run_block_combos_writes_accessions_concatenates_and_materializes_reps() {
        let root = temp_dir("blocks");
        let mut job = Job::new(3, 4096, &root).unwrap();
        job.set_round_count(2);
        let a = root.join("a.faa");
        let b = root.join("b.faa");
        let volumes = VolumedFile(vec![Volume::new(&a, 0, 2, 2), Volume::new(&b, 2, 4, 2)]);
        let mut backend = Backend::default();
        backend.records.lock().unwrap().insert(
            a,
            vec![
                SequenceRecord {
                    id: "a0 title".into(),
                    sequence: letters("ARND"),
                },
                SequenceRecord {
                    id: "a1".into(),
                    sequence: letters("CQ"),
                },
            ],
        );
        backend.records.lock().unwrap().insert(
            b,
            vec![
                SequenceRecord {
                    id: "b0".into(),
                    sequence: letters("EG"),
                },
                SequenceRecord {
                    id: "b1".into(),
                    sequence: letters("HI"),
                },
            ],
        );
        let config = MultinodeConfig {
            root_dir: root.clone(),
            threads: 1,
            ..MultinodeConfig::default()
        };
        let reps = run_block_combos(&job, &volumes, &config, &mut backend).unwrap();
        assert_eq!(backend.requests.len(), 3);
        assert_eq!(
            fs::read_to_string(root.join("accessions/0.txt")).unwrap(),
            "a0\na1\n"
        );
        assert_eq!(
            fs::read_to_string(job.base_dir().join("alignments.tsv")).unwrap(),
            "0\t0\n0\t1\n1\t1\n"
        );
        assert!(fs::read_to_string(job.base_dir().join("rep_ids"))
            .unwrap()
            .contains("3\n"));
        assert!(fs::read_to_string(reps).unwrap().contains("\t2\t0\t2"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn two_workers_claim_each_combo_once_and_share_the_concat_barrier() {
        let root = temp_dir("distributed");
        let mut first_job = Job::new(3, 4096, &root).unwrap();
        let mut second_job = Job::new(3, 4096, &root).unwrap();
        assert_eq!((first_job.worker_id(), second_job.worker_id()), (0, 1));
        first_job.set_round_count(1);
        second_job.set_round_count(1);
        let a = root.join("a.faa");
        let b = root.join("b.faa");
        let volumes = VolumedFile(vec![Volume::new(&a, 0, 2, 2), Volume::new(&b, 2, 4, 2)]);
        let make_backend = || {
            let backend = Backend::default();
            backend.records.lock().unwrap().insert(
                a.clone(),
                vec![
                    SequenceRecord {
                        id: "a0".into(),
                        sequence: letters("AR"),
                    },
                    SequenceRecord {
                        id: "a1".into(),
                        sequence: letters("ND"),
                    },
                ],
            );
            backend.records.lock().unwrap().insert(
                b.clone(),
                vec![
                    SequenceRecord {
                        id: "b0".into(),
                        sequence: letters("CQ"),
                    },
                    SequenceRecord {
                        id: "b1".into(),
                        sequence: letters("EG"),
                    },
                ],
            );
            backend
        };
        let config = MultinodeConfig {
            root_dir: root.clone(),
            threads: 1,
            ..MultinodeConfig::default()
        };
        let (first, second) = std::thread::scope(|scope| {
            let volumes1 = volumes.clone();
            let config1 = config.clone();
            let h1 = scope.spawn(move || {
                let mut backend = make_backend();
                let result = run_block_combos(&first_job, &volumes1, &config1, &mut backend);
                (result, backend.requests.len())
            });
            let h2 = scope.spawn(move || {
                let mut backend = make_backend();
                let result = run_block_combos(&second_job, &volumes, &config, &mut backend);
                (result, backend.requests.len())
            });
            (h1.join().unwrap(), h2.join().unwrap())
        });
        first.0.unwrap();
        second.0.unwrap();
        assert_eq!(first.1 + second.1, combos(2) as usize);
        assert_eq!(
            fs::read_to_string(root.join("round0/alignments.tsv")).unwrap(),
            "0\t0\n0\t1\n1\t1\n"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn multinode_prepares_blocks_runs_round_and_merges_accessions() {
        let root = temp_dir("workflow");
        let source = root.join("source.faa");
        let input = VolumedFile(vec![Volume::new(&source, 0, 2, 2)]);
        let mut backend = Backend::default();
        backend.records.lock().unwrap().insert(
            source,
            vec![
                SequenceRecord {
                    id: "alpha description".into(),
                    sequence: letters("ARND"),
                },
                SequenceRecord {
                    id: "beta".into(),
                    sequence: letters("CQEG"),
                },
            ],
        );
        // Prepared input paths are read again during the block-combo phase.
        backend.records.lock().unwrap().insert(
            root.join("input0.faa"),
            vec![
                SequenceRecord {
                    id: "0".into(),
                    sequence: letters("ARND"),
                },
                SequenceRecord {
                    id: "1".into(),
                    sequence: letters("CQEG"),
                },
            ],
        );
        let output = root.join("clusters.out");
        let mut config = MultinodeConfig {
            output_file: output.clone(),
            root_dir: root.clone(),
            memory_limit: 1_000_000,
            approx_min_id: 90.0,
            ..MultinodeConfig::default()
        };
        assert_eq!(multinode(&mut config, &input, &mut backend).unwrap(), 2);
        assert_eq!(config.file_buffer_size, 64 * 1024);
        assert_eq!(fs::read_to_string(output).unwrap(), "0\t0\n1\t1\n");
        assert!(jobless_blocks(&root).len() >= 16);
        fs::remove_dir_all(root).unwrap();
    }

    fn jobless_blocks(root: &Path) -> Vec<u8> {
        fs::read(root.join("round0/blocks")).unwrap()
    }
}
