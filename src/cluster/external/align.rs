//! External-clustering alignment stage from `cluster/external/align.cpp`.

use std::collections::HashMap;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use super::external::{ClusterStats, Edge, Job, PairEntryShort};
use super::pair_table::{Bucket, RadixedTable, RADIX_BITS, RADIX_COUNT};
use super::seed_table::{SequenceVolumeReader, Volume, VolumedFile};
use crate::basic::statistics::Statistics;
use crate::basic::value::Letter;
use crate::dp::swipe::{
    bin as swipe_bin, swipe, targets, CarryOver, DpTarget, Flags, HspValues, Params,
};
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::hash::hash64;
use crate::util::io::CompressedBuffer;

#[derive(Debug, Clone, PartialEq)]
pub struct AlignConfig {
    pub threads: usize,
    pub approx_min_identity: f64,
    pub member_cover: f64,
    pub mutual_cover: Option<f64>,
    pub alignment_output: Option<PathBuf>,
    pub cutoff_score_8bit: i32,
    pub max_swipe_dp: i64,
    pub approx_backtrace: bool,
    pub max_evalue: f64,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            approx_min_identity: 0.0,
            member_cover: 80.0,
            mutual_cover: None,
            alignment_output: None,
            cutoff_score_8bit: i8::MAX as i32,
            max_swipe_dp: 1_000_000,
            approx_backtrace: false,
            max_evalue: f64::MAX,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct ChunkSequences {
    sequences: HashMap<u64, Vec<Letter>>,
    letters: usize,
    volumes: usize,
}

impl ChunkSequences {
    pub fn from_volumes<R: SequenceVolumeReader>(
        volumes: &VolumedFile,
        reader: &mut R,
    ) -> Result<Self, String> {
        let mut sequences = HashMap::new();
        let mut letters = 0usize;
        for volume in &volumes.0 {
            for sequence in reader.read_volume(&volume.path)? {
                let oid = parse_oid(&sequence.id)?;
                letters += sequence.sequence.len();
                sequences.insert(oid, sequence.sequence);
            }
        }
        Ok(Self {
            sequences,
            letters,
            volumes: volumes.0.len(),
        })
    }

    pub fn get(&self, oid: u64) -> Result<&[Letter], String> {
        self.sequences
            .get(&oid)
            .map(Vec::as_slice)
            .ok_or_else(|| "ChunkSeqs".to_string())
    }

    pub fn oids(&self) -> usize {
        self.sequences.len()
    }

    pub fn letters(&self) -> usize {
        self.letters
    }

    pub fn volumes(&self) -> usize {
        self.volumes
    }
}

fn parse_oid(id: &str) -> Result<u64, String> {
    let token = id.split_whitespace().next().unwrap_or("");
    let digits = token
        .strip_prefix('+')
        .unwrap_or(token)
        .chars()
        .take_while(char::is_ascii_digit)
        .collect::<String>();
    if digits.is_empty() {
        return Err(format!("Invalid sequence OID: {id}"));
    }
    digits
        .parse::<u64>()
        .map_err(|_| format!("Invalid sequence OID: {id}"))
}

#[derive(Debug)]
struct EdgeFileArray {
    base_dir: PathBuf,
    worker_id: i64,
    edges: Vec<Vec<Edge>>,
}

impl EdgeFileArray {
    fn new(base_dir: impl AsRef<Path>, worker_id: i64) -> Result<Self, String> {
        let base_dir = base_dir.as_ref().to_path_buf();
        for radix in 0..RADIX_COUNT {
            fs::create_dir_all(base_dir.join(radix.to_string()))
                .map_err(|error| error.to_string())?;
        }
        Ok(Self {
            base_dir,
            worker_id,
            edges: vec![Vec::new(); RADIX_COUNT],
        })
    }

    fn write_msb(&mut self, edge: Edge) {
        let radix = (hash64(edge.member_oid) >> (64 - RADIX_BITS)) as usize;
        self.edges[radix].push(edge);
    }

    fn finish(&self) -> Result<RadixedTable, String> {
        let mut buckets = Vec::with_capacity(RADIX_COUNT);
        for radix in 0..RADIX_COUNT {
            let directory = self.base_dir.join(radix.to_string());
            let manifest = directory.join("bucket.tsv");
            let edges = &self.edges[radix];
            if edges.is_empty() {
                fs::write(&manifest, []).map_err(|error| error.to_string())?;
            } else {
                let volume = directory.join(format!("worker_{}_volume_0", self.worker_id));
                let mut compressed = CompressedBuffer::new();
                for edge in edges {
                    compressed
                        .write(&encode_edge(*edge))
                        .map_err(|error| error.to_string())?;
                }
                compressed.finish().map_err(|error| error.to_string())?;
                fs::write(&volume, compressed.data()).map_err(|error| error.to_string())?;
                fs::write(
                    &manifest,
                    format!("{}\t{}\n", volume.display(), edges.len()),
                )
                .map_err(|error| error.to_string())?;
            }
            buckets.push(Bucket::new(manifest, None));
        }
        Ok(RadixedTable(buckets))
    }
}

fn encode_edge(edge: Edge) -> [u8; 24] {
    let mut bytes = [0u8; 24];
    bytes[0..8].copy_from_slice(&edge.rep_oid.to_ne_bytes());
    bytes[8..16].copy_from_slice(&edge.member_oid.to_ne_bytes());
    bytes[16..20].copy_from_slice(&edge.rep_len.to_ne_bytes());
    bytes[20..24].copy_from_slice(&edge.member_len.to_ne_bytes());
    bytes
}

/// C++ file-local `align_rep` using the shared full-matrix SWIPE engine.
pub fn align_rep(
    chunk_sequences: &ChunkSequences,
    pairs: &[PairEntryShort],
    score_matrix: &ScoreMatrix,
    config: &AlignConfig,
    output: &mut Vec<Edge>,
    alignment_output: Option<&mut dyn Write>,
    stats: &mut ClusterStats,
) -> Result<(), String> {
    let Some(first) = pairs.first() else {
        return Ok(());
    };
    if pairs.iter().any(|pair| pair.rep_oid != first.rep_oid) {
        return Err("Pair group contains multiple representative OIDs.".to_string());
    }
    let representative = chunk_sequences.get(first.rep_oid)?;
    let mut dp_targets = targets();
    for (index, pair) in pairs.iter().enumerate() {
        let member = chunk_sequences.get(pair.member_oid)?;
        let bin = swipe_bin(
            HspValues::COORDS,
            representative.len() as i32,
            0,
            0,
            representative.len() as i64 * member.len() as i64,
            0,
            0,
            config.cutoff_score_8bit,
            config.max_swipe_dp,
            config.approx_backtrace,
        );
        dp_targets[bin].push_back(DpTarget::full(
            member.to_vec(),
            member.len() as i32,
            index as i64,
            CarryOver::default(),
        ));
        stats.extensions_computed += 1;
    }
    let swipe_statistics = Arc::new(Mutex::new(Statistics::new()));
    let mut params = Params::new(representative, score_matrix);
    params.flags = Flags::FULL_MATRIX;
    params.v = HspValues::COORDS;
    params.cutoff_score_8bit = config.cutoff_score_8bit;
    params.max_swipe_dp = config.max_swipe_dp;
    params.approx_backtrace = config.approx_backtrace;
    params.max_evalue = config.max_evalue;
    params.statistics = Some(swipe_statistics);
    let hsps = swipe(&dp_targets, &mut params);
    stats.hits_evalue_filtered += hsps.len() as u64;
    let mut alignment_output = alignment_output;
    for hsp in hsps {
        let pair = pairs
            .get(hsp.swipe_target as usize)
            .ok_or_else(|| "Invalid SWIPE target index.".to_string())?;
        let member = chunk_sequences.get(pair.member_oid)?;
        let query_cover = hsp.query_cover_percent(representative.len() as u32);
        let subject_cover = hsp.subject_cover_percent(member.len() as u32);
        if let Some(writer) = alignment_output.as_deref_mut() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                pair.rep_oid, pair.member_oid, query_cover, subject_cover, hsp.evalue
            )
            .map_err(|error| error.to_string())?;
        }
        if hsp.approx_id_percent(representative, member) < config.approx_min_identity {
            continue;
        }
        let mut passed = false;
        if let Some(mutual_cover) = config.mutual_cover {
            if subject_cover >= mutual_cover && query_cover >= mutual_cover {
                let (rep_oid, member_oid, rep_len, member_len) = if pair.rep_oid <= pair.member_oid
                {
                    (
                        pair.rep_oid,
                        pair.member_oid,
                        representative.len() as u32,
                        member.len() as u32,
                    )
                } else {
                    (
                        pair.member_oid,
                        pair.rep_oid,
                        member.len() as u32,
                        representative.len() as u32,
                    )
                };
                output.push(Edge {
                    rep_oid,
                    member_oid,
                    rep_len,
                    member_len,
                });
                passed = true;
            }
        } else {
            if subject_cover >= config.member_cover {
                output.push(Edge {
                    rep_oid: pair.rep_oid,
                    member_oid: pair.member_oid,
                    rep_len: representative.len() as u32,
                    member_len: member.len() as u32,
                });
                passed = true;
            }
            if query_cover >= config.member_cover {
                output.push(Edge {
                    rep_oid: pair.member_oid,
                    member_oid: pair.rep_oid,
                    rep_len: member.len() as u32,
                    member_len: representative.len() as u32,
                });
                passed = true;
            }
        }
        stats.hits_filtered += u64::from(passed);
    }
    Ok(())
}

fn read_volumed_file(path: &Path) -> Result<VolumedFile, String> {
    let text = fs::read_to_string(path).map_err(|error| error.to_string())?;
    let mut volumes = Vec::new();
    for line in text.lines().filter(|line| !line.trim().is_empty()) {
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.is_empty() {
            continue;
        }
        let records = fields
            .last()
            .and_then(|value| value.parse::<u64>().ok())
            .unwrap_or(0);
        volumes.push(Volume::new(fields[0], 0, 0, records));
    }
    Ok(VolumedFile(volumes))
}

fn read_pair_batches(path: &Path) -> Result<Vec<Vec<PairEntryShort>>, String> {
    let bytes = fs::read(path).map_err(|error| error.to_string())?;
    let mut offset = 0usize;
    let mut batches = Vec::new();
    while offset < bytes.len() {
        if bytes.len() - offset < std::mem::size_of::<usize>() {
            return Err("Short pair batch size.".to_string());
        }
        let mut size_bytes = [0u8; std::mem::size_of::<usize>()];
        size_bytes.copy_from_slice(&bytes[offset..offset + std::mem::size_of::<usize>()]);
        offset += std::mem::size_of::<usize>();
        let size = usize::from_ne_bytes(size_bytes);
        let byte_count = size
            .checked_mul(16)
            .ok_or_else(|| "Pair batch size overflow.".to_string())?;
        if bytes.len() - offset < byte_count {
            return Err("Short PairEntryShort batch.".to_string());
        }
        let mut batch = Vec::with_capacity(size);
        for record in bytes[offset..offset + byte_count].chunks_exact(16) {
            batch.push(PairEntryShort {
                rep_oid: u64::from_ne_bytes(record[0..8].try_into().unwrap()),
                member_oid: u64::from_ne_bytes(record[8..16].try_into().unwrap()),
            });
        }
        offset += byte_count;
        batches.push(batch);
    }
    Ok(batches)
}

/// C++ `External::align`; cross-process work claiming maps to deterministic
/// chunk iteration while retaining the exact files, grouping and radix output.
pub fn align<R: SequenceVolumeReader>(
    job: &Job,
    chunk_count: usize,
    _db_size: i64,
    score_matrix: &ScoreMatrix,
    config: &AlignConfig,
    reader: &mut R,
) -> Result<(RadixedTable, ClusterStats), String> {
    if config.threads == 0 {
        return Err("Alignment thread count must be positive.".to_string());
    }
    let chunks_path = job.base_dir(None).join("chunks");
    let alignment_path = job.base_dir(None).join("alignments");
    let mut output_files = EdgeFileArray::new(&alignment_path, job.worker_id())?;
    let mut all_stats = ClusterStats::default();
    for chunk in 0..chunk_count {
        let chunk_path = chunks_path.join(chunk.to_string());
        let volumes = read_volumed_file(&chunk_path.join("bucket.tsv"))?;
        let chunk_sequences = ChunkSequences::from_volumes(&volumes, reader)?;
        job.log(&format!(
            "Computing alignments. Chunk={}/{} Volumes={} Sequences={} Letters={}",
            chunk + 1,
            chunk_count,
            chunk_sequences.volumes(),
            chunk_sequences.oids(),
            chunk_sequences.letters()
        ))?;
        let pairs_path = chunk_path.join("pairs");
        let batches = read_pair_batches(&pairs_path)?;
        let mut alignment_output = if let Some(path) = &config.alignment_output {
            Some(
                OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(path)
                    .map_err(|_| {
                        format!("Failed to open alignment output file: {}", path.display())
                    })?,
            )
        } else {
            None
        };
        for batch in batches {
            let mut begin = 0usize;
            while begin < batch.len() {
                let rep_oid = batch[begin].rep_oid;
                let mut end = begin + 1;
                while end < batch.len() && batch[end].rep_oid == rep_oid {
                    end += 1;
                }
                let mut edges = Vec::new();
                align_rep(
                    &chunk_sequences,
                    &batch[begin..end],
                    score_matrix,
                    config,
                    &mut edges,
                    alignment_output.as_mut().map(|file| file as &mut dyn Write),
                    &mut all_stats,
                )?;
                for edge in edges {
                    output_files.write_msb(edge);
                }
                begin = end;
            }
        }
        fs::remove_file(&pairs_path).map_err(|error| error.to_string())?;
    }
    job.log_stats(&all_stats)?;
    Ok((output_files.finish()?, all_stats))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::AMINO_ACID_ALPHABET;
    use crate::cluster::external::seed_table::SeedSequence;
    use crate::stats::score_matrix::ScoreMatrix;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

    fn sequence(text: &str) -> Vec<Letter> {
        text.bytes()
            .map(|byte| {
                AMINO_ACID_ALPHABET
                    .iter()
                    .position(|&value| value == byte)
                    .unwrap() as Letter
            })
            .collect()
    }

    #[test]
    fn align_rep_emits_both_unidirectional_edges_and_exact_stats() {
        let mut sequences = ChunkSequences::default();
        sequences.sequences.insert(4, sequence("ARNDCQEGH"));
        sequences.sequences.insert(9, sequence("ARNDCQEGH"));
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1_000_000_000).unwrap();
        let config = AlignConfig {
            member_cover: 80.0,
            ..AlignConfig::default()
        };
        let mut edges = Vec::new();
        let mut alignment = Vec::new();
        let mut stats = ClusterStats::default();
        align_rep(
            &sequences,
            &[PairEntryShort {
                rep_oid: 4,
                member_oid: 9,
            }],
            &matrix,
            &config,
            &mut edges,
            Some(&mut alignment),
            &mut stats,
        )
        .unwrap();
        assert_eq!(stats.extensions_computed, 1);
        assert_eq!(stats.hits_evalue_filtered, 1);
        assert_eq!(stats.hits_filtered, 1);
        assert_eq!(edges.len(), 2);
        assert_eq!((edges[0].rep_oid, edges[0].member_oid), (4, 9));
        assert_eq!((edges[1].rep_oid, edges[1].member_oid), (9, 4));
        let line = String::from_utf8(alignment).unwrap();
        let fields = line.trim_end().split('\t').collect::<Vec<_>>();
        assert_eq!(&fields[..2], &["4", "9"]);
        assert_eq!(fields.len(), 5);
        assert!(fields[2].parse::<f64>().unwrap() >= 80.0);
        assert!(fields[3].parse::<f64>().unwrap() >= 80.0);
    }

    #[test]
    fn align_rep_mutual_cover_canonicalizes_oid_order() {
        let mut sequences = ChunkSequences::default();
        sequences.sequences.insert(8, sequence("ARNDCQEGH"));
        sequences.sequences.insert(2, sequence("ARNDCQEGH"));
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1_000_000_000).unwrap();
        let config = AlignConfig {
            mutual_cover: Some(80.0),
            ..AlignConfig::default()
        };
        let mut edges = Vec::new();
        let mut stats = ClusterStats::default();
        align_rep(
            &sequences,
            &[PairEntryShort {
                rep_oid: 8,
                member_oid: 2,
            }],
            &matrix,
            &config,
            &mut edges,
            None,
            &mut stats,
        )
        .unwrap();
        assert_eq!(edges.len(), 1);
        assert_eq!((edges[0].rep_oid, edges[0].member_oid), (2, 8));
        assert_eq!(stats.hits_filtered, 1);
    }

    #[test]
    fn pair_batches_reject_truncated_native_records() {
        let path = std::env::temp_dir().join(format!(
            "diamond-align-pairs-{}-{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("test")
        ));
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&1usize.to_ne_bytes());
        bytes.extend_from_slice(&3u64.to_ne_bytes());
        fs::write(&path, bytes).unwrap();
        assert_eq!(
            read_pair_batches(&path).unwrap_err(),
            "Short PairEntryShort batch."
        );
        fs::remove_file(path).unwrap();
    }

    #[derive(Default)]
    struct FixtureReader;

    impl SequenceVolumeReader for FixtureReader {
        fn read_volume(&mut self, _path: &Path) -> Result<Vec<SeedSequence>, String> {
            Ok(vec![
                SeedSequence {
                    id: "4 representative".to_string(),
                    sequence: sequence("ARNDCQEGH"),
                },
                SeedSequence {
                    id: "9 member".to_string(),
                    sequence: sequence("ARNDCQEGH"),
                },
            ])
        }
    }

    #[test]
    fn align_consumes_chunk_pairs_and_emits_radix_manifests() {
        let root = std::env::temp_dir().join(format!(
            "diamond-external-align-{}-{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        let job = Job::new(10, 1, 1 << 20, &root, 3).unwrap();
        let chunk = job.base_dir(None).join("chunks/0");
        fs::create_dir_all(&chunk).unwrap();
        fs::write(chunk.join("bucket.tsv"), "fixture.fa\t2\n").unwrap();
        let mut pairs = Vec::new();
        pairs.extend_from_slice(&1usize.to_ne_bytes());
        pairs.extend_from_slice(&4u64.to_ne_bytes());
        pairs.extend_from_slice(&9u64.to_ne_bytes());
        fs::write(chunk.join("pairs"), pairs).unwrap();
        let matrix = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 1_000_000_000).unwrap();
        let (table, stats) = align(
            &job,
            1,
            10,
            &matrix,
            &AlignConfig::default(),
            &mut FixtureReader,
        )
        .unwrap();

        assert_eq!(table.len(), RADIX_COUNT);
        assert_eq!(stats.extensions_computed, 1);
        assert_eq!(stats.hits_filtered, 1);
        assert!(!chunk.join("pairs").exists());
        assert!(table.iter().all(|bucket| bucket.path.exists()));
        fs::remove_dir_all(root).unwrap();
    }
}
