//! Pair-table construction from `diamond/src/cluster/external/pair_table.cpp`
//! and its inline `build_pair_table.h` algorithms.

use std::cmp::Ordering;
use std::fs;
use std::ops::{Deref, DerefMut};
use std::path::{Path, PathBuf};

use super::external::{Job, PairEntry};
use crate::util::hash::hash64;
use crate::util::io::{CompressedBuffer, InputFile};

pub const RADIX_BITS: u32 = 8;
pub const RADIX_COUNT: usize = 1 << RADIX_BITS;

/// Packed C++ `SeedEntry` (8-byte seed, 8-byte signed OID, 4-byte length).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SeedEntry {
    pub seed: u64,
    pub oid: i64,
    pub len: i32,
}

impl SeedEntry {
    pub const ENCODED_SIZE: usize = 20;

    pub const fn new(seed: u64, oid: i64, len: i32) -> Self {
        Self { seed, oid, len }
    }

    pub fn key(self) -> u64 {
        hash64(self.seed)
    }

    /// C++ `SeedEntry::Key::operator()` used for equal-seed merging.
    pub const fn merge_key(self) -> u64 {
        self.seed
    }

    pub fn serialize(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.seed.to_ne_bytes());
        out.extend_from_slice(&self.oid.to_ne_bytes());
        out.extend_from_slice(&self.len.to_ne_bytes());
    }

    pub fn deserialize(bytes: &[u8]) -> Result<Self, String> {
        if bytes.len() < Self::ENCODED_SIZE {
            return Err("Short SeedEntry".to_owned());
        }
        Ok(Self {
            seed: u64::from_ne_bytes(bytes[0..8].try_into().unwrap()),
            oid: i64::from_ne_bytes(bytes[8..16].try_into().unwrap()),
            len: i32::from_ne_bytes(bytes[16..20].try_into().unwrap()),
        })
    }
}

impl Ord for SeedEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.seed
            .cmp(&other.seed)
            .then_with(|| other.len.cmp(&self.len))
            .then_with(|| self.oid.cmp(&other.oid))
    }
}

impl PartialOrd for SeedEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bucket {
    pub path: PathBuf,
    pub records: Option<u64>,
}

impl Bucket {
    pub fn new(path: impl Into<PathBuf>, records: Option<u64>) -> Self {
        Self {
            path: path.into(),
            records,
        }
    }

    pub fn record_count(&self) -> Result<u64, String> {
        self.records
            .ok_or_else(|| format!("Record count not set for bucket: {}", self.path.display()))
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RadixedTable(pub Vec<Bucket>);

impl Deref for RadixedTable {
    type Target = [Bucket];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for RadixedTable {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl RadixedTable {
    pub fn from_manifest(path: impl AsRef<Path>) -> Result<Self, String> {
        let text = fs::read_to_string(path).map_err(|error| error.to_string())?;
        let mut buckets = Vec::new();
        let mut fields = text.split_whitespace();
        while let Some(path) = fields.next() {
            let records = fields
                .next()
                .ok_or_else(|| "Format error in RadixedTable".to_owned())?
                .parse::<u64>()
                .map_err(|_| "Format error in RadixedTable".to_owned())?;
            buckets.push(Bucket::new(path, Some(records)));
        }
        Ok(Self(buckets))
    }

    pub fn max_buckets(&self, mem_limit: u64, record_size: usize) -> Result<u64, String> {
        let mut counts = self
            .iter()
            .map(Bucket::record_count)
            .collect::<Result<Vec<_>, _>>()?;
        counts.sort_unstable_by(|a, b| b.cmp(a));
        let mut sum = 0_u64;
        for (index, count) in counts.iter().enumerate() {
            sum = sum.saturating_add(count.saturating_mul(record_size as u64));
            if sum >= mem_limit {
                return Ok(if index > 0 { index as u64 } else { 1 });
            }
        }
        Ok(counts.len() as u64)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairTableConfig {
    pub threads: usize,
    pub mutual_cover_present: bool,
    pub min_length_ratio: f64,
}

impl Default for PairTableConfig {
    fn default() -> Self {
        Self {
            threads: 1,
            mutual_cover_present: false,
            min_length_ratio: 0.0,
        }
    }
}

/// FileArray-compatible output owner. Each bucket is emitted as a zlib volume
/// plus the upstream tab-separated `bucket.tsv` volume manifest.
#[derive(Debug)]
pub struct PairFileArray {
    base_dir: PathBuf,
    worker_id: i64,
    pairs: Vec<Vec<PairEntry>>,
}

impl PairFileArray {
    pub fn new(base_dir: impl AsRef<Path>, worker_id: i64) -> Result<Self, String> {
        let base_dir = base_dir.as_ref().to_path_buf();
        for radix in 0..RADIX_COUNT {
            fs::create_dir_all(base_dir.join(radix.to_string()))
                .map_err(|error| error.to_string())?;
        }
        Ok(Self {
            base_dir,
            worker_id,
            pairs: vec![Vec::new(); RADIX_COUNT],
        })
    }

    pub fn write_msb(&mut self, entry: PairEntry) {
        let radix = (hash64(entry.rep_oid) >> (64 - RADIX_BITS)) as usize;
        self.pairs[radix].push(entry);
    }

    pub fn entries(&self, radix: usize) -> &[PairEntry] {
        &self.pairs[radix]
    }

    pub fn finish(&self) -> Result<RadixedTable, String> {
        let mut buckets = Vec::with_capacity(RADIX_COUNT);
        for radix in 0..RADIX_COUNT {
            let directory = self.base_dir.join(radix.to_string());
            let manifest = directory.join("bucket.tsv");
            let entries = &self.pairs[radix];
            if entries.is_empty() {
                fs::write(&manifest, []).map_err(|error| error.to_string())?;
            } else {
                let volume = directory.join(format!("worker_{}_volume_0", self.worker_id));
                let mut compressed = CompressedBuffer::new();
                for entry in entries {
                    compressed
                        .write(&encode_pair_entry(*entry))
                        .map_err(|error| error.to_string())?;
                }
                compressed.finish().map_err(|error| error.to_string())?;
                fs::write(&volume, compressed.data()).map_err(|error| error.to_string())?;
                fs::write(
                    &manifest,
                    format!("{}\t{}\n", volume.display(), entries.len()),
                )
                .map_err(|error| error.to_string())?;
            }
            buckets.push(Bucket::new(manifest, None));
        }
        Ok(RadixedTable(buckets))
    }
}

fn encode_pair_entry(entry: PairEntry) -> [u8; 24] {
    let mut bytes = [0_u8; 24];
    bytes[0..8].copy_from_slice(&entry.rep_oid.to_ne_bytes());
    bytes[8..16].copy_from_slice(&entry.member_oid.to_ne_bytes());
    bytes[16..20].copy_from_slice(&entry.rep_len.to_ne_bytes());
    bytes[20..24].copy_from_slice(&entry.member_len.to_ne_bytes());
    bytes
}

pub fn get_pairs_uni_cov(entries: &[SeedEntry], output: &mut PairFileArray) {
    let Some(rep) = entries.first() else {
        return;
    };
    for member in &entries[1..] {
        if rep.oid != member.oid {
            output.write_msb(PairEntry {
                rep_oid: rep.oid as u64,
                member_oid: member.oid as u64,
                rep_len: rep.len as u32,
                member_len: member.len as u32,
            });
        }
    }
}

pub fn get_pairs_mutual_cov(
    entries: &[SeedEntry],
    min_length_ratio: f64,
    output: &mut PairFileArray,
) {
    let size = entries.len();
    let mut i = 0;
    let mut j = 0;
    while i < size {
        let query_len = entries[i].len;
        let mut j1 = j;
        while j1 < size {
            let target_len = entries[j1].len;
            if f64::from(target_len) / f64::from(query_len) < min_length_ratio {
                break;
            }
            j1 += 1;
        }
        let query_position = i + (j1 - j) / 2;
        let rep = entries[query_position];
        for member in &entries[j..j1] {
            if rep.oid != member.oid {
                output.write_msb(PairEntry {
                    rep_oid: rep.oid as u64,
                    member_oid: member.oid as u64,
                    rep_len: rep.len as u32,
                    member_len: member.len as u32,
                });
            }
        }
        j = j1;
        if j == size {
            break;
        }
        let target_len = entries[j].len;
        while i < size {
            let query_len = entries[i].len;
            if f64::from(target_len) / f64::from(query_len) >= min_length_ratio {
                break;
            }
            i += 1;
        }
    }
}

/// C++ `External::build_pair_table`.
pub fn build_pair_table(
    job: &Job,
    seed_table: &RadixedTable,
    shape: i32,
    max_oid: i64,
    output_files: &mut PairFileArray,
    config: PairTableConfig,
) -> Result<RadixedTable, String> {
    let max_concurrent = seed_table.max_buckets(job.mem_limit, SeedEntry::ENCODED_SIZE)?;
    let concurrent_buckets = usize::try_from(max_concurrent)
        .unwrap_or(usize::MAX)
        .min(config.threads)
        .max(1);
    let bucket_workers = config.threads.div_ceil(concurrent_buckets);
    job.log(&format!(
        "Building pair table. Concurrent buckets={concurrent_buckets} Workers per bucket={bucket_workers}"
    ))?;
    let _shift = bit_length(max_oid as u64).saturating_sub(RADIX_BITS);
    let _queue_path = job
        .base_dir(None)
        .join(format!("seed_table_{shape}"))
        .join("build_pair_table_queue");

    for (bucket_index, bucket) in seed_table.iter().enumerate() {
        let (mut entries, volume_paths) = read_seed_bucket(bucket)?;
        job.log(&format!(
            "Building pair table. Bucket={}/{} Records={} Size={}",
            bucket_index + 1,
            seed_table.len(),
            entries.len(),
            entries.len() * SeedEntry::ENCODED_SIZE
        ))?;
        entries.sort_unstable();
        let mut begin = 0;
        while begin < entries.len() {
            let seed = entries[begin].merge_key();
            let mut end = begin + 1;
            while end < entries.len() && entries[end].merge_key() == seed {
                end += 1;
            }
            if config.mutual_cover_present {
                get_pairs_mutual_cov(&entries[begin..end], config.min_length_ratio, output_files);
            } else {
                get_pairs_uni_cov(&entries[begin..end], output_files);
            }
            begin = end;
        }
        for path in volume_paths {
            fs::remove_file(path).map_err(|error| error.to_string())?;
        }
        fs::remove_file(&bucket.path).map_err(|error| error.to_string())?;
        if let Some(parent) = bucket.path.parent() {
            let _ = fs::remove_dir(parent);
        }
    }
    output_files.finish()
}

fn bit_length(value: u64) -> u32 {
    u64::BITS - value.leading_zeros()
}

fn read_seed_bucket(bucket: &Bucket) -> Result<(Vec<SeedEntry>, Vec<PathBuf>), String> {
    let manifest = fs::read_to_string(&bucket.path).map_err(|error| error.to_string())?;
    let mut entries = Vec::new();
    let mut paths = Vec::new();
    for line in manifest.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 2 {
            return Err("Format error in VolumedFile".to_owned());
        }
        let path = PathBuf::from(fields[0]);
        let count = fields[1]
            .parse::<usize>()
            .map_err(|_| "Format error in VolumedFile".to_owned())?;
        let mut input = InputFile::new(
            path.to_str()
                .ok_or_else(|| "Non-UTF-8 seed volume path".to_owned())?,
            0,
        )
        .map_err(|error| error.to_string())?;
        for _ in 0..count {
            entries.push(SeedEntry {
                seed: input.read_u64().map_err(|error| error.to_string())?,
                oid: input.read_i64().map_err(|error| error.to_string())?,
                len: input.read_i32().map_err(|error| error.to_string())?,
            });
        }
        paths.push(path);
    }
    Ok((entries, paths))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

    fn temp_dir(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "diamond_pair_table_{name}_{}_{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, AtomicOrdering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    fn all_pairs(output: &PairFileArray) -> Vec<PairEntry> {
        let mut pairs: Vec<_> = output.pairs.iter().flatten().copied().collect();
        pairs.sort_unstable();
        pairs
    }

    #[test]
    fn seed_layout_and_order_match_packed_header_type() {
        let entry = SeedEntry::new(0x0102_0304_0506_0708, -9, 1234);
        let mut bytes = Vec::new();
        entry.serialize(&mut bytes);
        assert_eq!(bytes.len(), SeedEntry::ENCODED_SIZE);
        assert_eq!(SeedEntry::deserialize(&bytes).unwrap(), entry);

        let mut entries = [
            SeedEntry::new(7, 2, 50),
            SeedEntry::new(7, 1, 100),
            SeedEntry::new(7, 0, 100),
            SeedEntry::new(6, 9, 1),
        ];
        entries.sort_unstable();
        assert_eq!(
            entries.iter().map(|entry| entry.oid).collect::<Vec<_>>(),
            [9, 0, 1, 2]
        );
    }

    #[test]
    fn unidirectional_uses_longest_lowest_oid_rep_and_skips_self() {
        let root = temp_dir("uni");
        let mut output = PairFileArray::new(&root, 0).unwrap();
        let entries = [
            SeedEntry::new(1, 4, 100),
            SeedEntry::new(1, 4, 90),
            SeedEntry::new(1, 8, 80),
        ];
        get_pairs_uni_cov(&entries, &mut output);
        assert_eq!(
            all_pairs(&output),
            [PairEntry {
                rep_oid: 4,
                member_oid: 8,
                rep_len: 100,
                member_len: 80,
            }]
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn mutual_cover_windows_choose_middle_representatives() {
        let root = temp_dir("mutual");
        let mut output = PairFileArray::new(&root, 0).unwrap();
        let entries = [
            SeedEntry::new(1, 0, 100),
            SeedEntry::new(1, 1, 90),
            SeedEntry::new(1, 2, 80),
            SeedEntry::new(1, 3, 40),
        ];
        get_pairs_mutual_cov(&entries, 0.8, &mut output);
        assert_eq!(
            all_pairs(&output),
            [
                PairEntry {
                    rep_oid: 1,
                    member_oid: 0,
                    rep_len: 90,
                    member_len: 100,
                },
                PairEntry {
                    rep_oid: 1,
                    member_oid: 2,
                    rep_len: 90,
                    member_len: 80,
                },
            ]
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn build_reads_compressed_seed_volume_routes_pairs_and_removes_input() {
        let root = temp_dir("build");
        let seed_dir = root.join("seed_bucket");
        fs::create_dir_all(&seed_dir).unwrap();
        let volume = seed_dir.join("volume0");
        let manifest = seed_dir.join("bucket.tsv");
        let seeds = [
            SeedEntry::new(11, 2, 80),
            SeedEntry::new(10, 3, 70),
            SeedEntry::new(10, 1, 90),
            SeedEntry::new(11, 0, 100),
        ];
        let mut compressed = CompressedBuffer::new();
        for seed in seeds {
            let mut bytes = Vec::new();
            seed.serialize(&mut bytes);
            compressed.write(&bytes).unwrap();
        }
        compressed.finish().unwrap();
        fs::write(&volume, compressed.data()).unwrap();
        fs::write(
            &manifest,
            format!("{}\t{}\n", volume.display(), seeds.len()),
        )
        .unwrap();

        let job = Job::new(3, 1, 1 << 20, root.join("job"), 5).unwrap();
        let table = RadixedTable(vec![Bucket::new(&manifest, Some(seeds.len() as u64))]);
        let pair_dir = root.join("pairs");
        let mut output = PairFileArray::new(&pair_dir, job.worker_id()).unwrap();
        let result =
            build_pair_table(&job, &table, 0, 3, &mut output, PairTableConfig::default()).unwrap();

        assert_eq!(result.len(), RADIX_COUNT);
        assert_eq!(
            all_pairs(&output),
            [
                PairEntry {
                    rep_oid: 0,
                    member_oid: 2,
                    rep_len: 100,
                    member_len: 80,
                },
                PairEntry {
                    rep_oid: 1,
                    member_oid: 3,
                    rep_len: 90,
                    member_len: 70,
                },
            ]
        );
        assert!(!volume.exists());
        assert!(!manifest.exists());
        for pair in all_pairs(&output) {
            let radix = (hash64(pair.rep_oid) >> 56) as usize;
            assert!(output.entries(radix).contains(&pair));
            let line = fs::read_to_string(&result[radix].path).unwrap();
            assert!(line.ends_with("\t1\n"));
            let pair_volume = line.split_whitespace().next().unwrap();
            let mut input = InputFile::new(pair_volume, 0).unwrap();
            assert_eq!(input.read_u64().unwrap(), pair.rep_oid);
            assert_eq!(input.read_u64().unwrap(), pair.member_oid);
            assert_eq!(input.read_u32().unwrap(), pair.rep_len);
            assert_eq!(input.read_u32().unwrap(), pair.member_len);
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn max_buckets_matches_memory_cutoff_and_missing_counts_error() {
        let table = RadixedTable(vec![
            Bucket::new("a", Some(10)),
            Bucket::new("b", Some(5)),
            Bucket::new("c", Some(1)),
        ]);
        assert_eq!(table.max_buckets(250, 20).unwrap(), 1);
        assert_eq!(table.max_buckets(201, 20).unwrap(), 1);
        assert_eq!(table.max_buckets(301, 20).unwrap(), 2);
        assert!(RadixedTable(vec![Bucket::new("x", None)])
            .max_buckets(1, 20)
            .unwrap_err()
            .contains("Record count not set"));
    }
}
