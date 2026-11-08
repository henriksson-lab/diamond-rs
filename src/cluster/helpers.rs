//! Translation of `diamond/src/cluster/helpers.cpp`.

use std::fmt::Display;
use std::io::Write;
use std::path::Path;

use crate::basic::value::{BlockId, OId, SuperBlockId};
use crate::util::algo::Edge;
use crate::util::data_structures::{make_flat_array, FlatArray};
use crate::util::tsv::{Config as TsvConfig, File as TsvFile, Flags, Type};

pub const HEADER_LINE: &str = "centroid\tmember";
pub const DEFAULT_MEMBER_COVER: f64 = 80.0;

/// Database operations consumed by the upstream `SequenceFile` helpers.
pub trait ClusterHelperDatabase {
    fn sequence_count(&self) -> usize;
    fn accession_to_oid(&self, accession: &str) -> Result<Vec<OId>, String>;
    fn accession(&self, oid: OId) -> Result<&str, String>;
    fn accession_to_oid_bytes(&self, accession: &[u8]) -> Result<Vec<OId>, String> {
        self.accession_to_oid(std::str::from_utf8(accession).map_err(|error| error.to_string())?)
    }
    fn accession_bytes(&self, oid: OId) -> Result<&[u8], String> {
        Ok(self.accession(oid)?.as_bytes())
    }
}

fn first_oid<D: ClusterHelperDatabase>(
    db: &D,
    accession: &[u8],
    errors: Option<&crate::util::io::IoErrorOwner>,
) -> Result<OId, String> {
    if let Some(oid) = db.accession_to_oid_bytes(accession)?.first().copied() {
        return Ok(oid);
    }
    let mut message = b"Accession not found in database: ".to_vec();
    message.extend_from_slice(accession);
    let error = crate::util::io::RawDiagnosticError(message).into_io();
    if let Some(errors) = errors {
        errors.capture(Err(error))
    } else {
        Err(error.to_string())
    }
}

fn open_cluster_lines(
    path: &Path,
    simple_header: bool,
    missing_header: &str,
    buffer_size: usize,
    source_profile: bool,
    errors: Option<&crate::util::io::IoErrorOwner>,
) -> Result<
    (
        i64,
        crate::util::io::text_input_file::ByteTextInputFile<crate::util::io::InputFileStream>,
    ),
    String,
> {
    // Both Source read overloads count the decoded lines before opening the
    // parser, including before checking the optional header.
    let count = if source_profile {
        if let Some(errors) = errors {
            errors.capture(
                crate::util::tsv::tsv::count_lines_with_source_buffer_path_io(path, buffer_size),
            )?
        } else {
            crate::util::tsv::tsv::count_lines_with_source_buffer_path_io(path, buffer_size)
                .map_err(|error| error.to_string())?
        }
    } else {
        crate::util::tsv::tsv::count_lines_with_buffer_path(path, buffer_size)?
    };
    let input = if source_profile {
        crate::util::io::text_input_file::ByteTextInputFile::from_source_path_bytes_with_buffer_size(path, buffer_size)
    } else {
        crate::util::io::text_input_file::ByteTextInputFile::from_file_path_bytes_with_buffer_size(path, buffer_size)
    }
    .map_err(std::io::Error::other);
    let mut input = if let Some(errors) = errors {
        errors.capture(input)?
    } else {
        input.map_err(|error| error.to_string())?
    };
    if simple_header {
        let read = input.getline().map_err(std::io::Error::other);
        if let Some(errors) = errors {
            errors.capture(read)?;
        } else {
            read.map_err(|error| error.to_string())?;
        }
        if input.line != HEADER_LINE.as_bytes() {
            return Err(missing_header.to_owned());
        }
    }
    Ok((count, input))
}

fn visit_cluster_lines(
    input: &mut crate::util::io::text_input_file::ByteTextInputFile<
        crate::util::io::InputFileStream,
    >,
    errors: Option<&crate::util::io::IoErrorOwner>,
    mut visit: impl FnMut(&[u8], &[u8]) -> Result<(), String>,
) -> Result<usize, String> {
    let mut mappings = 0usize;
    loop {
        let loaded = input.getline().map_err(std::io::Error::other);
        if let Some(errors) = errors {
            errors.capture(loaded)?;
        } else {
            loaded.map_err(|error| error.to_string())?;
        }
        if input.eof() && input.line.is_empty() {
            break;
        }
        // This is Source's older C-string Tokenizer, not TSV CharTokenizer:
        // a delimiter at the end supplies an empty second token.
        let line = crate::util::sequence::c_string_bytes(&input.line);
        let mut fields = line.splitn(3, |&byte| byte == b'\t');
        let centroid = fields.next().unwrap_or(&[]);
        let member = fields
            .next()
            .ok_or_else(|| "Tokenizer Exception".to_owned())?;
        visit(centroid, member)?;
        mappings += 1;
        if mappings % 1_000_000 == 0 {
            crate::util::log_stream::log_stream()
                .lock()
                .map_err(|_| "Log stream lock poisoned".to_owned())?
                .write(format!("#Entries: {mappings}\n"))
                .map_err(|error| error.to_string())?;
        }
    }
    Ok(mappings)
}

/// `read(..., CentroidSorted)` overload.
pub fn read_centroid_sorted<Int, D>(
    file_name: &Path,
    db: &D,
    simple_header: bool,
) -> Result<(FlatArray<Int>, Vec<Int>), String>
where
    Int: Copy + Ord + std::hash::Hash + TryFrom<OId>,
    <Int as TryFrom<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    read_centroid_sorted_profile(file_name, db, simple_header, 8192, false, None)
}

pub fn read_centroid_sorted_with_source_buffer<Int, D>(
    file_name: &Path,
    db: &D,
    simple_header: bool,
    buffer_size: usize,
) -> Result<(FlatArray<Int>, Vec<Int>), String>
where
    Int: Copy + Ord + std::hash::Hash + TryFrom<OId>,
    <Int as TryFrom<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    read_centroid_sorted_profile(file_name, db, simple_header, buffer_size, true, None)
}

pub fn read_centroid_sorted_with_source_buffer_owned<Int, D>(
    file_name: &Path,
    db: &D,
    simple_header: bool,
    buffer_size: usize,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<(FlatArray<Int>, Vec<Int>), String>
where
    Int: Copy + Ord + std::hash::Hash + TryFrom<OId>,
    <Int as TryFrom<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    read_centroid_sorted_profile(
        file_name,
        db,
        simple_header,
        buffer_size,
        true,
        Some(errors),
    )
}

fn read_centroid_sorted_profile<Int, D>(
    file_name: &Path,
    db: &D,
    simple_header: bool,
    buffer_size: usize,
    source_profile: bool,
    errors: Option<&crate::util::io::IoErrorOwner>,
) -> Result<(FlatArray<Int>, Vec<Int>), String>
where
    Int: Copy + Ord + std::hash::Hash + TryFrom<OId>,
    <Int as TryFrom<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    let (count, mut input) = open_cluster_lines(
        file_name,
        simple_header,
        "Clusters file is missing header line.",
        buffer_size,
        source_profile,
        errors,
    )?;
    let mut pairs = Vec::new();
    pairs
        .try_reserve(usize::try_from(count).map_err(|error| error.to_string())?)
        .map_err(|error| error.to_string())?;
    visit_cluster_lines(&mut input, errors, |centroid, member| {
        pairs.push((
            Int::try_from(first_oid(db, centroid, errors)?)
                .map_err(|error| format!("{error:?}"))?,
            Int::try_from(first_oid(db, member, errors)?).map_err(|error| format!("{error:?}"))?,
        ));
        Ok(())
    })?;
    let closed = input.into_inner().close().map_err(std::io::Error::other);
    if let Some(errors) = errors {
        errors.capture(closed)?;
    } else {
        closed.map_err(|error| error.to_string())?;
    }
    Ok(make_flat_array(&mut pairs))
}

/// Mapping-vector `read` overload.
pub fn read_mapping<Int, D>(
    file_name: &Path,
    db: &D,
    simple_header: bool,
) -> Result<Vec<Int>, String>
where
    Int: Copy + Default + TryFrom<OId>,
    <Int as TryFrom<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    read_mapping_with_buffer(file_name, db, simple_header, 8192)
}

pub fn read_mapping_with_buffer<Int, D>(
    file_name: &Path,
    db: &D,
    simple_header: bool,
    buffer_size: usize,
) -> Result<Vec<Int>, String>
where
    Int: Copy + Default + TryFrom<OId>,
    <Int as TryFrom<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    read_mapping_profile(file_name, db, simple_header, buffer_size, false, None)
}

pub fn read_mapping_with_source_buffer<Int, D>(
    file_name: &Path,
    db: &D,
    simple_header: bool,
    buffer_size: usize,
) -> Result<Vec<Int>, String>
where
    Int: Copy + Default + TryFrom<OId>,
    <Int as TryFrom<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    read_mapping_profile(file_name, db, simple_header, buffer_size, true, None)
}

pub fn read_mapping_with_source_buffer_owned<Int, D>(
    file_name: &Path,
    db: &D,
    simple_header: bool,
    buffer_size: usize,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<Vec<Int>, String>
where
    Int: Copy + Default + TryFrom<OId>,
    <Int as TryFrom<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    read_mapping_profile(
        file_name,
        db,
        simple_header,
        buffer_size,
        true,
        Some(errors),
    )
}

fn read_mapping_profile<Int, D>(
    file_name: &Path,
    db: &D,
    simple_header: bool,
    buffer_size: usize,
    source_profile: bool,
    errors: Option<&crate::util::io::IoErrorOwner>,
) -> Result<Vec<Int>, String>
where
    Int: Copy + Default + TryFrom<OId>,
    <Int as TryFrom<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    let (_, mut input) = open_cluster_lines(
        file_name,
        simple_header,
        "Clustering input file is missing header line.",
        buffer_size,
        source_profile,
        errors,
    )?;
    let mut mapping = vec![Int::default(); db.sequence_count()];
    let mappings = visit_cluster_lines(&mut input, errors, |centroid, member| {
        let centroid = Int::try_from(first_oid(db, centroid, errors)?)
            .map_err(|error| format!("{error:?}"))?;
        let member = first_oid(db, member, errors)? as usize;
        let entry = mapping
            .get_mut(member)
            .ok_or_else(|| "Clustering member OID is outside the database".to_owned())?;
        *entry = centroid;
        Ok(())
    })?;
    let closed = input.into_inner().close().map_err(std::io::Error::other);
    if let Some(errors) = errors {
        errors.capture(closed)?;
    } else {
        closed.map_err(|error| error.to_string())?;
    }
    if mappings != db.sequence_count() {
        return Err("Invalid/incomplete clustering.".to_owned());
    }
    Ok(mapping)
}

pub fn member2centroid_mapping<Int>(clusters: &FlatArray<Int>, centroids: &[Int]) -> Vec<Int>
where
    Int: Copy + Default + TryInto<usize>,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
{
    let data_size: usize = clusters.data_size().try_into().unwrap();
    let mut mapping = vec![Int::default(); data_size];
    for (i, &centroid) in centroids.iter().enumerate() {
        for &member in clusters.range(i as u64) {
            mapping[member.try_into().unwrap()] = centroid;
        }
    }
    mapping
}

pub fn cluster_sorted<Int>(mapping: &[Int]) -> (FlatArray<Int>, Vec<Int>)
where
    Int: Copy + Ord + std::hash::Hash + TryFrom<usize>,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut pairs = mapping
        .iter()
        .copied()
        .enumerate()
        .map(|(member, centroid)| (centroid, Int::try_from(member).unwrap()))
        .collect::<Vec<_>>();
    make_flat_array(&mut pairs)
}

/// Join-style output path. Input records are `(centroid_oid, member_oid)` and
/// output is ordered by member OID, matching the two C++ sort/join passes.
pub fn output<Int, D>(out: &mut TsvFile, db: &D, oid_pairs: &[(Int, Int)]) -> Result<(), String>
where
    Int: Copy + Ord + TryInto<OId>,
    <Int as TryInto<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    let mut pairs = oid_pairs.to_vec();
    pairs.sort_by_key(|&(_, member)| member);
    for (centroid, member) in pairs {
        let centroid =
            db.accession_bytes(centroid.try_into().map_err(|error| format!("{error:?}"))?)?;
        let member =
            db.accession_bytes(member.try_into().map_err(|error| format!("{error:?}"))?)?;
        out.write_record_bytes(&[centroid, member])?;
    }
    Ok(())
}

pub fn output_mem_clusters<Int, D>(
    out: &mut TsvFile,
    db: &D,
    clusters: &FlatArray<Int>,
    centroids: &[Int],
    oid_output: bool,
) -> Result<(), String>
where
    Int: Copy + Display + TryInto<OId>,
    <Int as TryInto<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    for (i, &centroid_oid) in centroids.iter().enumerate() {
        for &member_oid in clusters.range(i as u64) {
            let centroid;
            let member;
            if oid_output {
                centroid = centroid_oid.to_string().into_bytes();
                member = member_oid.to_string().into_bytes();
            } else {
                centroid = db
                    .accession_bytes(
                        centroid_oid
                            .try_into()
                            .map_err(|error| format!("{error:?}"))?,
                    )?
                    .to_owned();
                member = db
                    .accession_bytes(
                        member_oid
                            .try_into()
                            .map_err(|error| format!("{error:?}"))?,
                    )?
                    .to_owned();
            }
            out.write_record_bytes(&[&centroid, &member])?;
        }
    }
    Ok(())
}

pub fn output_mem_mapping<Int, D>(
    out: &mut TsvFile,
    db: &D,
    mapping: &[Int],
    oid_output: bool,
) -> Result<(), String>
where
    Int: Copy + Ord + std::hash::Hash + Display + TryFrom<usize> + TryInto<OId>,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
    <Int as TryInto<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    let (clusters, centroids) = cluster_sorted(mapping);
    output_mem_clusters(out, db, &clusters, &centroids, oid_output)
}

/// Concrete Source output_mem file boundary. Retain the O(N) sorted numeric
/// clustering representation, but stage only one accession row as text.
pub fn output_mem_mapping_stream<Int, D, W>(
    out: &mut W,
    db: &D,
    mapping: &[Int],
    oid_output: bool,
) -> Result<(), String>
where
    Int: Copy + Ord + std::hash::Hash + Display + TryFrom<usize> + TryInto<OId>,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
    <Int as TryInto<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
    W: Write + ?Sized,
{
    let (clusters, centroids) = cluster_sorted(mapping);
    let mut line = Vec::new();
    for (index, &centroid_oid) in centroids.iter().enumerate() {
        let centroid = if oid_output {
            centroid_oid.to_string().into_bytes()
        } else {
            db.accession_bytes(
                centroid_oid
                    .try_into()
                    .map_err(|error| format!("{error:?}"))?,
            )?
            .to_owned()
        };
        for &member_oid in clusters.range(index as u64) {
            let member = if oid_output {
                member_oid.to_string().into_bytes()
            } else {
                db.accession_bytes(
                    member_oid
                        .try_into()
                        .map_err(|error| format!("{error:?}"))?,
                )?
                .to_owned()
            };
            line.clear();
            line.extend_from_slice(&centroid);
            line.push(b'\t');
            line.extend_from_slice(&member);
            line.push(b'\n');
            out.write_all(&line).map_err(|error| error.to_string())?;
        }
    }
    Ok(())
}

pub fn output_mem_pairs<Int, D>(
    out: &mut TsvFile,
    db: &D,
    centroid_member_pairs: &mut [(Int, Int)],
    oid_output: bool,
) -> Result<(), String>
where
    Int: Copy + Ord + std::hash::Hash + Display + TryInto<OId>,
    <Int as TryInto<OId>>::Error: std::fmt::Debug,
    D: ClusterHelperDatabase,
{
    let (clusters, centroids) = make_flat_array(centroid_member_pairs);
    output_mem_clusters(out, db, &clusters, &centroids, oid_output)
}

/// `output_mem(File&, SequenceFile&, File&)` overload. The C++ dispatcher
/// selects 32- or 64-bit records from the database size; Rust decodes the
/// signed TSV storage then follows the same selection.
pub fn output_mem_file<D: ClusterHelperDatabase>(
    out: &mut TsvFile,
    db: &D,
    oid_to_centroid_oid: &mut TsvFile,
    oid_output: bool,
) -> Result<(), String> {
    let mut records = Vec::<(i64, i64)>::new();
    oid_to_centroid_oid.read_typed(&mut records)?;
    if db.sequence_count() > i32::MAX as usize {
        let mut pairs = records
            .into_iter()
            .map(|(centroid, member)| (centroid as u64, member as u64))
            .collect::<Vec<_>>();
        output_mem_pairs(out, db, &mut pairs, oid_output)
    } else {
        let mut pairs = records
            .into_iter()
            .map(|(centroid, member)| (centroid as u32, member as u32))
            .collect::<Vec<_>>();
        output_mem_pairs(out, db, &mut pairs, oid_output)
    }
}

pub fn split<Int>(mapping: &[Int]) -> (Vec<Int>, Vec<Int>)
where
    Int: Copy + PartialEq + TryFrom<usize>,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
{
    let mut centroids = Vec::new();
    let mut members = Vec::with_capacity(mapping.len());
    for (i, &centroid) in mapping.iter().enumerate() {
        let oid = Int::try_from(i).unwrap();
        if centroid == oid {
            centroids.push(oid);
        } else {
            members.push(oid);
        }
    }
    (centroids, members)
}

pub fn member_counts(mapping: &[SuperBlockId]) -> Vec<SuperBlockId> {
    let mut counts = vec![0; mapping.len()];
    for &centroid in mapping {
        counts[centroid as usize] += 1;
    }
    counts
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClusterCommand {
    DeepClust,
    LinClust,
    ClusterReassign,
    Other,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ThresholdConfig {
    pub member_cover: Option<f64>,
    pub mutual_cover: Option<f64>,
    pub approx_min_id: Option<f64>,
    pub soft_masking: Option<String>,
    pub masking: Option<String>,
    pub diag_filter_id: Option<f64>,
    pub diag_filter_cov: Option<f64>,
    pub command: ClusterCommand,
}

pub fn init_thresholds(config: &mut ThresholdConfig) -> Result<(), String> {
    if config.member_cover.is_some() && config.mutual_cover.is_some() {
        return Err("--member-cover and --mutual-cover are mutually exclusive.".to_owned());
    }
    if config.mutual_cover.is_none() {
        config.member_cover.get_or_insert(DEFAULT_MEMBER_COVER);
    }
    config.approx_min_id.get_or_insert(match config.command {
        ClusterCommand::DeepClust => 0.0,
        ClusterCommand::LinClust => 90.0,
        _ => 50.0,
    });
    config
        .soft_masking
        .get_or_insert_with(|| "tantan".to_owned());
    config.masking.get_or_insert_with(|| "0".to_owned());
    if config.approx_min_id.unwrap() < 90.0 || config.mutual_cover.is_some() {
        return Ok(());
    }
    config
        .diag_filter_id
        .get_or_insert(config.approx_min_id.unwrap() - 10.0);
    let member_cover = config.member_cover.unwrap();
    config
        .diag_filter_cov
        .get_or_insert(if member_cover > 50.0 {
            member_cover - 10.0
        } else {
            0.0
        });
    Ok(())
}

pub fn open_out_tsv(file_name: &str, simple_header: bool) -> Result<TsvFile, String> {
    let mut file = TsvFile::new(
        vec![Type::String, Type::String],
        file_name,
        Flags::WRITE,
        TsvConfig::default(),
    )?;
    if simple_header {
        file.write_record_strings(&["centroid", "member"])?;
    }
    Ok(file)
}

pub fn len_sorted_clust(edges: &FlatArray<Edge<SuperBlockId>, SuperBlockId>) -> Vec<BlockId> {
    let mut clustering = vec![BlockId::MAX; edges.size() as usize];
    for i in 0..edges.size() {
        if clustering[i as usize] != BlockId::MAX {
            continue;
        }
        clustering[i as usize] = i as BlockId;
        for edge in edges.range(i) {
            if clustering[edge.node2 as usize] == BlockId::MAX {
                clustering[edge.node2 as usize] = i as BlockId;
            }
        }
    }
    clustering
}

pub fn output_edges<D: ClusterHelperDatabase>(
    file_name: &Path,
    db: &D,
    edges: &[Edge<SuperBlockId>],
) -> Result<(), String> {
    let mut out = std::fs::File::create(file_name).map_err(|error| error.to_string())?;
    for edge in edges {
        let node1 = db.accession_bytes(edge.node1 as OId)?;
        let node2 = db.accession_bytes(edge.node2 as OId)?;
        out.write_all(node1)
            .and_then(|()| out.write_all(b"\t"))
            .and_then(|()| out.write_all(node2))
            .and_then(|()| out.write_all(b"\n"))
            .map_err(|error| error.to_string())?;
    }
    Ok(())
}

pub fn round_value(
    parameters: &[String],
    name: &str,
    round: usize,
    round_count: usize,
) -> Result<f64, String> {
    if parameters.is_empty() || round_count == 0 || round >= round_count - 1 {
        return Ok(0.0);
    }
    if parameters.len() >= round_count {
        return Err(format!("Too many values provided for {name}"));
    }
    let mut values = Vec::with_capacity(round_count - 1);
    for parameter in parameters {
        let value = parameter
            .trim_start()
            .parse::<f64>()
            .map_err(|_| format!("Invalid value provided for {name}: {parameter}"))?;
        values.push(value);
    }
    values.splice(
        0..0,
        std::iter::repeat_n(values[0], round_count - 1 - values.len()),
    );
    Ok(values[round])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    struct Db {
        names: Vec<String>,
        ids: HashMap<String, OId>,
    }

    impl Db {
        fn new(names: &[&str]) -> Self {
            Self {
                names: names.iter().map(|s| (*s).to_owned()).collect(),
                ids: names
                    .iter()
                    .enumerate()
                    .map(|(i, s)| ((*s).to_owned(), i as OId))
                    .collect(),
            }
        }
    }

    impl ClusterHelperDatabase for Db {
        fn sequence_count(&self) -> usize {
            self.names.len()
        }
        fn accession_to_oid(&self, accession: &str) -> Result<Vec<OId>, String> {
            self.ids
                .get(accession)
                .copied()
                .map(|x| vec![x])
                .ok_or_else(|| "missing".to_owned())
        }
        fn accession(&self, oid: OId) -> Result<&str, String> {
            self.names
                .get(oid as usize)
                .map(String::as_str)
                .ok_or_else(|| "missing".to_owned())
        }
    }

    fn temp_path(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("diamond-rs-helpers-{}-{name}", std::process::id()))
    }

    #[test]
    fn read_overloads_and_mapping_helpers() {
        let db = Db::new(&["a", "b", "c", "d"]);
        let path = temp_path("read.tsv");
        std::fs::write(&path, "centroid\tmember\na\ta\na\tb\nc\tc\nc\td\n").unwrap();
        let mapping = read_mapping::<u32, _>(&path, &db, true).unwrap();
        assert_eq!(mapping, vec![0, 0, 2, 2]);
        let (clusters, centroids) = read_centroid_sorted::<u32, _>(&path, &db, true).unwrap();
        assert_eq!(centroids, vec![0, 2]);
        assert_eq!(clusters.range(0), &[0, 1]);
        assert_eq!(member2centroid_mapping(&clusters, &centroids), mapping);
        assert_eq!(split(&mapping), (vec![0, 2], vec![1, 3]));
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn output_overloads_preserve_their_distinct_orders() {
        let db = Db::new(&["a", "b", "c", "d"]);
        let mut joined = open_out_tsv("", false).unwrap();
        output(&mut joined, &db, &[(2u32, 3), (0, 1), (0, 0), (2, 2)]).unwrap();
        let got = (0..joined.table().size())
            .map(|i| {
                let r = joined.table().record(i);
                (r.get(0), r.get(1))
            })
            .collect::<Vec<_>>();
        assert_eq!(
            got,
            vec![
                ("a".into(), "a".into()),
                ("a".into(), "b".into()),
                ("c".into(), "c".into()),
                ("c".into(), "d".into())
            ]
        );

        let mut memory = open_out_tsv("", true).unwrap();
        output_mem_mapping(&mut memory, &db, &[0u32, 0, 2, 2], true).unwrap();
        assert_eq!(memory.table().size(), 5);
        assert_eq!(memory.table().record(0).get(0), "centroid");
        assert_eq!(memory.table().record(1).get(0), "0");
    }

    #[test]
    fn thresholds_rounds_counts_and_length_clustering() {
        let mut cfg = ThresholdConfig {
            member_cover: None,
            mutual_cover: None,
            approx_min_id: None,
            soft_masking: None,
            masking: None,
            diag_filter_id: None,
            diag_filter_cov: None,
            command: ClusterCommand::LinClust,
        };
        init_thresholds(&mut cfg).unwrap();
        assert_eq!(cfg.member_cover, Some(80.0));
        assert_eq!(cfg.approx_min_id, Some(90.0));
        assert_eq!(cfg.diag_filter_id, Some(80.0));
        assert_eq!(cfg.diag_filter_cov, Some(70.0));
        assert_eq!(member_counts(&[0, 0, 2, 2]), vec![2, 0, 2, 0]);

        let p = vec!["90".to_owned(), "80".to_owned()];
        assert_eq!(round_value(&p, "--x", 0, 4).unwrap(), 90.0);
        assert_eq!(round_value(&p, "--x", 1, 4).unwrap(), 90.0);
        assert_eq!(round_value(&p, "--x", 2, 4).unwrap(), 80.0);
        assert_eq!(round_value(&p, "--x", 3, 4).unwrap(), 0.0);

        let edges = FlatArray::from_limits_data(
            vec![0u32, 2, 2, 3, 3],
            vec![
                Edge::new(0, 1, 1.0),
                Edge::new(0, 2, 1.0),
                Edge::new(2, 3, 1.0),
            ],
        );
        // Row 2 is not traversed because node 2 was already claimed by row 0.
        assert_eq!(len_sorted_clust(&edges), vec![0, 0, 0, 3]);
    }

    #[test]
    fn exact_errors_file_overload_and_edge_output() {
        let db = Db::new(&["a", "b"]);
        let bad_header = temp_path("bad-header.tsv");
        std::fs::write(&bad_header, "wrong\na\ta\n").unwrap();
        assert_eq!(
            read_mapping::<u32, _>(&bad_header, &db, true).unwrap_err(),
            "Clustering input file is missing header line."
        );
        std::fs::remove_file(bad_header).unwrap();

        assert_eq!(
            round_value(&["x".to_owned()], "--round-id", 0, 2).unwrap_err(),
            "Invalid value provided for --round-id: x"
        );
        assert_eq!(
            round_value(&["1".into(), "2".into()], "--x", 2, 3).unwrap(),
            0.0
        );
        assert_eq!(
            round_value(&["1".into(), "2".into()], "--x", 0, 2).unwrap_err(),
            "Too many values provided for --x"
        );

        let mut records =
            TsvFile::from_lines(vec![Type::Int64, Type::Int64], ["0\t0", "0\t1"]).unwrap();
        let mut out = open_out_tsv("", false).unwrap();
        output_mem_file(&mut out, &db, &mut records, false).unwrap();
        assert_eq!(out.table().record(1).get(1), "b");

        let edge_path = temp_path("edges.tsv");
        output_edges(&edge_path, &db, &[Edge::new(0, 1, 3.5)]).unwrap();
        assert_eq!(std::fs::read_to_string(&edge_path).unwrap(), "a\tb\n");
        std::fs::remove_file(edge_path).unwrap();
    }
}
