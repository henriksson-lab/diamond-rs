//! Command-level greedy vertex cover.
//!
//! This mirrors `diamond/src/tools/greedy_vertex_cover.cpp`. The graph
//! algorithm itself remains in [`crate::util::algo::greedy_vertex_cover`];
//! this module implements the source file's mapping/edge parsing, directional
//! coverage rules, graph construction, and output generation.

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

use crate::util::algo::{greedy_vertex_cover as compute_vertex_cover, Edge};
use crate::util::data_structures::FlatArray;

/// Default used by `Cluster::DEFAULT_MEMBER_COVER` in the C++ command.
pub const DEFAULT_MEMBER_COVER: f64 = 80.0;

/// Explicit replacement for the process-wide C++ `config` object.
#[derive(Debug, Clone, PartialEq)]
pub struct GreedyVertexCoverConfig<T = String> {
    pub database: PathBuf,
    pub edges: PathBuf,
    pub query_or_target_cover: f64,
    pub member_cover: Option<f64>,
    pub edge_format: T,
    pub symmetric: bool,
    pub strict_gvc: bool,
    pub no_gvc_reassign: bool,
    pub connected_component_depth: Vec<T>,
    pub centroid_out: Option<PathBuf>,
    pub output_file: Option<PathBuf>,
}

impl<T: AsRef<[u8]>> GreedyVertexCoverConfig<T> {
    fn coverage_cutoff(&self) -> f64 {
        let member_cover = self.member_cover.unwrap_or(DEFAULT_MEMBER_COVER);
        // Match `std::max(a, b)` including its left-operand NaN behavior.
        if self.query_or_target_cover < member_cover {
            member_cover
        } else {
            self.query_or_target_cover
        }
    }

    fn connected_component_depth(&self) -> u64 {
        self.connected_component_depth
            .first()
            .map_or(0, |value| c_atoi(value.as_ref()) as u64)
    }
}

/// Generated command output, retained independently of optional output files.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GreedyVertexCoverOutput<T = String> {
    pub assignments: Vec<(T, T)>,
    pub centroids: Vec<T>,
}

pub type ByteGreedyVertexCoverOutput = GreedyVertexCoverOutput<Vec<u8>>;
impl<T> GreedyVertexCoverOutput<T> {
    pub fn centroid_count(&self) -> usize {
        self.centroids.len()
    }
}
impl GreedyVertexCoverOutput {
    fn assignments_text(&self) -> String {
        let mut output = String::new();
        for (centroid, member) in &self.assignments {
            output.push_str(centroid);
            output.push('\t');
            output.push_str(member);
            output.push('\n');
        }
        output
    }

    fn centroids_text(&self) -> String {
        let mut output = self.centroids.join("\n");
        if !output.is_empty() {
            output.push('\n');
        }
        output
    }
}

/// Translate C `atoi`'s command-line behavior.
fn c_atoi(value: &[u8]) -> i32 {
    let mut value = value[..value
        .iter()
        .position(|&byte| byte == 0)
        .unwrap_or(value.len())]
        .to_vec();
    value.push(0);
    unsafe { libc::atoi(value.as_ptr().cast()) }
}

/// Source LineIterator yields interior empty rows and excludes the final
/// zero-length field after a trailing LF. Only mapping getline trims CR.
fn byte_lines(bytes: &[u8]) -> impl Iterator<Item = &[u8]> {
    let end = bytes.len() - usize::from(bytes.last() == Some(&b'\n'));
    bytes[..end]
        .split(|&byte| byte == b'\n')
        .filter(move |_| !bytes.is_empty())
}
fn mapping_id(line: &[u8]) -> std::io::Result<Vec<u8>> {
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    mapping_token_id(line)
}
fn mapping_token_id(line: &[u8]) -> std::io::Result<Vec<u8>> {
    crate::util::string::ByteTokenizer::new(line, b'\t')
        .read_bytes()
        .map_err(std::io::Error::other)
}
fn add_mapping(map: &mut HashMap<Vec<u8>, u64>, accession: Vec<u8>) {
    let next = map.len() as u64;
    map.entry(accession).or_insert(next);
}
fn read_edge(line: &[u8], triplets: bool) -> std::io::Result<(Vec<u8>, Vec<u8>, f64, f64, f64)> {
    let mut fields = crate::util::string::ByteTokenizer::new(line, b'\t');
    let query = fields.read_bytes().map_err(std::io::Error::other)?;
    let target = fields.read_bytes().map_err(std::io::Error::other)?;
    let (qcov, tcov) = if triplets {
        (0.0, 0.0)
    } else {
        (
            fields.read_f64().map_err(std::io::Error::other)?,
            fields.read_f64().map_err(std::io::Error::other)?,
        )
    };
    let evalue = fields.read_f64().map_err(std::io::Error::other)?;
    Ok((query, target, qcov, tcov, evalue))
}
fn lookup(
    map: &HashMap<Vec<u8>, u64>,
    accession: &[u8],
    role: &[u8],
    line: usize,
) -> std::io::Result<u64> {
    map.get(accession).copied().ok_or_else(|| {
        // Retained fallible library convenience diagnostic. The actual Source
        // worker aborts without printing this what() string.
        let mut message = b"Unknown ".to_vec();
        message.extend_from_slice(role);
        message.extend_from_slice(b" accession `");
        message.extend_from_slice(accession);
        message.extend_from_slice(format!("` on edge line {line}.").as_bytes());
        crate::util::io::RawDiagnosticError(message).into_io()
    })
}
#[derive(Clone, Copy)]
struct EdgePolicy {
    triplets: bool,
    symmetric: bool,
    coverage: f64,
}
impl EdgePolicy {
    fn from_config<T: AsRef<[u8]>>(config: &GreedyVertexCoverConfig<T>) -> Self {
        Self {
            triplets: config.edge_format.as_ref() == b"triplet",
            symmetric: config.symmetric,
            coverage: config.coverage_cutoff(),
        }
    }
}
fn append_edge(
    line: &[u8],
    line_number: usize,
    map: &HashMap<Vec<u8>, u64>,
    policy: EdgePolicy,
    edges: &mut Vec<Edge<u64>>,
) -> std::io::Result<()> {
    let triplets = policy.triplets;
    let (query, target, qcov, tcov, evalue) = read_edge(line, triplets)?;
    let coverage = policy.coverage;
    if triplets || tcov >= coverage || qcov >= coverage {
        let query = lookup(map, &query, b"query", line_number)?;
        let target = lookup(map, &target, b"target", line_number)?;
        if query == target {
            return Ok(());
        }
        if triplets {
            edges.push(Edge::new(target, query, evalue));
            if policy.symmetric {
                edges.push(Edge::new(query, target, evalue));
            }
        } else {
            if tcov >= coverage {
                edges.push(Edge::new(query, target, evalue));
            }
            if qcov >= coverage {
                edges.push(Edge::new(target, query, evalue));
            }
        }
    }
    Ok(())
}

fn make_edge_array(mut edges: Vec<Edge<u64>>, node_count: usize) -> FlatArray<Edge<u64>, u64> {
    edges.sort_unstable_by(|left, right| {
        left.node1
            .cmp(&right.node1)
            .then_with(|| left.node2.cmp(&right.node2))
    });
    let mut limits = Vec::with_capacity(node_count + 1);
    limits.push(0);
    let mut cursor = 0usize;
    for node in 0..node_count {
        while cursor < edges.len() && edges[cursor].node1 == node as u64 {
            cursor += 1;
        }
        limits.push(cursor as u64);
    }
    FlatArray::from_limits_data(limits, edges)
}

/// Run the source workflow on in-memory text. This is the deterministic core
/// used by the file-backed command and focused tests.
pub fn greedy_vertex_cover(
    mapping_text: &str,
    edges_text: &str,
    config: &GreedyVertexCoverConfig,
) -> Result<GreedyVertexCoverOutput, String> {
    checked_output(
        greedy_vertex_cover_bytes(mapping_text.as_bytes(), edges_text.as_bytes(), config)
            .map_err(|error| error.to_string())?,
    )
    .map_err(|error| error.to_string())
}

pub fn greedy_vertex_cover_bytes<T: AsRef<[u8]>>(
    mapping: &[u8],
    edge_input: &[u8],
    config: &GreedyVertexCoverConfig<T>,
) -> std::io::Result<ByteGreedyVertexCoverOutput> {
    validate_mode(config)?;
    let mut map = HashMap::new();
    for line in byte_lines(mapping) {
        add_mapping(&mut map, mapping_id(line)?);
    }
    let mut edges = Vec::new();
    let policy = EdgePolicy::from_config(config);
    for (line, bytes) in byte_lines(edge_input).enumerate() {
        append_edge(bytes, line + 1, &map, policy, &mut edges)?;
    }
    Ok(finish_graph(map, edges, config))
}
fn validate_mode<T: AsRef<[u8]>>(config: &GreedyVertexCoverConfig<T>) -> std::io::Result<()> {
    if config.edge_format.as_ref() != b"triplet" && config.symmetric {
        return Err(std::io::Error::other(
            "--symmetric requires triplet edge format",
        ));
    }
    Ok(())
}
fn checked_output(output: ByteGreedyVertexCoverOutput) -> std::io::Result<GreedyVertexCoverOutput> {
    Ok(GreedyVertexCoverOutput {
        assignments: output
            .assignments
            .into_iter()
            .map(|(a, b)| {
                Ok((
                    String::from_utf8(a).map_err(std::io::Error::other)?,
                    String::from_utf8(b).map_err(std::io::Error::other)?,
                ))
            })
            .collect::<std::io::Result<_>>()?,
        centroids: output
            .centroids
            .into_iter()
            .map(|value| String::from_utf8(value).map_err(std::io::Error::other))
            .collect::<std::io::Result<_>>()?,
    })
}
fn cover<T: AsRef<[u8]>>(
    edge_array: &mut FlatArray<Edge<u64>, u64>,
    config: &GreedyVertexCoverConfig<T>,
) -> Vec<u64> {
    compute_vertex_cover(
        edge_array,
        None,
        !config.strict_gvc,
        !config.no_gvc_reassign,
        config.connected_component_depth(),
    )
}
fn reverse_mapping(accession_to_oid: HashMap<Vec<u8>, u64>) -> Vec<Vec<u8>> {
    let mut accessions = vec![Vec::new(); accession_to_oid.len()];
    for (accession, oid) in accession_to_oid {
        accessions[oid as usize] = accession;
    }
    accessions
}
fn emit_assignments(
    accessions: &[Vec<u8>],
    representatives: &[u64],
    mut row: impl FnMut(&[u8], &[u8], bool),
) -> ByteGreedyVertexCoverOutput {
    let mut assignments = Vec::with_capacity(accessions.len());
    let mut centroids = Vec::new();
    for (member, &representative) in representatives.iter().enumerate() {
        let is_centroid = representative as usize == member;
        if is_centroid {
            centroids.push(accessions[member].clone());
        }
        let centroid = &accessions[representative as usize];
        let member = &accessions[member];
        row(centroid, member, is_centroid);
        assignments.push((centroid.clone(), member.clone()));
    }
    GreedyVertexCoverOutput {
        assignments,
        centroids,
    }
}
fn finish_graph<T: AsRef<[u8]>>(
    accession_to_oid: HashMap<Vec<u8>, u64>,
    edges: Vec<Edge<u64>>,
    config: &GreedyVertexCoverConfig<T>,
) -> ByteGreedyVertexCoverOutput {
    let mut edge_array = make_edge_array(edges, accession_to_oid.len());
    let representatives = cover(&mut edge_array, config);
    emit_assignments(
        &reverse_mapping(accession_to_oid),
        &representatives,
        |_, _, _| {},
    )
}
fn source_message(text: String) {
    if let Ok(mut stream) = crate::util::log_stream::message_stream().lock() {
        let _ = stream.write(text);
    }
}

/// File-backed adapter for [`greedy_vertex_cover`].
pub fn greedy_vertex_cover_from_files(
    config: &GreedyVertexCoverConfig,
) -> Result<GreedyVertexCoverOutput, String> {
    let mapping_text = fs::read_to_string(&config.database).map_err(|error| {
        format!(
            "Failed to read mapping file {}: {error}",
            config.database.display()
        )
    })?;
    let edges_text = fs::read_to_string(&config.edges).map_err(|error| {
        format!(
            "Failed to read edge file {}: {error}",
            config.edges.display()
        )
    })?;
    let output = greedy_vertex_cover(&mapping_text, &edges_text, config)?;
    if let Some(path) = &config.centroid_out {
        fs::write(path, output.centroids_text()).map_err(|error| {
            format!("Failed to write centroid file {}: {error}", path.display())
        })?;
    }
    if let Some(path) = &config.output_file {
        fs::write(path, output.assignments_text())
            .map_err(|error| format!("Failed to write output file {}: {error}", path.display()))?;
    }
    Ok(output)
}

/// Source command I/O, retaining the actual failed OS operation through the
/// typed error boundary and decoding the same gzip/zlib build profile.
/// The public String-returning in-memory/plain-file APIs remain independent.
pub fn greedy_vertex_cover_from_source_files<T: AsRef<[u8]>>(
    config: &GreedyVertexCoverConfig<T>,
    buffer_size: usize,
) -> std::io::Result<GreedyVertexCoverOutput> {
    checked_output(greedy_vertex_cover_from_source_files_bytes(
        config,
        buffer_size,
    )?)
}
pub fn greedy_vertex_cover_from_source_files_bytes<T: AsRef<[u8]>>(
    config: &GreedyVertexCoverConfig<T>,
    buffer_size: usize,
) -> std::io::Result<ByteGreedyVertexCoverOutput> {
    greedy_vertex_cover_from_source_files_with_threads(
        config,
        buffer_size,
        crate::util::system::hardware_concurrency() as usize,
    )
}

/// Actual Source command path: the mapping table is the retained ID owner;
/// edge text lives only in bounded producer/consumer buffers on each pass.
pub fn greedy_vertex_cover_from_source_files_with_threads<T: AsRef<[u8]>>(
    config: &GreedyVertexCoverConfig<T>,
    buffer_size: usize,
    threads: usize,
) -> std::io::Result<ByteGreedyVertexCoverOutput> {
    use std::io;
    use std::sync::{
        atomic::{AtomicU64, Ordering},
        Mutex,
    };
    validate_mode(config)?;
    source_message(format!(
        "Coverage cutoff: {}%\n",
        crate::output::format::format_default_float(config.coverage_cutoff())
    ));
    let mut timer = crate::util::log_stream::TaskTimer::with_msg("Reading mapping file", 1);
    let count = source_operation(
        crate::util::tsv::tsv::count_lines_with_source_buffer_path_io(
            &config.database,
            buffer_size,
        ),
    )?;
    let mut accession_to_oid =
        HashMap::<Vec<u8>, u64>::with_capacity(usize::try_from(count).map_err(io::Error::other)?);
    let mut mapping = source_operation(
        crate::util::io::ByteTextInputFile::from_source_path_bytes_with_buffer_size(
            &config.database,
            buffer_size,
        )
        .map_err(io::Error::other),
    )?;
    loop {
        source_operation(mapping.getline().map_err(io::Error::other))?;
        if mapping.eof() && mapping.line.is_empty() {
            break;
        }
        add_mapping(&mut accession_to_oid, mapping_token_id(&mapping.line)?);
    }
    source_operation(mapping.into_inner().close().map_err(io::Error::other))?;
    timer.finish();
    source_message(format!("#OIds: {}\n", accession_to_oid.len()));
    timer.go_string("Counting input lines");
    let policy = EdgePolicy::from_config(config);
    let triplets = policy.triplets;
    let coverage = policy.coverage;
    let symmetric = policy.symmetric;
    let lines = AtomicU64::new(0);
    read_edge_pass(&config.edges, buffer_size, threads, &|bytes| {
        let mut n = 0u64;
        for line in byte_lines(bytes) {
            if triplets {
                n += if symmetric { 2 } else { 1 };
            } else {
                let mut fields = crate::util::string::ByteTokenizer::new(line, b'\t');
                fields.read_bytes().map_err(io::Error::other)?;
                fields.read_bytes().map_err(io::Error::other)?;
                let qcov = fields.read_f64().map_err(io::Error::other)?;
                let tcov = fields.read_f64().map_err(io::Error::other)?;
                n += u64::from(qcov >= coverage) + u64::from(tcov >= coverage);
            }
        }
        lines.fetch_add(n, Ordering::Relaxed);
        Ok(())
    })?;
    timer.finish();
    source_message(format!("#Lines: {}\n", lines.load(Ordering::Relaxed)));
    timer.go_string("Allocating memory");
    let mut edges = Vec::<Edge<u64>>::new();
    edges
        .try_reserve_exact(
            usize::try_from(lines.load(Ordering::Relaxed)).map_err(io::Error::other)?,
        )
        .map_err(io::Error::other)?;
    timer.go_string("Reading input lines");
    let edges = Mutex::new(edges);
    read_edge_pass(&config.edges, buffer_size, threads, &|bytes| {
        let mut chunk_edges = Vec::new();
        for (i, line) in byte_lines(bytes).enumerate() {
            append_edge(line, i + 1, &accession_to_oid, policy, &mut chunk_edges)?;
        }
        edges
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .extend(chunk_edges);
        Ok(())
    })?;
    let edges = edges
        .into_inner()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    timer.finish();
    let _ = crate::util::system::log_rss();
    timer.go_string("Making flat array");
    let mut edge_array = make_edge_array(edges, accession_to_oid.len());
    timer.finish();
    let _ = crate::util::system::log_rss();
    let representatives = cover(&mut edge_array, config);
    timer.go_string("Building reverse mapping");
    let accessions = reverse_mapping(accession_to_oid);
    timer.go_string("Generating output");
    // Both exception-disabled ofstreams are constructed before the first row.
    let mut centroids = SourceOfstream::new(config.centroid_out.as_deref());
    let mut assignments = SourceOfstream::new(config.output_file.as_deref());
    let output = emit_assignments(
        &accessions,
        &representatives,
        |centroid, member, is_centroid| {
            if is_centroid {
                centroids.row(&[member]);
            }
            assignments.row(&[centroid, b"\t", member]);
        },
    );
    centroids.close();
    assignments.close();
    timer.finish();
    Ok(output)
}

/// Source std::ofstream: text-mode CRT writes, sticky failbit, no exceptions
/// or perror from an unsuccessful open/write/endl flush/close.
struct SourceOfstream {
    file: Option<crate::util::io::cstdio_file::OwnedCstdio>,
    failed: bool,
}
impl SourceOfstream {
    fn new(path: Option<&std::path::Path>) -> Self {
        let file = path
            .filter(|path| !path.as_os_str().is_empty())
            .and_then(|path| {
                crate::util::io::cstdio_file::OwnedCstdio::open_file_role_path(path, "w").ok()
            });
        Self {
            file,
            failed: false,
        }
    }
    fn row(&mut self, fields: &[&[u8]]) {
        if self.failed {
            return;
        }
        let Some(file) = self.file.as_mut() else {
            return;
        };
        for field in fields
            .iter()
            .copied()
            .chain(std::iter::once(b"\n".as_slice()))
        {
            if !file.write_text_output(field).unwrap_or(false) {
                self.failed = true;
                return;
            }
        }
        self.failed = !file.flush_text_output().unwrap_or(false);
    }
    fn close(&mut self) {
        if let Some(mut file) = self.file.take() {
            let _ = file.close_file_role();
        }
    }
}
impl Drop for SourceOfstream {
    fn drop(&mut self) {
        self.close();
    }
}

/// Source already emitted perror at the low-level operation before a local
/// TaskTimer could unwind/flush cout. Keep that point and mark the retained
/// typed cause so the outer CLI reporter does not emit it a second time.
fn source_operation<T>(result: std::io::Result<T>) -> std::io::Result<T> {
    result.map_err(|error| match crate::util::io::emit_source_perror(&error) {
        Ok(true) => crate::util::io::PerrorEmittedError::wrap(error),
        Ok(false) => error,
        Err(conversion) => std::io::Error::other(conversion),
    })
}

fn abort_source_worker(error: &std::io::Error) -> ! {
    // Actual Source low-level I/O has already called perror before throwing.
    // Tokenizer/out_of_range exceptions produce no text: main's installed
    // terminate handler is just std::abort(), not the library's default one.
    if let Err(codec) = crate::util::io::emit_source_perror(error) {
        crate::util::io::abort_source_diagnostic_error(&codec);
    }
    crate::util::io::abort_source_process()
}

fn read_edge_pass<F>(
    path: &std::path::Path,
    buffer_size: usize,
    threads: usize,
    callback: &F,
) -> std::io::Result<()>
where
    F: Fn(&[u8]) -> std::io::Result<()> + Sync,
{
    use std::io;
    use std::sync::{mpsc, Mutex};
    let mut input = source_operation(
        crate::util::io::InputFile::with_source_profile_path(path, 0, buffer_size)
            .map_err(io::Error::other),
    )?;
    let consumers = threads.saturating_sub(1).max(1);
    let (send, receive) = mpsc::sync_channel::<Vec<u8>>(consumers);
    let receive = Mutex::new(receive);
    std::thread::scope(|scope| {
        let input = &mut input;
        scope.spawn(move || loop {
            let mut bytes = vec![0; crate::util::tsv::read_text_mt::READ_TEXT_MT_SIZE];
            let n = input
                .read_raw(&mut bytes)
                .map_err(io::Error::other)
                .unwrap_or_else(|error| abort_source_worker(&error));
            bytes.truncate(n);
            if n == crate::util::tsv::read_text_mt::READ_TEXT_MT_SIZE {
                input
                    .read_to(&mut bytes, b'\n')
                    .map_err(io::Error::other)
                    .unwrap_or_else(|error| abort_source_worker(&error));
            }
            let n = bytes.len();
            if n > 0 && send.send(bytes).is_err() {
                crate::util::io::abort_source_process();
            }
            if n < crate::util::tsv::read_text_mt::READ_TEXT_MT_SIZE {
                break;
            }
        });
        for _ in 0..consumers {
            let receive = &receive;
            scope.spawn(move || loop {
                let bytes = receive
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner())
                    .recv();
                let Ok(bytes) = bytes else {
                    break;
                };
                if let Err(error) = callback(&bytes) {
                    abort_source_worker(&error);
                }
            });
        }
        // The sender must be owned by the reader so consumers observe its
        // destruction after the final buffer, as Source's stop notification.
    });
    source_operation(input.close().map_err(io::Error::other))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> GreedyVertexCoverConfig {
        GreedyVertexCoverConfig {
            database: PathBuf::new(),
            edges: PathBuf::new(),
            query_or_target_cover: 0.0,
            member_cover: None,
            edge_format: "default".to_owned(),
            symmetric: false,
            strict_gvc: false,
            no_gvc_reassign: false,
            connected_component_depth: Vec::new(),
            centroid_out: None,
            output_file: None,
        }
    }

    #[test]
    fn byte_core_keeps_source_cstring_ids_and_crt_hex_scores() {
        let config = config();
        let output = greedy_vertex_cover_bytes(
            b"A\xff\0discarded\tignored\nA\xff\nB\x80\r\n",
            b"A\xff\tB\x80\t0x1.68p+6\t10\t1e-5\0ignored\n",
            &config,
        )
        .unwrap();
        assert_eq!(output.centroids, vec![b"B\x80".to_vec()]);
        assert_eq!(
            output.assignments,
            vec![
                (b"B\x80".to_vec(), b"A\xff".to_vec()),
                (b"B\x80".to_vec(), b"B\x80".to_vec())
            ]
        );
        // The edge LineIterator preserves CR, unlike mapping getline.
        assert!(greedy_vertex_cover_bytes(b"A\nB\n", b"A\tB\t90\t10\t1\r\n", &config).is_err());
    }

    #[test]
    fn symmetric_triplets_make_one_star_cluster() {
        let mut config = config();
        config.edge_format = "triplet".to_owned();
        config.symmetric = true;
        let output = greedy_vertex_cover("A\nB\nC\n", "A\tB\t1e-5\nA\tC\t2e-5\n", &config).unwrap();
        assert_eq!(output.centroids, ["A"]);
        assert_eq!(
            output.assignments,
            [
                ("A".into(), "A".into()),
                ("A".into(), "B".into()),
                ("A".into(), "C".into())
            ]
        );
    }

    #[test]
    fn default_edges_follow_both_coverage_directions() {
        let config = config();
        let target_covered =
            greedy_vertex_cover("A\nB\n", "A\tB\t10\t90\t1e-5\n", &config).unwrap();
        assert_eq!(target_covered.centroids, ["A"]);

        let query_covered = greedy_vertex_cover("A\nB\n", "A\tB\t90\t10\t1e-5\n", &config).unwrap();
        assert_eq!(query_covered.centroids, ["B"]);
    }

    #[test]
    fn validates_mode_and_accessions_and_keeps_first_mapping_occurrence() {
        let mut config = config();
        config.symmetric = true;
        assert_eq!(
            greedy_vertex_cover("A\n", "", &config).unwrap_err(),
            "--symmetric requires triplet edge format"
        );

        config.symmetric = false;
        let error = greedy_vertex_cover("A\n", "A\tB\t90\t90\t1\n", &config).unwrap_err();
        assert!(error.contains("Unknown target accession `B`"));

        let output = greedy_vertex_cover("A\tfirst\nA\tduplicate\nB\n", "", &config).unwrap();
        assert_eq!(output.assignments.len(), 2);
        assert_eq!(output.centroids, ["A", "B"]);
    }

    #[test]
    fn file_command_writes_source_ordered_outputs() {
        let root = std::env::temp_dir().join(format!(
            "diamond-gvc-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ));
        fs::create_dir_all(&root).unwrap();
        let mut config = config();
        config.database = root.join("mapping.tsv");
        config.edges = root.join("edges.tsv");
        config.output_file = Some(root.join("clusters.tsv"));
        config.centroid_out = Some(root.join("centroids.txt"));
        fs::write(&config.database, "A\nB\n").unwrap();
        fs::write(&config.edges, "A\tB\t10\t90\t1e-5\n").unwrap();

        let output = greedy_vertex_cover_from_files(&config).unwrap();
        assert_eq!(output.centroid_count(), 1);
        assert_eq!(
            fs::read_to_string(config.output_file.as_ref().unwrap()).unwrap(),
            "A\tA\nA\tB\n"
        );
        assert_eq!(
            fs::read_to_string(config.centroid_out.as_ref().unwrap()).unwrap(),
            "A\n"
        );

        fs::remove_dir_all(root).unwrap();
    }
}
