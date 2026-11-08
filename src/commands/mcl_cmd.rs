//! Concrete optional WITH_MCL command: native search, graph checkpoints and
//! component-at-a-time Markov expansion.
use super::blastp::BlastpConfig;
use crate::mcl::clustering_format::ClusteringFilter;
use crate::mcl::mcl::{GraphHandle, Mcl, MclConfig, MclStatistics, MASK_INVERSE};
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq)]
pub struct MclOptions {
    pub similarity: Vec<u8>,
    pub threshold: Option<f64>,
    pub graph_file: PathBuf,
    pub restart: bool,
    pub statistics: bool,
    pub component_chunk_size: u32,
    /// Preserve the Source distinction between an absent and explicit -k.
    pub max_target_seqs: Option<i64>,
    pub config: MclConfig,
}

impl Default for MclOptions {
    fn default() -> Self {
        Self {
            similarity: Vec::new(),
            threshold: None,
            graph_file: PathBuf::new(),
            restart: false,
            statistics: false,
            component_chunk_size: 1,
            max_target_seqs: None,
            config: MclConfig::default(),
        }
    }
}

/// Source MCL uses exception-disabled cout/ofstream, not OutputFile. A failed
/// open/write sets sticky stream state; identifier reads still continue.
struct MclOutput {
    stream: Option<Box<dyn Write>>,
    failed: bool,
}
impl MclOutput {
    fn new(path: Option<&Path>) -> Self {
        let stream: Option<Box<dyn Write>> = match path.filter(|path| !path.as_os_str().is_empty())
        {
            Some(path) => File::create(path)
                .ok()
                .map(|file| Box::new(file) as Box<dyn Write>),
            None => Some(Box::new(crate::util::io::console::stdout())),
        };
        Self {
            stream,
            failed: false,
        }
    }
    fn row(&mut self, id: &[u8], cluster: Option<u64>, label: char) {
        if self.failed {
            return;
        }
        if let Some(output) = self.stream.as_mut() {
            let result = output.write_all(id).and_then(|()| match cluster {
                Some(cluster) => writeln!(output, "\t{cluster}\t{label}"),
                None => writeln!(output, "\t-1\tu"),
            });
            self.failed = result.and_then(|()| output.flush()).is_err();
        }
    }
}

pub fn run(search: &BlastpConfig, options: &MclOptions) -> io::Result<()> {
    if search.database.as_os_str().is_empty() {
        return Err(io::Error::other(
            "Missing parameter: database file (--db/-d)",
        ));
    }
    if options.component_chunk_size == 0 {
        return Err(io::Error::other(
            "MCL component chunk size must be positive",
        ));
    }
    let mut search = search.clone();
    let expression = if options.similarity.is_empty() {
        b"normalized_bitscore_global".as_slice()
    } else {
        options.similarity.as_slice()
    };
    let minimum = if options.similarity.is_empty() {
        options
            .threshold
            .or(Some(50.0))
            .filter(|value| *value != 0.0)
    } else {
        options.threshold
    };
    if !options.restart && !options.similarity.is_empty() && minimum.is_none() {
        crate::util::io::console::source_eprintln!("WARNING: It is recommended to set a threshold value for the clustering similarity measure (option `--cluster-threshold`).");
    }
    if !options.restart {
        search.tuning.output.clustering_filter = Some(Arc::new(
            ClusteringFilter::new_bytes(expression, minimum).map_err(io::Error::other)?,
        ));
        search.outfmt = vec![b"clus".to_vec(), expression.to_vec()];
        search.no_self_hits = false;
        search.max_target_seqs = options.max_target_seqs.unwrap_or(i64::MAX);
        search = super::blastp::normalize_output_config(&search)?.into_owned();
    } else {
        // Source's main validates expression syntax independently. Restart
        // itself skips search-format setup and all search-only normalization.
        super::blastp::preflight_file_buffer_size(search.file_buffer_size)?;
    }
    let checkpoint =
        (!options.graph_file.as_os_str().is_empty()).then_some(options.graph_file.as_path());
    let (mut graph, mut identifiers) = super::search_blocks::run_mcl_search(
        &search,
        checkpoint,
        options.restart,
        options.config.symmetric,
    )?;
    let components = graph.components();
    let threads = if search.threads > 0 {
        search.threads as usize
    } else {
        crate::util::system::hardware_concurrency() as usize
    };
    let fixed_chunk = options.component_chunk_size as usize;
    // Source assigns each worker a strided cluster-ID lane and adapts later
    // component claims to the total number of nodes in its previous batch.
    // Source's zero-worker/zero-progress boundary is undefined; retain one
    // worker and a positive claim for nonempty graphs.
    let workers = if components.is_empty() {
        0
    } else {
        threads.min(components.len() / fixed_chunk).max(1)
    };
    let next_component = AtomicUsize::new(workers * fixed_chunk);
    let threads_done = AtomicUsize::new(0);
    let excess_threads = threads.saturating_sub(workers);
    let max_job_size = components
        .iter()
        .take(fixed_chunk)
        .map(Vec::len)
        .sum::<usize>();
    let worker_results = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for worker in 0..workers {
            let graph = &graph;
            let components = &components;
            let next_component = &next_component;
            let threads_done = &threads_done;
            handles.push(scope.spawn(move || {
                let started = std::time::Instant::now();
                let mut results = Vec::new();
                let mut cluster_id = worker as u64;
                let mut begin = worker * fixed_chunk;
                let mut claim = fixed_chunk;
                while begin < components.len() {
                    let end = begin.saturating_add(claim).min(components.len());
                    let mut nodes = 0usize;
                    for index in begin..end {
                        let component = &components[index];
                        nodes += component.len();
                        let edges = if component.len() > 1 {
                            graph.component(component)?
                        } else {
                            Vec::new()
                        };
                        let handle = GraphHandle {
                            sequence_count: component.len(),
                            symmetric: options.config.symmetric,
                            edges,
                            similarity: expression.to_vec(),
                            threshold: minimum,
                        };
                        let ids = vec![String::new(); component.len()];
                        let get_threads = || {
                            let done = threads_done.load(Ordering::Relaxed);
                            let remaining = workers.saturating_sub(done).max(1);
                            1 + (done + excess_threads) / remaining
                                + usize::from(worker < (done + excess_threads) % remaining)
                        };
                        let result = Mcl::default()
                            .run_component_with_threads(&ids, &handle, options.config, get_threads)
                            .map_err(io::Error::other)?;
                        let assignments = result
                            .encoded_assignments
                            .into_iter()
                            .map(|value| {
                                if value & MASK_INVERSE == 0 {
                                    value
                                } else {
                                    (value & MASK_INVERSE)
                                        | (cluster_id + (value & !MASK_INVERSE) * workers as u64)
                                }
                            })
                            .collect::<Vec<_>>();
                        cluster_id += result.statistics.clusters as u64 * workers as u64;
                        results.push((index, assignments, result.statistics));
                    }
                    claim = (end - begin)
                        .saturating_mul(max_job_size / nodes.max(1))
                        .max(1);
                    begin = next_component.fetch_add(claim, Ordering::Relaxed);
                }
                threads_done.fetch_add(1, Ordering::Relaxed);
                Ok::<_, io::Error>((results, started.elapsed().as_secs_f32()))
            }));
        }
        handles
            .into_iter()
            .map(|handle| {
                handle
                    .join()
                    .map_err(|_| io::Error::other("MCL clustering worker panicked"))?
            })
            .collect::<io::Result<Vec<_>>>()
    })?;
    search_message!(
        "Jobs per thread: {}",
        worker_results
            .iter()
            .map(|(results, _)| format!(" {:>8}", results.len()))
            .collect::<String>()
    );
    search_message!(
        "Time per thread: {}",
        worker_results
            .iter()
            .map(|(_, seconds)| format!(" {:>8}", seconds))
            .collect::<String>()
    );
    let mut results = worker_results
        .into_iter()
        .flat_map(|(results, _)| results)
        .collect::<Vec<_>>();
    results.sort_by_key(|(index, _, _)| *index);
    let mut encoded = vec![0u64; identifiers.len()];
    let mut stats = MclStatistics {
        components: components.len(),
        elements: graph.edge_count() as usize,
        ..Default::default()
    };
    for (index, assignments, local) in results {
        for (&oid, value) in components[index].iter().zip(assignments) {
            encoded[oid] = value;
        }
        stats.clusters += local.clusters;
        stats.singleton_clusters += local.singleton_clusters;
        stats.non_singleton_components += local.non_singleton_components;
        stats.dense_calculations += local.dense_calculations;
        stats.sparse_calculations += local.sparse_calculations;
        stats.failed_to_converge += local.failed_to_converge;
        stats.component_sizes.extend(local.component_sizes);
        stats.sparsities.extend(local.sparsities);
        stats.average_neighbors.extend(local.average_neighbors);
    }
    if options.statistics {
        print_statistics(&stats, &components, options.config, workers)?;
    }
    search_message!(
        "Clusters found {} ({} incl. singletons)",
        stats.clusters - stats.singleton_clusters,
        stats.clusters
    );
    search_message!(
        "\t number of failed calculations {}",
        stats.failed_to_converge
    );
    search_message!(
        "\t number of dense calculations {}",
        stats.dense_calculations
    );
    search_message!(
        "\t number of sparse calculations {}",
        stats.sparse_calculations
    );
    let mut output = MclOutput::new(search.output.as_deref());
    identifiers.begin_output()?;
    for value in encoded {
        let title = identifiers.next_title()?;
        let label = match value & MASK_INVERSE {
            crate::mcl::mcl::MASK_SINGLE_NODE => 's',
            crate::mcl::mcl::MASK_ATTRACTOR_NODE => 'a',
            crate::mcl::mcl::MASK_NORMAL_NODE => 'n',
            _ => 'u',
        };
        // Source seqid scans a C string as bytes; an EXTRA physical cursor
        // can expose non-UTF8 title bytes, which must not become U+FFFD.
        let end = title
            .iter()
            .position(|byte| {
                *byte == 0
                    || crate::util::sequence::ID_DELIMITERS
                        .as_bytes()
                        .contains(byte)
            })
            .unwrap_or(title.len());
        let id = &title[..end];
        output.row(
            id,
            (label != 'u').then_some((value & !MASK_INVERSE) + 1),
            label,
        );
    }
    Ok(())
}

fn print_statistics(
    stats: &MclStatistics,
    components: &[Vec<usize>],
    config: MclConfig,
    threads: usize,
) -> io::Result<()> {
    search_message!("Number of DIAMOND hits:          {}", stats.elements);
    search_message!(
        "Number of independet components: {} ({} incl. singletons)",
        stats.non_singleton_components,
        stats.components
    );
    let median = |values: &[f32]| {
        if values.is_empty() {
            0.0
        } else if values.len() % 2 == 0 {
            (values[values.len() / 2 - 1] + values[values.len() / 2]) / 2.0
        } else {
            values[values.len() / 2]
        }
    };
    let sizes = components
        .iter()
        .map(|component| component.len() as f32)
        .collect::<Vec<_>>();
    let non_single = sizes
        .iter()
        .copied()
        .filter(|&size| size > 1.0)
        .collect::<Vec<_>>();
    search_message!("Component size information:");
    search_message!("\tmax. : {:>12}", sizes.first().copied().unwrap_or(0.0));
    search_message!(
        "\tmed. : {:>12} ({} incl. singletons)",
        median(&non_single),
        median(&sizes)
    );
    search_message!(
        "\tmin. : {:>12} ({} incl. singletons)",
        non_single.last().copied().unwrap_or(0.0),
        sizes.last().copied().unwrap_or(0.0)
    );
    let mut memories = Vec::new();
    for ((&size, &sparsity), &neighbors) in stats
        .component_sizes
        .iter()
        .zip(&stats.sparsities)
        .zip(&stats.average_neighbors)
    {
        memories.push(if sparsity >= config.sparsity_switch {
            size as f32 * (1.0 + neighbors.powf(config.expansion)) * 12.0
        } else {
            4.0 * (size * size) as f32
        });
    }
    memories.sort_by(|a, b| b.total_cmp(a));
    for (name, values) in [
        ("Sparsity", &stats.sparsities),
        ("Average number of neighbors", &stats.average_neighbors),
    ] {
        let mut sorted = values.clone();
        sorted.sort_by(|a, b| b.total_cmp(a));
        search_message!("{name} of components (excluding singletons):");
        search_message!(
            "\tmax. : {:>12} - at max. size: {:>12}",
            sorted.first().copied().unwrap_or(0.0),
            values.first().copied().unwrap_or(0.0)
        );
        search_message!(
            "\tmed. : {:>12} - at med. size: {:>12}",
            median(&sorted),
            median(values)
        );
        search_message!(
            "\tmin. : {:>12} - at min. size: {:>12}",
            sorted.last().copied().unwrap_or(0.0),
            values.last().copied().unwrap_or(0.0)
        );
    }
    let memory = stats.elements as f32 * 12.0
        + memories
            .into_iter()
            .take(threads)
            .map(|value| 3.0 * value)
            .sum::<f32>();
    search_message!("Rough memory requirements: {memory}");
    Ok(())
}
