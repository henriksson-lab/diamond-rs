//! MCL's self-database search uses the ordinary streamed search pipeline.
use super::{ProteinQueries, QueryInput, ReferenceSource};
use crate::commands::blastp::{self, BlastpConfig, NativeSearchWriter};
use crate::mcl::sparse_matrix_stream::SparseMatrixStream;
use std::io::{self, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

/// The actual Source database owner survives search/checkpoint processing.
/// Identifiers are decoded one record at a time only after output opens.
pub(crate) struct MclIdentifiers {
    references: ReferenceSource,
    count: usize,
    cursor: usize,
}
impl MclIdentifiers {
    pub(crate) fn len(&self) -> usize {
        self.count
    }
    pub(crate) fn begin_output(&mut self) -> io::Result<()> {
        // Source MCL explicitly calls init_seq_access after output opens.
        // DMND must seek the fixed sizeof(Header1)+8+sizeof(Header2) cursor,
        // not an OID's indexed position or the constructor's current cursor.
        match &mut self.references {
            ReferenceSource::Database { file, .. } => file.init_source_sequential()?,
            ReferenceSource::Blast { file, .. } => file.init_seq_access(),
            ReferenceSource::Fasta(file) => file.cursor = 0,
            ReferenceSource::Memory { cursor, .. } => *cursor = 0,
        }
        self.cursor = 0;
        Ok(())
    }
    pub(crate) fn next_title(&mut self) -> io::Result<Vec<u8>> {
        if self.cursor >= self.count {
            return Err(io::Error::other("MCL identifier cursor outside database"));
        }
        // Follow the real sequential owner, including Source read_seq's
        // delimiter/EOF semantics; an indexed record lookup changes them.
        let title = match &mut self.references {
            ReferenceSource::Database { file, .. } => file.read_source_sequential_bytes()?.1,
            ReferenceSource::Blast { file, .. } => {
                file.read_seq_bytes().map_err(io::Error::other)?.1
            }
            ReferenceSource::Fasta(file) => {
                let record = file.record(file.cursor as u64)?;
                file.cursor += 1;
                record.id
            }
            ReferenceSource::Memory { records, cursor } => {
                let title = records
                    .get(*cursor)
                    .ok_or_else(|| io::Error::other("MCL identifier cursor outside database"))?
                    .id
                    .clone();
                *cursor += 1;
                title
            }
        };
        self.cursor += 1;
        Ok(title)
    }
}

struct GraphConsumer(Arc<Mutex<SparseMatrixStream>>);
impl Write for GraphConsumer {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0
            .lock()
            .map_err(|_| io::Error::other("MCL graph lock poisoned"))?
            .consume(bytes)?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

pub(crate) fn run(
    config: &BlastpConfig,
    graph_file: Option<&Path>,
    restart: bool,
    symmetric: bool,
) -> io::Result<(SparseMatrixStream, MclIdentifiers)> {
    if restart {
        run_inner(config, graph_file, restart, symmetric)
    } else {
        crate::basic::statistics::with_search_statistics_prepared(
            !config.tuning.multiprocessing.enabled
                || (!config.tuning.multiprocessing.initialize
                    && !config.tuning.multiprocessing.recover),
            || crate::commands::source_config::begin_search_temp(config),
            || run_inner(config, graph_file, restart, symmetric),
        )
    }
}

fn run_inner(
    config: &BlastpConfig,
    graph_file: Option<&Path>,
    restart: bool,
    symmetric: bool,
) -> io::Result<(SparseMatrixStream, MclIdentifiers)> {
    let mut reader_config = config.clone();
    if restart {
        reader_config.outfmt = vec!["6".into()];
        reader_config.tuning.output.clustering_filter = None;
    }
    let references = ReferenceSource::open(&reader_config)?;
    let count = references.totals().0;
    let count = usize::try_from(count).map_err(io::Error::other)?;
    if restart {
        let path = graph_file.ok_or_else(|| io::Error::other("Missing MCL graph file"))?;
        let graph = SparseMatrixStream::open(path, symmetric)?;
        if graph.sequence_count() != count {
            return Err(io::Error::other("Graph and sequence counts differ"));
        }
        return Ok((
            graph,
            MclIdentifiers {
                references,
                count,
                cursor: 0,
            },
        ));
    }
    // Source MPP sends joined Clus records to per-query chunk files instead
    // of the supplied master graph consumer. Its graph consequently stays
    // empty and the MCL phase emits singleton assignments. Reuse that exact
    // lifecycle rather than inventing a graph merge for this combination.
    let tmpdir = crate::util::io::temp_file::source_temp_directory(
        &config.tmpdir,
        config.output.as_deref().unwrap_or(""),
    );
    let graph = Arc::new(Mutex::new(SparseMatrixStream::create(
        count,
        symmetric,
        graph_file,
        &tmpdir,
        config.block_size.unwrap_or(0.0),
    )?));
    let writer = NativeSearchWriter::new_with_consumer(
        config,
        count as u64,
        Box::new(GraphConsumer(graph.clone())),
    )?;
    let queries = ReferenceSource::open(config)?;
    let budget = blastp::native_block_budget(config)?
        .ok_or_else(|| io::Error::other("Missing MCL search block budget"))?;
    let references = super::run_protein_with_writer_retaining_references(
        config,
        references,
        QueryInput::Protein(ProteinQueries::Database(Box::new(queries))),
        budget,
        Some(writer),
    )?;
    let mut graph = Arc::try_unwrap(graph)
        .map_err(|_| io::Error::other("MCL graph consumer is still active"))?
        .into_inner()
        .map_err(|_| io::Error::other("MCL graph lock poisoned"))?;
    graph.finish()?;
    Ok((
        graph,
        MclIdentifiers {
            references,
            count,
            cursor: 0,
        },
    ))
}
