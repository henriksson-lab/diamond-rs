//! Source clustering search budgets and the ordinary intermediate block join.
//!
//! This module is nested in search_blocks so graph searches share the same
//! dictionary, disk streams, and HSP reconstruction as file-based searches.
#[path = "cluster_double_indexed.rs"]
mod native_owned;

use std::io::{self, Write};
use std::ops::Range;
use std::path::PathBuf;

use super::{merge_query, IntermediateFile, ReferenceDictionary};
use crate::commands::blastp::{
    self, BlastpConfig, ClusterSearchOptions, InMemorySearchEdge, SearchBlockContext,
};
use crate::commands::search_output::SearchSubjects;
use crate::data::fasta::ByteFastaRecord as FastaRecord;
use crate::output::intermediate::OutputFormat;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ClusterOwnerKind {
    Fasta,
    Dmnd,
}

/// Source Search::run receives an already opened SequenceFile. Its original
/// constructor flags and physical OIDs survive compact native record blocks.
#[derive(Clone, Debug)]
pub(crate) struct ClusterSearchOwner {
    pub kind: ClusterOwnerKind,
    pub file_name: Option<PathBuf>,
    pub sequence_count: usize,
    /// Physical owner letters, independent of filtered/scoring overrides.
    pub letters: u64,
    pub active_filter: bool,
    /// Present for original owners opened with OID_TO_ACC_MAPPING; includes
    /// rows excluded by the current cascaded DbFilter.
    pub accessions: Option<Vec<Vec<u8>>>,
    pub original_oids: Vec<u64>,
}

impl ClusterSearchOwner {
    fn fasta(records: &[FastaRecord]) -> Self {
        Self {
            kind: ClusterOwnerKind::Fasta,
            file_name: None,
            sequence_count: records.len(),
            letters: records
                .iter()
                .map(|record| record.sequence.len() as u64)
                .sum(),
            active_filter: false,
            accessions: None,
            original_oids: Vec::new(),
        }
    }

    pub(super) fn taxonomy_access(&self, oid: u64) -> io::Error {
        match self.kind {
            ClusterOwnerKind::Fasta => io::Error::other("std::exception"),
            ClusterOwnerKind::Dmnd => io::Error::other(format!(
                "Source provided DIAMOND database taxonomy lookup at OID {oid} dereferences an unloaded taxonomy table; its owner was opened without TAXON_MAPPING."
            )),
        }
    }

    fn validate_filters(&self, config: &BlastpConfig) -> io::Result<()> {
        if config.tuning.multiprocessing.enabled && self.kind == ClusterOwnerKind::Fasta {
            return Err(io::Error::other(
                "Multiprocessing mode is not compatible with FASTA databases.",
            ));
        }
        let taxon_filter = !config.taxonlist.is_empty() || !config.taxon_exclude.is_empty();
        if taxon_filter {
            if !config.taxonlist.is_empty() && !config.taxon_exclude.is_empty() {
                return Err(io::Error::other(
                    "Options --taxonlist and --taxon-exclude are mutually exclusive.",
                ));
            }
            crate::data::sequence_file::parse_taxonomy_filter(
                if config.taxonlist.is_empty() {
                    &config.taxon_exclude
                } else {
                    &config.taxonlist
                },
                ',',
            )
            .map_err(io::Error::other)?;
            // A user taxonomy filter replaces, rather than intersects, a
            // supplied cascaded DbFilter. Its lookup visits the full owner.
            if self.sequence_count != 0 {
                return Err(self.taxonomy_access(0));
            }
        }
        if !config.seqidlist.as_os_str().is_empty() {
            if taxon_filter {
                return Err(io::Error::other(
                    "--seqidlist is not compatible with taxonomy filtering.",
                ));
            }
            return Err(io::Error::other(match self.kind {
                ClusterOwnerKind::Fasta => {
                    "The FASTA database format does not support filtering by accession."
                }
                ClusterOwnerKind::Dmnd => {
                    "The .dmnd database format does not support filtering by accession."
                }
            }));
        }
        Ok(())
    }
}

fn write_auxiliary_queries(
    output: &mut crate::commands::query_output::QueryOutputs,
    records: &[FastaRecord],
    aligned: &[bool],
) -> io::Result<()> {
    for pass in [false, true] {
        for (record, &is_aligned) in records.iter().zip(aligned) {
            if is_aligned == pass {
                output.write_query_bytes(
                    is_aligned,
                    &record.id,
                    &record.sequence,
                    record.quality.as_deref().filter(|value| !value.is_empty()),
                )?;
            }
        }
    }
    Ok(())
}

pub(crate) fn run(
    config: &BlastpConfig,
    database: &[FastaRecord],
    queries: Option<&[FastaRecord]>,
    chunk_size: f64,
    database_letters: u64,
    options: ClusterSearchOptions<'_>,
) -> io::Result<Vec<InMemorySearchEdge>> {
    run_with_owner(
        config,
        database,
        queries,
        chunk_size,
        database_letters,
        options,
        &ClusterSearchOwner::fasta(database),
    )
}

pub(crate) fn run_with_owner(
    config: &BlastpConfig,
    database: &[FastaRecord],
    queries: Option<&[FastaRecord]>,
    chunk_size: f64,
    database_letters: u64,
    options: ClusterSearchOptions<'_>,
    owner: &ClusterSearchOwner,
) -> io::Result<Vec<InMemorySearchEdge>> {
    crate::basic::statistics::with_search_statistics_prepared(
        !config.tuning.multiprocessing.enabled
            || (!config.tuning.multiprocessing.initialize
                && !config.tuning.multiprocessing.recover),
        || crate::commands::source_config::begin_search_temp(config),
        || {
            run_with_owner_inner(
                config,
                database,
                queries,
                chunk_size,
                database_letters,
                options,
                owner,
            )
        },
    )
}

fn run_with_owner_inner(
    config: &BlastpConfig,
    database: &[FastaRecord],
    queries: Option<&[FastaRecord]>,
    chunk_size: f64,
    database_letters: u64,
    options: ClusterSearchOptions<'_>,
    owner: &ClusterSearchOwner,
) -> io::Result<Vec<InMemorySearchEdge>> {
    if rayon::current_thread_index().is_none() {
        return blastp::install_search_threads(config.threads, move || {
            run_with_owner_inner(
                config,
                database,
                queries,
                chunk_size,
                database_letters,
                options,
                owner,
            )
        });
    }
    let letters = chunk_size * 1e9;
    if !letters.is_finite() || letters < 1.0 {
        return Err(io::Error::other("Invalid block size."));
    }
    let budget = letters as u64;
    let input = if options.self_search {
        database
    } else {
        queries.unwrap_or(database)
    };
    let search_setup = crate::run::config::Config::from_input(&crate::run::config::ConfigInput {
        self_search: config.tuning.self_search,
        iterate: config.iterate.clone(),
        sensitivity: config.sensitivity,
        lin_stage1_query: options.linear_query || config.tuning.setup.lin_stage1_query,
        lin_stage1_target: options.linear_target || config.tuning.setup.lin_stage1_target,
        target_indexed: config.target_indexed,
        global_ranking_targets: config.global_ranking_targets,
        frame_shift: config.frame_shift,
        multiprocessing: config.tuning.multiprocessing.enabled,
        taxonlist: config.taxonlist.clone(),
        taxon_exclude: config.taxon_exclude.clone(),
        query_cover: config.query_cover,
        subject_cover: config.subject_cover,
        min_length_ratio: config.tuning.min_length_ratio,
        algo: config.tuning.algorithm.source_request(
            config.query_indexed,
            config.contiguous_seed_mode,
            config.swipe_all,
        ),
        ..Default::default()
    })
    .map_err(io::Error::other)?;
    if let Some(message) = &search_setup.iteration_message {
        search_message!("{message}");
    }
    let min_length_ratio = search_setup.min_length_ratio;
    let rounds = search_setup.sensitivity;
    owner.validate_filters(config)?;
    let mut owner_config = config.clone();
    let mut index_owner = None;
    if let Some(path) = &owner.file_name {
        owner_config.database = path.to_path_buf();
    } else if config.target_indexed {
        // A provided in-memory/subdatabase owner is a real temporary FASTA
        // SequenceFile in Source. Retain its filename through index loading.
        let directory = crate::util::io::temp_file::source_temp_directory_path(
            &config.tmpdir,
            config.output.as_deref().unwrap_or(std::path::Path::new("")),
        );
        let temporary = crate::util::io::TempFileData::init_in(false, directory, true)
            .map_err(io::Error::other)?;
        let path = std::path::PathBuf::from(&temporary.name);
        drop(temporary);
        let guard = super::TemporaryDatabase(path.clone(), config.tuning.no_unlink);
        let mut output = std::io::BufWriter::with_capacity(
            config.file_buffer_size,
            std::fs::File::create(&path)
                .map_err(|error| crate::util::io::FileOpenError::into_io(&path, error))?,
        );
        crate::data::fasta::write_fasta_bytes(
            &mut output,
            database,
            crate::basic::value::AMINO_ACID_ALPHABET,
            160,
        )?;
        output.flush()?;
        owner_config.database = path.to_path_buf();
        index_owner = Some(guard);
    }
    let config = &owner_config;
    // Keep the physical owner alive until every seed-index access completes.
    let _index_owner = index_owner;
    if config.tuning.multiprocessing.enabled {
        let mut reader_config = config.clone();
        // The provided original owner has no TAXON_MAPPING flag. Do not
        // upgrade its constructor merely because a later join requests taxids.
        reader_config.taxon_k = 0;
        reader_config.range_culling = false;
        let mut references = super::ReferenceSource::open(&reader_config)?;
        let mut queries = super::ReferenceSource::open(&reader_config)?;
        if owner.active_filter {
            let mut filter = crate::data::sequence_file::DbFilter::new(owner.sequence_count);
            for (&oid, record) in owner.original_oids.iter().zip(database) {
                if let Some(selected) = filter.oid_filter.get_mut(oid as usize) {
                    *selected = true;
                    filter.letter_count += record.sequence.len() as u64;
                }
            }
            for source in [&mut references, &mut queries] {
                if let super::ReferenceSource::Database { filter: active, .. } = source {
                    *active = Some(filter.clone());
                }
            }
        }
        super::multiprocessing::run_cluster(
            config,
            references,
            super::QueryInput::Protein(super::ProteinQueries::Database(Box::new(queries))),
            budget,
            options,
            owner,
        )?;
        return Ok(Vec::new());
    }
    native_owned::run(
        config,
        database,
        input,
        budget,
        database_letters,
        options,
        owner,
        rounds,
        min_length_ratio,
    )
}
