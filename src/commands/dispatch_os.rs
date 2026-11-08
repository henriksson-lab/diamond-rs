//! Canonical Source argv parsing: option cells retain OS paths or narrow bytes.
use super::*;
use crate::util::command_line_parser::ReadOption;
use std::ffi::{OsStr, OsString};
use std::io;
use std::path::{Path, PathBuf};

pub(super) struct Options {
    pub(super) values: crate::commands::argument_registry::OsValues,
    pub(super) command: String,
}
fn raw_error(prefix: &[u8], value: &[u8], suffix: &[u8]) -> io::Error {
    let mut message = prefix.to_vec();
    message.extend_from_slice(value);
    message.extend_from_slice(suffix);
    crate::util::io::RawDiagnosticError(message).into_io()
}
impl Options {
    pub(super) fn new(args: &[OsString]) -> io::Result<Self> {
        let values = crate::commands::argument_registry::values_os(args)?;
        let command = super::super::argument_registry::command_name_os(args)?.unwrap_or("");
        Ok(Self {
            values,
            command: command.to_owned(),
        })
    }
    pub(super) fn has(&self, id: &str) -> bool {
        self.values.contains_key(id)
    }
    pub(super) fn number<T: ReadOption>(&self, id: &str, mut default: T) -> io::Result<T> {
        if let Some(values) = self.values.get(id) {
            T::read_option_os(&mut default, values)?;
        }
        Ok(default)
    }
    pub(super) fn path(&self, id: &str) -> Option<PathBuf> {
        self.values
            .get(id)
            .and_then(|v| v.first())
            .map(PathBuf::from)
    }
    pub(super) fn paths(&self, id: &str) -> Vec<PathBuf> {
        self.values
            .get(id)
            .into_iter()
            .flatten()
            .map(PathBuf::from)
            .collect()
    }
    pub(super) fn raw(&self, id: &str) -> io::Result<Option<Vec<u8>>> {
        self.values
            .get(id)
            .and_then(|v| v.first())
            .map(|v| {
                crate::util::io::source_os_bytes(v)
                    .map(|v| v.into_owned())
                    .map_err(io::Error::other)
            })
            .transpose()
    }
    pub(super) fn raw_vector(&self, id: &str) -> io::Result<Option<Vec<Vec<u8>>>> {
        self.values
            .get(id)
            .map(|values| {
                values
                    .iter()
                    .map(|v| {
                        crate::util::io::source_os_bytes(v)
                            .map(|v| v.into_owned())
                            .map_err(io::Error::other)
                    })
                    .collect()
            })
            .transpose()
    }
    // Explicit checked compatibility for callers whose semantic owner still
    // accepts UTF-8. Native paths and unconstrained text never use this hook.
    pub(super) fn text(&self, id: &str) -> io::Result<Option<String>> {
        self.raw(id)?
            .map(|v| String::from_utf8(v).map_err(io::Error::other))
            .transpose()
    }
    fn enum_text(&self, id: &str, allowed: &[&str], permitted: &str) -> io::Result<Option<String>> {
        let Some(value) = self.raw(id)? else {
            return Ok(None);
        };
        if let Some(text) = allowed.iter().find(|s| s.as_bytes() == value) {
            return Ok(Some((*text).to_owned()));
        }
        Err(raw_error(
            b"Invalid value for string field: ",
            &value,
            format!(". Permitted values: {permitted}").as_bytes(),
        ))
    }
    fn algorithm(&self) -> io::Result<Option<String>> {
        self.enum_text("algo", &["", "0", "1", "ctg"], ", 0, 1, ctg")
    }
    fn extension(&self) -> io::Result<Option<crate::search::sensitivity::ExtensionMode>> {
        let value = self.enum_text(
            "ext",
            &["", "banded-fast", "banded-slow", "full", "global", "none"],
            "banded-fast, banded-slow, full, global, none",
        )?;
        value
            .filter(|v| !v.is_empty())
            .map(|s| {
                crate::search::sensitivity::ExtensionMode::from_string(&s).map_err(io::Error::other)
            })
            .transpose()
    }
    fn sensitivity(&self) -> io::Result<Sensitivity> {
        let mut selected = None;
        for (id, level) in [
            ("faster", Sensitivity::Faster),
            ("fast", Sensitivity::Fast),
            ("mid-sensitive", Sensitivity::MidSensitive),
            ("sensitive", Sensitivity::Sensitive),
            ("more-sensitive", Sensitivity::MoreSensitive),
            ("very-sensitive", Sensitivity::VerySensitive),
            ("ultra-sensitive", Sensitivity::UltraSensitive),
            ("shapes-30x10", Sensitivity::Shapes30x10),
            ("shapes-6x10", Sensitivity::Shapes6x10),
            ("linclust-20", Sensitivity::Linclust20),
        ] {
            if self.has(id) {
                if selected.is_some() {
                    return Err(io::Error::other(
                        "Sensitivity switches are mutually exclusive.",
                    ));
                }
                selected = Some(level);
            }
        }
        Ok(selected.unwrap_or(Sensitivity::Default))
    }
}

pub(crate) fn parse_search_string_option_os(
    args: &[OsString],
    name: &str,
) -> io::Result<Option<String>> {
    Options::new(args)?.text(name)
}
pub(crate) fn parse_signed_search_option_os(
    args: &[OsString],
    name: &str,
    default: i32,
) -> io::Result<i32> {
    Options::new(args)?.number(name, default)
}
pub(crate) fn parse_float_search_option_os(
    args: &[OsString],
    name: &str,
    default: f64,
) -> io::Result<f64> {
    Options::new(args)?.number(name, default)
}
pub(crate) fn parse_native_search_tuning_os(
    args: &[OsString],
) -> io::Result<crate::commands::blastp::SearchTuning> {
    parse_tuning(&Options::new(args)?, true)
}
pub(crate) fn parse_source_config_tuning_os(
    args: &[OsString],
) -> io::Result<crate::commands::blastp::SearchTuning> {
    parse_tuning(&Options::new(args)?, false)
}

fn cpp_atoll(value: &[u8]) -> i64 {
    let mut v = value.to_vec();
    v.push(0);
    unsafe { libc::atoll(v.as_ptr().cast()) as i64 }
}

fn parse_tuning(
    options: &Options,
    resolve: bool,
) -> io::Result<crate::commands::blastp::SearchTuning> {
    let mut tuning = crate::commands::blastp::SearchTuning::default();
    tuning.algorithm.requested = match if resolve { options.algorithm()? } else { None }.as_deref()
    {
        Some("0") => crate::run::config::Algo::DoubleIndexed,
        Some("1") => crate::run::config::Algo::QueryIndexed,
        Some("ctg") => crate::run::config::Algo::CtgSeed,
        _ => crate::run::config::Algo::Auto,
    };
    tuning.output.qnum_offset = options
        .raw("qnum-offset")?
        .map_or(0, |value| cpp_atoll(&value));
    tuning.output.snum_offset = options
        .raw("snum-offset")?
        .map_or(0, |value| cpp_atoll(&value));
    tuning.output.xml_blord_format = options.has("xml-blord-format");
    tuning.output.no_parse_seqids = options.has("no-parse-seqids");
    tuning.matrix_file = options.path("custom-matrix").unwrap_or_default();
    tuning.stop_match_score = options.number("stop-match-score", 1)?;
    tuning.matrix_scale = options.number("cbs-matrix-scale", 1)?;
    tuning.cutoff_score_8bit = options.number("cutoff-score-8bit", 240)?;
    tuning.max_swipe_dp = options
        .raw("max-swipe-dp")?
        .map_or(1_000_000, |value| cpp_atoll(&value));
    tuning.approx_backtrace = options.has("approx-backtrace");
    tuning.cbs_window = options.number("cbs-window", 40)?;
    tuning.max_motif_len = options.number("max-motif-len", 30)?;
    tuning.tantan_min_mask_prob = options.number("tantan-minMaskProb", 0.9)?;
    tuning.db_size = options
        .raw("dbsize")?
        .map_or(0, |value| cpp_atoll(&value) as u64);
    tuning.no_ranking = options.has("no-ranking");
    tuning.no_reextend = options.has("no-reextend");
    tuning.mapany = options.has("mapany");
    tuning.self_search = options.has("self");
    tuning.add_self_aln = options.has("add-self-aln");
    tuning.target_hard_cap = options
        .raw("target_hard_cap")?
        .map_or(0, |value| cpp_atoll(&value));
    tuning.ranking_score_drop_factor = options.number("score-drop-factor", 0.95)?;
    tuning.ranking_cutoff_bitscore = options.number("ranking-cutoff-bitscore", 25.0)?;
    tuning.load_balancing = options.number("load-balancing", 0)?;
    tuning.log_query = options.has("log-query");
    tuning.log_extend = options.has("log-extend");
    tuning.target_parallel_verbosity = options.number("target-parallel-verbosity", u32::MAX)?;
    tuning.threads_align = options.number("palign", 0)?;
    tuning.band_bin = options.number("band-bin", 24)?;
    tuning.col_bin = options.number("col-bin", 400)?;
    tuning.min_task_trace_pts = options
        .raw("min_task_trace_pts")?
        .map_or(1024, |value| cpp_atoll(&value));
    tuning.trace_pt_fetch_size = options
        .raw("trace-pt-fetch-size")?
        .map_or(10_000_000_000, |value| cpp_atoll(&value) as u64);
    tuning.hamming_ext = options.has("hamming-ext");
    tuning.diag_filter_id = if options.has("diag-filter-id") {
        Some(options.number("diag-filter-id", 0.0)?)
    } else {
        None
    };
    tuning.diag_filter_cov = if options.has("diag-filter-cov") {
        Some(options.number("diag-filter-cov", 0.0)?)
    } else {
        None
    };
    tuning.mutual_cover = if options.has("mutual-cover") {
        Some(options.number("mutual-cover", 0.0)?)
    } else {
        None
    };
    tuning.chaining.max_gap = options.number("chaining-maxgap", 2000)?;
    tuning.chaining.range_cover = options
        .raw("chaining-range-cover")?
        .map_or(8, |value| cpp_atoll(&value) as usize);
    tuning.chaining.max_nodes = options
        .raw("chaining-maxnodes")?
        .map_or(0, |value| cpp_atoll(&value) as usize);
    tuning.chaining.length_cap = options.number("chaining-len-cap", 2.0)?;
    tuning.chaining.min_nodes = options
        .raw("chaining-min-nodes")?
        .map_or(200, |value| cpp_atoll(&value) as usize);
    tuning.chaining.no_merge_hsps = options.has("no_chaining_merge_hsps");
    tuning.chaining.stacked_hsp_ratio = options.number("chaining-stacked-hsp-ratio", 0.5)?;

    tuning.hash_join = crate::util::algo::HashJoinConfig {
        split_size: options.number("join-split-size", 100_000)? as usize,
        split_key_len: options.number("join-split-key-len", 17)?,
        radix_bits: options.number("radix-bits", 8)?,
        ht_factor: options.number("join-ht-factor", 1.3)?,
        swap: options.has("hash_join_swap"),
        radix_cluster_buffered: options.has("radix-cluster-buffered"),
    };
    tuning.kmer_ranking = options.has("kmer-ranking");
    tuning.min_length_ratio = options.number("min-len-ratio", 0.0)?;
    tuning.rank_ratio = options.number("rank-ratio", -1.0)?;
    tuning.rank_factor = options.number("rank-factor", -1.0)?;
    tuning.heartbeat = options.has("heartbeat");
    tuning.load_threads = options.number("load-threads", 16)?;
    tuning.minichunk = options
        .raw("minichunk")?
        .map_or(262_144, |value| cpp_atoll(&value) as usize);
    tuning.unaligned_targets = options.path("unaligned-targets").unwrap_or_default();
    tuning.setup.lin_stage1_target = options.has("linsearch");
    tuning.setup.lin_stage1_query = options.has("lin-stage1");
    tuning.setup.lin_stage1_combo = options.has("lin-combo");
    tuning.setup.linclust_banded_ext = options.has("linclust-banded-ext");
    tuning.tile_size = options.number("tile-size", 1024)? as usize;
    tuning.swipe_task_size = options
        .raw("swipe-task-size")?
        .map_or(100_000_000, |value| cpp_atoll(&value));
    tuning.setup.anchored_swipe = options.has("anchored-swipe");
    tuning.padding = options.raw("band")?.map(|value| cpp_atoll(&value) as i32);
    tuning.ignore_warnings = options.has("ignore-warnings");
    tuning.no_unlink = options.has("no-unlink");
    tuning.no_block_size_limit = options.has("no-block-size-limit");
    tuning.hit_membuf = options.has("hit-membuf");
    tuning.multiprocessing.enabled = options.has("multiprocessing");
    tuning.multiprocessing.initialize = options.has("mp-init");
    tuning.multiprocessing.recover = options.has("mp-recover");
    tuning.multiprocessing.self_compare = options.has("mp-self");
    tuning.multiprocessing.query_chunk = i64::from(options.number("mp-query-chunk", -1)?);
    tuning.setup.aln_out = options.path("aln-out").unwrap_or_default();
    tuning.output.cluster_threshold = if options.has("cluster-threshold") {
        Some(options.number("cluster-threshold", 0.0)?)
    } else {
        None
    };
    tuning.setup.parallel_tmpdir = options.path("parallel-tmpdir").unwrap_or_default();
    tuning.multiprocessing.directory = tuning.setup.parallel_tmpdir.clone().into();
    tuning.setup.shape_mask = options.raw_vector("shape-mask")?.unwrap_or_default();
    tuning.setup.shapes = options.number("shapes", 0u32)?;
    tuning.setup.seed_cut = options.number("seed-cut", 0.0)?;
    tuning.setup.freq_masking = options.has("freq-masking");
    tuning.setup.freq_sd = Some(options.number("freq-sd", 0.0)?);
    tuning.setup.min_identities = Some(options.number("id2", 0)?);
    tuning.setup.query_bins = options.raw("bin")?.map(|value| cpp_atoll(&value) as u32);
    tuning.setup.minimizer_window = Some(options.number("minimizer-window", 0)?);
    tuning.setup.sketch_size = Some(options.number("sketch-size", 0)?);
    tuning.setup.ungapped_evalue = Some(options.number("ungapped-evalue", -1.0)?);
    tuning.setup.ungapped_evalue_short = Some(options.number("ungapped-evalue-short", -1.0)?);
    tuning.setup.gapped_filter_evalue = Some(options.number("gapped-filter-evalue", -1.0)?);
    tuning.setup.gapped_filter_diag_bit_score = options.number("gapped-filter-diag-score", 12.0)?;
    tuning.short_query_ungapped_bitscore = options.number("short-query-ungapped-bitscore", 25.0)?;
    tuning.short_query_max_len = options.number("short-query-max-len", 60)?;
    tuning.store_query_quality = options.has("store-query-quality");
    tuning.ungapped_window = options.number("ungapped-window", 48)?;
    tuning.left_most_interval = options.number("left-most-interval", 32)?;
    tuning.gapped_filter_evalue1 = options.number("gapped-filter-evalue1", 2000.0)?;
    tuning.gapped_filter_window = options.number("gapped-filter-window", 200)?;
    tuning.swipe_chunk_size = options.number("swipe-chunk-size", 256)? as u32;
    tuning.transcript_len_estimate = options.number("transcript-len-estimate", 1.0)?;
    tuning.narrow_band_cov = options.number("narrow-band-cov", 0.0)?;
    tuning.narrow_band_factor = options.number("narrow-band-factor", 0.0)?;
    tuning.min_band_overlap = options.number("min-band-overlap", 0.0)?;
    tuning.cbs_err_tolerance = options.number("cbs-err-tolerance", 1e-8)?;
    tuning.cbs_it_limit = options.number("cbs-it-limit", 2000)?;
    tuning.anchor_window = options.number("anchor-window", 12)?;
    tuning.anchor_score = options.number("anchor-score", 1.0)?;
    tuning.matrix_adjust_thresholds = crate::stats::cbs::CbsThresholds::new(
        options.number("query-match-distance-threshold", -1.0)?,
        options.number("length-ratio-threshold", -1.0)?,
        options.number("cbs-angle", -1.0)?,
    );
    if resolve {
        tuning
            .algorithm
            .begin_search(options.sensitivity()?, tuning.algorithm.requested);
    }
    Ok(tuning)
}

fn append_path(path: &mut PathBuf, suffix: &str) {
    if !path.as_os_str().is_empty()
        && !path
            .as_os_str()
            .as_encoded_bytes()
            .ends_with(suffix.as_bytes())
    {
        path.as_mut_os_string().push(suffix);
    }
}
fn output_options(options: &Options) -> io::Result<(Option<PathBuf>, Vec<Vec<u8>>)> {
    let mut output = options.path("out");
    let mut format = options.raw_vector("outfmt")?.unwrap_or_default();
    if let Some(mut daa) = options.path("daa") {
        if output.is_some() {
            return Err(io::Error::other(
                "Options --daa and --out cannot be used together.",
            ));
        }
        if format
            .first()
            .is_some_and(|s| !s.eq_ignore_ascii_case(b"daa"))
        {
            return Err(io::Error::other("Invalid parameter: --daa/-a. Output file is specified with the --out/-o parameter."));
        }
        if !options.has("no-auto-append") {
            append_path(&mut daa, ".daa");
        }
        output = Some(daa);
        if format.is_empty() {
            format.push(b"daa".to_vec());
        }
    } else if format
        .first()
        .is_some_and(|s| s == b"100" || s.eq_ignore_ascii_case(b"daa"))
    {
        let path = output
            .as_mut()
            .ok_or_else(|| io::Error::other("DAA output requires -a/--daa or -o/--out."))?;
        if !options.has("no-auto-append") {
            append_path(path, ".daa");
        }
    }
    let compression = options.raw("compress")?.unwrap_or_default();
    if !compression.is_empty() && format.first().is_some_and(|s| s == b"100" || s == b"daa") {
        return Err(io::Error::other(
            "Compression is not supported for DAA format.",
        ));
    }
    if !options.has("no-auto-append") {
        if let Some(path) = output.as_mut() {
            if compression == b"1" {
                append_path(path, ".gz");
            } else if compression == b"zstd" {
                append_path(path, ".zst");
            }
        }
    }
    Ok((output, format))
}
fn masking(options: &Options) -> io::Result<crate::masking::MaskingMode> {
    let value = options
        .raw("masking")?
        .unwrap_or_else(|| b"tantan".to_vec());
    let mut lower = value.clone();
    lower.make_ascii_lowercase();
    match lower.as_slice() {
        b"0" | b"none" => Ok(crate::masking::MaskingMode::None),
        b"1" | b"tantan" => Ok(crate::masking::MaskingMode::Tantan),
        b"seg" => Ok(crate::masking::MaskingMode::BlastSeg),
        _ => Err(raw_error(
            b"Invalid value for string field: ",
            &value,
            b". Permitted values: none, tantan, seg",
        )),
    }
}
fn soft_masking(options: &Options) -> io::Result<crate::masking::MaskingAlgo> {
    let value = options.raw("soft-masking")?.unwrap_or_default();
    match value.as_slice() {
        b"" | b"0" | b"none" => Ok(crate::masking::MaskingAlgo::None),
        b"tantan" => Ok(crate::masking::MaskingAlgo::Tantan),
        b"seg" => Ok(crate::masking::MaskingAlgo::Seg),
        _ => Err(raw_error(
            b"Invalid value for string field: ",
            &value,
            b". Permitted values: 0, none, seg, tantan",
        )),
    }
}
fn memory_limit(options: &Options) -> io::Result<u64> {
    match options.raw("memory-limit")? {
        None => Ok(DEFAULT_NATIVE_MEMORY_LIMIT),
        Some(value) => {
            let size = crate::util::string::interpret_number_bytes(&value)?;
            u64::try_from(size).map_err(|_| {
                raw_error(
                    b"invalid --memory-limit '",
                    &value,
                    b"': value is out of range",
                )
            })
        }
    }
}
pub(crate) fn parse_native_blastp_config_os(
    args: &[OsString],
    invocation: &[u8],
) -> io::Result<BlastpConfig> {
    let options = Options::new(args)?;
    blastp_config(&options, invocation)
}
/// Project an already Source-validated workflow configuration into its internal
/// protein edge search. Source mutates these fields on the existing Config;
/// internal overrides are not a second user invocation or workflow permission.
pub(crate) fn parse_native_workflow_edge_config_os(
    args: &[OsString],
    invocation: &[u8],
) -> io::Result<BlastpConfig> {
    let mut options = Options::new(args)?;
    options.command = "blastp".to_owned();
    options
        .values
        .insert("outfmt", vec![OsString::from("edge")]);
    blastp_config(&options, invocation)
}
fn blastp_config(options: &Options, invocation: &[u8]) -> io::Result<BlastpConfig> {
    let dna = options.command == "blastn";
    let (output, outfmt) = output_options(options)?;
    let buffer = options.number("file-buffer-size", 67_108_864i64)?;
    let file_buffer_size = preflight_file_buffer_size(buffer)
        .ok_or_else(|| crate::util::io::BadAllocError::default().into_io())?;
    let global_ranking_targets = options.number("global-ranking", 0i64)?;
    let frame_shift = options.number("frameshift", 0i32)?;
    let extension_mode = if dna { None } else { options.extension()? };
    let toppercent = if options.has("top") {
        Some(options.number("top", 0.0)?)
    } else if options.has("long-reads") {
        Some(10.0)
    } else {
        None
    };
    let index_chunks = options.number("index-chunks", 0u32)?;
    let algorithm = options.algorithm()?;
    let cbs = options.number("comp-based-stats", 1u32)?;
    let comp_based_stats = if dna {
        CbsMode::parse(&cbs.to_string())
    } else {
        CbsMode::from_search_code(cbs, options.command == "blastx", global_ranking_targets)
            .map_err(io::Error::other)?
    };
    let target_indexed = options.has("target-indexed");
    if target_indexed && index_chunks != 1 {
        print_banner();
        return Err(io::Error::other("--target-indexed requires -c1."));
    }
    Ok(BlastpConfig {
        query_files: options.paths("query"),
        database: options.path("db").unwrap_or_default(),
        output,
        matrix: options
            .raw("matrix")?
            .unwrap_or_else(|| b"blosum62".to_vec()),
        gap_open: options.number("gapopen", -1)?,
        gap_extend: options.number("gapextend", -1)?,
        max_evalue: options.number("evalue", 0.001)?,
        min_bit_score: options.number("min-score", 0.0)?,
        max_target_seqs: match options.number("max-target-seqs", 25i64)? {
            0 => i64::MAX,
            value => value,
        },
        max_hsps: options.number("max-hsps", 1u32)?,
        report_unaligned: options.number("unal", u32::MAX)? as i32,
        include_lineage: options.has("include-lineage"),
        output_header: options.raw_vector("header")?,
        compression: options.raw("compress")?.unwrap_or_default(),
        file_buffer_size,
        sam: crate::commands::blastp::SamOptions {
            invocation: invocation.to_vec(),
            all_subject_titles: options.has("salltitles"),
            query_len: options.has("sam-query-len"),
        },
        archive: crate::commands::blastp::ArchiveOptions {
            all_subject_titles: options.has("salltitles"),
            all_subject_ids: options.has("sallseqid"),
        },
        tuning: parse_tuning(options, true)?,
        scoring_matrix: None,
        auxiliary_output: crate::commands::query_output::AuxiliaryOutputOptions {
            unaligned: options.path("un").unwrap_or_default(),
            aligned: options.path("al").unwrap_or_default(),
            unaligned_format: options.raw("unfmt")?.unwrap_or_else(|| b"fasta".to_vec()),
            aligned_format: options.raw("alfmt")?.unwrap_or_else(|| b"fasta".to_vec()),
        },
        taxonlist: options.raw("taxonlist")?.unwrap_or_default(),
        taxon_exclude: options.raw("taxon-exclude")?.unwrap_or_default(),
        seqidlist: options.path("seqidlist").unwrap_or_default(),
        skip_missing_seqids: options.has("skip-missing-seqids"),
        range_culling: options.has("range-culling") || options.has("long-reads"),
        range_cover: options.number("range-cover", 50.0)?,
        culling_overlap: options.number("culling-overlap", 50.0)?,
        taxon_k: options.number("taxon-k", 0u64)?,
        query_or_subject_cover: options.number("query-or-subject-cover", 0.0)?,
        reference_taxon_rank_ids: None,
        ext_chunk_size: options.number("ext-chunk-size", 0)?,
        toppercent,
        global_ranking_targets,
        frame_shift,
        extension_mode,
        iterate: options.raw_vector("iterate")?,
        round_linear_target: false,
        round_query_skip: None,
        min_id: options.number("id", 0.0)?,
        approx_min_id: options.number("approx-id", 0.0)?,
        soft_masking: if dna {
            crate::masking::MaskingAlgo::None
        } else {
            soft_masking(options)?
        },
        threads: options.number("threads", 0i32)?,
        index_chunks,
        block_size: if options.has("block-size") {
            Some(options.number("block-size", 0.0)?)
        } else {
            None
        },
        query_indexed: algorithm.as_deref() == Some("1"),
        contiguous_seed_mode: algorithm.as_deref() == Some("ctg"),
        swipe_all: options.has("swipe"),
        outfmt,
        sensitivity: options.sensitivity()?,
        masking: if dna {
            crate::masking::MaskingMode::None
        } else {
            masking(options)?
        },
        motif_masking: options.raw("motif-masking")?.unwrap_or_default(),
        min_query_len: options.number("min-query-len", 0i32)?,
        query_cover: options.number("query-cover", 0.0)?,
        subject_cover: options.number("subject-cover", 0.0)?,
        comp_based_stats,
        no_self_hits: options.has("no-self-hits"),
        ungapped_xdrop_bits: options.number("xdrop", 12.3)?,
        memory_limit: Some(memory_limit(options)?),
        tmpdir: options.path("tmpdir").unwrap_or_default(),
        target_indexed,
        translated_query_layout: None,
    })
}

fn blastx_config(options: &Options, invocation: &[u8]) -> io::Result<BlastxConfig> {
    let protein = blastp_config(options, invocation)?;
    let query_gencode = options.number("query-gencode", 1u32)?;
    crate::basic::translate::Translator::init(query_gencode).map_err(io::Error::other)?;
    let strand = options.raw("strand")?.unwrap_or_else(|| b"both".to_vec());
    if !matches!(strand.as_slice(), b"both" | b"plus" | b"minus") {
        return Err(io::Error::other("Invalid value for parameter --strand"));
    }
    Ok(BlastxConfig {
        query_files: protein.query_files,
        database: protein.database,
        output: protein.output,
        matrix: protein.matrix,
        gap_open: protein.gap_open,
        gap_extend: protein.gap_extend,
        max_evalue: protein.max_evalue,
        min_bit_score: protein.min_bit_score,
        max_target_seqs: protein.max_target_seqs,
        max_hsps: protein.max_hsps,
        report_unaligned: protein.report_unaligned,
        include_lineage: protein.include_lineage,
        output_header: protein.output_header,
        compression: protein.compression,
        file_buffer_size: protein.file_buffer_size,
        sam: protein.sam,
        archive: protein.archive,
        tuning: protein.tuning,
        scoring_matrix: protein.scoring_matrix,
        auxiliary_output: protein.auxiliary_output,
        taxonlist: protein.taxonlist,
        taxon_exclude: protein.taxon_exclude,
        seqidlist: protein.seqidlist,
        skip_missing_seqids: protein.skip_missing_seqids,
        range_culling: protein.range_culling,
        range_cover: protein.range_cover,
        culling_overlap: protein.culling_overlap,
        taxon_k: protein.taxon_k,
        query_or_subject_cover: protein.query_or_subject_cover,
        reference_taxon_rank_ids: protein.reference_taxon_rank_ids,
        ext_chunk_size: protein.ext_chunk_size,
        toppercent: protein.toppercent,
        global_ranking_targets: protein.global_ranking_targets,
        frame_shift: protein.frame_shift,
        extension_mode: protein.extension_mode,
        iterate: protein.iterate,
        min_id: protein.min_id,
        approx_min_id: protein.approx_min_id,
        soft_masking: protein.soft_masking,
        threads: options.number("threads", 0i32)?,
        index_chunks: protein.index_chunks,
        block_size: protein.block_size,
        query_indexed: protein.query_indexed,
        contiguous_seed_mode: protein.contiguous_seed_mode,
        min_query_len: protein.min_query_len,
        swipe_all: protein.swipe_all,
        outfmt: protein.outfmt,
        sensitivity: protein.sensitivity,
        query_gencode: query_gencode,
        strand: strand,
        min_orf: if options.has("min-orf") {
            Some(options.number("min-orf", 0u32)?)
        } else {
            None
        },
        masking: protein.masking,
        motif_masking: protein.motif_masking,
        query_cover: protein.query_cover,
        subject_cover: protein.subject_cover,
        comp_based_stats: protein.comp_based_stats,
        no_self_hits: protein.no_self_hits,
        ungapped_xdrop_bits: protein.ungapped_xdrop_bits,
        memory_limit: protein.memory_limit,
        tmpdir: protein.tmpdir,
        target_indexed: protein.target_indexed,
    })
}

pub(super) fn view_cli(options: &Options, argc: usize) -> io::Result<NativeViewCli> {
    let mut parsed = NativeViewCli::default();
    parsed.daa_file = options.path("daa").unwrap_or_default();
    parsed.output = options.path("out").unwrap_or_default();
    parsed.outfmt = options.raw_vector("outfmt")?.unwrap_or_default();
    parsed.header = options.raw_vector("header")?.unwrap_or_default();
    parsed.max_target_seqs = match options.number("max-target-seqs", 25i64)? {
        0 => i64::MAX,
        n => n,
    };
    parsed.toppercent = if options.has("top") {
        Some(options.number("top", 0.0)?)
    } else {
        None
    };
    if let Some(top) = parsed.toppercent {
        if !(0.0..=100.0).contains(&top) {
            return Err(io::Error::other(
                "Allowed value range for --top is between 0.0 and 100.0",
            ));
        }
        if options.has("max-target-seqs") {
            return Err(io::Error::other(
                "--top and -k/--max-target-seqs are mutually exclusive.",
            ));
        }
        if top == 100.0 {
            parsed.toppercent = None;
            parsed.max_target_seqs = i64::MAX;
        }
    }
    if options.has("header") && parsed.header.is_empty() {
        parsed.header.push(Vec::new());
    }
    parsed.forward_only = options.has("forwardonly");
    parsed.threads = options.number("threads", 0i32)?.max(0) as usize;
    parsed.qnum_offset = options.number("qnum-offset", 0u64)?;
    parsed.snum_offset = options.number("snum-offset", 0u64)?;
    parsed.file_buffer_size =
        preflight_file_buffer_size(options.number("file-buffer-size", 67_108_864i64)?)
            .ok_or_else(|| crate::util::io::BadAllocError::default().into_io())?;
    parsed.verbosity = if options.has("log") {
        3
    } else if options.has("quiet") {
        0
    } else if options.has("verbose") {
        2
    } else if parsed.output.as_os_str().is_empty() && argc != 2 {
        0
    } else {
        1
    };
    let runtime = &mut parsed.runtime;
    runtime.source_extra = cfg!(feature = "extra");
    if cfg!(feature = "extra") {
        runtime.compression = options.raw("compress")?.unwrap_or_default();
        runtime.database = options.path("db").unwrap_or_default();
        runtime.taxonomy_names = options.path("taxonnames").unwrap_or_default();
        runtime.xml_blord_format = options.has("xml-blord-format");
        runtime.no_parse_seqids = options.has("no-parse-seqids");
        runtime.extension = options.raw("ext")?.unwrap_or_default();
        runtime.multiprocessing = options.has("multiprocessing");
        runtime.query_range_culling = options.has("range-culling");
        runtime.frame_shift = options.number("frameshift", 0i32)?;
        runtime.global_ranking = options.number("global-ranking", 0i64)?;
        runtime.masking = options.raw("masking")?.unwrap_or_default();
        runtime.top_present = options.has("top");
        runtime.include_lineage = options.has("include-lineage");
        runtime.min_bit_score = options.number("min-score", 0.0)?;
        let input = &mut runtime.search_input;
        input.self_search = options.has("self");
        input.target_indexed = options.has("target-indexed");
        input.iterate = options.raw_vector("iterate")?;
        input.lin_stage1_query = options.has("lin-stage1");
        input.lin_stage1_target = options.has("linsearch");
        input.unaligned = options.path("un").unwrap_or_default();
        input.aligned_file = options.path("al").unwrap_or_default();
        input.taxonlist = options.raw("taxonlist")?.unwrap_or_default();
        input.taxon_exclude = options.raw("taxon-exclude")?.unwrap_or_default();
        input.freq_masking = options.has("freq-masking");
        input.seed_cut = options.number("seed-cut", 0.0)?;
        input.freq_sd = options.number("freq-sd", 0.0)?;
        input.minimizer_window = options.number("minimizer-window", 0i32)?;
        input.min_length_ratio = options.number("min-len-ratio", 0.0)?;
        input.query_cover = options.number("query-cover", 0.0)?;
        input.subject_cover = options.number("subject-cover", 0.0)?;
    }
    Ok(parsed)
}

pub(super) fn prepare(
    original_args: &[OsString],
    session: &mut super::super::source_config::Session,
) -> NativeResult<super::native_main::Preparation> {
    use super::native_main::{NativeCommand, Preparation};
    use crate::run::main::RunCommand;
    if original_args.len() < 2 {
        print_usage();
        return Ok(Preparation::Finished(0));
    }
    crate::commands::argument_registry::validate_os(original_args)?;
    if original_args
        .iter()
        .skip(2)
        .any(|v| v == OsStr::new("--legacy"))
    {
        return Err(NativeDispatchError::Untranslated);
    }
    let mut normalized = original_args.to_vec();
    if let Some(command) = crate::commands::argument_registry::command_name_os(original_args)? {
        normalized[1] = OsString::from(command);
    }
    let args = normalized.as_slice();
    let options = Options::new(args)?;
    if matches!(options.command.as_str(), "blastp" | "blastx") && options.has("help") {
        print_usage();
        return Ok(Preparation::Finished(0));
    }
    let startup = super::super::source_config::prepare_os_in(original_args, session)?;
    if options.command == "blastp" && options.number("frameshift", 0i32)? != 0 {
        return Err(io::Error::other(
            "Frameshift alignments are only supported for translated searches.",
        )
        .into());
    }
    if matches!(options.command.as_str(), "blastp" | "blastx")
        && options.number("frameshift", 0i32)? > 0
        && options.extension()? == Some(crate::search::sensitivity::ExtensionMode::Full)
    {
        return Err(io::Error::other(
            "Frameshift alignment does not support full matrix extension.",
        )
        .into());
    }

    #[cfg(feature = "extra")]
    if let Some(extra) = super::super::extra::prepare_os(args, original_args, startup.clone())? {
        let kind = extra.command();
        return super::native_main::ready_os(NativeCommand::Extra(extra), kind, args, startup);
    }
    let mut invocation = Vec::new();
    for (i, arg) in original_args.iter().enumerate() {
        if i > 0 {
            invocation.push(b' ');
        }
        invocation
            .extend_from_slice(&crate::util::io::source_os_bytes(arg).map_err(io::Error::other)?);
    }
    let workflow = |kind| NativeCommand::Workflow {
        kind,
        args: original_args.to_vec(),
        startup: startup.clone(),
    };
    let (command, kind) = match options.command.as_str() {
        "help" => (NativeCommand::Help, RunCommand::Help),
        "version" => (NativeCommand::Version, RunCommand::Version),
        "dbinfo" => {
            let database = options.path("db").unwrap_or_default();
            if database.as_os_str().is_empty() {
                return Err(io::Error::other("Missing parameter: database file (--db/-d)").into());
            }
            (
                NativeCommand::DbInfo {
                    database,
                    file_buffer_size: options.number("file-buffer-size", 67_108_864i64)? as usize,
                },
                RunCommand::DbInfo,
            )
        }
        "makedb" => {
            let database = options.path("db").unwrap_or_default();
            if database.as_os_str().is_empty() {
                return Err(io::Error::other("Missing parameter: database file (--db/-d)").into());
            }
            if options.number("block-size", 0.0f64)? != 0.0 {
                return Err(io::Error::other("Invalid option: --block-size/-b. Block size is set for the alignment commands.").into());
            }
            let sequence_type = match options
                .enum_text("dbtype", &["prot", "nucl"], "nucl, prot")?
                .as_deref()
                .unwrap_or("prot")
            {
                "nucl" => crate::basic::value::SequenceType::Nucleotide,
                _ => crate::basic::value::SequenceType::AminoAcid,
            };
            (
                NativeCommand::MakeDb {
                    inputs: options.paths("in"),
                    database,
                    threads: startup.config.threads,
                    options: crate::commands::makedb::MakedbOptions {
                        taxonmap: options.path("taxonmap"),
                        taxonnodes: options.path("taxonnodes"),
                        taxonnames: options.path("taxonnames"),
                        no_parse_seqids: options.has("no-parse-seqids"),
                    },
                    runtime: crate::commands::makedb::RuntimeOptions {
                        temporary_directory: options.path("tmpdir").unwrap_or_default(),
                        file_buffer_size: Some(4 * 1024 * 1024),
                        sequence_type,
                        no_auto_append: options.has("no-auto-append"),
                        mask: options.raw("masking")?.as_deref() != Some(b"0"),
                        matrix: options
                            .raw("matrix")?
                            .unwrap_or_else(|| b"blosum62".to_vec()),
                        gap_open: options.number("gapopen", -1)?,
                        gap_extend: options.number("gapextend", -1)?,
                        cbs: options.number("comp-based-stats", 1u32)?,
                        tuning: parse_tuning(&options, true)?,
                    },
                },
                RunCommand::MakeDb,
            )
        }
        "blastp" => (
            NativeCommand::BlastP(Box::new(blastp_config(&options, &invocation)?)),
            RunCommand::BlastP,
        ),
        "blastx" => (
            NativeCommand::BlastX(Box::new(blastx_config(&options, &invocation)?)),
            RunCommand::BlastX,
        ),
        "getseq" => (workflow(RunCommand::GetSeq), RunCommand::GetSeq),
        "makeidx" => {
            let mut config = crate::data::index::MakeIndexConfig::new(
                options.path("db").unwrap_or_default(),
                options.sensitivity()?,
            );
            config.shape_mask = options.raw_vector("shape-mask")?.unwrap_or_default();
            config.shapes = options.number("shapes", 0u32)?;
            config.threads = startup.config.threads.max(0) as usize;
            config.file_buffer_size = options.number("file-buffer-size", 67_108_864i64)? as usize;
            (NativeCommand::MakeIndex(config), RunCommand::MakeIdx)
        }
        "view" => {
            let parsed = view_cli(&options, args.len())?;
            (
                NativeCommand::View(
                    ViewConfig {
                        daa_file: parsed.daa_file,
                        output: parsed.output,
                        outfmt: parsed.outfmt,
                        max_target_seqs: parsed.max_target_seqs,
                        toppercent: parsed.toppercent,
                        forward_only: parsed.forward_only,
                        report_unaligned: !cfg!(feature = "extra")
                            && options.number("unal", 0i32)? != 0,
                        sam_qlen_field: options.has("sam-query-len"),
                        invocation: invocation.clone(),
                        header: parsed.header,
                        threads: parsed.threads,
                        qnum_offset: parsed.qnum_offset,
                        snum_offset: parsed.snum_offset,
                        verbosity: parsed.verbosity,
                        file_buffer_size: parsed.file_buffer_size,
                    },
                    parsed.runtime,
                ),
                RunCommand::View,
            )
        }
        "cluster" | "linclust" | "deepclust" => {
            let kind = match options.command.as_str() {
                "linclust" => RunCommand::LinClust,
                "deepclust" => RunCommand::DeepClust,
                _ => RunCommand::Cluster,
            };
            #[cfg(feature = "mcl")]
            if options.raw("cluster-algo")?.as_deref() == Some(b"mcl") {
                return super::native_main::ready_os(
                    NativeCommand::Mcl(
                        Box::new(blastp_config(&options, &invocation)?),
                        mcl_options(&options)?,
                    ),
                    kind,
                    args,
                    startup,
                );
            }
            (
                NativeCommand::Cluster(Box::new(cluster_config(&options, args)?)),
                kind,
            )
        }
        "merge-daa" => {
            let buffer =
                preflight_file_buffer_size(options.number("file-buffer-size", 67_108_864i64)?)
                    .ok_or_else(|| crate::util::io::BadAllocError::default().into_io())?;
            (
                NativeCommand::MergeDaa {
                    inputs: options.paths("in"),
                    output: options.path("out").unwrap_or_default(),
                    buffer,
                },
                RunCommand::MergeDaa,
            )
        }
        "roc" => (NativeCommand::Roc, RunCommand::Roc),
        "prepdb" => (NativeCommand::PrepDb, RunCommand::PrepDb),
        "greedy-vertex-cover" => (
            workflow(RunCommand::GreedyVertexCover),
            RunCommand::GreedyVertexCover,
        ),
        "realign" => (
            workflow(RunCommand::ClusterRealign),
            RunCommand::ClusterRealign,
        ),
        "reassign" => (
            workflow(RunCommand::ClusterReassign),
            RunCommand::ClusterReassign,
        ),
        "recluster" => (workflow(RunCommand::Recluster), RunCommand::Recluster),
        "benchmark" => (workflow(RunCommand::Benchmark), RunCommand::Benchmark),
        "test" => (
            workflow(RunCommand::RegressionTest),
            RunCommand::RegressionTest,
        ),
        _ => return Err(NativeDispatchError::Untranslated),
    };
    super::native_main::ready_os(command, kind, args, startup)
}

#[cfg(feature = "mcl")]
pub(super) fn mcl_options(options: &Options) -> io::Result<super::super::mcl_cmd::MclOptions> {
    Ok(super::super::mcl_cmd::MclOptions {
        similarity: options.raw("cluster-similarity")?.unwrap_or_default(),
        threshold: if options.has("cluster-threshold") {
            Some(options.number("cluster-threshold", 0.0)?)
        } else {
            None
        },
        graph_file: options.path("cluster-graph-file").unwrap_or_default(),
        restart: options.has("cluster-restart"),
        statistics: options.has("mcl-stats"),
        component_chunk_size: options.number("mcl-chunk-size", 1u32)?,
        max_target_seqs: if options.has("max-target-seqs") {
            Some(options.number("max-target-seqs", 0i64)?)
        } else {
            None
        },
        config: crate::mcl::mcl::MclConfig {
            expansion: options.number("mcl-expansion", 2u32)? as f32,
            inflation: options.number("mcl-inflation", 2.0f64)? as f32,
            max_iter: options.number("mcl-max-iterations", 100u32)?,
            sparsity_switch: options.number("mcl-sparsity-switch", 0.8f64)? as f32,
            symmetric: !options.has("mcl-nonsymmetric"),
        },
    })
}
pub(super) fn cluster_config(options: &Options, args: &[OsString]) -> io::Result<ClusterConfig> {
    let workflow = match options.command.as_str() {
        "linclust" => ClusterWorkflow::LinClust,
        "deepclust" => ClusterWorkflow::DeepClust,
        _ => ClusterWorkflow::Cascaded,
    };
    let parallel = options.path("parallel-tmpdir").unwrap_or_default();
    let external = workflow == ClusterWorkflow::LinClust && !parallel.as_os_str().is_empty();
    if !parallel.as_os_str().is_empty() && workflow != ClusterWorkflow::LinClust {
        return Err(io::Error::other(
            "Option is not permitted for this workflow: --parallel-tmpdir",
        ));
    }
    if parallel.as_os_str().is_empty()
        && options
            .path("aln-out")
            .is_some_and(|p| !p.as_os_str().is_empty())
    {
        return Err(io::Error::other(
            "Alignment output to file is not supported without --parallel-tmpdir.",
        ));
    }
    let graph_algo = match options.raw("graph-algo")?.as_deref().unwrap_or(b"gvc") {
        b"gvc" => GraphAlgo::GreedyVertexCover,
        b"len" => GraphAlgo::LengthSorted,
        value => return Err(raw_error(b"Invalid value for --graph-algo: ", value, b"")),
    };
    let header = options.raw_vector("header")?;
    let simple_header = match header.as_deref() {
        None => false,
        Some([]) => true,
        Some([one]) if one == b"0" => false,
        Some([one]) if one == b"simple" => true,
        Some([one]) if one == b"verbose" => {
            return Err(io::Error::other(
                "Verbose header format is not supported for cluster workflow.",
            ))
        }
        Some(values) => {
            let value = values.join(&b' ');
            return Err(raw_error(b"Invalid header format: ", &value, b""));
        }
    };
    let mask_value = options.raw("masking")?.unwrap_or_else(|| b"0".to_vec());
    let masking = if external {
        crate::masking::MaskingMode::None
    } else {
        let mut lower = mask_value.clone();
        lower.make_ascii_lowercase();
        match lower.as_slice() {
            b"0" | b"none" => crate::masking::MaskingMode::None,
            b"1" | b"tantan" => crate::masking::MaskingMode::Tantan,
            b"seg" => crate::masking::MaskingMode::BlastSeg,
            _ => {
                return Err(raw_error(
                    b"Invalid value for string field: ",
                    &mask_value,
                    b". Permitted values: none, tantan, seg",
                ))
            }
        }
    };
    let soft_value = options
        .raw("soft-masking")?
        .unwrap_or_else(|| b"tantan".to_vec());
    let soft_masking = match soft_value.as_slice() {
        b"0" | b"none" => crate::masking::MaskingAlgo::None,
        b"seg" => crate::masking::MaskingAlgo::Seg,
        b"tantan" => crate::masking::MaskingAlgo::Tantan,
        _ if external => crate::masking::MaskingAlgo::None,
        _ => {
            return Err(raw_error(
                b"Invalid value for string field: ",
                &soft_value,
                b". Permitted values: 0, none, seg, tantan",
            ))
        }
    };
    let comp_based_stats = options.number("comp-based-stats", 1i32)?;
    if !(0..=5).contains(&comp_based_stats) {
        return Err(io::Error::other(
            "Invalid value for --comp-based-stats. Permitted values: 0, 1, 2, 3, 4, 5.",
        ));
    }
    let extension_mode = options.raw("ext")?.unwrap_or_default();
    if !external && !extension_mode.is_empty() {
        options.extension()?;
    }
    let motif_masking = options.raw("motif-masking")?.unwrap_or_default();
    if !external && !matches!(motif_masking.as_slice(), b"" | b"0" | b"1") {
        return Err(io::Error::other(
            "Permitted values for --motif-masking: 0, 1",
        ));
    }
    let mut tuning = parse_tuning(options, true)?;
    tuning.setup.parallel_tmpdir = parallel.clone();
    tuning.setup.soft_masking = soft_value;
    if !tuning.matrix_file.as_os_str().is_empty() && comp_based_stats > 1 {
        return Err(io::Error::other("This value for --comp-based-stats is not supported when using a custom scoring matrix."));
    }
    let file_buffer_size =
        preflight_file_buffer_size(options.number("file-buffer-size", 67_108_864i64)?)
            .ok_or_else(|| crate::util::io::BadAllocError::default().into_io())?;
    Ok(ClusterConfig {
        search_args: Some(args.to_vec()),
        tuning,
        file_buffer_size,
        tmpdir: options.path("tmpdir").unwrap_or_default(),
        workflow,
        parallel_tmpdir: parallel,
        linclust_chunk_size: crate::util::string::interpret_number_bytes(
            &options
                .raw("linclust-chunk-size")?
                .unwrap_or_else(|| b"10G".to_vec()),
        )? as u64,
        database: options.path("db").unwrap_or_default(),
        output: options.path("out").unwrap_or_default(),
        threads: options.number("threads", 0i32)?,
        matrix: options
            .raw("matrix")?
            .unwrap_or_else(|| b"blosum62".to_vec()),
        gap_open: options.number("gapopen", -1i32)?,
        gap_extend: options.number("gapextend", -1i32)?,
        ext_chunk_size: options.number("ext-chunk-size", 0i64)?,
        max_evalue: options.number("evalue", 0.001)?,
        comp_based_stats,
        masking,
        soft_masking,
        motif_masking,
        extension_mode,
        member_cover: if options.has("member-cover") {
            Some(options.number("member-cover", 0.0)?)
        } else {
            None
        },
        mutual_cover: if options.has("mutual-cover") {
            Some(options.number("mutual-cover", 0.0)?)
        } else {
            None
        },
        approx_id: if options.has("approx-id") {
            Some(options.number("approx-id", 0.0)?)
        } else {
            None
        },
        cluster_steps: options.raw_vector("cluster-steps")?.unwrap_or_default(),
        connected_component_depth: options
            .raw_vector("connected-component-depth")?
            .unwrap_or_default(),
        round_coverage: options.raw_vector("round-coverage")?.unwrap_or_default(),
        round_approx_id: options.raw_vector("round-approx-id")?.unwrap_or_default(),
        kmer_ranking: options.has("kmer-ranking"),
        no_reassign: options.has("no-reassign"),
        weighted_gvc: options.has("weighted-gvc"),
        strict_gvc: options.has("strict-gvc"),
        graph_algo,
        memory_limit: memory_limit(options)?,
        oid_output: options.has("oid-output"),
        simple_header,
        alignment_output: options.path("aln-out"),
    })
}

#[cfg(test)]
mod raw_tests {
    use super::*;
    #[cfg(unix)]
    #[test]
    fn native_factory_keeps_os_paths_and_late_raw_controls() {
        use std::os::unix::ffi::{OsStrExt, OsStringExt};
        let args = vec![
            OsString::from("diamond"),
            OsString::from("blastp"),
            OsString::from("--query"),
            OsString::from_vec(b"query\xff.fa".to_vec()),
            OsString::from("--db"),
            OsString::from_vec(b"db\xfe".to_vec()),
            OsString::from("--matrix"),
            OsString::from_vec(b"name\xfd".to_vec()),
            OsString::from("--taxonlist"),
            OsString::from_vec(b"7\xfc".to_vec()),
            OsString::from("--outfmt"),
            OsString::from_vec(b"format\xfb".to_vec()),
            OsString::from("--threads"),
            OsString::from_vec(b"3\xfa".to_vec()),
        ];
        // This factory owns parsed configuration. Source constructors and
        // runtime owners validate matrix/format controls at their own phase.
        let config = parse_native_blastp_config_os(&args, b"original\xf9 invocation").unwrap();
        assert_eq!(
            config.query_files[0].as_os_str().as_bytes(),
            b"query\xff.fa"
        );
        assert_eq!(config.database.as_os_str().as_bytes(), b"db\xfe");
        assert_eq!(config.matrix, b"name\xfd");
        assert_eq!(config.taxonlist, b"7\xfc");
        assert_eq!(config.outfmt, vec![b"format\xfb".to_vec()]);
        assert_eq!(config.threads, 3);
        assert_eq!(config.sam.invocation, b"original\xf9 invocation");
    }
    #[test]
    fn workflow_edge_projection_validates_users_before_internal_overrides() {
        for workflow in ["reassign", "recluster", "cluster"] {
            let args =
                ["diamond", workflow, "--db", "original-db", "--threads", "1"].map(OsString::from);
            let invocation = b"original workflow invocation";
            let config = parse_native_workflow_edge_config_os(&args, invocation).unwrap();
            assert_eq!(config.outfmt, vec![b"edge".to_vec()]);
            assert_eq!(config.database, PathBuf::from("original-db"));
            assert_eq!(config.threads, 1);
            assert_eq!(config.sam.invocation, invocation);

            let mut user_format = args.to_vec();
            user_format.extend([OsString::from("--outfmt"), OsString::from("xml")]);
            let result = parse_native_workflow_edge_config_os(&user_format, invocation);
            if cfg!(feature = "extra") {
                assert_eq!(result.unwrap().outfmt, vec![b"edge".to_vec()]);
            } else {
                assert_eq!(
                    result.err().unwrap().to_string(),
                    "Option is not permitted for this workflow: outfmt"
                );
            }
        }
    }
    #[test]
    fn final_disabled_short_id_shadowing_is_used_by_real_factories() {
        let args = ["diamond", "blastp", "-w10"].map(OsString::from);
        let result = parse_native_search_tuning_os(&args);
        if cfg!(feature = "extra") {
            result.unwrap();
        } else {
            let error = result.unwrap_err();
            assert_eq!(
                crate::util::io::RawDiagnosticError::find(&error),
                Some(b"Invalid option: w".as_slice())
            );
        }
    }
}
