//! Testable translation of `diamond/src/basic/config.cpp`.

use crate::config::Sensitivity;
use crate::util::io::Compressor;
use std::ffi::OsString;
use std::path::{Path, PathBuf};

/// Source `EnumTraits<Sensitivity>::from_string`, including the sorted map's
/// permitted-value diagnostic. Keep workflow parsers on this shared registry.
pub fn sensitivity_from_string(value: &str) -> Result<Sensitivity, String> {
    sensitivity_from_bytes(value.as_bytes()).map_err(|error| error.to_string())
}

pub fn sensitivity_from_bytes(value: &[u8]) -> std::io::Result<Sensitivity> {
    match value {
        b"default" => Ok(Sensitivity::Default),
        b"fast" => Ok(Sensitivity::Fast),
        b"faster" => Ok(Sensitivity::Faster),
        b"linclust-20" => Ok(Sensitivity::Linclust20),
        b"linclust-40" => Ok(Sensitivity::Linclust40),
        b"mid-sensitive" => Ok(Sensitivity::MidSensitive),
        b"more-sensitive" => Ok(Sensitivity::MoreSensitive),
        b"sensitive" => Ok(Sensitivity::Sensitive),
        b"shapes-30x10" => Ok(Sensitivity::Shapes30x10),
        b"shapes-6x10" => Ok(Sensitivity::Shapes6x10),
        b"ultra-sensitive" => Ok(Sensitivity::UltraSensitive),
        b"very-sensitive" => Ok(Sensitivity::VerySensitive),
        _ => {
            let mut message = b"Invalid value for string field: ".to_vec();
            message.extend_from_slice(value);
            message.extend_from_slice(b". Permitted values: default, fast, faster, linclust-20, linclust-40, mid-sensitive, more-sensitive, sensitive, shapes-30x10, shapes-6x10, ultra-sensitive, very-sensitive");
            Err(crate::util::io::RawDiagnosticError(message).into_io())
        }
    }
}

pub const REGISTERED_COMMANDS: &str = "makedb prepdb blastp blastx cluster linclust realign recluster reassign view merge-daa help version getseq dbinfo test makeidx greedy-vertex-cover roc benchmark deepclust random-seqs sort dbstat mask fastq2fasta info seed-stat smith-waterman simulate-seqs split upgma upgmamc reverse compute-medoids mutate roc-id find-shapes hashseqs listseeds blastn wc cut model-seqs make-seed-table";

/// Exact registration order. DNA-only names occupy one contiguous portion and
/// are filtered by [`registered_options`].
pub const REGISTERED_OPTIONS: &str = "threads verbose log quiet tmpdir db out header in taxonmap taxonnodes taxonnames comp-based-stats masking soft-masking no-block-size-limit gapopen gapextend matrix custom-matrix evalue motif-masking approx-id ext max-target-seqs top faster fast mid-sensitive linclust-20 shapes-6x10 shapes-30x10 sensitive more-sensitive very-sensitive ultra-sensitive shapes query strand un al unfmt alfmt unal max-hsps range-culling compress min-score id query-cover subject-cover swipe iterate global-ranking block-size index-chunks frameshift long-reads query-gencode salltitles sallseqid no-self-hits taxonlist taxon-exclude seqidlist skip-missing-seqids outfmt qnum-offset snum-offset include-lineage cluster-steps kmer-ranking round-coverage round-approx-id aln-out memory-limit member-cover mutual-cover connected-component-depth no-reassign centroid-out edges edge-format symmetric clusters file-buffer-size no-unlink ignore-warnings no-parse-seqids parallel-tmpdir bin ext-chunk-size no-ranking dbsize no-auto-append tantan-minMaskProb oid-output swipe-task-size anchored-swipe query-match-distance-threshold length-ratio-threshold cbs-angle linclust-banded-ext linclust-chunk-size hit-membuf algo min-orf min-query-len load-threads minichunk seed-cut freq-masking freq-sd sketch-size tile-size id2 linsearch lin-stage1 lin-combo xdrop ungapped-evalue ungapped-evalue-short short-query-ungapped-bitscore gapped-filter-evalue band shape-mask multiprocessing mp-init mp-recover mp-query-chunk culling-overlap taxon-k range-cover xml-blord-format sam-query-len stop-match-score target-indexed unaligned-targets cut-bar check-multi-target roc-file family-map family-map-query query-parallel-limit log-evalue-scale bootstrap heartbeat mp-self zdrop repetition-cutoff extension chaining-out align-long-reads best-hsp-only chain-pen-gap-scale chain-pen-skip-scale penalty reward chain-align-cutoff min-chain-score max-overlap-extension zdrop-extension zdrop-global band-extension band-global query-or-subject-cover daa forwardonly seq window ungapped-score hit-band hit-score gapped-xdrop rank-ratio2 rank-ratio lambda K match1 match2 seed-freq space-penalty reverse neighborhood-score seed-weight superblock load-balancing log-query log-subject palign score-ratio fetch-size target-fetch-size rank-factor transcript-len-estimate family-counts radix-cluster-buffered join-split-size join-split-key-len radix-bits join-ht-factor sort-join simple-freq freq-treshold use-dataset-field store-query-quality swipe-chunk-size hard-masked cbs-window no-dict upgma-edge-limit tree upgma-dist upgma-input log-extend chaining-maxgap tantan-maxRepeatOffset tantan-ungapped chaining-range-cover no-swipe-realign chaining-maxnodes cutoff-score-8bit min-band-overlap min-realign-overhang ungapped-window gapped-filter-diag-score gapped-filter-window output-hits no-logfile band-bin col-bin self trace-pt-fetch-size short-query-max-len gapped-filter-evalue1 ext-yield full-sw-len relaxed-evalue-factor type raw chaining-len-cap chaining-min-nodes fast-tsv target-parallel-verbosity query-memory memory-intervals seed-hit-density chunk-size-multiplier score-drop-factor left-most-interval ranking-cutoff-bitscore no-forward-fp no-ref-masking target-bias output-fp family-cap cbs-matrix-scale query-count cbs-err-tolerance cbs-it-limit hash_join_swap deque_bucket_size max-swipe-dp no-reextend no-reorder file1 file2 key2 motif-mask-file max-motif-len chaining-stacked-hsp-ratio minimizer-window min_task_trace_pts oid-list bootstrap-block centroid-factor timeout resume target_hard_cap mapany neighbors reassign-overlap reassign-ratio reassign-max add-self-aln weighted-gvc hamming-ext diag-filter-id diag-filter-cov strict-gvc dbtype cluster-similarity cluster-threshold cluster-graph-file cluster-restart mcl-expansion mcl-inflation mcl-chunk-size mcl-max-iterations mcl-sparsity-switch mcl-nonsymmetric mcl-stats cluster-algo approx-backtrace narrow-band-cov narrow-band-factor anchor-window anchor-score classic-band no_8bit_extension no_chaining_merge_hsps graph-algo tsv-read-size min-len-ratio max-indirection promiscuous-seed-ratio";

const DNA: &[&str] = &[
    "zdrop",
    "repetition-cutoff",
    "extension",
    "chaining-out",
    "align-long-reads",
    "best-hsp-only",
    "chain-pen-gap-scale",
    "chain-pen-skip-scale",
    "penalty",
    "reward",
    "chain-align-cutoff",
    "min-chain-score",
    "max-overlap-extension",
    "zdrop-extension",
    "zdrop-global",
    "band-extension",
    "band-global",
];

pub fn registered_commands(extra: bool) -> Vec<&'static str> {
    let v: Vec<_> = REGISTERED_COMMANDS.split_whitespace().collect();
    if extra {
        v
    } else {
        v[..21].to_vec()
    }
}
pub fn registered_options(with_dna: bool) -> Vec<&'static str> {
    REGISTERED_OPTIONS
        .split_whitespace()
        .filter(|x| with_dna || !DNA.contains(x))
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Command {
    MakeDb,
    BlastP,
    BlastX,
    BlastN,
    View,
    DbInfo,
    Version,
    GetSeq,
    Fastq2Fasta,
    RegressionTest,
    Cluster,
    DeepClust,
    LinClust,
    Recluster,
    ClusterReassign,
    ClusterRealign,
    Benchmark,
    Mask,
    ListSeeds,
    ComputeMedoids,
    ModelSeqs,
    MakeSeedTable,
    Other,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Algo {
    Auto,
    DoubleIndexed,
    QueryIndexed,
    CtgSeed,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ConfigInput {
    pub argc: usize,
    pub argv: Vec<OsString>,
    pub check_io: bool,
    pub command: Command,
    pub debug_log: bool,
    pub quiet: bool,
    pub verbose: bool,
    pub database: PathBuf,
    pub output_file: PathBuf,
    pub daa_file: PathBuf,
    pub output_format: Vec<Vec<u8>>,
    pub compression: Vec<u8>,
    pub no_auto_append: bool,
    pub chunk_size: f64,
    pub top_present: bool,
    pub max_target_present: bool,
    pub no_self_hits: bool,
    pub long_reads: bool,
    pub range_culling: bool,
    pub top: Option<f64>,
    pub frame_shift: i32,
    pub global_ranking: i64,
    pub taxon_k: u64,
    pub multiprocessing: bool,
    pub mp_init: bool,
    pub mp_recover: bool,
    pub comp_based_stats: u32,
    pub translated_cbs_supported: bool,
    pub matrix_file: PathBuf,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub sensitivities: Vec<Sensitivity>,
    pub algo: Algo,
    pub query_gencode: u32,
    pub strand: Vec<u8>,
    pub unfmt: Vec<u8>,
    pub alfmt: Vec<u8>,
    pub aligned_file: PathBuf,
    pub query_file: Vec<PathBuf>,
    pub store_query_quality: bool,
    pub target_indexed: bool,
    pub lowmem: u32,
    pub swipe_all: bool,
    pub parallel_tmpdir: PathBuf,
    pub hit_membuf: bool,
    pub tmpdir: PathBuf,
    pub threads: Option<i32>,
    pub hardware_threads: i32,
    pub now: i64,
}
impl Default for ConfigInput {
    fn default() -> Self {
        Self {
            argc: 1,
            argv: vec!["diamond".into()],
            check_io: true,
            command: Command::Other,
            debug_log: false,
            quiet: false,
            verbose: false,
            database: PathBuf::new(),
            output_file: PathBuf::new(),
            daa_file: PathBuf::new(),
            output_format: vec![],
            compression: Vec::new(),
            no_auto_append: false,
            chunk_size: 0.0,
            top_present: false,
            max_target_present: false,
            no_self_hits: false,
            long_reads: false,
            range_culling: false,
            top: None,
            frame_shift: 0,
            global_ranking: 0,
            taxon_k: 0,
            multiprocessing: false,
            mp_init: false,
            mp_recover: false,
            comp_based_stats: 1,
            translated_cbs_supported: true,
            matrix_file: PathBuf::new(),
            gap_open: -1,
            gap_extend: -1,
            sensitivities: vec![],
            algo: Algo::Auto,
            query_gencode: 1,
            strand: b"both".to_vec(),
            unfmt: b"fasta".to_vec(),
            alfmt: b"fasta".to_vec(),
            aligned_file: PathBuf::new(),
            query_file: vec![],
            store_query_quality: false,
            target_indexed: false,
            lowmem: 4,
            swipe_all: false,
            parallel_tmpdir: PathBuf::new(),
            hit_membuf: false,
            tmpdir: PathBuf::new(),
            threads: None,
            hardware_threads: 1,
            now: 0,
        }
    }
}

pub trait ConfigBackend {
    fn invocation(&mut self, args: &[OsString]) -> Result<Vec<u8>, String> {
        invocation_bytes(args)
    }
    fn auxiliary_formats(&mut self, _: &ConfigInput) -> Result<(), String> {
        Ok(())
    }

    fn presentation(&mut self, _: &ConfigInput, _: u8) -> Result<(), String> {
        Ok(())
    }
    fn streams(&mut self, _: &ConfigInput, _: u8) -> Result<(), String> {
        Ok(())
    }
    fn completed(&mut self, _: &Config, _: &ConfigInput) -> Result<(), String> {
        Ok(())
    }
    fn cbs(&mut self, _: u32) -> Result<(), String>;
    fn scoring(&mut self, _: &ConfigInput) -> Result<(), String>;
    fn masking(&mut self) -> Result<(), String>;
    fn enums(&mut self) -> Result<(), String> {
        Ok(())
    }
    fn translator(&mut self, _: u32) -> Result<(), String>;
    fn temp(&mut self, _: &str) -> Result<(), String>;
    fn mkdir(&mut self, _: &str) -> Result<(), String>;
    fn temp_path(&mut self, path: &Path) -> Result<(), String> {
        self.temp(
            path.to_str()
                .ok_or("Non-UTF-8 path requires a path-aware ConfigBackend")?,
        )
    }
    fn mkdir_path(&mut self, path: &Path) -> Result<(), String> {
        self.mkdir(
            path.to_str()
                .ok_or("Non-UTF-8 path requires a path-aware ConfigBackend")?,
        )
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Config {
    pub command: Command,
    pub verbosity: u8,
    pub database: PathBuf,
    pub output_file: PathBuf,
    pub daa_file: PathBuf,
    pub range_culling: bool,
    pub top: Option<f64>,
    pub frame_shift: i32,
    pub sensitivity: Sensitivity,
    pub algo: Algo,
    pub store_query_quality: bool,
    pub tmpdir: PathBuf,
    pub threads: i32,
    pub trace_pt_membuf: bool,
    pub invocation: Vec<u8>,
    pub old_warning: bool,
}

impl Config {
    pub fn set_sens(&mut self, s: Sensitivity) -> Result<(), String> {
        if self.sensitivity != Sensitivity::Default {
            return Err("Sensitivity switches are mutually exclusive.".into());
        }
        self.sensitivity = s;
        Ok(())
    }
    pub fn single_query_file<'a>(&self, i: &'a ConfigInput) -> &'a str {
        self.single_query_path(i)
            .to_str()
            .expect("Non-UTF-8 query path; use single_query_path")
    }
    pub fn single_query_path<'a>(&self, i: &'a ConfigInput) -> &'a Path {
        i.query_file
            .first()
            .map(PathBuf::as_path)
            .unwrap_or_else(|| Path::new(""))
    }
    pub fn compressor(&self, i: &ConfigInput) -> Result<Compressor, String> {
        self.compressor_io(i).map_err(|error| error.to_string())
    }
    pub fn compressor_io(&self, i: &ConfigInput) -> std::io::Result<Compressor> {
        compressor_bytes(&i.compression)
    }
    pub fn new<B: ConfigBackend>(i: &ConfigInput, b: &mut B) -> Result<Self, String> {
        Self::new_with_extra(i, b, false)
    }
    /// Source EXTRA constructor. Its CBS::COUNT is still 6: the last
    /// explicitly assigned enum entry is 5, after the experimental entry 6.
    pub fn new_extra<B: ConfigBackend>(i: &ConfigInput, b: &mut B) -> Result<Self, String> {
        Self::new_with_extra(i, b, true)
    }
    fn new_with_extra<B: ConfigBackend>(
        i: &ConfigInput,
        b: &mut B,
        _extra: bool,
    ) -> Result<Self, String> {
        let verbosity = if i.debug_log {
            3
        } else if i.quiet {
            0
        } else if i.verbose {
            2
        } else if (matches!(
            i.command,
            Command::View | Command::BlastP | Command::BlastX | Command::BlastN
        ) && i.output_file.as_os_str().is_empty()
            && i.argc != 2)
            || matches!(
                i.command,
                Command::Version | Command::GetSeq | Command::Fastq2Fasta | Command::RegressionTest
            )
        {
            0
        } else {
            1
        };
        b.presentation(i, verbosity)?;
        if i.top_present && i.max_target_present {
            return Err("--top and -k/--max-target-seqs are mutually exclusive.".into());
        }
        if i.command == Command::BlastX && i.no_self_hits {
            return Err("--no-self-hits option is not supported in blastx mode.".into());
        }
        let (mut range, mut top, mut frame) = (i.range_culling, i.top, i.frame_shift);
        if i.long_reads {
            range = true;
            if top.is_none() {
                top = Some(10.0)
            }
            if frame == 0 {
                frame = 15
            }
        }
        if i.global_ranking > 0
            && (range
                || i.taxon_k != 0
                || i.multiprocessing
                || i.mp_init
                || i.mp_recover
                || i.comp_based_stats >= 2
                || frame > 0)
        {
            return Err("Global ranking is not supported in this mode.".into());
        }
        if i.comp_based_stats >= 6 {
            return Err(
                "Invalid value for --comp-based-stats. Permitted values: 0, 1, 2, 3, 4, 5.".into(),
            );
        }
        b.cbs(i.comp_based_stats)?;
        if i.command == Command::BlastX && !i.translated_cbs_supported {
            return Err(
                "This mode of composition based stats is not supported for translated searches."
                    .into(),
            );
        }
        let (mut db, mut out, mut daa) = (
            i.database.clone(),
            i.output_file.clone(),
            i.daa_file.clone(),
        );
        if i.check_io {
            if i.command == Command::MakeDb {
                if db.as_os_str().is_empty() {
                    return Err("Missing parameter: database file (--db/-d)".into());
                }
                if i.chunk_size != 0.0 {
                    return Err("Invalid option: --block-size/-b. Block size is set for the alignment commands.".into());
                }
            }
            if matches!(
                i.command,
                Command::BlastP | Command::BlastX | Command::BlastN
            ) {
                if db.as_os_str().is_empty() {
                    return Err("Missing parameter: database file (--db/-d)".into());
                }
                if !daa.as_os_str().is_empty() {
                    if !out.as_os_str().is_empty() {
                        return Err("Options --daa and --out cannot be used together.".into());
                    }
                    if i.output_format.first().is_some_and(|x| x != b"daa") {
                        return Err("Invalid parameter: --daa/-a. Output file is specified with the --out/-o parameter.".into());
                    }
                    out = daa.clone()
                }
                if !daa.as_os_str().is_empty()
                    || matches!(
                        i.output_format.first().map(Vec::as_slice),
                        Some(b"daa" | b"100")
                    )
                {
                    if !i.compression.is_empty() {
                        return Err("Compression is not supported for DAA format.".into());
                    }
                    if !i.no_auto_append {
                        append(&mut out, ".daa")
                    }
                }
            }
            if i.command == Command::DbInfo && db.as_os_str().is_empty() {
                return Err("Missing parameter: database file (--db/-d)".into());
            }
        }
        b.streams(i, verbosity)?;
        if !i.no_auto_append {
            if i.command == Command::MakeDb {
                append(&mut db, ".dmnd")
            }
            if i.command == Command::View {
                append(&mut daa, ".daa")
            }
            if i.compression == b"1" {
                append(&mut out, ".gz")
            }
            if i.compression == b"zstd" {
                append(&mut out, ".zst")
            }
        }
        let mut c = Self {
            command: i.command,
            verbosity,
            database: db,
            output_file: out,
            daa_file: daa,
            range_culling: range,
            top,
            frame_shift: frame,
            sensitivity: Sensitivity::Default,
            algo: i.algo,
            store_query_quality: i.store_query_quality,
            tmpdir: i.tmpdir.clone(),
            threads: i
                .threads
                .filter(|threads| *threads != 0)
                .unwrap_or(i.hardware_threads),
            trace_pt_membuf: i.hit_membuf,
            invocation: b.invocation(&i.argv)?,
            old_warning: i.command != Command::Version && i.now - 1_772_912_406 > 15_552_000,
        };
        if scoring(i.command) {
            if frame != 0 && i.command == Command::BlastP {
                return Err(
                    "Frameshift alignments are only supported for translated searches.".into(),
                );
            }
            if range && frame == 0 {
                return Err("Query range culling is only supported in frameshift alignment mode (option -F).".into());
            }
            if !i.matrix_file.as_os_str().is_empty() {
                if i.gap_open == -1 || i.gap_extend == -1 {
                    return Err("Custom scoring matrices require setting the --gapopen and --gapextend options.".into());
                }
                if matches!(
                    i.output_format.first().map(Vec::as_slice),
                    Some(b"daa" | b"100")
                ) {
                    return Err(
                        "Custom scoring matrices are not supported for the DAA format.".into(),
                    );
                }
                if i.comp_based_stats > 1 {
                    return Err("This value for --comp-based-stats is not supported when using a custom scoring matrix.".into());
                }
            }
            let mut scoring_input = i.clone();
            scoring_input.frame_shift = frame;
            b.scoring(&scoring_input)?;
            b.masking()?
        }
        if temp(i.command) {
            c.tmpdir =
                crate::util::io::temp_file::source_temp_directory_path(&c.tmpdir, &c.output_file);
            b.temp_path(&c.tmpdir)?
        }
        for &s in &i.sensitivities {
            c.set_sens(s)?
        }
        b.enums()?;
        b.translator(i.query_gencode)?;
        if !matches!(i.strand.as_slice(), b"both" | b"minus" | b"plus") {
            return Err("Invalid value for parameter --strand".into());
        }
        if i.unfmt == b"fastq" || i.alfmt == b"fastq" {
            c.store_query_quality = true
        }
        b.auxiliary_formats(i)?;
        if i.command == Command::BlastX {
            if i.query_file.len() > 2 {
                return Err("A maximum of 2 query files is supported in blastx mode.".into());
            }
        } else if i.query_file.len() > 1 {
            return Err("--query/-q has more than one argument.".into());
        }
        if i.target_indexed && i.lowmem != 1 {
            return Err("--target-indexed requires -c1.".into());
        }
        if i.swipe_all {
            c.algo = Algo::DoubleIndexed
        }
        if range && i.taxon_k != 0 {
            return Err("--taxon-k is not supported for --range-culling mode.".into());
        }
        if i.multiprocessing && i.parallel_tmpdir.as_os_str().is_empty() {
            return Err("--multiprocessing requires setting --parallel-tmpdir".into());
        }
        if i.multiprocessing {
            b.mkdir_path(&i.parallel_tmpdir)?
        }
        b.completed(&c, i)?;
        Ok(c)
    }
}
fn scoring(c: Command) -> bool {
    matches!(
        c,
        Command::BlastP
            | Command::BlastX
            | Command::Benchmark
            | Command::Mask
            | Command::MakeDb
            | Command::Cluster
            | Command::DeepClust
            | Command::LinClust
            | Command::RegressionTest
            | Command::ClusterReassign
            | Command::ClusterRealign
            | Command::ListSeeds
            | Command::ModelSeqs
            | Command::MakeSeedTable
            | Command::Recluster
            | Command::ComputeMedoids
    )
}
fn temp(c: Command) -> bool {
    matches!(
        c,
        Command::BlastP
            | Command::BlastX
            | Command::BlastN
            | Command::Benchmark
            | Command::Mask
            | Command::Cluster
            | Command::RegressionTest
            | Command::ClusterReassign
            | Command::Recluster
            | Command::ComputeMedoids
            | Command::DeepClust
            | Command::LinClust
    )
}
fn append(path: &mut PathBuf, suffix: &str) {
    if !path.as_os_str().is_empty()
        && !path
            .as_os_str()
            .as_encoded_bytes()
            .ends_with(suffix.as_bytes())
    {
        path.as_mut_os_string().push(suffix);
    }
}

pub(crate) fn invocation_bytes(args: &[OsString]) -> Result<Vec<u8>, String> {
    invocation_bytes_io(args).map_err(|e| e.to_string())
}
pub(crate) fn invocation_bytes_io(args: &[OsString]) -> std::io::Result<Vec<u8>> {
    let mut out = Vec::new();
    for (index, arg) in args.iter().enumerate() {
        if index != 0 {
            out.push(b' ');
        }
        out.extend_from_slice(
            &crate::util::io::source_os_bytes(arg.as_os_str()).map_err(std::io::Error::other)?,
        );
    }
    Ok(out)
}

pub fn compressor(s: &str) -> Result<Compressor, String> {
    compressor_bytes(s.as_bytes()).map_err(|error| error.to_string())
}
pub fn compressor_bytes(s: &[u8]) -> std::io::Result<Compressor> {
    match s {
        b"" | b"0" => Ok(Compressor::None),
        b"1" => Ok(Compressor::Zlib),
        b"zstd" => Ok(Compressor::Zstd),
        value => {
            let mut message = b"Invalid compression algorithm: ".to_vec();
            message.extend_from_slice(value);
            Err(crate::util::io::RawDiagnosticError(message).into_io())
        }
    }
}
pub fn set_string_option<T: Copy + Default>(
    s: &str,
    n: &str,
    v: &[(&str, T)],
) -> Result<T, String> {
    set_bytes_option(s.as_bytes(), n, v).map_err(|e| e.to_string())
}
pub fn set_bytes_option<T: Copy + Default>(
    s: &[u8],
    n: &str,
    v: &[(&str, T)],
) -> std::io::Result<T> {
    if s.is_empty() {
        return Ok(T::default());
    }
    v.iter()
        .find(|x| x.0.as_bytes() == s)
        .map(|x| x.1)
        .ok_or_else(|| {
            let message = format!(
                "Invalid argument for option {n}. Allowed values are:{}",
                v.iter().map(|x| format!(" {}", x.0)).collect::<String>()
            );
            crate::util::io::RawDiagnosticError(message.into_bytes()).into_io()
        })
}

pub fn block_size(m: i64, d: i64, s: Sensitivity, l: bool, t: i32, no_limit: bool) -> (f64, i32) {
    let (x, c) = crate::config::block_size(m, d, s, l, t);
    if !no_limit {
        return (x, c);
    }
    let tr = crate::search::sensitivity::get_traits(s);
    let r = crate::basic::reduction::Reduction::default_reduction();
    let w = crate::basic::shape::Shape::from_code(
        crate::search::sensitivity::get_shape_codes(s)[0],
        &r,
    )
    .weight;
    let mut spl = if tr.sketch_size > 0 {
        tr.sketch_size as f64 / 200.0
    } else {
        1.0
    } / c as f64;
    if tr.minimizer_window > 0 {
        spl /= tr.minimizer_window as f64 / 2.0
    }
    let bits = crate::search::sensitivity::seedp_bits(w, t, c, &r);
    let join = 1.0 + t as f64 / (crate::basic::seed::seedp_count(bits) / c as u64) as f64;
    (((m as f64 / 1e9) / (18.0 * join * spl + 2.0)).max(0.001), c)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[derive(Default)]
    struct Backend(Vec<String>);
    impl ConfigBackend for Backend {
        fn cbs(&mut self, n: u32) -> Result<(), String> {
            self.0.push(format!("cbs:{n}"));
            Ok(())
        }
        fn scoring(&mut self, _: &ConfigInput) -> Result<(), String> {
            self.0.push("score".into());
            Ok(())
        }
        fn masking(&mut self) -> Result<(), String> {
            self.0.push("mask".into());
            Ok(())
        }
        fn translator(&mut self, n: u32) -> Result<(), String> {
            self.0.push(format!("translator:{n}"));
            Ok(())
        }
        fn temp(&mut self, p: &str) -> Result<(), String> {
            self.0.push(format!("temp:{p}"));
            Ok(())
        }
        fn mkdir(&mut self, p: &str) -> Result<(), String> {
            self.0.push(format!("mkdir:{p}"));
            Ok(())
        }
    }
    #[test]
    fn registration_stream_is_complete() {
        assert_eq!(registered_commands(true).len(), 45);
        assert_eq!(registered_commands(false).len(), 21);
        assert_eq!(registered_options(true).len(), 316);
        assert_eq!(
            registered_options(false)
                .iter()
                .filter(|&&x| x == "freq-sd")
                .count(),
            1
        );
        assert!(!registered_options(false).contains(&"zdrop"));
    }
    #[test]
    fn constructor_transforms_and_initializes() {
        let i = ConfigInput {
            argc: 4,
            command: Command::BlastX,
            database: "db".into(),
            daa_file: "hits".into(),
            long_reads: true,
            hardware_threads: 8,
            now: 1_772_912_406 + 181 * 24 * 3600,
            ..Default::default()
        };
        let mut b = Backend::default();
        let c = Config::new(&i, &mut b).unwrap();
        assert_eq!(c.output_file, Path::new("hits.daa"));
        assert_eq!(
            (c.range_culling, c.top, c.frame_shift),
            (true, Some(10.0), 15)
        );
        assert_eq!(c.threads, 8);
        assert!(c.old_warning);
        assert_eq!(b.0, ["cbs:1", "score", "mask", "temp:", "translator:1"]);
    }
    #[test]
    fn validations_and_helpers_match() {
        let i = ConfigInput {
            command: Command::BlastP,
            database: "db".into(),
            top_present: true,
            max_target_present: true,
            ..Default::default()
        };
        assert_eq!(
            Config::new(&i, &mut Backend::default()).unwrap_err(),
            "--top and -k/--max-target-seqs are mutually exclusive."
        );
        let i = ConfigInput {
            command: Command::BlastP,
            database: "db".into(),
            sensitivities: vec![Sensitivity::Fast, Sensitivity::Sensitive],
            ..Default::default()
        };
        assert_eq!(
            Config::new(&i, &mut Backend::default()).unwrap_err(),
            "Sensitivity switches are mutually exclusive."
        );
        assert_eq!(compressor("zstd").unwrap(), Compressor::Zstd);
        assert_eq!(
            set_string_option("two", "mode", &[("one", 1), ("two", 2)]).unwrap(),
            2
        );
    }
}
