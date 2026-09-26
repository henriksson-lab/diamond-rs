//! Testable translation of `diamond/src/basic/config.cpp`.

use crate::config::Sensitivity;
use crate::util::io::Compressor;

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
    pub argv: Vec<String>,
    pub check_io: bool,
    pub command: Command,
    pub debug_log: bool,
    pub quiet: bool,
    pub verbose: bool,
    pub database: String,
    pub output_file: String,
    pub daa_file: String,
    pub output_format: Vec<String>,
    pub compression: String,
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
    pub matrix_file: String,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub sensitivities: Vec<Sensitivity>,
    pub algo: Algo,
    pub query_gencode: u32,
    pub strand: String,
    pub unfmt: String,
    pub alfmt: String,
    pub query_file: Vec<String>,
    pub store_query_quality: bool,
    pub target_indexed: bool,
    pub lowmem: u32,
    pub swipe_all: bool,
    pub parallel_tmpdir: String,
    pub hit_membuf: bool,
    pub tmpdir: String,
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
            database: String::new(),
            output_file: String::new(),
            daa_file: String::new(),
            output_format: vec![],
            compression: String::new(),
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
            matrix_file: String::new(),
            gap_open: -1,
            gap_extend: -1,
            sensitivities: vec![],
            algo: Algo::Auto,
            query_gencode: 1,
            strand: "both".into(),
            unfmt: "fasta".into(),
            alfmt: "fasta".into(),
            query_file: vec![],
            store_query_quality: false,
            target_indexed: false,
            lowmem: 4,
            swipe_all: false,
            parallel_tmpdir: String::new(),
            hit_membuf: false,
            tmpdir: String::new(),
            threads: None,
            hardware_threads: 1,
            now: 0,
        }
    }
}

pub trait ConfigBackend {
    fn cbs(&mut self, _: u32) -> Result<(), String>;
    fn scoring(&mut self, _: &ConfigInput) -> Result<(), String>;
    fn masking(&mut self) -> Result<(), String>;
    fn translator(&mut self, _: u32) -> Result<(), String>;
    fn temp(&mut self, _: &str) -> Result<(), String>;
    fn mkdir(&mut self, _: &str) -> Result<(), String>;
}

#[derive(Debug, Clone, PartialEq)]
pub struct Config {
    pub command: Command,
    pub verbosity: u8,
    pub database: String,
    pub output_file: String,
    pub daa_file: String,
    pub range_culling: bool,
    pub top: Option<f64>,
    pub frame_shift: i32,
    pub sensitivity: Sensitivity,
    pub algo: Algo,
    pub store_query_quality: bool,
    pub tmpdir: String,
    pub threads: i32,
    pub trace_pt_membuf: bool,
    pub invocation: String,
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
        i.query_file.first().map(String::as_str).unwrap_or("")
    }
    pub fn compressor(&self, i: &ConfigInput) -> Result<Compressor, String> {
        compressor(&i.compression)
    }
    pub fn new<B: ConfigBackend>(i: &ConfigInput, b: &mut B) -> Result<Self, String> {
        let verbosity = if i.debug_log {
            3
        } else if i.quiet {
            0
        } else if i.verbose {
            2
        } else if (matches!(
            i.command,
            Command::View | Command::BlastP | Command::BlastX | Command::BlastN
        ) && i.output_file.is_empty()
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
                if db.is_empty() {
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
                if db.is_empty() {
                    return Err("Missing parameter: database file (--db/-d)".into());
                }
                if !daa.is_empty() {
                    if !out.is_empty() {
                        return Err("Options --daa and --out cannot be used together.".into());
                    }
                    if i.output_format.first().is_some_and(|x| x != "daa") {
                        return Err("Invalid parameter: --daa/-a. Output file is specified with the --out/-o parameter.".into());
                    }
                    out = daa.clone()
                }
                if !daa.is_empty()
                    || matches!(
                        i.output_format.first().map(String::as_str),
                        Some("daa" | "100")
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
            if i.command == Command::DbInfo && db.is_empty() {
                return Err("Missing parameter: database file (--db/-d)".into());
            }
        }
        if !i.no_auto_append {
            if i.command == Command::MakeDb {
                append(&mut db, ".dmnd")
            }
            if i.command == Command::View {
                append(&mut daa, ".daa")
            }
            if i.compression == "1" {
                append(&mut out, ".gz")
            }
            if i.compression == "zstd" {
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
            threads: i.threads.unwrap_or(i.hardware_threads),
            trace_pt_membuf: i.hit_membuf,
            invocation: i.argv.join(" "),
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
            if !i.matrix_file.is_empty() {
                if i.gap_open == -1 || i.gap_extend == -1 {
                    return Err("Custom scoring matrices require setting the --gapopen and --gapextend options.".into());
                }
                if matches!(
                    i.output_format.first().map(String::as_str),
                    Some("daa" | "100")
                ) {
                    return Err(
                        "Custom scoring matrices are not supported for the DAA format.".into(),
                    );
                }
                if i.comp_based_stats > 1 {
                    return Err("This value for --comp-based-stats is not supported when using a custom scoring matrix.".into());
                }
            }
            b.scoring(i)?;
            b.masking()?
        }
        if temp(i.command) {
            if c.tmpdir.is_empty() {
                c.tmpdir = c
                    .output_file
                    .rsplit_once('/')
                    .map(|x| x.0.into())
                    .unwrap_or_default()
            }
            b.temp(&c.tmpdir)?
        }
        for &s in &i.sensitivities {
            c.set_sens(s)?
        }
        b.translator(i.query_gencode)?;
        if !matches!(i.strand.as_str(), "both" | "minus" | "plus") {
            return Err("Invalid value for parameter --strand".into());
        }
        if i.unfmt == "fastq" || i.alfmt == "fastq" {
            c.store_query_quality = true
        }
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
        if i.multiprocessing && i.parallel_tmpdir.is_empty() {
            return Err("--multiprocessing requires setting --parallel-tmpdir".into());
        }
        if i.multiprocessing {
            b.mkdir(&i.parallel_tmpdir)?
        }
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
            | Command::Recluster
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
            | Command::DeepClust
            | Command::LinClust
    )
}
fn append(s: &mut String, e: &str) {
    if !s.is_empty() && !s.ends_with(e) {
        s.push_str(e)
    }
}
pub fn compressor(s: &str) -> Result<Compressor, String> {
    match s {
        "" | "0" => Ok(Compressor::None),
        "1" => Ok(Compressor::Zlib),
        "zstd" => Ok(Compressor::Zstd),
        x => Err(format!("Invalid compression algorithm: {x}")),
    }
}
pub fn set_string_option<T: Copy + Default>(
    s: &str,
    n: &str,
    v: &[(&str, T)],
) -> Result<T, String> {
    if s.is_empty() {
        return Ok(T::default());
    }
    v.iter().find(|x| x.0 == s).map(|x| x.1).ok_or_else(|| {
        format!(
            "Invalid argument for option {n}. Allowed values are:{}",
            v.iter().map(|x| format!(" {}", x.0)).collect::<String>()
        )
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
        assert_eq!(c.output_file, "hits.daa");
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
