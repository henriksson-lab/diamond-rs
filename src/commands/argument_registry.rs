//! Baseline DIAMOND command/option parse protocol.
//!
//! Static registration metadata follows the verified developer inventory in
//! `translation/command_options.json`. Nothing is loaded from that directory at
//! runtime: the published crate contains this complete Rust table.
//!
//! Validation checks registration, workflow membership and parameter counts.
//! It deliberately does not parse values, normalize arguments, execute commands,
//! or establish that registered behavior has a native implementation.
//!
//! Source commit: `1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`.
//! Inventory SHA-256: `1e1dbba30aa55c9c7efbc817f04cc3abf1baac8d4199dd3cffcc2284b20e35a5`.

use std::cell::RefCell;
use std::ffi::{OsStr, OsString};
use std::io;
use std::rc::Rc;

use crate::util::command_line_parser::{CommandLineParser, ReadOption};

#[path = "argument_descriptions.rs"]
mod descriptions;

// Registration feature bits retain all four audited source variants.
// Cargo `extra` selects Source EXTRA registration. Registration does not
// establish that each optional command has a native implementation.
// Bit 0: baseline; bit 1: WITH_DNA; bit 2: EXTRA; bit 3: EXTRA + WITH_DNA.
const BASELINE: u8 = 1;

pub(crate) const fn registration_profile() -> u8 {
    if cfg!(feature = "dna") {
        8
    } else if cfg!(feature = "extra") {
        4
    } else {
        BASELINE
    }
}

/// Source documentation uses the original group memberships even in EXTRA,
/// where option storage bypasses workflow permission checks.
pub(crate) fn documentation(command: &str) -> String {
    let profile = registration_profile();
    let mut parser = CommandLineParser::new();
    for (index, rule) in COMMANDS.iter().enumerate() {
        if rule.variants & profile != 0 {
            parser.add_command(
                rule.name,
                descriptions::command(rule.name),
                index as u32 + 1,
            );
        }
    }
    let groups: Vec<_> = GROUPS
        .iter()
        .map(|rule| {
            parser.add_group(
                rule.title,
                command_codes(rule.commands),
                rule.disabled_variants & profile != 0,
            )
        })
        .collect();
    let value = Rc::new(RefCell::new(Scalar));
    for rule in OPTIONS.iter().filter(|rule| rule.variants & profile != 0) {
        parser.add_option(
            groups[rule.group],
            rule.name,
            rule.short,
            &descriptions::option(rule.name),
            value.clone(),
            Scalar,
            rule.minimum,
        );
    }
    let canonical = command.strip_prefix("--").unwrap_or(command);
    if matches!(canonical, "help" | "-h") {
        parser.print_help()
    } else if let Some(index) = COMMANDS.iter().position(|rule| rule.name == canonical) {
        parser.print_documentation(index as u32 + 1)
    } else {
        String::new()
    }
}

#[derive(Clone, Copy)]
enum Arity {
    Scalar,
    Flag,
    Vector,
}
struct CommandRule {
    name: &'static str,
    variants: u8,
}
struct GroupRule {
    title: &'static str,
    commands: u64,
    disabled_variants: u8,
}
struct OptionRule {
    name: &'static str,
    short: char,
    group: usize,
    arity: Arity,
    minimum: i32,
    variants: u8,
}

// Codes are private registration indices, not a public C++ enum/ABI mapping.
// Source command order and shared group memberships are retained exactly.
const fn command(name: &'static str, variants: u8) -> CommandRule {
    CommandRule { name, variants }
}
const fn group(title: &'static str, commands: u64, disabled_variants: u8) -> GroupRule {
    GroupRule {
        title,
        commands,
        disabled_variants,
    }
}
const fn option(
    name: &'static str,
    short: char,
    group: usize,
    arity: Arity,
    minimum: i32,
    variants: u8,
) -> OptionRule {
    OptionRule {
        name,
        short,
        group,
        arity,
        minimum,
        variants,
    }
}

const COMMANDS: &[CommandRule] = &[
    command("makedb", 15),
    command("prepdb", 15),
    command("blastp", 15),
    command("blastx", 15),
    command("cluster", 15),
    command("linclust", 15),
    command("realign", 15),
    command("recluster", 15),
    command("reassign", 15),
    command("view", 15),
    command("merge-daa", 15),
    command("help", 15),
    command("version", 15),
    command("getseq", 15),
    command("dbinfo", 15),
    command("test", 15),
    command("makeidx", 15),
    command("greedy-vertex-cover", 15),
    command("roc", 15),
    command("benchmark", 15),
    command("deepclust", 15),
    command("random-seqs", 12),
    command("sort", 12),
    command("dbstat", 12),
    command("mask", 12),
    command("fastq2fasta", 12),
    command("info", 12),
    command("seed-stat", 12),
    command("smith-waterman", 12),
    command("simulate-seqs", 12),
    command("split", 12),
    command("upgma", 12),
    command("upgmamc", 12),
    command("reverse", 12),
    command("compute-medoids", 12),
    command("mutate", 12),
    command("roc-id", 12),
    command("find-shapes", 12),
    command("hashseqs", 12),
    command("listseeds", 12),
    command("blastn", 12),
    command("wc", 12),
    command("cut", 12),
    command("model-seqs", 12),
    command("make-seed-table", 12),
];

// Group command lists are identical in all four audited source variants.
// EXTRA disables workflow permission checks upstream; this module uses the
// baseline parser checks. The hidden group remains disabled in baseline/DNA.
const GROUPS: &[GroupRule] = &[
    group("General options", 0x1367fd, 0),          // general
    group("General/database options", 0x1361fd, 0), // general_db
    group("General/output", 0x1227fc, 0),           // general_out
    group("General/output2", 0x1203fc, 0),          // general_out2
    group("Makedb options", 0x401, 0),              // makedb_opt
    group("Makedb/taxon options", 0x1, 0),          // makedb_tax_opt
    group("Aligner/Clustering/Realign options", 0x1001fc, 0), // align_clust_realign
    group("Aligner/Clustering options", 0x1001bc, 0), // align_clust
    group("Aligner/view options", 0x20c, 0),        // aligner_view
    group("Aligner/sens options", 0x1000c, 0),      // aligner_sens
    group("Aligner options", 0xc, 0),               // aligner
    group("Output format options", 0x24c, 0),       // format
    group("Taxon output format options", 0xc, 0),   // taxon_format
    group("Clustering options", 0x1000b0, 0),       // cluster_opt
    group("Memory options", 0x1201f0, 0),           // memory_opt
    group("Clustering/reassign options", 0x1201b0, 0), // cluster_reassign_opt
    group("GVC options", 0x20000, 0),               // gvc_opt
    group("Cluster input options", 0x1c0, 0),       // realign_opt
    group("Advanced/general", 0x1000010813d, 0),    // advanced_gen
    group("Advanced options aln/cluster", 0x100001081bc, 0), // advanced_aln_cluster
    group("Advanced options", 0x1000000800c, 0),    // advanced
    group("View/Align options", 0x20c, 0),          // view_align_options
    group("View options", 0x200, 0),                // view_options
    group("Getseq options", 0x2000, 0),             // getseq_options
    group("", 0xc, 0),                              // deprecated_options
    group("", 0x0, 3),                              // hidden_options
];

// All 316 source names, in constructor order. Baseline registers 299 of them:
// 161 enabled plus 138 disabled entries that must still shadow short aliases.
// The other 17 are WITH_DNA registrations and are absent from baseline maps.
const OPTIONS: &[OptionRule] = &[
    option("threads", 'p', 0, Arity::Scalar, 1, 15),
    option("verbose", 'v', 0, Arity::Flag, 0, 15),
    option("log", '\0', 0, Arity::Flag, 0, 15),
    option("quiet", '\0', 0, Arity::Flag, 0, 15),
    option("tmpdir", 't', 0, Arity::Scalar, 1, 15),
    option("db", 'd', 1, Arity::Scalar, 1, 15),
    option("out", 'o', 2, Arity::Scalar, 1, 15),
    option("header", '\0', 3, Arity::Vector, 0, 15),
    option("in", '\0', 4, Arity::Vector, 1, 15),
    option("taxonmap", '\0', 5, Arity::Scalar, 1, 15),
    option("taxonnodes", '\0', 5, Arity::Scalar, 1, 15),
    option("taxonnames", '\0', 5, Arity::Scalar, 1, 15),
    option("comp-based-stats", '\0', 6, Arity::Scalar, 1, 15),
    option("masking", '\0', 6, Arity::Scalar, 1, 15),
    option("soft-masking", '\0', 6, Arity::Scalar, 1, 15),
    option("no-block-size-limit", '\0', 6, Arity::Flag, 0, 15),
    option("gapopen", '\0', 7, Arity::Scalar, 1, 15),
    option("gapextend", '\0', 7, Arity::Scalar, 1, 15),
    option("matrix", '\0', 7, Arity::Scalar, 1, 15),
    option("custom-matrix", '\0', 7, Arity::Scalar, 1, 15),
    option("evalue", 'e', 7, Arity::Scalar, 1, 15),
    option("motif-masking", '\0', 7, Arity::Scalar, 1, 15),
    option("approx-id", '\0', 7, Arity::Scalar, 1, 15),
    option("ext", '\0', 7, Arity::Scalar, 1, 15),
    option("max-target-seqs", 'k', 8, Arity::Scalar, 1, 15),
    option("top", '\0', 8, Arity::Scalar, 1, 15),
    option("faster", '\0', 9, Arity::Flag, 0, 15),
    option("fast", '\0', 9, Arity::Flag, 0, 15),
    option("mid-sensitive", '\0', 9, Arity::Flag, 0, 15),
    option("linclust-20", '\0', 9, Arity::Flag, 0, 15),
    option("shapes-6x10", '\0', 9, Arity::Flag, 0, 15),
    option("shapes-30x10", '\0', 9, Arity::Flag, 0, 15),
    option("sensitive", '\0', 9, Arity::Flag, 0, 15),
    option("more-sensitive", '\0', 9, Arity::Flag, 0, 15),
    option("very-sensitive", '\0', 9, Arity::Flag, 0, 15),
    option("ultra-sensitive", '\0', 9, Arity::Flag, 0, 15),
    option("shapes", 's', 9, Arity::Scalar, 1, 15),
    option("query", 'q', 10, Arity::Vector, 1, 15),
    option("strand", '\0', 10, Arity::Scalar, 1, 15),
    option("un", '\0', 10, Arity::Scalar, 1, 15),
    option("al", '\0', 10, Arity::Scalar, 1, 15),
    option("unfmt", '\0', 10, Arity::Scalar, 1, 15),
    option("alfmt", '\0', 10, Arity::Scalar, 1, 15),
    option("unal", '\0', 10, Arity::Scalar, 1, 15),
    option("max-hsps", '\0', 10, Arity::Scalar, 1, 15),
    option("range-culling", '\0', 10, Arity::Flag, 0, 15),
    option("compress", '\0', 10, Arity::Scalar, 1, 15),
    option("min-score", '\0', 10, Arity::Scalar, 1, 15),
    option("id", '\0', 10, Arity::Scalar, 1, 15),
    option("query-cover", '\0', 10, Arity::Scalar, 1, 15),
    option("subject-cover", '\0', 10, Arity::Scalar, 1, 15),
    option("swipe", '\0', 10, Arity::Flag, 0, 15),
    option("iterate", '\0', 10, Arity::Vector, 0, 15),
    option("global-ranking", 'g', 10, Arity::Scalar, 1, 15),
    option("block-size", 'b', 10, Arity::Scalar, 1, 15),
    option("index-chunks", 'c', 10, Arity::Scalar, 1, 15),
    option("frameshift", 'F', 10, Arity::Scalar, 1, 15),
    option("long-reads", '\0', 10, Arity::Flag, 0, 15),
    option("query-gencode", '\0', 10, Arity::Scalar, 1, 15),
    option("salltitles", '\0', 10, Arity::Flag, 0, 15),
    option("sallseqid", '\0', 10, Arity::Flag, 0, 15),
    option("no-self-hits", '\0', 10, Arity::Flag, 0, 15),
    option("taxonlist", '\0', 10, Arity::Scalar, 1, 15),
    option("taxon-exclude", '\0', 10, Arity::Scalar, 1, 15),
    option("seqidlist", '\0', 10, Arity::Scalar, 1, 15),
    option("skip-missing-seqids", '\0', 10, Arity::Flag, 0, 15),
    option("outfmt", 'f', 11, Arity::Vector, 1, 15),
    option("qnum-offset", '\0', 11, Arity::Scalar, 1, 15),
    option("snum-offset", '\0', 11, Arity::Scalar, 1, 15),
    option("include-lineage", '\0', 12, Arity::Flag, 0, 15),
    option("cluster-steps", '\0', 13, Arity::Vector, 1, 15),
    option("kmer-ranking", '\0', 13, Arity::Flag, 0, 15),
    option("round-coverage", '\0', 13, Arity::Vector, 1, 15),
    option("round-approx-id", '\0', 13, Arity::Vector, 1, 15),
    option("aln-out", '\0', 13, Arity::Scalar, 1, 15),
    option("memory-limit", 'M', 14, Arity::Scalar, 1, 15),
    option("member-cover", '\0', 15, Arity::Scalar, 1, 15),
    option("mutual-cover", '\0', 15, Arity::Scalar, 1, 15),
    option("connected-component-depth", '\0', 15, Arity::Vector, 1, 15),
    option("no-reassign", '\0', 15, Arity::Flag, 0, 15),
    option("centroid-out", '\0', 16, Arity::Scalar, 1, 15),
    option("edges", '\0', 16, Arity::Scalar, 1, 15),
    option("edge-format", '\0', 16, Arity::Scalar, 1, 15),
    option("symmetric", '\0', 16, Arity::Flag, 0, 15),
    option("clusters", '\0', 17, Arity::Scalar, 1, 15),
    option("file-buffer-size", '\0', 18, Arity::Scalar, 1, 15),
    option("no-unlink", '\0', 18, Arity::Flag, 0, 15),
    option("ignore-warnings", '\0', 18, Arity::Flag, 0, 15),
    option("no-parse-seqids", '\0', 18, Arity::Flag, 0, 15),
    option("parallel-tmpdir", '\0', 19, Arity::Scalar, 1, 15),
    option("bin", '\0', 19, Arity::Scalar, 1, 15),
    option("ext-chunk-size", '\0', 19, Arity::Scalar, 1, 15),
    option("no-ranking", '\0', 19, Arity::Flag, 0, 15),
    option("dbsize", '\0', 19, Arity::Scalar, 1, 15),
    option("no-auto-append", '\0', 19, Arity::Flag, 0, 15),
    option("tantan-minMaskProb", '\0', 19, Arity::Scalar, 1, 15),
    option("oid-output", '\0', 19, Arity::Flag, 0, 15),
    option("swipe-task-size", '\0', 19, Arity::Scalar, 1, 15),
    option("anchored-swipe", '\0', 19, Arity::Flag, 0, 15),
    option(
        "query-match-distance-threshold",
        '\0',
        19,
        Arity::Scalar,
        1,
        15,
    ),
    option("length-ratio-threshold", '\0', 19, Arity::Scalar, 1, 15),
    option("cbs-angle", '\0', 19, Arity::Scalar, 1, 15),
    option("linclust-banded-ext", '\0', 19, Arity::Flag, 0, 15),
    option("linclust-chunk-size", '\0', 19, Arity::Scalar, 1, 15),
    option("hit-membuf", '\0', 19, Arity::Flag, 0, 15),
    option("algo", '\0', 20, Arity::Scalar, 1, 15),
    option("min-orf", 'l', 20, Arity::Scalar, 1, 15),
    option("min-query-len", '\0', 20, Arity::Scalar, 1, 15),
    option("load-threads", '\0', 20, Arity::Scalar, 1, 15),
    option("minichunk", '\0', 20, Arity::Scalar, 1, 15),
    option("seed-cut", '\0', 20, Arity::Scalar, 1, 15),
    option("freq-masking", '\0', 20, Arity::Flag, 0, 15),
    option("freq-sd", '\0', 20, Arity::Scalar, 1, 15),
    option("sketch-size", '\0', 20, Arity::Scalar, 1, 15),
    option("tile-size", '\0', 20, Arity::Scalar, 1, 15),
    option("id2", '\0', 20, Arity::Scalar, 1, 15),
    option("linsearch", '\0', 20, Arity::Flag, 0, 15),
    option("lin-stage1", '\0', 20, Arity::Flag, 0, 15),
    option("lin-combo", '\0', 20, Arity::Flag, 0, 15),
    option("xdrop", 'x', 20, Arity::Scalar, 1, 15),
    option("ungapped-evalue", '\0', 20, Arity::Scalar, 1, 15),
    option("ungapped-evalue-short", '\0', 20, Arity::Scalar, 1, 15),
    option(
        "short-query-ungapped-bitscore",
        '\0',
        20,
        Arity::Scalar,
        1,
        15,
    ),
    option("gapped-filter-evalue", '\0', 20, Arity::Scalar, 1, 15),
    option("band", '\0', 20, Arity::Scalar, 1, 15),
    option("shape-mask", '\0', 20, Arity::Vector, 1, 15),
    option("multiprocessing", '\0', 20, Arity::Flag, 0, 15),
    option("mp-init", '\0', 20, Arity::Flag, 0, 15),
    option("mp-recover", '\0', 20, Arity::Flag, 0, 15),
    option("mp-query-chunk", '\0', 20, Arity::Scalar, 1, 15),
    option("culling-overlap", '\0', 20, Arity::Scalar, 1, 15),
    option("taxon-k", '\0', 20, Arity::Scalar, 1, 15),
    option("range-cover", '\0', 20, Arity::Scalar, 1, 15),
    option("xml-blord-format", '\0', 20, Arity::Flag, 0, 15),
    option("sam-query-len", '\0', 20, Arity::Flag, 0, 15),
    option("stop-match-score", '\0', 20, Arity::Scalar, 1, 15),
    option("target-indexed", '\0', 20, Arity::Flag, 0, 15),
    option("unaligned-targets", '\0', 20, Arity::Scalar, 1, 15),
    option("cut-bar", '\0', 20, Arity::Flag, 0, 15),
    option("check-multi-target", '\0', 20, Arity::Flag, 0, 15),
    option("roc-file", '\0', 20, Arity::Scalar, 1, 15),
    option("family-map", '\0', 20, Arity::Scalar, 1, 15),
    option("family-map-query", '\0', 20, Arity::Scalar, 1, 15),
    option("query-parallel-limit", '\0', 20, Arity::Scalar, 1, 15),
    option("log-evalue-scale", '\0', 20, Arity::Scalar, 1, 15),
    option("bootstrap", '\0', 20, Arity::Flag, 0, 15),
    option("heartbeat", '\0', 20, Arity::Flag, 0, 15),
    option("mp-self", '\0', 20, Arity::Flag, 0, 15),
    option("query-or-subject-cover", '\0', 20, Arity::Scalar, 1, 15),
    option("daa", 'a', 21, Arity::Scalar, 1, 15),
    option("forwardonly", '\0', 22, Arity::Flag, 0, 15),
    option("seq", '\0', 23, Arity::Vector, 1, 15),
    option("window", 'w', 24, Arity::Scalar, 1, 15),
    option("ungapped-score", '\0', 24, Arity::Scalar, 1, 15),
    option("hit-band", '\0', 24, Arity::Scalar, 1, 15),
    option("hit-score", '\0', 24, Arity::Scalar, 1, 15),
    option("gapped-xdrop", 'X', 24, Arity::Scalar, 1, 15),
    option("rank-ratio2", '\0', 24, Arity::Scalar, 1, 15),
    option("rank-ratio", '\0', 24, Arity::Scalar, 1, 15),
    option("lambda", '\0', 24, Arity::Scalar, 1, 15),
    option("K", '\0', 24, Arity::Scalar, 1, 15),
    option("match1", '\0', 25, Arity::Scalar, 1, 15),
    option("match2", '\0', 25, Arity::Scalar, 1, 15),
    option("seed-freq", '\0', 25, Arity::Scalar, 1, 15),
    option("space-penalty", '\0', 25, Arity::Scalar, 1, 15),
    option("reverse", '\0', 25, Arity::Flag, 0, 15),
    option("neighborhood-score", '\0', 25, Arity::Scalar, 1, 15),
    option("seed-weight", 'w', 25, Arity::Scalar, 1, 15),
    option("superblock", '\0', 25, Arity::Scalar, 1, 15),
    option("load-balancing", '\0', 25, Arity::Scalar, 1, 15),
    option("log-query", '\0', 25, Arity::Flag, 0, 15),
    option("log-subject", '\0', 25, Arity::Flag, 0, 15),
    option("palign", '\0', 25, Arity::Scalar, 1, 15),
    option("score-ratio", '\0', 25, Arity::Scalar, 1, 15),
    option("fetch-size", '\0', 25, Arity::Scalar, 1, 15),
    option("target-fetch-size", '\0', 25, Arity::Scalar, 1, 15),
    option("rank-factor", '\0', 25, Arity::Scalar, 1, 15),
    option("transcript-len-estimate", '\0', 25, Arity::Scalar, 1, 15),
    option("family-counts", '\0', 25, Arity::Scalar, 1, 15),
    option("radix-cluster-buffered", '\0', 25, Arity::Flag, 0, 15),
    option("join-split-size", '\0', 25, Arity::Scalar, 1, 15),
    option("join-split-key-len", '\0', 25, Arity::Scalar, 1, 15),
    option("radix-bits", '\0', 25, Arity::Scalar, 1, 15),
    option("join-ht-factor", '\0', 25, Arity::Scalar, 1, 15),
    option("sort-join", '\0', 25, Arity::Flag, 0, 15),
    option("simple-freq", '\0', 25, Arity::Flag, 0, 15),
    option("freq-treshold", '\0', 25, Arity::Scalar, 1, 15),
    option("use-dataset-field", '\0', 25, Arity::Flag, 0, 15),
    option("store-query-quality", '\0', 25, Arity::Flag, 0, 15),
    option("swipe-chunk-size", '\0', 25, Arity::Scalar, 1, 15),
    option("hard-masked", '\0', 25, Arity::Flag, 0, 15),
    option("cbs-window", '\0', 25, Arity::Scalar, 1, 15),
    option("no-dict", '\0', 25, Arity::Flag, 0, 15),
    option("upgma-edge-limit", '\0', 25, Arity::Scalar, 1, 15),
    option("tree", '\0', 25, Arity::Scalar, 1, 15),
    option("upgma-dist", '\0', 25, Arity::Scalar, 1, 15),
    option("upgma-input", '\0', 25, Arity::Scalar, 1, 15),
    option("log-extend", '\0', 25, Arity::Flag, 0, 15),
    option("chaining-maxgap", '\0', 25, Arity::Scalar, 1, 15),
    option("tantan-maxRepeatOffset", '\0', 25, Arity::Scalar, 1, 15),
    option("tantan-ungapped", '\0', 25, Arity::Flag, 0, 15),
    option("chaining-range-cover", '\0', 25, Arity::Scalar, 1, 15),
    option("no-swipe-realign", '\0', 25, Arity::Flag, 0, 15),
    option("chaining-maxnodes", '\0', 25, Arity::Scalar, 1, 15),
    option("cutoff-score-8bit", '\0', 25, Arity::Scalar, 1, 15),
    option("min-band-overlap", '\0', 25, Arity::Scalar, 1, 15),
    option("min-realign-overhang", '\0', 25, Arity::Scalar, 1, 15),
    option("ungapped-window", '\0', 25, Arity::Scalar, 1, 15),
    option("gapped-filter-diag-score", '\0', 25, Arity::Scalar, 1, 15),
    option("gapped-filter-window", '\0', 25, Arity::Scalar, 1, 15),
    option("output-hits", '\0', 25, Arity::Flag, 0, 15),
    option("no-logfile", '\0', 25, Arity::Flag, 0, 15),
    option("band-bin", '\0', 25, Arity::Scalar, 1, 15),
    option("col-bin", '\0', 25, Arity::Scalar, 1, 15),
    option("self", '\0', 25, Arity::Flag, 0, 15),
    option("trace-pt-fetch-size", '\0', 25, Arity::Scalar, 1, 15),
    option("short-query-max-len", '\0', 25, Arity::Scalar, 1, 15),
    option("gapped-filter-evalue1", '\0', 25, Arity::Scalar, 1, 15),
    option("ext-yield", '\0', 25, Arity::Scalar, 1, 15),
    option("full-sw-len", '\0', 25, Arity::Scalar, 1, 15),
    option("relaxed-evalue-factor", '\0', 25, Arity::Scalar, 1, 15),
    option("type", '\0', 25, Arity::Scalar, 1, 15),
    option("raw", '\0', 25, Arity::Flag, 0, 15),
    option("chaining-len-cap", '\0', 25, Arity::Scalar, 1, 15),
    option("chaining-min-nodes", '\0', 25, Arity::Scalar, 1, 15),
    option("fast-tsv", '\0', 25, Arity::Flag, 0, 15),
    option("target-parallel-verbosity", '\0', 25, Arity::Scalar, 1, 15),
    option("query-memory", '\0', 25, Arity::Flag, 0, 15),
    option("memory-intervals", '\0', 25, Arity::Scalar, 1, 15),
    option("seed-hit-density", '\0', 25, Arity::Scalar, 1, 15),
    option("chunk-size-multiplier", '\0', 25, Arity::Scalar, 1, 15),
    option("score-drop-factor", '\0', 25, Arity::Scalar, 1, 15),
    option("left-most-interval", '\0', 25, Arity::Scalar, 1, 15),
    option("ranking-cutoff-bitscore", '\0', 25, Arity::Scalar, 1, 15),
    option("no-forward-fp", '\0', 25, Arity::Flag, 0, 15),
    option("no-ref-masking", '\0', 25, Arity::Flag, 0, 15),
    option("target-bias", '\0', 25, Arity::Flag, 0, 15),
    option("output-fp", '\0', 25, Arity::Flag, 0, 15),
    option("family-cap", '\0', 25, Arity::Scalar, 1, 15),
    option("cbs-matrix-scale", '\0', 25, Arity::Scalar, 1, 15),
    option("query-count", '\0', 25, Arity::Scalar, 1, 15),
    option("cbs-err-tolerance", '\0', 25, Arity::Scalar, 1, 15),
    option("cbs-it-limit", '\0', 25, Arity::Scalar, 1, 15),
    option("hash_join_swap", '\0', 25, Arity::Flag, 0, 15),
    option("deque_bucket_size", '\0', 25, Arity::Scalar, 1, 15),
    option("max-swipe-dp", '\0', 25, Arity::Scalar, 1, 15),
    option("no-reextend", '\0', 25, Arity::Flag, 0, 15),
    option("no-reorder", '\0', 25, Arity::Flag, 0, 15),
    option("file1", '\0', 25, Arity::Scalar, 1, 15),
    option("file2", '\0', 25, Arity::Scalar, 1, 15),
    option("key2", '\0', 25, Arity::Scalar, 1, 15),
    option("motif-mask-file", '\0', 25, Arity::Scalar, 1, 15),
    option("max-motif-len", '\0', 25, Arity::Scalar, 1, 15),
    option("chaining-stacked-hsp-ratio", '\0', 25, Arity::Scalar, 1, 15),
    option("minimizer-window", '\0', 25, Arity::Scalar, 1, 15),
    option("min_task_trace_pts", '\0', 25, Arity::Scalar, 1, 15),
    option("oid-list", '\0', 25, Arity::Scalar, 1, 15),
    option("bootstrap-block", '\0', 25, Arity::Scalar, 1, 15),
    option("centroid-factor", '\0', 25, Arity::Scalar, 1, 15),
    option("timeout", '\0', 25, Arity::Scalar, 1, 15),
    option("resume", '\0', 25, Arity::Scalar, 1, 15),
    option("target_hard_cap", '\0', 25, Arity::Scalar, 1, 15),
    option("mapany", '\0', 25, Arity::Flag, 0, 15),
    option("neighbors", '\0', 25, Arity::Scalar, 1, 15),
    option("reassign-overlap", '\0', 25, Arity::Scalar, 1, 15),
    option("reassign-ratio", '\0', 25, Arity::Scalar, 1, 15),
    option("reassign-max", '\0', 25, Arity::Scalar, 1, 15),
    option("add-self-aln", '\0', 25, Arity::Flag, 0, 15),
    option("weighted-gvc", '\0', 25, Arity::Flag, 0, 15),
    option("hamming-ext", '\0', 25, Arity::Flag, 0, 15),
    option("diag-filter-id", '\0', 25, Arity::Scalar, 1, 15),
    option("diag-filter-cov", '\0', 25, Arity::Scalar, 1, 15),
    option("strict-gvc", '\0', 25, Arity::Flag, 0, 15),
    option("dbtype", '\0', 25, Arity::Scalar, 1, 15),
    option("cluster-similarity", '\0', 25, Arity::Scalar, 1, 15),
    option("cluster-threshold", '\0', 25, Arity::Scalar, 1, 15),
    option("cluster-graph-file", '\0', 25, Arity::Scalar, 1, 15),
    option("cluster-restart", '\0', 25, Arity::Flag, 0, 15),
    option("mcl-expansion", '\0', 25, Arity::Scalar, 1, 15),
    option("mcl-inflation", '\0', 25, Arity::Scalar, 1, 15),
    option("mcl-chunk-size", '\0', 25, Arity::Scalar, 1, 15),
    option("mcl-max-iterations", '\0', 25, Arity::Scalar, 1, 15),
    option("mcl-sparsity-switch", '\0', 25, Arity::Scalar, 1, 15),
    option("mcl-nonsymmetric", '\0', 25, Arity::Flag, 0, 15),
    option("mcl-stats", '\0', 25, Arity::Flag, 0, 15),
    option("cluster-algo", '\0', 25, Arity::Scalar, 1, 15),
    option("approx-backtrace", '\0', 25, Arity::Flag, 0, 15),
    option("narrow-band-cov", '\0', 25, Arity::Scalar, 1, 15),
    option("narrow-band-factor", '\0', 25, Arity::Scalar, 1, 15),
    option("anchor-window", '\0', 25, Arity::Scalar, 1, 15),
    option("anchor-score", '\0', 25, Arity::Scalar, 1, 15),
    option("classic-band", '\0', 25, Arity::Flag, 0, 15),
    option("no_8bit_extension", '\0', 25, Arity::Flag, 0, 15),
    option("no_chaining_merge_hsps", '\0', 25, Arity::Flag, 0, 15),
    option("graph-algo", '\0', 25, Arity::Scalar, 1, 15),
    option("tsv-read-size", '\0', 25, Arity::Scalar, 1, 15),
    option("min-len-ratio", '\0', 25, Arity::Scalar, 1, 15),
    option("max-indirection", '\0', 25, Arity::Scalar, 1, 15),
    option("promiscuous-seed-ratio", '\0', 25, Arity::Scalar, 1, 15),
    option("zdrop", 'z', 20, Arity::Scalar, 1, 10),
    option("repetition-cutoff", '\0', 20, Arity::Scalar, 1, 10),
    option("extension", '\0', 20, Arity::Scalar, 1, 10),
    option("chaining-out", '\0', 20, Arity::Flag, 0, 10),
    option("align-long-reads", '\0', 20, Arity::Flag, 0, 10),
    option("best-hsp-only", '\0', 20, Arity::Flag, 0, 10),
    option("chain-pen-gap-scale", '\0', 20, Arity::Scalar, 1, 10),
    option("chain-pen-skip-scale", '\0', 20, Arity::Scalar, 1, 10),
    option("penalty", '\0', 20, Arity::Scalar, 1, 10),
    option("reward", '\0', 20, Arity::Scalar, 1, 10),
    option("chain-align-cutoff", '\0', 20, Arity::Scalar, 1, 10),
    option("min-chain-score", '\0', 20, Arity::Scalar, 1, 10),
    option("max-overlap-extension", '\0', 20, Arity::Scalar, 1, 10),
    option("zdrop-extension", '\0', 20, Arity::Scalar, 1, 10),
    option("zdrop-global", '\0', 20, Arity::Scalar, 1, 10),
    option("band-extension", '\0', 20, Arity::Scalar, 1, 10),
    option("band-global", '\0', 20, Arity::Scalar, 1, 10),
];

// Shared no-op stores retain only type-specific arity rules. Values remain
// untouched in the caller's argv and are read by its actual command parser.
#[derive(Clone, Default)]
struct Scalar;
impl ReadOption for Scalar {
    fn read_option(_: &mut Self, _: &[String]) -> Result<(), String> {
        Ok(())
    }
}
#[derive(Clone, Default)]
struct Flag;
impl ReadOption for Flag {
    fn read_option(_: &mut Self, _: &[String]) -> Result<(), String> {
        Ok(())
    }
    fn check_pcount(values: &[String], _: i32) -> bool {
        values.is_empty()
    }
}
#[derive(Clone, Default)]
struct Vector;
impl ReadOption for Vector {
    fn read_option(_: &mut Self, _: &[String]) -> Result<(), String> {
        Ok(())
    }
    fn check_pcount(values: &[String], minimum: i32) -> bool {
        values.len() >= minimum as usize
    }
}

fn command_codes(mask: u64) -> Vec<u32> {
    (0..COMMANDS.len())
        .filter(|&index| mask & (1u64 << index) != 0)
        .map(|index| index as u32 + 1)
        .collect()
}

/// Validate the baseline source argument protocol, preserving original argv.
///
/// Native extensions are intentionally limited to search memory limits,
/// search --help/-h, and -V/-h command conveniences. Exact --legacy bypasses
/// baseline validation for the explicit optional conformance adapter.
/// The caller retains its existing early no-arguments usage convenience; this
/// helper itself returns the translated source syntax error when no command
/// is present. Accepted options still require an actual native command backend.
pub fn validate(args: &[String]) -> Result<(), String> {
    let args: Vec<OsString> = args.iter().map(OsString::from).collect();
    validate_os(&args).map_err(|error| error.to_string())
}

pub(crate) type OsValues = std::collections::BTreeMap<&'static str, Vec<OsString>>;

/// A registry cell captures the actual option values after the shared parser
/// resolves aliases and arity. Its OS representation is authoritative.
#[derive(Clone)]
struct RegistryCell<const ARITY: u8> {
    key: &'static str,
    values: Rc<RefCell<OsValues>>,
}
impl<const ARITY: u8> RegistryCell<ARITY> {
    fn new(key: &'static str, values: Rc<RefCell<OsValues>>) -> Self {
        Self { key, values }
    }
}
impl<const ARITY: u8> ReadOption for RegistryCell<ARITY> {
    fn read_option(dst: &mut Self, values: &[String]) -> Result<(), String> {
        let values: Vec<OsString> = values.iter().map(OsString::from).collect();
        Self::read_option_os(dst, &values).map_err(|error| error.to_string())
    }
    fn read_option_os(dst: &mut Self, values: &[OsString]) -> io::Result<()> {
        dst.values.borrow_mut().insert(dst.key, values.to_vec());
        Ok(())
    }
    fn check_pcount_len(count: usize, minimum: i32) -> bool {
        match ARITY {
            0 => count == 1,
            1 => count == 0,
            _ => count >= minimum as usize,
        }
    }
}

pub(crate) fn validate_os(args: &[OsString]) -> io::Result<()> {
    values_os(args).map(|_| ())
}

/// Source option ownership, last occurrence and short-ID shadowing, using the
/// same parser for raw paths and text values without a preliminary UTF-8 gate.
pub(crate) fn values_os(args: &[OsString]) -> io::Result<OsValues> {
    if args
        .iter()
        .skip(2)
        .any(|argument| argument == OsStr::new("--legacy"))
    {
        return Ok(OsValues::new());
    }
    let mut parser = CommandLineParser::new();
    let profile = registration_profile();
    let extra = cfg!(feature = "extra");
    let all_commands = (1u64 << COMMANDS.len()) - 1;
    for (index, command) in COMMANDS.iter().enumerate() {
        if command.variants & profile != 0 {
            parser.add_command(command.name, "", index as u32 + 1);
        }
    }
    let groups: Vec<_> = GROUPS
        .iter()
        .map(|group| {
            parser.add_group(
                group.title,
                command_codes(if extra { all_commands } else { group.commands }),
                group.disabled_variants & profile != 0,
            )
        })
        .collect();
    let search_mask = COMMANDS
        .iter()
        .enumerate()
        .filter(|(_, command)| matches!(command.name, "blastp" | "blastx"))
        .fold(0u64, |mask, (index, _)| mask | 1u64 << index);
    // Extend this one option's permissions, not its whole source group.
    let memory_rule = OPTIONS
        .iter()
        .find(|rule| rule.name == "memory-limit")
        .unwrap();
    let memory_group = parser.add_group(
        GROUPS[memory_rule.group].title,
        command_codes(if extra {
            all_commands
        } else {
            GROUPS[memory_rule.group].commands | search_mask
        }),
        false,
    );
    let help_group = parser.add_group("Native search help", command_codes(search_mask), false);
    let captured = Rc::new(RefCell::new(OsValues::new()));
    for rule in OPTIONS.iter().filter(|rule| rule.variants & profile != 0) {
        let group = if rule.name == "memory-limit" {
            memory_group
        } else {
            groups[rule.group]
        };
        match rule.arity {
            Arity::Scalar => {
                parser.add_option(
                    group,
                    rule.name,
                    rule.short,
                    "",
                    Rc::new(RefCell::new(RegistryCell::<0>::new(
                        rule.name,
                        captured.clone(),
                    ))),
                    RegistryCell::<0>::new(rule.name, captured.clone()),
                    rule.minimum,
                );
            }
            Arity::Flag => {
                parser.add_option(
                    group,
                    rule.name,
                    rule.short,
                    "",
                    Rc::new(RefCell::new(RegistryCell::<1>::new(
                        rule.name,
                        captured.clone(),
                    ))),
                    RegistryCell::<1>::new(rule.name, captured.clone()),
                    rule.minimum,
                );
            }
            Arity::Vector => {
                parser.add_option(
                    group,
                    rule.name,
                    rule.short,
                    "",
                    Rc::new(RefCell::new(RegistryCell::<2>::new(
                        rule.name,
                        captured.clone(),
                    ))),
                    RegistryCell::<2>::new(rule.name, captured.clone()),
                    rule.minimum,
                );
            }
        }
    }
    // This native convenience must not change Source's unknown-option
    // diagnostics in unrelated commands merely by existing in the registry.
    if args.get(1).is_some_and(|command| {
        matches!(
            command
                .to_str()
                .and_then(|name| name.strip_prefix("--").or(Some(name)))
                .unwrap_or(""),
            "blastp" | "blastx"
        )
    }) {
        let cell = RegistryCell::<1>::new("help", captured.clone());
        parser.add_option(
            help_group,
            "help",
            'h',
            "",
            Rc::new(RefCell::new(cell.clone())),
            cell,
            0,
        );
    }
    let mut argv: Vec<&OsStr> = args.iter().map(OsString::as_os_str).collect();
    match argv.get(1).and_then(|arg| arg.to_str()) {
        Some("-V") => argv[1] = OsStr::new("version"),
        Some("-h") => argv[1] = OsStr::new("help"),
        _ => (),
    }
    parser.store_os(&argv, &mut 0)?;
    let values = std::mem::take(&mut *captured.borrow_mut());
    Ok(values)
}

/// Checked String convenience over the actual OS-value registry core.
pub(crate) fn values(
    args: &[String],
) -> Result<std::collections::BTreeMap<&'static str, Vec<String>>, String> {
    let args: Vec<OsString> = args.iter().map(OsString::from).collect();
    values_os(&args)
        .map_err(|error| error.to_string())?
        .into_iter()
        .map(|(key, values)| {
            values
                .into_iter()
                .map(|value| {
                    value
                        .into_string()
                        .map_err(|_| "Option value is not UTF-8; use values_os".to_owned())
                })
                .collect::<Result<Vec<_>, _>>()
                .map(|values| (key, values))
        })
        .collect()
}

#[cfg(feature = "extra")]
pub(crate) fn extra_values(
    args: &[String],
) -> Result<std::collections::BTreeMap<&'static str, Vec<String>>, String> {
    values(args)
}

#[cfg(feature = "extra")]
pub(crate) fn extra_values_os(args: &[OsString]) -> io::Result<OsValues> {
    values_os(args)
}

#[cfg(feature = "extra")]
pub(crate) fn is_extra_command(name: &str) -> bool {
    COMMANDS
        .iter()
        .any(|command| command.name == name && command.variants & BASELINE == 0)
}
