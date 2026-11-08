//! Search controls exercise native seed eligibility and DP filtering against
//! frozen outputs from the unchanged upstream executable.
use diamond::commands::blastp::{self, BlastpConfig, TranslatedQueryLayout};
use diamond::data::{block::Block, fasta};
use diamond::masking::{MaskingAlgo, MaskingMode};
use std::path::{Path, PathBuf};
use std::process::Command;

const FIELDS: &[&str] = &[
    "qseqid",
    "sseqid",
    "approx_pident",
    "bitscore",
    "qstart",
    "qend",
    "sstart",
    "send",
    "cigar",
    "qseq",
    "sseq",
];
fn directory(name: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "diamond-search-controls-{name}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&path).unwrap();
    path
}
fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/search_controls")
}
fn make_database(dir: &Path) -> PathBuf {
    let path = dir.join("reference.dmnd");
    let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .arg("makedb")
        .arg("--in")
        .arg(fixtures().join("reference.faa"))
        .arg("-d")
        .arg(&path)
        .args(["-p", "1"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    path
}
fn compare(actual: &[u8], expected: &[u8], label: &str) {
    let at = actual
        .iter()
        .zip(expected)
        .position(|(a, b)| a != b)
        .unwrap_or(actual.len().min(expected.len()));
    assert!(
        actual == expected,
        "{label}: first difference {at}; lengths {} and {}",
        actual.len(),
        expected.len()
    );
}
fn fields(kind: &str) -> &[&str] {
    match kind {
        "detailed" => FIELDS,
        "score" => &["qseqid", "sseqid", "bitscore"],
        "none" => &[],
        _ => panic!("unknown fields"),
    }
}
fn flags(row: &[&str]) -> Vec<String> {
    let mut args: Vec<String> = [
        "--comp-based-stats",
        "0",
        "--masking",
        row[6],
        "--soft-masking",
        row[4],
        "--approx-id",
        row[5],
        "--algo",
        row[2],
        "-c",
        "1",
        "-b",
        row[3],
        "--max-hsps",
        "1",
        "-f",
        row[8],
    ]
    .into_iter()
    .map(str::to_owned)
    .collect();
    args.extend(fields(row[9]).iter().map(|s| (*s).to_owned()));
    if row[7] != "default" {
        args.extend(["--motif-masking".into(), row[7].into()]);
    }
    args.extend(row[10].split_whitespace().map(str::to_owned));
    args
}
fn run_case(row: &[&str], dir: &Path, database: &Path, threads: &str, memory: &str) {
    let output = dir.join("output");
    let mut command = Command::new(env!("CARGO_BIN_EXE_diamond"));
    command
        .arg(row[1])
        .arg("-q")
        .arg(fixtures().join(if row[1] == "blastp" {
            "query.faa"
        } else {
            "query.fna"
        }))
        .arg("-d")
        .arg(database)
        .args(flags(row))
        .args(["-p", threads, "--memory-limit", memory])
        .arg("-o")
        .arg(&output);
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{}: {}",
        row[0],
        String::from_utf8_lossy(&result.stderr)
    );
    let actual = std::fs::read(if row[8] == "100" {
        output.with_extension("daa")
    } else {
        output
    })
    .unwrap();
    compare(
        &actual,
        &std::fs::read(fixtures().join(format!("{}.out", row[0]))).unwrap(),
        row[0],
    );
}
#[test]
fn native_seed_masking_and_approximate_identity_match_source() {
    let dir = directory("cli");
    let db = make_database(&dir);
    let text = std::fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    let mut checked = 0;
    for line in text.lines() {
        let row: Vec<_> = line.split('\t').collect();
        assert_eq!(row.len(), 11);
        for (threads, memory) in [("1", "16G"), ("4", "0G")] {
            run_case(&row, &dir, &db, threads, memory);
            checked += 1;
        }
    }
    assert!(checked >= 120);
    std::fs::remove_dir_all(dir).unwrap();
}
fn config(dir: &Path, db: &Path, row: &[&str], threads: i32) -> BlastpConfig {
    BlastpConfig {
        query_files: Vec::new(),
        database: db.to_path_buf(),
        output: Some(dir.join("library.out")),
        matrix: "blosum62".into(),
        gap_open: 11,
        gap_extend: 1,
        max_evalue: 0.001,
        min_bit_score: 0.0,
        max_target_seqs: 25,
        max_hsps: 1,
        report_unaligned: -1,
        include_lineage: false,
        output_header: None,
        compression: Vec::new(),
        file_buffer_size: 67_108_864,
        sam: Default::default(),
        archive: Default::default(),
        tuning: Default::default(),
        scoring_matrix: None,
        auxiliary_output: Default::default(),
        taxonlist: Vec::new(),
        taxon_exclude: Vec::new(),
        seqidlist: PathBuf::new(),
        skip_missing_seqids: false,
        range_culling: false,
        range_cover: 50.0,
        culling_overlap: 50.0,
        taxon_k: 0,
        query_or_subject_cover: 0.0,
        reference_taxon_rank_ids: None,
        ext_chunk_size: 0,
        toppercent: None,
        global_ranking_targets: 0,
        frame_shift: 0,
        extension_mode: None,
        iterate: None,
        round_linear_target: false,
        round_query_skip: None,
        min_id: 0.0,
        approx_min_id: row[5].parse().unwrap(),
        soft_masking: match row[4] {
            "none" | "0" => MaskingAlgo::None,
            "seg" => MaskingAlgo::Seg,
            "tantan" => MaskingAlgo::Tantan,
            _ => panic!("mask"),
        },
        threads,
        index_chunks: 1,
        block_size: Some(row[3].parse().unwrap()),
        query_indexed: row[2] == "1",
        contiguous_seed_mode: false,
        swipe_all: false,
        outfmt: std::iter::once("6")
            .chain(fields(row[9]).iter().copied())
            .map(|value| value.as_bytes().to_vec())
            .collect(),
        sensitivity: diamond::config::Sensitivity::Default,
        masking: MaskingMode::None,
        motif_masking: Vec::new(),
        min_query_len: 0,
        query_cover: 0.0,
        subject_cover: 0.0,
        comp_based_stats: diamond::stats::cbs::CbsMode::Disabled,
        no_self_hits: false,
        ungapped_xdrop_bits: 12.3,
        memory_limit: (threads == 4).then_some(0),
        tmpdir: dir.to_path_buf(),
        target_indexed: false,
        translated_query_layout: None,
    }
}
#[test]
fn supplied_query_controls_preserve_the_native_masking_lifecycle() {
    let dir = directory("library");
    let db = make_database(&dir);
    let text = std::fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    let mut checked = 0;
    for line in text.lines() {
        let row: Vec<_> = line.split('\t').collect();
        if row[3] != "1"
            || row[5] != "75"
            || !["seg", "tantan"].contains(&row[4])
            || row[6] != "0"
            || row[8] != "6"
            || !row[10].is_empty()
        {
            continue;
        }
        for threads in [1, 4] {
            let mut cfg = config(&dir, &db, &row, threads);
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads as usize)
                .build()
                .unwrap();
            if row[1] == "blastp" {
                let queries = fasta::read_fasta_file(
                    &fixtures().join("query.faa"),
                    diamond::basic::value::SequenceType::AminoAcid,
                )
                .unwrap();
                pool.install(|| blastp::run_with_query_records(&cfg, queries))
                    .unwrap();
            } else {
                let queries = fasta::read_fasta_file(
                    &fixtures().join("query.fna"),
                    diamond::basic::value::SequenceType::Nucleotide,
                )
                .unwrap();
                let mut block = Block::new();
                for (id, query) in queries.iter().enumerate() {
                    block
                        .push_back_translated_with_genetic_code(
                            &query.sequence,
                            Some(&query.id),
                            id as u64,
                            0b111111,
                            0,
                            1,
                        )
                        .unwrap();
                }
                cfg.translated_query_layout = Some(TranslatedQueryLayout {
                    mate_indices: vec![None; queries.len()],
                });
                pool.install(|| blastp::run_with_translated_query_block(&cfg, block))
                    .unwrap();
            }
            compare(
                &std::fs::read(cfg.output.unwrap()).unwrap(),
                &std::fs::read(fixtures().join(format!("{}.out", row[0]))).unwrap(),
                row[0],
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 16);
    std::fs::remove_dir_all(dir).unwrap();
}
#[test]
fn search_control_errors_and_impossible_typed_buffer_are_recoverable() {
    let dir = directory("errors");
    let db = make_database(&dir);
    for (options, error) in [
        (
            vec!["--soft-masking", "1"],
            "Invalid value for string field: 1. Permitted values: 0, none, seg, tantan",
        ),
        (
            vec!["--soft-masking", "SEG"],
            "Invalid value for string field: SEG. Permitted values: 0, none, seg, tantan",
        ),
        (
            vec!["--soft-masking", "motif"],
            "Invalid value for string field: motif. Permitted values: 0, none, seg, tantan",
        ),
        (
            vec!["--soft-masking"],
            "Invalid parameter count for option '--soft-masking'",
        ),
        (
            vec!["--soft-masking", "seg", "tantan"],
            "Invalid parameter count for option '--soft-masking'",
        ),
        (
            vec!["--soft-masking=seg"],
            "Invalid option: soft-masking=seg",
        ),
        (
            vec!["--approx-id"],
            "Invalid parameter count for option '--approx-id'",
        ),
        (
            vec!["--approx-id", "75", "90"],
            "Invalid parameter count for option '--approx-id'",
        ),
        (vec!["--approx-id=75"], "Invalid option: approx-id=75"),
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("blastp")
            .arg("-q")
            .arg(fixtures().join("query.faa"))
            .arg("-d")
            .arg(&db)
            .args(options)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(
            String::from_utf8_lossy(&result.stderr).ends_with(&format!("Error: {error}\n")),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    let row = [
        "typed", "blastp", "0", "1", "none", "0", "0", "default", "6", "score", "",
    ];
    let mut cfg = config(&dir, &db, &row, 1);
    cfg.file_buffer_size = usize::MAX;
    let error = blastp::run_with_query_records(&cfg, Vec::new()).unwrap_err();
    assert!(error
        .to_string()
        .contains("Failed to allocate sufficient memory"));
    std::fs::remove_dir_all(dir).unwrap();
}

/// Source SEG can index its alphabet tables with negative SEED_MASK bytes when
/// a previously empty masking table is recomputed. These cases compare the
/// defined native result to two agreeing successful source variants separately
/// from direct conformance with the original commands, which terminate by signal.
#[test]
fn empty_seg_mask_tables_keep_defined_results_across_reference_blocks() {
    let dir = directory("source-ub");
    let db = make_database(&dir);
    let table = std::fs::read_to_string(fixtures().join("compat-cases.tsv")).unwrap();
    let mut checked = 0;
    for line in table.lines() {
        let row: Vec<_> = line.split('\t').collect();
        for (threads, memory) in [("1", "16G"), ("4", "0G")] {
            run_case(&row, &dir, &db, threads, memory);
            checked += 1;
        }
    }
    assert_eq!(checked, 12);
    std::fs::remove_dir_all(dir).unwrap();
}
