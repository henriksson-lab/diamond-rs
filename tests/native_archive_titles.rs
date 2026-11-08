//! Archive title policies preserve the source dictionary and real separators.
#[path = "support/source_startup.rs"]
mod source_startup;

use diamond::commands::blastp::{self, ArchiveOptions, BlastpConfig, TranslatedQueryLayout};
use diamond::data::{block::Block, daa::DaaFile, fasta};
use std::path::{Path, PathBuf};
use std::process::Command;

const FIELDS: &[&str] = &[
    "qseqid",
    "sseqid",
    "sallseqid",
    "stitle",
    "salltitles",
    "bitscore",
];
fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/archive_titles")
}
fn directory(kind: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "diamond-archive-titles-{kind}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&d).unwrap();
    d
}
fn database(d: &Path) -> PathBuf {
    let path = d.join("reference.dmnd");
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
fn expected_metadata(bytes: Vec<u8>, format: &str, invocation: &str, database: &str) -> Vec<u8> {
    if format == "edge" || format == "100" {
        return bytes;
    }
    let mut text = String::from_utf8(bytes).unwrap();
    if format == "sam" {
        let marker = "\tCL:";
        let start = text.find(marker).unwrap() + marker.len();
        let end = start + text[start..].find('\n').unwrap();
        text.replace_range(start..end, invocation);
    }
    if format == "xml" {
        let marker = "<BlastOutput_db>";
        let start = text.find(marker).unwrap() + marker.len();
        let end = start + text[start..].find("</BlastOutput_db>").unwrap();
        text.replace_range(start..end, database);
    }
    text.into_bytes()
}
fn view(archive: &Path, output: &Path, expected: &[u8], label: &str) {
    let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .arg("view")
        .arg("-a")
        .arg(archive)
        .args(["-f", "6"])
        .args(FIELDS)
        .arg("-o")
        .arg(output)
        .args(["-p", "1"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{label}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    compare(&std::fs::read(output).unwrap(), expected, label);
}
#[test]
fn native_archive_titles_and_other_format_flags_match_source() {
    let dir = directory("cli");
    let db = database(&dir);
    let fixture = fixtures();
    let mut checked = 0;
    let mut viewed = 0;
    for line in std::fs::read_to_string(fixture.join("cases.tsv"))
        .unwrap()
        .lines()
    {
        let c: Vec<_> = line.split('\t').collect();
        assert_eq!(c.len(), 6);
        for (threads, memory) in [("1", "16G"), ("4", "0G")] {
            let output = dir.join("output");
            let mut command = Command::new(env!("CARGO_BIN_EXE_diamond"));
            command
                .arg(c[1])
                .arg("-q")
                .arg(fixture.join(if c[1] == "blastp" {
                    "query.faa"
                } else {
                    "query.fna"
                }))
                .arg("-d")
                .arg(&db)
                .args([
                    "--masking",
                    "0",
                    "--comp-based-stats",
                    "1",
                    "--max-hsps",
                    "0",
                    "-b",
                    c[2],
                    "-f",
                    c[4],
                ]);
            if ["6", "104"].contains(&c[4]) {
                command.args(FIELDS);
            }
            if ["titles", "both"].contains(&c[3]) {
                command.arg("--salltitles");
            }
            if ["ids", "both"].contains(&c[3]) {
                command.arg("--sallseqid");
            }
            if c[5] != "default" {
                command.args(["--unal", c[5]]);
            }
            command
                .args(["-p", threads, "--memory-limit", memory])
                .arg("-o")
                .arg(&output);
            let invocation = std::iter::once(command.get_program())
                .chain(command.get_args())
                .map(|s| s.to_string_lossy())
                .collect::<Vec<_>>()
                .join(" ");
            let result = command.output().unwrap();
            assert!(
                result.status.success(),
                "{}: {}",
                c[0],
                String::from_utf8_lossy(&result.stderr)
            );
            let path = if c[4] == "100" {
                output.with_extension("daa")
            } else {
                output
            };
            let actual = std::fs::read(&path).unwrap();
            let expected = expected_metadata(
                std::fs::read(fixture.join(format!("{}.out", c[0]))).unwrap(),
                c[4],
                &invocation,
                &db.to_string_lossy(),
            );
            compare(&actual, &expected, c[0]);
            checked += 1;
            if c[4] == "100" {
                view(
                    &path,
                    &dir.join("view.tsv"),
                    &std::fs::read(fixture.join(format!("{}.view", c[0]))).unwrap(),
                    c[0],
                );
                viewed += 1;
            }
        }
    }
    assert_eq!(checked, 156);
    assert_eq!(viewed, 48);
    std::fs::remove_dir_all(dir).unwrap();
}
fn config(dir: &Path, db: &Path, c: &[&str], threads: i32) -> BlastpConfig {
    BlastpConfig {
        query_files: Vec::new(),
        database: db.to_path_buf(),
        output: Some(dir.join("library.daa")),
        matrix: "blosum62".into(),
        gap_open: 11,
        gap_extend: 1,
        max_evalue: 0.001,
        min_bit_score: 0.0,
        max_target_seqs: 25,
        max_hsps: 0,
        report_unaligned: if c[5] == "default" {
            -1
        } else {
            c[5].parse().unwrap()
        },
        include_lineage: false,
        output_header: None,
        compression: Vec::new(),
        file_buffer_size: 67_108_864,
        sam: Default::default(),
        archive: ArchiveOptions {
            all_subject_titles: ["titles", "both"].contains(&c[3]),
            all_subject_ids: ["ids", "both"].contains(&c[3]),
        },
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
        approx_min_id: 0.0,
        soft_masking: diamond::masking::MaskingAlgo::None,
        threads,
        index_chunks: 0,
        block_size: Some(c[2].parse().unwrap()),
        query_indexed: false,
        contiguous_seed_mode: false,
        swipe_all: false,
        outfmt: vec!["100".into()],
        sensitivity: diamond::config::Sensitivity::Default,
        masking: diamond::masking::MaskingMode::None,
        motif_masking: Vec::new(),
        min_query_len: 0,
        query_cover: 0.0,
        subject_cover: 0.0,
        comp_based_stats: diamond::stats::cbs::CbsMode::Hauser,
        no_self_hits: false,
        ungapped_xdrop_bits: 12.3,
        memory_limit: (threads == 4).then_some(0),
        tmpdir: dir.to_path_buf(),
        target_indexed: false,
        translated_query_layout: None,
    }
}
#[test]
fn supplied_query_apis_preserve_archive_dictionary_policy_and_unaligned_behavior() {
    let dir = directory("library");
    let db = database(&dir);
    let fixture = fixtures();
    let mut checked = 0;
    for line in std::fs::read_to_string(fixture.join("cases.tsv"))
        .unwrap()
        .lines()
    {
        let c: Vec<_> = line.split('\t').collect();
        if c[4] != "100" {
            continue;
        }
        let threads = if c[2] == "1" { 1 } else { 4 };
        let mut cfg = config(&dir, &db, &c, threads);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads as usize)
            .build()
            .unwrap();
        if c[1] == "blastp" {
            let queries = fasta::read_fasta_file(
                &fixture.join("query.faa"),
                diamond::basic::value::SequenceType::AminoAcid,
            )
            .unwrap();
            if c[2] == "1" {
                pool.install(|| blastp::run_with_query_records(&cfg, queries))
                    .unwrap();
            } else {
                let references = fasta::read_fasta_file(
                    &fixture.join("reference.faa"),
                    diamond::basic::value::SequenceType::AminoAcid,
                )
                .unwrap();
                pool.install(|| blastp::run_with_records(&cfg, references, queries))
                    .unwrap();
            }
        } else {
            let queries = fasta::read_fasta_file(
                &fixture.join("query.fna"),
                diamond::basic::value::SequenceType::Nucleotide,
            )
            .unwrap();
            let mut block = Block::new();
            for (id, q) in queries.iter().enumerate() {
                block
                    .push_back_translated_with_genetic_code(
                        &q.sequence,
                        Some(&q.id),
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
        let output = Path::new(cfg.output.as_ref().unwrap());
        compare(
            &std::fs::read(output).unwrap(),
            &std::fs::read(fixture.join(format!("{}.out", c[0]))).unwrap(),
            c[0],
        );
        let daa = DaaFile::open(output).unwrap();
        assert_eq!(daa.db_seqs(), 3);
        assert_eq!(daa.query_records(), if c[1] == "blastp" { 1 } else { 2 });
        if ["titles", "both"].contains(&c[3]) {
            assert!((0..daa.db_seqs_used()).any(|id| daa.ref_name(id as usize).contains('\u{1}')));
            assert!((0..daa.db_seqs_used()).any(|id| daa.ref_name(id as usize).contains('\u{2}')));
        }
        checked += 1;
    }
    assert_eq!(checked, 24);
    std::fs::remove_dir_all(dir).unwrap();
}
#[test]
fn archive_flag_protocol_matches_source() {
    let dir = directory("errors");
    let db = database(&dir);
    for (flags, error) in [
        (
            vec!["--salltitles", "1"],
            "Invalid parameter count for option '--salltitles'",
        ),
        (vec!["--salltitles=1"], "Invalid option: salltitles=1"),
        (
            vec!["--sallseqid", "1"],
            "Invalid parameter count for option '--sallseqid'",
        ),
        (vec!["--sallseqid=1"], "Invalid option: sallseqid=1"),
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("blastp")
            .arg("-q")
            .arg(fixtures().join("query.faa"))
            .arg("-d")
            .arg(&db)
            .args(flags)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(
            String::from_utf8_lossy(&result.stderr).ends_with(
                &source_startup::source_console_text(&format!("Error: {error}\n"))
            ),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    std::fs::remove_dir_all(dir).unwrap();
}
