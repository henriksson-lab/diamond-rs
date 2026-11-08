//! The supplied-query API uses the same block loader and merger as file search.
use diamond::commands::blastp::{self, BlastpConfig, TranslatedQueryLayout};
use diamond::data::{block::Block, fasta};
use std::path::{Path, PathBuf};

fn config(directory: &Path, format: &str, threads: i32) -> BlastpConfig {
    BlastpConfig {
        query_files: Vec::new(),
        database: Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/search_output/reference.faa"),
        output: Some(directory.join("output")),
        matrix: "blosum62".into(),
        gap_open: 11,
        gap_extend: 1,
        max_evalue: 0.001,
        min_bit_score: 0.0,
        max_target_seqs: 25,
        max_hsps: 0,
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
        approx_min_id: 0.0,
        soft_masking: diamond::masking::MaskingAlgo::None,
        threads,
        index_chunks: 0,
        block_size: Some(0.0000001),
        query_indexed: false,
        contiguous_seed_mode: false,
        swipe_all: false,
        outfmt: vec![format.into()],
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
        tmpdir: directory.to_path_buf(),
        target_indexed: false,
        translated_query_layout: None,
    }
}
fn directory(name: &str) -> PathBuf {
    let path =
        std::env::temp_dir().join(format!("diamond-block-api-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&path).unwrap();
    path
}

#[test]
fn supplied_translated_contexts_match_five_frozen_source_formats() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let records = fasta::read_fasta_file(
        &root.join("tests/data/search_output/query.fna"),
        diamond::basic::value::SequenceType::Nucleotide,
    )
    .unwrap();
    let directory = directory("formats");
    for threads in [1, 4] {
        for format in ["6", "0", "xml", "paf", "100"] {
            let mut config = config(&directory, format, threads);
            config.translated_query_layout = Some(TranslatedQueryLayout {
                mate_indices: vec![None; records.len()],
            });
            let mut block = Block::new();
            for (id, record) in records.iter().enumerate() {
                block
                    .push_back_translated_with_genetic_code(
                        &record.sequence,
                        Some(&record.id),
                        id as u64,
                        0b111111,
                        0,
                        1,
                    )
                    .unwrap();
            }
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads as usize)
                .build()
                .unwrap()
                .install(|| blastp::run_with_translated_query_block(&config, block))
                .unwrap();
            let mut actual = std::fs::read(config.output.as_ref().unwrap()).unwrap();
            if format == "xml" {
                actual = String::from_utf8(actual)
                    .unwrap()
                    .replace(
                        &format!(
                            "<BlastOutput_db>{}</BlastOutput_db>",
                            config.database.display()
                        ),
                        "<BlastOutput_db>@DATABASE@</BlastOutput_db>",
                    )
                    .into_bytes();
            }
            let suffix = if format == "100" { "daa" } else { "out" };
            let expected = std::fs::read(
                root.join(format!("tests/data/search_blocks/blastx.{format}.{suffix}")),
            )
            .unwrap();
            assert!(
                actual == expected,
                "{format}/p{threads}: source bytes differ"
            );
        }
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn supplied_translated_blocks_validate_mate_layout_before_database_io() {
    let directory = directory("layout");
    let mut config = config(&directory, "6", 1);
    config.database = "database-that-does-not-exist".into();
    config.translated_query_layout = Some(TranslatedQueryLayout {
        mate_indices: vec![Some(9)],
    });
    let error = blastp::run_with_translated_query_block(&config, Block::new()).unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
    assert!(error
        .to_string()
        .contains("Invalid translated query block layout"));
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn supplied_paired_contexts_retain_genetic_code_orfs_and_mate_source_sequences() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let records = fasta::read_fasta_file(
        &root.join("tests/data/search_output/query.fna"),
        diamond::basic::value::SequenceType::Nucleotide,
    )
    .unwrap();
    let directory = directory("paired");
    let expected =
        std::fs::read(root.join("tests/data/search_blocks/library.paired-gencode.tsv")).unwrap();
    for threads in [1, 4] {
        let mut config = config(&directory, "6", threads);
        config.outfmt = [
            "6",
            "qseqid",
            "sseqid",
            "qstart",
            "qend",
            "qframe",
            "full_qseq_mate",
        ]
        .map(|value| value.as_bytes().to_vec())
        .to_vec();
        let mut block = Block::new();
        for (left, right) in records.iter().zip(records.iter().rev()) {
            for record in [left, right] {
                let oid = block.source_seqs().len() as u64;
                block
                    .push_back_translated_with_genetic_code(
                        &record.sequence,
                        Some(&record.id),
                        oid,
                        0b111111,
                        1,
                        4,
                    )
                    .unwrap();
            }
        }
        config.translated_query_layout = Some(TranslatedQueryLayout {
            mate_indices: (0..records.len() * 2).map(|id| Some(id ^ 1)).collect(),
        });
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads as usize)
            .build()
            .unwrap()
            .install(|| blastp::run_with_translated_query_block(&config, block))
            .unwrap();
        assert_eq!(
            std::fs::read(config.output.as_ref().unwrap()).unwrap(),
            expected,
            "p{threads}"
        );
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn supplied_translated_contexts_use_native_sam_callbacks() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let inputs = root.join("tests/data/tabular_context");
    let records = fasta::read_fasta_file(
        &inputs.join("query.fna"),
        diamond::basic::value::SequenceType::Nucleotide,
    )
    .unwrap();
    let directory = directory("sam");
    let frozen = std::fs::read_to_string(root.join("tests/data/sam_search/case-041.sam")).unwrap();
    for threads in [1, 4] {
        let mut config = config(&directory, "sam", threads);
        config.database = inputs.join("reference.faa");
        config.sam = blastp::SamOptions {
            invocation: "supplied translated query block".into(),
            all_subject_titles: true,
            query_len: true,
        };
        config.translated_query_layout = Some(TranslatedQueryLayout {
            mate_indices: vec![None; records.len()],
        });
        let mut block = Block::new();
        for (id, record) in records.iter().enumerate() {
            block
                .push_back_translated_with_genetic_code(
                    &record.sequence,
                    Some(&record.id),
                    id as u64,
                    0b111111,
                    0,
                    1,
                )
                .unwrap();
        }
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads as usize)
            .build()
            .unwrap()
            .install(|| blastp::run_with_translated_query_block(&config, block))
            .unwrap();
        let (prefix, after) = frozen.split_once("\tCL:").unwrap();
        let (_, suffix) = after.split_once('\n').unwrap();
        let expected = format!(
            "{prefix}\tCL:{}\n{suffix}",
            std::str::from_utf8(&config.sam.invocation).unwrap()
        );
        assert_eq!(
            std::fs::read(config.output.unwrap()).unwrap(),
            expected.as_bytes()
        );
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn supplied_protein_and_translated_queries_preserve_binary_edge_records() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let inputs = root.join("tests/data/tabular_context");
    let directory = directory("edge");
    for threads in [1, 4] {
        for (translated, fixture) in [(false, "case-003"), (true, "case-027")] {
            let mut config = config(&directory, "edge", threads);
            config.database = inputs.join("reference.faa");
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads as usize)
                .build()
                .unwrap();
            if translated {
                let records = fasta::read_fasta_file(
                    &inputs.join("query.fna"),
                    diamond::basic::value::SequenceType::Nucleotide,
                )
                .unwrap();
                config.translated_query_layout = Some(TranslatedQueryLayout {
                    mate_indices: vec![None; records.len()],
                });
                let mut block = Block::new();
                for (id, record) in records.iter().enumerate() {
                    block
                        .push_back_translated_with_genetic_code(
                            &record.sequence,
                            Some(&record.id),
                            id as u64,
                            0b111111,
                            0,
                            1,
                        )
                        .unwrap();
                }
                pool.install(|| blastp::run_with_translated_query_block(&config, block))
                    .unwrap();
            } else {
                let records = fasta::read_fasta_amino_acid(
                    std::fs::File::open(inputs.join("query.faa")).unwrap(),
                )
                .unwrap();
                pool.install(|| blastp::run_with_query_records(&config, records))
                    .unwrap();
            }
            assert_eq!(
                std::fs::read(config.output.unwrap()).unwrap(),
                std::fs::read(root.join(format!("tests/data/edge_null/{fixture}.out"))).unwrap()
            );
        }
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn supplied_queries_use_taxonomy_defaults_and_explicit_overrides() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/taxon_search");
    let inputs = root.join("tests/data/tabular_context");
    let directory = directory("taxonomy");
    let database = directory.join("weak.dmnd");
    let result = std::process::Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args(["makedb", "--in"])
        .arg(fixtures.join("weak.faa"))
        .arg("-d")
        .arg(&database)
        .arg("--taxonmap")
        .arg(fixtures.join("weak-mapping.tsv"))
        .arg("--taxonnodes")
        .arg(fixtures.join("root-nodes.dmp"))
        .arg("--taxonnames")
        .arg(fixtures.join("root-names.dmp"))
        .arg("-p1")
        .output()
        .unwrap();
    assert!(result.status.success());
    for threads in [1, 4] {
        for (translated, command) in [(false, "blastp"), (true, "blastx")] {
            for mode in ["default", "all", "score"] {
                let mut config = config(&directory, "102", threads);
                config.database = database.clone();
                config.include_lineage = true;
                config.max_hsps = 1;
                if mode == "all" {
                    config.toppercent = Some(100.0);
                }
                if mode == "score" {
                    config.min_bit_score = 80.0;
                }
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads as usize)
                    .build()
                    .unwrap();
                if translated {
                    let records = fasta::read_fasta_file(
                        &inputs.join("query.fna"),
                        diamond::basic::value::SequenceType::Nucleotide,
                    )
                    .unwrap();
                    config.translated_query_layout = Some(TranslatedQueryLayout {
                        mate_indices: vec![None; records.len()],
                    });
                    let mut block = Block::new();
                    for (id, record) in records.iter().enumerate() {
                        block
                            .push_back_translated_with_genetic_code(
                                &record.sequence,
                                Some(&record.id),
                                id as u64,
                                0b111111,
                                0,
                                1,
                            )
                            .unwrap();
                    }
                    pool.install(|| blastp::run_with_translated_query_block(&config, block))
                        .unwrap();
                } else {
                    let records = fasta::read_fasta_amino_acid(
                        std::fs::File::open(inputs.join("query.faa")).unwrap(),
                    )
                    .unwrap();
                    pool.install(|| blastp::run_with_query_records(&config, records))
                        .unwrap();
                }
                assert_eq!(
                    std::fs::read(config.output.unwrap()).unwrap(),
                    std::fs::read(fixtures.join(format!("weak.{command}.b0.0000001.{mode}.tsv")))
                        .unwrap()
                );
            }
        }
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn supplied_queries_match_complete_gzip_headers_with_fixed_invocation() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/output_options");
    let inputs = root.join("tests/data/tabular_context");
    let directory = directory("gzip-headers");
    let invocations = std::fs::read_to_string(fixtures.join("api-invocations.tsv")).unwrap();
    for threads in [1, 4] {
        for (translated, command) in [(false, "blastp"), (true, "blastx")] {
            for format in ["6", "sam"] {
                let label = format!("api.{command}.{format}");
                let invocation = invocations
                    .lines()
                    .find_map(|row| {
                        let (k, v) = row.split_once('\t').unwrap();
                        (k == label).then_some(v)
                    })
                    .unwrap();
                let mut config = config(&directory, format, threads);
                config.database = inputs.join("reference.faa");
                config.output_header = Some(vec!["verbose".into()]);
                config.compression = "1".into();
                config.file_buffer_size = 7;
                config.sam.invocation = invocation.into();
                if format == "6" {
                    config.outfmt.extend(
                        [
                            "qseqid",
                            "sseqid",
                            "qseq",
                            "sseq",
                            "cigar",
                            "full_qseq",
                            "full_qqual",
                            "hspnum",
                            "bitscore",
                        ]
                        .into_iter()
                        .map(|value| value.as_bytes().to_vec()),
                    );
                }
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads as usize)
                    .build()
                    .unwrap();
                if translated {
                    let records = fasta::read_fasta_file(
                        &inputs.join("query.fna"),
                        diamond::basic::value::SequenceType::Nucleotide,
                    )
                    .unwrap();
                    config.translated_query_layout = Some(TranslatedQueryLayout {
                        mate_indices: vec![None; records.len()],
                    });
                    let mut block = Block::new();
                    for (id, record) in records.iter().enumerate() {
                        block
                            .push_back_translated_with_genetic_code(
                                &record.sequence,
                                Some(&record.id),
                                id as u64,
                                0b111111,
                                0,
                                1,
                            )
                            .unwrap();
                    }
                    pool.install(|| blastp::run_with_translated_query_block(&config, block))
                        .unwrap();
                } else {
                    let records = fasta::read_fasta_amino_acid(
                        std::fs::File::open(inputs.join("query.faa")).unwrap(),
                    )
                    .unwrap();
                    pool.install(|| blastp::run_with_query_records(&config, records))
                        .unwrap();
                }
                let actual = std::fs::read(config.output.unwrap()).unwrap();
                let expected = std::fs::read(fixtures.join(format!("{label}.gzip"))).unwrap();
                assert!(
                    actual == expected,
                    "{label}/p{threads}: complete gzip differs"
                );
            }
        }
    }
    std::fs::remove_dir_all(directory).unwrap();
}
