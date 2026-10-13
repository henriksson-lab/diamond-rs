use std::env;
use std::path::PathBuf;

use diamond::cluster::cascaded::GraphAlgo;
use diamond::commands::blastp::BlastpConfig;
use diamond::commands::blastx::BlastxConfig;
use diamond::commands::cluster_cmd::{ClusterConfig, ClusterWorkflow};
use diamond::commands::view::ViewConfig;
use diamond::config::Sensitivity;
use diamond::output::format::FieldId;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        print_usage();
        return;
    }

    // Reject getseq options that standard upstream either disables or limits
    // to another workflow instead of silently accepting them on the native
    // route.
    if args[1] == "getseq" {
        if let Some(error) = native_getseq_option_error(&args) {
            eprintln!("Error: {error}");
            std::process::exit(1);
        }
    }

    // Match the upstream command-line parser's workflow checks before the
    // native view route can accidentally accept and ignore search-only flags.
    let native_view = if args[1] == "view" && !has_flag(&args, "--legacy") {
        match parse_native_view_cli(&args) {
            Ok(config) => config,
            Err(error) => {
                eprintln!("Error: {error}");
                std::process::exit(1);
            }
        }
    } else {
        NativeViewCli::default()
    };

    if args[1] == "view" && args.len() == 2 && !has_flag(&args, "--legacy") {
        print_view_documentation();
    }

    match args[1].as_str() {
        "version" => {
            println!("diamond version 2.1.24");
        }
        "help" | "--help" | "-h" => {
            print_usage();
        }
        "dbinfo" => {
            if let Some(db_path) = get_arg(&args, &["-d", "--db"]) {
                print_banner();
                run_or_exit(diamond::commands::dbinfo::run(&db_path));
            } else {
                eprintln!("Error: -d/--db argument required");
                std::process::exit(1);
            }
        }
        "makedb" => {
            // Upstream assigns each repeated `--in` occurrence to the option
            // vector, so the last occurrence replaces earlier ones. Multiple
            // values in that final occurrence are still rejected by makedb.
            let input_files = get_last_arg_group(&args, &["--in"]);
            if let Some(db) = get_arg(&args, &["-d", "--db"]) {
                if !input_files.is_empty() {
                    print_banner();
                    // Upstream defaults to `std::thread::hardware_concurrency()`;
                    // zero lets Rayon select the host's available parallelism.
                    let threads = parse_arg_or(&args, &["-p", "--threads"], 0);
                    run_or_exit(diamond::commands::makedb::run_with_options(
                        &input_files,
                        &db,
                        threads,
                        &diamond::commands::makedb::MakedbOptions {
                            taxonmap: get_arg(&args, &["--taxonmap"]),
                            taxonnodes: get_arg(&args, &["--taxonnodes"]),
                            taxonnames: get_arg(&args, &["--taxonnames"]),
                            no_parse_seqids: has_flag(&args, "--no-parse-seqids"),
                        },
                    ));
                } else {
                    eprintln!("Error: --in argument required");
                    std::process::exit(1);
                }
            } else {
                eprintln!("Error: -d/--db argument required");
                std::process::exit(1);
            }
        }
        "getseq" => {
            if let Some(db) = get_arg(&args, &["-d", "--db"]) {
                let seq_numbers = get_all_args(&args, &["--seq"]);
                let seq_numbers = (!seq_numbers.is_empty()).then(|| seq_numbers.join(","));
                let output = get_arg(&args, &["-o", "--out"]).map(PathBuf::from);
                run_or_exit(diamond::commands::getseq::run_with_output(
                    &db,
                    seq_numbers.as_deref(),
                    output.as_deref(),
                ));
            } else {
                eprintln!("Error: -d/--db argument required");
                std::process::exit(1);
            }
        }
        "makeidx" => {
            if let Some(database) = get_arg(&args, &["-d", "--db"]) {
                print_banner();
                let mut config =
                    diamond::data::index::MakeIndexConfig::new(database, parse_sensitivity(&args));
                config.shape_mask = get_all_args(&args, &["--shape-mask"]);
                config.shapes = parse_arg_or(&args, &["--shapes"], 0);
                config.threads = parse_arg_or(&args, &["-p", "--threads"], 0);
                match diamond::data::index::make_index(&config) {
                    Ok(path) => eprintln!("Wrote seed index: {}", path.display()),
                    Err(e) => run_or_exit(Err(e)),
                }
            } else {
                eprintln!("Error: -d/--db argument required");
                std::process::exit(1);
            }
        }
        "blastp" | "blastx" if has_flag(&args, "--help") || has_flag(&args, "-h") => {
            print_usage();
        }
        "blastp" if !has_flag(&args, "--legacy") && !route_blastp_to_legacy(&args) => {
            // Native Rust blastp pipeline
            print_banner();
            let (output, outfmt) = parse_native_search_output(&args).unwrap_or_else(|error| {
                eprintln!("Error: {error}");
                std::process::exit(1);
            });
            let config = BlastpConfig {
                query_files: get_all_args(&args, &["-q", "--query"]),
                database: get_arg(&args, &["-d", "--db"]).unwrap_or_default(),
                output,
                matrix: get_arg(&args, &["--matrix"]).unwrap_or_else(|| "blosum62".into()),
                gap_open: parse_arg_or(&args, &["--gapopen"], -1),
                gap_extend: parse_arg_or(&args, &["--gapextend"], -1),
                max_evalue: parse_arg_or(&args, &["-e", "--evalue"], 0.001),
                max_target_seqs: parse_max_target_seqs(&args),
                ext_chunk_size: parse_arg_or(&args, &["--ext-chunk-size"], 0),
                toppercent: get_arg(&args, &["--top"]).and_then(|s| s.parse().ok()),
                global_ranking_targets: parse_arg_or(&args, &["--global-ranking"], 0),
                min_id: parse_arg_or(&args, &["--id"], 0.0),
                // C++ defaults `--threads` to `std::thread::hardware_concurrency()`
                // (`config.cpp:783`). Match that — `0` here is interpreted by
                // `rayon::ThreadPoolBuilder` as "all cores".
                threads: parse_arg_or(&args, &["-p", "--threads"], 0),
                outfmt,
                sensitivity: parse_sensitivity(&args),
                masking: diamond::masking::MaskingMode::parse(
                    &get_arg(&args, &["--masking"]).unwrap_or_else(|| "tantan".into()),
                )
                .unwrap_or_else(|e| {
                    eprintln!("Error: {e}");
                    std::process::exit(1);
                }),
                motif_masking: get_arg(&args, &["--motif-masking"]).unwrap_or_default(),
                min_query_len: parse_arg_or(&args, &["--min-query-len"], 0usize),
                query_cover: parse_arg_or(&args, &["--query-cover"], 0.0),
                subject_cover: parse_arg_or(&args, &["--subject-cover"], 0.0),
                comp_based_stats: diamond::stats::cbs::CbsMode::parse(
                    &get_arg(&args, &["--comp-based-stats"]).unwrap_or_else(|| "1".into()),
                ),
                no_self_hits: has_flag(&args, "--no-self-hits"),
                ungapped_xdrop_bits: parse_arg_or(&args, &["-x", "--xdrop"], 12.3),
                memory_limit: parse_memory_limit_arg(&args),
                tmpdir: get_arg(&args, &["-t", "--tmpdir"])
                    .map(PathBuf::from)
                    .unwrap_or_default(),
                target_indexed: parse_target_indexed(&args),
                translated_query_layout: None,
            };
            run_or_exit(diamond::commands::blastp::run(&config));
        }
        "blastx" if !has_flag(&args, "--legacy") && !route_blastx_to_legacy(&args) => {
            // Upstream rejects this combination during configuration rather
            // than silently treating translated protein contexts as ordinary
            // self-search queries.
            if has_flag(&args, "--no-self-hits") {
                print_banner();
                eprintln!("Error: --no-self-hits option is not supported in blastx mode.");
                std::process::exit(1);
            }
            print_banner();
            let (output, outfmt) = parse_native_search_output(&args).unwrap_or_else(|error| {
                eprintln!("Error: {error}");
                std::process::exit(1);
            });
            let config = BlastxConfig {
                query_files: get_all_args(&args, &["-q", "--query"]),
                database: get_arg(&args, &["-d", "--db"]).unwrap_or_default(),
                output,
                matrix: get_arg(&args, &["--matrix"]).unwrap_or_else(|| "blosum62".into()),
                gap_open: parse_arg_or(&args, &["--gapopen"], -1),
                gap_extend: parse_arg_or(&args, &["--gapextend"], -1),
                max_evalue: parse_arg_or(&args, &["-e", "--evalue"], 0.001),
                max_target_seqs: parse_max_target_seqs(&args),
                ext_chunk_size: parse_arg_or(&args, &["--ext-chunk-size"], 0),
                toppercent: get_arg(&args, &["--top"]).and_then(|s| s.parse().ok()),
                global_ranking_targets: parse_arg_or(&args, &["--global-ranking"], 0),
                min_id: parse_arg_or(&args, &["--id"], 0.0),
                // C++ defaults `--threads` to `std::thread::hardware_concurrency()`.
                threads: parse_arg_or(&args, &["-p", "--threads"], 0),
                outfmt,
                sensitivity: parse_sensitivity(&args),
                query_gencode: parse_arg_or(&args, &["--query-gencode"], 1),
                strand: get_arg(&args, &["--strand"]).unwrap_or_else(|| "both".into()),
                min_orf: get_arg(&args, &["--min-orf"]).and_then(|s| s.parse().ok()),
                // Forward `--masking` / `--motif-masking` to the inner blastp
                // pipeline. Without this blastx silently overrode the user's
                // choice to `Tantan`/empty.
                masking: diamond::masking::MaskingMode::parse(
                    &get_arg(&args, &["--masking"]).unwrap_or_else(|| "tantan".into()),
                )
                .unwrap_or_else(|e| {
                    eprintln!("Error: {e}");
                    std::process::exit(1);
                }),
                motif_masking: get_arg(&args, &["--motif-masking"]).unwrap_or_default(),
                query_cover: parse_arg_or(&args, &["--query-cover"], 0.0),
                subject_cover: parse_arg_or(&args, &["--subject-cover"], 0.0),
                comp_based_stats: diamond::stats::cbs::CbsMode::parse(
                    &get_arg(&args, &["--comp-based-stats"]).unwrap_or_else(|| "1".into()),
                ),
                no_self_hits: has_flag(&args, "--no-self-hits"),
                ungapped_xdrop_bits: parse_arg_or(&args, &["-x", "--xdrop"], 12.3),
                memory_limit: parse_memory_limit_arg(&args),
                tmpdir: get_arg(&args, &["-t", "--tmpdir"])
                    .map(PathBuf::from)
                    .unwrap_or_default(),
                target_indexed: parse_target_indexed(&args),
            };
            run_or_exit(diamond::commands::blastx::run(&config));
        }
        "cluster" | "linclust" | "deepclust" if !has_flag(&args, "--legacy") => {
            print_banner();
            if let Err(error) = validate_native_cluster_options(&args) {
                eprintln!("Error: {error}");
                std::process::exit(1);
            }
            let config = match parse_native_cluster_config(&args) {
                Ok(config) => config,
                Err(error) => {
                    eprintln!("Error: {error}");
                    std::process::exit(1);
                }
            };
            run_or_exit(diamond::commands::cluster_cmd::run(&config));
        }
        "merge-daa" if !has_flag(&args, "--legacy") => {
            print_banner();
            let input_files = get_all_args(&args, &["--in"]);
            let output = get_arg(&args, &["-o", "--out"]);
            if input_files.is_empty() {
                eprintln!("Error: --in argument required");
                std::process::exit(1);
            }
            if let Some(output) = output {
                run_or_exit(diamond::commands::merge_daa::run(&input_files, &output));
            } else {
                eprintln!("Error: -o/--out argument required");
                std::process::exit(1);
            }
        }
        "view"
            if !has_flag(&args, "--legacy")
                && native_view
                    .outfmt
                    .first()
                    .is_some_and(|f| f == "100" || f == "daa") =>
        {
            present_native_view_startup(&args, &native_view);
            let daa_file = native_view.daa_file.clone();
            let output = native_view.output.clone();
            if daa_file.is_empty() {
                eprintln!("Error: The view command requires a DAA (option -a) input file.");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_daa(&ViewConfig {
                daa_file,
                output,
                outfmt: native_view.outfmt.clone(),
                max_target_seqs: native_view.max_target_seqs,
                toppercent: native_view.toppercent,
                forward_only: native_view.forward_only,
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: native_view.header.clone(),
                threads: native_view.threads,
                qnum_offset: native_view.qnum_offset,
                snum_offset: native_view.snum_offset,
                verbosity: native_view.verbosity,
            }));
        }
        "view"
            if !has_flag(&args, "--legacy") && {
                let outfmt = &native_view.outfmt;
                outfmt.is_empty()
                    || outfmt
                        .first()
                        .is_some_and(|f| f == "6" || f == "tab" || f == "104" || f == "json-flat")
            } =>
        {
            present_native_view_startup(&args, &native_view);
            let daa_file = native_view.daa_file.clone();
            let output = native_view.output.clone();
            if daa_file.is_empty() {
                eprintln!("Error: The view command requires a DAA (option -a) input file.");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_tabular(&ViewConfig {
                daa_file,
                output,
                outfmt: native_view.outfmt.clone(),
                max_target_seqs: native_view.max_target_seqs,
                toppercent: native_view.toppercent,
                forward_only: native_view.forward_only,
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: native_view.header.clone(),
                threads: native_view.threads,
                qnum_offset: native_view.qnum_offset,
                snum_offset: native_view.snum_offset,
                verbosity: native_view.verbosity,
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && native_view
                    .outfmt
                    .first()
                    .is_some_and(|f| f == "103" || f == "paf") =>
        {
            present_native_view_startup(&args, &native_view);
            let daa_file = native_view.daa_file.clone();
            let output = native_view.output.clone();
            if daa_file.is_empty() {
                eprintln!("Error: The view command requires a DAA (option -a) input file.");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_paf(&ViewConfig {
                daa_file,
                output,
                outfmt: native_view.outfmt.clone(),
                max_target_seqs: native_view.max_target_seqs,
                toppercent: native_view.toppercent,
                forward_only: native_view.forward_only,
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: native_view.header.clone(),
                threads: native_view.threads,
                qnum_offset: native_view.qnum_offset,
                snum_offset: native_view.snum_offset,
                verbosity: native_view.verbosity,
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && native_view
                    .outfmt
                    .first()
                    .is_some_and(|f| f == "101" || f == "sam") =>
        {
            present_native_view_startup(&args, &native_view);
            let daa_file = native_view.daa_file.clone();
            let output = native_view.output.clone();
            if daa_file.is_empty() {
                eprintln!("Error: The view command requires a DAA (option -a) input file.");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_sam(&ViewConfig {
                daa_file,
                output,
                outfmt: native_view.outfmt.clone(),
                max_target_seqs: native_view.max_target_seqs,
                toppercent: native_view.toppercent,
                forward_only: native_view.forward_only,
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: native_view.header.clone(),
                threads: native_view.threads,
                qnum_offset: native_view.qnum_offset,
                snum_offset: native_view.snum_offset,
                verbosity: native_view.verbosity,
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && native_view.outfmt.first().is_some_and(|f| f == "0") =>
        {
            present_native_view_startup(&args, &native_view);
            let daa_file = native_view.daa_file.clone();
            let output = native_view.output.clone();
            if daa_file.is_empty() {
                eprintln!("Error: The view command requires a DAA (option -a) input file.");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_pairwise(&ViewConfig {
                daa_file,
                output,
                outfmt: native_view.outfmt.clone(),
                max_target_seqs: native_view.max_target_seqs,
                toppercent: native_view.toppercent,
                forward_only: native_view.forward_only,
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: native_view.header.clone(),
                threads: native_view.threads,
                qnum_offset: native_view.qnum_offset,
                snum_offset: native_view.snum_offset,
                verbosity: native_view.verbosity,
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && native_view
                    .outfmt
                    .first()
                    .is_some_and(|f| f == "5" || f == "xml") =>
        {
            present_native_view_startup(&args, &native_view);
            let daa_file = native_view.daa_file.clone();
            let output = native_view.output.clone();
            if daa_file.is_empty() {
                eprintln!("Error: The view command requires a DAA (option -a) input file.");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_xml(&ViewConfig {
                daa_file,
                output,
                outfmt: native_view.outfmt.clone(),
                max_target_seqs: native_view.max_target_seqs,
                toppercent: native_view.toppercent,
                forward_only: native_view.forward_only,
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: native_view.header.clone(),
                threads: native_view.threads,
                qnum_offset: native_view.qnum_offset,
                snum_offset: native_view.snum_offset,
                verbosity: native_view.verbosity,
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && native_view.outfmt.first().is_some_and(|f| f == "null") =>
        {
            present_native_view_startup(&args, &native_view);
            let daa_file = native_view.daa_file.clone();
            let output = native_view.output.clone();
            if daa_file.is_empty() {
                eprintln!("Error: The view command requires a DAA (option -a) input file.");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_null(&ViewConfig {
                daa_file,
                output,
                outfmt: native_view.outfmt.clone(),
                max_target_seqs: native_view.max_target_seqs,
                toppercent: native_view.toppercent,
                forward_only: native_view.forward_only,
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: native_view.header.clone(),
                threads: native_view.threads,
                qnum_offset: native_view.qnum_offset,
                snum_offset: native_view.snum_offset,
                verbosity: native_view.verbosity,
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && native_view.outfmt.first().is_some_and(|f| f == "edge") =>
        {
            present_native_view_startup(&args, &native_view);
            let daa_file = native_view.daa_file.clone();
            let output = native_view.output.clone();
            if daa_file.is_empty() {
                eprintln!("Error: The view command requires a DAA (option -a) input file.");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_edge(&ViewConfig {
                daa_file,
                output,
                outfmt: native_view.outfmt.clone(),
                max_target_seqs: native_view.max_target_seqs,
                toppercent: native_view.toppercent,
                forward_only: native_view.forward_only,
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: native_view.header.clone(),
                threads: native_view.threads,
                qnum_offset: native_view.qnum_offset,
                snum_offset: native_view.snum_offset,
                verbosity: native_view.verbosity,
            }));
        }
        "test" => {
            print_banner();
            run_or_exit(diamond::commands::test_cmd::run());
        }
        _ => {
            // Fall back to C++ FFI for full compatibility
            // Filter out --legacy flag which is not known to C++
            let filtered: Vec<&str> = args
                .iter()
                .map(|s| s.as_str())
                .filter(|s| *s != "--legacy")
                .collect();
            let code = run_legacy(&filtered);
            std::process::exit(code);
        }
    }
}

fn print_banner() {
    eprintln!(
        "diamond v2.1.24.178 (C) Max Planck Society for the Advancement of Science, Benjamin J. Buchfink, University of Tuebingen"
    );
    eprintln!("Documentation, support and updates available at http://www.diamondsearch.org");
    eprintln!("Please cite: http://dx.doi.org/10.1038/s41592-021-01101-x Nature Methods (2021)");
    eprintln!();
}

fn present_native_view_startup(args: &[String], config: &NativeViewCli) {
    if config.verbosity == 0 {
        return;
    }
    print_banner();
    if config.verbosity == 3 {
        eprintln!("{}", args.join(" "));
    }
    if args.len() != 2 {
        let threads = if config.threads == 0 {
            std::thread::available_parallelism().map_or(1, usize::from)
        } else {
            config.threads
        };
        eprintln!("#CPU threads: {threads}");
    }
    if config.verbosity >= 2 && cfg!(debug_assertions) {
        eprintln!("Assertions enabled.");
    }
    if config.verbosity == 3 {
        // These are the translated equivalents of the compile/configuration
        // diagnostics emitted by the pinned upstream build.
        eprintln!("TRACE_PT_MEMBUF=OFF");
        eprintln!("MAX_SHAPE_LEN=19 SEQ_MASK STRICT_BAND");
    }
}

fn print_usage() {
    println!("diamond v2.1.24 — Rust port");
    println!();
    println!("Commands:");
    println!("  makedb     Build DIAMOND database from FASTA");
    println!("  blastp     Protein-protein alignment");
    println!("  blastx     Translated DNA-protein alignment");
    println!("  view       View DAA file");
    println!("  dbinfo     Print database info");
    println!("  getseq     Retrieve sequences from database");
    println!("  makeidx    Build persistent seed index");
    println!("  cluster    Cluster sequences");
    println!("  merge-daa  Merge DAA files");
    println!("  version    Show version");
    println!("  test       Run regression tests");
    println!();
    println!("Use 'diamond COMMAND --help' for command-specific options.");
    println!(
        "Native blastp/blastx: --memory-limit 16G is the default soft RSS ceiling; 0G forces disk mode."
    );
    #[cfg(all(feature = "ffi", diamond_ffi_available, not(windows)))]
    println!("Add --legacy to blastp/blastx to use C++ FFI backend.");
    #[cfg(not(all(feature = "ffi", diamond_ffi_available, not(windows))))]
    println!("C++ FFI fallback is not available in this build.");
}

/// Print the workflow-specific documentation emitted by upstream for
/// `diamond view` with no options. Keep this generated from the translated
/// output-field metadata so adding a field cannot leave the help stale.
fn print_view_documentation() {
    fn option(name: &str, description: &str) {
        println!("{name:<25}{description}");
    }

    println!("Options:");
    option("--threads", "number of CPU threads");
    option("--verbose", "verbose console output");
    option("--log", "enable debug log");
    option("--quiet", "disable console output");
    option("--tmpdir", "directory for temporary files");
    option("--out", "output file");
    option(
        "--header",
        "Use header lines in tabular output format (0/simple/verbose).",
    );
    option(
        "--max-target-seqs",
        "maximum number of target sequences to report alignments for (default=25)",
    );
    option(
        "--top",
        "report alignments within this percentage range of top alignment score (overrides --max-target-seqs)",
    );

    let base_fields = [
        "qseqid",
        "qlen",
        "sseqid",
        "sallseqid",
        "slen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "qseq",
        "sseq",
        "evalue",
        "bitscore",
        "score",
        "length",
        "pident",
        "nident",
        "mismatch",
        "positive",
        "gapopen",
        "gaps",
        "ppos",
        "qframe",
        "btop",
        "staxids",
        "sscinames",
        "sskingdoms",
        "stitle",
        "salltitles",
        "qcovhsp",
        "qtitle",
        "full_sseq",
        "qqual",
        "qnum",
        "snum",
        "scovhsp",
        "full_qqual",
        "full_qseq",
        "qseq_gapped",
        "sseq_gapped",
        "qstrand",
        "cigar",
        "skingdoms",
        "sphylums",
        "full_qseq_mate",
        "qseq_translated",
        "hspnum",
        "normalized_bitscore",
        "normalized_nident",
        "approx_pident",
        "corrected_bitscore",
        "slineages",
    ];
    let mut fields = base_fields
        .iter()
        .map(|key| {
            let field = FieldId::from_name(key)
                .and_then(FieldId::output_field)
                .expect("documented output field must be registered");
            (field.key.to_string(), field.description.to_string())
        })
        .collect::<Vec<_>>();
    for rank in diamond::data::taxonomy_nodes::Rank::NAMES.iter().skip(1) {
        fields.push((
            format!("s{}", rank.replace(' ', "_")),
            format!("Unique subject {rank}(s), separated by a ';'"),
        ));
    }
    let width = fields.iter().map(|(key, _)| key.len()).max().unwrap_or(0);
    let mut format_description = String::from(
        "output format\n\
\t0   = BLAST pairwise\n\
\t5   = BLAST XML\n\
\t6   = BLAST tabular\n\
\t100 = DIAMOND alignment archive (DAA)\n\
\t101 = SAM\n\
\t102 = Taxonomic classification\n\
\t103 = PAF\n\
\t104 = JSON (flat)\n\n\
\tValues 6 and 104 may be followed by a space-separated list of these keywords:\n\n",
    );
    for (key, description) in fields {
        format_description.push_str(&format!("\t{key:>width$} means {description}\n"));
    }
    format_description.push_str(
        "\n\tDefault: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\n",
    );
    option("--outfmt", &format_description);
    option(
        "--qnum-offset",
        "offset added to query ordinal id (qnum field)",
    );
    option(
        "--snum-offset",
        "offset added to subject ordinal id (snum field)",
    );
    option("--daa", "DIAMOND alignment archive (DAA) file");
    option("--forwardonly", "only show alignments of forward strand");
    println!();
}

#[cfg(all(feature = "ffi", diamond_ffi_available, not(windows)))]
fn run_legacy(args: &[&str]) -> i32 {
    diamond::ffi::run(args)
}

#[cfg(not(all(feature = "ffi", diamond_ffi_available, not(windows))))]
fn run_legacy(_args: &[&str]) -> i32 {
    eprintln!("Error: C++ FFI fallback is not available in this build.");
    eprintln!(
        "A non-Windows source checkout with the adjacent `diamond/` tree is required to enable it."
    );
    1
}

fn run_or_exit(result: std::io::Result<()>) {
    if let Err(e) = result {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

fn get_arg(args: &[String], flags: &[&str]) -> Option<String> {
    for (i, arg) in args.iter().enumerate() {
        if flags.contains(&arg.as_str()) {
            return args.get(i + 1).cloned();
        }
        for flag in flags {
            if flag.starts_with("--") {
                let Some(value) = arg.strip_prefix(&format!("{flag}=")) else {
                    continue;
                };
                return Some(value.to_string());
            }
            if flag.starts_with('-') && flag.len() == 2 {
                let Some(value) = arg.strip_prefix(flag) else {
                    continue;
                };
                if !value.is_empty() {
                    return Some(value.to_string());
                }
            }
        }
    }
    None
}

fn get_all_args(args: &[String], flags: &[&str]) -> Vec<String> {
    let mut values = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let mut joined = None;
        if flags.contains(&args[i].as_str()) {
            i += 1;
        } else {
            for flag in flags {
                if flag.starts_with("--") {
                    joined = args[i]
                        .strip_prefix(&format!("{flag}="))
                        .map(ToOwned::to_owned);
                } else if flag.starts_with('-') && flag.len() == 2 {
                    joined = args[i]
                        .strip_prefix(flag)
                        .filter(|value| !value.is_empty())
                        .map(ToOwned::to_owned);
                }
                if joined.is_some() {
                    break;
                }
            }
            if joined.is_none() {
                i += 1;
                continue;
            }
            i += 1;
        }
        if let Some(value) = joined {
            values.push(value);
        }
        while i < args.len() && !args[i].starts_with('-') {
            values.push(args[i].clone());
            i += 1;
        }
    }
    values
}

/// Return the values belonging to the last occurrence of an option.
///
/// DIAMOND's parser assigns repeated vector-valued options instead of
/// appending them. This differs from [`get_all_args`], which is appropriate
/// for options such as repeated query inputs and output-field lists.
fn get_last_arg_group(args: &[String], flags: &[&str]) -> Vec<String> {
    let mut last = None;
    let mut i = 0;
    while i < args.len() {
        let mut joined = None;
        if flags.contains(&args[i].as_str()) {
            i += 1;
        } else {
            for flag in flags {
                if flag.starts_with("--") {
                    joined = args[i]
                        .strip_prefix(&format!("{flag}="))
                        .map(ToOwned::to_owned);
                } else if flag.starts_with('-') && flag.len() == 2 {
                    joined = args[i]
                        .strip_prefix(flag)
                        .filter(|value| !value.is_empty())
                        .map(ToOwned::to_owned);
                }
                if joined.is_some() {
                    break;
                }
            }
            if joined.is_none() {
                i += 1;
                continue;
            }
            i += 1;
        }

        let mut values = Vec::new();
        if let Some(value) = joined {
            values.push(value);
        }
        while i < args.len() && !args[i].starts_with('-') {
            values.push(args[i].clone());
            i += 1;
        }
        last = Some(values);
    }
    last.unwrap_or_default()
}

/// Parse a vector-valued option using DIAMOND's last-occurrence-wins rule,
/// while rejecting an occurrence without any values as the upstream parser
/// does for options whose minimum argument count is one.
fn get_last_arg_group_checked(args: &[String], flags: &[&str]) -> Result<Vec<String>, String> {
    let values = get_last_arg_group(args, flags);
    let present = args.iter().any(|arg| {
        flags.iter().any(|flag| {
            arg == flag
                || (flag.starts_with("--") && arg.starts_with(&format!("{flag}=")))
                || (flag.starts_with('-')
                    && flag.len() == 2
                    && arg.starts_with(flag)
                    && arg.len() > flag.len())
        })
    });
    if present && values.is_empty() {
        Err(format!("Option {} requires at least one value.", flags[0]))
    } else {
        Ok(values)
    }
}

fn parse_native_search_output(args: &[String]) -> Result<(Option<String>, Vec<String>), String> {
    let daa_file = get_arg(args, &["-a", "--daa"]);
    let mut output = get_arg(args, &["-o", "--out"]);
    let mut outfmt = get_all_args(args, &["-f", "--outfmt"]);

    if let Some(mut daa_file) = daa_file {
        if output.is_some() {
            return Err("Options --daa and --out cannot be used together.".to_string());
        }
        if outfmt
            .first()
            .is_some_and(|format| !format.eq_ignore_ascii_case("daa"))
        {
            return Err(
                "Invalid parameter: --daa/-a. Output file is specified with the --out/-o parameter."
                    .to_string(),
            );
        }
        if !daa_file.ends_with(".daa") {
            daa_file.push_str(".daa");
        }
        output = Some(daa_file);
        if outfmt.is_empty() {
            outfmt.push("daa".to_string());
        }
    } else if outfmt
        .first()
        .is_some_and(|format| format == "100" || format.eq_ignore_ascii_case("daa"))
    {
        let path = output
            .as_mut()
            .ok_or_else(|| "DAA output requires -a/--daa or -o/--out.".to_string())?;
        if !path.ends_with(".daa") {
            path.push_str(".daa");
        }
    }

    Ok((output, outfmt))
}

fn has_flag(args: &[String], flag: &str) -> bool {
    args.iter().any(|a| a == flag)
}

/// Whether an option occurs in separated (`--foo value`), long joined
/// (`--foo=value`), or short joined (`-b2`) form. This deliberately checks
/// presence rather than successfully parsing a value: an option unsupported by
/// the native pipeline must never be silently ignored just because it is
/// malformed.
fn has_option(args: &[String], flags: &[&str]) -> bool {
    args.iter().any(|arg| {
        flags.iter().any(|flag| {
            arg == flag
                || (flag.starts_with("--")
                    && arg
                        .strip_prefix(flag)
                        .is_some_and(|suffix| suffix.starts_with('=')))
                || (flag.starts_with('-')
                    && flag.len() == 2
                    && arg
                        .strip_prefix(flag)
                        .is_some_and(|suffix| !suffix.is_empty()))
        })
    })
}

#[derive(Debug, Clone, PartialEq)]
struct NativeViewCli {
    daa_file: String,
    output: String,
    outfmt: Vec<String>,
    max_target_seqs: i64,
    toppercent: Option<f64>,
    forward_only: bool,
    header: Vec<String>,
    threads: usize,
    qnum_offset: u64,
    snum_offset: u64,
    verbosity: u8,
}

impl Default for NativeViewCli {
    fn default() -> Self {
        Self {
            daa_file: String::new(),
            output: String::new(),
            outfmt: Vec::new(),
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            header: Vec::new(),
            threads: 0,
            qnum_offset: 0,
            snum_offset: 0,
            verbosity: 0,
        }
    }
}

/// Match upstream's `atoll`-based parsing for signed ordinal offsets, but make
/// its result deterministic across Rust targets: leading whitespace and a
/// sign are accepted, parsing stops at the first non-digit, no digits means
/// zero, and overflow saturates to the signed 64-bit limits.
fn parse_cpp_atoll(value: &str) -> i64 {
    let bytes = value.as_bytes();
    let mut index = bytes
        .iter()
        .position(|byte| !byte.is_ascii_whitespace())
        .unwrap_or(bytes.len());
    let negative = match bytes.get(index) {
        Some(b'-') => {
            index += 1;
            true
        }
        Some(b'+') => {
            index += 1;
            false
        }
        _ => false,
    };
    let limit = if negative {
        (i64::MAX as u64) + 1
    } else {
        i64::MAX as u64
    };
    let mut magnitude = 0u64;
    let mut saw_digit = false;
    while let Some(byte) = bytes.get(index) {
        if !byte.is_ascii_digit() {
            break;
        }
        saw_digit = true;
        magnitude = magnitude
            .saturating_mul(10)
            .saturating_add(u64::from(byte - b'0'))
            .min(limit);
        index += 1;
    }
    if !saw_digit {
        return 0;
    }
    if negative {
        if magnitude == (i64::MAX as u64) + 1 {
            i64::MIN
        } else {
            -(magnitude as i64)
        }
    } else {
        magnitude as i64
    }
}

/// Match `atof`, which is what upstream's option reader uses. Parsing stops at
/// the first character outside the C floating grammar and a non-number is 0.
fn parse_cpp_atof(value: &str) -> f64 {
    let trimmed = value.trim_start();
    let (negative, unsigned) = match trimmed.as_bytes().first() {
        Some(b'-') => (true, &trimmed[1..]),
        Some(b'+') => (false, &trimmed[1..]),
        _ => (false, trimmed),
    };
    let lower = unsigned.to_ascii_lowercase();
    if lower.starts_with("infinity") || lower.starts_with("inf") {
        return if negative {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }
    if lower.starts_with("nan") {
        return if negative { -f64::NAN } else { f64::NAN };
    }

    // C strtod/atof also accepts hexadecimal floating constants.
    if lower.starts_with("0x") {
        let bytes = lower.as_bytes();
        let mut i = 2usize;
        let mut mantissa = 0.0f64;
        let mut digits = 0usize;
        while let Some(digit) = bytes.get(i).and_then(|b| (*b as char).to_digit(16)) {
            mantissa = mantissa * 16.0 + f64::from(digit);
            digits += 1;
            i += 1;
        }
        if bytes.get(i) == Some(&b'.') {
            i += 1;
            let mut place = 1.0 / 16.0;
            while let Some(digit) = bytes.get(i).and_then(|b| (*b as char).to_digit(16)) {
                mantissa += f64::from(digit) * place;
                place /= 16.0;
                digits += 1;
                i += 1;
            }
        }
        if digits != 0 {
            let mut exponent = 0i64;
            if matches!(bytes.get(i), Some(b'p' | b'P')) {
                let exponent_text = &lower[i + 1..];
                let exponent_len = exponent_text
                    .as_bytes()
                    .iter()
                    .enumerate()
                    .take_while(|(j, b)| {
                        b.is_ascii_digit() || (*j == 0 && matches!(b, b'+' | b'-'))
                    })
                    .count();
                if exponent_text[..exponent_len]
                    .bytes()
                    .any(|byte| byte.is_ascii_digit())
                {
                    exponent = parse_cpp_atoll(&exponent_text[..exponent_len]);
                }
            }
            let magnitude = if exponent > i64::from(i32::MAX) {
                f64::INFINITY
            } else if exponent < i64::from(i32::MIN) {
                0.0
            } else {
                mantissa * 2.0f64.powi(exponent as i32)
            };
            return if negative { -magnitude } else { magnitude };
        }
    }

    let bytes = unsigned.as_bytes();
    let mut i = 0usize;
    let mut digits = 0usize;
    while bytes.get(i).is_some_and(u8::is_ascii_digit) {
        i += 1;
        digits += 1;
    }
    if bytes.get(i) == Some(&b'.') {
        i += 1;
        while bytes.get(i).is_some_and(u8::is_ascii_digit) {
            i += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return 0.0;
    }
    if matches!(bytes.get(i), Some(b'e' | b'E')) {
        let exponent_start = i;
        i += 1;
        if matches!(bytes.get(i), Some(b'+' | b'-')) {
            i += 1;
        }
        let exponent_digits = i;
        while bytes.get(i).is_some_and(u8::is_ascii_digit) {
            i += 1;
        }
        if i == exponent_digits {
            i = exponent_start;
        }
    }
    let prefix = &unsigned[..i];
    let magnitude = prefix.parse::<f64>().unwrap_or(0.0);
    if negative {
        -magnitude
    } else {
        magnitude
    }
}

fn is_cpp_option_start(value: &str) -> bool {
    let bytes = value.as_bytes();
    bytes.len() > 1 && bytes[0] == b'-' && !bytes[1].is_ascii_digit()
}

fn unsupported_view_option(id: &str) -> bool {
    matches!(
        id,
        "query"
            | "q"
            | "db"
            | "d"
            | "evalue"
            | "e"
            | "matrix"
            | "gapopen"
            | "gapextend"
            | "unal"
            | "compress"
            | "sam-query-len"
            | "masking"
            | "comp-based-stats"
            | "id"
            | "query-cover"
            | "subject-cover"
            | "block-size"
            | "b"
            | "index-chunks"
            | "c"
            | "sensitive"
            | "very-sensitive"
            | "ultra-sensitive"
    )
}

/// Parse the native view surface with the same grouping rules as DIAMOND's
/// CommandLineParser: a switch owns every following token until the next
/// non-numeric `-` token, scalar options require exactly one value, vector
/// options are replaced by their last occurrence, and joined long options are
/// invalid while joined short options are accepted.
fn parse_native_view_cli(args: &[String]) -> Result<NativeViewCli, String> {
    let mut parsed = NativeViewCli::default();
    let mut top_present = false;
    let mut max_target_present = false;
    let mut verbose = false;
    let mut quiet = false;
    let mut debug_log = false;
    let mut i = 2usize;

    while i < args.len() {
        let option = &args[i];
        if !is_cpp_option_start(option) {
            return Err("Command line options must begin with - or --.".to_string());
        }
        let (id, attached, short) = if let Some(id) = option.strip_prefix("--") {
            (id, None, false)
        } else {
            let id_end = option
                .char_indices()
                .nth(2)
                .map_or(option.len(), |(n, _)| n);
            (
                &option[1..id_end],
                (option.len() > id_end).then(|| &option[id_end..]),
                true,
            )
        };

        let kind = match id {
            "threads" | "p" | "tmpdir" | "t" | "out" | "o" | "max-target-seqs" | "k" | "top"
            | "daa" | "a" | "qnum-offset" | "snum-offset" => 1,
            "outfmt" | "f" => 2,
            "header" => 3,
            "verbose" | "v" | "log" | "quiet" | "forwardonly" => 0,
            _ if unsupported_view_option(id) => {
                return Err(format!("Option is not permitted for this workflow: {id}"));
            }
            _ => return Err(format!("Invalid option: {id}")),
        };

        let mut values = Vec::new();
        if let Some(value) = attached {
            values.push(value.to_string());
        }
        i += 1;
        while i < args.len() && !is_cpp_option_start(&args[i]) {
            values.push(args[i].clone());
            i += 1;
        }

        let display = match id {
            "p" | "threads" => "-p/--threads",
            "t" | "tmpdir" => "-t/--tmpdir",
            "o" | "out" => "-o/--out",
            "k" | "max-target-seqs" => "-k/--max-target-seqs",
            "a" | "daa" => "-a/--daa",
            "f" | "outfmt" => "-f/--outfmt",
            _ if short => option.as_str(),
            _ => option.as_str(),
        };
        let valid_count = match kind {
            0 => values.is_empty(),
            1 => values.len() == 1,
            2 => !values.is_empty(),
            3 => true,
            _ => unreachable!(),
        };
        if !valid_count {
            return Err(format!("Invalid parameter count for option '{display}'"));
        }

        match id {
            "threads" | "p" => {
                parsed.threads = parse_cpp_atoll(&values[0]).max(0) as usize;
            }
            "out" | "o" => parsed.output = values[0].clone(),
            "max-target-seqs" | "k" => {
                max_target_present = true;
                parsed.max_target_seqs = parse_cpp_atoll(&values[0]);
                if parsed.max_target_seqs == 0 {
                    parsed.max_target_seqs = i64::MAX;
                }
            }
            "top" => {
                top_present = true;
                parsed.toppercent = Some(parse_cpp_atof(&values[0]));
            }
            "daa" | "a" => parsed.daa_file = values[0].clone(),
            "outfmt" | "f" => parsed.outfmt = values,
            "header" => {
                // Preserve presence: for view/tabular, a bare --header means
                // the verbose header, whereas absence means no header.
                parsed.header = if values.is_empty() {
                    vec![String::new()]
                } else {
                    values
                };
            }
            "qnum-offset" => parsed.qnum_offset = parse_cpp_atoll(&values[0]) as u64,
            "snum-offset" => parsed.snum_offset = parse_cpp_atoll(&values[0]) as u64,
            "forwardonly" => parsed.forward_only = true,
            "verbose" | "v" => verbose = true,
            "quiet" => quiet = true,
            "log" => debug_log = true,
            "tmpdir" | "t" => {}
            _ => unreachable!(),
        }
    }

    if top_present && max_target_present {
        return Err("--top and -k/--max-target-seqs are mutually exclusive.".to_string());
    }
    if let Some(format) = parsed.outfmt.first() {
        if format == "102" {
            return Err("Taxonomy features are not supported for the DAA format.".to_string());
        }
        if !matches!(
            format.as_str(),
            "0" | "5"
                | "xml"
                | "6"
                | "tab"
                | "100"
                | "daa"
                | "101"
                | "sam"
                | "103"
                | "paf"
                | "104"
                | "json-flat"
                | "null"
                | "edge"
        ) {
            return Err(format!(
                "Invalid output format: {format}\nAllowed values: 0,5,xml,6,tab,100,daa,101,sam,102,103,104,paf"
            ));
        }
    }
    if let Some(top) = parsed.toppercent {
        if top < 0.0 || top > 100.0 {
            return Err("Allowed value range for --top is between 0.0 and 100.0".to_string());
        }
        if top == 100.0 {
            parsed.toppercent = None;
            parsed.max_target_seqs = i64::MAX;
        }
    }
    parsed.verbosity = if debug_log {
        3
    } else if quiet {
        0
    } else if verbose {
        2
    } else if parsed.output.is_empty() && args.len() != 2 {
        0
    } else {
        1
    };
    Ok(parsed)
}

fn validate_native_cluster_options(args: &[String]) -> Result<(), String> {
    if has_option(args, &["--parallel-tmpdir"]) {
        if args.get(1).is_some_and(|command| command == "linclust") {
            return Err(
                "Native linclust --parallel-tmpdir is not translated; use --legacy.".to_string(),
            );
        }
        return Err("Option is not permitted for this workflow: --parallel-tmpdir".to_string());
    }
    if has_option(args, &["--aln-out"]) {
        return Err(
            "Alignment output to file is not supported without --parallel-tmpdir.".to_string(),
        );
    }
    Ok(())
}

fn parse_optional_cluster_number<T: std::str::FromStr>(
    args: &[String],
    flags: &[&str],
) -> Result<Option<T>, String> {
    match get_arg(args, flags) {
        Some(value) => value
            .parse()
            .map(Some)
            .map_err(|_| format!("Invalid value for {}: {value}", flags.last().unwrap())),
        None if has_option(args, flags) => Err(format!("Option {} requires a value.", flags[0])),
        None => Ok(None),
    }
}

fn parse_optional_cluster_string(
    args: &[String],
    flags: &[&str],
) -> Result<Option<String>, String> {
    match get_arg(args, flags) {
        Some(value) => Ok(Some(value)),
        None if has_option(args, flags) => Err(format!("Option {} requires a value.", flags[0])),
        None => Ok(None),
    }
}

fn parse_native_cluster_config(args: &[String]) -> Result<ClusterConfig, String> {
    let workflow = match args.get(1).map(String::as_str) {
        Some("linclust") => ClusterWorkflow::LinClust,
        Some("deepclust") => ClusterWorkflow::DeepClust,
        _ => ClusterWorkflow::Cascaded,
    };
    let memory_limit = match get_arg(args, &["-M", "--memory-limit"]) {
        Some(value) => parse_byte_size(&value)
            .map_err(|error| format!("invalid --memory-limit '{value}': {error}"))?,
        None if has_option(args, &["-M", "--memory-limit"]) => {
            return Err("Option --memory-limit requires a value.".to_string());
        }
        None => DEFAULT_NATIVE_MEMORY_LIMIT,
    };
    let graph_algo = match get_arg(args, &["--graph-algo"]).as_deref() {
        None => GraphAlgo::GreedyVertexCover,
        Some("gvc") => GraphAlgo::GreedyVertexCover,
        Some("len") => GraphAlgo::LengthSorted,
        Some(value) => return Err(format!("Invalid value for --graph-algo: {value}")),
    };
    let header = get_last_arg_group(args, &["--header"]);
    let simple_header = if !has_option(args, &["--header"]) {
        false
    } else if header.is_empty() {
        true
    } else if header.len() > 1 {
        return Err(format!("Invalid header format: {}", header.join(" ")));
    } else {
        match header[0].as_str() {
            "0" => false,
            "simple" => true,
            "verbose" => {
                return Err(
                    "Verbose header format is not supported for cluster workflow.".to_string(),
                );
            }
            value => return Err(format!("Invalid header format: {value}")),
        }
    };
    let matrix = match get_arg(args, &["--matrix"]) {
        Some(value) => value,
        None if has_option(args, &["--matrix"]) => {
            return Err("Option --matrix requires a value.".to_string());
        }
        None => "blosum62".to_string(),
    };
    let masking = diamond::masking::MaskingMode::parse(
        &parse_optional_cluster_string(args, &["--masking"])?.unwrap_or_else(|| "0".to_string()),
    )?;
    let soft_masking_value = parse_optional_cluster_string(args, &["--soft-masking"])?
        .unwrap_or_else(|| "tantan".to_string());
    let soft_masking = match soft_masking_value.as_str() {
        "0" | "none" => diamond::masking::MaskingAlgo::None,
        "tantan" => diamond::masking::MaskingAlgo::Tantan,
        "seg" => diamond::masking::MaskingAlgo::Seg,
        value => {
            return Err(format!(
                "Invalid value for string field: {value}. Permitted values: 0, none, seg, tantan"
            ));
        }
    };
    let comp_based_stats =
        parse_optional_cluster_number(args, &["--comp-based-stats"])?.unwrap_or(1);
    if !(0..=5).contains(&comp_based_stats) {
        return Err(
            "Invalid value for --comp-based-stats. Permitted values: 0, 1, 2, 3, 4, 5.".to_string(),
        );
    }
    let extension_mode = parse_optional_cluster_string(args, &["--ext"])?.unwrap_or_default();
    if !extension_mode.is_empty() {
        diamond::search::sensitivity::ExtensionMode::from_string(&extension_mode)?;
    }
    let motif_masking =
        parse_optional_cluster_string(args, &["--motif-masking"])?.unwrap_or_default();
    if !motif_masking.is_empty() && motif_masking != "0" && motif_masking != "1" {
        return Err("Permitted values for --motif-masking: 0, 1".to_string());
    }
    Ok(ClusterConfig {
        workflow,
        database: get_arg(args, &["-d", "--db"]).unwrap_or_default(),
        output: get_arg(args, &["-o", "--out"]).unwrap_or_default(),
        threads: parse_optional_cluster_number(args, &["-p", "--threads"])?.unwrap_or(0),
        matrix,
        gap_open: parse_optional_cluster_number(args, &["--gapopen"])?.unwrap_or(-1),
        gap_extend: parse_optional_cluster_number(args, &["--gapextend"])?.unwrap_or(-1),
        max_evalue: parse_optional_cluster_number(args, &["-e", "--evalue"])?.unwrap_or(0.001),
        comp_based_stats,
        masking,
        soft_masking,
        motif_masking,
        extension_mode,
        member_cover: parse_optional_cluster_number(args, &["--member-cover"])?,
        mutual_cover: parse_optional_cluster_number(args, &["--mutual-cover"])?,
        approx_id: parse_optional_cluster_number(args, &["--approx-id"])?,
        cluster_steps: get_last_arg_group_checked(args, &["--cluster-steps"])?,
        connected_component_depth: get_last_arg_group_checked(
            args,
            &["--connected-component-depth"],
        )?,
        round_coverage: get_last_arg_group_checked(args, &["--round-coverage"])?,
        round_approx_id: get_last_arg_group_checked(args, &["--round-approx-id"])?,
        kmer_ranking: has_flag(args, "--kmer-ranking"),
        no_reassign: has_flag(args, "--no-reassign"),
        weighted_gvc: has_flag(args, "--weighted-gvc"),
        strict_gvc: has_flag(args, "--strict-gvc"),
        graph_algo,
        memory_limit,
        oid_output: has_flag(args, "--oid-output"),
        simple_header,
        alignment_output: get_arg(args, &["--aln-out"]),
    })
}

/// Manual native-search parsing must not turn an unrecognised option into a
/// no-op. Unknown switches are sent to the compatibility parser, which either
/// handles them through the optional C++ backend or rejects them explicitly.
fn has_unknown_search_option(args: &[String], blastx: bool) -> bool {
    const COMMON: &[&str] = &[
        "--query",
        "--db",
        "--out",
        "--daa",
        "--matrix",
        "--gapopen",
        "--gapextend",
        "--evalue",
        "--max-target-seqs",
        "--ext-chunk-size",
        "--top",
        "--global-ranking",
        "--id",
        "--threads",
        "--outfmt",
        "--faster",
        "--fast",
        "--mid-sensitive",
        "--sensitive",
        "--more-sensitive",
        "--very-sensitive",
        "--ultra-sensitive",
        "--masking",
        "--motif-masking",
        "--min-query-len",
        "--query-cover",
        "--subject-cover",
        "--comp-based-stats",
        "--no-self-hits",
        "--xdrop",
        "--memory-limit",
        "--target-indexed",
        "--tmpdir",
        "--legacy",
        "--help",
        // Recognised incompatibilities are included here because the explicit
        // routing checks below decide whether a value is natively supported.
        "--max-hsps",
        "--block-size",
        "--index-chunks",
        "--algo",
        "--soft-masking",
        "--approx-id",
        "--file-buffer-size",
        "--query-parallel-limit",
        "--min-score",
        "--header",
        "--compress",
        "--iterate",
    ];
    const BLASTX: &[&str] = &[
        "--query-gencode",
        "--strand",
        "--min-orf",
        "--frameshift",
        "--swipe",
    ];
    const SHORT_COMMON: &[&str] = &[
        "-q", "-d", "-o", "-a", "-e", "-k", "-p", "-f", "-x", "-t", "-b", "-c", "-h",
    ];

    args.iter().skip(2).any(|arg| {
        if let Some(long) = arg.strip_prefix("--") {
            let name = &arg[..2 + long.find('=').unwrap_or(long.len())];
            !COMMON.contains(&name) && !(blastx && BLASTX.contains(&name))
        } else if arg.starts_with('-') && arg.as_bytes().get(1).is_some_and(u8::is_ascii_alphabetic)
        {
            let name = &arg[..2];
            (name == "-h" && arg.len() != 2)
                || (!SHORT_COMMON.contains(&name) && !(blastx && name == "-F"))
        } else {
            false
        }
    })
}

fn route_common_search_to_legacy(args: &[String], query_translated: bool) -> bool {
    // The native reader decodes NCBI volumes correctly, but the pinned C++
    // search path corrupts the first appended raw chunk through undefined
    // StringSet perimeter handling. Until that compatibility policy is
    // resolved, keep the whole command on the compatibility path instead of
    // presenting the intended Rust result as byte parity with upstream.
    let ncbi_blast_database = get_arg(args, &["-d", "--db"]).is_some_and(|database| {
        let path = std::path::Path::new(&database);
        let mut pin = path.as_os_str().to_owned();
        pin.push(".pin");
        let mut pal = path.as_os_str().to_owned();
        pal.push(".pal");
        path.extension()
            .is_some_and(|extension| extension == "pin" || extension == "pal")
            || std::path::PathBuf::from(pin).exists()
            || std::path::PathBuf::from(pal).exists()
    });
    let requested_output = get_all_args(args, &["-f", "--outfmt"]);
    let unsupported_output = requested_output.first().is_some_and(|format| {
        format != "6" && format != "tab" && format != "100" && format != "daa"
    }) || (!requested_output
        .first()
        .is_some_and(|format| format == "100" || format == "daa")
        && requested_output
            .iter()
            .skip(1)
            .any(|name| match FieldId::from_name(name) {
                Some(
                    FieldId::QSeqId
                    | FieldId::SSeqId
                    | FieldId::PIdent
                    | FieldId::Length
                    | FieldId::Mismatch
                    | FieldId::GapOpen
                    | FieldId::QStart
                    | FieldId::QEnd
                    | FieldId::SStart
                    | FieldId::SEnd
                    | FieldId::EValue
                    | FieldId::BitScore
                    | FieldId::Score
                    | FieldId::NIdent
                    | FieldId::Positive
                    | FieldId::Gaps
                    | FieldId::PPos
                    | FieldId::QLen
                    | FieldId::SLen
                    | FieldId::QFrame
                    | FieldId::QCovHsp
                    | FieldId::SCovHsp
                    | FieldId::STaxIds
                    | FieldId::SSciNames
                    | FieldId::SSKingdoms
                    | FieldId::SKingdoms
                    | FieldId::SPhylums
                    | FieldId::SLineages,
                )
                | None => false,
                Some(FieldId::FullQSeqMate) if query_translated => false,
                Some(_) => true,
            }));

    let target_indexed = has_flag(args, "--target-indexed");

    ncbi_blast_database
        || unsupported_output
        || get_arg(args, &["--masking"]).is_some_and(|m| m.eq_ignore_ascii_case("seg"))
        || get_arg(args, &["--max-hsps"]).is_some_and(|m| m != "1")
        || get_arg(args, &["--comp-based-stats"]).is_some_and(|m| m != "0" && m != "1")
        || get_arg(args, &["--global-ranking"]).is_some_and(|m| m != "0")
        || has_option(
            args,
            &[
                "-b",
                "--block-size",
                "--algo",
                "--soft-masking",
                "--approx-id",
                "--file-buffer-size",
                "--query-parallel-limit",
                "--min-score",
                "--header",
                "--compress",
                "--iterate",
            ],
        )
        || (has_option(args, &["-c", "--index-chunks"]) && !target_indexed)
}

fn route_blastp_to_legacy(args: &[String]) -> bool {
    has_unknown_search_option(args, false) || route_common_search_to_legacy(args, false)
}

fn route_blastx_to_legacy(args: &[String]) -> bool {
    has_unknown_search_option(args, true)
        || route_common_search_to_legacy(args, true)
        || get_arg(args, &["-F", "--frameshift"]).is_some_and(|f| f != "0")
        || has_flag(args, "--swipe")
        || has_option(args, &["--min-query-len"])
}

fn parse_arg_or<T: std::str::FromStr>(args: &[String], flags: &[&str], default: T) -> T {
    get_arg(args, flags)
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn parse_max_target_seqs(args: &[String]) -> i64 {
    match parse_arg_or(args, &["-k", "--max-target-seqs"], 25_i64) {
        // C++ normalizes the documented unlimited sentinel before search or
        // view code reaches `Config::output_range`.
        0 => i64::MAX,
        value => value,
    }
}

fn native_getseq_option_error(args: &[String]) -> Option<String> {
    args.iter().skip(2).find_map(|arg| {
        let error = match arg.as_str() {
            // Hidden and disabled in standard upstream builds.
            "--reverse" | "--hard-masked" | "--oid-list" | "--chunk-size" => {
                format!("Invalid option: {}", &arg[2..])
            }
            // Registered for other workflows, but not getseq.
            "--query" | "--block-size" => {
                format!("Option is not permitted for this workflow: {}", &arg[2..])
            }
            "-q" | "-b" => {
                format!("Option is not permitted for this workflow: {}", &arg[1..])
            }
            value if value.starts_with("-q") || value.starts_with("-b") => {
                format!(
                    "Option is not permitted for this workflow: {}",
                    &value[1..2]
                )
            }
            value
                if value.starts_with("--reverse=")
                    || value.starts_with("--hard-masked=")
                    || value.starts_with("--oid-list=")
                    || value.starts_with("--chunk-size=")
                    || value.starts_with("--query=")
                    || value.starts_with("--block-size=") =>
            {
                format!("Invalid option: {}", &value[2..])
            }
            _ => return None,
        };
        Some(error)
    })
}

fn parse_target_indexed(args: &[String]) -> bool {
    if !has_flag(args, "--target-indexed") {
        return false;
    }
    if get_arg(args, &["-c", "--index-chunks"]).as_deref() != Some("1") {
        print_banner();
        eprintln!("Error: --target-indexed requires -c1.");
        std::process::exit(1);
    }
    true
}

const DEFAULT_NATIVE_MEMORY_LIMIT: u64 = 16_000_000_000;

fn parse_memory_limit_arg(args: &[String]) -> Option<u64> {
    Some(
        get_arg(args, &["-M", "--memory-limit"]).map_or(DEFAULT_NATIVE_MEMORY_LIMIT, |value| {
            parse_byte_size(&value).unwrap_or_else(|error| {
                eprintln!("Error: invalid --memory-limit '{value}': {error}");
                std::process::exit(1);
            })
        }),
    )
}

fn parse_byte_size(value: &str) -> Result<u64, &'static str> {
    let bytes = diamond::util::string::interpret_number(value)
        .map_err(|_| "use K, M, G, or T (for example 4G)")?;
    u64::try_from(bytes).map_err(|_| "value is out of range")
}

fn parse_sensitivity(args: &[String]) -> Sensitivity {
    if has_flag(args, "--ultra-sensitive") {
        Sensitivity::UltraSensitive
    } else if has_flag(args, "--very-sensitive") {
        Sensitivity::VerySensitive
    } else if has_flag(args, "--more-sensitive") {
        Sensitivity::MoreSensitive
    } else if has_flag(args, "--sensitive") {
        Sensitivity::Sensitive
    } else if has_flag(args, "--mid-sensitive") {
        Sensitivity::MidSensitive
    } else if has_flag(args, "--fast") {
        Sensitivity::Fast
    } else if has_flag(args, "--faster") {
        Sensitivity::Faster
    } else {
        Sensitivity::Default
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_all_args_collects_outfmt_field_list() {
        let args = vec![
            "diamond".to_string(),
            "view".to_string(),
            "-f".to_string(),
            "6".to_string(),
            "qseqid".to_string(),
            "sseqid".to_string(),
            "--top".to_string(),
            "10".to_string(),
            "--outfmt".to_string(),
            "104".to_string(),
            "score".to_string(),
        ];

        assert_eq!(
            get_all_args(&args, &["-f", "--outfmt"]),
            ["6", "qseqid", "sseqid", "104", "score"]
        );
    }

    #[test]
    fn native_getseq_rejects_ignored_options_like_standard_upstream() {
        assert_eq!(
            native_getseq_option_error(&args(&["diamond", "getseq", "--hard-masked", "--reverse"])),
            Some("Invalid option: hard-masked".to_string())
        );
        assert_eq!(
            native_getseq_option_error(&args(&["diamond", "getseq", "--query", "ids.tsv"])),
            Some("Option is not permitted for this workflow: query".to_string())
        );
        assert_eq!(
            native_getseq_option_error(&args(&["diamond", "getseq", "-b2"])),
            Some("Option is not permitted for this workflow: b".to_string())
        );
        assert_eq!(
            native_getseq_option_error(&args(&["diamond", "getseq", "--oid-list=ids.tsv"])),
            Some("Invalid option: oid-list=ids.tsv".to_string())
        );
        assert_eq!(
            native_getseq_option_error(&args(&["diamond", "getseq", "--seq", "1"])),
            None
        );
    }

    #[test]
    fn test_get_all_args_accepts_joined_first_value() {
        assert_eq!(
            get_all_args(
                &args(&["diamond", "blastp", "--outfmt=6", "qseqid", "sseqid"]),
                &["-f", "--outfmt"]
            ),
            ["6", "qseqid", "sseqid"]
        );
        assert_eq!(
            get_all_args(
                &args(&["diamond", "blastp", "-f6", "qseqid"]),
                &["-f", "--outfmt"]
            ),
            ["6", "qseqid"]
        );
    }

    #[test]
    fn last_arg_group_replaces_repeated_occurrences() {
        assert_eq!(
            get_last_arg_group(
                &args(&[
                    "diamond",
                    "makedb",
                    "--in",
                    "first.faa",
                    "--in",
                    "second.faa"
                ]),
                &["--in"]
            ),
            ["second.faa"]
        );
        assert_eq!(
            get_last_arg_group(
                &args(&["diamond", "makedb", "--in", "one.faa", "two.faa"]),
                &["--in"]
            ),
            ["one.faa", "two.faa"]
        );
        assert_eq!(
            get_last_arg_group(
                &args(&["diamond", "makedb", "--in=first.faa", "--in=last.faa"]),
                &["--in"]
            ),
            ["last.faa"]
        );
    }

    #[test]
    fn checked_last_arg_group_rejects_missing_values() {
        assert_eq!(
            get_last_arg_group_checked(
                &args(&["diamond", "cluster", "--cluster-steps"]),
                &["--cluster-steps"]
            ),
            Err("Option --cluster-steps requires at least one value.".to_string())
        );
        assert_eq!(
            get_last_arg_group_checked(
                &args(&[
                    "diamond",
                    "cluster",
                    "--cluster-steps",
                    "fast",
                    "--cluster-steps",
                    "default"
                ]),
                &["--cluster-steps"]
            ),
            Ok(vec!["default".to_string()])
        );
    }

    #[test]
    fn native_cluster_rejects_untranslated_alignment_output_routes() {
        assert_eq!(
            validate_native_cluster_options(&args(&[
                "diamond",
                "cluster",
                "--aln-out",
                "edges.tsv"
            ])),
            Err("Alignment output to file is not supported without --parallel-tmpdir.".to_string())
        );
        assert_eq!(
            validate_native_cluster_options(&args(&[
                "diamond",
                "cluster",
                "--parallel-tmpdir=/tmp/work"
            ])),
            Err("Option is not permitted for this workflow: --parallel-tmpdir".to_string())
        );
        assert_eq!(
            validate_native_cluster_options(&args(&[
                "diamond",
                "linclust",
                "--parallel-tmpdir",
                "/tmp/work"
            ])),
            Err("Native linclust --parallel-tmpdir is not translated; use --legacy.".to_string())
        );
    }

    #[test]
    fn native_cluster_parser_forwards_vector_threshold_and_memory_controls() {
        let parsed = parse_native_cluster_config(&args(&[
            "diamond",
            "cluster",
            "-d",
            "input.dmnd",
            "-o",
            "clusters.tsv",
            "--member-cover",
            "75",
            "--mutual-cover=80",
            "--cluster-steps",
            "fast",
            "default",
            "--cluster-steps=faster_lin",
            "sensitive",
            "--connected-component-depth",
            "2",
            "--round-coverage",
            "81",
            "82",
            "--round-approx-id",
            "55",
            "60",
            "--kmer-ranking",
            "--no-reassign",
            "--weighted-gvc",
            "--strict-gvc",
            "--graph-algo=len",
            "--matrix=pam250",
            "--gapopen",
            "15",
            "--gapextend=3",
            "--evalue=1e-20",
            "--comp-based-stats=4",
            "--masking=tantan",
            "--soft-masking=seg",
            "--motif-masking=0",
            "--ext=banded-slow",
            "--oid-output",
            "--header=simple",
            "-M3G",
        ]))
        .unwrap();
        assert_eq!(parsed.member_cover, Some(75.0));
        assert_eq!(parsed.mutual_cover, Some(80.0));
        assert_eq!(parsed.cluster_steps, ["faster_lin", "sensitive"]);
        assert_eq!(parsed.connected_component_depth, ["2"]);
        assert_eq!(parsed.round_coverage, ["81", "82"]);
        assert_eq!(parsed.round_approx_id, ["55", "60"]);
        assert!(parsed.kmer_ranking && parsed.no_reassign);
        assert!(parsed.weighted_gvc && parsed.strict_gvc);
        assert_eq!(parsed.graph_algo, GraphAlgo::LengthSorted);
        assert_eq!(parsed.memory_limit, 3_000_000_000);
        assert_eq!(parsed.matrix, "pam250");
        assert_eq!(parsed.gap_open, 15);
        assert_eq!(parsed.gap_extend, 3);
        assert_eq!(parsed.max_evalue, 1e-20);
        assert_eq!(parsed.comp_based_stats, 4);
        assert_eq!(parsed.masking, diamond::masking::MaskingMode::Tantan);
        assert_eq!(parsed.soft_masking, diamond::masking::MaskingAlgo::Seg);
        assert_eq!(parsed.motif_masking, "0");
        assert_eq!(parsed.extension_mode, "banded-slow");
        assert!(parsed.oid_output);
        assert!(parsed.simple_header);
    }

    #[test]
    fn native_cluster_search_control_defaults_match_init_thresholds() {
        let parsed = parse_native_cluster_config(&args(&["diamond", "cluster"])).unwrap();
        assert_eq!(parsed.max_evalue, 0.001);
        assert_eq!(parsed.comp_based_stats, 1);
        assert_eq!(parsed.masking, diamond::masking::MaskingMode::None);
        assert_eq!(parsed.soft_masking, diamond::masking::MaskingAlgo::Tantan);
        assert!(parsed.motif_masking.is_empty());
        assert!(parsed.extension_mode.is_empty());
        assert!(parse_native_cluster_config(&args(&[
            "diamond",
            "cluster",
            "--comp-based-stats=6"
        ]))
        .is_err());
        assert!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--ext=invalid"])).is_err()
        );
    }

    #[test]
    fn native_cluster_header_matches_upstream_forms_and_errors() {
        assert!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--header"]))
                .unwrap()
                .simple_header
        );
        assert!(
            !parse_native_cluster_config(&args(&["diamond", "cluster", "--header", "0"]))
                .unwrap()
                .simple_header
        );
        assert_eq!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--header", "verbose"]))
                .unwrap_err(),
            "Verbose header format is not supported for cluster workflow."
        );
    }

    #[test]
    fn native_cluster_scoring_controls_require_values() {
        assert_eq!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--matrix"])).unwrap_err(),
            "Option --matrix requires a value."
        );
        assert_eq!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--gapopen", "bad"]))
                .unwrap_err(),
            "Invalid value for --gapopen: bad"
        );
        assert_eq!(
            parse_native_cluster_config(&args(&["diamond", "cluster", "--gapextend"])).unwrap_err(),
            "Option --gapextend requires a value."
        );
    }

    #[test]
    fn native_search_output_parses_daa_aliases_and_suffixes() {
        assert_eq!(
            parse_native_search_output(&args(&["diamond", "blastp", "-a", "hits"])),
            Ok((Some("hits.daa".to_string()), vec!["daa".to_string()]))
        );
        assert_eq!(
            parse_native_search_output(&args(&[
                "diamond", "blastx", "--outfmt", "100", "--out", "hits"
            ])),
            Ok((Some("hits.daa".to_string()), vec!["100".to_string()]))
        );
        assert_eq!(
            parse_native_search_output(&args(&[
                "diamond",
                "blastp",
                "--daa=hits.daa",
                "--outfmt",
                "daa"
            ])),
            Ok((Some("hits.daa".to_string()), vec!["daa".to_string()]))
        );
    }

    #[test]
    fn native_search_output_rejects_conflicting_or_invalid_daa_options() {
        assert!(parse_native_search_output(&args(&[
            "diamond", "blastp", "-a", "hits", "-o", "other"
        ]))
        .is_err());
        assert!(parse_native_search_output(&args(&[
            "diamond", "blastp", "-a", "hits", "--outfmt", "6"
        ]))
        .is_err());
        assert!(
            parse_native_search_output(&args(&["diamond", "blastp", "--outfmt", "100"])).is_err()
        );
    }

    #[test]
    fn parse_memory_sizes() {
        assert_eq!(parse_byte_size("512M"), Ok(512_000_000));
        assert_eq!(parse_byte_size("1.5G"), Ok(1_500_000_000));
        assert!(parse_byte_size("lots").is_err());
        assert!(parse_byte_size("512").is_err());

        let args = vec!["diamond".to_string(), "blastp".to_string()];
        assert_eq!(parse_memory_limit_arg(&args), Some(16_000_000_000));
        let args = vec![
            "diamond".to_string(),
            "blastp".to_string(),
            "--memory-limit".to_string(),
            "0G".to_string(),
        ];
        assert_eq!(parse_memory_limit_arg(&args), Some(0));
    }

    #[test]
    fn zero_max_target_sequences_means_unlimited() {
        assert_eq!(parse_max_target_seqs(&args(&["diamond", "view"])), 25);
        assert_eq!(
            parse_max_target_seqs(&args(&["diamond", "view", "-k", "0"])),
            i64::MAX
        );
        assert_eq!(
            parse_max_target_seqs(&args(&["diamond", "blastp", "--max-target-seqs=7"])),
            7
        );
    }

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_string()).collect()
    }

    #[test]
    fn native_view_rejects_search_only_options_like_upstream() {
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--unal", "0"])).unwrap_err(),
            "Option is not permitted for this workflow: unal"
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--compress", "1"])).unwrap_err(),
            "Option is not permitted for this workflow: compress"
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--wat"])).unwrap_err(),
            "Invalid option: wat"
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--top=10"])).unwrap_err(),
            "Invalid option: top=10"
        );
    }

    #[test]
    fn native_view_offsets_match_upstream_atoll_and_unsigned_conversion() {
        assert_eq!(parse_cpp_atoll("nope"), 0);
        assert_eq!(parse_cpp_atoll("  +12tail"), 12);
        assert_eq!(parse_cpp_atoll("9223372036854775808"), i64::MAX);
        assert_eq!(parse_cpp_atoll("-9223372036854775809"), i64::MIN);
        assert_eq!(
            parse_native_view_cli(&args(&[
                "diamond",
                "view",
                "--qnum-offset",
                "-1",
                "--snum-offset",
                "20",
            ]))
            .map(|config| (config.qnum_offset, config.snum_offset)),
            Ok((u64::MAX, 20))
        );
        // Like the C++ option parser, the last repeated occurrence wins.
        assert_eq!(
            parse_native_view_cli(&args(&[
                "diamond",
                "view",
                "--qnum-offset",
                "3",
                "--qnum-offset",
                "-2",
            ]))
            .map(|config| (config.qnum_offset, config.snum_offset)),
            Ok((u64::MAX - 1, 0))
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--qnum-offset"])).map(|_| ()),
            Err("Invalid parameter count for option '--qnum-offset'".to_string())
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--qnum-offset=-1",])).map(|_| ()),
            Err("Invalid option: qnum-offset=-1".to_string())
        );
    }

    #[test]
    fn native_view_top_and_vector_options_follow_cpp_rules() {
        let top_zero =
            parse_native_view_cli(&args(&["diamond", "view", "--top", "nonsense"])).unwrap();
        assert_eq!(top_zero.toppercent, Some(0.0));

        let unlimited = parse_native_view_cli(&args(&["diamond", "view", "--top", "100"])).unwrap();
        assert_eq!(unlimited.toppercent, None);
        assert_eq!(unlimited.max_target_seqs, i64::MAX);

        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--top", "-1"])).unwrap_err(),
            "Allowed value range for --top is between 0.0 and 100.0"
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "--top", "10", "-k5"])).unwrap_err(),
            "--top and -k/--max-target-seqs are mutually exclusive."
        );

        let repeated = parse_native_view_cli(&args(&[
            "diamond", "view", "-f", "6", "qseqid", "-f104", "score", "--header",
        ]))
        .unwrap();
        assert_eq!(repeated.outfmt, ["104", "score"]);
        assert_eq!(repeated.header, [""]);

        assert!(
            parse_native_view_cli(&args(&["diamond", "view", "-f", "wat"]))
                .unwrap_err()
                .starts_with("Invalid output format: wat\nAllowed values:")
        );
        assert_eq!(
            parse_native_view_cli(&args(&["diamond", "view", "-f102"])).unwrap_err(),
            "Taxonomy features are not supported for the DAA format."
        );
    }

    #[test]
    fn native_route_never_silently_ignores_unsupported_options() {
        for option in [
            ["--block-size", "1"].as_slice(),
            ["--block-size=1"].as_slice(),
            ["-b1"].as_slice(),
            ["--index-chunks", "1"].as_slice(),
            ["-c1"].as_slice(),
            ["--algo", "0"].as_slice(),
            ["--soft-masking", "0"].as_slice(),
            ["--approx-id", "80"].as_slice(),
            ["--file-buffer-size", "1024"].as_slice(),
            ["--query-parallel-limit", "2"].as_slice(),
            ["--min-score", "50"].as_slice(),
            ["--header", "1"].as_slice(),
            ["--compress", "1"].as_slice(),
            ["--iterate"].as_slice(),
        ] {
            let mut command = args(&["diamond", "blastp"]);
            command.extend(option.iter().map(|value| (*value).to_string()));
            assert!(route_blastp_to_legacy(&command), "did not route {option:?}");
        }

        assert!(route_blastp_to_legacy(&args(&[
            "diamond", "blastp", "--outfmt", "0"
        ])));
        assert!(route_blastp_to_legacy(&args(&[
            "diamond", "blastp", "--outfmt", "101"
        ])));
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond", "blastp", "--outfmt", "6", "qseqid", "sseqid"
        ])));
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--outfmt=6",
            "qlen",
            "nident",
            "score"
        ])));
        assert!(route_blastp_to_legacy(&args(&[
            "diamond", "blastp", "--outfmt", "6", "qseq"
        ])));
        // Native DMND databases carry the same mapping, nodes, ranks, and
        // scientific-name payloads as upstream DatabaseFile.
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--outfmt",
            "6",
            "staxids",
            "sscinames",
            "sskingdoms",
            "skingdoms",
            "sphylums",
            "slineages"
        ])));
        assert!(route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--outfmt",
            "6",
            "full_qseq_mate"
        ])));
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--outfmt",
            "6",
            "not_a_field"
        ])));
        // Pinned DIAMOND does not register the BLAST accession aliases in its
        // tabular field map and rejects them as invalid output fields.
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond", "blastp", "--outfmt", "6", "qacc"
        ])));
        assert!(route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--definitely-not-a-diamond-option"
        ])));
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "-qquery.faa"
        ])));
        assert!(!route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--target-indexed",
            "-c1"
        ])));
    }

    #[test]
    fn blastx_specific_unsupported_modes_route_to_legacy() {
        assert!(route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "--frameshift",
            "15"
        ])));
        assert!(route_blastx_to_legacy(&args(&[
            "diamond", "blastx", "--swipe"
        ])));
        assert!(route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "--min-query-len",
            "30"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "--query-gencode",
            "2",
            "--strand",
            "plus",
            "--frameshift",
            "0",
            "--outfmt",
            "6"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond", "blastx", "--outfmt", "6", "qseqid", "qlen", "sseqid"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "--outfmt",
            "6",
            "qseqid",
            "full_qseq_mate",
            "sseqid"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "--no-self-hits"
        ])));
        assert!(!route_blastx_to_legacy(&args(&[
            "diamond",
            "blastx",
            "-q",
            "query.fna",
            "-d",
            "db.faa",
            "--no-self-hits",
            "-o",
            "out.tsv"
        ])));
    }
}
