use std::env;
use std::path::PathBuf;

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

    // Upstream's command-line parser rejects `--unal` for the view workflow
    // based on option presence, even when its value is zero. DAA archives do
    // not contain records for queries with no reported alignments, so native
    // view cannot implement this option faithfully either.
    if args[1] == "view" && !has_flag(&args, "--legacy") && has_option(&args, &["--unal"]) {
        print_banner();
        eprintln!("Error: Option is not permitted for this workflow: unal");
        std::process::exit(1);
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
            let input_files = get_all_args(&args, &["--in"]);
            if let Some(db) = get_arg(&args, &["-d", "--db"]) {
                if !input_files.is_empty() {
                    print_banner();
                    let threads = parse_arg_or(&args, &["-p", "--threads"], 1);
                    run_or_exit(diamond::commands::makedb::run(&input_files, &db, threads));
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
                let seq_ids = get_arg(&args, &["--seq"]);
                let output = get_arg(&args, &["-o", "--out"]).map(PathBuf::from);
                run_or_exit(diamond::commands::getseq::run_with_output(
                    &db,
                    seq_ids.as_deref(),
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
                max_target_seqs: parse_arg_or(&args, &["-k", "--max-target-seqs"], 25),
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
                translated_query_layout: None,
            };
            run_or_exit(diamond::commands::blastp::run(&config));
        }
        "blastx" if !has_flag(&args, "--legacy") && !route_blastx_to_legacy(&args) => {
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
                max_target_seqs: parse_arg_or(&args, &["-k", "--max-target-seqs"], 25),
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
            };
            run_or_exit(diamond::commands::blastx::run(&config));
        }
        "cluster" | "linclust" | "deepclust" if !has_flag(&args, "--legacy") => {
            print_banner();
            let workflow = match args[1].as_str() {
                "linclust" => ClusterWorkflow::LinClust,
                "deepclust" => ClusterWorkflow::DeepClust,
                _ => ClusterWorkflow::Cascaded,
            };
            let config = ClusterConfig {
                workflow,
                database: get_arg(&args, &["-d", "--db"]).unwrap_or_default(),
                output: get_arg(&args, &["-o", "--out"]).unwrap_or_default(),
                threads: parse_arg_or(&args, &["-p", "--threads"], 0),
                member_cover: parse_arg_or(&args, &["--member-cover"], 80.0),
                approx_id: get_arg(&args, &["--approx-id"]).and_then(|s| s.parse().ok()),
                cluster_steps: get_all_args(&args, &["--cluster-steps"]),
                alignment_output: get_arg(&args, &["--aln-out"]),
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
                && get_all_args(&args, &["-f", "--outfmt"])
                    .first()
                    .is_some_and(|f| f == "100" || f == "daa") =>
        {
            print_banner();
            let daa_file = get_arg(&args, &["-a", "--daa"]).unwrap_or_default();
            let output = get_arg(&args, &["-o", "--out"]).unwrap_or_default();
            if daa_file.is_empty() {
                eprintln!("Error: -a/--daa argument required");
                std::process::exit(1);
            }
            if output.is_empty() {
                eprintln!("Error: -o/--out argument required for native DAA view output");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_daa(&ViewConfig {
                daa_file,
                output,
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
                max_target_seqs: parse_arg_or(&args, &["-k", "--max-target-seqs"], 25),
                toppercent: get_arg(&args, &["--top"]).and_then(|s| s.parse().ok()),
                forward_only: has_flag(&args, "--forwardonly"),
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: get_all_args(&args, &["--header"]),
            }));
        }
        "view"
            if !has_flag(&args, "--legacy") && {
                let outfmt = get_all_args(&args, &["-f", "--outfmt"]);
                outfmt.is_empty()
                    || outfmt
                        .first()
                        .is_some_and(|f| f == "6" || f == "tab" || f == "104" || f == "json-flat")
            } =>
        {
            print_banner();
            let daa_file = get_arg(&args, &["-a", "--daa"]).unwrap_or_default();
            let output = get_arg(&args, &["-o", "--out"]).unwrap_or_default();
            if daa_file.is_empty() {
                eprintln!("Error: -a/--daa argument required");
                std::process::exit(1);
            }
            if output.is_empty() {
                eprintln!("Error: -o/--out argument required for native tabular view output");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_tabular(&ViewConfig {
                daa_file,
                output,
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
                max_target_seqs: parse_arg_or(&args, &["-k", "--max-target-seqs"], 25),
                toppercent: get_arg(&args, &["--top"]).and_then(|s| s.parse().ok()),
                forward_only: has_flag(&args, "--forwardonly"),
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: get_all_args(&args, &["--header"]),
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && get_all_args(&args, &["-f", "--outfmt"])
                    .first()
                    .is_some_and(|f| f == "103" || f == "paf") =>
        {
            print_banner();
            let daa_file = get_arg(&args, &["-a", "--daa"]).unwrap_or_default();
            let output = get_arg(&args, &["-o", "--out"]).unwrap_or_default();
            if daa_file.is_empty() {
                eprintln!("Error: -a/--daa argument required");
                std::process::exit(1);
            }
            if output.is_empty() {
                eprintln!("Error: -o/--out argument required for native PAF view output");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_paf(&ViewConfig {
                daa_file,
                output,
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
                max_target_seqs: parse_arg_or(&args, &["-k", "--max-target-seqs"], 25),
                toppercent: get_arg(&args, &["--top"]).and_then(|s| s.parse().ok()),
                forward_only: has_flag(&args, "--forwardonly"),
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: get_all_args(&args, &["--header"]),
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && get_all_args(&args, &["-f", "--outfmt"])
                    .first()
                    .is_some_and(|f| f == "101" || f == "sam") =>
        {
            print_banner();
            let daa_file = get_arg(&args, &["-a", "--daa"]).unwrap_or_default();
            let output = get_arg(&args, &["-o", "--out"]).unwrap_or_default();
            if daa_file.is_empty() {
                eprintln!("Error: -a/--daa argument required");
                std::process::exit(1);
            }
            if output.is_empty() {
                eprintln!("Error: -o/--out argument required for native SAM view output");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_sam(&ViewConfig {
                daa_file,
                output,
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
                max_target_seqs: parse_arg_or(&args, &["-k", "--max-target-seqs"], 25),
                toppercent: get_arg(&args, &["--top"]).and_then(|s| s.parse().ok()),
                forward_only: has_flag(&args, "--forwardonly"),
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: get_all_args(&args, &["--header"]),
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && get_all_args(&args, &["-f", "--outfmt"])
                    .first()
                    .is_some_and(|f| f == "0") =>
        {
            print_banner();
            let daa_file = get_arg(&args, &["-a", "--daa"]).unwrap_or_default();
            let output = get_arg(&args, &["-o", "--out"]).unwrap_or_default();
            if daa_file.is_empty() {
                eprintln!("Error: -a/--daa argument required");
                std::process::exit(1);
            }
            if output.is_empty() {
                eprintln!("Error: -o/--out argument required for native pairwise view output");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_pairwise(&ViewConfig {
                daa_file,
                output,
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
                max_target_seqs: parse_arg_or(&args, &["-k", "--max-target-seqs"], 25),
                toppercent: get_arg(&args, &["--top"]).and_then(|s| s.parse().ok()),
                forward_only: has_flag(&args, "--forwardonly"),
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: get_all_args(&args, &["--header"]),
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && get_all_args(&args, &["-f", "--outfmt"])
                    .first()
                    .is_some_and(|f| f == "5" || f == "xml") =>
        {
            print_banner();
            let daa_file = get_arg(&args, &["-a", "--daa"]).unwrap_or_default();
            let output = get_arg(&args, &["-o", "--out"]).unwrap_or_default();
            if daa_file.is_empty() {
                eprintln!("Error: -a/--daa argument required");
                std::process::exit(1);
            }
            if output.is_empty() {
                eprintln!("Error: -o/--out argument required for native XML view output");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_xml(&ViewConfig {
                daa_file,
                output,
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
                max_target_seqs: parse_arg_or(&args, &["-k", "--max-target-seqs"], 25),
                toppercent: get_arg(&args, &["--top"]).and_then(|s| s.parse().ok()),
                forward_only: has_flag(&args, "--forwardonly"),
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: get_all_args(&args, &["--header"]),
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && get_all_args(&args, &["-f", "--outfmt"])
                    .first()
                    .is_some_and(|f| f == "null") =>
        {
            print_banner();
            let daa_file = get_arg(&args, &["-a", "--daa"]).unwrap_or_default();
            let output = get_arg(&args, &["-o", "--out"]).unwrap_or_default();
            if daa_file.is_empty() {
                eprintln!("Error: -a/--daa argument required");
                std::process::exit(1);
            }
            if output.is_empty() {
                eprintln!("Error: -o/--out argument required for native null view output");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_null(&ViewConfig {
                daa_file,
                output,
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
                max_target_seqs: parse_arg_or(&args, &["-k", "--max-target-seqs"], 25),
                toppercent: get_arg(&args, &["--top"]).and_then(|s| s.parse().ok()),
                forward_only: has_flag(&args, "--forwardonly"),
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: get_all_args(&args, &["--header"]),
            }));
        }
        "view"
            if !has_flag(&args, "--legacy")
                && get_all_args(&args, &["-f", "--outfmt"])
                    .first()
                    .is_some_and(|f| f == "edge") =>
        {
            print_banner();
            let daa_file = get_arg(&args, &["-a", "--daa"]).unwrap_or_default();
            let output = get_arg(&args, &["-o", "--out"]).unwrap_or_default();
            if daa_file.is_empty() {
                eprintln!("Error: -a/--daa argument required");
                std::process::exit(1);
            }
            if output.is_empty() {
                eprintln!("Error: -o/--out argument required for native edge view output");
                std::process::exit(1);
            }
            run_or_exit(diamond::commands::view::run_edge(&ViewConfig {
                daa_file,
                output,
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
                max_target_seqs: parse_arg_or(&args, &["-k", "--max-target-seqs"], 25),
                toppercent: get_arg(&args, &["--top"]).and_then(|s| s.parse().ok()),
                forward_only: has_flag(&args, "--forwardonly"),
                report_unaligned: parse_arg_or(&args, &["--unal"], 0) != 0,
                sam_qlen_field: has_flag(&args, "--sam-query-len"),
                invocation: args.join(" "),
                header: get_all_args(&args, &["--header"]),
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
    eprintln!("diamond v2.1.24.178 (C) Max Planck Society for the Advancement of Science, Benjamin J. Buchfink, University of Tuebingen");
    eprintln!("Documentation, support and updates available at http://www.diamondsearch.org");
    eprintln!("Please cite: http://dx.doi.org/10.1038/s41592-021-01101-x Nature Methods (2021)");
    eprintln!();
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
    #[cfg(all(feature = "ffi", not(windows)))]
    println!("Add --legacy to blastp/blastx to use C++ FFI backend.");
    #[cfg(not(all(feature = "ffi", not(windows))))]
    println!("C++ FFI fallback is not available in this build.");
}

#[cfg(all(feature = "ffi", not(windows)))]
fn run_legacy(args: &[&str]) -> i32 {
    diamond::ffi::run(args)
}

#[cfg(not(all(feature = "ffi", not(windows))))]
fn run_legacy(_args: &[&str]) -> i32 {
    eprintln!("Error: C++ FFI fallback is not available in this build.");
    eprintln!("Rebuild on a non-Windows target with `--features ffi` to enable it for testing.");
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

fn route_common_search_to_legacy(args: &[String]) -> bool {
    let requested_output = get_all_args(args, &["-f", "--outfmt"]);
    let unsupported_output = requested_output.first().is_some_and(|format| {
        format != "6" && format != "tab" && format != "100" && format != "daa"
    }) || (!requested_output
        .first()
        .is_some_and(|format| format == "100" || format == "daa")
        && requested_output.iter().skip(1).any(|name| {
            !matches!(
                FieldId::from_name(name),
                Some(
                    FieldId::QSeqId
                        | FieldId::QAcc
                        | FieldId::QAccVer
                        | FieldId::SSeqId
                        | FieldId::SAcc
                        | FieldId::SAccVer
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
                )
            )
        }));

    unsupported_output
        || get_arg(args, &["--masking"]).is_some_and(|m| m.eq_ignore_ascii_case("seg"))
        || get_arg(args, &["--max-hsps"]).is_some_and(|m| m != "1")
        || get_arg(args, &["--comp-based-stats"]).is_some_and(|m| m != "0" && m != "1")
        || get_arg(args, &["--global-ranking"]).is_some_and(|m| m != "0")
        || has_option(
            args,
            &[
                "-b",
                "--block-size",
                "-c",
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
            ],
        )
}

fn route_blastp_to_legacy(args: &[String]) -> bool {
    has_unknown_search_option(args, false) || route_common_search_to_legacy(args)
}

fn route_blastx_to_legacy(args: &[String]) -> bool {
    has_unknown_search_option(args, true)
        || route_common_search_to_legacy(args)
        || get_arg(args, &["-F", "--frameshift"]).is_some_and(|f| f != "0")
        || has_flag(args, "--swipe")
        || has_option(args, &["--min-query-len"])
}

fn parse_arg_or<T: std::str::FromStr>(args: &[String], flags: &[&str], default: T) -> T {
    get_arg(args, flags)
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

const DEFAULT_NATIVE_MEMORY_LIMIT: usize = 16_000_000_000;

fn parse_memory_limit_arg(args: &[String]) -> Option<usize> {
    Some(
        get_arg(args, &["--memory-limit"]).map_or(DEFAULT_NATIVE_MEMORY_LIMIT, |value| {
            parse_byte_size(&value).unwrap_or_else(|error| {
                eprintln!("Error: invalid --memory-limit '{value}': {error}");
                std::process::exit(1);
            })
        }),
    )
}

fn parse_byte_size(value: &str) -> Result<usize, &'static str> {
    let bytes = diamond::util::string::interpret_number(value)
        .map_err(|_| "use K, M, G, or T (for example 4G)")?;
    usize::try_from(bytes).map_err(|_| "value is out of range")
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

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_string()).collect()
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
        assert!(route_blastp_to_legacy(&args(&[
            "diamond",
            "blastp",
            "--outfmt",
            "6",
            "not_a_field"
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
    }
}
