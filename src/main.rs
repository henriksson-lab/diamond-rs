use std::env;

use diamond::commands::blastp::BlastpConfig;
use diamond::commands::blastx::BlastxConfig;
use diamond::commands::cluster_cmd::ClusterConfig;
use diamond::commands::view::ViewConfig;
use diamond::config::Sensitivity;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        print_usage();
        return;
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
                run_or_exit(diamond::commands::getseq::run(&db, seq_ids.as_deref()));
            } else {
                eprintln!("Error: -d/--db argument required");
                std::process::exit(1);
            }
        }
        "blastp" if !has_flag(&args, "--legacy") && !route_blastp_to_legacy(&args) => {
            // Native Rust blastp pipeline
            print_banner();
            let config = BlastpConfig {
                query_files: get_all_args(&args, &["-q", "--query"]),
                database: get_arg(&args, &["-d", "--db"]).unwrap_or_default(),
                output: get_arg(&args, &["-o", "--out"]),
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
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
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
            };
            run_or_exit(diamond::commands::blastp::run(&config));
        }
        "blastx" if !has_flag(&args, "--legacy") && !route_blastx_to_legacy(&args) => {
            print_banner();
            let config = BlastxConfig {
                query_files: get_all_args(&args, &["-q", "--query"]),
                database: get_arg(&args, &["-d", "--db"]).unwrap_or_default(),
                output: get_arg(&args, &["-o", "--out"]),
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
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
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
            };
            run_or_exit(diamond::commands::blastx::run(&config));
        }
        "cluster" | "linclust" | "deepclust" if !has_flag(&args, "--legacy") => {
            print_banner();
            let config = ClusterConfig {
                database: get_arg(&args, &["-d", "--db"]).unwrap_or_default(),
                output: get_arg(&args, &["-o", "--out"]).unwrap_or_else(|| "clusters.tsv".into()),
                threads: parse_arg_or(&args, &["-p", "--threads"], 1),
                member_cover: parse_arg_or(&args, &["--member-cover"], 80.0),
                approx_id: parse_arg_or(&args, &["--approx-id"], 0.0),
                sensitivity: parse_sensitivity(&args),
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
    println!("  cluster    Cluster sequences");
    println!("  merge-daa  Merge DAA files");
    println!("  version    Show version");
    println!("  test       Run regression tests");
    println!();
    println!("Use 'diamond COMMAND --help' for command-specific options.");
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
        if flags.contains(&args[i].as_str()) {
            i += 1;
            while i < args.len() && !args[i].starts_with('-') {
                values.push(args[i].clone());
                i += 1;
            }
            continue;
        }
        i += 1;
    }
    values
}

fn has_flag(args: &[String], flag: &str) -> bool {
    args.iter().any(|a| a == flag)
}

fn route_blastp_to_legacy(args: &[String]) -> bool {
    get_all_args(args, &["-f", "--outfmt"])
        .first()
        .is_some_and(|f| f == "100" || f == "daa")
        || get_arg(args, &["--masking"]).is_some_and(|m| m.eq_ignore_ascii_case("seg"))
        || get_arg(args, &["--max-hsps"]).is_some_and(|m| m != "1")
        || get_arg(args, &["--comp-based-stats"]).is_some_and(|m| m != "0" && m != "1")
        || get_arg(args, &["--global-ranking"]).is_some_and(|m| m != "0")
}

fn route_blastx_to_legacy(args: &[String]) -> bool {
    route_blastp_to_legacy(args)
        || get_arg(args, &["-F", "--frameshift"]).is_some_and(|f| f != "0")
        || has_flag(args, "--swipe")
}

fn parse_arg_or<T: std::str::FromStr>(args: &[String], flags: &[&str], default: T) -> T {
    get_arg(args, flags)
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
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
}
