use std::env;

use diamond::commands::blastp::BlastpConfig;
use diamond::commands::blastx::BlastxConfig;
use diamond::commands::cluster_cmd::ClusterConfig;
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
        "blastp" if !has_flag(&args, "--legacy") => {
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
                min_id: parse_arg_or(&args, &["--id"], 0.0),
                threads: parse_arg_or(&args, &["-p", "--threads"], 1),
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
                sensitivity: parse_sensitivity(&args),
            };
            run_or_exit(diamond::commands::blastp::run(&config));
        }
        "blastx" if !has_flag(&args, "--legacy") => {
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
                min_id: parse_arg_or(&args, &["--id"], 0.0),
                threads: parse_arg_or(&args, &["-p", "--threads"], 1),
                outfmt: get_all_args(&args, &["-f", "--outfmt"]),
                sensitivity: parse_sensitivity(&args),
                query_gencode: parse_arg_or(&args, &["--query-gencode"], 1),
                strand: get_arg(&args, &["--strand"]).unwrap_or_else(|| "both".into()),
                min_orf: get_arg(&args, &["--min-orf"]).and_then(|s| s.parse().ok()),
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
        _ => {
            // Fall back to C++ FFI for full compatibility
            // Filter out --legacy flag which is not known to C++
            let filtered: Vec<&str> = args.iter()
                .map(|s| s.as_str())
                .filter(|s| *s != "--legacy")
                .collect();
            let code = diamond::ffi::run(&filtered);
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
    println!("  version    Show version");
    println!("  test       Run regression tests");
    println!();
    println!("Use 'diamond COMMAND --help' for command-specific options.");
    println!("Add --legacy to blastp/blastx to use C++ FFI backend.");
}

fn run_or_exit(result: std::io::Result<()>) {
    if let Err(e) = result {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

fn get_arg(args: &[String], flags: &[&str]) -> Option<String> {
    args.iter()
        .position(|a| flags.contains(&a.as_str()))
        .and_then(|i| args.get(i + 1))
        .cloned()
}

fn get_all_args(args: &[String], flags: &[&str]) -> Vec<String> {
    let mut values = Vec::new();
    let mut i = 0;
    while i < args.len() {
        if flags.contains(&args[i].as_str()) {
            if let Some(val) = args.get(i + 1) {
                values.push(val.clone());
                i += 2;
                continue;
            }
        }
        i += 1;
    }
    values
}

fn has_flag(args: &[String], flag: &str) -> bool {
    args.iter().any(|a| a == flag)
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
