use std::io;
use std::path::Path;
use std::time::Instant;

use crate::basic::value::SequenceType;
use crate::data::db_builder;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MakedbOptions {
    pub taxonmap: Option<String>,
    pub taxonnodes: Option<String>,
    pub taxonnames: Option<String>,
    pub no_parse_seqids: bool,
}

/// Run the makedb command — build a DIAMOND database from FASTA file(s).
pub fn run(input_files: &[String], database: &str, threads: i32) -> io::Result<()> {
    run_with_options(input_files, database, threads, &MakedbOptions::default())
}

pub fn run_with_options(
    input_files: &[String],
    database: &str,
    threads: i32,
    options: &MakedbOptions,
) -> io::Result<()> {
    let start = Instant::now();

    // C++ rejects multiple values in the effective `--in` occurrence with
    // "Too many arguments provided for option --in."
    // (`diamond/src/data/dmnd/dmnd.cpp:261-262`). The CLI parser separately
    // mirrors upstream assignment semantics where a repeated occurrence
    // replaces the earlier one.
    if input_files.len() > 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Too many arguments provided for option --in.",
        ));
    }
    if input_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Missing required option --in.",
        ));
    }

    // Configure thread pool
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads as usize)
            .build_global()
            .ok(); // ignore if already initialized
    }

    let input_refs: Vec<&str> = input_files.iter().map(|s| s.as_str()).collect();

    // Append `.dmnd` to the output path unless it already ends in `.dmnd`.
    // C++ `auto_append_extension(database, ".dmnd")` (`basic/config.cpp:770`)
    // appends regardless of whether some OTHER extension is present (e.g.
    // `nr.foo` → `nr.foo.dmnd`). The prior Rust check only appended when no
    // extension was present at all, so `--db nr.foo` would write `nr.foo`
    // and produce a DB at the wrong path vs C++.
    let db_path = if Path::new(database).extension().is_some_and(|e| e == "dmnd") {
        database.to_string()
    } else {
        format!("{}.dmnd", database)
    };

    eprintln!("Database file: {}", db_path);
    eprintln!("Input files: {}", input_files.join(", "));

    let stats = db_builder::build_db_with_options(
        &input_refs,
        &db_path,
        SequenceType::AminoAcid,
        &db_builder::BuildOptions {
            taxonmap: options.taxonmap.clone(),
            taxonnodes: options.taxonnodes.clone(),
            taxonnames: options.taxonnames.clone(),
            no_parse_seqids: options.no_parse_seqids,
        },
    )?;

    let elapsed = start.elapsed();
    eprintln!(
        "Database sequences: {}, letters: {}, time: {:.1}s",
        stats.sequences,
        stats.letters,
        elapsed.as_secs_f64()
    );

    Ok(())
}
