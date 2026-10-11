use std::io;
use std::path::Path;
use std::time::Instant;

use crate::basic::value::SequenceType;
use crate::data::db_builder;

/// Run the makedb command — build a DIAMOND database from FASTA file(s).
pub fn run(input_files: &[String], database: &str, threads: i32) -> io::Result<()> {
    let start = Instant::now();

    // C++ rejects multiple `--in` values with "Too many arguments provided
    // for option --in." (`diamond/src/data/dmnd/dmnd.cpp:261-262`). Rust used
    // to silently concatenate; mirror C++.
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

    let stats = db_builder::build_db(&input_refs, &db_path, SequenceType::AminoAcid)?;

    let elapsed = start.elapsed();
    eprintln!(
        "Database sequences: {}, letters: {}, time: {:.1}s",
        stats.sequences,
        stats.letters,
        elapsed.as_secs_f64()
    );

    Ok(())
}
