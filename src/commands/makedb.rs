use std::io;
use std::path::Path;
use std::time::Instant;

use crate::basic::value::SequenceType;
use crate::data::db_builder;

/// Run the makedb command — build a DIAMOND database from FASTA file(s).
pub fn run(input_files: &[String], database: &str, threads: i32) -> io::Result<()> {
    let start = Instant::now();

    // Configure thread pool
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads as usize)
            .build_global()
            .ok(); // ignore if already initialized
    }

    let input_refs: Vec<&str> = input_files.iter().map(|s| s.as_str()).collect();

    // Append .dmnd extension if not present
    let db_path = if Path::new(database).extension().is_none() {
        format!("{}.dmnd", database)
    } else {
        database.to_string()
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::dmnd::{ReferenceHeader, MAGIC_NUMBER};
    use std::fs;

    #[test]
    fn test_makedb_from_test_fasta() {
        let input = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/1.faa"
        );
        let tmp = std::env::temp_dir().join("test_makedb.dmnd");
        let tmp_str = tmp.to_string_lossy().to_string();

        let result = run(&[input.to_string()], &tmp_str, 1);
        assert!(result.is_ok(), "makedb failed: {:?}", result.err());

        // Verify the output is a valid DIAMOND database
        let mut f = fs::File::open(&tmp).unwrap();
        let header = ReferenceHeader::read_from(&mut f).unwrap();
        assert_eq!(header.magic_number, MAGIC_NUMBER);
        assert!(header.sequences > 0);
        assert!(header.letters > 0);

        let _ = fs::remove_file(&tmp);
    }
}
