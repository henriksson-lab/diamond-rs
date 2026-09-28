use std::io::{self, Write};
use std::path::Path;

use crate::basic::value::AMINO_ACID_ALPHABET;
use crate::data::dmnd_reader;

/// Run the getseq command — retrieve sequences from a DIAMOND database.
pub fn run(database: &str, seq_ids: Option<&str>) -> io::Result<()> {
    run_with_output(database, seq_ids, None)
}

/// Retrieve sequences and optionally write them to `-o/--out`. Keeping stdout
/// as the `None` case preserves the library API while matching the CLI's
/// upstream output-file behavior.
pub fn run_with_output(
    database: &str,
    seq_ids: Option<&str>,
    output: Option<&Path>,
) -> io::Result<()> {
    let (_header, records) = dmnd_reader::read_dmnd_auto(database)?;

    let mut writer: Box<dyn Write> = if let Some(path) = output {
        Box::new(io::BufWriter::new(std::fs::File::create(path)?))
    } else {
        Box::new(io::BufWriter::new(io::stdout()))
    };

    if let Some(ids_str) = seq_ids {
        // Output specific sequences. Match against the BLAST seqid (truncated
        // at first whitespace / FASTA_HEADER_SEP) rather than the full
        // record.id, mirroring C++ `Util::Seq::seqid` (`sequence/sequence.cpp:74`)
        // used by `get_seq` (`sequence_file.cpp:404`). The full record.id
        // typically includes a description after the accession, so a literal
        // `record.id == "sp|P12345|FOO"` never matches in real databases.
        let requested: Vec<&str> = ids_str.split(',').collect();
        for record in &records {
            let seqid = crate::util::sequence::seqid(&record.id);
            if requested.iter().any(|&id| seqid == id || record.id == id) {
                write_fasta_record(&mut writer, &record.id, &record.sequence)?;
            }
        }
    } else {
        // Output all sequences
        for record in &records {
            write_fasta_record(&mut writer, &record.id, &record.sequence)?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn write_fasta_record<W: Write>(writer: &mut W, id: &str, sequence: &[i8]) -> io::Result<()> {
    writeln!(writer, ">{}", id)?;
    // C++ getseq writes 80-char lines (sequence_file.cpp:get_seq).
    for chunk in sequence.chunks(80) {
        for &letter in chunk {
            let idx = (letter & 0x1F) as usize;
            if idx < AMINO_ACID_ALPHABET.len() {
                writer.write_all(&[AMINO_ACID_ALPHABET[idx]])?;
            } else {
                writer.write_all(b"X")?;
            }
        }
        writeln!(writer)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_getseq() {
        let db = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/data.dmnd");
        // Should succeed without error
        let result = run(db, Some("d1ivsa4"));
        assert!(result.is_ok());
    }

    #[test]
    fn getseq_writes_requested_output_file() {
        let db = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/data.dmnd");
        let output = std::env::temp_dir().join(format!(
            "diamond-getseq-output-test-{}.faa",
            std::process::id()
        ));
        run_with_output(db, Some("d1ivsa4"), Some(&output)).unwrap();
        let fasta = std::fs::read_to_string(&output).unwrap();
        assert!(fasta.starts_with(">d1ivsa4"));
        assert!(fasta.lines().count() > 1);
        std::fs::remove_file(output).unwrap();
    }
}
