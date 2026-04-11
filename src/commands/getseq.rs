use std::io::{self, Write};

use crate::basic::value::AMINO_ACID_ALPHABET;
use crate::data::dmnd_reader;

/// Run the getseq command — retrieve sequences from a DIAMOND database.
pub fn run(database: &str, seq_ids: Option<&str>) -> io::Result<()> {
    let (header, records) = dmnd_reader::read_dmnd_auto(database)?;
    eprintln!("Database: {} sequences", header.sequences);

    let stdout = io::stdout();
    let mut writer = io::BufWriter::new(stdout.lock());

    if let Some(ids_str) = seq_ids {
        // Output specific sequences
        let requested: Vec<&str> = ids_str.split(',').collect();
        for record in &records {
            if requested.iter().any(|&id| record.id == id) {
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

fn write_fasta_record<W: Write>(
    writer: &mut W,
    id: &str,
    sequence: &[i8],
) -> io::Result<()> {
    writeln!(writer, ">{}", id)?;
    for chunk in sequence.chunks(60) {
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
        let db = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/data.dmnd"
        );
        // Should succeed without error
        let result = run(db, Some("d1ivsa4"));
        assert!(result.is_ok());
    }
}
