use std::io::{self, Write};

use crate::basic::value::AMINO_ACID_ALPHABET;
use crate::data::dmnd_reader;

/// Run the getseq command — retrieve sequences from a DIAMOND database.
pub fn run(database: &str, seq_ids: Option<&str>) -> io::Result<()> {
    let (_header, records) = dmnd_reader::read_dmnd_auto(database)?;

    let stdout = io::stdout();
    let mut writer = io::BufWriter::new(stdout.lock());

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
}
