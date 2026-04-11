use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use crate::basic::value::{
    CharRepresentation, Letter, SequenceType, AMINO_ACID_ALPHABET, MASK_LETTER,
    NUCLEOTIDE_ALPHABET,
};

/// A parsed FASTA/FASTQ record.
#[derive(Debug, Clone)]
pub struct FastaRecord {
    /// Sequence identifier (header line without '>').
    pub id: String,
    /// Encoded sequence as Letter values.
    pub sequence: Vec<Letter>,
}

/// Sequence file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqFileFormat {
    Fasta,
    Fastq,
}

/// Read FASTA records from a reader, encoding letters as amino acids.
pub fn read_fasta_amino_acid<R: Read>(reader: R) -> io::Result<Vec<FastaRecord>> {
    let cr = CharRepresentation::new(AMINO_ACID_ALPHABET, MASK_LETTER, b"UO-");
    read_fasta_with_encoding(reader, &cr)
}

/// Read FASTA records from a reader, encoding letters as nucleotides.
pub fn read_fasta_nucleotide<R: Read>(reader: R) -> io::Result<Vec<FastaRecord>> {
    let cr = CharRepresentation::new(NUCLEOTIDE_ALPHABET, 4, b"MRWSYKVHDBX");
    read_fasta_with_encoding(reader, &cr)
}

/// Read FASTA records from a file path.
pub fn read_fasta_file(
    path: &Path,
    seq_type: SequenceType,
) -> io::Result<Vec<FastaRecord>> {
    let file = std::fs::File::open(path)?;

    // Check for gzip
    let reader: Box<dyn Read> = if path
        .extension()
        .is_some_and(|ext| ext == "gz")
    {
        Box::new(flate2::read::GzDecoder::new(file))
    } else {
        Box::new(file)
    };

    match seq_type {
        SequenceType::AminoAcid => read_fasta_amino_acid(reader),
        SequenceType::Nucleotide => read_fasta_nucleotide(reader),
    }
}

fn read_fasta_with_encoding<R: Read>(
    reader: R,
    encoding: &CharRepresentation,
) -> io::Result<Vec<FastaRecord>> {
    let buf_reader = BufReader::new(reader);
    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq: Vec<Letter> = Vec::new();
    let mut is_fastq = false;
    let mut fastq_state = 0u8; // 0=header, 1=seq, 2=plus, 3=qual

    for line in buf_reader.lines() {
        let line = line?;
        let line = line.trim_end();

        if line.is_empty() {
            continue;
        }

        if !is_fastq && line.starts_with('>') {
            // FASTA header
            if let Some(id) = current_id.take() {
                records.push(FastaRecord {
                    id,
                    sequence: std::mem::take(&mut current_seq),
                });
            }
            // Parse ID: take everything after '>' up to first delimiter
            // Delimiters match C++: whitespace + \x01 (SOH, used in NCBI FASTA)
            let header = &line[1..];
            let id = header
                .split(|c: char| c.is_whitespace() || c == '\x01')
                .next()
                .unwrap_or("")
                .to_string();
            current_id = Some(id);
        } else if line.starts_with('@') && (records.is_empty() && current_id.is_none() || fastq_state == 0)
        {
            // FASTQ header
            is_fastq = true;
            if let Some(id) = current_id.take() {
                records.push(FastaRecord {
                    id,
                    sequence: std::mem::take(&mut current_seq),
                });
            }
            let header = &line[1..];
            let id = header
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_id = Some(id);
            fastq_state = 1;
        } else if is_fastq {
            match fastq_state {
                1 => {
                    // Sequence line
                    for ch in line.bytes() {
                        if let Ok(letter) = encoding.convert(ch) {
                            current_seq.push(letter);
                        }
                    }
                    fastq_state = 2;
                }
                2 => {
                    // Plus line
                    fastq_state = 3;
                }
                3 => {
                    // Quality line - skip
                    fastq_state = 0;
                }
                _ => {}
            }
        } else {
            // FASTA sequence line
            for ch in line.bytes() {
                if let Ok(letter) = encoding.convert(ch) {
                    current_seq.push(letter);
                }
            }
        }
    }

    // Push last record
    if let Some(id) = current_id {
        records.push(FastaRecord {
            id,
            sequence: current_seq,
        });
    }

    Ok(records)
}

/// Write FASTA records to a writer.
pub fn write_fasta<W: io::Write>(
    writer: &mut W,
    records: &[FastaRecord],
    alphabet: &[u8],
    line_width: usize,
) -> io::Result<()> {
    for record in records {
        writeln!(writer, ">{}", record.id)?;
        for chunk in record.sequence.chunks(line_width) {
            for &letter in chunk {
                writer.write_all(&[alphabet[letter as usize]])?;
            }
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::AMINO_ACID_ALPHABET;

    #[test]
    fn test_read_fasta_simple() {
        let input = b">seq1\nARND\n>seq2\nCQEG\n";
        let records = read_fasta_amino_acid(&input[..]).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].sequence, vec![0, 1, 2, 3]);
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].sequence, vec![4, 5, 6, 7]);
    }

    #[test]
    fn test_read_fasta_multiline() {
        let input = b">seq1 description here\nARND\nCQEG\n";
        let records = read_fasta_amino_acid(&input[..]).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].sequence.len(), 8);
    }

    #[test]
    fn test_read_fasta_lowercase() {
        let input = b">seq1\narnd\n";
        let records = read_fasta_amino_acid(&input[..]).unwrap();
        assert_eq!(records[0].sequence, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_read_fasta_masked_chars() {
        let input = b">seq1\nAUOA\n";
        let records = read_fasta_amino_acid(&input[..]).unwrap();
        // U and O should map to MASK_LETTER (23)
        assert_eq!(records[0].sequence, vec![0, MASK_LETTER, MASK_LETTER, 0]);
    }

    #[test]
    fn test_read_nucleotide() {
        let input = b">seq1\nACGTN\n";
        let records = read_fasta_nucleotide(&input[..]).unwrap();
        assert_eq!(records[0].sequence, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn test_read_real_test_file() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/1.faa");
        let records = read_fasta_file(Path::new(path), SequenceType::AminoAcid).unwrap();
        assert!(!records.is_empty());
        assert!(!records[0].sequence.is_empty());
    }

    #[test]
    fn test_write_fasta() {
        let records = vec![FastaRecord {
            id: "seq1".to_string(),
            sequence: vec![0, 1, 2, 3],
        }];
        let mut buf = Vec::new();
        write_fasta(&mut buf, &records, AMINO_ACID_ALPHABET, 80).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), ">seq1\nARND\n");
    }
}
