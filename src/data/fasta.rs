use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use crate::basic::value::{
    CharRepresentation, Letter, SequenceType, AMINO_ACID_ALPHABET, MASK_LETTER, NUCLEOTIDE_ALPHABET,
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

/// Matches C++ `read_fasta(reader, f)`.
pub fn read_fasta<R, F>(reader: &mut R, mut f: F) -> io::Result<()>
where
    R: BufRead,
    F: FnMut(&str, &str, u64) -> io::Result<()>,
{
    let mut line = String::new();
    let mut pos = 0u64;
    let mut start = 0u64;
    let mut id = String::new();
    let mut seq = String::new();

    let mut n = reader.read_line(&mut line)?;
    if n == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTA format error: file does not start with '>'",
        ));
    }
    let mut current = line.trim_end_matches('\n').to_string();
    if current.is_empty() || !current.starts_with('>') {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTA format error: file does not start with '>'",
        ));
    }

    loop {
        if current.ends_with('\r') {
            current.pop();
        }
        if !current.is_empty() {
            if let Some(rest) = current.strip_prefix('>') {
                if !id.is_empty() {
                    f(&id, &seq, start)?;
                }
                start = pos;
                id.clear();
                id.push_str(rest);
                if id.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("FASTA format error: empty id at file offset {pos}"),
                    ));
                }
                seq.clear();
            } else {
                seq.push_str(&current);
            }
        }

        pos += n as u64;
        line.clear();
        n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        current = line.trim_end_matches('\n').to_string();
    }

    if !id.is_empty() {
        f(&id, &seq, start)?;
    }
    Ok(())
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
pub fn read_fasta_file(path: &Path, seq_type: SequenceType) -> io::Result<Vec<FastaRecord>> {
    use std::io::Seek;
    let mut file = std::fs::File::open(path)?;

    // C++ `detect_compressor` (`diamond/src/util/io/input_file.cpp:52-61`)
    // sniffs the gzip magic `\x1F\x8B` from the first two bytes of the file,
    // regardless of file extension. A gzip-compressed query named
    // `query.fasta` (no `.gz`) would decompress fine in C++ but used to be
    // read as binary garbage here. Peek the first two bytes to detect gzip
    // magic before falling back to the extension check.
    let mut magic = [0u8; 2];
    let is_gzip = match file.read_exact(&mut magic) {
        Ok(()) => {
            // Rewind so the actual reader sees the whole file.
            file.seek(std::io::SeekFrom::Start(0))?;
            magic == [0x1F, 0x8B]
        }
        Err(_) => {
            // File shorter than 2 bytes — fall back to extension check.
            file.seek(std::io::SeekFrom::Start(0))?;
            path.extension().is_some_and(|ext| ext == "gz")
        }
    };

    // Use `MultiGzDecoder` so concatenated gzip members (e.g.
    // `cat a.fa.gz b.fa.gz > combined.fa.gz`) decompress past the first
    // member — `GzDecoder` stops at the first member's EOF and silently
    // truncates the rest. C++ uses zlib which decompresses all members by
    // default.
    let reader: Box<dyn Read> = if is_gzip {
        Box::new(flate2::read::MultiGzDecoder::new(file))
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
    let mut buf_reader = BufReader::new(reader);
    let mut records = Vec::new();
    {
        let buf = buf_reader.fill_buf()?;
        if buf.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Error detecting input file format. Input file seems to be empty.",
            ));
        }
        if buf.first().copied() == Some(b'>') {
            read_fasta(&mut buf_reader, |id, seq, _start| {
                if seq.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Missing fields in input line",
                    ));
                }
                let mut sequence = Vec::with_capacity(seq.len());
                for ch in seq.bytes() {
                    sequence.push(
                        encoding.convert(ch).map_err(|e| {
                            io::Error::new(io::ErrorKind::InvalidData, e.to_string())
                        })?,
                    );
                }
                records.push(FastaRecord {
                    id: id.to_string(),
                    sequence,
                });
                Ok(())
            })?;
            return Ok(records);
        }
        if buf.first().copied() != Some(b'@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Error detecting input file format. First line must begin with '>' (FASTA) or '@' (FASTQ).",
            ));
        }
    }

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
        } else if line.starts_with('@')
            && (records.is_empty() && current_id.is_none() || fastq_state == 0)
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
            let id = header.split_whitespace().next().unwrap_or("").to_string();
            current_id = Some(id);
            fastq_state = 1;
        } else if is_fastq {
            match fastq_state {
                1 => {
                    // Sequence line(s) — continue until we hit a '+' line.
                    // Bail on invalid characters like C++ does (its
                    // `CharRepresentation::operator()` throws on unknown
                    // chars). The previous `if let Ok(letter)` silently
                    // dropped invalid bytes — that's a faithfulness gap that
                    // also masks real data corruption from the user.
                    if line.starts_with('+') {
                        fastq_state = 3; // skip quality line next
                    } else {
                        for ch in line.bytes() {
                            let letter = encoding.convert(ch).map_err(|e| {
                                io::Error::new(io::ErrorKind::InvalidData, e.to_string())
                            })?;
                            current_seq.push(letter);
                        }
                    }
                }
                2 => {
                    // Plus line (unused — handled in state 1)
                    fastq_state = 3;
                }
                3 => {
                    // Quality line - skip
                    fastq_state = 0;
                }
                _ => {}
            }
        } else {
            // FASTA sequence line — error on invalid chars (matches C++).
            for ch in line.bytes() {
                let letter = encoding
                    .convert(ch)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
                current_seq.push(letter);
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
        assert_eq!(records[0].id, "seq1 description here");
        assert_eq!(records[0].sequence.len(), 8);
    }

    #[test]
    fn test_read_fasta_callback_matches_cpp_titles_offsets_and_blank_lines() {
        let input = b">seq1 description\r\nARND\n\n>seq2\nCQ\n";
        let mut rows = Vec::new();
        let mut reader = BufReader::new(&input[..]);
        read_fasta(&mut reader, |id, seq, start| {
            rows.push((id.to_string(), seq.to_string(), start));
            Ok(())
        })
        .unwrap();
        assert_eq!(
            rows,
            vec![
                ("seq1 description".to_string(), "ARND".to_string(), 0),
                ("seq2".to_string(), "CQ".to_string(), 25)
            ]
        );
    }

    #[test]
    fn test_read_fasta_callback_rejects_cpp_format_errors() {
        let mut no_header = BufReader::new(&b"seq\n"[..]);
        let err = read_fasta(&mut no_header, |_, _, _| Ok(())).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            err.to_string(),
            "FASTA format error: file does not start with '>'"
        );

        let mut empty_id = BufReader::new(&b">\nARND\n"[..]);
        let err = read_fasta(&mut empty_id, |_, _, _| Ok(())).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            err.to_string(),
            "FASTA format error: empty id at file offset 0"
        );
    }

    #[test]
    fn test_read_fasta_with_encoding_detects_empty_input_like_cpp() {
        let err = read_fasta_amino_acid(&b""[..]).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            err.to_string(),
            "Error detecting input file format. Input file seems to be empty."
        );
    }

    #[test]
    fn test_read_fasta_with_encoding_rejects_unknown_format_like_cpp() {
        let err = read_fasta_amino_acid(&b"seq\n"[..]).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            err.to_string(),
            "Error detecting input file format. First line must begin with '>' (FASTA) or '@' (FASTQ)."
        );
    }

    #[test]
    fn test_read_fasta_with_encoding_rejects_empty_sequence_like_cpp() {
        let err = read_fasta_amino_acid(&b">seq\n"[..]).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(err.to_string(), "Missing fields in input line");
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
