use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use crate::basic::value::{
    CharRepresentation, Letter, SequenceType, AMINO_ACID_ALPHABET, MASK_LETTER, NUCLEOTIDE_ALPHABET,
};

pub mod fasta_file;
pub use fasta_file::{FastaEntry, FastaFile, FastaFileConfig};

/// A parsed FASTA/FASTQ record.
#[derive(Debug, Clone)]
pub struct FastaRecord<T = String> {
    /// FASTQ quality values; absent for FASTA/database records.
    pub quality: Option<T>,
    /// Sequence identifier (header line without '>').
    pub id: T,
    /// Encoded sequence as Letter values.
    pub sequence: Vec<Letter>,
}

/// Actual Source strings are arbitrary bytes. Native CLI owners use this
/// specialization; the default String record remains a checked Rust API.
pub type ByteFastaRecord = FastaRecord<Vec<u8>>;
impl From<FastaRecord> for ByteFastaRecord {
    fn from(record: FastaRecord) -> Self {
        Self {
            id: record.id.into_bytes(),
            sequence: record.sequence,
            quality: record.quality.map(String::into_bytes),
        }
    }
}
impl TryFrom<ByteFastaRecord> for FastaRecord {
    type Error = std::string::FromUtf8Error;
    fn try_from(record: ByteFastaRecord) -> Result<Self, Self::Error> {
        Ok(Self {
            id: String::from_utf8(record.id)?,
            sequence: record.sequence,
            quality: record.quality.map(String::from_utf8).transpose()?,
        })
    }
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
    read_fasta_bytes(reader, |id, sequence, offset| {
        let id = std::str::from_utf8(id)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
        let sequence = std::str::from_utf8(sequence)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
        f(id, sequence, offset)
    })
}

pub fn read_fasta_bytes<R, F>(reader: &mut R, mut f: F) -> io::Result<()>
where
    R: BufRead,
    F: FnMut(&[u8], &[u8], u64) -> io::Result<()>,
{
    let mut line = Vec::new();
    let mut pos = 0u64;
    let mut start = 0u64;
    let mut id = Vec::new();
    let mut seq = Vec::new();

    let mut n = reader.read_until(b'\n', &mut line)?;
    if n == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTA format error: file does not start with '>'",
        ));
    }
    let mut current = line.strip_suffix(b"\n").unwrap_or(&line).to_vec();
    if current.is_empty() || !current.starts_with(b">") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "FASTA format error: file does not start with '>'",
        ));
    }

    loop {
        if current.ends_with(b"\r") {
            current.pop();
        }
        if !current.is_empty() {
            if let Some(rest) = current.strip_prefix(b">") {
                if !id.is_empty() {
                    f(&id, &seq, start)?;
                }
                start = pos;
                id.clear();
                id.extend_from_slice(rest);
                if id.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("FASTA format error: empty id at file offset {pos}"),
                    ));
                }
                seq.clear();
            } else {
                seq.extend_from_slice(&current);
            }
        }

        pos += n as u64;
        line.clear();
        n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        current = line.strip_suffix(b"\n").unwrap_or(&line).to_vec();
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
    let mut records = Vec::new();
    for_each_fasta_file(path, seq_type, |record| {
        records.push(record);
        Ok(())
    })?;
    Ok(records)
}

pub fn read_fasta_file_bytes(
    path: &Path,
    seq_type: SequenceType,
) -> io::Result<Vec<ByteFastaRecord>> {
    let mut records = Vec::new();
    for_each_fasta_file_with_buffer_size_bytes(path, seq_type, 8192, |record| {
        records.push(record);
        Ok(())
    })?;
    Ok(records)
}
pub fn for_each_fasta_file_with_buffer_size_bytes<F>(
    path: &Path,
    seq_type: SequenceType,
    buffer_size: usize,
    f: F,
) -> io::Result<()>
where
    F: FnMut(ByteFastaRecord) -> io::Result<()>,
{
    let input = prepare_source_sequence_reader(path, buffer_size)?;
    for_each_prepared_sequence_reader_bytes(input, seq_type, buffer_size, f)
}

/// FASTA/FASTQ autodetection for an already owned stream (including stdin).
/// This typed convenience API explicitly autodetects supported codecs, including
/// for borrowed cursors. Native file/stdin call sites use Source path policy.
pub fn read_fasta_reader<R: Read>(
    reader: R,
    seq_type: SequenceType,
) -> io::Result<Vec<FastaRecord>> {
    let mut records = Vec::new();
    for_each_fasta_decoded(
        reader,
        seq_type,
        8192,
        crate::util::io::std_stream::DecodeMode::NATIVE_AUTO,
        "",
        |record| {
            records.push(record);
            Ok(())
        },
    )?;
    Ok(records)
}

/// Open the actual Source InputFile owner and detect the decoded format once.
pub(crate) fn prepare_source_sequence_reader(
    path: &Path,
    buffer_size: usize,
) -> io::Result<crate::util::io::InputFile> {
    let name = path
        .to_str()
        .ok_or_else(|| io::Error::other("Non-UTF-8 sequence input path"))?;
    let input = crate::util::io::InputFile::with_source_profile(name, 0, buffer_size)
        .map_err(io::Error::other)?;
    prepare_decoded_sequence_reader(input)
}

/// Validate a retained decoded owner without installing another codec/buffer.
pub(crate) fn prepare_decoded_sequence_reader(
    mut input: crate::util::io::InputFile,
) -> io::Result<crate::util::io::InputFile> {
    let first = input.peek_bytes(1).map_err(io::Error::other)?;
    match first.first() {
        Some(b'>' | b'@') => Ok(input),
        None => Err(io::Error::other(
            "Error detecting input file format. Input file seems to be empty.",
        )),
        _ => Err(io::Error::other(
            "Error detecting input file format. First line must begin with '>' (FASTA) or '@' (FASTQ).",
        )),
    }
}

/// Parse an already decoded/prepared owner without codec probing or transport conversion.
pub(crate) fn for_each_prepared_sequence_reader<
    R: Read,
    F: FnMut(FastaRecord) -> io::Result<()>,
>(
    reader: R,
    seq_type: SequenceType,
    buffer_size: usize,
    mut f: F,
) -> io::Result<()> {
    for_each_prepared_sequence_reader_bytes(reader, seq_type, buffer_size, |record| {
        f(FastaRecord::try_from(record)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?)
    })
}

pub(crate) fn for_each_prepared_sequence_reader_bytes<
    R: Read,
    F: FnMut(ByteFastaRecord) -> io::Result<()>,
>(
    reader: R,
    seq_type: SequenceType,
    buffer_size: usize,
    f: F,
) -> io::Result<()> {
    let encoding = match seq_type {
        SequenceType::AminoAcid => {
            CharRepresentation::new(AMINO_ACID_ALPHABET, MASK_LETTER, b"UO-")
        }
        SequenceType::Nucleotide => CharRepresentation::new(NUCLEOTIDE_ALPHABET, 4, b"MRWSYKVHDBX"),
    };
    let mut reader = BufReader::with_capacity(buffer_size.max(1), reader);
    if reader.fill_buf()?.first() == Some(&b'@') {
        for_each_source_fastq(&mut reader, &encoding, f)
    } else {
        for_each_source_fasta(&mut reader, &encoding, f)
    }
}

fn for_each_source_fasta<R: BufRead, F: FnMut(ByteFastaRecord) -> io::Result<()>>(
    reader: &mut R,
    encoding: &CharRepresentation,
    mut emit: F,
) -> io::Result<()> {
    loop {
        let mut raw = Vec::new();
        loop {
            let n = reader.read_until(b'\n', &mut raw)?;
            if n == 0 {
                break;
            }
            if raw.last() == Some(&b'\n') {
                let next = reader.fill_buf()?.first().copied();
                if next.is_none() || next == Some(b'>') {
                    // Source consumes the separating LF, leaves the next '>'
                    // pending, and excludes that LF from this TSV record.
                    raw.pop();
                    break;
                }
            }
        }
        if raw.iter().all(u8::is_ascii_whitespace) {
            return Ok(());
        }
        if raw.first() != Some(&b'>') {
            return Err(io::Error::other("Missing fields in input line"));
        }
        let separator = raw
            .iter()
            .position(|&byte| byte == b'\n')
            .ok_or_else(|| io::Error::other("Missing fields in input line"))?;
        if separator + 1 == raw.len() {
            return Err(io::Error::other("Missing fields in input line"));
        }
        let mut id = raw[1..separator].to_vec();
        if id.last() == Some(&b'\r') {
            id.pop();
        }
        let sequence = raw[separator + 1..]
            .iter()
            .copied()
            .filter(|byte| !matches!(byte, b'\n' | b'\r'))
            .map(|byte| {
                encoding
                    .convert(byte)
                    .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
            })
            .collect::<io::Result<Vec<_>>>()?;
        emit(ByteFastaRecord {
            id,
            sequence,
            quality: None,
        })?;
        if reader.fill_buf()?.is_empty() {
            return Ok(());
        }
    }
}

/// Source first assembles a complete stream record, then applies its three-field
/// TSV schema. In particular, multiline sequence errors precede letter decoding
/// but follow the complete quality read (tokenizer_dyn.h:297 and :249).
fn for_each_source_fastq<R: BufRead, F: FnMut(ByteFastaRecord) -> io::Result<()>>(
    reader: &mut R,
    encoding: &CharRepresentation,
    mut emit: F,
) -> io::Result<()> {
    fn line<R: BufRead>(reader: &mut R, number: &mut u64) -> io::Result<(Vec<u8>, bool)> {
        let mut bytes = Vec::new();
        reader.read_until(b'\n', &mut bytes)?;
        let terminated = bytes.last() == Some(&b'\n');
        if terminated {
            bytes.pop();
        }
        if bytes.last() == Some(&b'\r') {
            bytes.pop();
        }
        *number += 1;
        Ok((bytes, !terminated))
    }
    let mut number = 0;
    loop {
        let (header, eof) = line(reader, &mut number)?;
        if eof && header.is_empty() {
            return Ok(());
        }
        let start = number;
        let malformed = || {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Malformed FASTQ record on line {start}"),
            )
        };
        if header.first() != Some(&b'@') || eof {
            return Err(malformed());
        }
        let mut raw = header;
        raw.push(b'\n');
        let mut len = 0;
        loop {
            let (bytes, eof) = line(reader, &mut number)?;
            if eof || bytes.is_empty() {
                return Err(malformed());
            }
            let separator = bytes.first() == Some(&b'+');
            if !separator {
                len += bytes.len();
            }
            raw.extend(bytes);
            raw.push(b'\n');
            if separator {
                break;
            }
        }
        let mut quality_len = 0;
        let eof = loop {
            let (bytes, eof) = line(reader, &mut number)?;
            quality_len += bytes.len();
            raw.extend(bytes);
            raw.push(b'\n');
            if eof || quality_len >= len {
                break eof;
            }
        };
        let (id, sequence, quality) = fasta_file::source_fastq_fields(&raw)?;
        let quality = quality
            .into_iter()
            .filter(|byte| !matches!(byte, b'\r' | b'\n'))
            .collect();
        let sequence = sequence
            .into_iter()
            .map(|byte| {
                encoding
                    .convert(byte)
                    .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
            })
            .collect::<io::Result<Vec<_>>>()?;
        emit(ByteFastaRecord {
            id,
            sequence,
            quality: Some(quality),
        })?;
        if eof {
            return Ok(());
        }
    }
}

/// Parse a forward-only input, detecting compression before record decoding.
pub fn for_each_fasta_reader<R: Read, F: FnMut(FastaRecord) -> io::Result<()>>(
    reader: R,
    seq_type: SequenceType,
    callback: F,
) -> io::Result<()> {
    for_each_fasta_decoded(
        reader,
        seq_type,
        8192,
        crate::util::io::std_stream::DecodeMode::NATIVE_AUTO,
        "",
        callback,
    )
}

/// Stream an opened sequence source with the configured staging capacity.
/// This supports non-seekable inputs such as standard input and named pipes.
pub fn for_each_fasta_reader_with_buffer_size<R: Read, F: FnMut(FastaRecord) -> io::Result<()>>(
    reader: R,
    path: &Path,
    seq_type: SequenceType,
    buffer_size: usize,
    callback: F,
) -> io::Result<()> {
    for_each_fasta_autodetect_with_buffer_size(reader, path, seq_type, buffer_size, callback)
}

/// Read a FASTA/FASTQ file one record at a time.
///
/// Unlike [`read_fasta_file`], this keeps only the parser's current record in
/// memory.  `makedb` uses it to reproduce `FastaFile::load_seqs`' bounded
/// block lifetime instead of retaining the complete decoded input.
pub fn for_each_fasta_file<F>(path: &Path, seq_type: SequenceType, f: F) -> io::Result<()>
where
    F: FnMut(FastaRecord) -> io::Result<()>,
{
    for_each_fasta_file_with_buffer_size(path, seq_type, 8192, f)
}

/// Search-side explicit input staging capacity, replacing source global config.
pub fn for_each_fasta_file_with_buffer_size<F>(
    path: &Path,
    seq_type: SequenceType,
    buffer_size: usize,
    f: F,
) -> io::Result<()>
where
    F: FnMut(FastaRecord) -> io::Result<()>,
{
    let input = prepare_source_sequence_reader(path, buffer_size)?;
    for_each_prepared_sequence_reader(input, seq_type, buffer_size, f)
}

#[cfg(test)]
fn read_fasta_autodetect<R: Read>(
    reader: R,
    path: &Path,
    seq_type: SequenceType,
) -> io::Result<Vec<FastaRecord>> {
    let mut records = Vec::new();
    for_each_fasta_decoded(
        reader,
        seq_type,
        8192,
        crate::util::io::std_stream::DecodeMode::NATIVE_AUTO,
        path.to_string_lossy(),
        |record| {
            records.push(record);
            Ok(())
        },
    )?;
    Ok(records)
}

fn for_each_fasta_autodetect_with_buffer_size<R: Read, F: FnMut(FastaRecord) -> io::Result<()>>(
    reader: R,
    path: &Path,
    seq_type: SequenceType,
    buffer_size: usize,
    f: F,
) -> io::Result<()> {
    let mode = crate::util::io::std_stream::source_path_mode(path)?;
    for_each_fasta_decoded(
        reader,
        seq_type,
        buffer_size,
        mode,
        path.to_string_lossy(),
        f,
    )
}

fn for_each_fasta_decoded<R: Read, F: FnMut(FastaRecord) -> io::Result<()>>(
    reader: R,
    seq_type: SequenceType,
    buffer_size: usize,
    mode: crate::util::io::std_stream::DecodeMode,
    label: impl Into<String>,
    f: F,
) -> io::Result<()> {
    let decoded =
        crate::util::io::std_stream::decode_source_format_reader(reader, buffer_size, mode, label)?;
    match seq_type {
        SequenceType::AminoAcid => for_each_fasta_with_encoding_buffered(
            decoded,
            &CharRepresentation::new(AMINO_ACID_ALPHABET, MASK_LETTER, b"UO-"),
            buffer_size,
            f,
        ),
        SequenceType::Nucleotide => for_each_fasta_with_encoding_buffered(
            decoded,
            &CharRepresentation::new(NUCLEOTIDE_ALPHABET, 4, b"MRWSYKVHDBX"),
            buffer_size,
            f,
        ),
    }
}

fn read_fasta_with_encoding<R: Read>(
    reader: R,
    encoding: &CharRepresentation,
) -> io::Result<Vec<FastaRecord>> {
    let mut records = Vec::new();
    for_each_fasta_with_encoding(reader, encoding, |record| {
        records.push(record);
        Ok(())
    })?;
    Ok(records)
}

fn for_each_fasta_with_encoding<R, F>(
    reader: R,
    encoding: &CharRepresentation,
    f: F,
) -> io::Result<()>
where
    R: Read,
    F: FnMut(FastaRecord) -> io::Result<()>,
{
    for_each_fasta_with_encoding_buffered(reader, encoding, 8192, f)
}

fn for_each_fasta_with_encoding_buffered<R: Read, F: FnMut(FastaRecord) -> io::Result<()>>(
    reader: R,
    encoding: &CharRepresentation,
    buffer_size: usize,
    mut f: F,
) -> io::Result<()> {
    let mut buf_reader = BufReader::with_capacity(buffer_size.max(1), reader);
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
                f(FastaRecord {
                    quality: None,
                    id: id.to_string(),
                    sequence,
                })
            })?;
            return Ok(());
        }
        if buf.first().copied() != Some(b'@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Error detecting input file format. First line must begin with '>' (FASTA) or '@' (FASTQ).",
            ));
        }
    }

    // FastqTokenizer reads one sequence line, then enough quality lines to
    // cover its length. EOF may leave a short quality string in the source.
    let mut id: Option<String> = None;
    let mut sequence = Vec::new();
    let mut quality = String::new();
    let mut state = 0;
    let mut header_line = 1;
    for (index, line) in buf_reader.lines().enumerate() {
        let line = line?;
        let line = line.trim_end_matches('\r');
        match state {
            0 => {
                if let Some(previous) = id.take() {
                    f(FastaRecord {
                        id: previous,
                        sequence: std::mem::take(&mut sequence),
                        quality: Some(std::mem::take(&mut quality)),
                    })?;
                }
                header_line = index + 1;
                id = Some(
                    line.strip_prefix('@')
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("Malformed FASTQ record on line {header_line}"),
                            )
                        })?
                        .to_owned(),
                );
                state = 1;
            }
            1 => {
                if line.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("Malformed FASTQ record on line {header_line}"),
                    ));
                }
                for letter in line.bytes() {
                    sequence.push(
                        encoding
                            .convert(letter)
                            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?,
                    );
                }
                state = 2;
            }
            2 => {
                if !line.starts_with('+') {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "Malformed FASTQ record",
                    ));
                }
                state = 3;
            }
            3 => {
                quality.push_str(line);
                if quality.len() >= sequence.len() {
                    state = 0;
                }
            }
            _ => unreachable!(),
        }
    }
    if state == 1 || state == 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Malformed FASTQ record on line {header_line}"),
        ));
    }
    if let Some(id) = id {
        f(FastaRecord {
            id,
            sequence,
            quality: Some(quality),
        })?;
    }
    Ok(())
}

/// Write FASTA records to a writer.
pub fn write_fasta<W: io::Write>(
    writer: &mut W,
    records: &[FastaRecord],
    alphabet: &[u8],
    line_width: usize,
) -> io::Result<()> {
    write_fasta_records(writer, records, alphabet, line_width)
}

pub fn write_fasta_bytes<W: io::Write>(
    writer: &mut W,
    records: &[ByteFastaRecord],
    alphabet: &[u8],
    line_width: usize,
) -> io::Result<()> {
    write_fasta_records(writer, records, alphabet, line_width)
}

fn write_fasta_records<W: io::Write, T: AsRef<[u8]>>(
    writer: &mut W,
    records: &[FastaRecord<T>],
    alphabet: &[u8],
    line_width: usize,
) -> io::Result<()> {
    for record in records {
        writer.write_all(b">")?;
        // FastaFile::write_seq passes id.c_str() to Util::Seq::format.
        writer.write_all(crate::util::sequence::c_string_bytes(record.id.as_ref()))?;
        writer.write_all(b"\n")?;
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
    #[test]
    fn raw_fasta_index_callback_counts_original_header_bytes() {
        let bytes = b">first\xff\r\nAR\n>second\x80\nNC\n";
        let mut reader = std::io::Cursor::new(bytes);
        let mut rows = Vec::new();
        super::read_fasta_bytes(&mut reader, |id, sequence, offset| {
            rows.push((id.to_vec(), sequence.to_vec(), offset));
            Ok(())
        })
        .unwrap();
        assert_eq!(
            rows,
            vec![
                (b"first\xff".to_vec(), b"AR".to_vec(), 0),
                (b"second\x80".to_vec(), b"NC".to_vec(), 12)
            ]
        );
    }

    use super::*;
    use crate::basic::value::AMINO_ACID_ALPHABET;

    struct NonSeekable<R>(R);

    impl<R: Read> Read for NonSeekable<R> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            self.0.read(buf)
        }
    }

    #[test]
    fn prepared_byte_records_preserve_titles_and_quality_without_utf8() {
        let mut records = Vec::new();
        for_each_prepared_sequence_reader_bytes(
            b">q\xff desc\0tail\nARND\n".as_slice(),
            SequenceType::AminoAcid,
            2,
            |record| {
                records.push(record);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(records[0].id, b"q\xff desc\0tail");
        assert!(FastaRecord::<String>::try_from(records.pop().unwrap()).is_err());
        for_each_prepared_sequence_reader_bytes(
            b"@q\xfe\nARND\n+\n!\xff!?\n".as_slice(),
            SequenceType::AminoAcid,
            2,
            |record| {
                records.push(record);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(records[0].id, b"q\xfe");
        assert_eq!(records[0].quality.as_deref(), Some(b"!\xff!?".as_slice()));
    }

    #[test]
    fn source_prepared_fastq_uses_record_then_schema_errors() {
        let mut records = Vec::new();
        let mut parse = |bytes: &[u8]| {
            for_each_prepared_sequence_reader(bytes, SequenceType::AminoAcid, 2, |record| {
                records.push(record);
                Ok(())
            })
        };
        assert_eq!(
            parse(b"@q\n?\nARND\n+\n!!!!!\n").unwrap_err().to_string(),
            "Malformed FASTQ record"
        );
        assert_eq!(
            parse(b"@q\n+\n!\n").unwrap_err().to_string(),
            "Missing fields in input line"
        );
        assert_eq!(
            parse(b"@q\nARND\n+").unwrap_err().to_string(),
            "Malformed FASTQ record on line 1"
        );
        parse(b"@q\n@ARND\n+\n!!!!!\n").unwrap();
        parse(b"@q\n@ARND\nQ\n+\n!!!!!!\n").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].sequence, [0, 1, 2, 3]);
        assert_eq!(records[0].quality.as_deref(), Some("!!!!!"));
        assert_eq!(records[1].quality.as_deref(), Some("Q"));
    }

    #[test]
    fn autodetect_reads_plain_nonseekable_input_without_consuming_prefix() {
        let input = NonSeekable(std::io::Cursor::new(b">pipe\nARND\n"));
        let records =
            read_fasta_autodetect(input, Path::new("/dev/stdin"), SequenceType::AminoAcid).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "pipe");
        assert_eq!(records[0].sequence, vec![0, 1, 2, 3]);
    }

    #[test]
    fn autodetect_reads_gzip_nonseekable_input_without_consuming_magic() {
        use std::io::Write;

        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), Default::default());
        encoder.write_all(b">pipe-gz\nCQEG\n").unwrap();
        let input = NonSeekable(std::io::Cursor::new(encoder.finish().unwrap()));
        let records =
            read_fasta_autodetect(input, Path::new("/dev/stdin"), SequenceType::AminoAcid).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "pipe-gz");
        assert_eq!(records[0].sequence, vec![4, 5, 6, 7]);
    }

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
    fn test_write_fasta() {
        let records = vec![FastaRecord {
            quality: None,
            id: "seq1".to_string(),
            sequence: vec![0, 1, 2, 3],
        }];
        let mut buf = Vec::new();
        write_fasta(&mut buf, &records, AMINO_ACID_ALPHABET, 80).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), ">seq1\nARND\n");
    }
}
