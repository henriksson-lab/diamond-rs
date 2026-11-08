//! Retained raw-byte tokenizer for Source FastaFile's TSV record boundary.
//! Mirrors tokenizer_dyn.h FastaTokenizer/FastqTokenizer and File::read_record;
//! it does not use the eager generic Rust TSV facade or cache sequence payloads.
use crate::data::fasta::SeqFileFormat;
use crate::util::io::{Deserializer, InputFile, InputFileStream};
use std::io;

pub(super) struct RecordReader {
    pub input: Deserializer<InputFileStream>,
    pub format: SeqFileFormat,
    pub eof: bool,
    pending: Option<u8>,
    line: u64,
    pub offset: u64,
    name: String,
    temporary: bool,
    unlinked: bool,
    closed: bool,
}

pub(super) struct RawRecord {
    pub id: Vec<u8>,
    pub sequence: Vec<u8>,
    pub quality: Option<Vec<u8>>,
    pub offset: u64,
}

impl RecordReader {
    pub fn open(input: InputFile, format: Option<SeqFileFormat>) -> io::Result<Self> {
        let name = input.file_name.clone();
        let temporary = input.temp_file;
        let unlinked = input.unlinked;
        let mut reader = Self {
            input: input.into_deserializer(),
            format: format.unwrap_or(SeqFileFormat::Fasta),
            eof: false,
            pending: None,
            line: 0,
            offset: 0,
            name,
            temporary,
            unlinked,
            closed: false,
        };
        if format.is_none() {
            let first = reader.input.peek_bytes(1).map_err(io::Error::other)?;
            reader.format = match first.first().copied() {
                Some(b'>') => SeqFileFormat::Fasta,
                Some(b'@') => SeqFileFormat::Fastq,
                None => return Err(io::Error::other("Error detecting input file format. Input file seems to be empty.")),
                _ => return Err(io::Error::other("Error detecting input file format. First line must begin with '>' (FASTA) or '@' (FASTQ).")),
            };
        }
        Ok(reader)
    }

    pub fn close(&mut self) -> io::Result<()> {
        if self.closed {
            return Ok(());
        }
        self.closed = true;
        self.input.close().map_err(io::Error::other)?;
        if self.temporary && !self.unlinked && std::fs::remove_file(&self.name).is_err() {
            crate::util::io::console::source_eprintln!(
                "Warning: Failed to delete temporary file {}",
                self.name
            );
        }
        Ok(())
    }

    pub fn rewind(&mut self) -> io::Result<()> {
        self.input.rewind().map_err(io::Error::other)?;
        self.pending = None;
        self.eof = false;
        self.line = 0;
        self.offset = 0;
        Ok(())
    }

    pub fn seek(&mut self, position: u64) -> io::Result<()> {
        self.input.seek(position as i64).map_err(io::Error::other)?;
        self.pending = None;
        self.eof = false;
        self.offset = position;
        Ok(())
    }

    fn getline(&mut self) -> io::Result<Vec<u8>> {
        let mut line = Vec::new();
        if let Some(byte) = self.pending.take() {
            line.push(byte);
        }
        let terminated = self
            .input
            .read_to(&mut line, b'\n')
            .map_err(io::Error::other)?;
        self.eof = !terminated;
        self.line = self.line.wrapping_add(1);
        self.offset += line.len() as u64 + u64::from(terminated);
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        Ok(line)
    }

    pub fn read(&mut self) -> io::Result<Option<RawRecord>> {
        if self.eof {
            return Ok(None);
        }
        let offset = self.offset;
        match self.format {
            SeqFileFormat::Fasta => {
                let mut raw = Vec::new();
                if let Some(byte) = self.pending.take() {
                    raw.push(byte);
                }
                let (_next, newlines) = self
                    .input
                    .read_to_record_start(&mut raw, b'\n', b'>')
                    .map_err(io::Error::other)?;
                let internal_lines = raw.iter().filter(|&&byte| byte == b'\n').count() as i64;
                self.offset += raw.len() as u64 + (newlines - internal_lines).max(0) as u64;
                // FastaTokenizer calls TextInputFile::read_to directly. It
                // never updates that class's eof_ or line_count (getline does).
                // Thus FastaFile::eof stays false even after physical EOF;
                // repeated failed reads still increment the owning OID.

                if raw.iter().all(u8::is_ascii_whitespace) {
                    return Ok(None);
                }
                // TSV schema requires TWO fields. A nonheader token occupies
                // only one field, including a mismatched paired input format.
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
                    .collect();
                Ok(Some(RawRecord {
                    id,
                    sequence,
                    quality: None,
                    offset,
                }))
            }
            SeqFileFormat::Fastq => {
                let header = self.getline()?;
                if self.eof && header.is_empty() {
                    return Ok(None);
                }
                let start_line = self.line;
                let malformed =
                    || io::Error::other(format!("Malformed FASTQ record on line {start_line}"));
                if header.first() != Some(&b'@') || self.eof {
                    return Err(malformed());
                }
                let mut raw = header;
                raw.push(b'\n');
                let mut len = 0;
                loop {
                    let line = self.getline()?;
                    if self.eof || line.is_empty() {
                        return Err(malformed());
                    }
                    let separator = line.first() == Some(&b'+');
                    if !separator {
                        len += line.len();
                    }
                    raw.extend(line);
                    raw.push(b'\n');
                    if separator {
                        break;
                    }
                }
                let mut quality_len = 0;
                loop {
                    let line = self.getline()?;
                    quality_len += line.len();
                    raw.extend(line);
                    raw.push(b'\n');
                    if self.eof || quality_len >= len {
                        break;
                    }
                }
                let (id, sequence, quality) = source_fastq_fields(&raw)?;
                Ok(Some(RawRecord {
                    id,
                    sequence,
                    quality: Some(
                        quality
                            .into_iter()
                            .filter(|byte| !matches!(byte, b'\n' | b'\r'))
                            .collect(),
                    ),
                    offset,
                }))
            }
        }
    }
}

/// Exact byte-level FastqTokenizer operator*/operator++ sequence used by the
/// three-string TSV schema. The increment runs even for its third field;
/// remaining tokens after a full schema are ignored by Source Table::push_back.
fn fastq_fields(raw: &[u8]) -> io::Result<Vec<Vec<u8>>> {
    let mut fields = Vec::new();
    let mut position = 0;
    while position < raw.len() && fields.len() < 3 {
        let marker = raw[position];
        let newline = raw[position..]
            .iter()
            .position(|&byte| byte == b'\n')
            .map_or(raw.len(), |offset| position + offset);
        let (start, end) = match marker {
            b'@' => (position + 1, newline),
            b'+' => {
                if newline == raw.len() {
                    (raw.len(), raw.len())
                } else {
                    (newline + 1, raw.len())
                }
            }
            _ => (position, newline),
        };
        let mut value = raw[start..end].to_vec();
        if value.last() == Some(&b'\r') {
            value.pop();
        }
        fields.push(value);
        match marker {
            b'+' => position = raw.len(),
            b'@' => position = newline.saturating_add(1),
            _ => {
                position = newline.saturating_add(1);
                if position < raw.len() && raw[position] != b'+' {
                    return Err(io::Error::other("Malformed FASTQ record"));
                }
            }
        }
    }
    if fields.len() < 3 {
        return Err(io::Error::other("Missing fields in input line"));
    }
    Ok(fields)
}

/// Shared Source schema boundary for the concrete prepared FASTQ reader.
pub(crate) fn source_fastq_fields(raw: &[u8]) -> io::Result<(Vec<u8>, Vec<u8>, Vec<u8>)> {
    let mut fields = fastq_fields(raw)?;
    let quality = fields.pop().unwrap();
    let sequence = fields.pop().unwrap();
    let id = fields.pop().unwrap();
    Ok((id, sequence, quality))
}
