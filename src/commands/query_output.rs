//! Streaming auxiliary query output from `data/queries.cpp` and
//! `run/double_indexed.cpp`. Classification belongs to the completed query
//! block: callers must supply the union of actual alignment results across
//! reference blocks and iterative rounds, before dropping its query owner.

use std::io;

use crate::basic::value::{
    Letter, SequenceType, ValueTraits, AMINO_ACID_ALPHABET, MASK_LETTER, NUCLEOTIDE_ALPHABET,
};
use crate::util::io::{Compressor, OutputFile};
use crate::util::sequence;
use crate::util::text_buffer::TextBuffer;

/// Source `unaligned`, `aligned_file`, `unfmt`, and `alfmt` options.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AuxiliaryOutputOptions {
    pub unaligned: std::path::PathBuf,
    pub aligned: std::path::PathBuf,
    pub unaligned_format: Vec<u8>,
    pub aligned_format: Vec<u8>,
}

impl Default for AuxiliaryOutputOptions {
    fn default() -> Self {
        Self {
            unaligned: std::path::PathBuf::new(),
            aligned: std::path::PathBuf::new(),
            unaligned_format: b"fasta".to_vec(),
            aligned_format: b"fasta".to_vec(),
        }
    }
}

impl AuxiliaryOutputOptions {
    pub fn enabled(&self) -> bool {
        !self.unaligned.as_os_str().is_empty() || !self.aligned.as_os_str().is_empty()
    }

    /// Source retains qualities whenever either format is FASTQ, including
    /// when that particular output path is disabled.
    pub fn needs_quality(&self) -> bool {
        self.unaligned_format == b"fastq" || self.aligned_format == b"fastq"
    }
}

/// One pair of streams for an entire search, not one pair per round/block.
/// Paths are used literally and outputs are uncompressed, as in Source;
/// the main search compressor and a `.gz` suffix do not change this.
pub(crate) struct QueryOutputs {
    options: AuxiliaryOutputOptions,
    unaligned: Option<OutputFile>,
    aligned: Option<OutputFile>,
    input_traits: ValueTraits,
    buffer: TextBuffer,
}

impl QueryOutputs {
    pub(crate) fn new(
        options: &AuxiliaryOutputOptions,
        input_type: SequenceType,
        file_buffer_size: usize,
    ) -> io::Result<Self> {
        let open = |path: &std::path::Path| {
            if path.as_os_str().is_empty() {
                Ok(None)
            } else {
                OutputFile::with_source_profile_path(path, Compressor::None, "wb", file_buffer_size)
                    .map(Some)
                    .map_err(io::Error::other)
            }
        };
        // Preserve Source's creation order, including when both paths name
        // the same file. Formats are intentionally validated only on write.
        let unaligned = open(&options.unaligned)?;
        let aligned = open(&options.aligned)?;
        let input_traits = match input_type {
            SequenceType::AminoAcid => {
                ValueTraits::new(AMINO_ACID_ALPHABET, MASK_LETTER, b"UO-", input_type)
            }
            SequenceType::Nucleotide => ValueTraits::new(NUCLEOTIDE_ALPHABET, 4, b"", input_type),
        };
        Ok(Self {
            options: options.clone(),
            unaligned,
            aligned,
            input_traits,
            buffer: TextBuffer::new(),
        })
    }

    /// Emit one source query in the completed block's order. For translated
    /// queries, `source` is original DNA, never translated/masked frames.
    /// Protein queries use the prepared/restored source Block sequence,
    /// retaining Source's hard masking. IDs and qualities come from that
    /// same query owner; paired mates are classified independently.
    pub(crate) fn write_query(
        &mut self,
        aligned: bool,
        title: &str,
        source: &[Letter],
        quality: Option<&[u8]>,
    ) -> io::Result<()> {
        self.write_query_bytes(aligned, title.as_bytes(), source, quality)
    }

    pub(crate) fn write_query_bytes(
        &mut self,
        aligned: bool,
        title: &[u8],
        source: &[Letter],
        quality: Option<&[u8]>,
    ) -> io::Result<()> {
        let (output, format) = if aligned {
            (&mut self.aligned, &self.options.aligned_format)
        } else {
            (&mut self.unaligned, &self.options.unaligned_format)
        };
        let Some(output) = output else {
            return Ok(());
        };
        self.buffer.clear();
        // FASTQ with absent qualities dereferences a null C string in
        // Source. Define that invalid-input case safely as an empty quality
        // line, matching the translated sequence formatter, without making
        // up quality values. Valid FASTQ input retains every original byte.
        sequence::format_bytes_with_format(
            source,
            title,
            quality,
            &mut self.buffer,
            format,
            &self.input_traits,
            160,
        )
        .map_err(io::Error::other)?;
        output
            .write_raw(self.buffer.data())
            .map_err(io::Error::other)
    }

    /// Call after finalizing the main search writer, as in Source.
    pub(crate) fn finish(mut self) -> io::Result<()> {
        if let Some(file) = &mut self.unaligned {
            file.close().map_err(io::Error::other)?;
        }
        if let Some(file) = &mut self.aligned {
            file.close().map_err(io::Error::other)?;
        }
        Ok(())
    }
}
