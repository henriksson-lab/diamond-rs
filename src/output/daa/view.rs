//! Translation of `diamond/src/output/daa/view.cpp`.

use std::io::{self, Seek, Write};

use crate::data::daa::{
    finish_daa_from_file, init_daa, view_query_daa, view_query_edge, view_query_null,
    view_query_paf, view_query_pairwise, view_query_sam, view_query_tabular, view_query_xml,
    DaaFile, DaaQueryRecord, DaaViewConfig,
};
use crate::output::format::{print_tabular_footer, print_tabular_header, Header, TabularFormat};
use crate::output::{pairwise, sam, xml};
use crate::stats::score_matrix::ScoreMatrix;

pub const VIEW_BUF_SIZE: usize = 32;

/// Explicit replacement for the polymorphic C++ `OutputFormat` plus the
/// format-specific globals consumed by its implementations.
pub enum DaaViewFormat<'a> {
    Daa,
    Tabular(&'a TabularFormat),
    Paf,
    Sam { qlen_field: bool },
    Pairwise,
    Xml,
    Null,
    Edge,
}

/// Remaining process-global settings read by C++ `view_daa`.
#[derive(Debug, Clone, Copy)]
pub struct DaaViewRuntimeConfig<'a> {
    pub invocation: &'a str,
    pub tabular_header: Header,
}

impl Default for DaaViewRuntimeConfig<'_> {
    fn default() -> Self {
        Self {
            invocation: "",
            tabular_header: Header::None,
        }
    }
}

/// Rust counterpart of C++ `ViewWriter`.
pub struct ViewWriter<W> {
    writer: W,
}

impl<W: Write> ViewWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// C++ `ViewWriter::operator()(TextBuffer&)`.
    pub fn consume(&mut self, buf: &mut Vec<u8>) -> io::Result<()> {
        self.writer.write_all(buf)?;
        buf.clear();
        Ok(())
    }

    pub fn writer_mut(&mut self) -> &mut W {
        &mut self.writer
    }

    pub fn into_inner(self) -> W {
        self.writer
    }
}

/// Rust counterpart of C++ `ViewFetcher`.
#[derive(Debug, Clone)]
pub struct ViewFetcher {
    pub batch_size: usize,
    finished: bool,
}

impl Default for ViewFetcher {
    fn default() -> Self {
        Self {
            batch_size: VIEW_BUF_SIZE,
            finished: false,
        }
    }
}

impl ViewFetcher {
    pub fn new() -> Self {
        Self::default()
    }

    /// C++ `ViewFetcher::operator()`, returning the records and their stable
    /// query numbers directly rather than exposing mutable member buffers.
    pub fn fetch_batch(&mut self, daa: &mut DaaFile) -> io::Result<Vec<(Vec<u8>, usize)>> {
        if self.finished {
            return Ok(Vec::new());
        }
        let mut records = Vec::with_capacity(self.batch_size);
        for _ in 0..self.batch_size {
            match daa.read_query_buffer()? {
                Some(record) => records.push(record),
                None => {
                    self.finished = true;
                    break;
                }
            }
        }
        Ok(records)
    }
}

/// C++ `view_query` with its cloned output format and global configuration
/// represented explicitly.
pub fn view_query(
    record: &DaaQueryRecord,
    daa: &DaaFile,
    out: &mut Vec<u8>,
    format: &DaaViewFormat<'_>,
    score_matrix: &ScoreMatrix,
    config: &DaaViewConfig,
) -> Result<(), String> {
    match format {
        DaaViewFormat::Daa => view_query_daa(record, daa, out, score_matrix, config),
        DaaViewFormat::Tabular(format) => {
            // A DAA contains records only for aligned queries. C++ constructs
            // `Output::Info` with `unaligned=false` here, irrespective of the
            // command's report-unaligned setting.
            view_query_tabular(record, daa, out, score_matrix, format, config, false)
        }
        DaaViewFormat::Paf => view_query_paf(record, daa, out, score_matrix, config),
        DaaViewFormat::Sam { qlen_field } => {
            view_query_sam(record, daa, out, score_matrix, config, *qlen_field)
        }
        DaaViewFormat::Pairwise => view_query_pairwise(record, daa, out, score_matrix, config),
        DaaViewFormat::Xml => view_query_xml(record, daa, out, score_matrix, config),
        DaaViewFormat::Null => view_query_null(record, daa, score_matrix, config),
        DaaViewFormat::Edge => view_query_edge(record, daa, out, score_matrix, config),
    }
}

/// C++ `view_worker`. Rust keeps fetching and output ordering deterministic;
/// callers can schedule workers externally without exposing DAA reader state.
pub fn view_worker<W: Write>(
    daa: &mut DaaFile,
    writer: &mut ViewWriter<W>,
    fetcher: &mut ViewFetcher,
    format: &DaaViewFormat<'_>,
    score_matrix: &ScoreMatrix,
    config: &DaaViewConfig,
) -> io::Result<()> {
    loop {
        let batch = fetcher.fetch_batch(daa)?;
        if batch.is_empty() {
            return Ok(());
        }
        for (buf, query_num) in batch {
            let record = DaaQueryRecord::from_buffer(daa, &buf, query_num)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
            let mut output = Vec::new();
            view_query(&record, daa, &mut output, format, score_matrix, config)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
            // This is intentionally the literal C++ ViewWriter behavior.  In
            // particular, `view_daa` does not use the format's query separator
            // (the search join path does), so flat JSON has no comma between
            // query buffers in pinned upstream DIAMOND 2.1.24.
            writer.consume(&mut output)?;
        }
    }
}

/// C++ `view_daa`, parameterized over its input file, output stream, format,
/// and formerly global settings.
pub fn view_daa<W: Write + Seek>(
    daa: &mut DaaFile,
    output: &mut W,
    format: &DaaViewFormat<'_>,
    config: &DaaViewConfig,
    runtime: &DaaViewRuntimeConfig<'_>,
) -> io::Result<()> {
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;

    if matches!(format, DaaViewFormat::Daa) {
        init_daa(output)?;
    }

    let first = daa.read_query_buffer()?;
    let first_record = first
        .as_ref()
        .map(|(buf, query_num)| {
            DaaQueryRecord::from_buffer(daa, buf, *query_num)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
        })
        .transpose()?;

    write_header(output, daa, format, runtime, first_record.as_ref())?;

    let mut writer = ViewWriter::new(output);
    if let Some(record) = first_record {
        let mut buf = Vec::new();
        view_query(&record, daa, &mut buf, format, &score_matrix, config)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
        writer.consume(&mut buf)?;
    }
    let mut fetcher = ViewFetcher::new();
    view_worker(
        daa,
        &mut writer,
        &mut fetcher,
        format,
        &score_matrix,
        config,
    )?;

    match format {
        DaaViewFormat::Daa => finish_daa_from_file(writer.writer_mut(), daa)?,
        DaaViewFormat::Tabular(format) => writer
            .writer_mut()
            .write_all(print_tabular_footer(format.is_json).as_bytes())?,
        DaaViewFormat::Pairwise => pairwise::print_footer(writer.writer_mut())?,
        DaaViewFormat::Xml => xml::print_footer(writer.writer_mut())?,
        DaaViewFormat::Paf
        | DaaViewFormat::Sam { .. }
        | DaaViewFormat::Null
        | DaaViewFormat::Edge => {}
    }
    writer.writer_mut().flush()
}

fn write_header<W: Write>(
    writer: &mut W,
    daa: &DaaFile,
    format: &DaaViewFormat<'_>,
    runtime: &DaaViewRuntimeConfig<'_>,
    first: Option<&DaaQueryRecord>,
) -> io::Result<()> {
    match format {
        DaaViewFormat::Daa | DaaViewFormat::Paf | DaaViewFormat::Null | DaaViewFormat::Edge => {
            Ok(())
        }
        DaaViewFormat::Tabular(format) => {
            let header = print_tabular_header(
                &format.fields,
                runtime.tabular_header,
                format.is_json,
                "",
                "",
            )
            .map_err(io::Error::other)?;
            writer.write_all(header.as_bytes())
        }
        DaaViewFormat::Sam { .. } => sam::print_header(writer, daa.mode(), runtime.invocation),
        DaaViewFormat::Pairwise => pairwise::print_header(writer),
        DaaViewFormat::Xml => xml::print_header(
            writer,
            daa.mode(),
            &daa.score_matrix(),
            daa.gap_open_penalty(),
            daa.gap_extension_penalty(),
            daa.evalue(),
            first.map_or("", |record| record.query_name.as_str()),
            first.map_or(0, |record| record.query_len() as u32),
            "",
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_transcript::{EditOperation, PackedOperation};
    use crate::basic::value::SequenceType;
    use crate::data::daa::{compute_flag, write_daa_query_record, DaaHeader1, DaaHeader2};
    use std::fs::File;
    use std::io::Cursor;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static FILE_ID: AtomicUsize = AtomicUsize::new(0);

    fn input_daa() -> (DaaFile, std::path::PathBuf) {
        let id = FILE_ID.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-daa-view-{}-{id}.daa",
            std::process::id()
        ));
        let mut query = Vec::new();
        let seek = write_daa_query_record(
            &mut query,
            "query0 description",
            &[0, 1, 2],
            SequenceType::AminoAcid,
        );
        let flag = compute_flag(30, 0, 0, false);
        query.extend_from_slice(&0u32.to_ne_bytes());
        query.extend_from_slice(&[flag, 30, 0, 0]);
        query.push(PackedOperation::from_op_count(EditOperation::Match, 3).code);
        query.push(PackedOperation::terminator().code);
        crate::data::daa::finish_daa_query_record(&mut query, seek);
        query.extend_from_slice(&0u32.to_ne_bytes());

        let mut header =
            DaaHeader2::with_params(1, 100, 11, 1, 2, -3, 0.7, 1.2, 0.001, "BLOSUM62", 2);
        header.db_seqs_used = 1;
        header.query_records = 1;
        header.block_type[..3].copy_from_slice(&[1, 2, 3]);
        header.block_size[0] = query.len() as u64;
        header.block_size[1] = 5;
        header.block_size[2] = 4;

        let mut file = File::create(&path).unwrap();
        DaaHeader1::new().write_to(&mut file).unwrap();
        header.write_to(&mut file).unwrap();
        file.write_all(&query).unwrap();
        file.write_all(b"ref0\0").unwrap();
        file.write_all(&100u32.to_ne_bytes()).unwrap();
        drop(file);
        (DaaFile::open(&path).unwrap(), path)
    }

    #[test]
    fn writer_consumes_and_clears_buffers() {
        let mut writer = ViewWriter::new(Vec::new());
        let mut buf = b"record".to_vec();
        writer.consume(&mut buf).unwrap();
        assert!(buf.is_empty());
        assert_eq!(writer.into_inner(), b"record");
    }

    #[test]
    fn view_daa_round_trips_archive_metadata_and_records() {
        let (mut input, input_path) = input_daa();
        let mut output = Cursor::new(Vec::new());
        view_daa(
            &mut input,
            &mut output,
            &DaaViewFormat::Daa,
            &DaaViewConfig::default(),
            &DaaViewRuntimeConfig::default(),
        )
        .unwrap();

        let output_path = input_path.with_extension("out.daa");
        std::fs::write(&output_path, output.into_inner()).unwrap();
        let mut viewed = DaaFile::open(&output_path).unwrap();
        assert_eq!(viewed.db_seqs(), 1);
        assert_eq!(viewed.db_seqs_used(), 1);
        assert_eq!(viewed.ref_name(0), "ref0");
        let (buf, query_num) = viewed.read_query_buffer().unwrap().unwrap();
        let record = DaaQueryRecord::from_buffer(&viewed, &buf, query_num).unwrap();
        assert_eq!(record.query_name, "query0");
        assert!(viewed.read_query_buffer().unwrap().is_none());
        let _ = std::fs::remove_file(input_path);
        let _ = std::fs::remove_file(output_path);
    }

    #[test]
    fn tabular_view_is_exact_and_never_emits_unaligned_placeholder() {
        let (mut input, input_path) = input_daa();
        let fields = ["tab", "qseqid", "sseqid", "score"];
        let format = TabularFormat::new(&fields, false, 0, false, "BLOSUM62").unwrap();
        let mut output = Cursor::new(Vec::new());
        view_daa(
            &mut input,
            &mut output,
            &DaaViewFormat::Tabular(&format),
            &DaaViewConfig::default(),
            &DaaViewRuntimeConfig::default(),
        )
        .unwrap();
        assert_eq!(output.into_inner(), b"query0\tref0\t30\n");
        let _ = std::fs::remove_file(input_path);
    }

    #[test]
    fn fetcher_preserves_query_numbers_and_batch_limit() {
        let (mut input, input_path) = input_daa();
        let mut fetcher = ViewFetcher {
            batch_size: 1,
            finished: false,
        };
        let batch = fetcher.fetch_batch(&mut input).unwrap();
        assert_eq!(batch.len(), 1);
        assert_eq!(batch[0].1, 0);
        assert!(fetcher.fetch_batch(&mut input).unwrap().is_empty());
        assert!(fetcher.fetch_batch(&mut input).unwrap().is_empty());
        let _ = std::fs::remove_file(input_path);
    }
}
