//! Translation of diamond/src/run/tools.cpp.
//!
//! The upstream FASTA/FASTQ reads and pairwise alignment call are commented
//! out, leaving four literal while(true) loops. Rust makes those intended
//! dependencies explicit through ToolsBackend so the command bodies keep
//! their formatting and processing semantics without becoming non-terminating.

use std::collections::BTreeSet;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::basic::value::{Letter, AMINO_ACID_ALPHABET, MASK_LETTER, NUCLEOTIDE_ALPHABET};
use crate::data::fasta::{ByteFastaRecord, FastaRecord};
use crate::masking::tantan::TantanMasker;
use crate::stats::matrices::BLOSUM62;
use crate::util::sequence::seqid_bytes;

#[derive(Debug, Clone, PartialEq)]
pub struct ToolsConfig {
    pub database: Option<PathBuf>,
    pub query_file: Option<PathBuf>,
    pub seq_no: Vec<Vec<u8>>,
    pub output_file: Option<PathBuf>,
    pub oid_list: Option<PathBuf>,
    pub reverse: bool,
    pub hardmasked: bool,
    pub chunk_size: f64,
    pub tantan_min_mask_prob: f32,
    pub threads: i32,
}

impl Default for ToolsConfig {
    fn default() -> Self {
        Self {
            database: None,
            query_file: None,
            seq_no: Vec::new(),
            output_file: None,
            oid_list: None,
            reverse: false,
            hardmasked: false,
            chunk_size: 0.0,
            tantan_min_mask_prob: 0.9,
            threads: 0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToolSequenceType {
    AminoAcid,
    Nucleotide,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairwiseOperation {
    Substitution,
    Deletion,
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PairwiseColumn {
    pub operation: PairwiseOperation,
    pub subject_pos: i32,
    pub query_pos: i32,
    pub query_char: char,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PairwiseSettings {
    pub matrix_name: &'static str,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub reward: i32,
    pub penalty: i32,
}

impl Default for PairwiseSettings {
    fn default() -> Self {
        Self {
            matrix_name: "DNA",
            gap_open: 5,
            gap_extend: 2,
            reward: 0,
            penalty: 1,
        }
    }
}

pub trait ToolsBackend {
    fn get_seq(&mut self, config: &ToolsConfig) -> Result<(), String>;
    fn database_sequences(&mut self, database: &str) -> Result<Vec<FastaRecord>, String>;
    fn read_sequences(
        &mut self,
        path: &str,
        sequence_type: ToolSequenceType,
    ) -> Result<Vec<FastaRecord>, String>;
    /// Canonical Source byte records; provided bridges preserve existing String backends.
    fn database_sequences_bytes(&mut self, database: &str) -> Result<Vec<ByteFastaRecord>, String> {
        self.database_sequences(database)
            .map(|records| records.into_iter().map(Into::into).collect())
    }
    fn read_sequences_bytes(
        &mut self,
        path: &str,
        sequence_type: ToolSequenceType,
    ) -> Result<Vec<ByteFastaRecord>, String> {
        self.read_sequences(path, sequence_type)
            .map(|records| records.into_iter().map(Into::into).collect())
    }
    fn database_sequences_path_bytes(
        &mut self,
        database: &Path,
    ) -> Result<Vec<ByteFastaRecord>, String> {
        self.database_sequences_bytes(
            database
                .to_str()
                .ok_or("Path is not UTF-8; implement database_sequences_path_bytes")?,
        )
    }
    fn read_sequences_path_bytes(
        &mut self,
        path: &Path,
        sequence_type: ToolSequenceType,
    ) -> Result<Vec<ByteFastaRecord>, String> {
        self.read_sequences_bytes(
            path.to_str()
                .ok_or("Path is not UTF-8; implement read_sequences_path_bytes")?,
            sequence_type,
        )
    }
    fn random_index(&mut self, sequence_count: usize) -> usize;
    fn mask_sequence(&mut self, sequence: &mut [Letter], probability: f32) -> Result<(), String> {
        TantanMasker::new(&BLOSUM62, probability).mask_hard(sequence);
        Ok(())
    }
    /// Source initializes its score matrix before opening pairwise input,
    /// including with zero worker threads. Native backends must preserve any
    /// constructor error at that point rather than delaying it until a pair.
    fn prepare_pairwise(&mut self, _settings: PairwiseSettings) -> Result<(), String> {
        Ok(())
    }
    fn pairwise_alignment(
        &mut self,
        reference: &[Letter],
        query: &[Letter],
        settings: PairwiseSettings,
    ) -> Result<Vec<PairwiseColumn>, String>;
}

fn required<'a>(value: &'a Option<PathBuf>, name: &str) -> Result<&'a Path, String> {
    value
        .as_deref()
        .ok_or_else(|| format!("Missing required parameter: {name}"))
}

fn first_seq_no(config: &ToolsConfig) -> Result<usize, String> {
    let value = config
        .seq_no
        .first()
        .ok_or_else(|| "Missing required parameter: seq-no".to_owned())?;
    let mut value = value.clone();
    value.push(0);
    Ok(unsafe { libc::atoi(value.as_ptr().cast()) } as usize)
}

fn sequence_text(sequence: &[Letter], alphabet: &[u8]) -> String {
    sequence
        .iter()
        .map(|&letter| {
            let value = letter as u8;
            let index = (value & 127) as usize;
            let character = alphabet.get(index).copied().unwrap_or(b'X');
            if value & 128 == 0 {
                character as char
            } else {
                character.to_ascii_lowercase() as char
            }
        })
        .collect()
}

fn reversed_sequence_text(sequence: &[Letter], alphabet: &[u8]) -> String {
    sequence
        .iter()
        .rev()
        .map(|&letter| {
            alphabet
                .get((letter as u8 & 127) as usize)
                .copied()
                .unwrap_or(b'X') as char
        })
        .collect()
}

pub fn get_seq<B: ToolsBackend>(backend: &mut B, config: &ToolsConfig) -> Result<(), String> {
    required(&config.database, "database")?;
    backend.get_seq(config)
}

pub fn random_seqs<B: ToolsBackend, W: Write, S: Write>(
    backend: &mut B,
    config: &ToolsConfig,
    out: &mut W,
    status: &mut S,
) -> Result<(), String> {
    let records = backend.database_sequences_path_bytes(required(&config.database, "database")?)?;
    writeln!(status, "Sequences = {}", records.len()).map_err(|error| error.to_string())?;
    let count = first_seq_no(config)?;
    if count > records.len() {
        return Err("Requested sequence count exceeds database size.".to_owned());
    }
    let mut selected = BTreeSet::new();
    while selected.len() < count {
        selected.insert(backend.random_index(records.len()) % records.len());
    }
    for (number, index) in selected.into_iter().enumerate() {
        writeln!(out, ">{number}").map_err(|error| error.to_string())?;
        // Preserve the active C++ reverse branch: its actual reverse-printing
        // statement is commented out, so it emits an empty sequence line.
        if !config.reverse {
            writeln!(
                out,
                "{}",
                sequence_text(&records[index].sequence, AMINO_ACID_ALPHABET)
            )
            .map_err(|error| error.to_string())?;
        } else {
            writeln!(out).map_err(|error| error.to_string())?;
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct MaskerSummary {
    pub masked_sequences: usize,
    pub total_sequences: usize,
    pub masked_letters: usize,
}

pub fn run_masker_records<W: Write, E: Write>(
    records: &[FastaRecord],
    out: &mut W,
    err: &mut E,
    elapsed_seconds: f64,
    min_mask_prob: f32,
) -> Result<MaskerSummary, String> {
    let summary = mask_records(records, out, min_mask_prob)?;
    write_masker_summary(err, summary, elapsed_seconds)?;
    Ok(summary)
}

fn mask_records<W: Write, N: AsRef<[u8]>>(
    records: &[FastaRecord<N>],
    out: &mut W,
    min_mask_prob: f32,
) -> Result<MaskerSummary, String> {
    let masker = TantanMasker::new(&BLOSUM62, min_mask_prob);
    mask_records_with(records, out, |sequence| {
        masker.mask_hard(sequence);
        Ok(())
    })
}
fn mask_records_with<W: Write, N: AsRef<[u8]>>(
    records: &[FastaRecord<N>],
    out: &mut W,
    mut mask: impl FnMut(&mut [Letter]) -> Result<(), String>,
) -> Result<MaskerSummary, String> {
    let mut summary = MaskerSummary::default();
    for record in records {
        out.write_all(b">").map_err(|error| error.to_string())?;
        out.write_all(record.id.as_ref())
            .map_err(|error| error.to_string())?;
        writeln!(out).map_err(|error| error.to_string())?;
        let mut sequence = record.sequence.clone();
        mask(&mut sequence)?;
        writeln!(out, "{}", sequence_text(&sequence, AMINO_ACID_ALPHABET))
            .map_err(|error| error.to_string())?;
        let masked = sequence
            .iter()
            .filter(|&&letter| letter == MASK_LETTER)
            .count();
        summary.masked_letters += masked;
        summary.masked_sequences += usize::from(masked > 0);
        summary.total_sequences += 1;
    }
    Ok(summary)
}

fn write_masker_summary<E: Write>(
    err: &mut E,
    summary: MaskerSummary,
    elapsed_seconds: f64,
) -> Result<(), String> {
    writeln!(
        err,
        "#Sequences: {}/{}, #Letters: {}, t={}",
        summary.masked_sequences, summary.total_sequences, summary.masked_letters, elapsed_seconds
    )
    .map_err(|error| error.to_string())?;
    Ok(())
}

pub fn run_masker<B: ToolsBackend, W: Write, E: Write>(
    backend: &mut B,
    config: &ToolsConfig,
    out: &mut W,
    err: &mut E,
) -> Result<MaskerSummary, String> {
    let started = Instant::now();
    let records = backend.read_sequences_path_bytes(
        required(&config.query_file, "query")?,
        ToolSequenceType::AminoAcid,
    )?;
    let summary = mask_records_with(records.as_slice(), out, |sequence| {
        backend.mask_sequence(sequence, config.tantan_min_mask_prob)
    })?;
    // Source reports elapsed whole milliseconds after masking and writing,
    // rather than timing only the input read before entering the masker.
    write_masker_summary(err, summary, started.elapsed().as_millis() as f64 / 1000.0)?;
    Ok(summary)
}

pub fn fastq2fasta<B: ToolsBackend, W: Write>(
    backend: &mut B,
    config: &ToolsConfig,
    out: &mut W,
) -> Result<(), String> {
    let records = backend.read_sequences_path_bytes(
        required(&config.query_file, "query")?,
        ToolSequenceType::Nucleotide,
    )?;
    let max = first_seq_no(config)?;
    for record in records.iter().take(max) {
        out.write_all(b">").map_err(|error| error.to_string())?;
        out.write_all(record.id.as_ref())
            .map_err(|error| error.to_string())?;
        writeln!(out).map_err(|error| error.to_string())?;
        writeln!(
            out,
            "{}",
            sequence_text(&record.sequence, NUCLEOTIDE_ALPHABET)
        )
        .map_err(|error| error.to_string())?;
    }
    Ok(())
}

pub fn architecture_flags() -> Vec<&'static str> {
    #[allow(unused_mut)]
    let mut flags = Vec::new();
    #[cfg(target_feature = "sse2")]
    flags.push("sse2");
    #[cfg(target_feature = "sse3")]
    flags.push("sse3");
    #[cfg(target_feature = "ssse3")]
    flags.push("ssse3");
    #[cfg(target_feature = "popcnt")]
    flags.push("popcnt");
    #[cfg(target_feature = "neon")]
    flags.push("neon");
    flags
}

pub fn info<W: Write>(out: &mut W) -> Result<(), String> {
    write!(out, "Architecture flags: ").map_err(|error| error.to_string())?;
    for flag in architecture_flags() {
        write!(out, "{flag} ").map_err(|error| error.to_string())?;
    }
    writeln!(out).map_err(|error| error.to_string())
}

pub fn pairwise_worker<B: ToolsBackend, W: Write>(
    backend: &mut B,
    records: &[FastaRecord],
    out: &mut W,
) -> Result<(), String> {
    pairwise_worker_raw(backend, records, out)
}

fn pairwise_worker_raw<B: ToolsBackend, W: Write, N: AsRef<[u8]>>(
    backend: &mut B,
    records: &[FastaRecord<N>],
    out: &mut W,
) -> Result<(), String> {
    for pair in records.chunks_exact(2) {
        let reference = &pair[0];
        let query = &pair[1];
        let reference_id = seqid_bytes(reference.id.as_ref());
        let query_id = seqid_bytes(query.id.as_ref());
        let columns = backend.pairwise_alignment(
            &reference.sequence,
            &query.sequence,
            PairwiseSettings::default(),
        )?;
        let mut buffer = Vec::new();
        for column in columns {
            if matches!(
                column.operation,
                PairwiseOperation::Substitution | PairwiseOperation::Deletion
            ) {
                buffer
                    .write_all(&reference_id)
                    .map_err(|error| error.to_string())?;
                buffer.write_all(b"\t").map_err(|error| error.to_string())?;
                buffer
                    .write_all(&query_id)
                    .map_err(|error| error.to_string())?;
                match column.operation {
                    PairwiseOperation::Substitution => writeln!(
                        buffer,
                        "\t{}\t{}\t{}",
                        column.subject_pos, column.query_pos, column.query_char
                    ),
                    PairwiseOperation::Deletion => {
                        writeln!(buffer, "\t{}\t-1\t-", column.subject_pos)
                    }
                    PairwiseOperation::Other => unreachable!(),
                }
                .map_err(|error| error.to_string())?;
            }
        }
        out.write_all(&buffer).map_err(|error| error.to_string())?;
    }
    Ok(())
}

pub fn pairwise<B: ToolsBackend, W: Write>(
    backend: &mut B,
    config: &ToolsConfig,
    out: &mut W,
) -> Result<(), String> {
    backend.prepare_pairwise(PairwiseSettings::default())?;
    let records = backend.read_sequences_path_bytes(
        required(&config.query_file, "query")?,
        ToolSequenceType::Nucleotide,
    )?;
    if config.threads <= 0 {
        return Ok(());
    }
    // C++ workers share the input stream and serialize complete per-pair
    // buffers. Sequential ownership preserves those record boundaries; only
    // the unspecified cross-thread ordering is made deterministic.
    pairwise_worker_raw(backend, &records, out)
}

pub fn reverse<B: ToolsBackend, W: Write>(
    backend: &mut B,
    config: &ToolsConfig,
    out: &mut W,
) -> Result<(), String> {
    let records = backend.read_sequences_path_bytes(
        required(&config.query_file, "query")?,
        ToolSequenceType::AminoAcid,
    )?;
    for record in records {
        out.write_all(b">\\").map_err(|error| error.to_string())?;
        out.write_all(&record.id)
            .map_err(|error| error.to_string())?;
        writeln!(out).map_err(|error| error.to_string())?;
        writeln!(
            out,
            "{}",
            reversed_sequence_text(&record.sequence, AMINO_ACID_ALPHABET)
        )
        .map_err(|error| error.to_string())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Default)]
    struct Backend {
        database: Vec<FastaRecord>,
        input: Vec<FastaRecord>,
        random: Vec<usize>,
        alignments: Vec<Vec<PairwiseColumn>>,
        got_get_seq: Option<PathBuf>,
        sequence_types: Vec<ToolSequenceType>,
    }

    impl ToolsBackend for Backend {
        fn get_seq(&mut self, config: &ToolsConfig) -> Result<(), String> {
            self.got_get_seq = config.database.clone();
            Ok(())
        }
        fn database_sequences(&mut self, _: &str) -> Result<Vec<FastaRecord>, String> {
            Ok(self.database.clone())
        }
        fn read_sequences(
            &mut self,
            _: &str,
            kind: ToolSequenceType,
        ) -> Result<Vec<FastaRecord>, String> {
            self.sequence_types.push(kind);
            Ok(self.input.clone())
        }
        fn random_index(&mut self, sequence_count: usize) -> usize {
            self.random.remove(0) % sequence_count
        }
        fn pairwise_alignment(
            &mut self,
            _: &[Letter],
            _: &[Letter],
            settings: PairwiseSettings,
        ) -> Result<Vec<PairwiseColumn>, String> {
            assert_eq!(settings, PairwiseSettings::default());
            Ok(self.alignments.remove(0))
        }
    }

    #[test]
    fn byte_backend_reaches_tool_headers_without_checked_string_reads() {
        struct RawBackend;
        impl ToolsBackend for RawBackend {
            fn get_seq(&mut self, _: &ToolsConfig) -> Result<(), String> {
                Ok(())
            }
            fn database_sequences(&mut self, _: &str) -> Result<Vec<FastaRecord>, String> {
                Err("String reader must not be called".into())
            }
            fn read_sequences(
                &mut self,
                _: &str,
                _: ToolSequenceType,
            ) -> Result<Vec<FastaRecord>, String> {
                Err("String reader must not be called".into())
            }
            fn read_sequences_bytes(
                &mut self,
                _: &str,
                _: ToolSequenceType,
            ) -> Result<Vec<ByteFastaRecord>, String> {
                Ok(vec![ByteFastaRecord {
                    id: b"raw\xff\0tail".to_vec(),
                    sequence: vec![0, 1],
                    quality: None,
                }])
            }
            fn random_index(&mut self, _: usize) -> usize {
                0
            }
            fn pairwise_alignment(
                &mut self,
                _: &[Letter],
                _: &[Letter],
                _: PairwiseSettings,
            ) -> Result<Vec<PairwiseColumn>, String> {
                Ok(vec![])
            }
        }
        let config = ToolsConfig {
            query_file: Some("unused".into()),
            seq_no: vec!["1".into()],
            ..Default::default()
        };
        let mut out = Vec::new();
        fastq2fasta(&mut RawBackend, &config, &mut out).unwrap();
        assert_eq!(out, b">raw\xff\0tail\nAC\n");
        out.clear();
        reverse(&mut RawBackend, &config, &mut out).unwrap();
        assert_eq!(out, b">\\raw\xff\0tail\nRA\n");
    }

    fn record(id: &str, sequence: &[Letter]) -> FastaRecord {
        FastaRecord {
            quality: None,
            id: id.to_owned(),
            sequence: sequence.to_vec(),
        }
    }

    #[test]
    fn get_seq_and_random_sequences_preserve_headers_order_and_reverse_quirk() {
        let mut backend = Backend {
            database: vec![record("a", &[0, 1]), record("b", &[2]), record("c", &[3])],
            random: vec![2, 0],
            ..Default::default()
        };
        let mut config = ToolsConfig {
            database: Some("db.dmnd".into()),
            seq_no: vec!["2junk".into()],
            ..Default::default()
        };
        get_seq(&mut backend, &config).unwrap();
        assert_eq!(backend.got_get_seq.as_deref(), Some(Path::new("db.dmnd")));
        let mut out = Vec::new();
        let mut status = Vec::new();
        random_seqs(&mut backend, &config, &mut out, &mut status).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">0\nAR\n>1\nD\n");
        assert_eq!(String::from_utf8(status).unwrap(), "Sequences = 3\n");
        backend.random = vec![1];
        config.seq_no = vec!["1".into()];
        config.reverse = true;
        let mut out = Vec::new();
        let mut status = Vec::new();
        random_seqs(&mut backend, &config, &mut out, &mut status).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">0\n\n");
    }

    #[test]
    fn masker_fastq_conversion_info_and_reverse_formats() {
        let records = vec![record("rep", &[0; 80])];
        let mut out = Vec::new();
        let mut err = Vec::new();
        let summary = run_masker_records(&records, &mut out, &mut err, 1.25, 0.9).unwrap();
        assert_eq!(summary.total_sequences, 1);
        assert!(String::from_utf8(out).unwrap().starts_with(">rep\n"));
        assert!(String::from_utf8(err).unwrap().ends_with("t=1.25\n"));

        let mut backend = Backend {
            input: vec![record("q", &[0, 1, 4])],
            ..Default::default()
        };
        let config = ToolsConfig {
            query_file: Some("q.fq".into()),
            seq_no: vec!["1".into()],
            ..Default::default()
        };
        let mut out = Vec::new();
        fastq2fasta(&mut backend, &config, &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">q\nACN\n");
        let mut out = Vec::new();
        reverse(&mut backend, &config, &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">\\q\nCRA\n");
        let mut out = Vec::new();
        info(&mut out).unwrap();
        assert!(String::from_utf8(out)
            .unwrap()
            .starts_with("Architecture flags: "));
    }

    #[test]
    fn pairwise_formats_only_substitutions_and_deletions() {
        let mut backend = Backend {
            input: vec![record("ref desc", &[0]), record("query desc", &[1])],
            alignments: vec![vec![
                PairwiseColumn {
                    operation: PairwiseOperation::Substitution,
                    subject_pos: 4,
                    query_pos: 7,
                    query_char: 'C',
                },
                PairwiseColumn {
                    operation: PairwiseOperation::Deletion,
                    subject_pos: 5,
                    query_pos: 8,
                    query_char: 'A',
                },
                PairwiseColumn {
                    operation: PairwiseOperation::Other,
                    subject_pos: 6,
                    query_pos: 9,
                    query_char: 'T',
                },
            ]],
            ..Default::default()
        };
        let config = ToolsConfig {
            query_file: Some("pairs.fa".into()),
            threads: 4,
            ..Default::default()
        };
        let mut out = Vec::new();
        pairwise(&mut backend, &config, &mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "ref\tquery\t4\t7\tC\nref\tquery\t5\t-1\t-\n"
        );
        assert_eq!(backend.sequence_types, vec![ToolSequenceType::Nucleotide]);
    }

    #[test]
    fn sequence_printing_matches_soft_mask_and_reverse_rules() {
        let sequence = [0, crate::basic::value::SEED_MASK | 1, 4];
        assert_eq!(sequence_text(&sequence, AMINO_ACID_ALPHABET), "ArC");
        assert_eq!(
            reversed_sequence_text(&sequence, AMINO_ACID_ALPHABET),
            "CRA"
        );
    }
}
