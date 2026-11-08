//! Persistent seed-index construction.
//!
//! Translation of `diamond/src/data/index.cpp`.

use std::io::{self, Write};
use std::path::{Path, PathBuf};

use crate::basic::reduction::Reduction;
use crate::basic::shape_config::ShapeConfig;
use crate::basic::value::SequenceType;
use crate::config::Sensitivity;
use crate::data::block::Block;
use crate::data::dmnd::{DatabaseOpenOptions, SourceDatabaseFile};
use crate::data::seed_set::{HashedSeedSet, Table, SEED_INDEX_MAGIC_NUMBER, SEED_INDEX_VERSION};
use crate::search::sensitivity::{get_shape_codes, get_traits, soft_masking_algo};

/// C++ `makeindex` refuses larger databases because it materializes the whole
/// reference block before building the persistent hash tables.
pub const MAX_INDEX_LETTERS: u64 = 100_000_000;

#[derive(Debug, Clone)]
pub struct MakeIndexConfig {
    pub database: PathBuf,
    pub sensitivity: Sensitivity,
    pub shape_mask: Vec<Vec<u8>>,
    pub shapes: u32,
    /// Zero selects Source's hardware concurrency, matching the CLI's
    /// upstream default.
    pub threads: usize,
    pub file_buffer_size: usize,
}

impl MakeIndexConfig {
    pub fn new(database: impl Into<PathBuf>, sensitivity: Sensitivity) -> Self {
        Self {
            database: database.into(),
            sensitivity,
            shape_mask: Vec::new(),
            shapes: 0,
            threads: 0,
            file_buffer_size: 67_108_864,
        }
    }
}

/// Match DatabaseFile's actual filename fallback: append only when the
/// suffixed file exists, retaining an absent caller path for its open error.
fn resolve_database_path(database: &Path) -> PathBuf {
    let mut suffixed = database.as_os_str().to_owned();
    suffixed.push(".dmnd");
    let suffixed = PathBuf::from(suffixed);
    if !database.as_os_str().as_encoded_bytes().ends_with(b".dmnd") && suffixed.exists() {
        suffixed
    } else {
        database.to_owned()
    }
}

/// Serialize the on-disk format emitted by C++ `makeindex`.
pub fn write_seed_index<W: Write>(
    mut writer: W,
    index: &HashedSeedSet,
    shape_count: usize,
) -> io::Result<()> {
    writer.write_all(&SEED_INDEX_MAGIC_NUMBER.to_ne_bytes())?;
    writer.write_all(&SEED_INDEX_VERSION.to_ne_bytes())?;
    writer.write_all(&(shape_count as u32).to_ne_bytes())?;

    for i in 0..shape_count {
        writer.write_all(&index.table(i).size().to_ne_bytes())?;
    }
    for i in 0..shape_count {
        let table = index.table(i);
        writer.write_all(&table.data()[..table.size() + Table::PADDING])?;
    }
    Ok(())
}

/// Build `<resolved database path>.seed_idx`.
///
/// This is the Rust counterpart of the single C++ `makeindex()` function in
/// `data/index.cpp`; configuration that was global in C++ is explicit here.
pub fn make_index(config: &MakeIndexConfig) -> io::Result<PathBuf> {
    if config.database.as_os_str().is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Missing parameter: database file (--db/-d).",
        ));
    }

    let database_path = resolve_database_path(&config.database);
    let mut database = SourceDatabaseFile::open(
        &database_path,
        DatabaseOpenOptions {
            file_buffer_size: Some(config.file_buffer_size),
            ..Default::default()
        },
    )?;
    let header = &database.ref_header;
    if header.letters > MAX_INDEX_LETTERS {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Indexing is only supported for databases of < 100000000 letters.",
        ));
    }

    let reduction = Reduction::default_reduction();
    let shape_codes: Vec<Vec<u8>> = if config.shape_mask.is_empty() {
        get_shape_codes(config.sensitivity)
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect()
    } else {
        config.shape_mask.clone()
    };
    let shapes = ShapeConfig::from_codes_bytes(&shape_codes, config.shapes, &reduction)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let mut block = Block::new();
    database.visit_sequences(MAX_INDEX_LETTERS, |oid, sequence| {
        block
            .push_back(sequence, None, None, oid, SequenceType::AminoAcid, 0, false)
            .map(|_| ())
            .map_err(io::Error::other)
    })?;

    // Source makeindex -> DatabaseFile::load_seqs -> load_twopass clears
    // persisted DMND SEG annotations before HashedSeedSet applies the selected
    // sensitivity's soft masking (sequence_file.cpp:267). Raw indexed readers
    // intentionally retain those bytes; this is the concrete loader boundary.
    for sequence_id in 0..block.seqs().len() {
        crate::masking::remove_bit_mask(block.seqs_mut().get_mut(sequence_id));
    }

    let traits = get_traits(config.sensitivity);
    let soft_masking = soft_masking_algo(&traits, "", false, false)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let threads = if config.threads == 0 {
        crate::util::system::hardware_concurrency() as usize
    } else {
        config.threads
    };
    block.soft_mask_with_threads(soft_masking, threads);
    let index = HashedSeedSet::from_block(&mut block, None, 0.0, soft_masking, &shapes, &reduction)
        .map_err(io::Error::other)?;

    let mut output_name = database_path.as_os_str().to_owned();
    output_name.push(".seed_idx");
    let output_path = PathBuf::from(output_name);
    let mut output = crate::util::io::OutputFile::with_source_profile_path(
        &output_path,
        crate::util::io::Compressor::None,
        "wb",
        config.file_buffer_size,
    )
    .map_err(io::Error::other)?;
    struct SourceWriter<'a>(&'a mut crate::util::io::OutputFile);
    impl Write for SourceWriter<'_> {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.0
                .write_raw(bytes)
                .map(|_| bytes.len())
                .map_err(io::Error::other)
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    write_seed_index(SourceWriter(&mut output), &index, shapes.count() as usize)?;
    output.close().map_err(io::Error::other)?;
    Ok(output_path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_missing_database_parameter() {
        let err = make_index(&MakeIndexConfig::new("", Sensitivity::Default)).unwrap_err();
        assert_eq!(
            err.to_string(),
            "Missing parameter: database file (--db/-d)."
        );
    }
}
