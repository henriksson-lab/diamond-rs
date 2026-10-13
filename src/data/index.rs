//! Persistent seed-index construction.
//!
//! Translation of `diamond/src/data/index.cpp`.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::basic::reduction::Reduction;
use crate::basic::shape_config::ShapeConfig;
use crate::basic::value::SequenceType;
use crate::config::Sensitivity;
use crate::data::block::Block;
use crate::data::dmnd::ReferenceHeader;
use crate::data::dmnd_reader::for_each_dmnd_record;
use crate::data::seed_set::{HashedSeedSet, Table, SEED_INDEX_MAGIC_NUMBER, SEED_INDEX_VERSION};
use crate::search::sensitivity::{get_shape_codes, get_traits, soft_masking_algo};

/// C++ `makeindex` refuses larger databases because it materializes the whole
/// reference block before building the persistent hash tables.
pub const MAX_INDEX_LETTERS: u64 = 100_000_000;

#[derive(Debug, Clone)]
pub struct MakeIndexConfig {
    pub database: String,
    pub sensitivity: Sensitivity,
    pub shape_mask: Vec<String>,
    pub shapes: u32,
    /// Zero selects the platform's available parallelism, matching the CLI's
    /// upstream default.
    pub threads: usize,
}

impl MakeIndexConfig {
    pub fn new(database: impl Into<String>, sensitivity: Sensitivity) -> Self {
        Self {
            database: database.into(),
            sensitivity,
            shape_mask: Vec::new(),
            shapes: 0,
            threads: 0,
        }
    }
}

/// Resolve the database exactly like C++ `auto_append_extension_if_exists`:
/// retain an existing path, otherwise append (rather than replace) `.dmnd`.
fn resolve_database_path(database: &str) -> PathBuf {
    let path = Path::new(database);
    if path.exists() {
        return path.to_path_buf();
    }
    let mut appended = path.as_os_str().to_owned();
    appended.push(".dmnd");
    PathBuf::from(appended)
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
    if config.database.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Missing parameter: database file (--db/-d).",
        ));
    }

    let database_path = resolve_database_path(&config.database);
    let mut header_file = File::open(&database_path)?;
    let header = ReferenceHeader::read_from(&mut header_file)?;
    header
        .validate()
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    if header.letters > MAX_INDEX_LETTERS {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Indexing is only supported for databases of < 100000000 letters.",
        ));
    }

    let reduction = Reduction::default_reduction();
    let shape_codes: Vec<String> = if config.shape_mask.is_empty() {
        get_shape_codes(config.sensitivity)
            .iter()
            .map(|s| (*s).to_string())
            .collect()
    } else {
        config.shape_mask.clone()
    };
    let shapes = ShapeConfig::from_codes(&shape_codes, config.shapes, &reduction)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let mut block = Block::new();
    for_each_dmnd_record(&database_path, |oid, id, sequence| {
        // DMND stores encoded letters as raw bytes. `Letter` is `i8`, so this
        // is a same-size, alignment-one borrowed view valid for the callback.
        let sequence =
            unsafe { std::slice::from_raw_parts(sequence.as_ptr().cast::<i8>(), sequence.len()) };
        block
            .push_back(
                sequence,
                Some(id),
                None,
                oid as u64,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .map(|_| ())
            .map_err(io::Error::other)
    })?;

    let traits = get_traits(config.sensitivity);
    let soft_masking = soft_masking_algo(&traits, "", false, false)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let threads = if config.threads == 0 {
        std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1)
    } else {
        config.threads
    };
    block.soft_mask_with_threads(soft_masking, threads);
    let index = HashedSeedSet::from_block(&mut block, None, 0.0, soft_masking, &shapes, &reduction)
        .map_err(io::Error::other)?;

    let mut output_name = database_path.as_os_str().to_owned();
    output_name.push(".seed_idx");
    let output_path = PathBuf::from(output_name);
    let mut output = BufWriter::new(File::create(&output_path)?);
    write_seed_index(&mut output, &index, shapes.count() as usize)?;
    output.flush()?;
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
