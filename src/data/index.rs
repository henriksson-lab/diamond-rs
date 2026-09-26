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
use crate::data::dmnd_reader::read_dmnd;
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
}

impl MakeIndexConfig {
    pub fn new(database: impl Into<String>, sensitivity: Sensitivity) -> Self {
        Self {
            database: database.into(),
            sensitivity,
            shape_mask: Vec::new(),
            shapes: 0,
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
    let (header, records) = read_dmnd(&database_path)?;
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
    for (oid, record) in records.iter().enumerate() {
        block
            .push_back(
                &record.sequence,
                Some(&record.id),
                None,
                oid as u64,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .map_err(io::Error::other)?;
    }

    let traits = get_traits(config.sensitivity);
    let soft_masking = soft_masking_algo(&traits, "", false, false)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
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

    #[test]
    fn writes_an_index_readable_by_hashed_seed_set() {
        let database = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/data.dmnd");
        let output = make_index(&MakeIndexConfig::new(database, Sensitivity::Default)).unwrap();

        let reduction = Reduction::default_reduction();
        let codes = get_shape_codes(Sensitivity::Default)
            .iter()
            .map(|s| (*s).to_string())
            .collect::<Vec<_>>();
        let shapes = ShapeConfig::from_codes(&codes, 0, &reduction).unwrap();
        let parsed = HashedSeedSet::from_index_file(&output, &shapes).unwrap();
        assert_eq!(parsed.table(0).size() > 0, true);
        assert_eq!(parsed.table(1).size() > 0, true);

        std::fs::remove_file(output).unwrap();
    }
}
