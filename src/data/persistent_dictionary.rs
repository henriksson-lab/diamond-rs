//! Source `SequenceFile` dictionary persistence for cross-process block joins.
//!
//! On-disk entries follow `data/sequence_file.cpp::write_dict_entry`: native
//! u32 OID, optional native u32 length, optional nul-terminated title, optional
//! sequence bytes, and optional native f64 self score. IDs are implicit entry
//! ordinals, reset for each multiprocessing (query, reference) dictionary.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use super::sequence_file::{ByteDictionaryEntry, DictionaryEntry, DICT_EMPTY};
use crate::basic::value::DictId;

#[derive(Clone, Copy, Debug, Default)]
pub struct DictionaryLayout {
    pub lengths: bool,
    pub titles: bool,
    pub sequences: bool,
    pub self_scores: bool,
}

pub fn dictionary_file_name(directory: &Path, query: usize, reference: usize) -> PathBuf {
    directory.join(format!("ref_dict_{query}_{reference}"))
}

pub struct DictionaryWriter {
    output: BufWriter<File>,
    layout: DictionaryLayout,
    block_ids: HashMap<usize, Vec<DictId>>,
    next_id: DictId,
}
impl DictionaryWriter {
    pub fn create(path: &Path, layout: DictionaryLayout, buffer_size: usize) -> io::Result<Self> {
        if layout.sequences && !layout.lengths {
            return Err(io::Error::other(
                "Sequence dictionaries require stored lengths.",
            ));
        }
        Ok(Self {
            output: BufWriter::with_capacity(buffer_size, File::create(path)?),
            layout,
            block_ids: HashMap::new(),
            next_id: 0,
        })
    }

    pub fn init_block(&mut self, block: usize, count: usize, persist: bool) {
        if !persist {
            self.block_ids.clear();
        }
        self.block_ids
            .entry(block)
            .or_insert_with(|| vec![DICT_EMPTY; count]);
    }

    /// Assign each loaded target once, writing its record at that ID's ordinal.
    pub fn dict_id(
        &mut self,
        block: usize,
        target: usize,
        entry: &DictionaryEntry,
    ) -> io::Result<DictId> {
        self.dict_id_impl(block, target, entry)
    }
    pub fn dict_id_bytes(
        &mut self,
        block: usize,
        target: usize,
        entry: &ByteDictionaryEntry,
    ) -> io::Result<DictId> {
        self.dict_id_impl(block, target, entry)
    }
    fn dict_id_impl<T: AsRef<[u8]>>(
        &mut self,
        block: usize,
        target: usize,
        entry: &DictionaryEntry<T>,
    ) -> io::Result<DictId> {
        let ids = self
            .block_ids
            .get_mut(&block)
            .ok_or_else(|| io::Error::other("Dictionary not initialized."))?;
        let id = ids
            .get_mut(target)
            .ok_or_else(|| io::Error::other("Dictionary not initialized."))?;
        if *id != DICT_EMPTY {
            return Ok(*id);
        }
        self.output.write_all(&(entry.oid as u32).to_ne_bytes())?;
        if self.layout.lengths {
            self.output.write_all(&(entry.len as u32).to_ne_bytes())?;
        }
        if self.layout.titles {
            // Source writes a C string, stopping at the first nul.
            let title = entry.title.as_ref();
            let end = title
                .iter()
                .position(|&byte| byte == 0)
                .unwrap_or(title.len());
            self.output.write_all(&title[..end])?;
            self.output.write_all(&[0])?;
        }
        if self.layout.sequences {
            let length = usize::try_from(entry.len)
                .map_err(|_| io::Error::other("Negative dictionary sequence length."))?;
            if length != entry.seq.len() {
                return Err(io::Error::other("Dictionary sequence length mismatch."));
            }
            for &letter in &entry.seq {
                self.output.write_all(&[letter as u8])?;
            }
        }
        if self.layout.self_scores {
            self.output.write_all(&entry.self_aln_score.to_ne_bytes())?;
        }
        *id = self.next_id;
        self.next_id = self
            .next_id
            .checked_add(1)
            .ok_or_else(|| io::Error::other("Dictionary ID overflow."))?;
        Ok(*id)
    }

    pub fn close_block(&mut self, persist: bool) -> io::Result<()> {
        self.output.flush()?;
        if !persist {
            self.block_ids.clear();
        }
        Ok(())
    }

    pub fn finish(mut self) -> io::Result<()> {
        self.output.flush()
    }
}

pub fn load_dictionary(
    path: &Path,
    layout: DictionaryLayout,
    buffer_size: usize,
) -> io::Result<Vec<DictionaryEntry>> {
    load_dictionary_bytes(path, layout, buffer_size)?
        .into_iter()
        .map(|entry| {
            DictionaryEntry::try_from(entry)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
        })
        .collect()
}

pub fn load_dictionary_bytes(
    path: &Path,
    layout: DictionaryLayout,
    buffer_size: usize,
) -> io::Result<Vec<ByteDictionaryEntry>> {
    if layout.sequences && !layout.lengths {
        return Err(io::Error::other(
            "Sequence dictionaries require stored lengths.",
        ));
    }
    let mut input = BufReader::with_capacity(buffer_size, File::open(path)?);
    let mut entries = Vec::new();
    loop {
        let mut oid = [0; 4];
        if input.read(&mut oid[..1])? == 0 {
            break;
        }
        input.read_exact(&mut oid[1..])?;
        let length = if layout.lengths {
            let mut value = [0; 4];
            input.read_exact(&mut value)?;
            u32::from_ne_bytes(value) as i32
        } else {
            0
        };
        let title = if layout.titles {
            let mut value = Vec::new();
            if input.read_until(0, &mut value)? == 0 || value.last() != Some(&0) {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Truncated dictionary title.",
                ));
            }
            value.pop();
            value
        } else {
            Vec::new()
        };
        let seq = if layout.sequences {
            let count = usize::try_from(length)
                .map_err(|_| io::Error::other("Negative dictionary sequence length."))?;
            let mut bytes = vec![0; count];
            input.read_exact(&mut bytes)?;
            bytes.into_iter().map(|letter| letter as i8).collect()
        } else {
            Vec::new()
        };
        let self_aln_score = if layout.self_scores {
            let mut bytes = [0; 8];
            input.read_exact(&mut bytes)?;
            f64::from_ne_bytes(bytes)
        } else {
            0.0
        };
        entries.push(DictionaryEntry {
            oid: u32::from_ne_bytes(oid) as u64,
            len: length,
            title,
            seq,
            self_aln_score,
        });
    }
    Ok(entries)
}

/// Compatibility spelling for callers describing the load as a read.
pub fn read_dictionary_bytes(
    path: &Path,
    layout: DictionaryLayout,
    buffer_size: usize,
) -> io::Result<Vec<ByteDictionaryEntry>> {
    load_dictionary_bytes(path, layout, buffer_size)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn raw_title_dictionary_preserves_bytes_and_source_c_string_extent() {
        let path = std::env::temp_dir().join(format!(
            "diamond-byte-dictionary-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ));
        let layout = DictionaryLayout {
            lengths: true,
            titles: true,
            sequences: true,
            self_scores: true,
        };
        let entry = ByteDictionaryEntry {
            oid: 7,
            len: 3,
            title: b"raw\xff\0ignored".to_vec(),
            seq: vec![0, 1, 2],
            self_aln_score: 12.5,
        };
        let mut writer = DictionaryWriter::create(&path, layout, 128).unwrap();
        writer.init_block(0, 1, false);
        assert_eq!(writer.dict_id_bytes(0, 0, &entry).unwrap(), 0);
        assert_eq!(writer.dict_id_bytes(0, 0, &entry).unwrap(), 0);
        writer.finish().unwrap();
        let mut expected = 7u32.to_ne_bytes().to_vec();
        expected.extend_from_slice(&3u32.to_ne_bytes());
        expected.extend_from_slice(b"raw\xff\0");
        expected.extend_from_slice(&[0, 1, 2]);
        expected.extend_from_slice(&12.5f64.to_ne_bytes());
        assert_eq!(std::fs::read(&path).unwrap(), expected);
        let loaded = load_dictionary_bytes(&path, layout, 128).unwrap();
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded[0].title, b"raw\xff");
        assert_eq!(loaded[0].seq, entry.seq);
        assert_eq!(loaded[0].self_aln_score, 12.5);
        assert_eq!(
            load_dictionary(&path, layout, 128).unwrap_err().kind(),
            io::ErrorKind::InvalidData
        );
        std::fs::remove_file(path).unwrap();
    }
}
