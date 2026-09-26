//! BLAST protein alias (`.pal`) parsing.
//!
//! Mirrors `diamond/src/data/blastdb/pal.cpp` and `pal.h`. Volume decoding
//! remains in [`super::volume`], matching the upstream dependency direction.

use super::volume::BlastVolume;
use crate::util::string::ends_with;
use crate::util::system::{absolute_path, exists, is_absolute_path, PATH_SEPARATOR};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Read;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Pal {
    pub volumes: Vec<String>,
    pub metadata: BTreeMap<String, String>,
    pub oid_index: Vec<u64>,
    pub sequence_count: u64,
    pub letters: u64,
    pub version: i32,
}

impl Pal {
    pub fn new(path: &str) -> Result<Self, String> {
        let supported_keys: BTreeSet<&str> = [
            "TITLE",
            "MEMB_BIT",
            "SEQIDLIST",
            "NSEQ",
            "LENGTH",
            "TAXIDLIST",
        ]
        .into_iter()
        .collect();
        let (db_dir, file) = absolute_path(path);
        let mut pal = Self::default();
        if !exists(&format!("{path}.pal")) && !ends_with(path, ".pal") {
            pal.volumes.push(format!("{db_dir}{PATH_SEPARATOR}{file}"));
        } else {
            let pal_path = if ends_with(path, ".pal") {
                path.to_string()
            } else {
                format!("{path}.pal")
            };
            let mut text = String::new();
            std::fs::File::open(&pal_path)
                .map_err(|_| format!("Unable to open PAL file: {pal_path}"))?
                .read_to_string(&mut text)
                .map_err(|e| e.to_string())?;
            for (line_number0, mut line) in text.lines().map(str::to_string).enumerate() {
                let line_number = line_number0 + 1;
                if let Some(comment) = line.find('#') {
                    line.truncate(comment);
                }
                line = trim(&line);
                if line.is_empty() {
                    continue;
                }

                let key_end = line.find([' ', '\t']).ok_or_else(|| {
                    format!("Error parsing PAL file: line {line_number} is missing a value: {line}")
                })?;
                let key = line[..key_end].to_string();
                let value = trim(&line[key_end + 1..]);
                if value.is_empty() {
                    return Err(format!(
                        "Error parsing PAL file: line {line_number} has an empty value: {line}"
                    ));
                }

                if key == "DBLIST" {
                    let vls = split_whitespace(&value);
                    if vls.is_empty() {
                        return Err(format!(
                            "Error parsing PAL file: DBLIST on line {line_number} does not list any volumes"
                        ));
                    }
                    pal.volumes.extend(vls);
                    for s in &mut pal.volumes {
                        if !is_absolute_path(s) && !s.is_empty() && !s.starts_with('"') {
                            *s = format!("{db_dir}{PATH_SEPARATOR}{s}");
                        }
                    }
                    continue;
                }

                if !supported_keys.contains(key.as_str()) {
                    return Err(format!(
                        "Error parsing PAL file: Unsupported PAL key '{key}' on line {line_number}"
                    ));
                }
                if pal.metadata.contains_key(&key) {
                    return Err(format!(
                        "Error parsing PAL file: Duplicate key '{key}' on line {line_number}"
                    ));
                }
                pal.metadata.insert(key, value);
            }
        }

        pal.oid_index.push(0);
        let mut it = 0usize;
        while it < pal.volumes.len() {
            let volume = pal.volumes[it].clone();
            if volume.len() >= 2 && volume.starts_with('"') && volume.ends_with('"') {
                let nested = volume[1..volume.len() - 1].to_string();
                pal.volumes.remove(it);
                let nested_path = if is_absolute_path(&nested) {
                    nested
                } else {
                    format!("{db_dir}{PATH_SEPARATOR}{nested}")
                };
                it += pal.recurse(&nested_path, it)?;
            } else {
                let vol = BlastVolume::new(&volume, 0, 0, 0, false)?;
                pal.sequence_count += u64::from(vol.index().num_oids);
                pal.oid_index
                    .push(u64::from(vol.index().num_oids) + pal.oid_index.last().unwrap());
                pal.letters += vol.index().total_length;
                pal.version = vol.index().version as i32;
                it += 1;
            }
        }

        if let Some(seqidlist) = pal.metadata.get_mut("SEQIDLIST") {
            if ends_with(seqidlist, ".bsl") {
                return Err(format!(
                    "Binary SEQIDLIST files(.bsl) are not supported, use text file instead : {seqidlist}"
                ));
            }
            if !is_absolute_path(seqidlist) {
                *seqidlist = format!("{db_dir}{PATH_SEPARATOR}{seqidlist}");
            }
        }
        if let Some(taxidlist) = pal.metadata.get_mut("TAXIDLIST") {
            if !is_absolute_path(taxidlist) {
                *taxidlist = format!("{db_dir}{PATH_SEPARATOR}{taxidlist}");
            }
        }
        assert!(pal.sequence_count > 0);
        Ok(pal)
    }

    fn recurse(&mut self, path: &str, volume_it: usize) -> Result<usize, String> {
        let pal = Self::new(path)?;
        let inserted = pal.volumes.len();
        self.volumes.splice(volume_it..volume_it, pal.volumes);

        for (key, value) in pal.metadata {
            if self.metadata.contains_key(&key) {
                if key == "TITLE" || key == "NSEQ" || key == "LENGTH" {
                    continue;
                }
                return Err(format!("Duplicate key '{key}' in nested PAL file: {path}"));
            }
            self.metadata.insert(key, value);
        }
        let base = *self.oid_index.last().unwrap();
        self.oid_index
            .extend(pal.oid_index.iter().skip(1).map(|oid| oid + base));
        self.sequence_count += pal.sequence_count;
        self.letters += pal.letters;
        self.version = pal.version;
        Ok(inserted)
    }

    pub fn volume(&self, oid: u64) -> i32 {
        assert!(oid < self.sequence_count);
        self.oid_index.partition_point(|&x| x <= oid) as i32 - 1
    }
}

fn trim(text: &str) -> String {
    text.trim_matches([' ', '\t', '\r', '\n']).to_string()
}

fn split_whitespace(text: &str) -> Vec<String> {
    text.split_whitespace().map(str::to_string).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_pal(contents: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-pal-errors-{}-{}.pal",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    fn trim_and_split_match_stream_whitespace_semantics() {
        assert_eq!(trim(" \t alpha beta\r\n"), "alpha beta");
        assert_eq!(
            split_whitespace(" one\t two\r\nthree "),
            ["one", "two", "three"]
        );
    }

    #[test]
    fn rejects_missing_unsupported_and_duplicate_values() {
        let cases = [
            ("TITLE\n", "line 1 is missing a value"),
            ("TITLE  \t\n", "line 1 is missing a value"),
            ("UNKNOWN value\n", "Unsupported PAL key 'UNKNOWN' on line 1"),
            ("TITLE one\nTITLE two\n", "Duplicate key 'TITLE' on line 2"),
        ];
        for (text, expected) in cases {
            let path = temp_pal(text);
            let error = Pal::new(path.to_str().unwrap()).unwrap_err();
            assert!(error.contains(expected), "{error:?}");
            std::fs::remove_file(path).unwrap();
        }
    }

    #[test]
    fn comments_are_removed_before_validation() {
        let path = temp_pal("# comment\nTITLE # missing after comment\n");
        let error = Pal::new(path.to_str().unwrap()).unwrap_err();
        assert!(error.contains("line 2 is missing a value"), "{error:?}");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn volume_lookup_uses_upper_bound_at_boundaries() {
        let pal = Pal {
            oid_index: vec![0, 2, 5, 9],
            sequence_count: 9,
            ..Pal::default()
        };
        assert_eq!(pal.volume(0), 0);
        assert_eq!(pal.volume(1), 0);
        assert_eq!(pal.volume(2), 1);
        assert_eq!(pal.volume(8), 2);
    }
}
