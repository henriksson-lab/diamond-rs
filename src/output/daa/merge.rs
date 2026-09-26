//! Mirrored implementation boundary for `diamond/src/output/daa/merge.cpp`.
//!
//! The DAA wire-format implementation remains in [`crate::data::daa`]; these
//! wrappers preserve the upstream source hierarchy and make command-global
//! configuration explicit.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use crate::data::daa::{
    copy_match_record_raw, finish_daa_from_refs, finish_daa_query_record, init_daa,
    write_daa_query_record, DaaFile, DaaQueryRecord,
};

/// Explicit counterpart of the C++ `config` fields consumed by `merge_daa`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DaaMergeConfig {
    pub input_files: Vec<PathBuf>,
    pub output_file: PathBuf,
}

impl DaaMergeConfig {
    pub fn new<I, O>(input_files: I, output_file: O) -> Self
    where
        I: IntoIterator,
        I::Item: Into<PathBuf>,
        O: Into<PathBuf>,
    {
        Self {
            input_files: input_files.into_iter().map(Into::into).collect(),
            output_file: output_file.into(),
        }
    }

    pub fn validate(&self) -> io::Result<()> {
        if self.input_files.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Missing parameter: input files (--in)",
            ));
        }
        if self.output_file.as_os_str().is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Missing parameter: output file (--out)",
            ));
        }
        Ok(())
    }
}

/// Mirrors C++ `build_mapping`.
pub fn build_mapping(
    accession_to_oid: &mut HashMap<String, u32>,
    sequence_ids: &mut Vec<String>,
    sequence_lengths: &mut Vec<u32>,
    file: &DaaFile,
) -> HashMap<u32, u32> {
    let mut mapping = HashMap::new();
    for input_oid in 0..file.db_seqs_used() as usize {
        let name = file.ref_name(input_oid).to_string();
        let next_oid = accession_to_oid.len() as u32;
        let output_oid = *accession_to_oid.entry(name.clone()).or_insert(next_oid);
        mapping.insert(input_oid as u32, output_oid);
        if output_oid == next_oid {
            sequence_ids.push(name);
            sequence_lengths.push(file.ref_len_at(input_oid));
        }
    }
    mapping
}

/// Mirrors C++ `write_file`, including raw match copying and subject-ID remap.
pub fn write_file<W: Write>(
    file: &mut DaaFile,
    output: &mut W,
    subject_map: &HashMap<u32, u32>,
) -> io::Result<i64> {
    let mut output_buffer = Vec::new();
    let mut last_query_number = None;
    while let Some((buffer, query_number)) = file.read_query_buffer()? {
        let record = DaaQueryRecord::from_buffer(file, &buffer, query_number)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
        let size_position = write_daa_query_record(
            &mut output_buffer,
            &record.query_name,
            record.query_source(),
            record.input_sequence_type(),
        );
        let mut matches = record.raw_begin();
        while matches.good() {
            copy_match_record_raw(&mut matches, &mut output_buffer, subject_map)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
        }
        finish_daa_query_record(&mut output_buffer, size_position);
        output.write_all(&output_buffer)?;
        output_buffer.clear();
        last_query_number = Some(query_number);
    }
    Ok(last_query_number.map_or(0, |query_number| query_number as i64 + 1))
}

/// Mirrors C++ `merge_daa`, with its global config made explicit.
pub fn merge_daa(config: &DaaMergeConfig) -> io::Result<i64> {
    config.validate()?;

    let mut files = Vec::with_capacity(config.input_files.len());
    let mut accession_to_oid = HashMap::new();
    let mut oid_maps = Vec::with_capacity(config.input_files.len());
    let mut sequence_ids = Vec::new();
    let mut sequence_lengths = Vec::new();

    for input_file in &config.input_files {
        let file = DaaFile::open(input_file)?;
        oid_maps.push(build_mapping(
            &mut accession_to_oid,
            &mut sequence_ids,
            &mut sequence_lengths,
            &file,
        ));
        files.push(file);
    }

    let mut output = File::create(&config.output_file)?;
    init_daa(&mut output)?;
    let mut query_count = 0;
    for (file, subject_map) in files.iter_mut().zip(&oid_maps) {
        query_count += write_file(file, &mut output, subject_map)?;
    }
    finish_daa_from_refs(
        &mut output,
        &files[0],
        &sequence_ids,
        &sequence_lengths,
        query_count,
    )?;
    Ok(query_count)
}

/// Compatibility entry point retaining the existing path-oriented API.
pub fn merge_daa_files<P: AsRef<Path>>(input_files: &[P], output_file: P) -> io::Result<i64> {
    crate::data::daa::merge_daa_files(input_files, output_file)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn explicit_config_matches_cpp_parameter_validation() {
        let missing_inputs = DaaMergeConfig::new(Vec::<PathBuf>::new(), "out.daa");
        let error = merge_daa(&missing_inputs).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
        assert_eq!(error.to_string(), "Missing parameter: input files (--in)");

        let missing_output = DaaMergeConfig::new(["in.daa"], "");
        let error = merge_daa(&missing_output).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
        assert_eq!(error.to_string(), "Missing parameter: output file (--out)");
    }

    #[test]
    fn config_preserves_input_order_and_paths() {
        let config = DaaMergeConfig::new(["first.daa", "second.daa"], "merged.daa");
        assert_eq!(
            config.input_files,
            [PathBuf::from("first.daa"), PathBuf::from("second.daa")]
        );
        assert_eq!(config.output_file, PathBuf::from("merged.daa"));
        assert!(config.validate().is_ok());
    }
}
