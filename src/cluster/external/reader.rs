//! Bounded SequenceFile::auto_create volume traversal shared by external stages.
use super::seed_table::{ByteSeedSequence, SeedSequence};
use crate::data::sequence_file::{
    detect_type_source_with_buffer, SequenceFileConfig, SequenceFileFlags as Flags,
    SequenceFileType,
};
use std::path::Path;

pub fn visit_volume(
    path: &Path,
    visit: &mut dyn FnMut(SeedSequence) -> Result<(), String>,
) -> Result<(), String> {
    visit_volume_owned(path, visit, &mut || Ok(()), &Default::default())
}

pub fn visit_volume_owned(
    path: &Path,
    visit: &mut dyn FnMut(SeedSequence) -> Result<(), String>,
    on_open: &mut dyn FnMut() -> Result<(), String>,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<(), String> {
    visit_volume_bytes_owned(
        path,
        &mut |record| {
            visit(SeedSequence::<String>::try_from(record).map_err(|error| error.to_string())?)
        },
        on_open,
        errors,
    )
}

pub fn visit_volume_bytes_owned(
    path: &Path,
    visit: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    on_open: &mut dyn FnMut() -> Result<(), String>,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<(), String> {
    visit_volume_bytes_while_owned(path, &mut || true, visit, on_open, errors)
}

/// Evaluate Source's loop condition before reading each record from the retained
/// owner. A false condition does not decode or validate a trailing record.
pub fn visit_volume_bytes_while_owned(
    path: &Path,
    should_read: &mut dyn FnMut() -> bool,
    visit: &mut dyn FnMut(ByteSeedSequence) -> Result<(), String>,
    on_open: &mut dyn FnMut() -> Result<(), String>,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<(), String> {
    let (kind, resolved) = errors.capture(detect_type_source_with_buffer(
        &[path.to_path_buf()],
        Flags::NONE,
        SequenceFileConfig::default(),
        65_536,
    ))?;
    match kind {
        SequenceFileType::Block => Err("Block-backed volume has no file owner.".to_owned()),
        SequenceFileType::Dmnd => {
            let mut file =
                errors.capture(crate::data::dmnd::DatabaseFile::open_with_options_io(
                    resolved,
                    crate::data::dmnd::DatabaseOpenOptions {
                        file_buffer_size: Some(65_536),
                        ..Default::default()
                    },
                ))?;
            on_open()?;
            for oid in 0..file.sequence_count() {
                if !should_read() {
                    break;
                }
                let (sequence, id) = errors.capture(file.record_bytes_io(oid))?;
                visit(ByteSeedSequence { id, sequence })?;
            }
            Ok(())
        }
        SequenceFileType::Blast => {
            let mut file = errors.capture(
                crate::data::blastdb::blastdb::BlastDB::new_path_with_io_owner(
                    &resolved,
                    Flags::SEQS | Flags::TITLES | Flags::FULL_TITLES,
                    errors,
                ),
            )?;
            on_open()?;
            for oid in 0..file.sequence_count() {
                if !should_read() {
                    break;
                }
                let (sequence, id) = errors.capture(file.record_at_bytes_io(oid))?;
                visit(ByteSeedSequence { id, sequence })?;
            }
            Ok(())
        }
        SequenceFileType::Fasta => {
            let mut file = errors.capture(
                crate::data::fasta::fasta_file::FastaFile::open_source_with_buffer(
                    &[resolved],
                    Default::default(),
                    65_536,
                ),
            )?;
            file.set_io_errors(errors.clone());
            on_open()?;
            while should_read() {
                let Some(record) = errors.capture(file.read_seq_bytes_io())? else {
                    break;
                };
                visit(ByteSeedSequence {
                    id: record.id,
                    sequence: record.sequence,
                })?;
            }
            errors.capture(file.close_io())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT: AtomicU64 = AtomicU64::new(0);

    #[test]
    fn actual_volume_reader_retains_non_utf8_fasta_title() {
        let path = std::env::temp_dir().join(format!(
            "diamond-raw-volume-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, b">q\xff description\nARND\n").unwrap();
        let mut ids = Vec::new();
        visit_volume_bytes_owned(
            &path,
            &mut |record| {
                ids.push(record.id);
                Ok(())
            },
            &mut || Ok(()),
            &Default::default(),
        )
        .unwrap();
        assert_eq!(ids, vec![b"q\xff description".to_vec()]);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn pre_read_condition_skips_trailing_invalid_sequence_but_keeps_constructor() {
        let path = std::env::temp_dir().join(format!(
            "diamond-volume-stop-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, b">first\nARND\n>invalid\n!\n").unwrap();
        let read = std::cell::Cell::new(false);
        let mut opened = false;
        visit_volume_bytes_while_owned(
            &path,
            &mut || !read.get(),
            &mut |record| {
                assert_eq!(record.id, b"first");
                read.set(true);
                Ok(())
            },
            &mut || {
                opened = true;
                Ok(())
            },
            &Default::default(),
        )
        .unwrap();
        assert!(opened && read.get());
        opened = false;
        visit_volume_bytes_while_owned(
            &path,
            &mut || false,
            &mut |_| panic!("false pre-read condition must not visit a record"),
            &mut || {
                opened = true;
                Ok(())
            },
            &Default::default(),
        )
        .unwrap();
        assert!(opened);
        // An unrestricted visit still validates and rejects the bad tail.
        assert!(visit_volume_bytes_owned(
            &path,
            &mut |_| Ok(()),
            &mut || Ok(()),
            &Default::default()
        )
        .is_err());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn failed_input_construction_precedes_the_output_open_callback() {
        let directory = std::env::temp_dir().join(format!(
            "diamond-external-reader-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&directory).unwrap();
        for (name, contents) in [
            ("missing", None),
            ("empty", Some(&b""[..])),
            ("invalid", Some(&b"\xff"[..])),
        ] {
            let path = directory.join(name);
            if let Some(contents) = contents {
                std::fs::write(&path, contents).unwrap();
            }
            let mut opened = false;
            let errors = crate::util::io::IoErrorOwner::default();
            let result = visit_volume_owned(
                &path,
                &mut |_| panic!("failed constructor must not visit records"),
                &mut || {
                    opened = true;
                    Ok(())
                },
                &errors,
            );
            assert!(!opened, "{name} opened its output before input validation");
            let error = errors.finish(result).unwrap_err();
            if name == "missing" {
                assert_eq!(
                    Path::new(&crate::util::io::file_open_error::find(&error).unwrap().path),
                    path
                );
            } else {
                assert!(error
                    .to_string()
                    .starts_with("Error detecting input file format."));
                assert_eq!(error.to_string().contains("empty"), name == "empty");
            }
        }
        std::fs::remove_dir_all(directory).unwrap();
    }
}
