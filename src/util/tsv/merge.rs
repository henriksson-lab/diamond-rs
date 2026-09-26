//! Ordered multi-file TSV merge.
//!
//! This mirrors `diamond/src/util/tsv/merge.cpp`. Each input must already be
//! sorted on the selected signed 64-bit key column.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use super::{Config, File, Flags};

/// Merge sorted files, ordering equal keys by their input-file index.
pub fn merge(files: &mut [&mut File], column: usize) -> Result<File, String> {
    // Upstream dereferences `begin` unconditionally. Retain the established
    // Rust API's defined error instead of reproducing undefined behavior.
    if files.is_empty() {
        return Err("merge with empty input".to_string());
    }

    let schema = files[0].schema();
    let mut out = File::new(schema, "", Flags::TEMP, Config::default())?;
    let mut queue = BinaryHeap::<Reverse<(i64, usize)>>::new();
    let mut tables = Vec::with_capacity(files.len());
    for (i, file) in files.iter_mut().enumerate() {
        let table = file.read_record();
        if !table.empty() {
            queue.push(Reverse((table.front().get_t::<i64>(column), i)));
        }
        tables.push(table);
    }

    while let Some(Reverse((_, file_index))) = queue.pop() {
        out.write(&tables[file_index].front())
            .map_err(|e| e.to_string())?;
        tables[file_index] = files[file_index].read_record();
        if !tables[file_index].empty() {
            queue.push(Reverse((
                tables[file_index].front().get_t::<i64>(column),
                file_index,
            )));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::text_buffer::TextBuffer;
    use crate::util::tsv::{Schema, Type};

    fn schema() -> Schema {
        vec![Type::Int64, Type::String]
    }

    fn as_text(file: &File) -> Vec<u8> {
        let mut text = TextBuffer::new();
        file.table().write(&mut text);
        text.data().to_vec()
    }

    #[test]
    fn equal_keys_are_ordered_by_input_index_before_next_key() {
        let mut first = File::from_lines(schema(), ["1\ta0", "1\ta1", "3\ta3"]).unwrap();
        let mut second = File::from_lines(schema(), ["1\tb0", "1\tb1", "2\tb2"]).unwrap();
        let mut third = File::from_lines(schema(), ["1\tc0", "4\tc4"]).unwrap();

        let merged = merge(&mut [&mut first, &mut second, &mut third], 0).unwrap();

        assert_eq!(merged.schema_ref(), schema());
        assert_eq!(
            as_text(&merged),
            b"1\ta0\n1\ta1\n1\tb0\n1\tb1\n1\tc0\n2\tb2\n3\ta3\n4\tc4\n"
        );
    }

    #[test]
    fn empty_files_are_skipped_but_first_schema_is_selected() {
        let first_schema = schema();
        let mut empty = File::from_lines(first_schema.clone(), std::iter::empty::<&str>()).unwrap();
        let mut populated = File::from_lines(schema(), ["2\tb", "5\te"]).unwrap();

        let merged = merge(&mut [&mut empty, &mut populated], 0).unwrap();

        assert_eq!(merged.schema_ref(), first_schema);
        assert_eq!(as_text(&merged), b"2\tb\n5\te\n");
    }

    #[test]
    fn all_empty_files_produce_empty_output_with_first_schema() {
        let first_schema = vec![Type::Int64];
        let mut first = File::from_lines(first_schema.clone(), std::iter::empty::<&str>()).unwrap();
        let mut second = File::from_lines(schema(), std::iter::empty::<&str>()).unwrap();

        let merged = merge(&mut [&mut first, &mut second], 0).unwrap();

        assert_eq!(merged.schema_ref(), first_schema);
        assert!(merged.table().empty());
    }

    #[test]
    fn empty_input_has_defined_compatibility_error() {
        let mut files: [&mut File; 0] = [];
        assert_eq!(merge(&mut files, 0).unwrap_err(), "merge with empty input");
    }

    #[test]
    fn mismatching_nonempty_schema_is_reported() {
        let mut first = File::from_lines(schema(), ["1\ta"]).unwrap();
        let mut second = File::from_lines(vec![Type::Int64], ["2"]).unwrap();

        assert_eq!(
            merge(&mut [&mut first, &mut second], 0).unwrap_err(),
            "Mismatching schema."
        );
    }
}
