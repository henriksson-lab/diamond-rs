//! Ordered two-way TSV joins.
//!
//! This mirrors `diamond/src/util/tsv/join.cpp`. Inputs must be sorted by
//! their respective key columns. Equal keys are consumed pairwise rather than
//! expanded into a Cartesian product.

use super::{Config, File, FileColumn, Flags, Schema, Type, Value};

/// Build the output schema by projecting columns from the two input schemas.
pub fn join_schema(schemas: &[Schema; 2], output_fields: &[FileColumn]) -> Schema {
    let mut schema = Vec::with_capacity(output_fields.len());
    for field in output_fields {
        schema.push(schemas[field.file][field.column]);
    }
    schema
}

/// Join two sorted files into an existing output file.
pub fn join(
    file1: &mut File,
    file2: &mut File,
    column1: usize,
    column2: usize,
    output_fields: &[FileColumn],
    out: &mut File,
) -> Result<(), String> {
    if output_fields.is_empty() {
        return Err("Join with empty output".to_string());
    }

    let mut tables = [file1.read_record(), file2.read_record()];
    while !tables[0].empty() && !tables[1].empty() {
        let keys = [
            tables[0].front().get_t::<i64>(column1),
            tables[1].front().get_t::<i64>(column2),
        ];
        if keys[0] < keys[1] {
            tables[0] = file1.read_record();
        } else if keys[1] < keys[0] {
            tables[1] = file2.read_record();
        } else {
            let mut values = Vec::with_capacity(output_fields.len());
            for field in output_fields {
                values.push(match tables[field.file].schema()[field.column] {
                    Type::String => {
                        Value::String(tables[field.file].front().get_t::<String>(field.column))
                    }
                    Type::Int64 => {
                        Value::Int64(tables[field.file].front().get_t::<i64>(field.column))
                    }
                });
            }
            out.write_record_values(&values)?;
            tables[0] = file1.read_record();
            tables[1] = file2.read_record();
        }
    }
    Ok(())
}

/// Join two sorted files into a temporary file with the projected schema.
pub fn joined(
    file1: &mut File,
    file2: &mut File,
    column1: usize,
    column2: usize,
    output_fields: &[FileColumn],
) -> Result<File, String> {
    let schemas = [file1.schema(), file2.schema()];
    let schema = join_schema(&schemas, output_fields);
    let mut out = File::new(schema, "", Flags::TEMP, Config::default())?;

    // Upstream currently calls the output overload with `file1` twice here.
    // Passing both supplied files is the intended behavior and preserves the
    // established Rust API rather than reproducing that apparent typo.
    join(file1, file2, column1, column2, output_fields, &mut out)?;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const LEFT_VALUE: FileColumn = FileColumn { file: 0, column: 1 };
    const RIGHT_VALUE: FileColumn = FileColumn { file: 1, column: 1 };

    fn schema() -> Schema {
        vec![Type::Int64, Type::String]
    }

    #[test]
    fn schema_projection_preserves_requested_order_and_duplicates() {
        let schemas = [
            vec![Type::Int64, Type::String],
            vec![Type::String, Type::Int64],
        ];
        let fields = [
            FileColumn { file: 1, column: 1 },
            FileColumn { file: 0, column: 1 },
            FileColumn { file: 1, column: 1 },
        ];

        assert_eq!(
            join_schema(&schemas, &fields),
            vec![Type::Int64, Type::String, Type::Int64]
        );
    }

    #[test]
    fn duplicate_keys_pair_one_for_one_and_keep_input_order() {
        let mut left = File::from_lines(schema(), ["1\ta", "1\tb", "2\tc", "4\td"]).unwrap();
        let mut right =
            File::from_lines(schema(), ["1\tx", "1\ty", "1\tz", "3\tq", "4\tw"]).unwrap();

        let out = joined(&mut left, &mut right, 0, 0, &[LEFT_VALUE, RIGHT_VALUE]).unwrap();

        assert_eq!(out.schema_ref(), &[Type::String, Type::String]);
        assert_eq!(out.table().size(), 3);
        assert_eq!(out.table().record(0).get_t::<String>(0), "a");
        assert_eq!(out.table().record(0).get_t::<String>(1), "x");
        assert_eq!(out.table().record(1).get_t::<String>(0), "b");
        assert_eq!(out.table().record(1).get_t::<String>(1), "y");
        assert_eq!(out.table().record(2).get_t::<String>(0), "d");
        assert_eq!(out.table().record(2).get_t::<String>(1), "w");
        let mut text = crate::util::text_buffer::TextBuffer::new();
        out.table().write(&mut text);
        assert_eq!(text.data(), b"a\tx\nb\ty\nd\tw\n");
    }

    #[test]
    fn empty_input_produces_empty_projected_output() {
        let mut left = File::from_lines(schema(), std::iter::empty::<&str>()).unwrap();
        let mut right = File::from_lines(schema(), ["1\tx"]).unwrap();

        let out = joined(&mut left, &mut right, 0, 0, &[RIGHT_VALUE, LEFT_VALUE]).unwrap();

        assert_eq!(out.schema_ref(), &[Type::String, Type::String]);
        assert!(out.table().empty());
    }

    #[test]
    fn empty_output_is_rejected_before_consuming_inputs() {
        let mut left = File::from_lines(schema(), ["1\ta"]).unwrap();
        let mut right = File::from_lines(schema(), ["1\tx"]).unwrap();
        let mut out = File::from_table(super::super::Table::new(vec![]), Flags::TEMP);

        assert_eq!(
            join(&mut left, &mut right, 0, 0, &[], &mut out).unwrap_err(),
            "Join with empty output"
        );
        assert!(!left.eof());
        assert!(!right.eof());
    }
}
