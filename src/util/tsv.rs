pub type Schema = Vec<Type>;
pub type RecordId = i64;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Type {
    String,
    Int64,
}

pub mod record;
pub use record::{FromTsvRecord, Record, RecordField, RecordIterator, TsvType, TsvValue};
pub mod tsv;
pub use tsv::{count_lines, count_lines_str};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InvalidType;

impl std::fmt::Display for InvalidType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Invalid type in schema.")
    }
}

impl std::error::Error for InvalidType {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SchemaMismatch;

impl std::fmt::Display for SchemaMismatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Mismatching schema.")
    }
}

impl std::error::Error for SchemaMismatch {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Value {
    String(String),
    Int64(i64),
}

pub mod table;
pub use table::Table;
pub mod read_text_mt;
pub use read_text_mt::{read_text_mt, READ_TEXT_MT_SIZE};
pub mod file;
pub use file::{Config, File, FileColumn, Flags};
pub mod join;
pub use join::{join, join_schema, joined};
pub mod merge;
pub use merge::{merge, merge as merge_files};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BuildHelper {
    pub buffer: Vec<u8>,
    pub limits: Vec<i64>,
    pub counts: Vec<i32>,
}

impl BuildHelper {
    pub fn new(size: i64) -> Self {
        const TOLERANCE: f64 = 1.1;
        let mut buffer = Vec::new();
        buffer.reserve((size as f64 * TOLERANCE) as usize);
        Self {
            buffer,
            limits: vec![0],
            counts: Vec::new(),
        }
    }

    pub fn add(&mut self, t: &Table) {
        self.buffer.extend_from_slice(&t.data);
        self.limits.extend_from_slice(&t.limits[1..]);
        self.counts.push(t.size() as i32);
    }

    pub fn get(mut self, schema: Schema) -> Table {
        let n = self.counts.len();
        if n == 0 {
            return Table::new(schema);
        }

        let mut offsets = Vec::with_capacity(n);
        offsets.push(0);
        let mut it = self.counts[0] as usize;
        for i in 1..n {
            offsets.push(self.limits[it] + offsets[i - 1]);
            it += self.counts[i] as usize;
        }

        let mut begin = self.counts[0] as usize + 1;
        for (i, d) in offsets.iter().enumerate().skip(1) {
            let end = begin + self.counts[i] as usize;
            for limit in &mut self.limits[begin..end] {
                *limit += *d;
            }
            begin += self.counts[i] as usize;
        }
        Table::from_parts(schema, self.buffer, self.limits)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::text_buffer::TextBuffer;

    fn as_string(buf: &TextBuffer) -> String {
        String::from_utf8(buf.data().to_vec()).unwrap()
    }

    #[test]
    fn test_record_and_table_write() {
        let schema = vec![Type::Int64, Type::String];
        let mut table = Table::new(schema.clone());
        table
            .write_record(&[Value::Int64(7), Value::String("abc".to_string())])
            .unwrap();
        table.push_back_line("9\tdef", None).unwrap();
        assert_eq!(table.size(), 2);
        assert_eq!(table.front().get_t::<i64>(0), 7);
        assert_eq!(table.record(1).get(0), "9");
        assert_eq!(table.record(1).get_t::<String>(1), "def");

        let mut out = TextBuffer::new();
        table.write(&mut out);
        assert_eq!(as_string(&out), "7\tabc\n9\tdef\n");
    }

    #[test]
    fn test_append_sort_shuffle_and_alloc() {
        let schema = vec![Type::String, Type::Int64];
        let mut a = Table::new(schema.clone());
        a.write_record(&[Value::String("b".to_string()), Value::Int64(2)])
            .unwrap();
        let mut b = Table::new(schema.clone());
        b.write_record(&[Value::String("a".to_string()), Value::Int64(1)])
            .unwrap();
        a.append(&b).unwrap();
        assert!(a.alloc_size() > 0);

        let sorted = a.sorted(1, 1).unwrap();
        assert_eq!(sorted.record(0).get_t::<String>(0), "a");
        assert_eq!(sorted.record(1).get_t::<i64>(1), 2);

        let shuffled = sorted.shuffle([1, 0]);
        assert_eq!(shuffled.record(0).get_t::<String>(0), "b");
    }

    #[test]
    fn test_table_map_writes_mapped_records() {
        let schema = vec![Type::Int64, Type::String];
        let mut table = Table::new(schema.clone());
        table
            .write_record(&[Value::Int64(1), Value::String("a".to_string())])
            .unwrap();
        table
            .write_record(&[Value::Int64(2), Value::String("b".to_string())])
            .unwrap();
        let mut out = File::from_table(Table::new(schema.clone()), Flags::TEMP);

        table
            .map(
                1,
                |record| {
                    let mut mapped = Table::new(schema.clone());
                    mapped
                        .write_record(&[
                            Value::Int64(record.get_t::<i64>(0) * 10),
                            Value::String(record.get_t::<String>(1)),
                        ])
                        .unwrap();
                    mapped
                },
                &mut out,
            )
            .unwrap();

        assert_eq!(out.table().record(0).get_t::<i64>(0), 10);
        assert_eq!(out.table().record(1).get_t::<String>(1), "b");
    }

    #[test]
    fn test_build_helper() {
        let schema = vec![Type::Int64, Type::String];
        let mut a = Table::new(schema.clone());
        a.write_record(&[Value::Int64(1), Value::String("x".to_string())])
            .unwrap();
        let mut b = Table::new(schema.clone());
        b.write_record(&[Value::Int64(2), Value::String("y".to_string())])
            .unwrap();

        let mut helper = BuildHelper::new(a.alloc_size() + b.alloc_size());
        helper.add(&a);
        helper.add(&b);
        let table = helper.get(schema);
        assert_eq!(table.size(), 2);
        assert_eq!(table.record(0).get_t::<String>(1), "x");
        assert_eq!(table.record(1).get_t::<i64>(0), 2);
    }

    #[test]
    fn test_errors() {
        let mut table = Table::new(vec![Type::Int64, Type::String]);
        assert!(table.push_back_line("1", None).is_err());
        assert_eq!(table.size(), 1);
        assert_eq!(
            table
                .write_record(&[
                    Value::Int64(1),
                    Value::String("x".to_string()),
                    Value::Int64(2)
                ])
                .unwrap_err(),
            "write_record with too many fields."
        );
        assert_eq!(
            table.write_record(&[Value::Int64(1)]).unwrap_err(),
            "Mismatching field count for Table::write_record"
        );
        assert_eq!(
            table
                .write_record(&[Value::String("bad".to_string()), Value::Int64(1)])
                .unwrap_err(),
            "Invalid type in schema"
        );
        assert!(table.sorted(1, 1).is_err());

        let other = Table::new(vec![Type::String]);
        assert_eq!(table.append(&other).unwrap_err(), SchemaMismatch);
    }

    #[test]
    fn test_file_read_write_sort_map() {
        let schema = vec![Type::Int64, Type::String];
        let mut file = File::new(schema.clone(), "", Flags::TEMP, Config::default()).unwrap();
        file.write_record_values(&[Value::Int64(2), Value::String("b".to_string())])
            .unwrap();
        file.write_record_values(&[Value::Int64(1), Value::String("a".to_string())])
            .unwrap();
        assert_eq!(file.read_record().front().get_t::<i64>(0), 2);
        assert_eq!(file.read_record().front().get_t::<String>(1), "a");
        assert!(file.read_record().empty());
        file.rewind();
        let sorted = file.sort(0, 1).unwrap();
        assert_eq!(sorted.table().front().get_t::<i64>(0), 1);
        let mapped = file.map(1, |r| {
            let mut t = Table::new(schema.clone());
            t.write_record(&[
                Value::Int64(r.get_t::<i64>(0) + 10),
                Value::String(r.get_t::<String>(1)),
            ])
            .unwrap();
            t
        });
        assert_eq!(mapped.table().record(0).get_t::<i64>(0), 12);

        let mut typed = Vec::<(i64, String)>::new();
        file.read_typed(&mut typed).unwrap();
        assert_eq!(typed, vec![(2, "b".to_string()), (1, "a".to_string())]);
        assert_eq!(
            file.read_typed::<(String, i64)>(&mut Vec::new())
                .unwrap_err(),
            "Template parameters do not match schema."
        );
    }

    #[test]
    fn test_file_new_reads_path_and_record_id_column() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-tsv-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, "11\talpha\r\n12\tbeta\r\n").unwrap();
        let name = path.to_string_lossy().into_owned();

        let mut file = File::new(
            vec![Type::Int64, Type::String],
            &name,
            Flags::READ,
            Config::default(),
        )
        .unwrap();
        assert_eq!(file.file_name(), name);
        assert_eq!(file.read_record().front().get_t::<i64>(0), 11);
        assert_eq!(file.read_record().front().get_t::<String>(1), "beta");
        assert!(file.read_record().empty());

        let mut record_id_file = File::new(
            vec![Type::Int64, Type::Int64, Type::String],
            &name,
            Flags::READ | Flags::RECORD_ID_COLUMN,
            Config::default(),
        )
        .unwrap();
        let first = record_id_file.read_record();
        assert_eq!(first.front().get_t::<i64>(0), 0);
        assert_eq!(first.front().get_t::<i64>(1), 11);
        let second = record_id_file.read_record();
        assert_eq!(second.front().get_t::<i64>(0), 1);
        assert_eq!(second.front().get_t::<String>(2), "beta");

        let mut out = File::new(
            vec![Type::Int64, Type::Int64, Type::String],
            "",
            Flags::TEMP | Flags::RECORD_ID_COLUMN,
            Config::default(),
        )
        .unwrap();
        out.write_record_values(&[Value::Int64(7), Value::String("x".to_string())])
            .unwrap();
        out.write_record_strings(&["8", "y"]).unwrap();
        assert_eq!(out.table().record(0).get_t::<i64>(0), 0);
        assert_eq!(out.table().record(1).get_t::<i64>(0), 1);
        assert_eq!(out.table().record(1).get_t::<String>(2), "y");

        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_read_text_mt_extends_full_chunks_to_record_boundary() {
        let mut data = vec![b'a'; READ_TEXT_MT_SIZE];
        data.extend_from_slice(b"tail\nb\n");
        let mut chunks = Vec::new();
        read_text_mt(&data, i64::MAX, 4, |chunk, bytes| {
            chunks.push((chunk, bytes.len(), bytes.last().copied()));
        });
        assert_eq!(
            chunks,
            vec![(0, READ_TEXT_MT_SIZE + 4, Some(b'l')), (1, 2, Some(b'\n'))]
        );

        let mut limited = Vec::new();
        read_text_mt(&data, 1, 1, |chunk, bytes| {
            limited.push((chunk, bytes.len()));
        });
        assert_eq!(limited, vec![(0, READ_TEXT_MT_SIZE + 4)]);
    }

    #[test]
    fn test_file_read_raw_chunks_and_read_tables() {
        let schema = vec![Type::Int64, Type::String];
        let mut file = File::from_lines(schema.clone(), ["1\ta", "2\tb"]).unwrap();
        let mut raw = Vec::new();
        file.read_raw_chunks(i64::MAX, 1, |chunk, bytes| {
            raw.push((chunk, String::from_utf8(bytes.to_vec()).unwrap()));
        });
        assert_eq!(raw, vec![(0, "1\ta\n2\tb\n".to_string())]);

        let mut tables = Vec::new();
        file.read_tables(1, |chunk, table| {
            tables.push((chunk, table.size(), table.record(1).get_t::<String>(1)));
        })
        .unwrap();
        assert_eq!(tables, vec![(0, 2, "b".to_string())]);
    }

    #[test]
    fn test_merge_join_and_count_lines() {
        let schema = vec![Type::Int64, Type::String];
        let mut a = File::from_lines(schema.clone(), ["1\ta", "3\tc"]).unwrap();
        let mut b = File::from_lines(schema.clone(), ["2\tb", "4\td"]).unwrap();
        let merged = merge_files(&mut [&mut a, &mut b], 0).unwrap();
        assert_eq!(merged.table().size(), 4);
        assert_eq!(merged.table().record(2).get_t::<i64>(0), 3);

        let mut j1 = File::from_lines(schema.clone(), ["1\ta", "3\tc"]).unwrap();
        let mut j2 = File::from_lines(schema, ["1\tx", "2\ty", "3\tz"]).unwrap();
        let out = joined(
            &mut j1,
            &mut j2,
            0,
            0,
            &[
                FileColumn { file: 0, column: 1 },
                FileColumn { file: 1, column: 1 },
            ],
        )
        .unwrap();
        assert_eq!(out.table().size(), 2);
        assert_eq!(out.table().record(0).get_t::<String>(0), "a");
        assert_eq!(out.table().record(1).get_t::<String>(1), "z");
        assert_eq!(count_lines_str("a\nb\n"), 2);
        assert_eq!(count_lines_str("a\nb"), 2);
        assert_eq!(count_lines_str("\n"), 1);
        assert_eq!(count_lines_str("\n\n"), 2);
        assert_eq!(count_lines_str(""), 0);
    }
}
