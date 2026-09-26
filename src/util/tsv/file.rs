//! Typed TSV file facade and chunked text reads.
//!
//! This mirrors `diamond/src/util/tsv/file.{h,cpp}`. The Rust facade keeps an
//! eagerly parsed table, while preserving the upstream cursor and output order.

use crate::util::text_buffer::TextBuffer;

use super::read_text_mt::{read_text_mt, READ_TEXT_MT_SIZE};
use super::{FromTsvRecord, Record, RecordId, Schema, SchemaMismatch, Table, Type, Value};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Flags(pub u32);

impl Flags {
    pub const READ: Self = Self(0);
    pub const READ_WRITE: Self = Self(1);
    pub const WRITE: Self = Self(1 << 1);
    pub const OVERWRITE: Self = Self(1 << 2);
    pub const RECORD_ID_COLUMN: Self = Self(1 << 3);
    pub const TEMP: Self = Self(1 << 4);

    pub fn any(self, other: Self) -> bool {
        (self.0 & other.0) != 0
    }

    pub fn all(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }
}

impl std::ops::BitOr for Flags {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl std::ops::BitOrAssign for Flags {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FileColumn {
    pub file: usize,
    pub column: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Config {
    pub line_delimiter: u8,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            line_delimiter: b'\n',
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct File {
    flags: Flags,
    schema: Schema,
    config: Config,
    table: Table,
    cursor: i64,
    write_buf: TextBuffer,
    record_id: RecordId,
    file_name: String,
    raw_text: Option<String>,
    closed: bool,
}

impl File {
    pub fn new(
        schema: Schema,
        file_name: &str,
        flags: Flags,
        config: Config,
    ) -> Result<Self, String> {
        let flags = get_flags(flags)?;
        if flags.all(Flags::WRITE | Flags::OVERWRITE) {
            return Err("Invalid File flags".to_string());
        }
        if flags.any(Flags::RECORD_ID_COLUMN) && schema.first() != Some(&Type::Int64) {
            return Err("Schema does not contain record_id column.".to_string());
        }
        let mut table = Table::new(schema.clone());
        let mut raw_text = None;
        if !flags.any(Flags::WRITE) || flags.any(Flags::READ_WRITE) {
            if !file_name.is_empty() && !flags.any(Flags::TEMP) {
                let text = std::fs::read_to_string(file_name).map_err(|e| e.to_string())?;
                for line in text.split_terminator(config.line_delimiter as char) {
                    let line = line.strip_suffix('\r').unwrap_or(line);
                    let record_id = flags
                        .any(Flags::RECORD_ID_COLUMN)
                        .then_some(table.size() as RecordId);
                    table.push_back_line(line, record_id)?;
                }
                raw_text = Some(text);
            }
        }
        Ok(Self {
            table,
            flags,
            schema,
            config,
            cursor: 0,
            write_buf: TextBuffer::new(),
            record_id: 0,
            file_name: file_name.to_string(),
            raw_text,
            closed: false,
        })
    }

    pub fn from_table(table: Table, flags: Flags) -> Self {
        let schema = table.schema.clone();
        Self {
            flags,
            schema,
            config: Config::default(),
            table,
            cursor: 0,
            write_buf: TextBuffer::new(),
            record_id: 0,
            file_name: String::new(),
            raw_text: None,
            closed: false,
        }
    }

    pub fn from_lines<'a, I>(schema: Schema, lines: I) -> Result<Self, String>
    where
        I: IntoIterator<Item = &'a str>,
    {
        Ok(Self::from_table(
            Table::from_lines(schema, lines)?,
            Flags::READ,
        ))
    }

    pub fn rewind(&mut self) {
        self.cursor = 0;
        self.record_id = 0;
    }

    /// Seek to a textual record boundary.
    ///
    /// The eagerly parsed representation cannot expose a cursor in the middle
    /// of a record. Positions accepted here are therefore the same positions
    /// at which the next upstream `read_record` begins, plus end-of-file.
    pub fn seek(&mut self, position: usize) -> Result<(), String> {
        let text = self.text_bytes();
        if position > text.len()
            || (position != 0
                && position != text.len()
                && text[position - 1] != self.config.line_delimiter)
        {
            return Err("TSV seek position is not a record boundary".to_string());
        }
        self.cursor = if position == text.len() {
            self.table.size()
        } else {
            text[..position]
                .iter()
                .filter(|&&byte| byte == self.config.line_delimiter)
                .count() as i64
        };
        Ok(())
    }

    pub fn eof(&self) -> bool {
        self.cursor >= self.table.size()
    }

    pub fn size(&self) -> i64 {
        self.text_bytes().len() as i64
    }

    pub fn schema(&self) -> Schema {
        self.schema.clone()
    }

    pub fn schema_ref(&self) -> &[Type] {
        &self.schema
    }

    pub fn read_table(&mut self, _threads: usize) -> Table {
        self.rewind();
        self.table.clone()
    }

    pub fn read_chunk(&mut self, max_size: i64, _threads: usize) -> Table {
        let mut out = Table::new(self.schema.clone());
        // `read_text_mt` always performs at least one 1 MiB raw read and only
        // checks `max_size` between full blocks. A record crossing that block
        // boundary is included in full.
        let blocks = (max_size.max(READ_TEXT_MT_SIZE as i64) / READ_TEXT_MT_SIZE as i64).max(1);
        let target = blocks * READ_TEXT_MT_SIZE as i64;
        let mut text_size = 0i64;
        while !self.eof() && text_size < target {
            let record = self.table.record(self.cursor);
            let mut text = TextBuffer::new();
            record.write(&mut text);
            text_size += text.size() as i64;
            out.push_back_record(&record).unwrap();
            self.cursor += 1;
        }
        out
    }

    pub fn read_raw_chunks<F>(&mut self, max_size: i64, threads: usize, callback: F)
    where
        F: FnMut(i64, &[u8]),
    {
        if let Some(text) = &self.raw_text {
            read_text_mt(text.as_bytes(), max_size, threads, callback);
        } else {
            let mut buf = TextBuffer::new();
            self.table.write(&mut buf);
            read_text_mt(buf.data(), max_size, threads, callback);
        }
    }

    pub fn read_tables<F>(&mut self, threads: usize, mut callback: F) -> Result<(), String>
    where
        F: FnMut(i64, &Table),
    {
        let schema = self.schema.clone();
        let line_delimiter = self.config.line_delimiter;
        let mut error = None;
        self.read_raw_chunks(i64::MAX, threads, |chunk, begin| {
            if error.is_some() {
                return;
            }
            let text = String::from_utf8_lossy(begin);
            let mut table = Table::new(schema.clone());
            if let Err(e) = table.append_lines(text.split_terminator(line_delimiter as char)) {
                error = Some(e);
                return;
            }
            callback(chunk, &table);
        });
        if let Some(e) = error {
            Err(e)
        } else {
            Ok(())
        }
    }

    pub fn read_record(&mut self) -> Table {
        let mut table = Table::new(self.schema.clone());
        if self.eof() {
            return table;
        }
        let record = self.table.record(self.cursor);
        table.push_back_record(&record).unwrap();
        self.cursor += 1;
        self.record_id += 1;
        table
    }

    pub fn write_record_values(&mut self, values: &[Value]) -> Result<(), String> {
        let stop = if self.flags.any(Flags::RECORD_ID_COLUMN) {
            self.schema.len().saturating_sub(1)
        } else {
            self.schema.len()
        };
        if values.len() > stop {
            return Err("write_record with too many fields.".to_string());
        }
        if values.len() < stop {
            return Err("write_record with insufficient field count.".to_string());
        }
        if self.flags.any(Flags::RECORD_ID_COLUMN) {
            let mut record = Vec::with_capacity(values.len() + 1);
            record.push(Value::Int64(self.record_id));
            record.extend_from_slice(values);
            self.record_id += 1;
            self.table.write_record(&record)?;
        } else {
            self.table.write_record(values)?;
        }
        self.raw_text = None;
        Ok(())
    }

    pub fn write_record_strings(&mut self, values: &[&str]) -> Result<(), String> {
        let stop = if self.flags.any(Flags::RECORD_ID_COLUMN) {
            self.schema.len().saturating_sub(1)
        } else {
            self.schema.len()
        };
        if values.len() > stop {
            return Err("write_record with too many fields.".to_string());
        }
        if values.len() < stop {
            return Err("write_record with insufficient field count.".to_string());
        }
        self.write_buf.clear();
        for (i, value) in values.iter().enumerate() {
            if i > 0 {
                self.write_buf.append_char('\t');
            }
            self.write_buf.append_str(value);
        }
        self.write_buf
            .append_char(self.config.line_delimiter as char);
        let line = String::from_utf8(self.write_buf.data().to_vec()).map_err(|e| e.to_string())?;
        let record_id = if self.flags.any(Flags::RECORD_ID_COLUMN) {
            let record_id = Some(self.record_id);
            self.record_id += 1;
            record_id
        } else {
            None
        };
        self.table.push_back_line(
            line.trim_end_matches(self.config.line_delimiter as char),
            record_id,
        )?;
        self.raw_text = None;
        Ok(())
    }

    pub fn read_typed<T>(&mut self, out: &mut Vec<T>) -> Result<(), String>
    where
        T: FromTsvRecord,
    {
        if self.schema.len() != T::TYPES.len()
            || self.schema.iter().zip(T::TYPES.iter()).any(|(a, b)| a != b)
        {
            return Err("Template parameters do not match schema.".to_string());
        }
        self.rewind();
        while !self.eof() {
            let record = self.table.record(self.cursor);
            out.push(T::from_record(&record));
            self.cursor += 1;
        }
        Ok(())
    }

    pub fn write(&mut self, record: &Record<'_>) -> Result<(), SchemaMismatch> {
        self.table.push_back_record(record)?;
        self.raw_text = None;
        Ok(())
    }

    pub fn write_table(&mut self, table: &Table) -> Result<(), SchemaMismatch> {
        for i in 0..table.size() {
            self.write(&table.record(i))?;
        }
        Ok(())
    }

    /// Typed-facade counterpart of C++ `write(const TextBuffer*)`.
    pub fn write_text(&mut self, text: &TextBuffer) -> Result<(), String> {
        let text = std::str::from_utf8(text.data()).map_err(|error| error.to_string())?;
        for line in text.split_terminator(self.config.line_delimiter as char) {
            let record_id = self
                .flags
                .any(Flags::RECORD_ID_COLUMN)
                .then_some(self.record_id);
            self.table.push_back_line(line, record_id)?;
            self.record_id += 1;
        }
        self.raw_text = None;
        Ok(())
    }

    pub fn map<F>(&mut self, _threads: usize, mut f: F) -> Self
    where
        F: FnMut(&Record<'_>) -> Table,
    {
        let mut out = Self::from_table(Table::new(self.schema.clone()), Flags::TEMP);
        for i in 0..self.table.size() {
            let mapped = f(&self.table.record(i));
            out.write_table(&mapped).unwrap();
        }
        out
    }

    pub fn sort(&mut self, column: usize, threads: usize) -> Result<Self, String> {
        Ok(Self::from_table(
            self.table.sorted(column, threads)?,
            Flags::TEMP,
        ))
    }

    pub fn close(&mut self) {
        self.closed = true;
    }

    pub fn is_closed(&self) -> bool {
        self.closed
    }

    pub fn file_name(&self) -> String {
        self.file_name.clone()
    }

    pub fn table(&self) -> &Table {
        &self.table
    }

    fn text_bytes(&self) -> Vec<u8> {
        if let Some(text) = &self.raw_text {
            text.as_bytes().to_vec()
        } else {
            let mut text = TextBuffer::new();
            self.table.write(&mut text);
            text.data().to_vec()
        }
    }
}

fn get_flags(mut flags: Flags) -> Result<Flags, String> {
    if flags.any(Flags::TEMP) {
        if flags.any(Flags::WRITE) {
            return Err("Write-only temp file.".to_string());
        }
        flags |= Flags::READ_WRITE | Flags::OVERWRITE;
    }
    Ok(flags)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn size_is_textual_file_size_not_binary_table_allocation() {
        let file = File::from_lines(vec![Type::Int64, Type::String], ["7\tabc", "9\tx"]).unwrap();
        assert_eq!(file.size(), 10);
    }

    #[test]
    fn flag_errors_and_close_state_match_contract() {
        assert_eq!(
            File::new(
                vec![Type::Int64],
                "",
                Flags::TEMP | Flags::WRITE,
                Config::default()
            )
            .unwrap_err(),
            "Write-only temp file."
        );
        assert_eq!(
            File::new(
                vec![Type::String],
                "",
                Flags::RECORD_ID_COLUMN,
                Config::default()
            )
            .unwrap_err(),
            "Schema does not contain record_id column."
        );
        assert_eq!(
            File::new(
                vec![Type::Int64],
                "",
                Flags::WRITE | Flags::OVERWRITE,
                Config::default()
            )
            .unwrap_err(),
            "Invalid File flags"
        );

        let mut file = File::from_lines(vec![Type::Int64], ["1"]).unwrap();
        assert!(!file.is_closed());
        file.close();
        assert!(file.is_closed());
    }

    #[test]
    fn rewind_resets_cursor_and_record_ids() {
        let schema = vec![Type::Int64, Type::Int64];
        let mut file = File::from_table(
            Table::new(schema.clone()),
            Flags::TEMP | Flags::RECORD_ID_COLUMN,
        );
        file.write_record_values(&[Value::Int64(4)]).unwrap();
        file.write_record_values(&[Value::Int64(5)]).unwrap();
        assert_eq!(file.read_record().front().get_t::<i64>(0), 0);
        assert_eq!(file.read_record().front().get_t::<i64>(0), 1);
        assert!(file.eof());
        file.rewind();
        assert_eq!(file.read_record().front().get_t::<i64>(0), 0);
    }

    #[test]
    fn seek_uses_text_offsets_and_does_not_reset_record_id() {
        let schema = vec![Type::Int64, Type::Int64];
        let mut file = File::new(
            schema,
            "",
            Flags::TEMP | Flags::RECORD_ID_COLUMN,
            Config::default(),
        )
        .unwrap();
        file.write_record_strings(&["10"]).unwrap();
        file.write_record_strings(&["20"]).unwrap();

        assert_eq!(file.read_record().front().get_t::<i64>(0), 0);
        // The first rendered record is `0\t10\n`, five bytes long.
        file.seek(5).unwrap();
        assert_eq!(file.read_record().front().get_t::<i64>(0), 1);
        file.write_record_values(&[Value::Int64(30)]).unwrap();
        assert_eq!(file.table().record(2).get_t::<i64>(0), 4);
        assert_eq!(
            file.seek(2).unwrap_err(),
            "TSV seek position is not a record boundary"
        );
    }

    #[test]
    fn small_max_size_still_reads_one_upstream_raw_block() {
        let mut file = File::from_lines(vec![Type::Int64], ["1", "2", "3"]).unwrap();
        let chunk = file.read_chunk(1, 4);
        assert_eq!(chunk.size(), 3);
        assert!(file.eof());
    }

    #[test]
    fn raw_text_write_is_parsed_in_order() {
        let mut file =
            File::from_lines(vec![Type::Int64, Type::String], std::iter::empty::<&str>()).unwrap();
        let mut text = TextBuffer::new();
        text.append_str("2\tb\n1\ta\n");

        file.write_text(&text).unwrap();

        assert_eq!(file.table().size(), 2);
        assert_eq!(file.table().record(0).get_t::<i64>(0), 2);
        assert_eq!(file.table().record(1).get_t::<String>(1), "a");
    }
}
