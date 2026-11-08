//! Typed TSV file facade and chunked text reads.
//!
//! Path constructors retain Source Input/Output owners and defer parsing.
//! Explicit in-memory constructors preserve the reusable table facade.

use crate::util::io::{Compressor, InputFile, IoErrorOwner, OutputFile, TempFile};
use crate::util::text_buffer::TextBuffer;
use std::cell::RefCell;
use std::io;
use std::path::{Path, PathBuf};

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

#[derive(Debug)]
pub struct File {
    flags: Flags,
    schema: Schema,
    config: Config,
    table: Table,
    cursor: i64,
    write_buf: TextBuffer,
    record_id: RecordId,
    file_name: PathBuf,
    raw_text: Option<Vec<u8>>,
    closed: bool,
    source: Option<RefCell<SourceOwner>>,
}

impl File {
    pub fn new(
        schema: Schema,
        file_name: &str,
        flags: Flags,
        config: Config,
    ) -> Result<Self, String> {
        Self::new_path_io(schema, Path::new(file_name), flags, config)
            .map_err(|error| error.to_string())
    }

    pub fn new_path(
        schema: Schema,
        file_name: &Path,
        flags: Flags,
        config: Config,
    ) -> Result<Self, String> {
        Self::new_path_io(schema, file_name, flags, config).map_err(|error| error.to_string())
    }

    pub fn new_path_io(
        schema: Schema,
        file_name: &Path,
        flags: Flags,
        config: Config,
    ) -> io::Result<Self> {
        Self::new_path_with_buffer_size_io(schema, file_name, flags, config, 67_108_864)
    }

    pub fn new_path_with_buffer_size_io(
        schema: Schema,
        file_name: &Path,
        flags: Flags,
        config: Config,
        buffer_size: usize,
    ) -> io::Result<Self> {
        let flags = get_flags(flags).map_err(io::Error::other)?;
        if flags.all(Flags::WRITE | Flags::OVERWRITE) {
            return Err(io::Error::other("Invalid File flags"));
        }
        if flags.any(Flags::RECORD_ID_COLUMN) && schema.first() != Some(&Type::Int64) {
            return Err(io::Error::other(
                "Schema does not contain record_id column.",
            ));
        }
        let source = SourceOwner::new(file_name, flags, buffer_size)?;
        let actual_name = source
            .input
            .as_ref()
            .map(InputFile::file_path)
            .or_else(|| source.output.as_ref().map(SourceOutput::path))
            .unwrap_or(file_name)
            .to_path_buf();
        Ok(Self {
            table: Table::new(schema.clone()),
            flags,
            schema,
            config,
            cursor: 0,
            write_buf: TextBuffer::new(),
            record_id: 0,
            file_name: actual_name,
            raw_text: None,
            closed: false,
            source: Some(RefCell::new(source)),
        })
    }

    /// Source File overload adopting an already constructed TextInputFile.
    /// Does not open again, inspect size, seek or change the supplied buffers.
    pub fn from_source_input_io(
        schema: Schema,
        input: InputFile,
        flags: Flags,
        config: Config,
    ) -> io::Result<Self> {
        let normalized = get_flags(flags).map_err(io::Error::other)?;
        if flags != Flags::READ {
            return Err(io::Error::other("File::File"));
        }
        let file_name = input.file_name.clone();
        Ok(Self {
            flags: normalized,
            schema: schema.clone(),
            config,
            table: Table::new(schema),
            cursor: 0,
            write_buf: TextBuffer::new(),
            record_id: 0,
            file_name,
            raw_text: None,
            closed: false,
            source: Some(RefCell::new(SourceOwner {
                input: Some(input),
                output: None,
                errors: IoErrorOwner::default().with_source_destructors(),
            })),
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
            file_name: PathBuf::new(),
            raw_text: None,
            closed: false,
            source: None,
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

    pub fn from_byte_lines<'a>(
        schema: Schema,
        lines: impl IntoIterator<Item = &'a [u8]>,
    ) -> io::Result<Self> {
        let mut table = Table::new(schema);
        for line in lines {
            table.push_back_raw_line_io(line, None)?;
        }
        Ok(Self::from_table(table, Flags::READ))
    }

    pub fn rewind_io(&mut self) -> io::Result<()> {
        if let Some(source) = &self.source {
            let mut source = source.borrow_mut();
            let result = source.rewind();
            source.errors.capture_io(result)?;
        }
        self.cursor = 0;
        self.record_id = 0;
        Ok(())
    }

    pub fn rewind(&mut self) {
        self.rewind_io().expect("TSV rewind failed; use rewind_io");
        self.cursor = 0;
        self.record_id = 0;
    }

    /// Seek to a textual record boundary.
    ///
    /// Source path owners forward arbitrary byte positions to InputFile.
    /// The explicit in-memory facade accepts only rendered record boundaries.
    pub fn seek_io(&mut self, position: usize) -> io::Result<()> {
        if let Some(source) = &self.source {
            let mut source = source.borrow_mut();
            let result = source
                .input_mut()?
                .seek(position as i64)
                .map(|_| ())
                .map_err(io::Error::other);
            return source.errors.capture_io(result);
        }
        self.seek_memory(position).map_err(io::Error::other)
    }

    pub fn seek(&mut self, position: usize) -> Result<(), String> {
        self.seek_io(position).map_err(|error| error.to_string())
    }

    fn seek_memory(&mut self, position: usize) -> Result<(), String> {
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
        if self.source.is_some() {
            // Source TextInputFile's eof_ is changed only by getline; its
            // default CharTokenizer never calls getline (read_record throws).
            return false;
        }
        self.cursor >= self.table.size()
    }

    pub fn size_io(&self) -> io::Result<i64> {
        if let Some(source) = &self.source {
            let mut source = source.borrow_mut();
            let result = source.input_mut()?.file_size().map_err(io::Error::other);
            return source.errors.capture_io(result);
        }
        Ok(self.text_bytes().len() as i64)
    }

    pub fn size(&self) -> i64 {
        self.size_io().expect("TSV size failed; use size_io")
    }

    pub fn schema(&self) -> Schema {
        self.schema.clone()
    }

    pub fn is_source_owner(&self) -> bool {
        self.source.is_some()
    }

    pub fn schema_ref(&self) -> &[Type] {
        &self.schema
    }

    pub fn read_table_io(&mut self, threads: usize) -> io::Result<Table> {
        self.rewind_io()?;
        if self.source.is_some() {
            self.read_chunk_io(i64::MAX, threads)
        } else {
            Ok(self.table.clone())
        }
    }
    pub fn read_table(&mut self, threads: usize) -> Table {
        self.read_table_io(threads)
            .expect("TSV read failed; use read_table_io")
    }

    pub fn read_chunk_io(&mut self, max_size: i64, threads: usize) -> io::Result<Table> {
        if let Some(source) = &self.source {
            let mut source = source.borrow_mut();
            let result = super::streaming::StreamingFile::read_input_io(
                source.input_mut()?,
                self.schema.clone(),
                max_size,
                threads,
            );
            return source.errors.capture_io(result);
        }
        Ok(self.read_chunk_memory(max_size, threads))
    }
    pub fn read_chunk(&mut self, max_size: i64, threads: usize) -> Table {
        self.read_chunk_io(max_size, threads)
            .expect("TSV chunk read failed; use read_chunk_io")
    }
    fn read_chunk_memory(&mut self, max_size: i64, _threads: usize) -> Table {
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
        self.read_raw_chunks_io(max_size, threads, callback)
            .expect("TSV raw read failed; use read_raw_chunks_io")
    }
    /// Ordered FnMut compatibility; Source raw callbacks use the parallel
    /// overload, while this service keeps a non-Send callback on its caller.
    pub fn read_raw_chunks_io<F>(
        &mut self,
        max_size: i64,
        threads: usize,
        mut callback: F,
    ) -> io::Result<()>
    where
        F: FnMut(i64, &[u8]),
    {
        if let Some(source) = &self.source {
            let mut source = source.borrow_mut();
            let result = source.read_raw_chunks(max_size, threads, &mut callback);
            return source.errors.capture_io(result);
        }
        if let Some(text) = &self.raw_text {
            read_text_mt(text, max_size, threads, callback);
        } else {
            let mut buf = TextBuffer::new();
            self.table.write(&mut buf);
            read_text_mt(buf.data(), max_size, threads, callback);
        }
        Ok(())
    }

    /// Canonical Source callback overload. Callbacks run concurrently in the
    /// consumer workers; chunk indices follow FIFO dequeue order.
    pub fn read_raw_chunks_parallel_io<F>(
        &mut self,
        max_size: i64,
        threads: usize,
        callback: F,
    ) -> io::Result<()>
    where
        F: Fn(i64, &[u8]) -> io::Result<()> + Sync,
    {
        if let Some(source) = &self.source {
            let mut source = source.borrow_mut();
            let result = super::streaming::read_input_parallel_io(
                source.input_mut()?,
                max_size,
                threads,
                true,
                callback,
            );
            return source.errors.capture_io(result);
        }
        let bytes = self.text_bytes();
        super::streaming::read_bytes_parallel_io(&bytes, max_size, threads, callback)
    }

    pub fn read_tables_parallel_io<F>(&mut self, threads: usize, callback: F) -> io::Result<()>
    where
        F: Fn(i64, &Table) -> io::Result<()> + Sync,
    {
        let schema = self.schema.clone();
        let delimiter = self.config.line_delimiter;
        self.read_raw_chunks_parallel_io(i64::MAX, threads, |index, bytes| {
            let mut table = Table::new(schema.clone());
            for line in raw_lines(bytes, delimiter) {
                table.push_back_raw_line_io(line, None)?;
            }
            callback(index, &table)
        })
    }

    /// Mutable-callback compatibility. Parsing and callbacks are ordered and
    /// serial; use read_tables_parallel_io for Source's concurrent callbacks.
    pub fn read_tables<F>(&mut self, threads: usize, callback: F) -> Result<(), String>
    where
        F: FnMut(i64, &Table),
    {
        self.read_tables_io(threads, callback)
            .map_err(|error| error.to_string())
    }
    pub fn read_tables_io<F>(&mut self, threads: usize, mut callback: F) -> io::Result<()>
    where
        F: FnMut(i64, &Table),
    {
        let schema = self.schema.clone();
        let line_delimiter = self.config.line_delimiter;
        let source_profile = self.source.is_some();
        let mut error = None;
        self.read_raw_chunks_io(i64::MAX, threads, |chunk, begin| {
            if error.is_some() {
                return;
            }
            let mut table = Table::new(schema.clone());
            for line in raw_lines(begin, line_delimiter) {
                if let Err(e) = table.push_back_raw_line_io(line, None) {
                    if source_profile {
                        crate::util::io::abort_source_process();
                    }
                    error = Some(e);
                    return;
                }
            }
            callback(chunk, &table);
        })?;
        error.map_or(Ok(()), Err)
    }

    pub fn read_record_io(&mut self) -> io::Result<Table> {
        if self.source.is_some() {
            return Err(io::Error::other("Invalid use of CharTokenizer"));
        }
        Ok(self.read_record_memory())
    }
    pub fn read_record(&mut self) -> Table {
        self.read_record_io()
            .expect("TSV read_record failed; use read_record_io")
    }
    fn read_record_memory(&mut self) -> Table {
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
        self.write_record_values_io(values)
            .map_err(|error| error.to_string())
    }
    pub fn write_record_values_io(&mut self, values: &[Value]) -> io::Result<()> {
        if self.source.is_some() {
            let stop = self
                .schema
                .len()
                .saturating_sub(usize::from(self.flags.any(Flags::RECORD_ID_COLUMN)));
            if values.len() != stop {
                return Err(io::Error::other(if values.len() > stop {
                    "write_record with too many fields."
                } else {
                    "write_record with insufficient field count."
                }));
            }
            let mut bytes = TextBuffer::new();
            for (i, value) in values.iter().enumerate() {
                if i > 0 {
                    bytes.write_raw(b"\t");
                }
                match value {
                    Value::String(value) => bytes.write_raw(value.as_bytes()),
                    Value::Bytes(value) => bytes.write_raw(value),
                    Value::Int64(value) => {
                        bytes.append_display(value);
                    }
                };
            }
            bytes.write_raw(b"\n");
            return self.source_write_io(bytes.data());
        }
        self.write_record_values_memory(values)
            .map_err(io::Error::other)
    }
    fn write_record_values_memory(&mut self, values: &[Value]) -> Result<(), String> {
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
        let bytes: Vec<&[u8]> = values.iter().map(|value| value.as_bytes()).collect();
        self.write_record_bytes(&bytes)
    }
    pub fn write_record_bytes(&mut self, values: &[&[u8]]) -> Result<(), String> {
        self.write_record_bytes_io(values)
            .map_err(|error| error.to_string())
    }
    pub fn write_record_bytes_io(&mut self, values: &[&[u8]]) -> io::Result<()> {
        if self.source.is_some() {
            let stop = self
                .schema
                .len()
                .saturating_sub(usize::from(self.flags.any(Flags::RECORD_ID_COLUMN)));
            if values.len() != stop {
                return Err(io::Error::other(if values.len() > stop {
                    "write_record with too many fields."
                } else {
                    "write_record with insufficient field count."
                }));
            }
            let mut bytes = Vec::new();
            for (i, value) in values.iter().enumerate() {
                if i > 0 {
                    bytes.push(b'\t');
                }
                bytes.extend_from_slice(value);
            }
            bytes.push(b'\n');
            return self.source_write_io(&bytes);
        }
        self.write_record_bytes_memory(values)
            .map_err(io::Error::other)
    }
    fn write_record_bytes_memory(&mut self, values: &[&[u8]]) -> Result<(), String> {
        let stop = self
            .schema
            .len()
            .saturating_sub(usize::from(self.flags.any(Flags::RECORD_ID_COLUMN)));
        if values.len() != stop {
            return Err(if values.len() > stop {
                "write_record with too many fields."
            } else {
                "write_record with insufficient field count."
            }
            .into());
        }
        self.write_buf.clear();
        for (i, value) in values.iter().enumerate() {
            if i > 0 {
                self.write_buf.append_char('\t');
            }
            self.write_buf.write_raw(value);
        }
        self.write_buf
            .append_char(self.config.line_delimiter as char);
        let mut line = self.write_buf.data();
        while line.last() == Some(&self.config.line_delimiter) {
            line = &line[..line.len() - 1];
        }
        let record_id = if self.flags.any(Flags::RECORD_ID_COLUMN) {
            let record_id = Some(self.record_id);
            self.record_id += 1;
            record_id
        } else {
            None
        };
        self.table.push_back_raw_line(line, record_id)?;
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
        self.rewind_io().map_err(|error| error.to_string())?;
        loop {
            let table = self.read_record_io().map_err(|error| error.to_string())?;
            if table.empty() {
                break;
            }
            out.push(T::from_record(&table.front()));
        }
        Ok(())
    }

    pub fn write_io(&mut self, record: &Record<'_>) -> io::Result<()> {
        if self.source.is_some() {
            let mut bytes = TextBuffer::new();
            record.write(&mut bytes);
            return self.source_write_io(bytes.data());
        }
        self.table
            .push_back_record(record)
            .map_err(io::Error::other)?;
        self.raw_text = None;
        Ok(())
    }
    pub fn write(&mut self, record: &Record<'_>) -> Result<(), SchemaMismatch> {
        if self.source.is_some() {
            self.write_io(record)
                .expect("TSV write failed; use write_io");
            return Ok(());
        }
        self.table.push_back_record(record)?;
        self.raw_text = None;
        Ok(())
    }

    pub fn write_table_io(&mut self, table: &Table) -> io::Result<()> {
        for i in 0..table.size() {
            self.write_io(&table.record(i))?;
        }
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
        self.write_text_bytes_io(text.data())
            .map_err(|error| error.to_string())
    }

    pub fn write_text_bytes_io(&mut self, text: &[u8]) -> io::Result<()> {
        if self.source.is_some() {
            return self.source_write_io(text);
        }
        for line in raw_lines(text, self.config.line_delimiter) {
            let record_id = self
                .flags
                .any(Flags::RECORD_ID_COLUMN)
                .then_some(self.record_id);
            self.table.push_back_raw_line_io(line, record_id)?;
            self.record_id += 1;
        }
        self.raw_text = None;
        Ok(())
    }

    /// Source map callback overload. Per-chunk transforms execute in the
    /// workers; the real output owner receives completed chunks in index order.
    pub fn map_parallel_io<F>(&mut self, threads: usize, f: F) -> io::Result<Self>
    where
        F: Fn(&Record<'_>) -> Table + Sync,
    {
        if self.source.is_some() {
            let out = Self::new_path_io(
                self.schema.clone(),
                Path::new(""),
                Flags::TEMP,
                Config::default(),
            )?;
            let output = out
                .source
                .as_ref()
                .unwrap()
                .borrow_mut()
                .output
                .take()
                .ok_or_else(|| io::Error::other("TSV file is not writable"))?;
            struct OrderedOutput {
                next: i64,
                pending: std::collections::BTreeMap<i64, Vec<u8>>,
                output: SourceOutput,
            }
            let ordered = std::sync::Mutex::new(OrderedOutput {
                next: 0,
                pending: Default::default(),
                output,
            });
            let result = self.read_tables_parallel_io(threads, |index, table| {
                let mut bytes = TextBuffer::new();
                for row in 0..table.size() {
                    f(&table.record(row)).write(&mut bytes);
                }
                let mut ordered = ordered
                    .lock()
                    .map_err(|_| io::Error::other("TSV output reorder queue poisoned"))?;
                ordered.pending.insert(index, bytes.data().to_vec());
                loop {
                    let next = ordered.next;
                    let Some(bytes) = ordered.pending.remove(&next) else {
                        break;
                    };
                    ordered.output.write(&bytes)?;
                    ordered.next += 1;
                }
                Ok(())
            });
            let ordered = ordered
                .into_inner()
                .map_err(|_| io::Error::other("TSV output reorder queue poisoned"))?;
            out.source.as_ref().unwrap().borrow_mut().output = Some(ordered.output);
            result?;
            return Ok(out);
        }
        struct OrderedTables {
            next: i64,
            pending: std::collections::BTreeMap<i64, Table>,
            table: Table,
        }
        let ordered = std::sync::Mutex::new(OrderedTables {
            next: 0,
            pending: Default::default(),
            table: Table::new(self.schema.clone()),
        });
        self.read_tables_parallel_io(threads, |index, table| {
            let mut mapped = Table::new(table.schema().to_vec());
            for row in 0..table.size() {
                mapped
                    .append(&f(&table.record(row)))
                    .map_err(io::Error::other)?;
            }
            let mut ordered = ordered
                .lock()
                .map_err(|_| io::Error::other("TSV output reorder queue poisoned"))?;
            ordered.pending.insert(index, mapped);
            loop {
                let next = ordered.next;
                let Some(table) = ordered.pending.remove(&next) else {
                    break;
                };
                ordered.table.append(&table).map_err(io::Error::other)?;
                ordered.next += 1;
            }
            Ok(())
        })?;
        Ok(Self::from_table(
            ordered
                .into_inner()
                .map_err(|_| io::Error::other("TSV output reorder queue poisoned"))?
                .table,
            Flags::TEMP,
        ))
    }

    /// FnMut compatibility facade; map_parallel_io represents Source's
    /// concurrent const-callable callback and ordered output owner.
    pub fn map_io<F>(&mut self, threads: usize, mut f: F) -> io::Result<Self>
    where
        F: FnMut(&Record<'_>) -> Table,
    {
        if self.source.is_some() {
            let mut out = Self::new_path_io(
                self.schema.clone(),
                Path::new(""),
                Flags::TEMP,
                Config::default(),
            )?;
            let mut error = None;
            self.read_tables_io(threads, |_, table| {
                if error.is_some() {
                    return;
                }
                for i in 0..table.size() {
                    if let Err(e) = out.write_source_table(&f(&table.record(i))) {
                        error = Some(e);
                        break;
                    }
                }
            })?;
            if let Some(error) = error {
                return Err(error);
            }
            return Ok(out);
        }
        let mut out = Self::from_table(Table::new(self.schema.clone()), Flags::TEMP);
        for i in 0..self.table.size() {
            out.write_table(&f(&self.table.record(i)))
                .map_err(io::Error::other)?;
        }
        Ok(out)
    }
    pub fn map<F>(&mut self, threads: usize, f: F) -> Self
    where
        F: FnMut(&Record<'_>) -> Table,
    {
        self.map_io(threads, f).expect("TSV map failed; use map_io")
    }

    pub fn sort_io(&mut self, column: usize, threads: usize) -> io::Result<Self> {
        if self.source.is_some() {
            // Source's 1,000,000 parameter is a raw-byte read limit, while
            // its loop continuation compares the resulting ROW count.
            self.rewind_io()?;
            let mut files = Vec::new();
            loop {
                let table = self.read_chunk_io(1_000_000, threads)?;
                if table.empty() {
                    break;
                }
                let mut out = Self::new_path_io(
                    self.schema.clone(),
                    Path::new(""),
                    Flags::TEMP,
                    Config::default(),
                )?;
                out.write_source_table(&table.sorted(column, threads).map_err(io::Error::other)?)?;
                out.rewind_io()?;
                let again = table.size() == 1_000_000;
                files.push(out);
                if !again {
                    break;
                }
            }
            let mut refs: Vec<_> = files.iter_mut().collect();
            return super::merge::merge_io(&mut refs, column);
        }
        Ok(Self::from_table(
            self.table
                .sorted(column, threads)
                .map_err(io::Error::other)?,
            Flags::TEMP,
        ))
    }
    pub fn sort(&mut self, column: usize, threads: usize) -> Result<Self, String> {
        self.sort_io(column, threads)
            .map_err(|error| error.to_string())
    }

    pub fn close_io(&mut self) -> io::Result<()> {
        if let Some(source) = &self.source {
            let mut source = source.borrow_mut();
            let result = source.close();
            source.errors.capture_io(result)?;
        }
        self.closed = true;
        Ok(())
    }
    pub fn close(&mut self) {
        self.close_io().expect("TSV close failed; use close_io");
    }

    fn source_write_io(&mut self, bytes: &[u8]) -> io::Result<()> {
        let mut source = self.source.as_ref().expect("Source TSV owner").borrow_mut();
        let result = source
            .output
            .as_mut()
            .ok_or_else(|| io::Error::other("TSV file is not writable"))?
            .write(bytes);
        source.errors.capture_io(result)
    }
    fn write_source_table(&mut self, table: &Table) -> io::Result<()> {
        let mut bytes = TextBuffer::new();
        table.write(&mut bytes);
        self.source_write_io(bytes.data())
    }

    pub fn is_closed(&self) -> bool {
        self.closed
    }

    pub fn file_path(&self) -> &Path {
        &self.file_name
    }

    pub fn try_file_name(&self) -> io::Result<String> {
        self.file_name
            .to_str()
            .map(str::to_owned)
            .ok_or_else(|| io::Error::other("Non-UTF-8 TSV path; use file_path"))
    }

    pub fn file_name(&self) -> String {
        self.try_file_name()
            .expect("Non-UTF-8 TSV path; use file_path")
    }

    /// The explicit in-memory representation. Source owners are read through
    /// read_table_io, never an eagerly populated duplicate table.
    pub fn table(&self) -> &Table {
        assert!(
            self.source.is_none(),
            "Path-backed TSV owner; use read_table_io"
        );
        &self.table
    }

    fn text_bytes(&self) -> Vec<u8> {
        if let Some(text) = &self.raw_text {
            text.clone()
        } else {
            let mut text = TextBuffer::new();
            self.table.write(&mut text);
            text.data().to_vec()
        }
    }
}
#[derive(Debug)]
enum SourceOutput {
    Named(OutputFile),
    Temp(TempFile),
}
impl SourceOutput {
    fn path(&self) -> &Path {
        match self {
            Self::Named(out) => out.file_path(),
            Self::Temp(out) => out.file_path(),
        }
    }
    fn write(&mut self, bytes: &[u8]) -> io::Result<()> {
        match self {
            Self::Named(out) => out.write_raw(bytes),
            Self::Temp(out) => out.write_raw(bytes),
        }
        .map_err(io::Error::other)
    }
    fn rewind(&mut self) -> io::Result<()> {
        match self {
            Self::Named(out) => out.rewind(),
            Self::Temp(out) => out.rewind(),
        }
        .map_err(io::Error::other)
    }
    fn close(&mut self) -> io::Result<()> {
        match self {
            Self::Named(out) => out.close(),
            Self::Temp(out) => out.close(),
        }
        .map_err(io::Error::other)
    }
}
#[derive(Debug)]
struct SourceOwner {
    input: Option<InputFile>,
    output: Option<SourceOutput>,
    errors: IoErrorOwner,
}
impl SourceOwner {
    fn new(path: &Path, flags: Flags, buffer_size: usize) -> io::Result<Self> {
        let mut output = if flags.any(Flags::WRITE | Flags::READ_WRITE) {
            Some(if flags.any(Flags::TEMP) {
                SourceOutput::Temp({
                    let data = TempFile::init(true).map_err(io::Error::other)?;
                    TempFile::from_temp_file_data_with_buffer_size(&data, buffer_size)
                        .map_err(io::Error::other)?
                })
            } else {
                SourceOutput::Named(
                    OutputFile::with_source_profile_path(
                        path,
                        Compressor::None,
                        if flags.any(Flags::OVERWRITE) {
                            "w+b"
                        } else if flags.any(Flags::READ_WRITE) {
                            "r+b"
                        } else {
                            "wb"
                        },
                        buffer_size,
                    )
                    .map_err(io::Error::other)?,
                )
            })
        } else {
            None
        };
        let input = if flags.any(Flags::READ_WRITE) {
            Some(
                match output.as_mut().expect("TSV read/write output") {
                    SourceOutput::Named(out) => {
                        InputFile::from_output_file_with_buffer_size(out, 0, buffer_size)
                    }
                    SourceOutput::Temp(out) => InputFile::from_temp_file_with_buffer_size(
                        out,
                        0,
                        Compressor::None,
                        buffer_size,
                    ),
                }
                .map_err(io::Error::other)?,
            )
        } else if !flags.any(Flags::WRITE) {
            Some(
                InputFile::with_source_profile_path(path, 0, buffer_size)
                    .map_err(io::Error::other)?,
            )
        } else {
            None
        };
        Ok(Self {
            input,
            output,
            errors: IoErrorOwner::default().with_source_destructors(),
        })
    }
    fn input_mut(&mut self) -> io::Result<&mut InputFile> {
        self.input
            .as_mut()
            .ok_or_else(|| io::Error::other("TSV file is not readable"))
    }
    fn rewind(&mut self) -> io::Result<()> {
        if let Some(output) = self.output.as_mut() {
            output.rewind()?;
        }
        if let Some(input) = self.input.as_mut() {
            input.rewind().map_err(io::Error::other)?;
        }
        Ok(())
    }
    fn read_raw_chunks(
        &mut self,
        max_size: i64,
        threads: usize,
        callback: &mut impl FnMut(i64, &[u8]),
    ) -> io::Result<()> {
        super::streaming::read_input_serial_io(self.input_mut()?, max_size, threads, true, callback)
    }
    fn close(&mut self) -> io::Result<()> {
        if self.input.as_ref().is_some_and(|input| input.temp_file) {
            if let Some(input) = self.input.as_mut() {
                input.close_and_delete().map_err(io::Error::other)?;
            }
        } else if let Some(output) = self.output.as_mut() {
            output.close()?;
        } else if let Some(input) = self.input.as_mut() {
            input.close().map_err(io::Error::other)?;
        }
        self.input = None;
        self.output = None;
        Ok(())
    }
}
impl Drop for SourceOwner {
    fn drop(&mut self) {
        if let Err(error) = self.close() {
            if let Err(error) = self.errors.emit_pending_perror() {
                crate::util::io::abort_source_diagnostic_error(&error);
            }
            if let Err(error) = crate::util::io::emit_source_perror(&error) {
                crate::util::io::abort_source_diagnostic_error(&error);
            }
            crate::util::io::abort_source_process();
        }
    }
}

// Byte equivalent of str::split_terminator: retain interior empty rows,
// omit only the empty row following a final delimiter.
fn raw_lines(bytes: &[u8], delimiter: u8) -> impl Iterator<Item = &[u8]> {
    let end = bytes.len() - usize::from(bytes.last() == Some(&delimiter));
    bytes[..end]
        .split(move |&byte| byte == delimiter)
        .filter(move |_| !bytes.is_empty())
        .map(|line| line.strip_suffix(b"\r").unwrap_or(line))
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
    fn source_path_owner_defers_parsing_and_accepts_byte_position_seek() {
        let path = std::env::temp_dir().join(format!(
            "diamond-tsv-owned-{}-{}.tsv",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, b"header\xff\ttail\n").unwrap();
        let mut file = File::new_path_io(
            vec![Type::String, Type::String],
            &path,
            Flags::READ,
            Config::default(),
        )
        .unwrap();
        assert_eq!(file.file_path(), path);
        assert_eq!(
            file.read_record_io().unwrap_err().to_string(),
            "Invalid use of CharTokenizer"
        );
        file.seek_io(2).unwrap();
        let table = file.read_chunk_io(i64::MAX, 2).unwrap();
        assert_eq!(table.size(), 1);
        assert_eq!(table.record(0).string_bytes(0), b"ader\xff");
        assert_eq!(table.record(0).string_bytes(1), b"tail");
        assert!(!file.eof());
        file.close_io().unwrap();
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn source_temp_writer_rewinds_retained_alias_and_does_not_parse_on_write() {
        let mut file = File::new_path_io(
            vec![Type::String],
            Path::new(""),
            Flags::TEMP,
            Config::default(),
        )
        .unwrap();
        file.write_record_bytes(&[b"id\xff"]).unwrap();
        // Serializer bytes are still pending, just as in Source TSV TempFile.
        assert_eq!(file.size_io().unwrap(), 0);
        let table = file.read_table_io(1).unwrap();
        assert_eq!(table.record(0).string_bytes(0), b"id\xff");
        file.close_io().unwrap();

        let mut unchecked = File::new_path_io(
            vec![Type::Int64],
            Path::new(""),
            Flags::TEMP,
            Config::default(),
        )
        .unwrap();
        unchecked
            .write_record_bytes(&[b"not-an-integer\xff"])
            .unwrap();
        unchecked.close_io().unwrap();
    }

    #[test]
    fn source_named_read_write_alias_close_deletes_without_flushing_pending_text() {
        let path = std::env::temp_dir().join(format!(
            "diamond-tsv-named-{}-{}.tsv",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut file = File::new_path_io(
            vec![Type::String],
            &path,
            Flags::READ_WRITE | Flags::OVERWRITE,
            Config::default(),
        )
        .unwrap();
        file.write_record_bytes(&[b"pending\xff"]).unwrap();
        assert_eq!(std::fs::metadata(&path).unwrap().len(), 0);
        file.close_io().unwrap();
        assert!(!path.exists());
    }

    #[test]
    fn raw_memory_read_tables_and_write_text_keep_byte_fields() {
        let mut file = File::from_byte_lines(
            vec![Type::String, Type::String],
            [b"q\xff\ts\xfe".as_slice()],
        )
        .unwrap();
        let mut rows = Vec::new();
        file.read_tables_io(1, |_, table| rows.push(table.record(0).get_t::<Vec<u8>>(0)))
            .unwrap();
        assert_eq!(rows, [b"q\xff".to_vec()]);
        file.write_text_bytes_io(b"x\xff\ty\xfe\n").unwrap();
        assert_eq!(file.table().record(1).string_bytes(1), b"y\xfe");
        assert!(file.table().record(1).try_get_string(1).is_err());
    }

    #[test]
    fn byte_record_write_preserves_non_utf8_through_table_serialization() {
        let mut file = File::from_table(Table::new(vec![Type::String, Type::String]), Flags::TEMP);
        file.write_record_bytes(&[b"q\xff", b"s\xfe"]).unwrap();
        assert_eq!(file.table().record(0).string_bytes(0), b"q\xff");
        let mut output = TextBuffer::new();
        file.table().write(&mut output);
        assert_eq!(output.data(), b"q\xff\ts\xfe\n");
    }

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
        let mut file = File::from_table(Table::new(schema), Flags::TEMP | Flags::RECORD_ID_COLUMN);
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
