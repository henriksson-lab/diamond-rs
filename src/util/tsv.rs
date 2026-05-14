use crate::util::string::convert_string_i64;
use crate::util::text_buffer::TextBuffer;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

pub type Schema = Vec<Type>;
pub type RecordId = i64;

const READ_TEXT_MT_SIZE: usize = 1 << 20;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Type {
    String,
    Int64,
}

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

pub trait TsvValue: Sized {
    fn interpret(ptr: &[u8]) -> Self;
    fn push_to(self, table: &mut Table);
}

impl TsvValue for i64 {
    fn interpret(ptr: &[u8]) -> Self {
        i64::from_ne_bytes(ptr[..8].try_into().unwrap())
    }

    fn push_to(self, table: &mut Table) {
        table.push_i64(self);
    }
}

impl TsvValue for i32 {
    fn interpret(ptr: &[u8]) -> Self {
        i32::from_ne_bytes(ptr[..4].try_into().unwrap())
    }

    fn push_to(self, table: &mut Table) {
        table.push_i32(self);
    }
}

impl TsvValue for String {
    fn interpret(ptr: &[u8]) -> Self {
        let len = i32::from_ne_bytes(ptr[..4].try_into().unwrap()) as usize;
        String::from_utf8(ptr[4..4 + len].to_vec()).unwrap()
    }

    fn push_to(self, table: &mut Table) {
        table.push_str(&self);
    }
}

impl TsvValue for &str {
    fn interpret(_: &[u8]) -> Self {
        panic!("&str cannot be read from a TSV record")
    }

    fn push_to(self, table: &mut Table) {
        table.push_str(self);
    }
}

pub trait FromTsvRecord: Sized {
    const TYPES: &'static [Type];

    fn from_record(record: &Record<'_>) -> Self;
}

impl FromTsvRecord for () {
    const TYPES: &'static [Type] = &[];

    fn from_record(_: &Record<'_>) -> Self {}
}

macro_rules! impl_from_tsv_record_tuple {
    ($($name:ident : $idx:tt),+) => {
        impl<$($name),+> FromTsvRecord for ($($name,)+)
        where
            $($name: TsvValue + TsvType,)+
        {
            const TYPES: &'static [Type] = &[$($name::TYPE,)+];

            fn from_record(record: &Record<'_>) -> Self {
                ($(record.get_t::<$name>($idx),)+)
            }
        }
    };
}

pub trait TsvType {
    const TYPE: Type;
}

impl TsvType for i64 {
    const TYPE: Type = Type::Int64;
}

impl TsvType for String {
    const TYPE: Type = Type::String;
}

impl_from_tsv_record_tuple!(A: 0);
impl_from_tsv_record_tuple!(A: 0, B: 1);
impl_from_tsv_record_tuple!(A: 0, B: 1, C: 2);
impl_from_tsv_record_tuple!(A: 0, B: 1, C: 2, D: 3);
impl_from_tsv_record_tuple!(A: 0, B: 1, C: 2, D: 3, E: 4);
impl_from_tsv_record_tuple!(A: 0, B: 1, C: 2, D: 3, E: 4, F: 5);

#[derive(Debug, Clone, Copy)]
pub struct Record<'a> {
    schema: &'a [Type],
    buf: &'a [u8],
}

impl<'a> Record<'a> {
    pub fn new(schema: &'a [Type], begin: &'a [u8], end: usize) -> Self {
        Self {
            schema,
            buf: &begin[..end],
        }
    }

    pub fn get_t<T: TsvValue>(&self, i: usize) -> T {
        let mut it = self.begin();
        for _ in 0..i {
            it.next();
        }
        it.get()
    }

    pub fn get(&self, i: usize) -> String {
        let mut it = self.begin();
        for _ in 0..i {
            it.next();
        }
        it.value_string()
    }

    pub fn begin(&self) -> RecordIterator<'a> {
        RecordIterator {
            schema: self.schema,
            buf: self.buf,
            idx: 0,
            ptr: 0,
        }
    }

    pub fn end(&self) -> RecordIterator<'a> {
        RecordIterator {
            schema: self.schema,
            buf: self.buf,
            idx: self.schema.len(),
            ptr: self.buf.len(),
        }
    }

    pub fn raw_size(&self) -> i64 {
        self.buf.len() as i64
    }

    pub fn write(&self, buf: &mut TextBuffer) {
        let mut n = 0usize;
        for field in self.begin() {
            if n > 0 {
                buf.append_char('\t');
            }
            match field.type_() {
                Type::Int64 => {
                    buf.append_display(field.get::<i64>());
                }
                Type::String => {
                    buf.append_str(&field.get::<String>());
                }
            }
            n += 1;
        }
        buf.append_char('\n');
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RecordField<'a> {
    type_: Type,
    ptr: &'a [u8],
}

impl<'a> RecordField<'a> {
    pub fn type_(&self) -> Type {
        self.type_
    }

    pub fn get<T: TsvValue>(&self) -> T {
        T::interpret(self.ptr)
    }

    pub fn value_string(&self) -> String {
        match self.type_ {
            Type::String => self.get::<String>(),
            Type::Int64 => self.get::<i64>().to_string(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RecordIterator<'a> {
    schema: &'a [Type],
    buf: &'a [u8],
    idx: usize,
    ptr: usize,
}

impl<'a> RecordIterator<'a> {
    pub fn type_(&self) -> Type {
        self.schema[self.idx]
    }

    pub fn get<T: TsvValue>(&self) -> T {
        T::interpret(&self.buf[self.ptr..])
    }

    pub fn value_string(&self) -> String {
        match self.type_() {
            Type::String => self.get::<String>(),
            Type::Int64 => self.get::<i64>().to_string(),
        }
    }
}

impl<'a> Iterator for RecordIterator<'a> {
    type Item = RecordField<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.schema.len() {
            return None;
        }
        let type_ = self.schema[self.idx];
        let start = self.ptr;
        self.ptr += match type_ {
            Type::String => {
                4 + i32::from_ne_bytes(self.buf[start..start + 4].try_into().unwrap()) as usize
            }
            Type::Int64 => 8,
        };
        self.idx += 1;
        Some(RecordField {
            type_,
            ptr: &self.buf[start..self.ptr],
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Table {
    schema: Schema,
    data: Vec<u8>,
    limits: Vec<i64>,
}

impl Table {
    pub fn new(schema: Schema) -> Self {
        Self {
            schema,
            data: Vec::new(),
            limits: vec![0],
        }
    }

    pub fn from_parts(schema: Schema, data: Vec<u8>, limits: Vec<i64>) -> Self {
        Self {
            schema,
            data,
            limits,
        }
    }

    pub fn from_lines<'a, I>(schema: Schema, lines: I) -> Result<Self, String>
    where
        I: IntoIterator<Item = &'a str>,
    {
        let mut table = Self::new(schema);
        for line in lines {
            table.push_back_line(line, None)?;
        }
        Ok(table)
    }

    pub fn schema(&self) -> &[Type] {
        &self.schema
    }

    pub fn size(&self) -> i64 {
        self.limits.len() as i64 - 1
    }

    pub fn empty(&self) -> bool {
        self.size() == 0
    }

    pub fn record(&self, i: i64) -> Record<'_> {
        let begin = self.limits[i as usize] as usize;
        let end = self.limits[i as usize + 1] as usize;
        Record {
            schema: &self.schema,
            buf: &self.data[begin..end],
        }
    }

    pub fn front(&self) -> Record<'_> {
        self.record(0)
    }

    pub fn push_back_record(&mut self, record: &Record<'_>) -> Result<(), SchemaMismatch> {
        if self.schema != record.schema {
            return Err(SchemaMismatch);
        }
        self.limits
            .push(self.limits.last().copied().unwrap() + record.raw_size());
        self.data.extend_from_slice(record.buf);
        Ok(())
    }

    pub fn push_back_line(
        &mut self,
        line: &str,
        record_id: Option<RecordId>,
    ) -> Result<(), String> {
        let mut i = 0usize;
        self.limits.push(*self.limits.last().unwrap());
        if let Some(record_id) = record_id {
            i += 1;
            self.push_i64(record_id);
        }
        for token in line.split('\t') {
            if i >= self.schema.len() {
                break;
            }
            match self.schema[i] {
                Type::String => self.push_str(token),
                Type::Int64 => self.push_i64(convert_string_i64(token)?),
            }
            i += 1;
        }
        if i < self.schema.len() {
            return Err("Missing fields in input line".to_string());
        }
        Ok(())
    }

    pub fn append(&mut self, table: &Table) -> Result<(), SchemaMismatch> {
        if self.schema != table.schema {
            return Err(SchemaMismatch);
        }
        let d = *self.limits.last().unwrap();
        self.data.extend_from_slice(&table.data);
        self.limits.reserve(table.limits.len().saturating_sub(1));
        for limit in table.limits.iter().skip(1) {
            self.limits.push(*limit + d);
        }
        Ok(())
    }

    pub fn append_lines<'a, I>(&mut self, lines: I) -> Result<(), String>
    where
        I: IntoIterator<Item = &'a str>,
    {
        for line in lines {
            self.push_back_line(line, None)?;
        }
        Ok(())
    }

    pub fn write(&self, buf: &mut TextBuffer) {
        for i in 0..self.size() {
            self.record(i).write(buf);
        }
    }

    pub fn write_record(&mut self, values: &[Value]) -> Result<(), String> {
        self.limits.push(*self.limits.last().unwrap());
        for i in 0..self.schema.len().min(values.len()) {
            match (self.schema[i], &values[i]) {
                (Type::String, Value::String(s)) => self.push_str(s),
                (Type::Int64, Value::Int64(x)) => self.push_i64(*x),
                _ => return Err("Invalid type in schema".to_string()),
            }
        }
        if values.len() > self.schema.len() {
            return Err("write_record with too many fields.".to_string());
        }
        if values.len() < self.schema.len() {
            return Err("Mismatching field count for Table::write_record".to_string());
        }
        Ok(())
    }

    pub fn sort(&mut self, col: usize, _threads: usize) -> Result<(), String> {
        *self = self.sorted(col, _threads)?;
        Ok(())
    }

    pub fn sorted(&self, col: usize, _threads: usize) -> Result<Self, String> {
        if col >= self.schema.len() || self.schema[col] != Type::Int64 {
            return Err("Invalid sort".to_string());
        }
        let mut v = Vec::with_capacity(self.size() as usize);
        for i in 0..self.size() {
            v.push((self.record(i).get_t::<i64>(col), i));
        }
        v.sort();
        Ok(self.shuffle(v.iter().map(|&(_, i)| i)))
    }

    pub fn map<F>(&self, _threads: usize, mut f: F, out: &mut File) -> Result<(), SchemaMismatch>
    where
        F: FnMut(&Record<'_>) -> Table,
    {
        for i in 0..self.size() {
            let mapped = f(&self.record(i));
            out.write_table(&mapped)?;
        }
        Ok(())
    }

    pub fn shuffle<I>(&self, iter: I) -> Self
    where
        I: IntoIterator<Item = i64>,
    {
        let mut t = Self::new(self.schema.clone());
        t.data.reserve(self.data.len());
        t.limits.reserve(self.limits.len());
        for i in iter {
            t.push_back_record(&self.record(i)).unwrap();
        }
        t
    }

    pub fn alloc_size(&self) -> i64 {
        (self.data.len() + self.limits.len() * std::mem::size_of::<i64>()) as i64
    }

    fn push_str(&mut self, s: &str) {
        self.push_i32(s.len() as i32);
        self.data.extend_from_slice(s.as_bytes());
        *self.limits.last_mut().unwrap() += s.len() as i64;
    }

    fn push_i32(&mut self, x: i32) {
        self.data.extend_from_slice(&x.to_ne_bytes());
        *self.limits.last_mut().unwrap() += 4;
    }

    fn push_i64(&mut self, x: i64) {
        self.data.extend_from_slice(&x.to_ne_bytes());
        *self.limits.last_mut().unwrap() += 8;
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BuildHelper {
    pub buffer: Vec<u8>,
    pub limits: Vec<i64>,
    pub counts: Vec<i32>,
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

    pub fn eof(&self) -> bool {
        self.cursor >= self.table.size()
    }

    pub fn size(&self) -> i64 {
        self.table.alloc_size()
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
        let mut size = 0i64;
        while !self.eof() && size < max_size {
            let record = self.table.record(self.cursor);
            size += record.raw_size();
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
        if self.schema.len() != T::TYPES.len() {
            return Err("Template parameters do not match schema.".to_string());
        }
        for (a, b) in self.schema.iter().zip(T::TYPES.iter()) {
            if a != b {
                return Err("Template parameters do not match schema.".to_string());
            }
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

    pub fn file_name(&self) -> String {
        self.file_name.clone()
    }

    pub fn table(&self) -> &Table {
        &self.table
    }
}

pub fn count_lines_str(s: &str) -> i64 {
    s.split_terminator('\n').count() as i64
}

pub fn count_lines(file_name: &str) -> Result<i64, String> {
    let text = std::fs::read_to_string(file_name).map_err(|e| e.to_string())?;
    Ok(count_lines_str(&text))
}

pub fn read_text_mt<F>(data: &[u8], max_size: i64, _threads: usize, mut callback: F)
where
    F: FnMut(i64, &[u8]),
{
    let mut next_chunk = 0i64;
    let mut total = 0i64;
    let mut pos = 0usize;
    loop {
        if pos >= data.len() {
            break;
        }
        let read_end = (pos + READ_TEXT_MT_SIZE).min(data.len());
        let mut end = read_end;
        if read_end - pos == READ_TEXT_MT_SIZE {
            while end < data.len() && data[end] != b'\n' {
                end += 1;
            }
        }
        let n = end - pos;
        total += n as i64;
        if n > 0 {
            callback(next_chunk, &data[pos..end]);
            next_chunk += 1;
        }
        if n < READ_TEXT_MT_SIZE || total + READ_TEXT_MT_SIZE as i64 > max_size {
            break;
        }
        pos = if end < data.len() && data[end] == b'\n' {
            end + 1
        } else {
            end
        };
    }
}

pub fn merge_files(files: &mut [&mut File], column: usize) -> Result<File, String> {
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
    while let Some(Reverse((_, f))) = queue.pop() {
        out.write(&tables[f].front()).map_err(|e| e.to_string())?;
        tables[f] = files[f].read_record();
        if !tables[f].empty() {
            queue.push(Reverse((tables[f].front().get_t::<i64>(column), f)));
        }
    }
    Ok(out)
}

pub fn join_schema(schemas: &[Schema; 2], output_fields: &[FileColumn]) -> Schema {
    let mut schema = Vec::with_capacity(output_fields.len());
    for field in output_fields {
        schema.push(schemas[field.file][field.column]);
    }
    schema
}

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
    join(file1, file2, column1, column2, output_fields, &mut out)?;
    Ok(out)
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
