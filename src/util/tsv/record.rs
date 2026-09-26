//! Binary-backed TSV record views.
//!
//! This mirrors `diamond/src/util/tsv/record.{h,cpp}`. Integer and string
//! length fields use native byte order, matching the upstream in-memory table
//! representation.

use crate::util::text_buffer::TextBuffer;

use super::{Table, Type};

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
    pub(super) schema: &'a [Type],
    pub(super) buf: &'a [u8],
}

impl<'a> Record<'a> {
    pub fn new(schema: &'a [Type], begin: &'a [u8], end: usize) -> Self {
        Self {
            schema,
            buf: &begin[..end],
        }
    }

    pub fn get_t<T: TsvValue>(&self, i: usize) -> T {
        let mut iterator = self.begin();
        for _ in 0..i {
            iterator.next();
        }
        iterator.get()
    }

    pub fn get(&self, i: usize) -> String {
        let mut iterator = self.begin();
        for _ in 0..i {
            iterator.next();
        }
        iterator.value_string()
    }

    pub fn begin(&self) -> RecordIterator<'a> {
        RecordIterator::new(self.schema, self.buf)
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
        for (idx, field) in self.begin().enumerate() {
            if idx > 0 {
                buf.append_char('\t');
            }
            match field.type_() {
                Type::Int64 => buf.append_display(field.get::<i64>()),
                Type::String => buf.append_str(&field.get::<String>()),
            };
        }
        buf.append_char('\n');
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RecordField<'a> {
    type_: Type,
    ptr: &'a [u8],
}

impl RecordField<'_> {
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

impl RecordIterator<'_> {
    /// Construct the iterator at the first schema field.
    pub fn new<'a>(schema: &'a [Type], buf: &'a [u8]) -> RecordIterator<'a> {
        RecordIterator {
            schema,
            buf,
            idx: 0,
            ptr: 0,
        }
    }

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

impl PartialEq for RecordIterator<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.schema.as_ptr() == other.schema.as_ptr()
            && self.schema.len() == other.schema.len()
            && self.idx == other.idx
    }
}

impl Eq for RecordIterator<'_> {}

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_binary_layout_iteration_and_text_are_exact() {
        let schema = [Type::String, Type::Int64, Type::String];
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&3i32.to_ne_bytes());
        bytes.extend_from_slice(b"a\0b");
        bytes.extend_from_slice(&(-42i64).to_ne_bytes());
        bytes.extend_from_slice(&2i32.to_ne_bytes());
        bytes.extend_from_slice("é".as_bytes());
        let record = Record::new(&schema, &bytes, bytes.len());

        assert_eq!(record.raw_size(), bytes.len() as i64);
        assert_eq!(record.get_t::<String>(0).as_bytes(), b"a\0b");
        assert_eq!(record.get_t::<i64>(1), -42);
        assert_eq!(record.get(1), "-42");
        assert_eq!(record.get_t::<String>(2), "é");

        let fields: Vec<_> = record.begin().collect();
        assert_eq!(fields.len(), 3);
        assert_eq!(fields[0].type_(), Type::String);
        assert_eq!(fields[1].type_(), Type::Int64);
        assert_eq!(fields[2].value_string(), "é");

        let mut iterator = record.begin();
        assert_ne!(iterator, record.end());
        iterator.by_ref().for_each(drop);
        assert_eq!(iterator, record.end());

        let mut text = TextBuffer::new();
        record.write(&mut text);
        assert_eq!(text.data(), "a\0b\t-42\té\n".as_bytes());
    }

    #[test]
    fn empty_schema_writes_one_record_terminator() {
        let record = Record::new(&[], &[], 0);
        assert_eq!(record.raw_size(), 0);
        assert!(record.begin().next().is_none());

        let mut text = TextBuffer::new();
        record.write(&mut text);
        assert_eq!(text.data(), b"\n");
    }
}
