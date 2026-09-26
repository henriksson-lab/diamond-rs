//! Binary-backed TSV table storage and transforms.
//!
//! This mirrors `diamond/src/util/tsv/table.{h,cpp}`. Records are concatenated
//! in `data`, while `limits` stores native signed 64-bit byte offsets.

use crate::util::string::convert_string_i64;
use crate::util::text_buffer::TextBuffer;

use super::{File, Record, RecordId, Schema, SchemaMismatch, Type, Value};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Table {
    pub(super) schema: Schema,
    pub(super) data: Vec<u8>,
    pub(super) limits: Vec<i64>,
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
        table.append_lines(lines)?;
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
        // C++ relies on this same-schema precondition without checking it.
        // Rust reports the violation before the record bytes can corrupt the
        // destination table's interpretation.
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
        let mut field = 0usize;
        self.limits.push(*self.limits.last().unwrap());
        if let Some(record_id) = record_id {
            field += 1;
            self.push_i64(record_id);
        }
        for token in line.split('\t') {
            if field >= self.schema.len() {
                break;
            }
            match self.schema[field] {
                Type::String => self.push_str(token),
                Type::Int64 => self.push_i64(convert_string_i64(token)?),
            }
            field += 1;
        }
        if field < self.schema.len() {
            return Err("Missing fields in input line".to_string());
        }
        Ok(())
    }

    pub fn append(&mut self, table: &Table) -> Result<(), SchemaMismatch> {
        if self.schema != table.schema {
            return Err(SchemaMismatch);
        }
        let offset = *self.limits.last().unwrap();
        self.data.extend_from_slice(&table.data);
        self.limits.reserve(table.limits.len().saturating_sub(1));
        for limit in table.limits.iter().skip(1) {
            self.limits.push(*limit + offset);
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

    pub fn sort(&mut self, col: usize, threads: usize) -> Result<(), String> {
        *self = self.sorted(col, threads)?;
        Ok(())
    }

    pub fn sorted(&self, col: usize, _threads: usize) -> Result<Self, String> {
        if col >= self.schema.len() || self.schema[col] != Type::Int64 {
            return Err("Invalid sort".to_string());
        }
        let mut order = Vec::with_capacity(self.size() as usize);
        for i in 0..self.size() {
            order.push((self.record(i).get_t::<i64>(col), i));
        }
        order.sort();
        Ok(self.shuffle(order.iter().map(|&(_, i)| i)))
    }

    pub fn map<F>(&self, _threads: usize, mut f: F, out: &mut File) -> Result<(), SchemaMismatch>
    where
        F: FnMut(&Record<'_>) -> Table,
    {
        // The upstream worker pool writes through an ordered reorder queue.
        // Sequential evaluation preserves the same observable record order.
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
        let mut table = Self::new(self.schema.clone());
        table.data.reserve(self.data.len());
        table.limits.reserve(self.limits.len());
        for i in iter {
            table.push_back_record(&self.record(i)).unwrap();
        }
        table
    }

    pub fn alloc_size(&self) -> i64 {
        (self.data.len() + self.limits.len() * std::mem::size_of::<i64>()) as i64
    }

    pub(super) fn push_str(&mut self, s: &str) {
        self.push_i32(s.len() as i32);
        self.data.extend_from_slice(s.as_bytes());
        *self.limits.last_mut().unwrap() += s.len() as i64;
    }

    pub(super) fn push_i32(&mut self, x: i32) {
        self.data.extend_from_slice(&x.to_ne_bytes());
        *self.limits.last_mut().unwrap() += 4;
    }

    pub(super) fn push_i64(&mut self, x: i64) {
        self.data.extend_from_slice(&x.to_ne_bytes());
        *self.limits.last_mut().unwrap() += 8;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn append_offsets_and_allocation_size_match_binary_layout() {
        let schema = vec![Type::Int64, Type::String];
        let mut left = Table::from_lines(schema.clone(), ["2\tb"]).unwrap();
        let right = Table::from_lines(schema, ["1\tx", "1\ta"]).unwrap();
        left.append(&right).unwrap();

        assert_eq!(left.size(), 3);
        assert_eq!(left.limits, [0, 13, 26, 39]);
        assert_eq!(left.alloc_size(), 39 + 4 * 8);
        assert_eq!(left.record(2).get_t::<String>(1), "a");
    }

    #[test]
    fn integer_sort_ties_use_original_record_index() {
        let schema = vec![Type::Int64, Type::String];
        let table = Table::from_lines(schema, ["2\tb", "1\tx", "1\ta"]).unwrap();

        let sorted = table.sorted(0, 4).unwrap();

        assert_eq!(sorted.record(0).get_t::<String>(1), "x");
        assert_eq!(sorted.record(1).get_t::<String>(1), "a");
        assert_eq!(sorted.record(2).get_t::<String>(1), "b");
        assert_eq!(table.sorted(1, 1).unwrap_err(), "Invalid sort");
    }

    #[test]
    fn record_id_extra_fields_and_shuffle_match_upstream() {
        let schema = vec![Type::Int64, Type::String];
        let mut table = Table::new(schema);
        table.push_back_line("alpha\tignored", Some(7)).unwrap();
        table.push_back_line("beta", Some(8)).unwrap();

        assert_eq!(table.record(0).get_t::<i64>(0), 7);
        assert_eq!(table.record(0).get_t::<String>(1), "alpha");
        let shuffled = table.shuffle([1, 0, 1]);
        assert_eq!(shuffled.size(), 3);
        assert_eq!(shuffled.record(0).get_t::<i64>(0), 8);
        assert_eq!(shuffled.record(1).get_t::<i64>(0), 7);
        assert_eq!(shuffled.record(2).get_t::<i64>(0), 8);
    }

    #[test]
    fn write_record_count_errors_leave_cpp_compatible_partial_record() {
        let mut too_many = Table::new(vec![Type::Int64]);
        assert_eq!(
            too_many
                .write_record(&[Value::Int64(4), Value::Int64(5)])
                .unwrap_err(),
            "write_record with too many fields."
        );
        assert_eq!(too_many.size(), 1);
        assert_eq!(too_many.front().get_t::<i64>(0), 4);

        let mut too_few = Table::new(vec![Type::Int64, Type::String]);
        assert_eq!(
            too_few.write_record(&[Value::Int64(4)]).unwrap_err(),
            "Mismatching field count for Table::write_record"
        );
        assert_eq!(too_few.size(), 1);
        assert_eq!(too_few.front().raw_size(), 8);
    }
}
