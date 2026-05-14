#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Finish;

/// Size-prefixed legacy record reader matching `dmnd/record_reader.h`.
#[derive(Debug, Clone)]
pub struct DynamicRecordReader<'a> {
    data: &'a [u8],
    pos: usize,
    size: u64,
}

impl<'a> DynamicRecordReader<'a> {
    pub fn new(data: &'a [u8]) -> Result<Self, String> {
        if data.len() < std::mem::size_of::<u64>() {
            return Err("DynamicRecordReader: missing record size".to_string());
        }
        let mut size = [0u8; 8];
        size.copy_from_slice(&data[..8]);
        Ok(Self {
            data,
            pos: 8,
            size: u64::from_le_bytes(size),
        })
    }

    pub fn read_u64(&mut self) -> Result<u64, String> {
        if self.size >= std::mem::size_of::<u64>() as u64 {
            let bytes = self.read_exact(std::mem::size_of::<u64>())?;
            let mut out = [0u8; 8];
            out.copy_from_slice(bytes);
            self.size -= std::mem::size_of::<u64>() as u64;
            Ok(u64::from_le_bytes(out))
        } else {
            Ok(0)
        }
    }

    pub fn read_ulong(&mut self) -> Result<u64, String> {
        self.read_u64()
    }

    pub fn read_i32(&mut self) -> Result<i32, String> {
        if self.size >= std::mem::size_of::<i32>() as u64 {
            let bytes = self.read_exact(std::mem::size_of::<i32>())?;
            let mut out = [0u8; 4];
            out.copy_from_slice(bytes);
            self.size -= std::mem::size_of::<i32>() as u64;
            Ok(i32::from_le_bytes(out))
        } else {
            Ok(0)
        }
    }

    pub fn read_bytes(&mut self, out: &mut [u8]) -> Result<&mut Self, String> {
        let s = out.len();
        if self.size >= s as u64 {
            out.copy_from_slice(self.read_exact(s)?);
            self.size -= s as u64;
        } else {
            out.fill(0);
        }
        Ok(self)
    }

    pub fn finish(&mut self, _finish: Finish) -> Result<(), String> {
        if self.size == 0 {
            return Ok(());
        }
        self.read_exact(self.size as usize)?;
        self.size = 0;
        Ok(())
    }

    pub fn remaining_size(&self) -> u64 {
        self.size
    }

    pub fn position(&self) -> usize {
        self.pos
    }

    fn read_exact(&mut self, n: usize) -> Result<&'a [u8], String> {
        let end = self.pos + n;
        if end > self.data.len() {
            return Err("DynamicRecordReader: unexpected end of input".to_string());
        }
        let out = &self.data[self.pos..end];
        self.pos = end;
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(payload: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(&(payload.len() as u64).to_le_bytes());
        out.extend_from_slice(payload);
        out
    }

    #[test]
    fn test_dynamic_record_reader_reads_present_fields() {
        let mut payload = Vec::new();
        payload.extend_from_slice(&7u64.to_le_bytes());
        payload.extend_from_slice(&(-3i32).to_le_bytes());
        payload.extend_from_slice(&[1, 2, 3]);
        payload.extend_from_slice(&[9, 9]);

        let data = record(&payload);
        let mut r = DynamicRecordReader::new(&data).unwrap();
        assert_eq!(r.read_u64().unwrap(), 7);
        assert_eq!(r.read_i32().unwrap(), -3);
        let mut bytes = [0u8; 3];
        r.read_bytes(&mut bytes).unwrap();
        assert_eq!(bytes, [1, 2, 3]);
        assert_eq!(r.remaining_size(), 2);
        r.finish(Finish).unwrap();
        assert_eq!(r.remaining_size(), 0);
        assert_eq!(r.position(), data.len());
    }

    #[test]
    fn test_dynamic_record_reader_blanks_absent_fields() {
        let data = record(&[1, 2, 3]);
        let mut r = DynamicRecordReader::new(&data).unwrap();
        assert_eq!(r.read_u64().unwrap(), 0);
        assert_eq!(r.remaining_size(), 3);
        assert_eq!(r.read_i32().unwrap(), 0);
        let mut bytes = [7u8; 4];
        r.read_bytes(&mut bytes).unwrap();
        assert_eq!(bytes, [0, 0, 0, 0]);
        assert_eq!(r.remaining_size(), 3);
        r.finish(Finish).unwrap();
        assert_eq!(r.remaining_size(), 0);
    }
}
