use crate::util::algo::read_varuint32;

/// Compact vector-of-vectors storage matching legacy `dmnd/compact_array.h`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompactArray {
    data: Vec<u8>,
    limits: CompactArrayLimits,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum CompactArrayLimits {
    U32(Vec<u32>),
    I64(Vec<i64>),
}

impl CompactArray {
    pub fn new(data: Vec<u8>, size: usize, data_size: usize) -> Result<Self, String> {
        let mut data = data;
        data.resize(data_size, 0);
        let limits = if data_size > u32::MAX as usize {
            CompactArrayLimits::I64(Self::init_i64(size as i64, data_size as i64, &data)?)
        } else {
            CompactArrayLimits::U32(Self::init_u32(size as i64, data_size as i64, &data)?)
        };
        Ok(Self { data, limits })
    }

    pub fn get(&self, i: usize) -> Result<Vec<i32>, String> {
        match &self.limits {
            CompactArrayLimits::U32(limits) => Self::read_vec(&self.data[limits[i] as usize..]),
            CompactArrayLimits::I64(limits) => Self::read_vec(&self.data[limits[i] as usize..]),
        }
    }

    pub fn size(&self) -> usize {
        match &self.limits {
            CompactArrayLimits::U32(limits) => limits.len() - 1,
            CompactArrayLimits::I64(limits) => limits.len() - 1,
        }
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    fn init_u32(size: i64, data_size: i64, data: &[u8]) -> Result<Vec<u32>, String> {
        let mut limits = Vec::with_capacity(size as usize + 1);
        limits.push(0);
        let mut ptr = 0usize;
        for _ in 0..size {
            ptr = Self::skip_vec(data, ptr)?;
            if ptr > u32::MAX as usize {
                return Err("Array size overflow.".to_string());
            }
            limits.push(ptr as u32);
        }
        if *limits.last().unwrap() as i64 != data_size {
            return Err("Error loading CompactArray.".to_string());
        }
        Ok(limits)
    }

    fn init_i64(size: i64, data_size: i64, data: &[u8]) -> Result<Vec<i64>, String> {
        let mut limits = Vec::with_capacity(size as usize + 1);
        limits.push(0);
        let mut ptr = 0usize;
        for _ in 0..size {
            ptr = Self::skip_vec(data, ptr)?;
            limits.push(ptr as i64);
        }
        if *limits.last().unwrap() != data_size {
            return Err("Error loading CompactArray.".to_string());
        }
        Ok(limits)
    }

    fn read_vec(ptr: &[u8]) -> Result<Vec<i32>, String> {
        let (n, mut offset) = read_varuint32(ptr)?;
        let mut out = Vec::with_capacity(n as usize);
        for _ in 0..n {
            let (v, consumed) = read_varuint32(&ptr[offset..])?;
            offset += consumed;
            out.push(v as i32);
        }
        Ok(out)
    }

    fn skip_vec(data: &[u8], ptr: usize) -> Result<usize, String> {
        let (n, consumed) = read_varuint32(&data[ptr..])?;
        let mut ptr = ptr + consumed;
        for _ in 0..n {
            let (_, consumed) = read_varuint32(&data[ptr..])?;
            ptr += consumed;
        }
        Ok(ptr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::algo::write_varuint32;

    fn push_vec(out: &mut Vec<u8>, values: &[u32]) {
        write_varuint32(values.len() as u32, out);
        for &value in values {
            write_varuint32(value, out);
        }
    }

    #[test]
    fn test_compact_array_roundtrip() {
        let mut data = Vec::new();
        push_vec(&mut data, &[1, 2, 130]);
        push_vec(&mut data, &[]);
        push_vec(&mut data, &[16_384, 1 << 28]);

        let array = CompactArray::new(data.clone(), 3, data.len()).unwrap();
        assert_eq!(array.size(), 3);
        assert_eq!(array.data(), data);
        assert_eq!(array.get(0).unwrap(), vec![1, 2, 130]);
        assert_eq!(array.get(1).unwrap(), Vec::<i32>::new());
        assert_eq!(array.get(2).unwrap(), vec![16_384, 1 << 28]);
    }

    #[test]
    fn test_compact_array_rejects_size_mismatch() {
        let mut data = Vec::new();
        push_vec(&mut data, &[1, 2, 3]);
        assert_eq!(
            CompactArray::new(data.clone(), 2, data.len())
                .unwrap_err()
                .as_str(),
            "Format error: Invalid varint encoding."
        );
        assert_eq!(
            CompactArray::new(data, 0, 1).unwrap_err().as_str(),
            "Error loading CompactArray."
        );
    }
}
