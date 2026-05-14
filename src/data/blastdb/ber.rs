use crate::util::io::{File, IoResult};

pub fn read_be32(file: &mut File) -> IoResult<u32> {
    let mut buf = [0u8; 4];
    file.read_exact(&mut buf)?;
    let mut value = 0u32;
    for b in buf {
        value = (value << 8) | u32::from(b);
    }
    Ok(value)
}

pub fn read_le64(file: &mut File) -> IoResult<u64> {
    let mut buf = [0u8; 8];
    file.read_exact(&mut buf)?;
    let mut value = 0u64;
    for b in buf.into_iter().rev() {
        value = (value << 8) | u64::from(b);
    }
    Ok(value)
}

pub fn decode_integer(data: &[u8]) -> i64 {
    if data.len() > std::mem::size_of::<i64>() {
        return 0;
    }
    if data.is_empty() {
        return 0;
    }
    let mut value = if (data[0] & 0x80) != 0 { -1i64 } else { 0 };
    for &byte in data {
        value = (value << 8) | i64::from(byte);
    }
    value
}

pub fn read_pascal_string(file: &mut File) -> IoResult<String> {
    let l = read_be32(file)? as usize;
    let mut s = vec![0u8; l];
    file.read_exact(&mut s)?;
    Ok(String::from_utf8_lossy(&s).into_owned())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_decode_integer_sign_extends_like_cpp() {
        assert_eq!(decode_integer(&[0x7f]), 127);
        assert_eq!(decode_integer(&[0x00, 0x80]), 128);
        assert_eq!(decode_integer(&[0xff]), -1);
        assert_eq!(decode_integer(&[0xff, 0x7f]), -129);
        assert_eq!(decode_integer(&[0xff; 9]), 0);
        assert_eq!(decode_integer(&[]), 0);
    }

    #[test]
    fn test_read_ber_endian_and_pascal_values() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-ber-{}-{}.bin",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            f.write_all(&[0x12, 0x34, 0x56, 0x78]).unwrap();
            f.write_all(&[0x08, 0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01])
                .unwrap();
            f.write_all(&[0x00, 0x00, 0x00, 0x04, b't', b'e', b's', b't'])
                .unwrap();
        }
        let mut file = File::open(path.to_str().unwrap(), "r").unwrap();
        assert_eq!(read_be32(&mut file).unwrap(), 0x1234_5678);
        assert_eq!(read_le64(&mut file).unwrap(), 0x0102_0304_0506_0708);
        assert_eq!(read_pascal_string(&mut file).unwrap(), "test");
        std::fs::remove_file(path).unwrap();
    }
}
