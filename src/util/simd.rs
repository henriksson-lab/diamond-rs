#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Arch {
    None,
    Generic,
    Sse4_1,
    Avx2,
    Avx512,
    Neon,
}

pub const SSSE3: i32 = 1;
pub const POPCNT: i32 = 2;
pub const SSE4_1: i32 = 4;
pub const AVX2: i32 = 8;
pub const AVX512: i32 = 16;
pub const NEON: i32 = 32;

pub fn cpuid(info: &mut [i32; 4], info_type: i32) {
    #[cfg(target_arch = "x86")]
    {
        #[cfg(target_env = "msvc")]
        let r = std::arch::x86::__cpuid_count(info_type as u32, 0);
        #[cfg(not(target_env = "msvc"))]
        let r = unsafe { std::arch::x86::__cpuid_count(info_type as u32, 0) };
        info[0] = r.eax as i32;
        info[1] = r.ebx as i32;
        info[2] = r.ecx as i32;
        info[3] = r.edx as i32;
    }
    #[cfg(target_arch = "x86_64")]
    {
        #[cfg(target_env = "msvc")]
        let r = std::arch::x86_64::__cpuid_count(info_type as u32, 0);
        #[cfg(not(target_env = "msvc"))]
        let r = unsafe { std::arch::x86_64::__cpuid_count(info_type as u32, 0) };
        info[0] = r.eax as i32;
        info[1] = r.ebx as i32;
        info[2] = r.ecx as i32;
        info[3] = r.edx as i32;
    }
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {
        let _ = info_type;
        *info = [0; 4];
    }
}

pub fn flags() -> i32 {
    let mut flags = 0;
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if std::is_x86_feature_detected!("ssse3") {
            flags |= SSSE3;
        }
        if std::is_x86_feature_detected!("popcnt") {
            flags |= POPCNT;
        }
        if std::is_x86_feature_detected!("sse4.1") {
            flags |= SSE4_1;
        }
        if std::is_x86_feature_detected!("avx2") {
            flags |= AVX2;
        }
        if std::is_x86_feature_detected!("avx512f") && std::is_x86_feature_detected!("avx512bw") {
            flags |= AVX512;
        }
    }
    #[cfg(any(target_arch = "aarch64", target_arch = "arm"))]
    {
        flags |= NEON;
    }
    flags
}

pub fn init_arch() -> Arch {
    let flags = flags();
    if flags & NEON != 0 {
        return Arch::Neon;
    }
    if flags & AVX512 != 0 {
        return Arch::Avx512;
    }
    if (flags & (SSSE3 | POPCNT | SSE4_1 | AVX2)) == (SSSE3 | POPCNT | SSE4_1 | AVX2) {
        return Arch::Avx2;
    }
    if (flags & (SSSE3 | POPCNT | SSE4_1)) == (SSSE3 | POPCNT | SSE4_1) {
        return Arch::Sse4_1;
    }
    Arch::Generic
}

pub fn arch() -> Arch {
    init_arch()
}

pub fn features() -> String {
    let flags = flags();
    let mut r = Vec::new();
    if flags & NEON != 0 {
        r.push("neon");
    }
    if flags & SSSE3 != 0 {
        r.push("ssse3");
    }
    if flags & POPCNT != 0 {
        r.push("popcnt");
    }
    if flags & SSE4_1 != 0 {
        r.push("sse4.1");
    }
    if flags & AVX2 != 0 {
        r.push("avx2");
    }
    if r.is_empty() {
        "None".to_string()
    } else {
        r.join(" ")
    }
}

pub fn transpose(data: &[&[i8]], n: usize, out: &mut [i8], width: usize) {
    assert!(width == 16 || width == 32);
    assert!(n <= width);
    assert!(out.len() >= width * width);
    out[..width * width].fill(0);
    let row_offset = width - n;
    for row in 0..n {
        assert!(data[row].len() >= width);
        for col in 0..width {
            out[col * width + row_offset + row] = data[row][col];
        }
    }
}

pub fn transpose_16(data: &[&[i8]], n: usize, out: &mut [i8]) {
    transpose(data, n, out, 16);
}

pub fn transpose_32(data: &[&[i8]], n: usize, out: &mut [i8]) {
    transpose(data, n, out, 32);
}

pub fn transpose_offset(data: &[&[i8]], n: usize, offset: isize, out: &mut [i8], width: usize) {
    assert!(offset >= 0);
    assert!(width == 16 || width == 32);
    assert!(n <= width);
    assert!(out.len() >= width * width);
    out[..width * width].fill(0);
    let start = offset as usize * width;
    let row_offset = width - n;
    for row in 0..n {
        assert!(data[row].len() >= start + width);
        for col in 0..width {
            out[col * width + row_offset + row] = data[row][start + col];
        }
    }
}

pub fn transpose_offset_i16(
    data: &[&[i16]],
    n: usize,
    offset: isize,
    out: &mut [i16],
    width: usize,
) {
    assert!(offset >= 0);
    assert!(width == 16);
    assert!(n <= width);
    assert!(out.len() >= width * width);
    out[..width * width].fill(0);
    let start = offset as usize * width;
    for row in 0..n {
        assert!(data[row].len() >= start + width);
        for col in 0..width {
            out[col * width + row] = data[row][start + col];
        }
    }
}

pub fn transpose_offset_8bit_i16(data: &[&[i16]], n: usize, offset: isize, out: &mut [i16]) {
    assert!(offset >= 0);
    assert!(n <= 16);
    assert!(out.len() >= 16 * 16);
    out[..16 * 16].fill(0);
    let start = offset as usize * 16;
    for row in 0..n {
        assert!(data[row].len() >= start + 16);
        for col in 0..16 {
            out[col * 16 + row] = data[row][start + col] & 0xff;
        }
    }
}

pub mod dispatch_arch {
    pub mod simd {
        use std::marker::PhantomData;

        #[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
        pub struct Vector<T> {
            _marker: PhantomData<T>,
        }

        #[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
        pub struct Traits<T> {
            _marker: PhantomData<T>,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(width: usize, n: usize) -> Vec<Vec<i8>> {
        (0..n)
            .map(|r| (0..width).map(|c| (r * 40 + c) as i8).collect())
            .collect()
    }

    #[test]
    fn test_features_and_arch_are_callable() {
        let mut info = [0; 4];
        cpuid(&mut info, 0);
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        assert!(info[0] >= 0);
        let _ = flags();
        let _ = arch();
        assert!(!features().is_empty());
    }

    #[test]
    fn test_transpose_16_full_and_partial() {
        let full_rows = rows(16, 16);
        let refs: Vec<_> = full_rows.iter().map(Vec::as_slice).collect();
        let mut out = vec![0i8; 16 * 16];
        transpose_16(&refs, 16, &mut out);
        assert_eq!(out[0], full_rows[0][0]);
        assert_eq!(out[1], full_rows[1][0]);
        assert_eq!(out[16], full_rows[0][1]);
        assert_eq!(out[16 * 15 + 15], full_rows[15][15]);

        let partial_rows = rows(16, 3);
        let refs: Vec<_> = partial_rows.iter().map(Vec::as_slice).collect();
        transpose_16(&refs, 3, &mut out);
        assert_eq!(out[0], 0);
        assert_eq!(out[13], partial_rows[0][0]);
        assert_eq!(out[14], partial_rows[1][0]);
        assert_eq!(out[15], partial_rows[2][0]);
        assert_eq!(out[16 + 13], partial_rows[0][1]);
    }

    #[test]
    fn test_transpose_32_and_offset() {
        let rows = rows(64, 32);
        let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
        let mut out = vec![0i8; 32 * 32];
        transpose_32(&refs, 32, &mut out);
        assert_eq!(out[31], rows[31][0]);
        assert_eq!(out[32], rows[0][1]);
        transpose_offset(&refs, 32, 1, &mut out, 32);
        assert_eq!(out[0], rows[0][32]);
        assert_eq!(out[32 * 31 + 31], rows[31][63]);
    }

    #[test]
    fn test_transpose_offset_i16() {
        let rows: Vec<Vec<i16>> = (0..16)
            .map(|r| (0..32).map(|c| (r * 100 + c) as i16).collect())
            .collect();
        let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
        let mut out = vec![0i16; 16 * 16];
        transpose_offset_i16(&refs, 16, 1, &mut out, 16);
        assert_eq!(out[0], rows[0][16]);
        assert_eq!(out[15], rows[15][16]);
        assert_eq!(out[16], rows[0][17]);
        transpose_offset_8bit_i16(&refs, 16, 1, &mut out);
        assert_eq!(out[0], rows[0][16] & 0xff);
    }
}
