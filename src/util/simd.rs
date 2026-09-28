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
        let r = std::arch::x86::__cpuid_count(info_type as u32, 0);
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
        let r = std::arch::x86_64::__cpuid_count(info_type as u32, 0);
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
    static ARCH: std::sync::OnceLock<Arch> = std::sync::OnceLock::new();
    *ARCH.get_or_init(init_arch)
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
    assert!(data.len() >= n);
    assert!(out.len() >= width * width);
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    unsafe {
        if width == 32 && std::arch::is_x86_feature_detected!("avx2") {
            let mut pointers = [std::ptr::null(); 32];
            for (pointer, row) in pointers.iter_mut().zip(data.iter().take(n)) {
                *pointer = row.as_ptr();
            }
            transpose_32_avx2(&pointers[..n], n, 0, out);
            return;
        }
        if width == 16 && std::arch::is_x86_feature_detected!("sse2") {
            transpose_16_sse2(data, n, out);
            return;
        }
    }
    #[cfg(target_arch = "aarch64")]
    unsafe {
        if width == 16 {
            transpose_16_neon(data, n, out);
            return;
        }
    }
    out[..width * width].fill(0);
    let row_offset = width - n;
    for row in 0..n {
        assert!(data[row].len() >= width);
        for col in 0..width {
            out[col * width + row_offset + row] = data[row][col];
        }
    }
}

#[inline]
#[cfg(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64"))]
fn bit_reverse(value: usize, bits: u32) -> usize {
    value.reverse_bits() >> (usize::BITS - bits)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "sse2")]
unsafe fn transpose_16_sse2(data: &[&[i8]], n: usize, out: &mut [i8]) {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    let zero = _mm_setzero_si128();
    let mut rows = [zero; 16];
    let row0 = 16 - n;
    for (dst, src) in rows[row0..].iter_mut().zip(data.iter().take(n)) {
        debug_assert!(src.len() >= 16);
        *dst = _mm_loadu_si128(src.as_ptr().cast());
    }
    for group in (0..16).step_by(2) {
        let a = rows[group];
        let b = rows[group + 1];
        rows[group] = _mm_unpacklo_epi8(a, b);
        rows[group + 1] = _mm_unpackhi_epi8(a, b);
    }
    for group in (0..16).step_by(4) {
        for offset in 0..2 {
            let a = rows[group + offset];
            let b = rows[group + offset + 2];
            rows[group + offset] = _mm_unpacklo_epi16(a, b);
            rows[group + offset + 2] = _mm_unpackhi_epi16(a, b);
        }
    }
    for group in (0..16).step_by(8) {
        for offset in 0..4 {
            let a = rows[group + offset];
            let b = rows[group + offset + 4];
            rows[group + offset] = _mm_unpacklo_epi32(a, b);
            rows[group + offset + 4] = _mm_unpackhi_epi32(a, b);
        }
    }
    for offset in 0..8 {
        let a = rows[offset];
        let b = rows[offset + 8];
        rows[offset] = _mm_unpacklo_epi64(a, b);
        rows[offset + 8] = _mm_unpackhi_epi64(a, b);
    }
    for column in 0..16 {
        _mm_storeu_si128(
            out.as_mut_ptr().add(column * 16).cast(),
            rows[bit_reverse(column, 4)],
        );
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn transpose_32_avx2(
    data: &[*const i8],
    n: usize,
    offset: usize,
    out: &mut [i8],
) {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    debug_assert!(n <= 32 && data.len() >= n && out.len() >= 32 * 32);
    let row0 = 32 - n;
    let mut rows = if n == 32 {
        macro_rules! load_row {
            ($i:expr) => {
                _mm256_loadu_si256(data.get_unchecked($i).add(offset).cast())
            };
        }
        [
            load_row!(0),
            load_row!(1),
            load_row!(2),
            load_row!(3),
            load_row!(4),
            load_row!(5),
            load_row!(6),
            load_row!(7),
            load_row!(8),
            load_row!(9),
            load_row!(10),
            load_row!(11),
            load_row!(12),
            load_row!(13),
            load_row!(14),
            load_row!(15),
            load_row!(16),
            load_row!(17),
            load_row!(18),
            load_row!(19),
            load_row!(20),
            load_row!(21),
            load_row!(22),
            load_row!(23),
            load_row!(24),
            load_row!(25),
            load_row!(26),
            load_row!(27),
            load_row!(28),
            load_row!(29),
            load_row!(30),
            load_row!(31),
        ]
    } else {
        let zero = _mm256_setzero_si256();
        let mut rows = [zero; 32];
        for (dst, &src) in rows[row0..].iter_mut().zip(data.iter().take(n)) {
            *dst = _mm256_loadu_si256(src.add(offset).cast());
        }
        rows
    };
    macro_rules! unpack8 {
        ($a:expr, $b:expr) => {{
            let x = rows[$a];
            let y = rows[$b];
            rows[$a] = _mm256_unpacklo_epi8(x, y);
            rows[$b] = _mm256_unpackhi_epi8(x, y);
        }};
    }
    macro_rules! unpack16 {
        ($a:expr, $b:expr) => {{
            let x = rows[$a];
            let y = rows[$b];
            rows[$a] = _mm256_unpacklo_epi16(x, y);
            rows[$b] = _mm256_unpackhi_epi16(x, y);
        }};
    }
    macro_rules! unpack32 {
        ($a:expr, $b:expr) => {{
            let x = rows[$a];
            let y = rows[$b];
            rows[$a] = _mm256_unpacklo_epi32(x, y);
            rows[$b] = _mm256_unpackhi_epi32(x, y);
        }};
    }
    macro_rules! unpack64 {
        ($a:expr, $b:expr) => {{
            let x = rows[$a];
            let y = rows[$b];
            rows[$a] = _mm256_unpacklo_epi64(x, y);
            rows[$b] = _mm256_unpackhi_epi64(x, y);
        }};
    }
    macro_rules! unpack128 {
        ($a:expr, $b:expr) => {{
            let x = rows[$a];
            let y = rows[$b];
            rows[$a] = _mm256_permute2x128_si256(x, y, 0x20);
            rows[$b] = _mm256_permute2x128_si256(x, y, 0x31);
        }};
    }
    macro_rules! pairs {
        ($op:ident; $(($a:expr, $b:expr)),+ $(,)?) => {
            $($op!($a, $b);)+
        };
    }

    pairs!(unpack8; (0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (10, 11), (12, 13), (14, 15),
        (16, 17), (18, 19), (20, 21), (22, 23), (24, 25), (26, 27), (28, 29), (30, 31));
    pairs!(unpack16; (0, 2), (1, 3), (4, 6), (5, 7), (8, 10), (9, 11), (12, 14), (13, 15),
        (16, 18), (17, 19), (20, 22), (21, 23), (24, 26), (25, 27), (28, 30), (29, 31));
    pairs!(unpack32; (0, 4), (2, 6), (1, 5), (3, 7), (8, 12), (10, 14), (9, 13), (11, 15),
        (16, 20), (18, 22), (17, 21), (19, 23), (24, 28), (26, 30), (25, 29), (27, 31));
    pairs!(unpack64; (0, 8), (4, 12), (2, 10), (6, 14), (1, 9), (5, 13), (3, 11), (7, 15),
        (16, 24), (20, 28), (18, 26), (22, 30), (17, 25), (21, 29), (19, 27), (23, 31));
    pairs!(unpack128; (0, 16), (8, 24), (4, 20), (12, 28), (2, 18), (10, 26), (6, 22), (14, 30),
        (1, 17), (9, 25), (5, 21), (13, 29), (3, 19), (11, 27), (7, 23), (15, 31));

    macro_rules! store_row {
        ($column:expr, $register:expr) => {
            _mm256_storeu_si256(out.as_mut_ptr().add($column * 32).cast(), rows[$register]);
        };
    }
    store_row!(0, 0);
    store_row!(1, 8);
    store_row!(2, 4);
    store_row!(3, 12);
    store_row!(4, 2);
    store_row!(5, 10);
    store_row!(6, 6);
    store_row!(7, 14);
    store_row!(8, 1);
    store_row!(9, 9);
    store_row!(10, 5);
    store_row!(11, 13);
    store_row!(12, 3);
    store_row!(13, 11);
    store_row!(14, 7);
    store_row!(15, 15);
    store_row!(16, 16);
    store_row!(17, 24);
    store_row!(18, 20);
    store_row!(19, 28);
    store_row!(20, 18);
    store_row!(21, 26);
    store_row!(22, 22);
    store_row!(23, 30);
    store_row!(24, 17);
    store_row!(25, 25);
    store_row!(26, 21);
    store_row!(27, 29);
    store_row!(28, 19);
    store_row!(29, 27);
    store_row!(30, 23);
    store_row!(31, 31);
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn transpose_16_neon(data: &[&[i8]], n: usize, out: &mut [i8]) {
    use std::arch::aarch64::*;

    let zero = vdupq_n_s8(0);
    let mut rows = [zero; 16];
    let row0 = 16 - n;
    for (dst, src) in rows[row0..].iter_mut().zip(data.iter().take(n)) {
        debug_assert!(src.len() >= 16);
        *dst = vld1q_s8(src.as_ptr());
    }
    for group in (0..16).step_by(2) {
        let pair = vtrnq_s8(rows[group], rows[group + 1]);
        rows[group] = pair.0;
        rows[group + 1] = pair.1;
    }
    for group in (0..16).step_by(4) {
        for offset in 0..2 {
            let pair = vtrnq_s16(
                vreinterpretq_s16_s8(rows[group + offset]),
                vreinterpretq_s16_s8(rows[group + offset + 2]),
            );
            rows[group + offset] = vreinterpretq_s8_s16(pair.0);
            rows[group + offset + 2] = vreinterpretq_s8_s16(pair.1);
        }
    }
    for group in (0..16).step_by(8) {
        for offset in 0..4 {
            let pair = vtrnq_s32(
                vreinterpretq_s32_s8(rows[group + offset]),
                vreinterpretq_s32_s8(rows[group + offset + 4]),
            );
            rows[group + offset] = vreinterpretq_s8_s32(pair.0);
            rows[group + offset + 4] = vreinterpretq_s8_s32(pair.1);
        }
    }
    for column in 0..8 {
        let low = vcombine_s8(
            vget_low_s8(rows[bit_reverse(column, 3)]),
            vget_low_s8(rows[bit_reverse(column, 3) + 8]),
        );
        let high = vcombine_s8(
            vget_high_s8(rows[bit_reverse(column, 3)]),
            vget_high_s8(rows[bit_reverse(column, 3) + 8]),
        );
        vst1q_s8(out.as_mut_ptr().add(column * 16), low);
        vst1q_s8(out.as_mut_ptr().add((column + 8) * 16), high);
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

    fn scalar_transpose(rows: &[Vec<i8>], width: usize) -> Vec<i8> {
        let mut out = vec![0; width * width];
        let row_offset = width - rows.len();
        for (row, input) in rows.iter().enumerate() {
            for column in 0..width {
                out[column * width + row_offset + row] = input[column];
            }
        }
        out
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
    fn transposes_match_scalar_for_every_batch_size() {
        for width in [16, 32] {
            for n in 0..=width {
                let rows = rows(width, n);
                let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
                let expected = scalar_transpose(&rows, width);
                let mut actual = vec![0x55; width * width];
                transpose(&refs, n, &mut actual, width);
                assert_eq!(actual, expected, "width={width}, n={n}");
            }
        }
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
