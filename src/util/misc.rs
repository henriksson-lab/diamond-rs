use crate::util::system::PATH_SEPARATOR;

pub fn div_up<T>(x: T, m: T) -> T
where
    T: Copy
        + From<u8>
        + std::ops::Add<Output = T>
        + std::ops::Sub<Output = T>
        + std::ops::Div<Output = T>,
{
    (x + (m - T::from(1))) / m
}

pub fn round_up<T>(x: T, m: T) -> T
where
    T: Copy
        + From<u8>
        + std::ops::Add<Output = T>
        + std::ops::Sub<Output = T>
        + std::ops::Div<Output = T>
        + std::ops::Mul<Output = T>,
{
    div_up(x, m) * m
}

pub fn extract_dir(s: &str) -> String {
    match s.rfind(PATH_SEPARATOR) {
        Some(pos) => s[..pos].to_string(),
        None => String::new(),
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sd {
    a: f64,
    q: f64,
    k: f64,
}

impl Default for Sd {
    fn default() -> Self {
        Self::new()
    }
}

impl Sd {
    pub fn new() -> Self {
        Self {
            a: 0.0,
            q: 0.0,
            k: 1.0,
        }
    }

    pub fn from_groups(groups: &[Sd]) -> Self {
        let mut out = Self {
            a: 0.0,
            q: 0.0,
            k: 0.0,
        };
        for group in groups {
            out.k += group.k;
            out.a += group.a * group.k;
            out.q += group.q;
        }
        out.a /= out.k;
        for group in groups {
            out.q += (group.a - out.a) * (group.a - out.a) * group.k;
        }
        out
    }

    pub fn add(&mut self, x: f64) {
        let d = x - self.a;
        self.q += (self.k - 1.0) / self.k * d * d;
        self.a += d / self.k;
        self.k += 1.0;
    }

    pub fn mean(&self) -> f64 {
        self.a
    }

    pub fn sd(&self) -> f64 {
        (self.q / (self.k - 1.0)).sqrt()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Matrix<T> {
    cols: i64,
    data: Vec<T>,
}

impl<T: Clone + Default> Matrix<T> {
    pub fn new(rows: i64, cols: i64) -> Self {
        Self {
            cols,
            data: vec![T::default(); (rows * cols) as usize],
        }
    }

    pub fn row(&self, i: i64) -> &[T] {
        let start = (i * self.cols) as usize;
        &self.data[start..start + self.cols as usize]
    }

    pub fn row_mut(&mut self, i: i64) -> &mut [T] {
        let start = (i * self.cols) as usize;
        &mut self.data[start..start + self.cols as usize]
    }
}

pub fn round_down_const<const N: i32>(x: i32) -> i32 {
    (x / N) * N
}

pub fn round_up_const<const N: i32>(x: i32) -> i32 {
    ((x + N - 1) / N) * N
}

pub fn print_char(c: u8) -> String {
    if c < 32 {
        format!("ASCII {}", c)
    } else {
        (c as char).to_string()
    }
}

pub fn percentage<T1, T2>(x: T2, y: T2) -> T1
where
    T1: From<f64>,
    T2: Copy + Into<f64>,
{
    T1::from(x.into() * 100.0 / y.into())
}

pub fn print_binary(mut x: u64) -> String {
    let mut out = String::with_capacity(64);
    for _ in 0..64 {
        out.push(if x & 1 != 0 { '1' } else { '0' });
        x >>= 1;
    }
    out
}

pub fn megabytes(x: usize) -> f64 {
    x as f64 / (1usize << 20) as f64
}

pub fn make_multiple<T>(x: T, m: T) -> T
where
    T: Copy
        + Eq
        + From<u8>
        + std::ops::Rem<Output = T>
        + std::ops::Sub<Output = T>
        + std::ops::Add<Output = T>
        + Bounded,
{
    if x % m == T::from(0) {
        return x;
    }
    let d = m - x % m;
    if T::max_value() - d < x {
        x
    } else {
        x + d
    }
}

pub trait Bounded: Copy + PartialOrd + std::ops::Sub<Output = Self> {
    fn max_value() -> Self;
}

macro_rules! impl_bounded {
    ($($t:ty),* $(,)?) => {
        $(
            impl Bounded for $t {
                fn max_value() -> Self {
                    <$t>::MAX
                }
            }
        )*
    };
}

impl_bounded!(u8, u16, u32, u64, usize, i8, i16, i32, i64, isize);

pub fn hex_print(x: &[u8]) -> String {
    let mut out = String::with_capacity(x.len() * 2);
    for b in x {
        out.push_str(&format!("{:02x}", b));
    }
    out
}

pub fn combine<T1: Clone, T2: Clone>(v1: &[T1], v2: &[T2]) -> Vec<(T1, T2)> {
    let mut r = Vec::with_capacity(v1.len());
    for i in 0..v1.len() {
        r.push((v1[i].clone(), v2[i].clone()));
    }
    r
}

#[derive(Debug, Clone)]
pub struct AsyncKeyMerger<'a, T, Key, GetKey> {
    slice: &'a [T],
    it: usize,
    end: usize,
    key: GetKey,
    _marker: std::marker::PhantomData<Key>,
}

impl<'a, T, Key, GetKey> AsyncKeyMerger<'a, T, Key, GetKey>
where
    Key: PartialEq,
    GetKey: Copy + Fn(&T) -> Key,
{
    pub fn new(slice: &'a [T], key: GetKey) -> Self {
        Self {
            slice,
            it: 0,
            end: slice.len(),
            key,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn increment(&mut self) -> (usize, usize) {
        if self.it == self.end {
            return (self.end, self.end);
        }
        let begin = self.it;
        self.it += 1;
        let key = (self.key)(&self.slice[begin]);
        while self.it != self.end && (self.key)(&self.slice[self.it]) == key {
            self.it += 1;
        }
        (begin, self.it)
    }
}

#[derive(Debug, Clone)]
pub struct KeyMergeIterator<'a, T, Key, GetKey> {
    slice: &'a [T],
    end: usize,
    begin: usize,
    key_end: usize,
    get_key: GetKey,
    key: Option<Key>,
    next_key: Option<Key>,
}

impl<'a, T, Key, GetKey> KeyMergeIterator<'a, T, Key, GetKey>
where
    Key: Copy + PartialEq,
    GetKey: Copy + Fn(&T) -> Key,
{
    pub fn new(slice: &'a [T], key: GetKey) -> Self {
        let mut out = Self {
            slice,
            end: slice.len(),
            begin: 0,
            key_end: 0,
            get_key: key,
            key: None,
            next_key: slice.first().map(key),
        };
        if !slice.is_empty() {
            out.increment();
        }
        out
    }

    pub fn increment(&mut self) {
        self.begin = self.key_end;
        if self.begin == self.end {
            return;
        }
        self.key = self.next_key;
        self.key_end += 1;
        while self.key_end != self.end {
            self.next_key = Some((self.get_key)(&self.slice[self.key_end]));
            if self.next_key != self.key {
                break;
            }
            self.key_end += 1;
        }
    }

    pub fn good(&self) -> bool {
        self.begin != self.end
    }

    pub fn begin(&self) -> usize {
        self.begin
    }

    pub fn end(&self) -> usize {
        self.key_end
    }

    pub fn range(&self) -> &'a [T] {
        &self.slice[self.begin..self.key_end]
    }

    pub fn key(&self) -> Key {
        self.key.unwrap()
    }

    pub fn count(&self) -> usize {
        self.key_end - self.begin
    }
}

pub fn merge_keys_by<T, Key, GetKey>(
    slice: &[T],
    key: GetKey,
) -> KeyMergeIterator<'_, T, Key, GetKey>
where
    Key: Copy + PartialEq,
    GetKey: Copy + Fn(&T) -> Key,
{
    KeyMergeIterator::new(slice, key)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexIterator {
    i: usize,
}

impl IndexIterator {
    pub fn new(i: usize) -> Self {
        Self { i }
    }
}

impl Iterator for IndexIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        let i = self.i;
        self.i += 1;
        Some(i)
    }
}

pub fn log2_approximate(x: f32) -> f32 {
    let mut bits = x.to_bits();
    let mut log_2 = ((bits >> 23) & 255) as f32 - 128.0;
    bits &= !(255 << 23);
    bits += 127 << 23;
    let z = f32::from_bits(bits);
    log_2 += (-0.34484843f32 * z + 2.02466578f32) * z - 0.67487759f32;
    log_2
}

pub fn to_upper_case(s: &str) -> String {
    s.chars().map(|c| c.to_ascii_uppercase()).collect()
}

pub fn to_lower_case(s: &str) -> String {
    s.chars().map(|c| c.to_ascii_lowercase()).collect()
}

pub fn safe_cast_i64_to_i32(value: i64) -> Result<i32, String> {
    i32::try_from(value).map_err(|_| "safe_cast: out of range (signed -> signed)".to_string())
}

pub fn safe_cast_i64_to_u32(value: i64) -> Result<u32, String> {
    if value < 0 {
        Err("safe_cast: negative value to unsigned".to_string())
    } else {
        u32::try_from(value).map_err(|_| "safe_cast: overflow (signed -> unsigned)".to_string())
    }
}

pub fn safe_cast_u64_to_i32(value: u64) -> Result<i32, String> {
    i32::try_from(value).map_err(|_| "safe_cast: overflow (unsigned -> signed)".to_string())
}

pub fn safe_cast_f64_to_i32(value: f64) -> Result<i32, String> {
    if !value.is_finite() {
        return Err("safe_cast: non-finite value (NaN/Inf)".to_string());
    }
    if value < i32::MIN as f64 || value > i32::MAX as f64 {
        return Err("safe_cast: out of range (float -> int)".to_string());
    }
    Ok(value as i32)
}

pub trait SafeCastFrom<From>: Sized {
    fn safe_cast_from(value: From) -> Result<Self, String>;
}

macro_rules! impl_safe_cast_int_to_int {
    ($to:ty, $from:ty) => {
        impl SafeCastFrom<$from> for $to {
            fn safe_cast_from(value: $from) -> Result<Self, String> {
                <$to>::try_from(value).map_err(|_| {
                    let from_signed = (<$from>::MIN as i128) < 0;
                    let to_signed = (<$to>::MIN as i128) < 0;
                    match (from_signed, to_signed) {
                        (true, true) => "safe_cast: out of range (signed -> signed)",
                        (true, false) if (value as i128) < 0 => {
                            "safe_cast: negative value to unsigned"
                        }
                        (true, false) => "safe_cast: overflow (signed -> unsigned)",
                        (false, true) => "safe_cast: overflow (unsigned -> signed)",
                        (false, false) => "safe_cast: overflow (unsigned -> unsigned)",
                    }
                    .to_string()
                })
            }
        }
    };
}

macro_rules! impl_safe_cast_float_to_int {
    ($to:ty, $from:ty) => {
        impl SafeCastFrom<$from> for $to {
            fn safe_cast_from(value: $from) -> Result<Self, String> {
                if !value.is_finite() {
                    return Err("safe_cast: non-finite value (NaN/Inf)".to_string());
                }
                if (<$to>::MIN as i128) >= 0 && value < 0.0 {
                    return Err("safe_cast: negative value to unsigned".to_string());
                }
                if value < <$to>::MIN as $from || value > <$to>::MAX as $from {
                    return Err("safe_cast: out of range (float -> int)".to_string());
                }
                Ok(value as $to)
            }
        }
    };
}

impl_safe_cast_int_to_int!(i8, i8);
impl_safe_cast_int_to_int!(i8, i16);
impl_safe_cast_int_to_int!(i8, i32);
impl_safe_cast_int_to_int!(i8, i64);
impl_safe_cast_int_to_int!(i8, u8);
impl_safe_cast_int_to_int!(i8, u16);
impl_safe_cast_int_to_int!(i8, u32);
impl_safe_cast_int_to_int!(i8, u64);
impl_safe_cast_int_to_int!(i16, i8);
impl_safe_cast_int_to_int!(i16, i16);
impl_safe_cast_int_to_int!(i16, i32);
impl_safe_cast_int_to_int!(i16, i64);
impl_safe_cast_int_to_int!(i16, u8);
impl_safe_cast_int_to_int!(i16, u16);
impl_safe_cast_int_to_int!(i16, u32);
impl_safe_cast_int_to_int!(i16, u64);
impl_safe_cast_int_to_int!(i32, i8);
impl_safe_cast_int_to_int!(i32, i16);
impl_safe_cast_int_to_int!(i32, i32);
impl_safe_cast_int_to_int!(i32, i64);
impl_safe_cast_int_to_int!(i32, u8);
impl_safe_cast_int_to_int!(i32, u16);
impl_safe_cast_int_to_int!(i32, u32);
impl_safe_cast_int_to_int!(i32, u64);
impl_safe_cast_int_to_int!(i64, i8);
impl_safe_cast_int_to_int!(i64, i16);
impl_safe_cast_int_to_int!(i64, i32);
impl_safe_cast_int_to_int!(i64, i64);
impl_safe_cast_int_to_int!(i64, u8);
impl_safe_cast_int_to_int!(i64, u16);
impl_safe_cast_int_to_int!(i64, u32);
impl_safe_cast_int_to_int!(i64, u64);
impl_safe_cast_int_to_int!(u8, i8);
impl_safe_cast_int_to_int!(u8, i16);
impl_safe_cast_int_to_int!(u8, i32);
impl_safe_cast_int_to_int!(u8, i64);
impl_safe_cast_int_to_int!(u8, u8);
impl_safe_cast_int_to_int!(u8, u16);
impl_safe_cast_int_to_int!(u8, u32);
impl_safe_cast_int_to_int!(u8, u64);
impl_safe_cast_int_to_int!(u16, i8);
impl_safe_cast_int_to_int!(u16, i16);
impl_safe_cast_int_to_int!(u16, i32);
impl_safe_cast_int_to_int!(u16, i64);
impl_safe_cast_int_to_int!(u16, u8);
impl_safe_cast_int_to_int!(u16, u16);
impl_safe_cast_int_to_int!(u16, u32);
impl_safe_cast_int_to_int!(u16, u64);
impl_safe_cast_int_to_int!(u32, i8);
impl_safe_cast_int_to_int!(u32, i16);
impl_safe_cast_int_to_int!(u32, i32);
impl_safe_cast_int_to_int!(u32, i64);
impl_safe_cast_int_to_int!(u32, u8);
impl_safe_cast_int_to_int!(u32, u16);
impl_safe_cast_int_to_int!(u32, u32);
impl_safe_cast_int_to_int!(u32, u64);
impl_safe_cast_int_to_int!(u64, i8);
impl_safe_cast_int_to_int!(u64, i16);
impl_safe_cast_int_to_int!(u64, i32);
impl_safe_cast_int_to_int!(u64, i64);
impl_safe_cast_int_to_int!(u64, u8);
impl_safe_cast_int_to_int!(u64, u16);
impl_safe_cast_int_to_int!(u64, u32);
impl_safe_cast_int_to_int!(u64, u64);

impl_safe_cast_float_to_int!(i8, f32);
impl_safe_cast_float_to_int!(i16, f32);
impl_safe_cast_float_to_int!(i32, f32);
impl_safe_cast_float_to_int!(i64, f32);
impl_safe_cast_float_to_int!(u8, f32);
impl_safe_cast_float_to_int!(u16, f32);
impl_safe_cast_float_to_int!(u32, f32);
impl_safe_cast_float_to_int!(u64, f32);
impl_safe_cast_float_to_int!(i8, f64);
impl_safe_cast_float_to_int!(i16, f64);
impl_safe_cast_float_to_int!(i32, f64);
impl_safe_cast_float_to_int!(i64, f64);
impl_safe_cast_float_to_int!(u8, f64);
impl_safe_cast_float_to_int!(u16, f64);
impl_safe_cast_float_to_int!(u32, f64);
impl_safe_cast_float_to_int!(u64, f64);

pub fn safe_cast<To, From>(value: From) -> Result<To, String>
where
    To: SafeCastFrom<From>,
{
    To::safe_cast_from(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rounding_dirs_and_case() {
        assert_eq!(div_up(10u32, 3), 4);
        assert_eq!(round_up(10u32, 3), 12);
        assert_eq!(round_down_const::<8>(31), 24);
        assert_eq!(round_up_const::<8>(31), 32);
        assert_eq!(extract_dir("/tmp/x"), "/tmp");
        assert_eq!(to_upper_case("aBc"), "ABC");
        assert_eq!(to_lower_case("aBc"), "abc");
    }

    #[test]
    fn test_sd_matrix_and_formatting() {
        let mut sd = Sd::new();
        sd.add(2.0);
        sd.add(4.0);
        assert!((sd.mean() - 3.0).abs() < 1e-12);
        assert!((sd.sd() - 1.0).abs() < 1e-12);
        let combined = Sd::from_groups(&[sd, sd]);
        assert!((combined.mean() - 3.0).abs() < 1e-12);
        let mut m = Matrix::<i32>::new(2, 3);
        m.row_mut(1)[2] = 7;
        assert_eq!(m.row(1)[2], 7);
        assert_eq!(print_char(7), "ASCII 7");
        assert_eq!(print_char(b'A'), "A");
        assert_eq!(hex_print(&[0, 15, 255]), "000fff");
    }

    #[test]
    fn test_binary_multiple_log2_and_casts() {
        assert_eq!(&print_binary(5)[..4], "1010");
        assert_eq!(megabytes(1 << 20), 1.0);
        assert_eq!(make_multiple(10u32, 4), 12);
        assert_eq!(combine(&[1, 2], &["a", "b"]), vec![(1, "a"), (2, "b")]);
        assert!((log2_approximate(8.0) - 3.0).abs() < 0.05);
        assert_eq!(safe_cast_i64_to_i32(7).unwrap(), 7);
        assert!(safe_cast_i64_to_i32(i64::MAX).is_err());
        assert!(safe_cast_i64_to_u32(-1).is_err());
        assert!(safe_cast_u64_to_i32(u64::MAX).is_err());
        assert_eq!(safe_cast_f64_to_i32(1.9).unwrap(), 1);
        assert!(safe_cast_f64_to_i32(f64::INFINITY).is_err());
        assert_eq!(safe_cast::<i32, i64>(7).unwrap(), 7);
        assert_eq!(
            safe_cast::<u32, i64>(-1).unwrap_err(),
            "safe_cast: negative value to unsigned"
        );
        assert_eq!(
            safe_cast::<i32, u64>(u64::MAX).unwrap_err(),
            "safe_cast: overflow (unsigned -> signed)"
        );
        assert_eq!(safe_cast::<u32, f64>(1.9).unwrap(), 1);
        assert_eq!(
            safe_cast::<u32, f64>(-1.0).unwrap_err(),
            "safe_cast: negative value to unsigned"
        );
        assert_eq!(
            safe_cast::<i32, f64>(f64::INFINITY).unwrap_err(),
            "safe_cast: non-finite value (NaN/Inf)"
        );
    }

    #[test]
    fn test_key_mergers() {
        let rows = [(1, 'a'), (1, 'b'), (2, 'c'), (4, 'd'), (4, 'e')];
        let mut merger = merge_keys_by(&rows, |x| x.0);
        assert!(merger.good());
        assert_eq!(merger.key(), 1);
        assert_eq!(merger.begin(), 0);
        assert_eq!(merger.end(), 2);
        assert_eq!(merger.count(), 2);
        assert_eq!(merger.range(), &[(1, 'a'), (1, 'b')]);
        merger.increment();
        assert_eq!(merger.key(), 2);
        assert_eq!(merger.range(), &[(2, 'c')]);
        merger.increment();
        assert_eq!(merger.key(), 4);
        assert_eq!(merger.range(), &[(4, 'd'), (4, 'e')]);
        merger.increment();
        assert!(!merger.good());

        let mut async_merger = AsyncKeyMerger::new(&rows, |x| x.0);
        assert_eq!(async_merger.increment(), (0, 2));
        assert_eq!(async_merger.increment(), (2, 3));
        assert_eq!(async_merger.increment(), (3, 5));
        assert_eq!(async_merger.increment(), (5, 5));
    }

    #[test]
    fn test_index_iterator() {
        let mut it = IndexIterator::new(3);
        assert_eq!(it.next(), Some(3));
        assert_eq!(it.next(), Some(4));
    }
}
