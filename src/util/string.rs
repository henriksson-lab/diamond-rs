use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TokenizerException {
    msg: String,
}

impl TokenizerException {
    pub fn new() -> Self {
        Self {
            msg: "Tokenizer Exception".to_string(),
        }
    }

    pub fn with_msg(msg: &str) -> Self {
        Self {
            msg: msg.to_string(),
        }
    }
}

impl std::fmt::Display for TokenizerException {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.msg)
    }
}

impl std::error::Error for TokenizerException {}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Skip;

pub trait Delimiter {
    fn scan(&self, s: &str, p: usize) -> Option<(usize, usize)>;
    fn next(&self, s: &str, p: usize) -> Result<Option<usize>, TokenizerException>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CharDelimiter {
    pub c: char,
}

impl CharDelimiter {
    pub fn new(c: char) -> Self {
        Self { c }
    }
}

impl Delimiter for CharDelimiter {
    fn scan(&self, s: &str, p: usize) -> Option<(usize, usize)> {
        s[p..].find(self.c).map(|q| {
            let begin = p + q;
            (begin, begin + self.c.len_utf8())
        })
    }

    fn next(&self, s: &str, p: usize) -> Result<Option<usize>, TokenizerException> {
        if s[p..].starts_with(self.c) {
            Ok(Some(p + self.c.len_utf8()))
        } else if p == s.len() {
            Ok(None)
        } else {
            Err(TokenizerException::new())
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StringDelimiter {
    pub s: String,
    pub len: usize,
}

impl StringDelimiter {
    pub fn new(s: &str) -> Self {
        Self {
            s: s.to_string(),
            len: s.len(),
        }
    }
}

impl Delimiter for StringDelimiter {
    fn scan(&self, s: &str, p: usize) -> Option<(usize, usize)> {
        s[p..].find(&self.s).map(|q| {
            let begin = p + q;
            (begin, begin + self.len)
        })
    }

    fn next(&self, s: &str, p: usize) -> Result<Option<usize>, TokenizerException> {
        if s[p..].starts_with(&self.s) {
            Ok(Some(p + self.len))
        } else if p == s.len() {
            Ok(None)
        } else {
            Err(TokenizerException::new())
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StringDelimiters {
    pub s: Vec<String>,
}

impl StringDelimiters {
    pub fn new(s: &[&str]) -> Self {
        Self {
            s: s.iter().map(|x| x.to_string()).collect(),
        }
    }
}

impl Delimiter for StringDelimiters {
    fn scan(&self, s: &str, p: usize) -> Option<(usize, usize)> {
        for delimiter in &self.s {
            if let Some(q) = s[p..].find(delimiter) {
                let begin = p + q;
                return Some((begin, begin + delimiter.len()));
            }
        }
        None
    }

    fn next(&self, s: &str, p: usize) -> Result<Option<usize>, TokenizerException> {
        for delimiter in &self.s {
            if s[p..].starts_with(delimiter) {
                return Ok(Some(p + delimiter.len()));
            }
        }
        if p == s.len() {
            Ok(None)
        } else {
            Err(TokenizerException::new())
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Tokenizer<'a, Delimiters> {
    s: &'a str,
    p: Option<usize>,
    delimiters: Delimiters,
}

impl<'a, Delimiters> Tokenizer<'a, Delimiters>
where
    Delimiters: Delimiter,
{
    pub fn new(s: &'a str, delimiters: Delimiters) -> Self {
        Self {
            s,
            p: Some(0),
            delimiters,
        }
    }

    pub fn ptr(&self) -> Option<usize> {
        self.p
    }

    pub fn set(&mut self, ptr: Option<usize>) {
        self.p = ptr;
    }

    pub fn read_skip(&mut self, _skip: Skip) -> Result<&mut Self, TokenizerException> {
        let p = self.p.ok_or_else(TokenizerException::new)?;
        self.p = self.delimiters.scan(self.s, p).map(|(_, next)| next);
        Ok(self)
    }

    pub fn read_string(&mut self) -> Result<String, TokenizerException> {
        let p = self.p.ok_or_else(TokenizerException::new)?;
        let out = if let Some((begin, next)) = self.delimiters.scan(self.s, p) {
            self.p = Some(next);
            self.s[p..begin].to_string()
        } else {
            self.p = None;
            self.s[p..].to_string()
        };
        Ok(out)
    }

    pub fn read_i64(&mut self) -> Result<i64, TokenizerException> {
        let (token, end) = self.read_integer_prefix()?;
        let value = token
            .parse::<i64>()
            .map_err(|_| TokenizerException::new())?;
        self.p = self.delimiters.next(self.s, end)?;
        Ok(value)
    }

    pub fn read_u64(&mut self) -> Result<u64, TokenizerException> {
        let (token, end) = self.read_integer_prefix()?;
        let value = token
            .parse::<i64>()
            .map(|x| x as u64)
            .map_err(|_| TokenizerException::new())?;
        self.p = self.delimiters.next(self.s, end)?;
        Ok(value)
    }

    pub fn read_i32(&mut self) -> Result<i32, TokenizerException> {
        self.read_i64().map(|x| x as i32)
    }

    pub fn read_u32(&mut self) -> Result<u32, TokenizerException> {
        self.read_i64().map(|x| x as u32)
    }

    pub fn read_f64(&mut self) -> Result<f64, TokenizerException> {
        let (token, end) = self.read_float_prefix("No token left", "Unable to parse double")?;
        let value = token
            .parse::<f64>()
            .map_err(|_| TokenizerException::with_msg("Unable to parse double"))?;
        self.p = self.delimiters.next(self.s, end)?;
        Ok(value)
    }

    pub fn read_f32(&mut self) -> Result<f32, TokenizerException> {
        let (token, end) = self.read_float_prefix("No token left", "Unable to parse float")?;
        let value = token
            .parse::<f32>()
            .map_err(|_| TokenizerException::with_msg("Unable to parse float"))?;
        self.p = self.delimiters.next(self.s, end)?;
        Ok(value)
    }

    pub fn good(&self) -> bool {
        self.p.is_some_and(|p| p < self.s.len())
    }

    pub fn skip_to(&mut self, c: char) {
        self.p = self
            .p
            .and_then(|p| self.s[p..].find(c).map(|q| p + q + c.len_utf8()));
    }

    pub fn getline(&self) -> String {
        let p = self.p.unwrap_or(self.s.len());
        let q = self.s[p..]
            .find('\n')
            .map(|i| p + i)
            .unwrap_or(self.s.len());
        self.s[p..q].to_string()
    }

    fn read_integer_prefix(&self) -> Result<(&str, usize), TokenizerException> {
        if !self.good() {
            return Err(TokenizerException::new());
        }
        let p = self.p.unwrap();
        let mut end = p;
        let bytes = self.s.as_bytes();
        if end < self.s.len() && (bytes[end] == b'+' || bytes[end] == b'-') {
            end += 1;
        }
        let digits_begin = end;
        while end < self.s.len() && bytes[end].is_ascii_digit() {
            end += 1;
        }
        if end == digits_begin {
            return Err(TokenizerException::new());
        }
        Ok((&self.s[p..end], end))
    }

    fn read_float_prefix(
        &self,
        empty_msg: &str,
        parse_msg: &str,
    ) -> Result<(&str, usize), TokenizerException> {
        if !self.good() {
            return Err(TokenizerException::with_msg(empty_msg));
        }
        let p = self.p.unwrap();
        let bytes = self.s.as_bytes();
        let mut end = p;
        if end < self.s.len() && (bytes[end] == b'+' || bytes[end] == b'-') {
            end += 1;
        }
        let mut seen_digit = false;
        while end < self.s.len() && bytes[end].is_ascii_digit() {
            end += 1;
            seen_digit = true;
        }
        if end < self.s.len() && bytes[end] == b'.' {
            end += 1;
            while end < self.s.len() && bytes[end].is_ascii_digit() {
                end += 1;
                seen_digit = true;
            }
        }
        if !seen_digit {
            return Err(TokenizerException::with_msg(parse_msg));
        }
        if end < self.s.len() && (bytes[end] == b'e' || bytes[end] == b'E') {
            let exp = end;
            end += 1;
            if end < self.s.len() && (bytes[end] == b'+' || bytes[end] == b'-') {
                end += 1;
            }
            let exp_digits = end;
            while end < self.s.len() && bytes[end].is_ascii_digit() {
                end += 1;
            }
            if end == exp_digits {
                end = exp;
            }
        }
        Ok((&self.s[p..end], end))
    }
}

pub fn ends_with(s: &str, t: &str) -> bool {
    if s.len() < t.len() {
        return false;
    }
    &s[s.len() - t.len()..] == t
}

pub fn rstrip(s: &str, t: &str) -> String {
    if s.len() < t.len() {
        return s.to_string();
    }
    if &s[s.len() - t.len()..] == t {
        s[..s.len() - t.len()].to_string()
    } else {
        s.to_string()
    }
}

pub fn max_len(s: &[&str]) -> usize {
    let mut l = 0;
    for item in s {
        l = l.max(item.len());
    }
    l
}

pub fn convert_size(mut size: usize) -> String {
    const SIZES: [&str; 6] = ["B", "KB", "MB", "GB", "TB", "PB"];
    let mut div = 0usize;
    let mut rem = 0usize;

    while size >= 1024 && div < SIZES.len() {
        rem = size % 1024;
        div += 1;
        size /= 1024;
    }

    format!("{:.1} {}", size as f64 + rem as f64 / 1024.0, SIZES[div])
}

pub fn join<T: std::fmt::Display>(sep: &str, values: &[T]) -> String {
    let mut out = String::new();
    let mut first = true;
    for value in values {
        if !first {
            out.push_str(sep);
        }
        first = false;
        out.push_str(&value.to_string());
    }
    out
}

pub fn format_double(x: f64) -> String {
    if x >= 100.0 {
        format!("{}", x.floor() as i64)
    } else {
        let i = (x * 10.0).round() as i64;
        format!("{}.{}", i / 10, i % 10)
    }
}

pub fn replace(s: &str, a: char, b: char) -> String {
    s.chars().map(|c| if c == a { b } else { c }).collect()
}

pub fn ratio_percentage_f64(x: f64, y: f64) -> String {
    format!("{:.0}/{:.0} ({:.2}%)", x, y, x / y * 100.0)
}

pub fn ratio_percentage_usize(x: usize, y: usize) -> String {
    ratio_percentage_f64(x as f64, y as f64)
}

pub trait RatioPercentageArg {
    fn to_f64(self) -> f64;
}

impl RatioPercentageArg for f64 {
    fn to_f64(self) -> f64 {
        self
    }
}

impl RatioPercentageArg for usize {
    fn to_f64(self) -> f64 {
        self as f64
    }
}

pub fn ratio_percentage<T>(x: T, y: T) -> String
where
    T: RatioPercentageArg + Copy,
{
    ratio_percentage_f64(x.to_f64(), y.to_f64())
}

pub fn interpret_number(s: &str) -> Result<i64, String> {
    let trimmed = s.trim_start();
    let mut split = trimmed.len();
    for (i, c) in trimmed.char_indices() {
        if !(c.is_ascii_digit() || c == '.' || c == '+' || c == '-') {
            split = i;
            break;
        }
    }
    let n: f64 = trimmed[..split]
        .parse()
        .map_err(|_| format!("Invalid number format: {}", s))?;
    let suffix = trimmed[split..].chars().next().ok_or_else(|| {
        format!(
            "Missing size specifier in number: {}. Permitted values: T, G, M, K",
            s
        )
    })?;
    let mult = match suffix {
        'T' | 't' => 1e12,
        'G' | 'g' => 1e9,
        'M' | 'm' => 1e6,
        'K' | 'k' => 1e3,
        _ => {
            return Err(format!(
                "Invalid size specifier ({}) in number: {}. Permitted suffixes: T, G, M, K",
                suffix, s
            ))
        }
    };
    if trimmed[split + suffix.len_utf8()..]
        .chars()
        .any(|c| !c.is_ascii_whitespace())
    {
        return Err(format!("Invalid number format: {}", s));
    }
    Ok((n * mult) as i64)
}

pub fn tokenize(str_: &str, delimiters: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut i = 0usize;
    let bytes = str_.as_bytes();
    while i < bytes.len() {
        while i < bytes.len() && delimiters.as_bytes().contains(&bytes[i]) {
            i += 1;
        }
        let start = i;
        while i < bytes.len() && !delimiters.as_bytes().contains(&bytes[i]) {
            i += 1;
        }
        if i > start {
            out.push(str_[start..i].to_string());
        }
    }
    if out.is_empty() {
        out.push(String::new());
    }
    out
}

pub fn parse_csv(s: &str) -> BTreeSet<i32> {
    let mut r = BTreeSet::new();
    let t = tokenize(s, ",");
    for item in t {
        if !item.is_empty() {
            r.insert(item.parse::<i32>().unwrap_or(0));
        }
    }
    r
}

pub fn convert_string_i64(s: &str) -> Result<i64, String> {
    let trimmed = s.trim_start_matches(|c: char| c.is_ascii_whitespace());
    let value = trimmed
        .parse::<i64>()
        .map_err(|_| format!("Error converting integer value: {}", s))?;
    if value == i64::MAX || value == i64::MIN {
        return Err(format!("Error converting integer value: {}", s));
    }
    Ok(value)
}

pub fn convert_string_i32(s: &str) -> Result<i32, String> {
    convert_string_i64(s).and_then(|value| {
        i32::try_from(value).map_err(|_| format!("Error converting integer value: {}", s))
    })
}

pub fn convert_string_u64(s: &str) -> Result<u64, String> {
    let trimmed = s.trim_start_matches(|c: char| c.is_ascii_whitespace());
    let value = if let Some(rest) = trimmed.strip_prefix('-') {
        let magnitude = rest
            .parse::<u64>()
            .map_err(|_| format!("Error converting integer value: {}", s))?;
        0u64.wrapping_sub(magnitude)
    } else if let Some(rest) = trimmed.strip_prefix('+') {
        rest.parse::<u64>()
            .map_err(|_| format!("Error converting integer value: {}", s))?
    } else {
        trimmed
            .parse::<u64>()
            .map_err(|_| format!("Error converting integer value: {}", s))?
    };
    if value == u64::MAX {
        return Err(format!("Error converting integer value: {}", s));
    }
    Ok(value)
}

pub fn convert_string_u32(s: &str) -> Result<u32, String> {
    convert_string_i64(s).and_then(|value| {
        u32::try_from(value).map_err(|_| format!("Error converting integer value: {}", s))
    })
}

pub trait ConvertString: Sized {
    fn convert_string(s: &str) -> Result<Self, String>;
}

impl ConvertString for i64 {
    fn convert_string(s: &str) -> Result<Self, String> {
        convert_string_i64(s)
    }
}

impl ConvertString for i32 {
    fn convert_string(s: &str) -> Result<Self, String> {
        convert_string_i32(s)
    }
}

impl ConvertString for u64 {
    fn convert_string(s: &str) -> Result<Self, String> {
        convert_string_u64(s)
    }
}

impl ConvertString for u32 {
    fn convert_string(s: &str) -> Result<Self, String> {
        convert_string_u32(s)
    }
}

pub fn convert_string<T>(s: &str) -> Result<T, String>
where
    T: ConvertString,
{
    T::convert_string(s)
}

pub fn format_number(number: f64) -> String {
    const SUFFIXES: [&str; 7] = ["", "K", "M", "G", "T", "P", "E"];

    if number == 0.0 {
        return "0".to_string();
    }

    let is_negative = number < 0.0;
    let abs_number = number.abs();

    let mut suffix_index = 0usize;
    if abs_number >= 1000.0 {
        let exp = (abs_number.log10() / 3.0) as usize;
        suffix_index = exp.min(SUFFIXES.len() - 1);
    }

    let mut divisor = 1000.0f64.powi(suffix_index as i32);
    let mut formatted_number = abs_number / divisor;
    let mut num_str = format!("{:.2}", formatted_number);
    while num_str.ends_with('0') {
        num_str.pop();
    }
    if num_str.ends_with('.') {
        num_str.pop();
    }

    if formatted_number >= 999.995 && suffix_index < SUFFIXES.len() - 1 {
        suffix_index += 1;
        divisor *= 1000.0;
        formatted_number = abs_number / divisor;
        num_str = format!("{:.2}", formatted_number);
        while num_str.ends_with('0') {
            num_str.pop();
        }
        if num_str.ends_with('.') {
            num_str.pop();
        }
    }

    format!(
        "{}{}{}",
        if is_negative { "-" } else { "" },
        num_str,
        SUFFIXES[suffix_index]
    )
}

pub trait FormatNumber {
    fn to_format_f64(self) -> f64;
}

impl FormatNumber for f64 {
    fn to_format_f64(self) -> f64 {
        self
    }
}

impl FormatNumber for u64 {
    fn to_format_f64(self) -> f64 {
        self as f64
    }
}

impl FormatNumber for u32 {
    fn to_format_f64(self) -> f64 {
        self as f64
    }
}

impl FormatNumber for u16 {
    fn to_format_f64(self) -> f64 {
        self as f64
    }
}

impl FormatNumber for u8 {
    fn to_format_f64(self) -> f64 {
        self as f64
    }
}

impl FormatNumber for usize {
    fn to_format_f64(self) -> f64 {
        self as f64
    }
}

impl FormatNumber for i64 {
    fn to_format_f64(self) -> f64 {
        self as f64
    }
}

impl FormatNumber for i32 {
    fn to_format_f64(self) -> f64 {
        self as f64
    }
}

impl FormatNumber for i16 {
    fn to_format_f64(self) -> f64 {
        self as f64
    }
}

impl FormatNumber for i8 {
    fn to_format_f64(self) -> f64 {
        self as f64
    }
}

impl FormatNumber for isize {
    fn to_format_f64(self) -> f64 {
        self as f64
    }
}

pub fn format<T>(number: T) -> String
where
    T: FormatNumber,
{
    format_number(number.to_format_f64())
}

pub fn trim_cr<'a>(begin: &'a str, end: usize) -> &'a str {
    let mut end = end;
    if end > 0 && begin.as_bytes()[end - 1] == b'\r' {
        end -= 1;
    }
    &begin[..end]
}

pub fn find_newline(begin: &str) -> (usize, i32) {
    match begin.as_bytes().iter().position(|&b| b == b'\n') {
        Some(p) => {
            if p > 0 && begin.as_bytes()[p - 1] == b'\r' {
                (p + 1, 2)
            } else {
                (p + 1, 1)
            }
        }
        None => (begin.len(), 0),
    }
}

pub fn is_whitespace(s: &str) -> bool {
    for c in s.chars() {
        if c != ' ' && c != '\n' && c != '\r' {
            return false;
        }
    }
    true
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CharTokenizer<'a> {
    delimiter: u8,
    ptr: usize,
    end: usize,
    s: &'a str,
}

impl<'a> CharTokenizer<'a> {
    pub fn new(delimiter: char) -> Self {
        Self {
            delimiter: delimiter as u8,
            ptr: 0,
            end: 0,
            s: "",
        }
    }

    pub fn reset(&mut self, s: &'a str) {
        self.ptr = 0;
        self.end = s.len();
        self.s = s;
    }

    pub fn good(&self) -> bool {
        self.ptr < self.end
    }

    pub fn ptr(&self) -> usize {
        self.ptr
    }

    pub fn current(&self) -> String {
        let q = self.s.as_bytes()[self.ptr..self.end]
            .iter()
            .position(|&b| b == self.delimiter)
            .map(|i| self.ptr + i)
            .unwrap_or(self.end);
        self.s[self.ptr..q].to_string()
    }

    pub fn increment(&mut self) -> &mut Self {
        self.ptr = self.s.as_bytes()[self.ptr..self.end]
            .iter()
            .position(|&b| b == self.delimiter)
            .map(|i| self.ptr + i)
            .unwrap_or(self.end);
        if self.ptr < self.end {
            self.ptr += 1;
        }
        self
    }

    pub fn read_record(&mut self) -> (bool, String) {
        if self.ptr >= self.end {
            return (false, String::new());
        }
        let b = self.ptr;
        let found = self.s.as_bytes()[self.ptr..self.end]
            .iter()
            .position(|&c| c == b'\n')
            .map(|i| self.ptr + i);
        self.ptr = found.unwrap_or(self.end);
        if found.is_none() {
            return (true, self.s[b..self.end].to_string());
        }
        if self.ptr > b && self.s.as_bytes()[self.ptr - 1] == b'\r' {
            self.ptr += usize::from(self.ptr < self.end);
            (true, self.s[b..self.ptr - 2].to_string())
        } else {
            self.ptr += usize::from(self.ptr < self.end);
            let end = self.ptr.saturating_sub(1).max(b);
            (true, self.s[b..end].to_string())
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TokenIterator<'a, const DELIMITER: u8> {
    ptr: usize,
    end: usize,
    s: &'a str,
}

impl<'a, const DELIMITER: u8> TokenIterator<'a, DELIMITER> {
    pub fn new(s: &'a str) -> Self {
        Self {
            ptr: 0,
            end: s.len(),
            s,
        }
    }

    pub fn good(&self) -> bool {
        self.ptr < self.end
    }

    pub fn current(&self) -> String {
        let q = self.s.as_bytes()[self.ptr..self.end]
            .iter()
            .position(|&b| b == DELIMITER)
            .map(|i| self.ptr + i)
            .unwrap_or(self.end);
        self.s[self.ptr..q].to_string()
    }

    pub fn increment(&mut self) -> &mut Self {
        self.ptr = self.s.as_bytes()[self.ptr..self.end]
            .iter()
            .position(|&b| b == DELIMITER)
            .map(|i| self.ptr + i)
            .unwrap_or(self.end);
        if self.ptr < self.end {
            self.ptr += 1;
        }
        self
    }

    pub fn ptr(&self) -> usize {
        self.ptr
    }
}

pub type TabIterator<'a> = TokenIterator<'a, b'\t'>;
pub type LineIterator<'a> = TokenIterator<'a, b'\n'>;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FastaTokenizer<'a> {
    ptr: usize,
    end: usize,
    s: &'a str,
}

impl<'a> FastaTokenizer<'a> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn reset(&mut self, s: &'a str) {
        self.ptr = 0;
        self.end = s.len();
        self.s = s;
    }

    pub fn good(&self) -> bool {
        self.ptr < self.end
    }

    pub fn current(&self) -> String {
        if self.s.as_bytes()[self.ptr] == b'>' {
            let i = self.s.as_bytes()[self.ptr..self.end]
                .iter()
                .position(|&b| b == b'\n')
                .map(|i| self.ptr + i)
                .unwrap_or(self.end);
            trim_cr(&self.s[self.ptr + 1..], i - self.ptr - 1).to_string()
        } else {
            let mut r = String::new();
            let mut p = self.ptr;
            while p < self.end {
                let i = self.s.as_bytes()[p..self.end]
                    .iter()
                    .position(|&b| b == b'\n')
                    .map(|i| p + i)
                    .unwrap_or(self.end);
                r.push_str(trim_cr(&self.s[p..], i - p));
                p = i + usize::from(i < self.end);
            }
            r
        }
    }

    pub fn increment(&mut self) -> &mut Self {
        if self.s.as_bytes()[self.ptr] == b'>' {
            self.ptr = self.s.as_bytes()[self.ptr..self.end]
                .iter()
                .position(|&b| b == b'\n')
                .map(|i| self.ptr + i + 1)
                .unwrap_or(self.end);
        } else {
            self.ptr = self.end;
        }
        self
    }

    pub fn read_record(&mut self) -> (bool, String) {
        if self.ptr >= self.end {
            return (false, String::new());
        }
        let p = self.ptr;
        self.ptr = self.s[self.ptr..self.end]
            .find("\n>")
            .map(|i| self.ptr + i)
            .unwrap_or(self.end);
        (true, self.s[p..self.ptr].to_string())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MalformedFastqRecord {
    line: i64,
}

impl MalformedFastqRecord {
    pub fn new(line: i64) -> Self {
        Self { line }
    }
}

impl std::fmt::Display for MalformedFastqRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.line == -1 {
            f.write_str("Malformed FASTQ record")
        } else {
            write!(f, "Malformed FASTQ record on line {}", self.line)
        }
    }
}

impl std::error::Error for MalformedFastqRecord {}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FastqTokenizer<'a> {
    ptr: usize,
    end: usize,
    s: &'a str,
}

impl<'a> FastqTokenizer<'a> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn reset(&mut self, s: &'a str) {
        self.ptr = 0;
        self.end = s.len();
        self.s = s;
    }

    pub fn current(&self) -> String {
        match self.s.as_bytes()[self.ptr] {
            b'@' => {
                let i = self.find_byte(self.ptr, b'\n');
                trim_cr(&self.s[self.ptr + 1..], i - self.ptr - 1).to_string()
            }
            b'+' => {
                let i = self.find_byte(self.ptr, b'\n');
                if i == self.end {
                    String::new()
                } else {
                    trim_cr(&self.s[i + 1..], self.end - i - 1).to_string()
                }
            }
            _ => {
                let i = self.find_byte(self.ptr, b'\n');
                trim_cr(&self.s[self.ptr..], i - self.ptr).to_string()
            }
        }
    }

    pub fn increment(&mut self) -> Result<&mut Self, MalformedFastqRecord> {
        match self.s.as_bytes()[self.ptr] {
            b'@' => self.ptr = (self.find_byte(self.ptr, b'\n') + 1).min(self.end),
            b'+' => self.ptr = self.end,
            _ => {
                self.ptr = (self.find_byte(self.ptr, b'\n') + 1).min(self.end);
                if self.ptr < self.end && self.s.as_bytes()[self.ptr] != b'+' {
                    return Err(MalformedFastqRecord::new(-1));
                }
            }
        }
        Ok(self)
    }

    pub fn read_record(&mut self) -> Result<(bool, String), MalformedFastqRecord> {
        if self.ptr == self.end {
            return Ok((false, String::new()));
        }
        if self.s.as_bytes()[self.ptr] != b'@' {
            return Err(MalformedFastqRecord::new(-1));
        }
        let begin = self.ptr;
        let (mut l, mut n) = find_newline(&self.s[self.ptr..self.end]);
        if n == 0 {
            return Err(MalformedFastqRecord::new(-1));
        }
        l += self.ptr;
        self.ptr = l;
        let mut len = 0i64;
        loop {
            let r = find_newline(&self.s[self.ptr..self.end]);
            l = self.ptr + r.0;
            n = r.1;
            if n == 0 {
                return Err(MalformedFastqRecord::new(-1));
            }
            len += (l - self.ptr - n as usize) as i64;
            self.ptr = l;
            if self.ptr == self.end {
                return Err(MalformedFastqRecord::new(-1));
            }
            if self.s.as_bytes()[self.ptr] == b'+' {
                let r = find_newline(&self.s[self.ptr..self.end]);
                l = self.ptr + r.0;
                n = r.1;
                if n == 0 {
                    return Err(MalformedFastqRecord::new(-1));
                }
                self.ptr = l;
                break;
            }
        }
        let mut qlen = 0i64;
        loop {
            let r = find_newline(&self.s[self.ptr..self.end]);
            l = self.ptr + r.0;
            n = r.1;
            qlen += (l - self.ptr - n as usize) as i64;
            self.ptr = l;
            if self.ptr == self.end || qlen >= len {
                break;
            }
        }
        Ok((true, self.s[begin..self.ptr].to_string()))
    }

    fn find_byte(&self, p: usize, c: u8) -> usize {
        self.s.as_bytes()[p..self.end]
            .iter()
            .position(|&b| b == c)
            .map(|i| p + i)
            .unwrap_or(self.end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_suffix_helpers() {
        assert!(ends_with("abc.dmnd", ".dmnd"));
        assert!(!ends_with("abc", ".dmnd"));
        assert_eq!(rstrip("abc.dmnd", ".dmnd"), "abc");
        assert_eq!(rstrip("abc", ".dmnd"), "abc");
        assert_eq!(max_len(&["a", "abcd", "xy"]), 4);
    }

    #[test]
    fn test_convert_size_join_and_format_double() {
        assert_eq!(convert_size(1), "1.0 B");
        assert_eq!(convert_size(1536), "1.5 KB");
        assert_eq!(join(",", &[1, 2, 3]), "1,2,3");
        assert_eq!(format_double(99.94), "99.9");
        assert_eq!(format_double(99.95), "100.0");
        assert_eq!(format_double(100.9), "100");
    }

    #[test]
    fn test_replace_ratio_interpret_number() {
        assert_eq!(replace("a-b-c", '-', '_'), "a_b_c");
        assert_eq!(ratio_percentage_usize(1, 4), "1/4 (25.00%)");
        assert_eq!(ratio_percentage(1usize, 4usize), "1/4 (25.00%)");
        assert_eq!(ratio_percentage(1.5f64, 3.0f64), "2/3 (50.00%)");
        assert_eq!(interpret_number("1.5K").unwrap(), 1500);
        assert_eq!(interpret_number("2M").unwrap(), 2_000_000);
        assert_eq!(interpret_number(" 2M \t").unwrap(), 2_000_000);
        assert!(interpret_number("10").is_err());
        assert!(interpret_number("10X").is_err());
    }

    #[test]
    fn test_tokenize_parse_csv_convert_string() {
        assert_eq!(tokenize("a,b;;c", ",;"), vec!["a", "b", "c"]);
        assert_eq!(tokenize(",,,", ","), vec![""]);
        assert_eq!(
            parse_csv("1,2,,2").into_iter().collect::<Vec<_>>(),
            vec![1, 2]
        );
        assert_eq!(convert_string_i64("-7").unwrap(), -7);
        assert_eq!(convert_string_i64(" +7").unwrap(), 7);
        assert!(convert_string_i64("7 ").is_err());
        assert_eq!(convert_string_i32("2147483647").unwrap(), i32::MAX);
        assert!(convert_string_i32("2147483648").is_err());
        assert_eq!(convert_string_u64("7").unwrap(), 7);
        assert_eq!(convert_string_u64(" -2").unwrap(), u64::MAX - 1);
        assert!(convert_string_u32("-1").is_err());
        assert_eq!(convert_string::<i64>("-8").unwrap(), -8);
        assert_eq!(convert_string::<u32>("8").unwrap(), 8);
    }

    #[test]
    fn test_format_number() {
        assert_eq!(format_number(0.0), "0");
        assert_eq!(format_number(999.0), "999");
        assert_eq!(format_number(1000.0), "1K");
        assert_eq!(format_number(1_234_567.0), "1.23M");
        assert_eq!(format_number(-1500.0), "-1.5K");
        assert_eq!(format_number(999_995.0), "1M");
        assert_eq!(format(1_234_567u64), "1.23M");
        assert_eq!(format(1_234_567u32), "1.23M");
        assert_eq!(format(-1500i64), "-1.5K");
        assert_eq!(format(-1500i32), "-1.5K");
    }

    #[test]
    fn test_dyn_tokenizer_helpers() {
        assert_eq!(trim_cr("abc\r\n", 4), "abc");
        assert_eq!(trim_cr("abc\n", 3), "abc");
        assert_eq!(find_newline("abc\nrest"), (4, 1));
        assert_eq!(find_newline("abc\r\nrest"), (5, 2));
        assert_eq!(find_newline("abc"), (3, 0));
        assert!(is_whitespace(" \r\n"));
        assert!(!is_whitespace(" \tx"));
    }

    #[test]
    fn test_delimiters_and_tokenizer() {
        let d = CharDelimiter::new(',');
        assert_eq!(d.scan("a,b", 0), Some((1, 2)));
        assert_eq!(d.next("a,", 2).unwrap(), None);
        assert!(d.next("a,b", 1).is_ok());
        assert!(d.next("a,b", 0).is_err());

        let d = StringDelimiter::new("::");
        assert_eq!(d.scan("a::b", 0), Some((1, 3)));
        assert_eq!(d.next("a::b", 1).unwrap(), Some(3));

        let d = StringDelimiters::new(&["::", "|"]);
        assert_eq!(d.scan("a|b::c", 0), Some((3, 5)));
        assert_eq!(d.next("::x", 0).unwrap(), Some(2));

        let mut t = Tokenizer::new("12,34,tail", CharDelimiter::new(','));
        assert_eq!(t.read_i32().unwrap(), 12);
        assert_eq!(t.read_u64().unwrap(), 34);
        assert_eq!(t.read_string().unwrap(), "tail");
        assert!(!t.good());

        let mut t = Tokenizer::new("x=2.5", CharDelimiter::new('='));
        assert_eq!(t.read_skip(Skip).unwrap().ptr(), Some(2));
        assert_eq!(t.read_f32().unwrap(), 2.5);
        assert_eq!(t.getline(), "");

        let mut t = Tokenizer::new("12abc,7", CharDelimiter::new(','));
        assert!(t.read_i64().is_err());
        assert_eq!(t.ptr(), Some(0));

        let mut t = Tokenizer::new("12,7", CharDelimiter::new(','));
        assert_eq!(t.read_i64().unwrap(), 12);
        assert_eq!(t.ptr(), Some(3));

        let mut t = Tokenizer::new("1.5e2,tail", CharDelimiter::new(','));
        assert_eq!(t.read_f64().unwrap(), 150.0);
        assert_eq!(t.read_string().unwrap(), "tail");

        let mut t = Tokenizer::new("1e,tail", CharDelimiter::new(','));
        assert!(t.read_f64().is_err());
        assert_eq!(t.ptr(), Some(0));

        let mut t = Tokenizer::new("abc\ndef", CharDelimiter::new(','));
        assert_eq!(t.getline(), "abc");
        t.skip_to('\n');
        assert_eq!(t.getline(), "def");
    }

    #[test]
    fn test_char_and_token_iterators() {
        let mut t = CharTokenizer::new('\t');
        t.reset("a\tb\nc");
        assert!(t.good());
        assert_eq!(t.current(), "a");
        t.increment();
        assert_eq!(t.current(), "b\nc");

        let mut lines = LineIterator::new("a\nb\n");
        assert_eq!(lines.current(), "a");
        lines.increment();
        assert_eq!(lines.current(), "b");
        lines.increment();
        assert!(!lines.good());

        let mut tab = TabIterator::new("x\ty");
        assert_eq!(tab.current(), "x");
        tab.increment();
        assert_eq!(tab.ptr(), 2);
        assert_eq!(tab.current(), "y");

        let mut r = CharTokenizer::new(',');
        r.reset("abc\r\ndef");
        assert_eq!(r.read_record(), (true, "abc".to_string()));
        assert_eq!(r.read_record(), (true, "def".to_string()));
        assert_eq!(r.read_record(), (false, String::new()));
    }

    #[test]
    fn test_fasta_tokenizer() {
        let mut t = FastaTokenizer::new();
        t.reset(">id\r\nAC\nGT\r\n");
        assert_eq!(t.current(), "id");
        t.increment();
        assert_eq!(t.current(), "ACGT");

        t.reset(">a\nAA\n>b\nBB\n");
        assert_eq!(t.read_record(), (true, ">a\nAA".to_string()));
        assert_eq!(t.read_record(), (true, String::new()));
    }

    #[test]
    fn test_fastq_tokenizer() {
        let mut t = FastqTokenizer::new();
        t.reset("@id\r\nAC\n+\n!!\n");
        assert_eq!(t.current(), "id");
        t.increment().unwrap();
        assert_eq!(t.current(), "AC");
        t.increment().unwrap();
        assert_eq!(t.current(), "!!\n");

        t.reset("@id\nAC\n+\n!!\n");
        assert_eq!(
            t.read_record().unwrap(),
            (true, "@id\nAC\n+\n!!\n".to_string())
        );
        assert_eq!(t.read_record().unwrap(), (false, String::new()));

        t.reset("bad\n");
        assert!(t.read_record().is_err());
        assert_eq!(
            MalformedFastqRecord::new(7).to_string(),
            "Malformed FASTQ record on line 7"
        );
    }
}
