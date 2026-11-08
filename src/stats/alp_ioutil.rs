#![allow(non_snake_case)]

use std::sync::Mutex;

/// Default C++ stream insertion for the primitive numeric values used by
/// ALP's Matrix/Vector templates. Custom value types retain their Display
/// implementation, corresponding to their own stream insertion operator.
pub(crate) fn format_value<T: std::fmt::Display + std::any::Any>(value: &T) -> String {
    // Inspect the primitive itself before Display can discard a NaN sign.
    // Any proves the concrete type safely; custom owned values keep their
    // own Display insertion, rather than being guessed from rendered text.
    let primitive = value as &dyn std::any::Any;
    let number = if let Some(value) = primitive.downcast_ref::<f32>() {
        Some((
            f64::from(*value),
            value.is_nan() && value.is_sign_negative(),
        ))
    } else {
        primitive
            .downcast_ref::<f64>()
            .map(|value| (*value, value.is_nan() && value.is_sign_negative()))
    };
    match number {
        Some((_, true)) => "-nan".to_owned(),
        Some((value, false)) => crate::output::format::format_default_float(value),
        None => value.to_string(),
    }
}

/// Formatted extraction for ALP's primitive values. Numeric extraction consumes
/// a prefix and leaves the next character in the stream, unlike token.parse().
pub(crate) fn read_value<T: std::str::FromStr>(input: &str, offset: &mut usize) -> Option<T> {
    let bytes = input.as_bytes();
    while *offset < bytes.len() && bytes[*offset].is_ascii_whitespace() {
        *offset += 1;
    }
    let start = *offset;
    if start == bytes.len() {
        return None;
    }
    let kind = std::any::type_name::<T>();
    let float = matches!(kind, "f32" | "f64");
    let integer = matches!(
        kind,
        "i16" | "i32" | "i64" | "i128" | "isize" | "u16" | "u32" | "u64" | "u128" | "usize"
    );
    if float || integer {
        if matches!(bytes.get(*offset), Some(b'+' | b'-')) {
            *offset += 1;
        }
        let digits = *offset;
        while bytes.get(*offset).is_some_and(u8::is_ascii_digit) {
            *offset += 1;
        }
        let mut has_digits = *offset > digits;
        if float && bytes.get(*offset) == Some(&b'.') {
            *offset += 1;
            let fraction = *offset;
            while bytes.get(*offset).is_some_and(u8::is_ascii_digit) {
                *offset += 1;
            }
            has_digits |= *offset > fraction;
        }
        if !has_digits {
            return None;
        }
        if float && matches!(bytes.get(*offset), Some(b'e' | b'E')) {
            *offset += 1;
            if matches!(bytes.get(*offset), Some(b'+' | b'-')) {
                *offset += 1;
            }
            let exponent = *offset;
            while bytes.get(*offset).is_some_and(u8::is_ascii_digit) {
                *offset += 1;
            }
            if *offset == exponent {
                return None;
            }
        }
        let token = &input[start..*offset];
        if float {
            let finite = if kind == "f32" {
                token.parse::<f32>().ok()?.is_finite()
            } else {
                token.parse::<f64>().ok()?.is_finite()
            };
            return finite.then(|| token.parse().ok()).flatten();
        }
        if kind.starts_with('u') && token.starts_with('-') {
            let bits = match kind {
                "u16" => 16,
                "u32" => 32,
                "u64" => 64,
                "usize" => usize::BITS,
                _ => 128,
            };
            let magnitude = token[1..].parse::<u128>().ok()?;
            let max = u128::MAX >> (128 - bits);
            if magnitude > max {
                return None;
            }
            return (0u128.wrapping_sub(magnitude) & max)
                .to_string()
                .parse()
                .ok();
        }
        token.parse().ok()
    } else {
        while *offset < bytes.len() && !bytes[*offset].is_ascii_whitespace() {
            *offset += 1;
        }
        input[start..*offset].parse().ok()
    }
}

pub(crate) fn read_values<T: std::str::FromStr>(input: &str) -> (Vec<T>, bool) {
    let mut offset = 0;
    let mut values = Vec::new();
    while let Some(value) = read_value(input, &mut offset) {
        values.push(value);
    }
    (values, offset == input.len())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    HUMAN,
    MACHINE,
    FORMAT,
}

const FORMAT_DEFAULT: Format = Format::HUMAN;
const TERMINATOR_DEFAULT: char = '!';

static FORMAT_STATE: Mutex<Format> = Mutex::new(FORMAT_DEFAULT);
static TERMINATOR_STATE: Mutex<char> = Mutex::new(TERMINATOR_DEFAULT);

pub fn getFormat() -> Format {
    *FORMAT_STATE.lock().unwrap()
}

pub fn setFormat(format_: Format) {
    *FORMAT_STATE.lock().unwrap() = format_;
}

pub fn clearFormat() -> Format {
    let mut format = FORMAT_STATE.lock().unwrap();
    *format = FORMAT_DEFAULT;
    *format
}

pub fn getTerminator() -> char {
    *TERMINATOR_STATE.lock().unwrap()
}

pub fn setTerminator(t_: char) {
    *TERMINATOR_STATE.lock().unwrap() = t_;
}

pub fn clearTerminator() -> char {
    let mut terminator = TERMINATOR_STATE.lock().unwrap();
    *terminator = TERMINATOR_DEFAULT;
    *terminator
}

pub fn abort() -> ! {
    std::process::exit(1)
}

pub fn abort_msg(s_: &str) -> ! {
    eprintln!("{s_}");
    abort()
}

pub fn getLine_lines<I>(in_: &mut I, str_: &mut String, t_: char) -> bool
where
    I: Iterator<Item = String>,
{
    assert!(t_ != '\0');

    str_.clear();
    for mut line in in_ {
        let first_non_ws = line
            .char_indices()
            .find(|&(_, c)| !c.is_ascii_whitespace())
            .map(|(_, c)| c);
        if first_non_ws.is_none() || first_non_ws == Some(t_) {
            continue;
        }

        if t_ != '\n' {
            if let Some(pos) = line.find(t_) {
                line.truncate(pos);
            }
        }
        *str_ = line;
        return true;
    }

    false
}

pub fn getLine(input: &str, offset: &mut usize, str_: &mut String, t_: char) -> bool {
    assert!(t_ != '\0');

    str_.clear();
    while *offset <= input.len() {
        if *offset == input.len() {
            return false;
        }

        let rest = &input[*offset..];
        let line_end = rest.find('\n').map(|p| *offset + p).unwrap_or(input.len());
        let mut line = input[*offset..line_end].to_string();
        *offset = if line_end < input.len() {
            line_end + 1
        } else {
            line_end
        };

        let first_non_ws = line
            .char_indices()
            .find(|&(_, c)| !c.is_ascii_whitespace())
            .map(|(_, c)| c);
        if first_non_ws.is_none() || first_non_ws == Some(t_) {
            continue;
        }

        if t_ != '\n' {
            if let Some(pos) = line.find(t_) {
                line.truncate(pos);
            }
        }
        *str_ = line;
        return true;
    }

    false
}

pub fn getString<T>(
    input: &str,
    offset: &mut usize,
    value_: &mut T,
    line_: &mut String,
    t_: char,
) -> bool
where
    T: std::str::FromStr,
{
    if !getLine(input, offset, line_, t_) {
        return false;
    }
    let Some(str_) = line_.split_whitespace().next() else {
        return false;
    };
    if let Ok(value) = str_.parse::<T>() {
        *value_ = value;
    }
    // The template parses through a separate stringstream and returns the
    // state of the original input stream, so conversion failure is not
    // propagated to the caller.
    true
}

pub fn getLine_stream(input: &str, offset: &mut usize, t_: char) -> Option<String> {
    let mut s = String::new();
    if getLine(input, offset, &mut s, t_) {
        Some(s)
    } else {
        None
    }
}

pub fn in_double(token: &str, x_: &mut f64) -> bool {
    let s = token.to_ascii_lowercase();
    if s == "1.#inf" || s == "nan" {
        *x_ = f64::INFINITY;
        true
    } else {
        let parsed = s
            .char_indices()
            .map(|(index, _)| index)
            .chain(std::iter::once(s.len()))
            .rev()
            .find_map(|end| s[..end].parse::<f64>().ok());
        match parsed {
            Some(x) => {
                *x_ = x;
                true
            }
            None => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static TEST_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn test_format_and_terminator_state() {
        let _guard = TEST_LOCK.lock().unwrap();
        assert_eq!(clearFormat(), Format::HUMAN);
        assert_eq!(getFormat(), Format::HUMAN);
        setFormat(Format::MACHINE);
        assert_eq!(getFormat(), Format::MACHINE);
        assert_eq!(clearFormat(), Format::HUMAN);

        assert_eq!(clearTerminator(), '!');
        setTerminator('#');
        assert_eq!(getTerminator(), '#');
        assert_eq!(clearTerminator(), '!');
    }

    #[test]
    fn test_get_line_skips_blank_and_comment_parts() {
        let mut offset = 0;
        let input = "   \n! comment\n  42 24 ! trailing\nnext\n";
        let mut line = String::new();
        assert!(getLine(input, &mut offset, &mut line, '!'));
        assert_eq!(line, "  42 24 ");
        assert!(getLine(input, &mut offset, &mut line, '!'));
        assert_eq!(line, "next");
        assert!(!getLine(input, &mut offset, &mut line, '!'));
    }

    #[test]
    fn test_get_string_and_double_input() {
        let mut offset = 0;
        let mut line = String::new();
        let mut value = 0i32;
        assert!(getString(
            "  17 rest\n",
            &mut offset,
            &mut value,
            &mut line,
            '!'
        ));
        assert_eq!(value, 17);
        assert_eq!(line, "  17 rest");

        let mut x = 0.0;
        assert!(in_double("1.#INF", &mut x));
        assert!(x.is_infinite());
        assert!(in_double("nan", &mut x));
        assert!(x.is_infinite());
        assert!(in_double("3.5", &mut x));
        assert_eq!(x, 3.5);
        assert!(!in_double("abc", &mut x));
        assert!(in_double("3.5trailing", &mut x));
        assert_eq!(x, 3.5);
    }

    #[test]
    fn get_string_does_not_propagate_secondary_conversion_failure() {
        let mut offset = 0;
        let mut line = String::new();
        let mut value = 41_i32;
        assert!(getString(
            "not-an-int rest\n",
            &mut offset,
            &mut value,
            &mut line,
            '!'
        ));
        assert_eq!(value, 41);
        assert_eq!(line, "not-an-int rest");
    }
}
