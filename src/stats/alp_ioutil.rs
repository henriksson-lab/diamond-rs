#![allow(non_snake_case)]

use std::sync::Mutex;

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
    panic!("IoUtil::abort")
}

pub fn abort_msg(s_: &str) -> ! {
    panic!("{}", s_)
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
            .find(|&(_, c)| !c.is_whitespace())
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
            .find(|&(_, c)| !c.is_whitespace())
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
    match str_.parse::<T>() {
        Ok(value) => {
            *value_ = value;
            true
        }
        Err(_) => false,
    }
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
        match s.parse::<f64>() {
            Ok(x) => {
                *x_ = x;
                true
            }
            Err(_) => false,
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
    }
}
