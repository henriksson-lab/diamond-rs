//! Hierarchy-preserving facade for `diamond/src/lib/alp/njn_ioutil.cpp`.
//!
//! The implementation predates the mirrored `lib/alp` tree and remains in
//! `stats::alp_ioutil` for compatibility with its matrix/vector consumers.

pub use crate::stats::alp_ioutil::Format;

pub const FORMAT: Format = Format::HUMAN;
pub const TERMINATOR: char = '!';

pub fn get_format() -> Format {
    crate::stats::alp_ioutil::getFormat()
}
pub fn set_format(format: Format) {
    crate::stats::alp_ioutil::setFormat(format)
}
pub fn clear_format() -> Format {
    crate::stats::alp_ioutil::clearFormat()
}
pub fn get_terminator() -> char {
    crate::stats::alp_ioutil::getTerminator()
}
pub fn set_terminator(value: char) {
    crate::stats::alp_ioutil::setTerminator(value)
}
pub fn clear_terminator() -> char {
    crate::stats::alp_ioutil::clearTerminator()
}
pub fn abort() -> ! {
    crate::stats::alp_ioutil::abort()
}
pub fn abort_message(message: &str) -> ! {
    crate::stats::alp_ioutil::abort_msg(message)
}
pub fn get_line(input: &str, offset: &mut usize, line: &mut String, terminator: char) -> bool {
    crate::stats::alp_ioutil::getLine(input, offset, line, terminator)
}
pub fn get_line_stream(input: &str, offset: &mut usize, terminator: char) -> Option<String> {
    crate::stats::alp_ioutil::getLine_stream(input, offset, terminator)
}
pub fn get_string<T: std::str::FromStr>(
    input: &str,
    offset: &mut usize,
    value: &mut T,
    line: &mut String,
    terminator: char,
) -> bool {
    crate::stats::alp_ioutil::getString(input, offset, value, line, terminator)
}
pub fn input_double(token: &str, value: &mut f64) -> bool {
    crate::stats::alp_ioutil::in_double(token, value)
}

/// Header `operator<<`/`operator>>` both only update the shared format state.
pub fn apply_format(format: Format) {
    set_format(format);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn facade_constants_and_format_manipulator_match_header() {
        assert_eq!(FORMAT, Format::HUMAN);
        assert_eq!(TERMINATOR, '!');
        apply_format(Format::MACHINE);
        assert_eq!(get_format(), Format::MACHINE);
        assert_eq!(clear_format(), FORMAT);
    }
}
