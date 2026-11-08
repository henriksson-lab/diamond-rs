//! Expression parser used by the MCL clustering output format.
//!
//! This mirrors the implementation embedded in
//! `diamond/src/contrib/mcl/recursive_parser.h`.

use super::clustering_variables::{UnknownVariableError, Variable, VariableRegistry};
use crate::align::hsp::HspContext;
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseError {
    message: crate::util::io::RawDiagnosticError,
}

impl ParseError {
    pub(crate) fn new(message: impl Into<String>) -> Self {
        Self {
            message: crate::util::io::RawDiagnosticError(message.into().into_bytes()),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.message, f)
    }
}

impl std::error::Error for ParseError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.message)
    }
}
impl ParseError {
    pub fn from_bytes(message: Vec<u8>) -> Self {
        Self {
            message: crate::util::io::RawDiagnosticError(message),
        }
    }
    pub fn into_io(self) -> std::io::Error {
        std::io::Error::other(self)
    }
}

impl From<UnknownVariableError> for ParseError {
    fn from(value: UnknownVariableError) -> Self {
        Self::new(value.to_string())
    }
}

pub struct RecursiveParser<'a, T = String> {
    context: Option<&'a HspContext<T>>,
    expression: &'a [u8],
    position: usize,
    variables: Vec<&'static Variable>,
    query_translated: bool,
}

impl<'a> RecursiveParser<'a> {
    pub fn clean_expression(expression: &str) -> String {
        String::from_utf8(Self::clean_expression_bytes(expression.as_bytes()))
            .expect("ASCII whitespace removal preserves UTF-8")
    }
    pub fn clean_expression_bytes(expression: &[u8]) -> Vec<u8> {
        expression
            .iter()
            .copied()
            .filter(|byte| !byte.is_ascii_whitespace())
            .collect()
    }
    pub fn new_bytes(context: Option<&'a HspContext>, expression: &'a [u8]) -> Self {
        Self::new_with_context_mode_bytes(context, expression, false)
    }
    pub fn new(context: Option<&'a HspContext>, expression: &'a str) -> Self {
        Self::new_with_context_mode(context, expression, false)
    }
    pub fn new_with_mode(
        context: Option<&'a HspContext>,
        expression: &'a str,
        translated: bool,
    ) -> Self {
        Self::new_with_context_mode(context, expression, translated)
    }
}

impl<'a, T> RecursiveParser<'a, T> {
    pub fn new_with_context(context: Option<&'a HspContext<T>>, expression: &'a str) -> Self {
        Self::new_with_context_mode(context, expression, false)
    }

    pub fn new_with_context_mode(
        context: Option<&'a HspContext<T>>,
        expression: &'a str,
        query_translated: bool,
    ) -> Self {
        Self::new_with_context_mode_bytes(context, expression.as_bytes(), query_translated)
    }
    pub fn new_with_context_mode_bytes(
        context: Option<&'a HspContext<T>>,
        expression: &'a [u8],
        query_translated: bool,
    ) -> Self {
        // Source parser receives c_str(), not the full stored std::string.
        let expression = crate::util::sequence::c_string_bytes(expression);
        Self {
            context,
            expression,
            position: 0,
            variables: Vec::new(),
            query_translated,
        }
    }

    fn peek(&self, ahead: usize) -> Option<u8> {
        self.expression.get(self.position + ahead).copied()
    }

    fn get(&mut self) -> Result<u8, ParseError> {
        let value = self
            .peek(0)
            .ok_or_else(|| ParseError::new("Unexpected end of clustering expression"))?;
        self.position += 1;
        Ok(value)
    }

    fn consume(&mut self, _expected: u8) -> Result<(), ParseError> {
        // Source grammar simply advances over these separator positions.
        // Retain a bounds error only at its genuine out-of-allocation edge.
        if self.position >= self.expression.len() {
            return Err(ParseError::new("Unexpected end of clustering expression"));
        }
        self.position += 1;
        Ok(())
    }

    fn starts_with(&self, text: &[u8]) -> bool {
        self.expression[self.position..].starts_with(text)
    }

    fn relation(&mut self) -> Result<fn(f64, f64) -> bool, ParseError> {
        let relations: &[(&[u8], fn(f64, f64) -> bool)] = &[
            (b">=", |a, b| a >= b),
            (b"<=", |a, b| a <= b),
            (b"==", |a, b| a == b),
            (b"=", |a, b| a == b),
            (b">", |a, b| a > b),
            (b"<", |a, b| a < b),
        ];
        for &(token, operation) in relations {
            if self.starts_with(token) {
                self.position += token.len();
                return Ok(operation);
            }
        }
        // The Source constructs a runtime_error here but never throws it.
        // Its fallback predicate is false and leaves the cursor untouched.
        Ok(|_, _| false)
    }

    fn variable_name(&mut self) -> Result<&'a str, ParseError> {
        let start = self.position;
        while matches!(self.peek(0), Some(b'a'..=b'z' | b'A'..=b'Z' | b'_')) {
            self.position += 1;
        }
        // The variable lexer admits ASCII letters/underscore only, including
        // an empty key for an unrecognized leading byte (Source registry error).
        Ok(std::str::from_utf8(&self.expression[start..self.position])
            .expect("ASCII variable lexer"))
    }

    fn integer(&mut self) -> Result<u32, ParseError> {
        let first = self.get()?;
        // Source integer() consumes the first byte unconditionally. This
        // also governs the fractional substring following a decimal point.
        let mut result = (i32::from(first as libc::c_char) - i32::from(b'0')) as u32;
        while let Some(digit @ b'0'..=b'9') = self.peek(0) {
            self.position += 1;
            result = result
                .wrapping_mul(10)
                .wrapping_add(u32::from(digit - b'0'));
        }
        Ok(result)
    }

    fn number(&mut self) -> Result<f64, ParseError> {
        let mut result = self.integer()? as f64;
        if self.peek(0) == Some(b'.') {
            self.position += 1;
            let start = self.position;
            let decimals = self.integer()?;
            // Preserve the source implementation literally: its divisor is
            // the number of decimal characters, not a power of ten.
            result += decimals as f64 / (self.position - start) as f64;
        }
        Ok(result)
    }

    fn factor(&mut self) -> Result<f64, ParseError> {
        match self.peek(0) {
            Some(b'0'..=b'9') => self.number(),
            Some(b'(') => {
                self.position += 1;
                let result = self.expression()?;
                self.consume(b')')?;
                Ok(result)
            }
            Some(b'-') => {
                self.position += 1;
                Ok(-self.factor()?)
            }
            _ if self.starts_with(b"max(") => {
                self.position += 4;
                let first = self.expression()?;
                self.consume(b',')?;
                let second = self.expression()?;
                self.consume(b')')?;
                // std::max returns its first operand on an unordered or
                // equal comparison; f64::max instead discards a NaN.
                Ok(if first < second { second } else { first })
            }
            _ if self.starts_with(b"min(") => {
                self.position += 4;
                let first = self.expression()?;
                self.consume(b',')?;
                let second = self.expression()?;
                self.consume(b')')?;
                Ok(if second < first { second } else { first })
            }
            _ if self.starts_with(b"exp(") => {
                self.position += 4;
                let result = self.expression()?.exp();
                self.consume(b')')?;
                Ok(result)
            }
            _ if self.starts_with(b"log(") => {
                self.position += 4;
                let result = self.expression()?.ln();
                self.consume(b')')?;
                Ok(result)
            }
            _ if self.starts_with(b"I(") => {
                self.position += 2;
                let first = self.expression()?;
                let relation = self.relation()?;
                let second = self.expression()?;
                self.consume(b')')?;
                Ok(u8::from(relation(first, second)) as f64)
            }
            Some(_) => {
                let name = self.variable_name()?;
                let variable = VariableRegistry::get(name)?;
                if let Some(context) = self.context {
                    Ok(variable.get(context, self.query_translated))
                } else {
                    self.variables.push(variable);
                    // C++ uses a harmless dummy value while discovering vars.
                    Ok(4.0)
                }
            }
            None => Err(ParseError::new("Unknown variable: ")),
        }
    }

    fn term(&mut self) -> Result<f64, ParseError> {
        let mut result = self.factor()?;
        while matches!(self.peek(0), Some(b'*' | b'/')) {
            if self.get()? == b'*' {
                result *= self.factor()?;
            } else {
                result /= self.factor()?;
            }
        }
        Ok(result)
    }

    fn expression(&mut self) -> Result<f64, ParseError> {
        let mut result = self.term()?;
        while matches!(self.peek(0), Some(b'+' | b'-')) {
            if self.get()? == b'+' {
                result += self.term()?;
            } else {
                result -= self.term()?;
            }
        }
        Ok(result)
    }

    pub fn evaluate(&mut self) -> Result<f64, ParseError> {
        self.expression()
    }

    pub fn variables(&self) -> &[&'static Variable] {
        &self.variables
    }
}

/// Actual Source syntax-only evaluation (the caller chooses the run phase).
pub fn evaluate_bytes(expression: &[u8]) -> std::io::Result<f64> {
    let cleaned = RecursiveParser::clean_expression_bytes(expression);
    RecursiveParser::new_bytes(None, &cleaned)
        .evaluate()
        .map_err(ParseError::into_io)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;

    fn context() -> HspContext {
        let mut hsp = Hsp::default();
        hsp.bit_score = 50.0;
        hsp.length = 10;
        hsp.identities = 8;
        HspContext::new(
            hsp,
            1,
            2,
            vec![vec![0; 20]],
            20,
            "q",
            3,
            40,
            "s",
            0,
            0,
            vec![0; 40],
            100.0,
            80.0,
        )
    }

    #[test]
    fn byte_expression_preserves_source_lexer_and_separator_behavior() {
        assert_eq!(evaluate_bytes(b" max(1;2] ").unwrap(), 2.0);
        assert_eq!(evaluate_bytes(b"1\0+invalid").unwrap(), 1.0);
        assert_eq!(evaluate_bytes(b"4294967296").unwrap(), 0.0);
        assert_eq!(evaluate_bytes(b"1.a").unwrap(), 50.0);
        let error = evaluate_bytes(b"\xff").unwrap_err();
        assert_eq!(
            crate::util::io::RawDiagnosticError::find(&error).unwrap(),
            b"Unknown variable: "
        );
        assert_eq!(
            evaluate_bytes(b"I(1!2)").unwrap_err().to_string(),
            "Unknown variable: "
        );
    }

    #[test]
    fn evaluates_source_grammar_and_relations() {
        let context = context();
        let mut parser =
            RecursiveParser::new(Some(&context), "max(bitscore/2,min(30,qlen))+I(pident>=80)");
        assert_eq!(parser.evaluate().unwrap(), 26.0);
        let mut parser = RecursiveParser::new(None, "exp(log(4))-(-2)");
        assert!((parser.evaluate().unwrap() - 6.0).abs() < 1e-12);
    }

    #[test]
    fn discovery_retains_duplicates_and_source_decimal_quirk() {
        let mut parser = RecursiveParser::new(None, "bitscore+bitscore+1.25");
        assert_eq!(parser.evaluate().unwrap(), 21.5);
        assert_eq!(
            parser
                .variables()
                .iter()
                .map(|variable| variable.get_name())
                .collect::<Vec<_>>(),
            ["bitscore", "bitscore"]
        );
    }

    #[test]
    fn cleaning_and_errors_are_explicit() {
        assert_eq!(RecursiveParser::clean_expression(" a +\n b\t"), "a+b");
        assert!(RecursiveParser::new(None, "unknown").evaluate().is_err());
        assert!(RecursiveParser::new(None, "I(1!2)").evaluate().is_err());
    }
}
