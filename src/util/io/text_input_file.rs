//! Translation of `diamond/src/util/io/text_input_file.{h,cpp}`.

use super::{
    Compressor, Deserializer, InputFile, InputFileStream, IoError, IoResult, OutputFile,
    StreamEntity, TempFile,
};

#[derive(Debug, Clone)]
pub struct TextInputFile<S: StreamEntity> {
    input: Deserializer<S>,
    pub line: String,
    pub line_count: usize,
    putback_line: bool,
    eof: bool,
    line_separator: u8,
}

impl<S: StreamEntity> TextInputFile<S> {
    pub const DEFAULT_LINE_SEPARATOR: u8 = b'\n';

    /// Generic compatibility constructor. The three C++ constructors are
    /// represented by `from_file_name`, `from_temp_file`, and
    /// `from_output_file` below.
    pub fn new(input: Deserializer<S>, line_separator: u8) -> Self {
        Self {
            input,
            line: String::new(),
            line_count: 0,
            putback_line: false,
            eof: false,
            line_separator,
        }
    }

    /// Rust spelling of the C++ constructor's default separator argument.
    pub fn with_default_separator(input: Deserializer<S>) -> Self {
        Self::new(input, Self::DEFAULT_LINE_SEPARATOR)
    }

    pub fn rewind(&mut self) -> IoResult<()> {
        self.input.rewind()?;
        self.line_count = 0;
        self.putback_line = false;
        self.eof = false;
        self.line.clear();
        Ok(())
    }

    pub fn eof(&self) -> bool {
        self.eof
    }

    pub fn getline(&mut self) -> IoResult<()> {
        if self.putback_line {
            self.putback_line = false;
            self.line_count = self.line_count.wrapping_add(1);
            return Ok(());
        }

        let mut line = Vec::new();
        self.eof = !self.input.read_to(&mut line, self.line_separator)?;
        self.line_count = self.line_count.wrapping_add(1);
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        self.line = String::from_utf8(line).map_err(|error| IoError::Other(error.to_string()))?;
        Ok(())
    }

    pub fn putback_line(&mut self) {
        self.putback_line = true;
        // C++ line_count is size_t, so its unchecked `--line_count` wraps.
        self.line_count = self.line_count.wrapping_sub(1);
    }

    /// C++ `operator bool()`.
    pub fn is_open(&self) -> bool {
        !self.eof()
    }

    pub fn into_inner(self) -> Deserializer<S> {
        self.input
    }
}

impl TextInputFile<InputFileStream> {
    pub fn from_file_name(file_name: &str, line_separator: u8) -> IoResult<Self> {
        let input = InputFile::new(file_name, 0)?;
        Ok(Self::new(input.into_deserializer(), line_separator))
    }

    pub fn from_file_name_default(file_name: &str) -> IoResult<Self> {
        Self::from_file_name(file_name, Self::DEFAULT_LINE_SEPARATOR)
    }

    pub fn from_temp_file(tmp_file: &mut TempFile, line_separator: u8) -> IoResult<Self> {
        let input = InputFile::from_temp_file(tmp_file, 0, Compressor::None)?;
        Ok(Self::new(input.into_deserializer(), line_separator))
    }

    pub fn from_temp_file_default(tmp_file: &mut TempFile) -> IoResult<Self> {
        Self::from_temp_file(tmp_file, Self::DEFAULT_LINE_SEPARATOR)
    }

    pub fn from_output_file(out_file: &mut OutputFile, line_separator: u8) -> IoResult<Self> {
        let input = InputFile::from_output_file(out_file, 0)?;
        Ok(Self::new(input.into_deserializer(), line_separator))
    }

    pub fn from_output_file_default(out_file: &mut OutputFile) -> IoResult<Self> {
        Self::from_output_file(out_file, Self::DEFAULT_LINE_SEPARATOR)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::io::{TempFileData, VecStream};

    fn text(bytes: &[u8], separator: u8) -> TextInputFile<VecStream> {
        TextInputFile::new(
            Deserializer::new(VecStream::from_vec(bytes.to_vec())),
            separator,
        )
    }

    #[test]
    fn getline_handles_crlf_custom_separator_and_exact_eof_boundary() {
        let mut input = text(b"one\r|two|tail", b'|');

        input.getline().unwrap();
        assert_eq!(input.line, "one");
        assert_eq!(input.line_count, 1);
        assert!(!input.eof());
        assert!(input.is_open());

        input.getline().unwrap();
        assert_eq!(input.line, "two");
        assert!(!input.eof());

        // C++ read_to returns false after copying an unterminated final line.
        input.getline().unwrap();
        assert_eq!(input.line, "tail");
        assert!(input.eof());
        assert!(!input.is_open());
        assert_eq!(input.line_count, 3);
    }

    #[test]
    fn terminated_input_reaches_eof_only_on_the_following_getline() {
        let mut input = TextInputFile::with_default_separator(Deserializer::new(
            VecStream::from_vec(b"a\n".to_vec()),
        ));
        input.getline().unwrap();
        assert_eq!(input.line, "a");
        assert!(!input.eof());
        input.getline().unwrap();
        assert_eq!(input.line, "");
        assert!(input.eof());
        assert_eq!(input.line_count, 2);
    }

    #[test]
    fn putback_line_replays_without_reading_and_matches_size_t_wrapping() {
        let mut input = text(b"a\nb\n", b'\n');
        input.getline().unwrap();
        input.putback_line();
        assert_eq!(input.line_count, 0);
        input.getline().unwrap();
        assert_eq!(input.line, "a");
        assert_eq!(input.line_count, 1);
        input.getline().unwrap();
        assert_eq!(input.line, "b");

        let mut empty = text(b"", b'\n');
        empty.putback_line();
        assert_eq!(empty.line_count, usize::MAX);
        empty.getline().unwrap();
        assert_eq!(empty.line_count, 0);
        assert_eq!(empty.line, "");
        assert!(!empty.eof());
    }

    #[test]
    fn rewind_clears_all_text_state_and_restarts_input() {
        let mut input = text(b"x\ny", b'\n');
        input.getline().unwrap();
        input.putback_line();
        input.rewind().unwrap();
        assert_eq!(input.line, "");
        assert_eq!(input.line_count, 0);
        assert!(!input.eof());
        input.getline().unwrap();
        assert_eq!(input.line, "x");
    }

    #[test]
    fn all_three_cpp_constructor_paths_read_exact_lines() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-text-input-{}-{}.tmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let name = path.to_string_lossy().into_owned();
        std::fs::write(&path, b"named\n").unwrap();

        let mut named = TextInputFile::from_file_name_default(&name).unwrap();
        named.getline().unwrap();
        assert_eq!(named.line, "named");

        let mut output = OutputFile::new(&name, Compressor::None, "w+b").unwrap();
        output.write_raw(b"output\n").unwrap();
        let mut from_output = TextInputFile::from_output_file_default(&mut output).unwrap();
        from_output.getline().unwrap();
        assert_eq!(from_output.line, "output");
        output.close().unwrap();

        let data = TempFileData::init(false).unwrap();
        let mut temp = TempFile::from_temp_file_data(&data).unwrap();
        temp.write_raw(b"temp\n").unwrap();
        let mut from_temp = TextInputFile::from_temp_file_default(&mut temp).unwrap();
        from_temp.getline().unwrap();
        assert_eq!(from_temp.line, "temp");
        temp.close().unwrap();

        std::fs::remove_file(path).unwrap();
        std::fs::remove_file(&data.name).unwrap();
    }

    #[test]
    fn invalid_utf8_is_reported_by_the_string_compatibility_api() {
        let mut input = text(&[0xff, b'\n'], b'\n');
        assert!(matches!(input.getline(), Err(IoError::Other(_))));
        assert_eq!(input.line_count, 1);
    }
}
