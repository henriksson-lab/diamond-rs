use std::io;
use std::path::{Path, PathBuf};

use crate::output::daa::merge::{merge_daa, DaaMergeConfig};

pub fn run(input_files: &[String], output_file: &str) -> io::Result<()> {
    run_with_buffer_size(input_files, output_file, 67_108_864)
}

pub fn run_with_buffer_size(
    input_files: &[String],
    output_file: &str,
    file_buffer_size: usize,
) -> io::Result<()> {
    run_with_buffer_size_paths(
        &input_files.iter().map(PathBuf::from).collect::<Vec<_>>(),
        Path::new(output_file),
        file_buffer_size,
    )
}

pub fn run_with_buffer_size_paths(
    input_files: &[PathBuf],
    output_file: &Path,
    file_buffer_size: usize,
) -> io::Result<()> {
    let mut config = DaaMergeConfig::new(input_files.iter().cloned(), output_file);
    config.file_buffer_size = file_buffer_size;
    merge_daa(&config).map(|_| ())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_daa_command_module_loads() {}

    #[test]
    fn run_rejects_missing_cpp_parameters() {
        let error = run(&[], "out.daa").unwrap_err();
        assert_eq!(error.to_string(), "Missing parameter: input files (--in)");

        let error = run(&["in.daa".to_string()], "").unwrap_err();
        assert_eq!(error.to_string(), "Missing parameter: output file (--out)");
    }
}
