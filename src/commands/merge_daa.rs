use std::io;

use crate::output::daa::merge::{merge_daa, DaaMergeConfig};

pub fn run(input_files: &[String], output_file: &str) -> io::Result<()> {
    let config = DaaMergeConfig::new(input_files.iter().cloned(), output_file);
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
