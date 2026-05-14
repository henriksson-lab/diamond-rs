use std::io;

use crate::data::daa::merge_daa_files;

pub fn run(input_files: &[String], output_file: &str) -> io::Result<()> {
    merge_daa_files(input_files, output_file.to_string()).map(|_| ())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_merge_daa_command_module_loads() {}
}
