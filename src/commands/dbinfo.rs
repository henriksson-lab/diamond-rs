use std::fs::File;
use std::io;
use std::path::Path;

use crate::data::dmnd::ReferenceHeader;

/// Run the dbinfo command — print information about a DIAMOND database file.
pub fn run(database_path: &str) -> io::Result<()> {
    let path = Path::new(database_path);

    // Add .dmnd extension if not present
    let path = if path.extension().is_none() {
        path.with_extension("dmnd")
    } else {
        path.to_path_buf()
    };

    let mut file = File::open(&path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Error opening database file {}: {}", path.display(), e),
        )
    })?;

    let header = ReferenceHeader::read_from(&mut file)?;
    header.validate().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let db_type = "Diamond database";
    println!("{:>22}  {db_type}", "Database type");
    println!("{:>22}  {}", "Database format version", header.db_version);
    println!("{:>22}  {}", "Diamond build", header.build);
    println!("{:>22}  {}", "Sequences", header.sequences);
    println!("{:>22}  {}", "Letters", header.letters);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dbinfo_test_db() {
        let db_path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/data.dmnd"
        );
        let result = run(db_path);
        assert!(result.is_ok());
    }
}
