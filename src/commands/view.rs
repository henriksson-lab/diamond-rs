use std::io::{self, BufReader};
use std::path::Path;

use crate::data::daa::{DaaHeader1, DaaHeader2};

/// Run the view command — show DAA file metadata.
///
/// Full DAA record parsing (alignment viewing) uses FFI via --legacy.
pub fn run(daa_file: &str) -> io::Result<()> {
    let path = Path::new(daa_file);
    let mut file = BufReader::new(std::fs::File::open(path)?);

    let header1 = DaaHeader1::read_from(&mut file)?;
    header1.validate().map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let header2 = DaaHeader2::read_from(&mut file)?;

    let mode_str = match header2.mode {
        2 => "blastp",
        3 => "blastx",
        4 => "blastn",
        _ => "unknown",
    };

    eprintln!("DAA file: {}", daa_file);
    eprintln!("Mode: {}", mode_str);
    eprintln!("Database sequences: {}", header2.db_seqs);
    eprintln!("Database letters: {}", header2.db_letters);
    eprintln!("Query records: {}", header2.query_records);
    eprintln!("Matrix: {}", header2.matrix_name());
    eprintln!("Gap penalties: {}/{}", header2.gap_open, header2.gap_extend);

    // For full DAA viewing, fall back to FFI
    eprintln!("Note: Full DAA record parsing requires --legacy mode");

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_view_stub() {
        // View command shows DAA metadata — full record parsing uses FFI
    }
}
