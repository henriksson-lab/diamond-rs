use std::collections::BTreeSet;
use std::io::{self, Write};
use std::path::Path;

use crate::basic::value::AMINO_ACID_ALPHABET;
use crate::data::dmnd_reader;

/// Run the getseq command — retrieve sequences from a DIAMOND database.
pub fn run(database: &str, seq_numbers: Option<&str>) -> io::Result<()> {
    run_with_output(database, seq_numbers, None)
}

/// Retrieve sequences and optionally write them to `-o/--out`. Keeping stdout
/// as the `None` case preserves the library API while matching the CLI's
/// upstream output-file behavior.
pub fn run_with_output(
    database: &str,
    seq_numbers: Option<&str>,
    output: Option<&Path>,
) -> io::Result<()> {
    let seq_numbers: Vec<String> = seq_numbers
        .map(|numbers| {
            numbers
                .split(|c: char| c == ',' || c.is_ascii_whitespace())
                .filter(|number| !number.is_empty())
                .map(ToOwned::to_owned)
                .collect()
        })
        .unwrap_or_default();
    run_with_output_numbers(database, &seq_numbers, output)
}

/// Retrieve sequences by DIAMOND's 1-based numeric sequence numbers.
///
/// The CLI accepts a space-separated list after `--seq`. The compatibility
/// wrapper above also accepts comma-separated numbers for existing Rust API
/// callers.
pub fn run_with_output_numbers(
    database: &str,
    seq_numbers: &[String],
    output: Option<&Path>,
) -> io::Result<()> {
    let mut writer: Box<dyn Write> = if let Some(path) = output {
        Box::new(io::BufWriter::new(std::fs::File::create(path)?))
    } else {
        Box::new(io::BufWriter::new(io::stdout()))
    };

    let requested = if seq_numbers.is_empty() {
        None
    } else {
        // Upstream stores `atoi(value) - 1` in a `set<size_t>` and compares it
        // with the zero-based database ordinal. Preserve those selection and
        // deduplication semantics rather than interpreting values as seqids.
        Some(
            seq_numbers
                .iter()
                .map(|number| c_atoi(number).wrapping_sub(1) as usize)
                .collect::<BTreeSet<_>>(),
        )
    };
    // C++ formats each FASTA record in a reusable TextBuffer and hands the
    // complete record to OutputFile in one call.  Do the same here: writing
    // each encoded residue through `dyn Write` adds a virtual call per letter
    // on this output-heavy command.
    let mut record = Vec::new();
    dmnd_reader::for_each_dmnd_record_auto(database, |ordinal, id, sequence| {
        if requested
            .as_ref()
            .is_none_or(|requested| requested.contains(&ordinal))
        {
            format_fasta_record(&mut record, id, sequence);
            writer.write_all(&record)?;
        }
        Ok(())
    })?;

    writer.flush()?;
    Ok(())
}

// Match the portion of C `atoi` that is defined and relevant to DIAMOND's
// sequence-number option: leading ASCII space, an optional sign, and the
// initial decimal digit run. Values outside `int` range are undefined for
// `atoi`; saturating keeps malformed CLI input deterministic and portable.
fn c_atoi(value: &str) -> i32 {
    let bytes = value.as_bytes();
    let mut i = bytes
        .iter()
        .position(|byte| !byte.is_ascii_whitespace())
        .unwrap_or(bytes.len());
    let negative = match bytes.get(i) {
        Some(b'-') => {
            i += 1;
            true
        }
        Some(b'+') => {
            i += 1;
            false
        }
        _ => false,
    };
    let mut magnitude = 0_i64;
    let mut any = false;
    while let Some(byte) = bytes.get(i) {
        if !byte.is_ascii_digit() {
            break;
        }
        any = true;
        magnitude = magnitude
            .saturating_mul(10)
            .saturating_add(i64::from(byte - b'0'));
        i += 1;
    }
    if !any {
        return 0;
    }
    let signed = if negative { -magnitude } else { magnitude };
    signed.clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32
}

fn format_fasta_record(record: &mut Vec<u8>, id: &str, sequence: &[u8]) {
    record.clear();
    let line_count = sequence.len().div_ceil(80);
    record.reserve(1 + id.len() + 1 + sequence.len() + line_count);
    record.push(b'>');
    record.extend_from_slice(id.as_bytes());
    record.push(b'\n');

    // C++ getseq writes 80-char lines (sequence_file.cpp:get_seq).
    for chunk in sequence.chunks(80) {
        for &letter in chunk {
            let idx = (letter & 0x1F) as usize;
            record.push(AMINO_ACID_ALPHABET.get(idx).copied().unwrap_or(b'X'));
        }
        record.push(b'\n');
    }
}

#[cfg(test)]
mod tests {
    use super::{c_atoi, format_fasta_record};

    #[test]
    fn atoi_matches_upstream_sequence_number_inputs() {
        assert_eq!(c_atoi("1"), 1);
        assert_eq!(c_atoi("  +3tail"), 3);
        assert_eq!(c_atoi("-2"), -2);
        assert_eq!(c_atoi("not-a-number"), 0);
    }

    #[test]
    fn reusable_record_buffer_formats_wrapping_mask_bits_and_unknown_letters() {
        let mut sequence = vec![0; 80];
        sequence.extend_from_slice(&[0x80 | 1, 31]);
        let mut record = b"stale bytes".to_vec();
        format_fasta_record(&mut record, "seq", &sequence);

        let mut expected = b">seq\n".to_vec();
        expected.extend_from_slice(&[b'A'; 80]);
        expected.extend_from_slice(b"\nRX\n");
        assert_eq!(record, expected);
    }
}
