//! Protein sequence (`.psq`) decoding translated from
//! `data/blastdb/psq.cpp` and the corresponding `BlastVolume` declarations in
//! `data/blastdb/volume.h`.

use crate::basic::value::{Letter, Loc, NCBI_TO_STD};

/// Decode one NCBI BLAST protein record.
///
/// BLAST records may carry a leading null sentinel and a trailing null
/// terminator. A null anywhere else is malformed. Residue bytes are converted
/// through DIAMOND's exact 28-entry `NCBI_TO_STD` table.
pub fn decode_protein_sequence(data: &[u8]) -> Result<Vec<Letter>, String> {
    let mut decoded = Vec::with_capacity(data.len());
    for (i, &aa) in data.iter().enumerate() {
        if aa == b'\0' {
            if i == 0 {
                continue;
            }
            if i == data.len() - 1 {
                break;
            }
            return Err("Unexpected null terminator in sequence data".to_string());
        }
        if usize::from(aa) >= NCBI_TO_STD.len() {
            return Err("Invalid amino acid code in sequence data".to_string());
        }
        decoded.push(NCBI_TO_STD[aa as usize]);
    }
    Ok(decoded)
}

/// C++ `BlastVolume::length`, separated from file ownership for reuse and
/// preserving the PSQ convention that each indexed span includes one sentinel.
pub fn length(sequence_index: &[u32], oid: usize) -> Loc {
    (sequence_index[oid + 1] - sequence_index[oid] - 1) as Loc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodes_every_ncbi_residue_code_exactly() {
        let encoded: Vec<u8> = (1..NCBI_TO_STD.len() as u8).collect();
        let expected: Vec<Letter> = encoded
            .iter()
            .map(|&code| NCBI_TO_STD[code as usize])
            .collect();
        assert_eq!(decode_protein_sequence(&encoded).unwrap(), expected);
    }

    #[test]
    fn applies_cpp_null_sentinel_rules() {
        assert_eq!(decode_protein_sequence(&[]).unwrap(), Vec::<Letter>::new());
        assert_eq!(decode_protein_sequence(&[0]).unwrap(), Vec::<Letter>::new());
        assert_eq!(
            decode_protein_sequence(&[0, 0]).unwrap(),
            Vec::<Letter>::new()
        );
        assert_eq!(decode_protein_sequence(&[0, 1, 2, 0]).unwrap(), vec![0, 20]);
        assert_eq!(decode_protein_sequence(&[1, 2, 0]).unwrap(), vec![0, 20]);
        assert_eq!(
            decode_protein_sequence(&[1, 0, 2]).unwrap_err(),
            "Unexpected null terminator in sequence data"
        );
    }

    #[test]
    fn rejects_codes_outside_the_ncbi_table() {
        for code in [NCBI_TO_STD.len() as u8, u8::MAX] {
            assert_eq!(
                decode_protein_sequence(&[code]).unwrap_err(),
                "Invalid amino acid code in sequence data"
            );
        }
    }

    #[test]
    fn indexed_length_excludes_one_psq_sentinel() {
        let offsets = [0, 4, 6, 12];
        assert_eq!(length(&offsets, 0), 3);
        assert_eq!(length(&offsets, 1), 1);
        assert_eq!(length(&offsets, 2), 5);
    }

    #[test]
    fn volume_module_preserves_compatibility_exports() {
        assert_eq!(
            crate::data::blastdb::volume::decode_protein_sequence(&[0, 1, 0]).unwrap(),
            vec![0]
        );
        assert_eq!(crate::data::blastdb::volume::length(&[10, 15], 0), 4);
    }
}
