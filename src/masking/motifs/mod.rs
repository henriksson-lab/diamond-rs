use std::collections::HashSet;

use crate::basic::value::{Letter, AMINO_ACID_ALPHABET, LETTER_MASK};

mod motif_data;

/// A motif table for identifying known low-complexity patterns.
pub struct MotifTable {
    table: HashSet<[u8; motif_data::MOTIF_LEN]>,
}

impl MotifTable {
    /// Initialize the motif table from the built-in motif data.
    pub fn new() -> Self {
        let mut table = HashSet::with_capacity(motif_data::MOTIFS.len());
        for &motif in motif_data::MOTIFS {
            table.insert(*motif);
        }
        MotifTable { table }
    }

    /// Check if a kmer (as amino acid characters) is in the motif table.
    pub fn contains(&self, kmer: &[u8; motif_data::MOTIF_LEN]) -> bool {
        self.table.contains(kmer)
    }

    /// Find all motif positions in a sequence.
    ///
    /// Returns a list of positions where motifs start.
    pub fn find_motifs(&self, seq: &[Letter]) -> Vec<usize> {
        let mut positions = Vec::new();
        if seq.len() < motif_data::MOTIF_LEN {
            return positions;
        }

        for i in 0..=(seq.len() - motif_data::MOTIF_LEN) {
            let mut kmer = [0u8; motif_data::MOTIF_LEN];
            let mut valid = true;
            for j in 0..motif_data::MOTIF_LEN {
                let l = (seq[i + j] & LETTER_MASK) as usize;
                if l >= AMINO_ACID_ALPHABET.len() {
                    valid = false;
                    break;
                }
                kmer[j] = AMINO_ACID_ALPHABET[l];
            }
            if valid && self.contains(&kmer) {
                positions.push(i);
            }
        }

        positions
    }

    /// Number of motifs in the table.
    pub fn len(&self) -> usize {
        self.table.len()
    }

    pub fn is_empty(&self) -> bool {
        self.table.is_empty()
    }
}

impl Default for MotifTable {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_motif_table_init() {
        let table = MotifTable::new();
        // 8000 motifs in source, 16 duplicates → 7984 unique
        assert_eq!(table.len(), 7984);
    }

    #[test]
    fn test_motif_lookup() {
        let table = MotifTable::new();
        assert!(table.contains(b"FRKYTAFT"));
        assert!(!table.contains(b"AAAAAAAA"));
    }

    #[test]
    fn test_find_motifs_in_sequence() {
        let table = MotifTable::new();
        // Create a sequence containing the motif FRKYTAFT
        // F=13 R=1 K=11 Y=18 T=16 A=0 F=13 T=16
        let seq: Vec<Letter> = vec![0, 0, 13, 1, 11, 18, 16, 0, 13, 16, 0, 0];
        let positions = table.find_motifs(&seq);
        assert!(
            positions.contains(&2),
            "Should find FRKYTAFT at position 2"
        );
    }
}
