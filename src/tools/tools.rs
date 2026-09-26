//! Standalone helpers from `diamond/src/tools/tools.cpp`.
//!
//! File/database globals and seed enumeration are explicit inputs. This keeps
//! the three tool algorithms reusable while preserving their formatting,
//! hashing, frequency, sorting, and chunk-boundary behavior.

use crate::basic::reduction::Reduction;
use crate::basic::value::{Letter, AMINO_ACID_ALPHABET};
use crate::util::hash::murmurhash3_x64_128;
use crate::util::sequence::seqid;
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceRecord {
    pub title: String,
    pub sequence: Vec<Letter>,
}

/// Safe equivalent of the local seed-enumeration callback in `list_seeds`.
#[derive(Debug, Default)]
pub struct SeedCollector {
    pub seeds: Vec<u64>,
}

impl SeedCollector {
    pub fn call(&mut self, seed: u64, _position: usize, _shape: u32, _sequence: usize) -> bool {
        self.seeds.push(seed);
        true
    }

    pub fn finish(&mut self) {}
}

/// Translation of `split`: format sequences as FASTA chunks, switching files
/// before a record whenever the current chunk has reached the letter limit.
pub fn split(records: &[SequenceRecord], chunk_letters: usize) -> Vec<Vec<u8>> {
    let mut chunks = vec![Vec::new()];
    let mut letters = 0_usize;
    for record in records {
        if letters >= chunk_letters {
            chunks.push(Vec::new());
            letters = 0;
        }
        let output = chunks.last_mut().expect("initial chunk exists");
        output.push(b'>');
        output.extend_from_slice(seqid(&record.title).as_bytes());
        output.push(b'\n');
        for &letter in &record.sequence {
            output.push(AMINO_ACID_ALPHABET[letter as usize]);
        }
        output.push(b'\n');
        letters += record.sequence.len();
    }
    chunks
}

/// Translation of `hash_seqs`, returning the exact tabular lines instead of
/// writing process-global stdout. The commented-out upstream FASTA parser is
/// represented by the finite record slice.
pub fn hash_seqs(records: &[SequenceRecord]) -> Vec<String> {
    records
        .iter()
        .map(|record| {
            let bytes: Vec<u8> = record.sequence.iter().map(|&letter| letter as u8).collect();
            let hash = murmurhash3_x64_128(&bytes, &[0; 16]);
            let hex: String = hash.iter().map(|byte| format!("{byte:02x}")).collect();
            format!("{}\t{hex}", seqid(&record.title))
        })
        .collect()
}

/// Mean reduced background frequency for an amino-acid string.
pub fn freq(sequence: &str, reduction: &Reduction) -> f64 {
    if sequence.is_empty() {
        return f64::NAN;
    }
    let mut sum = 0.0;
    for byte in sequence.bytes() {
        let letter = AMINO_ACID_ALPHABET
            .iter()
            .position(|&candidate| candidate.eq_ignore_ascii_case(&byte))
            .unwrap_or(23) as Letter;
        sum += reduction.freq(reduction.reduce(letter));
    }
    sum / sequence.len() as f64
}

/// Final counting/ranking stage of `list_seeds`. Masking and seed enumeration
/// remain explicit upstream inputs; duplicate seeds are counted and ordered by
/// descending `(count, seed)`, exactly like sorting C++ `pair<uint64_t,u64>`
/// and iterating it in reverse.
pub fn list_seeds(
    seeds: &[u64],
    reduction: &Reduction,
    shape_weight: usize,
    query_count: usize,
) -> Vec<(u64, String)> {
    let mut counts = BTreeMap::<u64, u64>::new();
    for &seed in seeds {
        *counts.entry(seed).or_default() += 1;
    }
    let mut counts: Vec<(u64, u64)> = counts
        .into_iter()
        .map(|(seed, count)| (count, seed))
        .collect();
    counts.sort_unstable();
    counts
        .into_iter()
        .rev()
        .take(query_count)
        .map(|(count, seed)| (count, reduction.decode_seed(seed, shape_weight)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity_reduction() -> Reduction {
        Reduction::new(
            "A R N D C Q E G H I L K M F P S T W Y V",
            AMINO_ACID_ALPHABET,
        )
    }

    #[test]
    fn split_uses_seqids_and_source_boundary_rule() {
        let records = vec![
            SequenceRecord {
                title: "id one text".into(),
                sequence: vec![0, 1, 2],
            },
            SequenceRecord {
                title: "id2".into(),
                sequence: vec![3, 4],
            },
            SequenceRecord {
                title: "id3".into(),
                sequence: vec![5],
            },
        ];
        let chunks = split(&records, 3);
        assert_eq!(
            chunks,
            vec![b">id\nARN\n".to_vec(), b">id2\nDC\n>id3\nQ\n".to_vec()]
        );
    }

    #[test]
    fn hashes_are_seeded_from_zero_for_each_sequence() {
        let records = vec![
            SequenceRecord {
                title: "a description".into(),
                sequence: vec![0, 1, 2],
            },
            SequenceRecord {
                title: "b".into(),
                sequence: vec![0, 1, 2],
            },
        ];
        let hashes = hash_seqs(&records);
        assert_eq!(
            hashes[0].split_once('\t').unwrap().1,
            hashes[1].split_once('\t').unwrap().1
        );
        assert!(hashes[0].starts_with("a\t"));
        assert_eq!(hashes[0].len(), 34);
    }

    #[test]
    fn seed_counts_sort_by_count_then_seed_descending() {
        let reduction = identity_reduction();
        let mut collector = SeedCollector::default();
        for seed in [1, 2, 1, 3, 2, 3, 3] {
            assert!(collector.call(seed, 0, 0, 0));
        }
        collector.finish();
        let result = list_seeds(&collector.seeds, &reduction, 2, 2);
        assert_eq!(result, vec![(3, "AD".into()), (2, "AN".into())]);
        assert!(freq("ARND", &reduction).is_finite());
        assert!(freq("", &reduction).is_nan());
    }
}
