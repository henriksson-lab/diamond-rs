use std::collections::HashSet;
use std::sync::OnceLock;

use crate::basic::value::{
    is_amino_acid, letter_mask, Letter, AMINO_ACID_ALPHABET, LETTER_MASK, MASK_LETTER, SEED_MASK,
    TRUE_AA,
};

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

/// Default cap from C++ `config.cpp:600`: `("max-motif-len", 0, "", max_motif_len, 30)`.
const MAX_MOTIF_LEN: usize = 30;

/// Abort threshold from C++ `mask_motifs`: if matched coverage reaches 50%
/// the whole sequence is left unmasked — heavily repetitive sequences keep
/// their letters.
const ABORT_RATIO_NUM: usize = 1;
const ABORT_RATIO_DEN: usize = 2;

/// Lookup table built at first use: AA `Letter` value (0..19) → cached so the
/// rolling kmer code in `mask_motifs` runs without HashSet overhead per byte.
fn motif_code_set() -> &'static HashSet<u64> {
    static SET: OnceLock<HashSet<u64>> = OnceLock::new();
    SET.get_or_init(|| {
        // Build a HashSet<u64> of motif codes in base-20 over the AA alphabet,
        // matching C++ `Kmer<8>(const char* s)` in `util/kmer/kmer.h:54`. The
        // motif strings are stored as AA character bytes (e.g. b'F'); convert
        // them to letter indices via the alphabet lookup.
        let mut char_to_letter = [u8::MAX; 256];
        for (i, &c) in AMINO_ACID_ALPHABET
            .iter()
            .enumerate()
            .take(TRUE_AA as usize)
        {
            char_to_letter[c as usize] = i as u8;
            char_to_letter[(c as char).to_ascii_lowercase() as usize] = i as u8;
        }
        let mut set = HashSet::with_capacity(motif_data::MOTIFS.len());
        for m in motif_data::MOTIFS {
            let mut code = 0u64;
            for &c in *m {
                let l = char_to_letter[c as usize];
                debug_assert!(l != u8::MAX);
                code = code * TRUE_AA as u64 + l as u64;
            }
            set.insert(code);
        }
        set
    })
}

/// A recorded motif mask — the byte that was overwritten with `MASK_LETTER`
/// at `pos` was originally `original`. Used by `restore_motifs` to undo the
/// hard-mask after seed enumeration, mirroring C++'s
/// `Block::remove_soft_masking` flow.
#[derive(Debug, Clone, Copy)]
pub struct MotifMaskEntry {
    pub pos: usize,
    pub original: Letter,
}

/// Apply motif masking in place. Returns a vector of saved positions so the
/// caller can restore the original letters after seed enumeration. The mask
/// is `MASK_LETTER` (hard mask) — seed enumeration drops seeds overlapping a
/// motif because `set_seed` rejects non-AA letters at spaced positions.
///
/// Mirrors C++ `mask_motifs` (`masking/masking.cpp:110`):
/// - rolling base-20 kmer over the AA letters, restarting on any non-AA letter
/// - look up each 8-mer code in the motif table; record `[pos, pos+8)` ranges
///   on a hit
/// - if total recorded coverage ≥ 50% of the sequence length, return early
///   without modifying the sequence
/// - for each surviving range whose length ≤ `MAX_MOTIF_LEN` (= 30), overwrite
///   those positions with `MASK_LETTER` — the same hard-mask the C++
///   `MaskingTable::add` applies via `std::fill`
pub fn mask_motifs(seq: &mut [Letter]) -> Vec<MotifMaskEntry> {
    let mut saved = Vec::new();
    let len = seq.len();
    if len < motif_data::MOTIF_LEN {
        return saved;
    }
    let table = motif_code_set();
    let modulus = (TRUE_AA as u64).pow(motif_data::MOTIF_LEN as u32 - 1);

    // C++ uses `Mask::Ranges` (`masking/def.h:73-81`) which merges adjacent
    // / overlapping ranges in `push_back`: a new `[begin, end)` either extends
    // the last range (if `begin <= back.end`) or starts a new one. This makes
    // both the 50%-coverage abort and the `max_motif_len` cap operate on
    // **merged** spans. Storing raw 8-mers as separate ranges would:
    //   - overcount coverage (3 overlapping 8-mers = 24 letters, vs ~10 merged)
    //     and prematurely trigger the 50% abort on dense low-complexity input
    //   - bypass `max_motif_len` entirely (raw len = 8 always ≤ 30)
    let mut hits: Vec<(usize, usize)> = Vec::new();
    let push_merged = |hits: &mut Vec<(usize, usize)>, begin: usize, end: usize| {
        if let Some(last) = hits.last_mut() {
            if begin <= last.1 {
                if end > last.1 {
                    last.1 = end;
                }
                return;
            }
        }
        hits.push((begin, end));
    };
    let mut code: u64 = 0;
    let mut filled: usize = 0;

    for i in 0..len {
        let l = letter_mask(seq[i]);
        if is_amino_acid(l) && (l as i32) < TRUE_AA {
            if filled == motif_data::MOTIF_LEN {
                code %= modulus;
            } else {
                filled += 1;
            }
            code = code * TRUE_AA as u64 + l as u64;
            if filled == motif_data::MOTIF_LEN && table.contains(&code) {
                let start = i + 1 - motif_data::MOTIF_LEN;
                push_merged(&mut hits, start, start + motif_data::MOTIF_LEN);
            }
        } else {
            code = 0;
            filled = 0;
        }
    }

    let total: usize = hits.iter().map(|(a, b)| b - a).sum();
    if total * ABORT_RATIO_DEN >= len * ABORT_RATIO_NUM {
        return saved;
    }

    // Collect the unique positions to mask first, then save each position's
    // original letter EXACTLY ONCE before overwriting. The previous
    // `if seq[p] != MASK_LETTER { save; mask; }` pattern looks safe but
    // silently leaks `MASK_LETTER` (= 23) into the post-restore sequence when
    // motifs overlap: the second overlapping motif sees `MASK_LETTER` already
    // there, skips the save, and on restore the position never gets put back.
    // For low-complexity sequences (e.g. proteoglycans) this affects hundreds
    // of positions and quietly inflates the Hauser CBS computation later
    // because lost letters that should have been "X" (mask) end up scoring as
    // real AAs in alignment.
    use std::collections::BTreeSet;
    let mut to_mask: BTreeSet<usize> = BTreeSet::new();
    for (a, b) in hits {
        // `MAX_MOTIF_LEN` applies to merged spans now — chains of overlapping
        // motifs whose merged span exceeds the cap are skipped entirely,
        // matching C++ `if (i->second - i->first <= config.max_motif_len)`
        // (`diamond/src/masking/masking.cpp:126`).
        if b - a <= MAX_MOTIF_LEN {
            for p in a..b {
                to_mask.insert(p);
            }
        }
    }
    for p in to_mask {
        saved.push(MotifMaskEntry {
            pos: p,
            original: seq[p],
        });
        seq[p] = MASK_LETTER;
    }
    saved
}

/// Restore letters overwritten by a prior `mask_motifs` call.
///
/// When `template_len > 0`, the SEED_MASK bit is OR'd over each motif span
/// extended by `template_len - 1` letters to the left — this mirrors C++
/// `MaskingTable::remove(... template_len, add_bit_mask=true)`
/// (`diamond/src/masking/masking.cpp:91-101`), which is invoked from
/// `enum_seeds.h:223` on the query side with `mask_seeds=true` and
/// `template_len = max(shape.length_)`. Downstream DP code zeros any cell
/// whose `(ql|sl) & SEED_MASK != 0`, so an HSP traversing the trailing edge
/// of a motif scores those positions as 0 rather than against the restored
/// real letters. Pass `template_len = 0` for the reference-side restore
/// (C++ uses `mask_seeds=false` there).
pub fn restore_motifs(seq: &mut [Letter], saved: &[MotifMaskEntry], template_len: i32) {
    for e in saved {
        seq[e.pos] = e.original;
    }
    if template_len <= 0 || saved.is_empty() {
        return;
    }
    // Walk the saved positions to recover contiguous motif spans (the merged
    // [begin, end) ranges in `mask_motifs`). Each maximal run of consecutive
    // positions is one C++ `MaskingTable` entry.
    let mut run_start = saved[0].pos;
    let mut run_end = saved[0].pos + 1;
    for e in &saved[1..] {
        if e.pos == run_end {
            run_end += 1;
        } else {
            apply_seed_mask_extension(seq, run_start, run_end, template_len);
            run_start = e.pos;
            run_end = e.pos + 1;
        }
    }
    apply_seed_mask_extension(seq, run_start, run_end, template_len);
}

#[inline]
fn apply_seed_mask_extension(seq: &mut [Letter], begin: usize, end: usize, template_len: i32) {
    let i0 = (begin as i32 - template_len + 1).max(0) as usize;
    let i1 = end.min(seq.len());
    for letter in &mut seq[i0..i1] {
        *letter |= SEED_MASK;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_motif_table_init() {
        let table = MotifTable::new();
        // 1000 active motifs (the C++ source has ~8000 but ~7000 are wrapped
        // in a /* ... */ comment and never inserted). No duplicates.
        assert_eq!(table.len(), 1000);
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
        // F=13 R=1 K=11 Y=18 T=16 A=0 F=13 T=16
        let seq: Vec<Letter> = vec![0, 0, 13, 1, 11, 18, 16, 0, 13, 16, 0, 0];
        let positions = table.find_motifs(&seq);
        assert!(positions.contains(&2), "Should find FRKYTAFT at position 2");
    }

    #[test]
    fn mask_motifs_hits_known_motif() {
        // FRKYTAFT at position 8 inside a 32-letter sequence (motif covers
        // 8/32 = 25%, well below the 50% abort threshold).
        let mut seq: Vec<Letter> = vec![0; 32];
        let motif: [Letter; 8] = [13, 1, 11, 18, 16, 0, 13, 16];
        seq[8..16].copy_from_slice(&motif);
        let original = seq.clone();
        let saved = mask_motifs(&mut seq);
        assert_eq!(saved.len(), 8);
        for i in 8..16 {
            assert_eq!(seq[i], MASK_LETTER);
        }
        // restore returns the sequence to its original state.
        restore_motifs(&mut seq, &saved, 0);
        assert_eq!(seq, original);
    }

    #[test]
    fn mask_motifs_bails_when_too_much_coverage() {
        // Repeat FRKYTAFT enough times that more than 50% of the sequence
        // would be masked; expect zero positions to be touched.
        let motif: Vec<Letter> = vec![13, 1, 11, 18, 16, 0, 13, 16];
        let mut seq: Vec<Letter> = motif.iter().cloned().cycle().take(16).collect();
        let original = seq.clone();
        let saved = mask_motifs(&mut seq);
        assert!(saved.is_empty());
        assert_eq!(seq, original);
    }

    #[test]
    fn mask_motifs_short_seq_noop() {
        let mut seq: Vec<Letter> = vec![0; 5];
        assert!(mask_motifs(&mut seq).is_empty());
    }

    #[test]
    fn mask_motifs_overlapping_restore_round_trip() {
        // Use real adjacent motifs from the table. The first three entries in
        // `MOTIF_STRINGS` are "FRKYTAFT", "KYTAFTIP", "RKYTAFTI" — three
        // sliding windows over "FRKYTAFTIP". Placing that 10-letter span means
        // all three 8-mers fall inside it and the union is positions
        // [start, start+10). Before the bugfix, only the FIRST motif's letters
        // were saved; the overlapping positions in the next two motifs would
        // silently NOT be restored, leaking MASK_LETTER into the post-restore
        // sequence. The C++ `n / len >= 0.5` abort uses the raw sum of motif
        // lengths (no overlap dedup), so we need a long enough background to
        // keep `3*8 / len < 0.5`.
        let mut letter_for = [0u8; 256];
        for (i, &c) in AMINO_ACID_ALPHABET.iter().enumerate() {
            letter_for[c as usize] = i as u8;
        }
        let composite = b"FRKYTAFTIP";
        // 80-letter background → 24/80 = 30%, well under the 50% abort.
        let mut seq: Vec<Letter> = vec![0; 80];
        for (i, &c) in composite.iter().enumerate() {
            seq[i + 8] = letter_for[c as usize] as Letter;
        }
        let original = seq.clone();
        let saved = mask_motifs(&mut seq);
        // Union of the three motif spans is exactly [8, 18).
        for i in 8..18 {
            assert_eq!(seq[i], MASK_LETTER, "position {} not masked", i);
        }
        // Restore must put every masked position back to its original letter.
        restore_motifs(&mut seq, &saved, 0);
        let diffs: Vec<usize> = seq
            .iter()
            .zip(original.iter())
            .enumerate()
            .filter(|(_, (a, b))| a != b)
            .map(|(i, _)| i)
            .collect();
        assert!(
            diffs.is_empty(),
            "restore_motifs left positions changed: {:?}",
            diffs
        );
    }

    #[test]
    fn restore_motifs_sets_seed_mask_with_template_len() {
        // Mirror C++ `MaskingTable::remove(... template_len, add_bit_mask=true)`:
        // SEED_MASK should be OR'd over `[motif_begin - template_len + 1,
        // motif_begin + motif_len)`.
        let mut letter_for = [0u8; 256];
        for (i, &c) in AMINO_ACID_ALPHABET.iter().enumerate() {
            letter_for[c as usize] = i as u8;
        }
        let motif = b"FRKYTAFT";
        let mut seq: Vec<Letter> = vec![0; 64];
        for (i, &c) in motif.iter().enumerate() {
            seq[16 + i] = letter_for[c as usize] as Letter;
        }
        let saved = mask_motifs(&mut seq);
        assert!(!saved.is_empty());

        let template_len: i32 = 12;
        restore_motifs(&mut seq, &saved, template_len);

        // Letters are restored
        for (i, &c) in motif.iter().enumerate() {
            assert_eq!(seq[16 + i] & LETTER_MASK, letter_for[c as usize] as Letter);
        }
        // SEED_MASK is set on `[16 - 11, 16 + 8)` = `[5, 24)`.
        for i in 5..24 {
            assert!(seq[i] & SEED_MASK != 0, "position {} missing SEED_MASK", i);
        }
        // SEED_MASK is NOT set outside that range.
        for i in 0..5 {
            assert_eq!(
                seq[i] & SEED_MASK,
                0,
                "position {} unexpectedly has SEED_MASK",
                i
            );
        }
        for i in 24..64 {
            assert_eq!(
                seq[i] & SEED_MASK,
                0,
                "position {} unexpectedly has SEED_MASK",
                i
            );
        }
    }
}
