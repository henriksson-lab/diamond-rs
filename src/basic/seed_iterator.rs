use std::collections::VecDeque;

use super::reduction::Reduction;
use super::seed::PackedSeed;
use super::sequence::Sequence;
use super::shape::Shape;
use super::value::{is_amino_acid, letter_mask, Letter, Loc};
use crate::util::hash::murmur_hash_u64;

pub struct SeedIterator<'a> {
    seq: &'a [Letter],
    ptr: usize,
    end: usize,
}

impl<'a> SeedIterator<'a> {
    pub fn new(seq: &'a [Letter], sh: &Shape) -> Self {
        let end = seq
            .len()
            .saturating_sub(sh.length as usize)
            .saturating_add(1);
        Self { seq, ptr: 0, end }
    }

    pub fn good(&self) -> bool {
        self.ptr < self.end
    }

    pub fn get(&mut self, sh: &Shape, reduction: &Reduction) -> Option<PackedSeed> {
        let seed = sh.set_seed_reduced(&self.seq[self.ptr..], reduction);
        self.ptr += 1;
        seed
    }

    pub fn increment(&mut self) -> &mut Self {
        self.ptr += 1;
        self
    }
}

#[derive(Clone)]
pub struct MinimizerIterator<'a> {
    seq: &'a [Letter],
    ptr: usize,
    begin: usize,
    end: usize,
    window: Loc,
    sh: Shape,
    reduction: &'a Reduction,
    seeds: VecDeque<u64>,
    hashes: VecDeque<u64>,
    pos: VecDeque<Loc>,
    min_idx: usize,
}

impl<'a> MinimizerIterator<'a> {
    pub fn new(seq: &'a [Letter], sh: &Shape, window: Loc, reduction: &'a Reduction) -> Self {
        let end = seq
            .len()
            .saturating_sub(sh.length as usize)
            .saturating_add(1);
        let mut it = Self {
            seq,
            ptr: 0,
            begin: 0,
            end,
            window,
            sh: *sh,
            reduction,
            seeds: VecDeque::new(),
            hashes: VecDeque::new(),
            pos: VecDeque::new(),
            min_idx: 0,
        };
        it.next();
        if it.good() {
            it.min_idx = it.get_min_idx();
        }
        it
    }

    pub fn good(&self) -> bool {
        self.seeds.len() as Loc == self.window
    }

    pub fn get(&self) -> u64 {
        self.seeds[self.min_idx]
    }

    pub fn increment(&mut self) -> &mut Self {
        let mut m = 0usize;
        let current = self.get();
        loop {
            self.seeds.pop_front();
            self.hashes.pop_front();
            self.pos.pop_front();
            self.next();
            if !self.good() {
                break;
            }
            m = self.get_min_idx();
            if self.seeds[m] != current {
                break;
            }
        }
        self.min_idx = m;
        self
    }

    pub fn pos(&self) -> Loc {
        self.pos[self.min_idx]
    }

    fn next(&mut self) {
        while (self.seeds.len() as Loc) < self.window && self.ptr < self.end {
            if let Some(s) = self
                .sh
                .set_seed_reduced(&self.seq[self.ptr..], self.reduction)
            {
                self.seeds.push_back(s);
                self.hashes.push_back(murmur_hash_u64(s));
                self.pos.push_back((self.ptr - self.begin) as Loc);
            }
            self.ptr += 1;
        }
    }

    fn get_min_idx(&self) -> usize {
        let mut j = 0usize;
        let mut s = self.hashes[0];
        for i in 1..self.hashes.len() {
            if self.hashes[i] < s {
                s = self.hashes[i];
                j = i;
            }
        }
        j
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SketchKmer {
    seed: u64,
    hash: u64,
    pos: Loc,
}

impl Ord for SketchKmer {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.hash.cmp(&other.hash)
    }
}

impl PartialOrd for SketchKmer {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

pub struct SketchIterator {
    data: Vec<SketchKmer>,
    it: usize,
}

impl SketchIterator {
    pub fn new(seq: &[Letter], sh: &Shape, n: Loc, reduction: &Reduction) -> Self {
        let mut v = Vec::new();
        let end = seq
            .len()
            .saturating_sub(sh.length as usize)
            .saturating_add(1);
        v.reserve(end);
        for p in 0..end {
            if let Some(s) = sh.set_seed_reduced(&seq[p..], reduction) {
                v.push(SketchKmer {
                    seed: s,
                    hash: murmur_hash_u64(s),
                    pos: p as Loc,
                });
            }
        }
        v.sort();
        v.truncate((n as usize).min(v.len()));
        Self { data: v, it: 0 }
    }

    pub fn good(&self) -> bool {
        self.it < self.data.len()
    }

    pub fn get(&self) -> u64 {
        self.data[self.it].seed
    }

    pub fn pos(&self) -> Loc {
        self.data[self.it].pos
    }

    pub fn increment(&mut self) -> &mut Self {
        self.it += 1;
        self
    }
}

pub struct HashedSeedIterator<'a, const B: u64> {
    long_mask: u64,
    seq: &'a [Letter],
    ptr: usize,
    end: usize,
    last: u64,
    reduction: &'a Reduction,
}

impl<'a, const B: u64> HashedSeedIterator<'a, B> {
    pub fn new(seq: &'a [Letter], len: Loc, sh: &Shape, reduction: &'a Reduction) -> Self {
        let end = len.max(0) as usize;
        let mut it = Self {
            long_mask: sh.long_mask,
            seq,
            ptr: 0,
            end,
            last: 0,
            reduction,
        };
        for _ in 0..sh.length {
            if it.ptr < it.end {
                it.last = (it.last << B) | reduction.reduce(letter_mask(it.seq[it.ptr])) as u64;
                it.ptr += 1;
            }
        }
        it
    }

    pub fn good(&self) -> bool {
        self.ptr <= self.end
    }

    pub fn get(&self) -> u64 {
        murmur_hash_u64(self.last & self.long_mask)
    }

    pub fn increment(&mut self) -> &mut Self {
        while self.ptr < self.end {
            self.last <<= B;
            let l = letter_mask(self.seq[self.ptr]);
            self.ptr += 1;
            if !is_amino_acid(l) {
                continue;
            }
            self.last |= self.reduction.reduce(l) as u64;
            return self;
        }
        self.ptr += 1;
        self
    }

    pub fn seq_ptr(&self, sh: &Shape) -> usize {
        self.ptr - sh.length as usize
    }
}

pub struct ContiguousSeedIterator<'a, const L: usize, const B: u64, const FILTER_MASKED: bool> {
    seq: &'a [Letter],
    ptr: usize,
    end: usize,
    last: u64,
    mask: u32,
    reduction: &'a Reduction,
}

impl<'a, const L: usize, const B: u64, const FILTER_MASKED: bool>
    ContiguousSeedIterator<'a, L, B, FILTER_MASKED>
{
    pub fn new(seq: Sequence<'a>, reduction: &'a Reduction) -> Self {
        let data = seq.data();
        let mut it = Self {
            seq: data,
            ptr: 0,
            end: data.len(),
            last: 0,
            mask: 0,
            reduction,
        };
        for _ in 0..L.saturating_sub(1) {
            let l = letter_mask(it.seq[it.ptr]);
            it.ptr += 1;
            it.last = (it.last << B) | reduction.reduce(l) as u64;
            if FILTER_MASKED && !is_amino_acid(l) {
                it.mask |= 1;
            }
            if FILTER_MASKED {
                it.mask <<= 1;
            }
        }
        it
    }

    pub fn good(&self) -> bool {
        self.ptr < self.end
    }

    pub fn get(&mut self) -> Option<u64> {
        loop {
            self.last <<= B;
            self.last &= (1u64 << (B * L as u64)) - 1;
            if FILTER_MASKED {
                self.mask <<= 1;
                self.mask &= (1u32 << L) - 1;
            }
            let l = letter_mask(self.seq[self.ptr]);
            self.ptr += 1;
            let r = self.reduction.reduce(l);
            self.last |= r as u64;
            if FILTER_MASKED && !is_amino_acid(l) {
                self.mask |= 1;
            }
            return if !FILTER_MASKED || self.mask == 0 {
                Some(self.last)
            } else {
                None
            };
        }
    }

    pub fn length() -> i32 {
        L as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::MASK_LETTER;

    #[test]
    fn test_seed_iterator_reduced() {
        let r = Reduction::default_reduction();
        let sh = Shape::from_code("101", &r);
        let seq = [0i8, 9, 1, 2];
        let mut it = SeedIterator::new(&seq, &sh);
        assert!(it.good());
        assert_eq!(it.get(&sh, &r), Some(1));
        assert_eq!(it.get(&sh, &r), Some(92));
        assert!(!it.good());
    }

    #[test]
    fn test_minimizer_iterator() {
        let r = Reduction::default_reduction();
        let sh = Shape::from_code("11", &r);
        let seq = [0i8, 1, 2, 3, 4];
        let mut it = MinimizerIterator::new(&seq, &sh, 3, &r);
        assert!(it.good());
        let first = it.get();
        let first_pos = it.pos();
        it.increment();
        if it.good() {
            assert!(it.pos() >= first_pos);
            assert_ne!(it.get(), first);
        }
    }

    #[test]
    fn test_sketch_iterator() {
        let r = Reduction::default_reduction();
        let sh = Shape::from_code("11", &r);
        let seq = [0i8, 1, 2, 3, 4];
        let mut it = SketchIterator::new(&seq, &sh, 2, &r);
        assert!(it.good());
        let a = (it.get(), it.pos());
        it.increment();
        assert!(it.good());
        let b = (it.get(), it.pos());
        assert_ne!(a, b);
        it.increment();
        assert!(!it.good());
    }

    #[test]
    fn test_hashed_seed_iterator() {
        let r = Reduction::default_reduction();
        let sh = Shape::from_code("111", &r);
        let seq = [0i8, 1, 2, 3];
        let mut it = HashedSeedIterator::<4>::new(&seq, seq.len() as Loc, &sh, &r);
        assert!(it.good());
        assert_eq!(it.seq_ptr(&sh), 0);
        let first = it.get();
        it.increment();
        assert!(it.good());
        assert_eq!(it.seq_ptr(&sh), 1);
        assert_ne!(it.get(), first);
        it.increment();
        assert!(!it.good());
    }

    #[test]
    fn test_contiguous_seed_iterator_filtered_and_unfiltered() {
        let r = Reduction::default_reduction();
        let seq = [0i8, 1, MASK_LETTER, 3];
        let mut unfiltered = ContiguousSeedIterator::<2, 4, false>::new(Sequence::new(&seq), &r);
        assert_eq!(ContiguousSeedIterator::<2, 4, false>::length(), 2);
        assert!(unfiltered.good());
        assert!(unfiltered.get().is_some());
        assert!(unfiltered.get().is_some());

        let mut filtered = ContiguousSeedIterator::<2, 4, true>::new(Sequence::new(&seq), &r);
        assert!(filtered.get().is_some());
        assert_eq!(filtered.get(), None);
        assert_eq!(filtered.get(), None);
    }
}
