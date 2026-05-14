use crate::basic::sequence::Sequence;
use crate::basic::value::{letter_mask, Letter, AMINO_ACID_ALPHABET, TRUE_AA};
use crate::util::math::power;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct IdentityReduction;

pub trait KmerReduction: Copy {
    fn bit_size(&self) -> i32;
    fn size(&self) -> u64;
    fn reduce(&self, x: u64) -> u64;
}

impl KmerReduction for IdentityReduction {
    fn bit_size(&self) -> i32 {
        5
    }

    fn size(&self) -> u64 {
        20
    }

    fn reduce(&self, x: u64) -> u64 {
        x
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Kmer<const K: usize> {
    pub code: u64,
}

impl<const K: usize> Kmer<K> {
    pub fn new() -> Self {
        Self { code: 0 }
    }

    pub fn from_ascii(s: &str) -> Self {
        assert_eq!(s.len(), K);
        let mut code = 0u64;
        for ch in s.bytes() {
            code = code * TRUE_AA as u64 + amino_acid_from_char(ch) as u64;
        }
        Self { code }
    }

    pub fn as_u64(self) -> u64 {
        self.code
    }
}

impl<const K: usize> From<Kmer<K>> for u64 {
    fn from(value: Kmer<K>) -> Self {
        value.code
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KmerIterator<'a, const K: usize, R = IdentityReduction> {
    reduction: R,
    seq: &'a [Letter],
    ptr: isize,
    end: isize,
    mod_: u64,
    kmer: Kmer<K>,
}

impl<'a, const K: usize> KmerIterator<'a, K, IdentityReduction> {
    pub fn new(seq: Sequence<'a>) -> Self {
        Self::with_reduction(seq, IdentityReduction)
    }
}

impl<'a, const K: usize, R> KmerIterator<'a, K, R>
where
    R: KmerReduction,
{
    pub fn with_reduction(seq: Sequence<'a>, reduction: R) -> Self {
        let mut out = Self {
            reduction,
            seq: seq.data(),
            ptr: -1,
            end: seq.data().len() as isize,
            mod_: power(reduction.size() as usize, K - 1) as u64,
            kmer: Kmer::new(),
        };
        out.inc(0, 1);
        out
    }

    pub fn get(&self) -> Kmer<K> {
        self.kmer
    }

    pub fn good(&self) -> bool {
        self.ptr < self.end
    }

    pub fn increment(&mut self) -> &mut Self {
        self.inc((K - 1) as u64, self.mod_);
        self
    }

    pub fn offset_from_start(&self) -> i32 {
        (self.ptr + 1 - K as isize) as i32
    }

    fn inc(&mut self, mut n: u64, mod_: u64) {
        self.kmer.code %= mod_;
        loop {
            self.ptr += 1;
            if self.ptr == self.end {
                return;
            }
            let l = letter_mask(self.seq[self.ptr as usize]) as u64;
            if l < TRUE_AA as u64 {
                self.kmer.code = self.kmer.code * self.reduction.size() + self.reduction.reduce(l);
                n += 1;
            } else {
                self.kmer.code = 0;
                n = 0;
            }
            if n >= K as u64 {
                break;
            }
        }
    }
}

fn amino_acid_from_char(ch: u8) -> i8 {
    for (i, &aa) in AMINO_ACID_ALPHABET.iter().enumerate() {
        if ch == aa || ch == aa.to_ascii_lowercase() {
            return i as i8;
        }
    }
    panic!("Invalid character in sequence: '{}'", ch as char);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::MASK_LETTER;

    #[test]
    fn test_identity_reduction() {
        let r = IdentityReduction;
        assert_eq!(r.bit_size(), 5);
        assert_eq!(r.size(), 20);
        assert_eq!(r.reduce(17), 17);
    }

    #[test]
    fn test_kmer_from_ascii() {
        let k = Kmer::<3>::from_ascii("ARN");
        assert_eq!(k.code, 0 * 20 * 20 + 1 * 20 + 2);
        assert_eq!(u64::from(k), 22);
    }

    #[test]
    fn test_kmer_iterator() {
        let data = vec![0, 1, 2, 3];
        let seq = Sequence::new(&data);
        let mut it = KmerIterator::<3>::new(seq);
        assert!(it.good());
        assert_eq!(it.get().code, 22);
        assert_eq!(it.offset_from_start(), 0);
        it.increment();
        assert!(it.good());
        assert_eq!(it.get().code, 1 * 20 * 20 + 2 * 20 + 3);
        assert_eq!(it.offset_from_start(), 1);
        it.increment();
        assert!(!it.good());
    }

    #[test]
    fn test_kmer_iterator_skips_masked_runs() {
        let data = vec![0, 1, MASK_LETTER, 2, 3, 4];
        let seq = Sequence::new(&data);
        let mut it = KmerIterator::<3>::new(seq);
        assert!(it.good());
        assert_eq!(it.get().code, 2 * 20 * 20 + 3 * 20 + 4);
        assert_eq!(it.offset_from_start(), 3);
        it.increment();
        assert!(!it.good());
    }
}
