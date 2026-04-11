use super::value::{Letter, LETTER_MASK};

/// Compute amino acid frequencies for a sequence.
pub fn amino_acid_frequencies(seq: &[Letter]) -> [f64; 20] {
    let mut counts = [0u32; 20];
    let mut total = 0u32;
    for &l in seq {
        let l = (l & LETTER_MASK) as usize;
        if l < 20 {
            counts[l] += 1;
            total += 1;
        }
    }
    let mut freq = [0.0f64; 20];
    if total > 0 {
        for i in 0..20 {
            freq[i] = counts[i] as f64 / total as f64;
        }
    }
    freq
}

/// GC content of a nucleotide sequence.
/// Nucleotide encoding: A=0, C=1, G=2, T=3, N=4
pub fn gc_content(seq: &[Letter]) -> f64 {
    let mut gc = 0u64;
    let mut total = 0u64;
    for &l in seq {
        let l = (l & LETTER_MASK) as u8;
        if l < 4 {
            total += 1;
            if l == 1 || l == 2 {
                // C or G
                gc += 1;
            }
        }
    }
    if total == 0 {
        0.0
    } else {
        gc as f64 / total as f64
    }
}

/// Reverse complement a nucleotide sequence.
/// A(0) <-> T(3), C(1) <-> G(2), N(4) stays N(4)
pub fn reverse_complement(seq: &[Letter]) -> Vec<Letter> {
    seq.iter()
        .rev()
        .map(|&l| {
            let l = l & LETTER_MASK;
            match l {
                0 => 3, // A -> T
                1 => 2, // C -> G
                2 => 1, // G -> C
                3 => 0, // T -> A
                _ => l, // N stays N
            }
        })
        .collect()
}

/// Convert a Letter sequence to a readable string.
pub fn to_string(seq: &[Letter], alphabet: &[u8]) -> String {
    seq.iter()
        .map(|&l| {
            let idx = (l & LETTER_MASK) as usize;
            if idx < alphabet.len() {
                alphabet[idx] as char
            } else {
                'X'
            }
        })
        .collect()
}

/// Count identities between two aligned sequences.
pub fn count_identities(a: &[Letter], b: &[Letter]) -> i32 {
    a.iter()
        .zip(b.iter())
        .filter(|(&la, &lb)| (la & LETTER_MASK) == (lb & LETTER_MASK))
        .count() as i32
}

/// Compute sliding window scores between two sequences.
pub fn window_scores(
    query: &[Letter],
    subject: &[Letter],
    window_size: usize,
    score_fn: impl Fn(Letter, Letter) -> i32,
) -> Vec<i32> {
    let len = query.len().min(subject.len());
    if len < window_size {
        return vec![];
    }

    let mut scores = Vec::with_capacity(len - window_size + 1);
    let mut current_score: i32 = 0;

    // Initial window
    for i in 0..window_size {
        current_score += score_fn(query[i] & LETTER_MASK, subject[i] & LETTER_MASK);
    }
    scores.push(current_score);

    // Slide window
    for i in window_size..len {
        current_score += score_fn(query[i] & LETTER_MASK, subject[i] & LETTER_MASK);
        current_score -= score_fn(
            query[i - window_size] & LETTER_MASK,
            subject[i - window_size] & LETTER_MASK,
        );
        scores.push(current_score);
    }

    scores
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_amino_acid_frequencies() {
        let seq = vec![0i8, 0, 1, 1, 1, 2]; // AARRN
        let freq = amino_acid_frequencies(&seq);
        assert!((freq[0] - 2.0 / 6.0).abs() < 0.001); // A
        assert!((freq[1] - 3.0 / 6.0).abs() < 0.001); // R
    }

    #[test]
    fn test_reverse_complement() {
        let seq = vec![0i8, 1, 2, 3]; // ACGT
        let rc = reverse_complement(&seq);
        assert_eq!(rc, vec![0, 1, 2, 3]); // ACGT is its own reverse complement
    }

    #[test]
    fn test_reverse_complement_asym() {
        let seq = vec![0i8, 0, 0, 0]; // AAAA
        let rc = reverse_complement(&seq);
        assert_eq!(rc, vec![3, 3, 3, 3]); // TTTT
    }

    #[test]
    fn test_gc_content() {
        let seq = vec![0i8, 1, 2, 3]; // ACGT - 50% GC
        assert!((gc_content(&seq) - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_count_identities() {
        let a = vec![0i8, 1, 2, 3, 4];
        let b = vec![0i8, 1, 3, 3, 5];
        assert_eq!(count_identities(&a, &b), 3); // positions 0,1,3 match
    }

    #[test]
    fn test_to_string() {
        let seq = vec![0i8, 1, 2, 3];
        assert_eq!(to_string(&seq, crate::basic::value::AMINO_ACID_ALPHABET), "ARND");
    }

    #[test]
    fn test_window_scores() {
        let q = vec![0i8, 1, 2, 3, 4];
        let s = vec![0i8, 1, 2, 3, 4];
        let scores = window_scores(&q, &s, 3, |a, b| if a == b { 1 } else { 0 });
        assert_eq!(scores.len(), 3);
        assert_eq!(scores[0], 3); // all 3 match
    }
}
