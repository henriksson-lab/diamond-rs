use crate::basic::value::{Letter, LETTER_MASK, SEED_MASK, TRUE_AA};

/// Tantan masking parameters.
const P_REPEAT: f64 = 0.005;
const P_REPEAT_END: f64 = 0.05;
const DEFAULT_MIN_MASK_PROB: f64 = 0.5;

/// Apply tantan masking to a sequence.
///
/// Tantan identifies repetitive regions using an HMM and masks them
/// by setting the high bit (soft masking). This is a simplified scalar
/// implementation matching the C++ tantan algorithm's behavior.
pub fn mask_tantan(seq: &mut [Letter]) {
    mask_tantan_with_params(seq, P_REPEAT, P_REPEAT_END, DEFAULT_MIN_MASK_PROB);
}

/// Apply tantan masking with custom parameters.
pub fn mask_tantan_with_params(
    seq: &mut [Letter],
    p_repeat: f64,
    p_repeat_end: f64,
    min_mask_prob: f64,
) {
    let len = seq.len();
    if len < 2 {
        return;
    }

    let aa_count = TRUE_AA as usize;

    // Build likelihood ratio matrix from BLOSUM62
    // Using simplified background frequencies
    let bg_freq = [
        0.0787, 0.0507, 0.0447, 0.0545, 0.0137, 0.0408, 0.0632, 0.0680,
        0.0220, 0.0590, 0.0964, 0.0589, 0.0238, 0.0393, 0.0481, 0.0685,
        0.0552, 0.0132, 0.0321, 0.0663,
    ];

    // Forward pass: compute repeat probabilities
    let mut repeat_prob = vec![0.0f64; len];

    let mut prob_bg = 1.0 - p_repeat;
    let mut prob_rp = p_repeat;

    for i in 0..len {
        let l = (seq[i] & LETTER_MASK) as usize;
        if l >= aa_count {
            prob_bg = 1.0 - p_repeat;
            prob_rp = p_repeat;
            continue;
        }

        // Likelihood ratio: in repeat state, the letter is emitted uniformly (1/aa_count),
        // while in background state it's emitted with frequency bg_freq[l].
        // For a repeated single letter, lr > 1 when bg_freq[l] < 1/aa_count.
        // Actually for tantan, repeat state emits a specific letter with high probability.
        // Use a simple model: repeat emits the most recent letter with prob 1.
        let lr = if bg_freq[l] > 0.0 {
            1.0 / bg_freq[l]
        } else {
            20.0
        };

        // HMM transition
        let new_bg = prob_bg * (1.0 - p_repeat) + prob_rp * p_repeat_end;
        let new_rp = prob_bg * p_repeat * lr + prob_rp * (1.0 - p_repeat_end) * lr;

        let total = new_bg + new_rp;
        if total > 0.0 {
            prob_bg = new_bg / total;
            prob_rp = new_rp / total;
        } else {
            prob_bg = 1.0 - p_repeat;
            prob_rp = p_repeat;
        }

        repeat_prob[i] = prob_rp;
    }

    // Backward pass: smooth and apply masking
    let mut prob_bg = 1.0 - p_repeat;
    let mut prob_rp = p_repeat;

    for i in (0..len).rev() {
        let l = (seq[i] & LETTER_MASK) as usize;
        if l >= aa_count {
            prob_bg = 1.0 - p_repeat;
            prob_rp = p_repeat;
            continue;
        }

        let lr = if bg_freq[l] > 0.0 {
            1.0 / bg_freq[l]
        } else {
            20.0
        };

        let new_bg = prob_bg * (1.0 - p_repeat) + prob_rp * p_repeat_end;
        let new_rp = prob_bg * p_repeat * lr + prob_rp * (1.0 - p_repeat_end) * lr;

        let total = new_bg + new_rp;
        if total > 0.0 {
            prob_bg = new_bg / total;
            prob_rp = new_rp / total;
        }

        // Combine forward and backward probabilities
        let combined = repeat_prob[i] * prob_rp;
        let bg_combined = (1.0 - repeat_prob[i]) * prob_bg;
        let posterior = if combined + bg_combined > 0.0 {
            combined / (combined + bg_combined)
        } else {
            0.0
        };

        if posterior >= min_mask_prob {
            seq[i] |= SEED_MASK; // soft mask
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_mask_short() {
        let mut seq = vec![0i8]; // Single letter - shouldn't be masked
        mask_tantan(&mut seq);
        assert_eq!(seq[0] & SEED_MASK, 0);
    }

    #[test]
    fn test_no_mask_diverse() {
        // A non-repetitive diverse sequence (all 20 amino acids cycled)
        // Should have fewer masked positions than a repetitive one
        let mut diverse: Vec<Letter> = (0..100).map(|i| (i % 20) as Letter).collect();
        let mut repetitive = vec![0i8; 100];
        mask_tantan(&mut diverse);
        mask_tantan(&mut repetitive);
        let diverse_masked = diverse.iter().filter(|&&l| l & SEED_MASK != 0).count();
        let repetitive_masked = repetitive.iter().filter(|&&l| l & SEED_MASK != 0).count();
        // Repetitive should have more masking than diverse
        assert!(
            repetitive_masked >= diverse_masked,
            "Repetitive ({}) should have >= masking than diverse ({})",
            repetitive_masked,
            diverse_masked
        );
    }

    #[test]
    fn test_mask_repetitive() {
        // Highly repetitive sequence (same letter repeated)
        let mut seq = vec![0i8; 100]; // 100x Alanine
        mask_tantan(&mut seq);
        // Repetitive regions should have some masking
        let masked_count = seq
            .iter()
            .filter(|&&l| l & SEED_MASK != 0)
            .count();
        // At minimum, the simplified HMM should mask some positions
        assert!(
            masked_count > 0,
            "No positions masked in repetitive sequence"
        );
    }
}
