use super::hsp::Match;

/// Target culling — remove redundant targets based on overlap.
///
/// When multiple targets align to the same query region, this module
/// decides which targets to keep based on score and coverage overlap.
pub struct TargetCulling {
    /// Overlap fraction threshold for considering targets redundant.
    pub overlap_threshold: f64,
}

impl TargetCulling {
    pub fn new(overlap_threshold: f64) -> Self {
        TargetCulling { overlap_threshold }
    }

    /// Apply culling to a list of matches, keeping only non-redundant ones.
    /// Matches should be pre-sorted by score (best first).
    pub fn cull(&self, matches: &mut Vec<Match>) {
        if matches.is_empty() {
            return;
        }

        // Keep track of which matches to retain
        let mut keep = vec![true; matches.len()];

        for i in 0..matches.len() {
            if !keep[i] {
                continue;
            }
            for j in (i + 1)..matches.len() {
                if !keep[j] {
                    continue;
                }
                // Check if match j is dominated by match i
                if self.is_dominated(&matches[i], &matches[j]) {
                    keep[j] = false;
                }
            }
        }

        let mut write_idx = 0;
        for (read_idx, &should_keep) in keep.iter().enumerate() {
            if should_keep {
                if write_idx != read_idx {
                    matches.swap(write_idx, read_idx);
                }
                write_idx += 1;
            }
        }
        matches.truncate(write_idx);
    }

    /// Check if target b is dominated by target a (a has higher score and overlapping HSPs).
    fn is_dominated(&self, a: &Match, b: &Match) -> bool {
        if a.top_score() <= b.top_score() {
            return false;
        }

        for hsp_b in &b.hsps {
            let mut is_covered = false;
            for hsp_a in &a.hsps {
                let overlap = hsp_a
                    .subject_range
                    .overlap_factor(&hsp_b.subject_range);
                if overlap >= self.overlap_threshold {
                    is_covered = true;
                    break;
                }
            }
            if !is_covered {
                return false;
            }
        }
        true
    }
}

/// Apply max-target-seqs filtering.
pub fn apply_max_target_seqs(matches: &mut Vec<Match>, max_targets: usize) {
    if matches.len() > max_targets {
        matches.truncate(max_targets);
    }
}

/// Apply top-percent filtering: keep matches within `top_percent` of the best score.
pub fn apply_top_percent(matches: &mut Vec<Match>, top_percent: f64) {
    if matches.is_empty() {
        return;
    }
    let top_score = matches[0].top_score() as f64;
    let cutoff = top_score * (1.0 - top_percent / 100.0);
    matches.retain(|m| m.top_score() as f64 >= cutoff);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;
    use crate::basic::value::Score;
    use crate::util::interval::Interval;

    fn make_match(block_id: u32, score: Score, subject_range: (i32, i32)) -> Match {
        let mut m = Match::new(block_id, block_id as u64);
        let mut hsp = Hsp::new();
        hsp.score = score;
        hsp.subject_range = Interval::new(subject_range.0, subject_range.1);
        m.hsps.push(hsp);
        m
    }

    #[test]
    fn test_max_target_seqs() {
        let mut matches = vec![
            make_match(0, 100, (0, 50)),
            make_match(1, 90, (0, 50)),
            make_match(2, 80, (0, 50)),
        ];
        apply_max_target_seqs(&mut matches, 2);
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn test_top_percent() {
        let mut matches = vec![
            make_match(0, 100, (0, 50)),
            make_match(1, 95, (0, 50)),
            make_match(2, 50, (0, 50)),
        ];
        apply_top_percent(&mut matches, 10.0);
        assert_eq!(matches.len(), 2);
    }
}
