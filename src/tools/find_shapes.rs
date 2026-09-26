//! Safe translation of `diamond/src/tools/find_shapes.cpp`.

use std::collections::{BTreeSet, HashMap};

pub type Pattern = u32;
pub const W: usize = 12;
pub const L: usize = 24;
pub const N: usize = 64;
pub const T: usize = 6;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShapeRound {
    pub alignments: usize,
    pub total: usize,
    pub pattern: Pattern,
    pub hits: usize,
    pub excluded: Vec<Pattern>,
}

fn process_window_with(p: Pattern, patterns: &mut BTreeSet<Pattern>, weight: usize, length: usize) {
    let positions: Vec<usize> = (1..length).filter(|&bit| p & (1 << bit) != 0).collect();
    let choose = weight.saturating_sub(1);
    fn combinations(
        positions: &[usize],
        choose: usize,
        start: usize,
        pattern: Pattern,
        output: &mut BTreeSet<Pattern>,
    ) {
        if choose == 0 {
            output.insert(pattern | 1);
            return;
        }
        if positions.len().saturating_sub(start) < choose {
            return;
        }
        for index in start..=positions.len() - choose {
            combinations(
                positions,
                choose - 1,
                index + 1,
                pattern | (1 << positions[index]),
                output,
            );
        }
    }
    combinations(&positions, choose, 0, 0, patterns);
}

pub fn process_window(p: Pattern, patterns: &mut BTreeSet<Pattern>) {
    process_window_with(p, patterns, W, L);
}

fn process_pattern_with(
    sequence: &str,
    patterns: &mut BTreeSet<Pattern>,
    weight: usize,
    length: usize,
    tolerance: usize,
) {
    let mask = (1_u32 << length) - 1;
    let mut pattern = 0_u32;
    for byte in sequence.bytes() {
        pattern = (pattern << 1) & mask;
        if byte == b'1' {
            pattern |= 1;
            let count = pattern.count_ones() as usize;
            if count >= weight && count <= weight + tolerance {
                process_window_with(pattern, patterns, weight, length);
            }
        }
    }
}

pub fn process_pattern(sequence: &str, patterns: &mut BTreeSet<Pattern>) {
    process_pattern_with(sequence, patterns, W, L, T);
}

fn is_excluded_with(
    lines: &[&str],
    exclude: &BTreeSet<Pattern>,
    weight: usize,
    length: usize,
) -> bool {
    if exclude.is_empty() {
        return false;
    }
    let mask = (1_u32 << length) - 1;
    for line in lines {
        let mut pattern = 0_u32;
        for byte in line.bytes() {
            pattern <<= 1;
            if byte == b'1' {
                pattern |= 1;
            }
            pattern &= mask;
            if pattern.count_ones() as usize >= weight
                && exclude
                    .iter()
                    .any(|excluded| pattern & excluded == *excluded)
            {
                return true;
            }
        }
    }
    false
}

pub fn is_excluded(lines: &[&str], exclude: &BTreeSet<Pattern>) -> bool {
    is_excluded_with(lines, exclude, W, L)
}

pub fn is_id(sequence: &str) -> bool {
    sequence.bytes().all(|byte| byte == b'1')
}

fn process_aln_with(
    lines: &[&str],
    exclude: &BTreeSet<Pattern>,
    counts: &mut HashMap<Pattern, usize>,
    weight: usize,
    length: usize,
    tolerance: usize,
) {
    let mut patterns = BTreeSet::new();
    for line in lines {
        process_pattern_with(line, &mut patterns, weight, length, tolerance);
    }
    if patterns.is_disjoint(exclude) {
        for pattern in patterns {
            *counts.entry(pattern).or_default() += 1;
        }
    }
}

pub fn process_aln(
    lines: &[&str],
    exclude: &BTreeSet<Pattern>,
    counts: &mut HashMap<Pattern, usize>,
) {
    process_aln_with(lines, exclude, counts, W, L, T);
}

pub fn as_string(pattern: Pattern) -> String {
    if pattern == 0 {
        return String::new();
    }
    let width = (u32::BITS - pattern.leading_zeros()) as usize;
    (0..width)
        .rev()
        .map(|bit| if pattern & (1 << bit) != 0 { '1' } else { '0' })
        .collect()
}

pub fn print_all(exclude: &BTreeSet<Pattern>) -> Vec<String> {
    exclude.iter().copied().map(as_string).collect()
}

fn find_shapes_with(
    lines: &[String],
    rounds: usize,
    weight: usize,
    length: usize,
    tolerance: usize,
) -> Vec<ShapeRound> {
    let mut exclude = BTreeSet::new();
    let mut output = Vec::new();
    for _ in 0..rounds {
        let mut counts = HashMap::new();
        let mut alignments = 0;
        let mut total = 0;
        for line in lines {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() > 1 || !is_id(fields[0]) {
                if !is_excluded_with(&fields, &exclude, weight, length) {
                    process_aln_with(&fields, &exclude, &mut counts, weight, length, tolerance);
                    alignments += 1;
                }
                total += 1;
            }
        }
        let Some((&pattern, &hits)) = counts
            .iter()
            .filter(|(_, count)| **count > 0)
            .min_by_key(|(pattern, count)| (usize::MAX - **count, **pattern))
        else {
            break;
        };
        exclude.insert(pattern);
        output.push(ShapeRound {
            alignments,
            total,
            pattern,
            hits,
            excluded: exclude.iter().copied().collect(),
        });
    }
    output
}

pub fn find_shapes(lines: &[String]) -> Vec<ShapeRound> {
    find_shapes_with(lines, N, W, L, T)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn window_enumeration_selects_exact_weight_and_formats_without_padding() {
        let mut patterns = BTreeSet::new();
        process_window_with(0b1111, &mut patterns, 3, 4);
        assert_eq!(patterns, BTreeSet::from([0b0111, 0b1011, 0b1101]));
        assert_eq!(as_string(0b1011), "1011");
        assert_eq!(as_string(0), "");
    }

    #[test]
    fn exclusion_scans_sliding_windows_and_ids_are_all_ones() {
        let exclude = BTreeSet::from([0b101]);
        assert!(is_excluded_with(&["00101"], &exclude, 2, 5));
        assert!(!is_excluded_with(&["11000"], &exclude, 2, 5));
        assert!(is_id("111"));
        assert!(!is_id("101"));
        assert!(is_id(""));
    }

    #[test]
    fn iterative_finder_uses_lowest_pattern_to_break_count_ties() {
        let lines = vec!["1111".to_string(), "1011\t1101".to_string()];
        let rounds = find_shapes_with(&lines, 2, 2, 4, 2);
        assert!(!rounds.is_empty());
        assert_eq!(rounds[0].total, 1);
        assert_eq!(rounds[0].alignments, 1);
        assert!(rounds[0].hits > 0);
        assert_eq!(rounds[0].excluded, vec![rounds[0].pattern]);
    }
}
