use crate::align::hsp::{Hsp, Match};
use crate::basic::packed_transcript::EditOperation;
use crate::basic::value::{Letter, LETTER_MASK};
use crate::basic::{reduction::Reduction, shape::Shape};
use crate::data::sequence_set::LetterStringSet;
use crate::dna::alignment::{compute_ksw_cigar, compute_wfa_extension, Cigar, DnaScoring, KswMode};
use crate::dna::build_score::Blastn_Score;
use crate::dna::dna_index::Index;
use crate::dna::extension_seed_matches::merge_and_extend_seeds;
use crate::dna::seed_set_dna::SeedMatch;
use crate::util::interval::Interval;

pub const KSW2_END_BONUS: i32 = 5;
pub const KSW2_BAND: i32 = 40;
pub const WFA_CUTOFF_STEPS: i32 = 10;
pub const KSW_FLAG_R: i32 = 1;
pub const KSW_FLAG_L: i32 = 2;
pub const KSW_FLAG_G: i32 = 3;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DnaExtensionAlgo {
    Ksw,
    Wfa,
}

impl DnaExtensionAlgo {
    pub fn to_string(self) -> &'static str {
        match self {
            Self::Ksw => "ksw",
            Self::Wfa => "wfa",
        }
    }

    pub fn from_string(s: &str) -> Option<Self> {
        match s {
            "ksw" => Some(Self::Ksw),
            "wfa" => Some(Self::Wfa),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExtendedSeed {
    pub i_min_extended: i32,
    pub i_max_extended: i32,
    pub j_min_extended: i32,
    pub j_max_extended: i32,
    pub length: i32,
}

impl ExtendedSeed {
    pub fn new(i_min: i32, i_max: i32, j_min: i32, j_max: i32) -> Self {
        Self {
            i_min_extended: i_min,
            i_max_extended: i_max,
            j_min_extended: j_min,
            j_max_extended: j_max,
            length: i_max - i_min,
        }
    }
}

pub fn intersection(hit: &SeedMatch, extended: &[ExtendedSeed]) -> bool {
    if extended.is_empty() {
        return false;
    }

    extended.iter().any(|s| {
        hit.i_start() >= s.i_min_extended
            && hit.i() <= s.i_max_extended
            && hit.j_start() >= s.j_min_extended
            && hit.j() <= s.j_max_extended
    })
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct CigarShort {
    pub cigar_data: Vec<(i32, char)>,
    score: i32,
    max_query: i32,
    max_target: i32,
}

impl CigarShort {
    pub fn score(&self) -> i32 {
        self.score
    }

    pub fn max_query(&self) -> i32 {
        self.max_query
    }

    pub fn max_target(&self) -> i32 {
        self.max_target
    }
}

fn scoring(score_builder: &Blastn_Score) -> DnaScoring {
    DnaScoring {
        reward: score_builder.reward(),
        penalty: score_builder.penalty().min(-score_builder.penalty()),
        gap_open: score_builder.gap_open(),
        gap_extend: score_builder.gap_extend(),
    }
}

fn consumed(cigar: &[(i32, char)]) -> (i32, i32) {
    cigar
        .iter()
        .fold((0, 0), |(query, target), &(length, op)| match op {
            'M' | '=' | 'X' => (query + length, target + length),
            'I' => (query + length, target),
            'D' => (query, target + length),
            _ => (query, target),
        })
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KswCigar(pub CigarShort);

impl KswCigar {
    pub fn new(
        target: &[Letter],
        query: &[Letter],
        score_builder: &Blastn_Score,
        flag: i32,
        ungapped_score: i32,
        band: i32,
    ) -> Self {
        let target_bytes = target.iter().map(|&x| x as u8).collect::<Vec<_>>();
        let query_bytes = query.iter().map(|&x| x as u8).collect::<Vec<_>>();
        let mut extension = Cigar::new();
        let mode = if flag == KSW_FLAG_L {
            KswMode::Left
        } else if flag == KSW_FLAG_G {
            KswMode::Global
        } else {
            KswMode::Right
        };
        compute_ksw_cigar(
            &target_bytes,
            &query_bytes,
            scoring(score_builder),
            mode,
            &mut extension,
            0,
            band,
        );
        let mut cigar_data = extension.get_cigar_data_const().to_vec();
        let (query_consumed, target_consumed) = consumed(&cigar_data);
        if flag == KSW_FLAG_L {
            cigar_data.reverse();
            cigar_data.push((ungapped_score, 'M'));
        }
        Self(CigarShort {
            cigar_data,
            score: extension.score
                + if flag == KSW_FLAG_L {
                    ungapped_score * score_builder.reward()
                } else {
                    0
                },
            max_query: query_consumed - 1,
            max_target: target_consumed - 1,
        })
    }

    pub fn into_inner(self) -> CigarShort {
        self.0
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WfaCigar(pub CigarShort);

impl WfaCigar {
    pub fn new(
        target: &[Letter],
        query: &[Letter],
        score_builder: &Blastn_Score,
        left: bool,
        anchor_score: i32,
        band: i32,
    ) -> Result<Self, String> {
        let encode = |sequence: &[Letter]| {
            sequence
                .iter()
                .map(|&letter| match letter & LETTER_MASK {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => 'N',
                })
                .collect::<String>()
        };
        let (_, encoded) = compute_wfa_extension(
            &encode(query),
            &encode(target),
            scoring(score_builder),
            band,
        );
        let mut steps = 0i32;
        let mut score = 0i32;
        let mut cigar_data = Vec::new();
        for op in encoded.chars() {
            if op.is_ascii_digit() {
                steps = steps * 10 + op.to_digit(10).unwrap() as i32;
                continue;
            }
            if steps <= 0 || !matches!(op, '=' | 'X' | 'I' | 'D') {
                return Err(format!("WFA Cigar_short: Invalid Cigar_short Symbol {op}"));
            }
            score += match op {
                '=' => steps * score_builder.reward(),
                'X' => steps * score_builder.penalty(),
                'I' | 'D' => -(score_builder.gap_open() + steps * score_builder.gap_extend()),
                _ => unreachable!(),
            };
            cigar_data.push((steps, op));
            steps = 0;
        }
        if steps != 0 {
            return Err("WFA Cigar_short: missing operation after length".to_owned());
        }
        let (query_consumed, target_consumed) = consumed(&cigar_data);
        if left {
            cigar_data.reverse();
            cigar_data.push((anchor_score, '='));
            score += anchor_score * score_builder.reward();
        }
        Ok(Self(CigarShort {
            cigar_data,
            score,
            max_query: query_consumed - 1,
            max_target: target_consumed - 1,
        }))
    }

    pub fn into_inner(self) -> CigarShort {
        self.0
    }
}

impl std::ops::Add for CigarShort {
    type Output = CigarShort;

    fn add(mut self, other: CigarShort) -> CigarShort {
        self.cigar_data.extend(other.cigar_data);
        self.score += other.score;
        self
    }
}

pub fn cigar_to_hsp_seed_match(
    target: &[Letter],
    query: &[Letter],
    hit: &SeedMatch,
    out: &mut Hsp,
    reverse: bool,
) {
    let mut pattern_pos = hit.i_start();
    let mut text_pos = hit.j_start();
    out.query_range.begin = pattern_pos;
    out.subject_range.begin = text_pos;

    for _ in 0..hit.ungapped_score() {
        out.push_match(target[text_pos as usize], query[pattern_pos as usize], true);
        pattern_pos += 1;
        text_pos += 1;
    }

    out.query_range.end = pattern_pos;
    out.subject_range.end = text_pos;
    out.transcript.push_terminator();
    out.target_seq = target.to_vec();
    out.query_source_range = out.query_range;
    out.subject_source_range = if reverse {
        Interval::new(out.subject_range.end, out.subject_range.begin)
    } else {
        Interval::new(out.subject_range.begin, out.subject_range.end)
    };
    out.frame = reverse as i32;
}

pub fn cigar_to_hsp(
    cigar: &CigarShort,
    target: &[Letter],
    query: &[Letter],
    pos_i: i32,
    pos_j: i32,
    out: &mut Hsp,
    reverse: bool,
) {
    let mut pattern_pos = pos_i - cigar.max_query() - 1;
    let mut text_pos = pos_j - cigar.max_target() - 1;
    out.query_range.begin = pattern_pos;
    out.subject_range.begin = text_pos;

    for operation in &cigar.cigar_data {
        match operation.1 {
            'M' | '=' | 'X' => {
                for _ in 0..operation.0 {
                    out.push_match(target[text_pos as usize], query[pattern_pos as usize], true);
                    pattern_pos += 1;
                    text_pos += 1;
                }
            }
            'D' => {
                let end = (text_pos + operation.0) as usize;
                out.push_gap(EditOperation::Deletion, operation.0, &target[..end]);
                text_pos += operation.0;
            }
            'I' => {
                let end = (pattern_pos + operation.0).max(0) as usize;
                out.push_gap(EditOperation::Insertion, operation.0, &query[..end]);
                pattern_pos += operation.0;
            }
            _ => {}
        }
    }

    out.query_range.end = pattern_pos;
    out.subject_range.end = text_pos;
    out.transcript.push_terminator();
    out.target_seq = target.to_vec();
    out.query_source_range = out.query_range;
    out.subject_source_range = if reverse {
        Interval::new(out.subject_range.end, out.subject_range.begin)
    } else {
        Interval::new(out.subject_range.begin, out.subject_range.end)
    };
    out.frame = reverse as i32 + 2;
}

pub type ChainingExtension<'a> = dyn Fn(&[Letter], &[Letter]) -> Vec<Match> + Send + Sync + 'a;
pub type ExtensionTiming<'a> = dyn Fn(usize, std::time::Duration) + Send + Sync + 'a;

/// Explicit counterpart of the C++ search configuration and process globals
/// consumed by this translation unit.
pub struct ExtensionConfig<'a> {
    pub targets: &'a LetterStringSet,
    pub index: &'a Index,
    pub seed_shape: &'a Shape,
    pub reduction: &'a Reduction,
    pub score_builder: &'a Blastn_Score,
    pub minimizer_window: i32,
    pub kmer_size: i32,
    pub max_evalue: f64,
    pub dna_extension: DnaExtensionAlgo,
    pub chaining_out: bool,
    pub align_long_reads: bool,
    pub chaining_extension: &'a ChainingExtension<'a>,
    pub timing: Option<&'a ExtensionTiming<'a>>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ExtensionStats;

pub fn target_extension(
    cfg: &ExtensionConfig<'_>,
    id: u32,
    query: &[Letter],
    hits: &[SeedMatch],
    reverse: bool,
) -> Result<Match, String> {
    let target = cfg.targets.get(id as usize);
    let mut extended = Vec::new();
    let mut result = Match::new_extension(id, target, None, 0, 0, f64::MAX);

    for hit in hits {
        if intersection(hit, &extended) {
            continue;
        }

        let mut out = Hsp::new();
        if hit.ungapped_score() == query.len() as i32 {
            out.score = hit.ungapped_score() * cfg.score_builder.reward();
            out.bit_score = cfg.score_builder.blast_bit_score(out.score);
            out.evalue = cfg
                .score_builder
                .blast_e_value(out.score, query.len() as i32);
            if out.evalue >= cfg.max_evalue {
                continue;
            }
            cigar_to_hsp_seed_match(target, query, hit, &mut out, reverse);
        } else {
            let query_right = &query[hit.i() as usize..];
            let target_right_end = target
                .len()
                .min(hit.j() as usize + query_right.len().saturating_mul(2));
            let target_right = &target[hit.j() as usize..target_right_end];
            let band_right = KSW2_BAND.min(query_right.len().min(target_right.len()) as i32 / 3);

            let query_left = query[..hit.i_start() as usize]
                .iter()
                .rev()
                .copied()
                .collect::<Vec<_>>();
            let target_left_begin = (hit.j_start() - (query_left.len() as i32 * 2)).max(0) as usize;
            let target_left = target[target_left_begin..hit.j_start() as usize]
                .iter()
                .rev()
                .copied()
                .collect::<Vec<_>>();
            let band_left = KSW2_BAND.min(query_left.len().min(target_left.len()) as i32 / 3);

            let extension = match cfg.dna_extension {
                DnaExtensionAlgo::Wfa => {
                    WfaCigar::new(
                        &target_left,
                        &query_left,
                        cfg.score_builder,
                        true,
                        hit.ungapped_score(),
                        crate::dna::alignment::WFA_BAND_EXTENSION,
                    )?
                    .into_inner()
                        + WfaCigar::new(
                            target_right,
                            query_right,
                            cfg.score_builder,
                            false,
                            0,
                            crate::dna::alignment::WFA_BAND_EXTENSION,
                        )?
                        .into_inner()
                }
                DnaExtensionAlgo::Ksw => {
                    KswCigar::new(
                        &target_left,
                        &query_left,
                        cfg.score_builder,
                        KSW_FLAG_L,
                        hit.ungapped_score(),
                        band_left,
                    )
                    .into_inner()
                        + KswCigar::new(
                            target_right,
                            query_right,
                            cfg.score_builder,
                            KSW_FLAG_R,
                            0,
                            band_right,
                        )
                        .into_inner()
                }
            };

            out.score = extension.score();
            out.bit_score = cfg.score_builder.blast_bit_score(out.score);
            out.evalue = cfg
                .score_builder
                .blast_e_value(out.score, query.len() as i32);
            if out.evalue >= cfg.max_evalue {
                continue;
            }
            cigar_to_hsp(
                &extension,
                target,
                query,
                hit.i_start(),
                hit.j_start(),
                &mut out,
                reverse,
            );
        }

        extended.push(ExtendedSeed::new(
            out.query_range.begin,
            out.query_range.end,
            out.subject_range.begin,
            out.subject_range.end,
        ));
        result.hsps.push(out);
    }
    Ok(result)
}

pub fn query_extension(
    cfg: &ExtensionConfig<'_>,
    query: &[Letter],
    reverse: bool,
) -> Result<Vec<Match>, String> {
    let mut seed_hits = crate::dna::seed_set_dna::seed_lookup(
        query,
        cfg.targets,
        cfg.index,
        cfg.minimizer_window,
        cfg.seed_shape,
        cfg.reduction,
    );
    let targets = (0..cfg.targets.size())
        .map(|id| cfg.targets.get(id as usize).to_vec())
        .collect::<Vec<_>>();
    let ungapped_started = std::time::Instant::now();
    seed_hits = merge_and_extend_seeds(&mut seed_hits, query, &targets, cfg.kmer_size);
    if let Some(timing) = cfg.timing {
        timing(1, ungapped_started.elapsed());
    }
    seed_hits.sort_by(|a, b| {
        a.id()
            .cmp(&b.id())
            .then_with(|| b.ungapped_score().cmp(&a.ungapped_score()))
    });

    let mut matches = Vec::new();
    let extension_started = std::time::Instant::now();
    let mut begin = 0;
    while begin < seed_hits.len() {
        let id = seed_hits[begin].id();
        let mut end = begin + 1;
        while end < seed_hits.len() && seed_hits[end].id() == id {
            end += 1;
        }
        let result = target_extension(cfg, id, query, &seed_hits[begin..end], reverse)?;
        if !result.hsps.is_empty() {
            matches.push(result);
        }
        begin = end;
    }
    if let Some(timing) = cfg.timing {
        timing(4, extension_started.elapsed());
    }
    Ok(matches)
}

pub fn extend(
    cfg: &ExtensionConfig<'_>,
    query: &[Letter],
) -> Result<(Vec<Match>, ExtensionStats), String> {
    let reverse_query = crate::basic::sequence_utils::reverse_complement(query);
    if cfg.chaining_out || cfg.align_long_reads {
        return Ok((
            (cfg.chaining_extension)(query, &reverse_query),
            ExtensionStats,
        ));
    }

    let mut matches = query_extension(cfg, query, false)?;
    matches.extend(query_extension(cfg, &reverse_query, true)?);
    Ok((matches, ExtensionStats))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_loc::PackedLoc;
    use crate::basic::packed_transcript::EditOperation;
    use crate::basic::seed::{seed_partition, seed_partition_offset, seedp_mask};
    use crate::basic::value::BlockId;
    use crate::data::seed_array::SeedArrayEntry;

    fn dna(s: &[u8]) -> Vec<Letter> {
        s.iter().map(|&x| x as Letter).collect()
    }

    fn seed(i: i32, id: BlockId, j: i32, score: i32) -> SeedMatch {
        let mut s = SeedMatch::new(i, id, j, score);
        s.score(score);
        s
    }

    fn score_builder() -> Blastn_Score {
        Blastn_Score::new_with_parameters(2, -3, 5, 1, 100_000, 10, 0.625, 0.41, 1.1)
    }

    fn build_index(
        targets: &LetterStringSet,
        shape: &Shape,
        reduction: &Reduction,
        seedp_bits: i32,
    ) -> Index {
        let mut entries = vec![Vec::new(); 1usize << seedp_bits];
        let mask = seedp_mask(seedp_bits);
        for id in 0..targets.size() {
            let target = targets.get(id as usize);
            for j in 0..=target.len().saturating_sub(shape.length as usize) {
                if let Some(seed) = shape.set_seed_reduced(&target[j..], reduction) {
                    let partition = seed_partition(seed, mask) as usize;
                    entries[partition].push(SeedArrayEntry::new(
                        seed_partition_offset(seed, seedp_bits as u64),
                        PackedLoc::new(targets.position(id, j as i32) as u64),
                    ));
                }
            }
        }
        Index::new(entries, seedp_bits, 0.0)
    }

    #[test]
    fn test_dna_extension_algo_parse_and_display() {
        assert_eq!(
            DnaExtensionAlgo::from_string("ksw"),
            Some(DnaExtensionAlgo::Ksw)
        );
        assert_eq!(
            DnaExtensionAlgo::from_string("wfa"),
            Some(DnaExtensionAlgo::Wfa)
        );
        assert_eq!(DnaExtensionAlgo::from_string("bad"), None);
        assert_eq!(DnaExtensionAlgo::Wfa.to_string(), "wfa");
    }

    #[test]
    fn test_extended_seed_and_intersection() {
        let ext = vec![ExtendedSeed::new(5, 20, 10, 25)];
        assert_eq!(ext[0].length, 15);
        assert!(intersection(&seed(15, 0, 20, 5), &ext));
        assert!(!intersection(&seed(25, 0, 30, 5), &ext));
        assert!(!intersection(&seed(15, 0, 20, 5), &[]));
    }

    #[test]
    fn test_cigar_short_add_and_accessors() {
        let mut left = CigarShort::default();
        left.cigar_data.push((2, '='));
        left.score = 4;
        left.max_query = 1;
        left.max_target = 2;

        let mut right = CigarShort::default();
        right.cigar_data.push((1, 'I'));
        right.score = -3;

        let combined = left + right;
        assert_eq!(combined.cigar_data, vec![(2, '='), (1, 'I')]);
        assert_eq!(combined.score(), 1);
        assert_eq!(combined.max_query(), 1);
        assert_eq!(combined.max_target(), 2);
    }

    #[test]
    fn test_cigar_to_hsp_seed_match() {
        let target = dna(b"AACCGG");
        let query = dna(b"AACAGG");
        let hit = seed(4, 0, 4, 4);
        let mut hsp = Hsp::new();

        cigar_to_hsp_seed_match(&target, &query, &hit, &mut hsp, false);

        assert_eq!(hsp.query_range, Interval::new(0, 4));
        assert_eq!(hsp.subject_range, Interval::new(0, 4));
        assert_eq!(hsp.subject_source_range, Interval::new(0, 4));
        assert_eq!(hsp.frame, 0);
        assert_eq!(hsp.length, 4);
        assert_eq!(hsp.identities, 3);
        assert_eq!(hsp.mismatches, 1);
    }

    #[test]
    fn test_cigar_to_hsp_short() {
        let target = dna(b"AACCGGTT");
        let query = dna(b"AACAGGTT");
        let mut cigar = CigarShort::default();
        cigar.cigar_data = vec![(3, '='), (1, 'I'), (2, 'M'), (1, 'D')];
        let mut hsp = Hsp::new();

        cigar_to_hsp(&cigar, &target, &query, 1, 1, &mut hsp, true);

        assert_eq!(hsp.query_range, Interval::new(0, 6));
        assert_eq!(hsp.subject_range, Interval::new(0, 6));
        assert_eq!(hsp.subject_source_range, Interval::new(6, 0));
        assert_eq!(hsp.frame, 3);
        assert_eq!(hsp.length, 7);
        assert_eq!(hsp.gaps, 2);
        let ops = hsp.transcript.iter().collect::<Vec<_>>();
        assert_eq!(ops[0].op, EditOperation::Match);
        assert!(hsp.transcript.data().last().unwrap().is_terminator());
    }

    #[test]
    fn test_ksw_and_wfa_cigar_constructors_include_left_anchor() {
        let score = score_builder();
        let target = [0, 1];
        let query = [0, 1];

        let ksw = KswCigar::new(&target, &query, &score, KSW_FLAG_L, 3, 40).into_inner();
        assert_eq!(ksw.cigar_data, vec![(2, 'M'), (3, 'M')]);
        assert_eq!((ksw.score(), ksw.max_query(), ksw.max_target()), (10, 1, 1));

        let wfa = WfaCigar::new(&target, &query, &score, true, 3, 20)
            .unwrap()
            .into_inner();
        assert_eq!(wfa.cigar_data, vec![(2, '='), (3, '=')]);
        assert_eq!((wfa.score(), wfa.max_query(), wfa.max_target()), (10, 1, 1));
    }

    #[test]
    fn test_target_extension_full_query_fast_path_and_overlap_suppression() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("11", &reduction);
        let mut targets = LetterStringSet::new();
        targets.push_back(&[0, 1, 2, 3]);
        let index = build_index(&targets, &shape, &reduction, 2);
        let score = score_builder();
        let chaining = |_: &[Letter], _: &[Letter]| Vec::new();
        let cfg = ExtensionConfig {
            targets: &targets,
            index: &index,
            seed_shape: &shape,
            reduction: &reduction,
            score_builder: &score,
            minimizer_window: 1,
            kmer_size: 2,
            max_evalue: f64::INFINITY,
            dna_extension: DnaExtensionAlgo::Ksw,
            chaining_out: false,
            align_long_reads: false,
            chaining_extension: &chaining,
            timing: None,
        };
        let query = [0, 1, 2, 3];
        let hits = [seed(4, 0, 4, 4), seed(4, 0, 4, 4)];

        let result = target_extension(&cfg, 0, &query, &hits, true).unwrap();
        assert_eq!(result.hsps.len(), 1);
        let hsp = &result.hsps[0];
        assert_eq!(hsp.score, 8);
        assert_eq!(hsp.query_range, Interval::new(0, 4));
        assert_eq!(hsp.subject_source_range, Interval::new(4, 0));
        assert_eq!(hsp.frame, 1);
    }

    #[test]
    fn test_query_extension_groups_targets_and_extend_dispatches_both_modes() {
        let reduction = Reduction::default_reduction();
        let shape = Shape::from_code("11", &reduction);
        let mut targets = LetterStringSet::new();
        targets.push_back(&[0, 1, 2, 3]);
        targets.push_back(&[3, 0, 1, 2, 3]);
        let index = build_index(&targets, &shape, &reduction, 2);
        let score = score_builder();
        let chaining = |query: &[Letter], reverse: &[Letter]| {
            assert_eq!(query, &[0, 1, 2, 3]);
            assert_eq!(reverse, &[0, 1, 2, 3]);
            vec![Match::new(9, 9)]
        };
        let timing_stages = std::sync::Mutex::new(Vec::new());
        let timing = |stage, _: std::time::Duration| timing_stages.lock().unwrap().push(stage);
        let mut cfg = ExtensionConfig {
            targets: &targets,
            index: &index,
            seed_shape: &shape,
            reduction: &reduction,
            score_builder: &score,
            minimizer_window: 1,
            kmer_size: 2,
            max_evalue: f64::INFINITY,
            dna_extension: DnaExtensionAlgo::Ksw,
            chaining_out: false,
            align_long_reads: false,
            chaining_extension: &chaining,
            timing: Some(&timing),
        };

        let forward = query_extension(&cfg, &[0, 1, 2, 3], false).unwrap();
        assert_eq!(
            forward
                .iter()
                .map(|m| m.target_block_id)
                .collect::<Vec<_>>(),
            vec![0, 1]
        );

        let (both_strands, _) = extend(&cfg, &[0, 1, 2, 3]).unwrap();
        assert!(both_strands.len() >= forward.len());
        assert!(both_strands.iter().any(|m| m.hsps[0].frame == 0));
        assert!(both_strands.iter().any(|m| m.hsps[0].frame == 1));
        assert!(timing_stages
            .lock()
            .unwrap()
            .windows(2)
            .any(|s| s == [1, 4]));

        cfg.align_long_reads = true;
        let (chained, _) = extend(&cfg, &[0, 1, 2, 3]).unwrap();
        assert_eq!(chained.len(), 1);
        assert_eq!(chained[0].target_block_id, 9);
    }
}
