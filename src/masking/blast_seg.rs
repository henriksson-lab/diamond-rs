//! Safe Rust translation of `diamond/src/lib/blast/blast_seg.cpp`.
//!
//! DIAMOND stores amino acids as numeric letters.  The original NCBI SEG
//! implementation treats values 0 through 19 as the standard alphabet and
//! every other value as a bogus residue.

pub const SEG_WINDOW: i32 = 10;
pub const SEG_LOCUT: f64 = 1.8;
pub const SEG_HICUT: f64 = 2.1;

#[derive(Debug, Clone, PartialEq)]
pub struct SegParameters {
    pub window: i32,
    pub locut: f64,
    pub hicut: f64,
    pub period: i32,
    pub hilenmin: i32,
    pub overlaps: bool,
    pub maxtrim: i32,
    pub maxbogus: i32,
}

impl SegParameters {
    /// Matches C++ `SegParametersNewAa`.
    pub fn new_aa() -> Self {
        Self {
            window: SEG_WINDOW,
            locut: SEG_LOCUT,
            hicut: SEG_HICUT,
            period: 1,
            hilenmin: 0,
            overlaps: false,
            maxtrim: 50,
            maxbogus: 2,
        }
    }
}

impl Default for SegParameters {
    fn default() -> Self {
        Self::new_aa()
    }
}

#[derive(Debug, Clone)]
struct Alpha {
    alphasize: usize,
    lnalphasize: f64,
}

#[derive(Debug, Clone)]
struct SSequence<'a> {
    seq: &'a [u8],
    start: usize,
    length: usize,
    /// Exclusive end of the parent window; bounds `s_ShiftWin1` exactly as
    /// the C++ `parent->length` link does.
    parent_end: usize,
    alpha: Alpha,
    bogus: i32,
    punctuation: bool,
    composition: Vec<i32>,
    state: Vec<i32>,
    entropy: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SSeg {
    begin: i32,
    end: i32,
}

#[inline]
fn is_bogus(letter: u8) -> bool {
    letter >= 20
}

/// Matches C++ `s_SSequenceNew`; Rust ownership supplies initialized storage.
fn s_ssequence_new(sequence: &[u8]) -> SSequence<'_> {
    SSequence {
        seq: sequence,
        start: 0,
        length: sequence.len(),
        parent_end: sequence.len(),
        alpha: s_aa20alpha_std(),
        bogus: 0,
        punctuation: false,
        composition: Vec::new(),
        state: Vec::new(),
        entropy: 0.0,
    }
}

/// C++ free routines map to Rust ownership drops.
fn s_alpha_free(_alpha: Alpha) {}
fn s_ssequence_free(seq: SSequence<'_>) {
    s_alpha_free(seq.alpha);
}
fn s_seg_free(_segs: Vec<SSeg>) {}

fn s_has_dash(win: &SSequence<'_>) -> bool {
    win.seq[win.start..win.start + win.length].contains(&b'-')
}

fn s_state_cmp(a: &i32, b: &i32) -> std::cmp::Ordering {
    b.cmp(a)
}

fn s_comp_on(win: &mut SSequence<'_>) {
    win.composition = vec![0; win.alpha.alphasize];
    win.bogus = 0;
    for &letter in &win.seq[win.start..win.start + win.length] {
        if is_bogus(letter) {
            win.bogus += 1;
        } else {
            win.composition[letter as usize] += 1;
        }
    }
}

fn s_state_on(win: &mut SSequence<'_>) {
    if win.composition.is_empty() {
        s_comp_on(win);
    }
    win.state = win
        .composition
        .iter()
        .copied()
        .filter(|&count| count != 0)
        .collect();
    win.state.sort_by(s_state_cmp);
    win.state.resize(win.alpha.alphasize + 1, 0);
}

fn s_open_win<'a>(parent: &SSequence<'a>, start: i32, length: i32) -> Option<SSequence<'a>> {
    if start < 0 || length < 0 || start as usize + length as usize > parent.length {
        return None;
    }
    let mut win = SSequence {
        seq: parent.seq,
        start: parent.start + start as usize,
        length: length as usize,
        parent_end: parent.start + parent.length,
        alpha: parent.alpha.clone(),
        bogus: 0,
        punctuation: false,
        composition: Vec::new(),
        state: Vec::new(),
        entropy: -2.0,
    };
    s_state_on(&mut win);
    Some(win)
}

fn s_entropy(state: &[i32]) -> f64 {
    let total: i32 = state.iter().copied().take_while(|&v| v != 0).sum();
    if total == 0 {
        return 0.0;
    }
    let total_f = f64::from(total);
    let ent: f64 = state
        .iter()
        .copied()
        .take_while(|&v| v != 0)
        .map(|v| f64::from(v) * (f64::from(v) / total_f).ln() / std::f64::consts::LN_2)
        .sum();
    (ent / total_f).abs()
}

fn s_decrement_sv(state: &mut [i32], class: i32) {
    for i in 0..state.len().saturating_sub(1) {
        if state[i] == 0 {
            break;
        }
        if state[i] == class && state[i + 1] < class {
            state[i] -= 1;
            break;
        }
    }
}

fn s_increment_sv(state: &mut [i32], class: i32) {
    if let Some(value) = state.iter_mut().find(|value| **value == class) {
        *value += 1;
    }
}

fn s_shift_win1(win: &mut SSequence<'_>) -> bool {
    if win.start + win.length >= win.parent_end {
        return false;
    }
    let outgoing = win.seq[win.start];
    if is_bogus(outgoing) {
        win.bogus -= 1;
    } else {
        let slot = outgoing as usize;
        let old = win.composition[slot];
        s_decrement_sv(&mut win.state, old);
        win.composition[slot] -= 1;
    }

    let incoming = win.seq[win.start + win.length];
    win.start += 1;
    if is_bogus(incoming) {
        win.bogus += 1;
    } else {
        let slot = incoming as usize;
        let old = win.composition[slot];
        s_increment_sv(&mut win.state, old);
        win.composition[slot] += 1;
    }
    if win.entropy > -2.0 {
        win.entropy = s_entropy(&win.state);
    }
    true
}

fn s_close_win(_win: SSequence<'_>) {}

fn s_entropy_on(win: &mut SSequence<'_>) {
    if win.state.is_empty() {
        s_state_on(win);
    }
    win.entropy = s_entropy(&win.state);
}

fn s_seq_entropy(seq: &SSequence<'_>, window: i32, maxbogus: i32) -> Option<Vec<f64>> {
    if window < 0 || window as usize > seq.length {
        return None;
    }
    let downset = (window + 1) / 2 - 1;
    let upset = window - downset;
    let mut entropy = vec![-1.0; seq.length];
    let mut win = s_open_win(seq, 0, window)?;
    s_entropy_on(&mut win);
    let first = downset;
    let last = seq.length as i32 - upset;
    for i in first..=last {
        if (!seq.punctuation || !s_has_dash(&win)) && win.bogus <= maxbogus {
            entropy[i as usize] = win.entropy;
        }
        s_shift_win1(&mut win);
    }
    s_close_win(win);
    Some(entropy)
}

fn s_find_low(i: i32, limit: i32, hicut: f64, entropy: &[f64]) -> i32 {
    let mut j = i;
    while j >= limit && entropy[j as usize] != -1.0 && entropy[j as usize] <= hicut {
        j -= 1;
    }
    j + 1
}

fn s_find_high(i: i32, limit: i32, hicut: f64, entropy: &[f64]) -> i32 {
    let mut j = i;
    while j <= limit && entropy[j as usize] != -1.0 && entropy[j as usize] <= hicut {
        j += 1;
    }
    j - 1
}

/// Matches the six-decimal table shipped by NCBI for table-sized inputs.
fn s_lnfact(n: u32) -> f64 {
    // The source table contains entries 0..=10_000, rounded to six decimals.
    if n <= 10_000 {
        static TABLE: std::sync::OnceLock<Vec<f64>> = std::sync::OnceLock::new();
        let table = TABLE.get_or_init(|| {
            let mut values = Vec::with_capacity(10_001);
            let mut sum = 0.0;
            values.push(0.0);
            for i in 1..=10_000 {
                if i > 1 {
                    sum += (i as f64).ln();
                }
                values.push((sum * 1_000_000.0).round() / 1_000_000.0);
            }
            values
        });
        table[n as usize]
    } else {
        let n = f64::from(n);
        (n + 0.5) * n.ln() - n + 0.918_938_533_2
    }
}

fn s_ln_perm(state: &[i32], window_length: i32) -> f64 {
    let mut answer = s_lnfact(window_length as u32);
    for &count in state.iter().take_while(|&&v| v != 0) {
        answer -= s_lnfact(count as u32);
    }
    answer
}

fn s_ln_ass(state: &[i32], alphasize: i32) -> f64 {
    let mut answer = s_lnfact(alphasize as u32);
    if state[0] == 0 {
        return answer;
    }
    let mut total = alphasize;
    let mut class = 1;
    let mut previous = state[0];
    for i in 1..=alphasize as usize {
        if i == alphasize as usize {
            answer -= s_lnfact(class);
            break;
        }
        let current = state[i];
        if current == previous {
            class += 1;
        } else {
            total -= class as i32;
            answer -= s_lnfact(class);
            if current == 0 {
                answer -= s_lnfact(total as u32);
                break;
            }
            class = 1;
        }
        previous = current;
    }
    answer
}

fn s_get_prob(state: &[i32], total: i32, alpha: &Alpha) -> f64 {
    s_ln_ass(state, alpha.alphasize as i32) + s_ln_perm(state, total)
        - f64::from(total) * alpha.lnalphasize
}

fn s_trim(seq: &SSequence<'_>, leftend: &mut i32, rightend: &mut i32, params: &SegParameters) {
    let mut lend = 0;
    let mut rend = seq.length as i32 - 1;
    let minlen = 1.max(seq.length as i32 - params.maxtrim);
    let mut minprob = 1.0;
    for len in (minlen + 1..=seq.length as i32).rev() {
        let mut win = s_open_win(seq, 0, len).expect("valid SEG trim window");
        let mut i = 0;
        loop {
            let prob = s_get_prob(&win.state, len, &win.alpha);
            if prob < minprob {
                minprob = prob;
                lend = i;
                rend = len + i - 1;
            }
            if !s_shift_win1(&mut win) {
                break;
            }
            i += 1;
        }
        s_close_win(win);
    }
    *leftend += lend;
    *rightend -= seq.length as i32 - rend - 1;
}

fn s_seg_seq(seq: &SSequence<'_>, params: &mut SegParameters, segs: &mut Vec<SSeg>, offset: i32) {
    if params.window <= 0 {
        return;
    }
    params.locut = params.locut.max(0.0);
    params.hicut = params.hicut.max(0.0);
    let window = params.window;
    let downset = (window + 1) / 2 - 1;
    let upset = window - downset;
    let Some(entropy) = s_seq_entropy(seq, window, params.maxbogus) else {
        return;
    };
    let first = downset;
    let last = seq.length as i32 - upset;
    let mut lowlim = first;
    let mut i = first;
    while i <= last {
        if entropy[i as usize] <= params.locut && entropy[i as usize] != -1.0 {
            let loi = s_find_low(i, lowlim, params.hicut, &entropy);
            let hii = s_find_high(i, last, params.hicut, &entropy);
            let mut leftend = loi - downset;
            let mut rightend = hii + upset - 1;
            let temp_seq = s_open_win(seq, leftend, rightend - leftend + 1)
                .expect("valid SEG candidate window");
            s_trim(&temp_seq, &mut leftend, &mut rightend, params);
            s_close_win(temp_seq);

            if i + upset - 1 < leftend {
                let lend = loi - downset;
                let rend = leftend - 1;
                if let Some(leftseq) = s_open_win(seq, lend, rend - lend + 1) {
                    let mut leftsegs = Vec::new();
                    s_seg_seq(&leftseq, params, &mut leftsegs, offset + lend);
                    // The C implementation links only the head returned by this
                    // recursive call before the existing list.
                    if let Some(head) = leftsegs.first().copied() {
                        segs.insert(0, head);
                    }
                }
            }
            segs.insert(
                0,
                SSeg {
                    begin: leftend + offset,
                    end: rightend + offset,
                },
            );
            i = hii.min(rightend + downset);
            lowlim = i + 1;
        }
        i += 1;
    }
}

fn s_merge_segs(sequence_length: usize, segs: &mut Vec<SSeg>) {
    if segs.is_empty() {
        return;
    }
    // `hilenmin` is hard-coded to zero in the upstream routine.  Therefore
    // only overlapping (not merely adjacent) reverse-ordered segments merge.
    let mut index = 0;
    while index + 1 < segs.len() {
        if segs[index].begin - segs[index + 1].end - 1 < 0 {
            segs[index].begin = segs[index].begin.min(segs[index + 1].begin);
            segs[index].end = segs[index].end.max(segs[index + 1].end);
            segs.remove(index + 1);
        } else {
            index += 1;
        }
    }
    let _ = sequence_length;
}

fn s_segs_to_blast_seq_loc(segs: &[SSeg], offset: i32) -> Vec<(i32, i32)> {
    // C++ traverses the reverse-ordered SSeg list and prepends each BlastSeqLoc.
    segs.iter()
        .rev()
        .map(|seg| (seg.begin + offset, seg.end + offset))
        .collect()
}

fn s_aa20alpha_std() -> Alpha {
    Alpha {
        alphasize: 20,
        lnalphasize: 2.995_732_273_553_991,
    }
}

fn s_seg_parameters_check(params: &mut SegParameters) {
    if params.window <= 0 {
        params.window = 12;
    }
    params.locut = params.locut.max(0.0);
    params.hicut = params.hicut.max(0.0);
    if params.locut > params.hicut {
        params.hicut = params.locut;
    }
    params.maxbogus = params.maxbogus.clamp(0, params.window);
    if params.period <= 0 {
        params.period = 1;
    }
    params.maxtrim = params.maxtrim.max(0);
}

/// Matches C++ `SeqBufferSeg`. Returned ranges are inclusive.
pub fn seq_buffer_seg(
    sequence: &[u8],
    offset: u32,
    params: Option<&SegParameters>,
) -> Vec<(i32, i32)> {
    let mut params = params.cloned().unwrap_or_default();
    s_seg_parameters_check(&mut params);
    let seqwin = s_ssequence_new(sequence);
    let mut segs = Vec::new();
    s_seg_seq(&seqwin, &mut params, &mut segs, 0);
    if params.overlaps {
        s_merge_segs(sequence.len(), &mut segs);
    }
    let locations = s_segs_to_blast_seq_loc(&segs, offset as i32);
    s_ssequence_free(seqwin);
    s_seg_free(segs);
    locations
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_blast_seg() {
        let p = SegParameters::new_aa();
        assert_eq!(p.window, 10);
        assert_eq!((p.locut, p.hicut), (1.8, 2.1));
        assert_eq!((p.period, p.hilenmin, p.maxtrim, p.maxbogus), (1, 0, 50, 2));
        assert!(!p.overlaps);
    }

    #[test]
    fn entropy_matches_simple_states() {
        assert_eq!(s_entropy(&[10, 0]), 0.0);
        assert!((s_entropy(&[5, 5, 0]) - 1.0).abs() < 1e-12);
        assert!((s_entropy(&[1, 1, 1, 1, 0]) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn low_complexity_homopolymer_is_masked() {
        assert_eq!(seq_buffer_seg(&[0; 20], 0, None), vec![(0, 19)]);
    }

    #[test]
    fn diverse_short_sequence_is_not_masked() {
        let sequence: Vec<u8> = (0..20).collect();
        assert!(seq_buffer_seg(&sequence, 0, None).is_empty());
    }

    #[test]
    fn offset_is_added_to_inclusive_ranges() {
        assert_eq!(seq_buffer_seg(&[3; 12], 7, None), vec![(7, 18)]);
    }

    #[test]
    fn parameter_check_matches_cpp_clamps() {
        let mut p = SegParameters {
            window: 0,
            locut: -1.0,
            hicut: -2.0,
            period: 0,
            hilenmin: 0,
            overlaps: false,
            maxtrim: -3,
            maxbogus: 99,
        };
        s_seg_parameters_check(&mut p);
        assert_eq!(p.window, 12);
        assert_eq!((p.locut, p.hicut), (0.0, 0.0));
        assert_eq!((p.period, p.maxtrim, p.maxbogus), (1, 0, 12));
    }
}
