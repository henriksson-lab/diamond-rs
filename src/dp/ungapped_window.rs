//! Window-based ungapped score used as the stage-2 filter.
//!
//! Ports C++ `dp/ungapped_align.cpp:ungapped_window` and the per-subject driver
//! in `dp/ungapped_simd.cpp`. The window starts 16 letters before the seed
//! anchor and runs forward until a delimiter or up to `2 * window` letters
//! total. Within the window the score accumulates with a zero floor and the
//! filter keeps the maximum running score.
//!
//! This is the score C++ uses for the seed cutoff check in stage 2. It is
//! deliberately less sensitive than full x-drop extension — using x-drop here
//! lets too many marginal targets through.
use crate::basic::value::{Letter, DELIMITER_LETTER, LETTER_MASK};
use crate::stats::score_matrix::ScoreMatrix;

/// Default window from C++ `config.cpp:555`: `("ungapped-window", 0, "", ungapped_window, 48)`.
pub const UNGAPPED_WINDOW: usize = 48;

/// C++ `Util::Seq::clip(seq, 2*window, window)`: clamp the forward window so
/// the slice extends from `start` up to the first DELIMITER on or after
/// `start + anchor`, or up to `len` letters, whichever is shorter.
fn clip_len(seq: &[Letter], start: isize, len: usize, anchor: usize) -> usize {
    let mut begin = start;
    let end = start + len as isize;
    let anchor_pos = start + anchor as isize;
    loop {
        let mut p: Option<isize> = None;
        let scan_from = begin.max(0) as usize;
        if scan_from >= seq.len() {
            break;
        }
        let scan_end = (end.min(seq.len() as isize)).max(0) as usize;
        if scan_from >= scan_end {
            break;
        }
        for i in scan_from..scan_end {
            if seq[i] == DELIMITER_LETTER {
                p = Some(i as isize);
                break;
            }
        }
        match p {
            None => return (end - begin).max(0) as usize,
            Some(pi) if pi >= anchor_pos => return (pi - begin).max(0) as usize,
            Some(pi) => begin = pi + 1,
        }
    }
    (end - begin).max(0) as usize
}

/// Run an ungapped accumulating walk of `window` letters. Score is clipped at 0
/// (local) and the maximum encountered along the walk is returned.
#[inline]
pub fn ungapped_window_score(
    query: &[Letter],
    subject: &[Letter],
    q_start: isize,
    s_start: isize,
    window: usize,
    sm: &ScoreMatrix,
) -> i32 {
    let mut score = 0i32;
    let mut st = 0i32;
    for n in 0..window {
        let qi = q_start + n as isize;
        let si = s_start + n as isize;
        if qi < 0 || si < 0 {
            continue;
        }
        let qu = qi as usize;
        let su = si as usize;
        if qu >= query.len() || su >= subject.len() {
            break;
        }
        let q = query[qu] & LETTER_MASK;
        let s = subject[su] & LETTER_MASK;
        st += sm.score(q, s);
        if st < 0 {
            st = 0;
        }
        if st > score {
            score = st;
        }
    }
    score
}

/// Filter-style ungapped score: matches the C++ stage-2 driver.
/// `q_pos`/`s_pos` are the seed anchor positions. The window is centred so it
/// starts `window` letters before the anchor and runs forward `2 * window`
/// letters, clipped by sequence boundaries and delimiters.
pub fn filter_ungapped_score(
    query: &[Letter],
    subject: &[Letter],
    q_pos: usize,
    s_pos: usize,
    window: usize,
    sm: &ScoreMatrix,
) -> i32 {
    let window_left = window as isize;
    let q_clip = clip_len(
        query,
        q_pos as isize - window_left,
        2 * window,
        window_left as usize,
    );
    if q_clip == 0 {
        return 0;
    }
    let q_start = (q_pos as isize - window_left).max(0);
    let s_start = s_pos as isize - window_left + (q_start - (q_pos as isize - window_left));
    ungapped_window_score(query, subject, q_start, s_start, q_clip, sm)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sm() -> ScoreMatrix {
        ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap()
    }

    #[test]
    fn identical_segment_scores_above_threshold() {
        let seq: Vec<Letter> = (0..40).map(|i| (i % 20) as Letter).collect();
        let s = filter_ungapped_score(&seq, &seq, 20, 20, UNGAPPED_WINDOW, &sm());
        assert!(s > 30, "got {}", s);
    }

    #[test]
    fn delimiter_clips_window() {
        // delimiter mid-walk truncates the window.
        let mut q: Vec<Letter> = (0..40).map(|i| (i % 20) as Letter).collect();
        let mut s = q.clone();
        q[25] = DELIMITER_LETTER;
        s[25] = DELIMITER_LETTER;
        let s_clipped = filter_ungapped_score(&q, &s, 20, 20, UNGAPPED_WINDOW, &sm());
        let s_full = filter_ungapped_score(
            &(0..40).map(|i| (i % 20) as Letter).collect::<Vec<_>>(),
            &(0..40).map(|i| (i % 20) as Letter).collect::<Vec<_>>(),
            20,
            20,
            UNGAPPED_WINDOW,
            &sm(),
        );
        assert!(s_clipped < s_full);
    }
}
