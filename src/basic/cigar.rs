use super::packed_transcript::{EditOperation, PackedTranscript};

/// Generate a CIGAR string from a PackedTranscript.
///
/// Uses standard CIGAR operations: M (match/mismatch), I (insertion), D (deletion).
pub fn to_cigar(transcript: &PackedTranscript) -> String {
    let mut cigar = String::new();
    let mut current_op: Option<char> = None;
    let mut current_count = 0u32;

    for combined in transcript.iter() {
        let cigar_op = match combined.op {
            EditOperation::Match => 'M',
            EditOperation::Substitution => 'M', // CIGAR uses M for both
            EditOperation::Insertion => 'I',
            EditOperation::Deletion => 'D',
            EditOperation::FrameshiftForward | EditOperation::FrameshiftReverse => 'F',
        };

        if Some(cigar_op) == current_op {
            current_count += combined.count;
        } else {
            if let Some(op) = current_op {
                cigar.push_str(&format!("{}{}", current_count, op));
            }
            current_op = Some(cigar_op);
            current_count = combined.count;
        }
    }

    if let Some(op) = current_op {
        cigar.push_str(&format!("{}{}", current_count, op));
    }

    cigar
}

/// Generate a BTOP (Blast Traceback Operations) string.
///
/// BTOP format: numbers for matches, pairs of letters for substitutions,
/// dashes for gaps.
pub fn to_btop(
    transcript: &PackedTranscript,
    query: &[i8],
    subject: &[i8],
    query_start: usize,
    subject_start: usize,
    alphabet: &[u8],
) -> String {
    let mut btop = String::new();
    let mut qi = query_start;
    let mut si = subject_start;

    for combined in transcript.iter() {
        match combined.op {
            EditOperation::Match => {
                btop.push_str(&combined.count.to_string());
                qi += combined.count as usize;
                si += combined.count as usize;
            }
            EditOperation::Substitution => {
                let ql = (query[qi] & 0x1F) as usize;
                let sl = (subject[si] & 0x1F) as usize;
                if ql < alphabet.len() && sl < alphabet.len() {
                    btop.push(alphabet[ql] as char);
                    btop.push(alphabet[sl] as char);
                }
                qi += 1;
                si += 1;
            }
            EditOperation::Insertion => {
                for _ in 0..combined.count {
                    btop.push('-');
                    if qi < query.len() {
                        let ql = (query[qi] & 0x1F) as usize;
                        if ql < alphabet.len() {
                            btop.push(alphabet[ql] as char);
                        }
                    }
                    qi += 1;
                }
            }
            EditOperation::Deletion => {
                if si < subject.len() {
                    let sl = (subject[si] & 0x1F) as usize;
                    if sl < alphabet.len() {
                        btop.push(alphabet[sl] as char);
                    }
                }
                btop.push('-');
                si += 1;
            }
            _ => {}
        }
    }

    btop
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_all_match() {
        let mut t = PackedTranscript::new();
        t.push_with_count(EditOperation::Match, 10);
        t.push_terminator();
        assert_eq!(to_cigar(&t), "10M");
    }

    #[test]
    fn test_cigar_mixed() {
        let mut t = PackedTranscript::new();
        t.push_with_count(EditOperation::Match, 5);
        t.push_with_letter(EditOperation::Substitution, 3);
        t.push_with_count(EditOperation::Match, 3);
        t.push_terminator();
        // Substitutions map to M in CIGAR, so 5+1+3 = 9M
        assert_eq!(to_cigar(&t), "9M");
    }

    #[test]
    fn test_cigar_with_gaps() {
        let mut t = PackedTranscript::new();
        t.push_with_count(EditOperation::Match, 5);
        t.push_with_count(EditOperation::Insertion, 2);
        t.push_with_count(EditOperation::Match, 3);
        t.push_terminator();
        assert_eq!(to_cigar(&t), "5M2I3M");
    }

    #[test]
    fn test_cigar_empty() {
        let t = PackedTranscript::new();
        assert_eq!(to_cigar(&t), "");
    }
}
