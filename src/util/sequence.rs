use crate::basic::translate::translate_6_frames;
use crate::basic::value::{
    letter_mask, Letter, Loc, Score, ValueTraits, AMINO_ACID_COUNT, MASK_LETTER, STOP_LETTER,
    TRUE_AA,
};
use crate::util::text_buffer::{find_first_of, TextBuffer};

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AccessionParsing {
    pub uniref_prefix: i64,
    pub gi_prefix: i64,
    pub prefix_before_pipe: i64,
    pub suffix_after_pipe: i64,
    pub suffix_after_dot: i64,
    pub pdb_suffix: i64,
}

pub const ID_DELIMITERS: &str = " \u{0007}\u{0008}\u{000c}\n\r\t\u{000b}\u{0001}";
pub const FASTA_HEADER_SEP: [&str; 2] = ["\u{0001}", " >"];
pub const TAB_ERR: &str = "Tabulator character in sequence title";
pub const SPACES_ERR: &str = "Leading spaces in sequence title";
pub const BLANK_ERR: &str = "Blank sequence title";

pub fn format(
    seq: &[Letter],
    id: &str,
    qual: Option<&[u8]>,
    out: &mut TextBuffer,
    format: &str,
    value_traits: &ValueTraits,
    wrap: Loc,
) -> Result<(), String> {
    if format == "fasta" {
        out.append_char('>').append_str(id).append_char('\n');
        let mut i = 0usize;
        let wrap = wrap as usize;
        while i < seq.len() {
            let end = (i + wrap).min(seq.len());
            for &letter in &seq[i..end] {
                out.push_back(value_traits.alphabet[letter_mask(letter) as usize]);
            }
            out.append_char('\n');
            i += wrap;
        }
    } else if format == "fastq" {
        out.append_char('@').append_str(id).append_char('\n');
        for &letter in seq {
            let l = letter as i32;
            if (l & 128) == 0 {
                out.push_back(value_traits.alphabet[l as usize]);
            } else {
                out.push_back(value_traits.alphabet[(l & 127) as usize].to_ascii_lowercase());
            }
        }
        out.append_str("\n+\n");
        out.write_raw(qual.unwrap_or(&[]));
        out.append_char('\n');
    } else {
        return Err("Invalid sequence file format".to_string());
    }
    Ok(())
}

pub fn clip(seq: &[Letter], anchor: i32) -> &[Letter] {
    let anchor = anchor as usize;
    let mut begin = 0usize;
    loop {
        match seq[begin..]
            .iter()
            .position(|&x| x == crate::basic::value::DELIMITER_LETTER)
        {
            None => return &seq[begin..],
            Some(p) => {
                let p = begin + p;
                if p >= anchor {
                    return &seq[begin..p];
                }
                begin = p + 1;
            }
        }
    }
}

pub fn seqid(title: &str) -> String {
    let end = title
        .find(|c| ID_DELIMITERS.contains(c))
        .unwrap_or(title.len());
    title[..end].to_string()
}

pub fn get_accession(title: &str, stat: &mut AccessionParsing) -> String {
    let mut t = title.to_string();
    if t.starts_with("UniRef") {
        let n = t.find('_').map(|i| i + 1).unwrap_or(0);
        t.drain(..n);
        stat.uniref_prefix += 1;
    } else if let Some(mut i) = t.find('|') {
        if t.starts_with("gi|") {
            let erase_to = t[i + 1..].find('|').map(|j| i + 1 + j + 1).unwrap_or(0);
            t.drain(..erase_to);
            i = t.find('|').unwrap_or(t.len());
            stat.gi_prefix += 1;
        }
        t.drain(..i + 1);
        stat.prefix_before_pipe += 1;
        if let Some(j) = t.find('|') {
            t.truncate(j);
            stat.suffix_after_pipe += 1;
        }
    }
    if let Some(i) = t.rfind('.') {
        t.truncate(i);
        stat.suffix_after_dot += 1;
    }
    t
}

pub fn accession_from_title(
    title: &str,
    parse_seqids: bool,
    stat: &mut AccessionParsing,
) -> Vec<String> {
    let mut out = Vec::new();
    for token in title.split('\u{0001}') {
        for part in token.split(" >") {
            if part.is_empty() {
                continue;
            }
            let id = seqid(part);
            out.push(if parse_seqids {
                get_accession(&id, stat)
            } else {
                id
            });
        }
    }
    out
}

pub fn all_seqids(s: &str) -> String {
    let mut out = Vec::new();
    for token in s.split('\u{0001}') {
        for part in token.split(" >") {
            if !part.is_empty() {
                out.push(seqid(part));
            }
        }
    }
    out.join("\u{0001}")
}

pub fn fix_title(s: &mut String) -> Option<&'static str> {
    let mut i = 0usize;
    let bytes = s.as_bytes();
    while i < bytes.len() && bytes[i] < 33 {
        i += 1;
    }
    let mut r = None;
    if i > 0 {
        s.drain(..i);
        r = Some(SPACES_ERR);
    }
    if s.is_empty() {
        s.push_str("N/A");
        return Some(BLANK_ERR);
    }
    let mut j = 0usize;
    while j < s.len() {
        if s.as_bytes()[j] == b'\t' {
            s.replace_range(j..j + 1, "\\t");
            r = Some(TAB_ERR);
            j += 2;
        } else {
            j += 1;
        }
    }
    r
}

pub fn get_title_def(s: &str) -> (String, String) {
    let i = find_first_of(s, ID_DELIMITERS);
    let title = s[..i].to_string();
    let def = if i >= s.len() {
        String::new()
    } else {
        s[i + 1..].to_string()
    };
    (title, def)
}

pub fn is_fully_masked(seq: &[Letter]) -> bool {
    let mut n = 0usize;
    for &letter in seq {
        if letter >= TRUE_AA as Letter {
            n += 1;
        }
    }
    n == seq.len()
}

pub fn translate(seq: &[Letter]) -> [Vec<Letter>; 6] {
    if seq.len() < 3 {
        return Default::default();
    }
    let dna: Vec<u8> = seq.iter().map(|&x| x as u8).collect();
    translate_6_frames(&dna).map(|frame| frame.into_iter().map(|x| x as Letter).collect())
}

pub fn find_orfs(seq: &mut [Letter], min_len: Loc) -> Loc {
    let mut begin = 0usize;
    let mut n = 0;
    while let Some(offset) = seq[begin..].iter().position(|&x| x == STOP_LETTER) {
        let it = begin + offset;
        let l = (it - begin) as Loc;
        if l < min_len {
            seq[begin..it].fill(MASK_LETTER);
        } else {
            n += l;
        }
        begin = it + 1;
    }
    let l = (seq.len() - begin) as Loc;
    if l < min_len {
        seq[begin..].fill(MASK_LETTER);
    } else {
        n += l;
    }
    n
}

pub fn looks_like_dna(seq: &[Letter], value_traits: &ValueTraits) -> Result<bool, String> {
    let mut count = [0i32; AMINO_ACID_COUNT];
    for &letter in seq {
        count[letter_mask(letter) as usize] += 1;
    }
    let a = value_traits.from_char.convert(b'A')? as usize;
    let c = value_traits.from_char.convert(b'C')? as usize;
    let g = value_traits.from_char.convert(b'G')? as usize;
    let t = value_traits.from_char.convert(b'T')? as usize;
    let n = value_traits.from_char.convert(b'N')? as usize;
    Ok((count[a] + count[c] + count[g] + count[t] + count[n]) as usize == seq.len())
}

pub fn window_scores<F>(
    seq1: &[Letter],
    seq2: &[Letter],
    window: Loc,
    score_matrix: F,
) -> Vec<Score>
where
    F: Fn(Letter, Letter) -> Score,
{
    assert_eq!(seq1.len(), seq2.len());
    let mut v = Vec::with_capacity(seq1.len());
    let mut s = 0;
    let l = seq1.len().min(window as usize);
    for i in 0..l {
        s += score_matrix(letter_mask(seq1[i]), letter_mask(seq2[i]));
        v.push(s);
    }
    let mut j = 0usize;
    for i in window as usize..seq1.len() {
        s += score_matrix(letter_mask(seq1[i]), letter_mask(seq2[i]));
        s -= score_matrix(letter_mask(seq1[j]), letter_mask(seq2[j]));
        v.push(s);
        j += 1;
    }
    v
}

pub fn from_string(s: &str, t: &ValueTraits, _line: i64) -> Result<Vec<Letter>, String> {
    let mut v = Vec::with_capacity(s.len());
    for c in s.bytes() {
        v.push(t.from_char.convert(c)?);
    }
    Ok(v)
}

pub fn from_string_into(
    s: &str,
    out: &mut Vec<Letter>,
    t: &ValueTraits,
    _line: i64,
) -> Result<(), String> {
    out.clear();
    out.reserve(s.len());
    for c in s.bytes() {
        if c == b'\n' || c == b'\r' {
            continue;
        }
        out.push(t.from_char.convert(c)?);
    }
    Ok(())
}

pub fn remove_newlines(s: &str) -> String {
    let mut r = String::with_capacity(s.len());
    for c in s.chars() {
        if c != '\n' && c != '\r' {
            r.push(c);
        }
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::{
        SequenceType, AMINO_ACID_ALPHABET, MASK_LETTER, NUCLEOTIDE_ALPHABET,
    };

    #[test]
    fn test_seqid_and_all_seqids() {
        assert_eq!(seqid("sp|P12345|ABC desc"), "sp|P12345|ABC");
        assert_eq!(seqid("abc\tdef"), "abc");
        assert_eq!(
            all_seqids("a desc\u{0001}b more >c tail"),
            "a\u{0001}b\u{0001}c"
        );
    }

    #[test]
    fn test_get_accession_cpp_rules() {
        let mut s = AccessionParsing::default();
        assert_eq!(get_accession("UniRef90_P12345.1", &mut s), "P12345");
        assert_eq!(s.uniref_prefix, 1);
        assert_eq!(s.suffix_after_dot, 1);
        assert_eq!(get_accession("gi|1|ref|XP_1.2|", &mut s), "XP_1");
        assert_eq!(s.gi_prefix, 1);
        assert_eq!(s.prefix_before_pipe, 1);
        assert_eq!(s.suffix_after_pipe, 1);
        assert_eq!(get_accession("sp|Q9ABC1|NAME", &mut s), "Q9ABC1");
    }

    #[test]
    fn test_accession_from_title() {
        let mut s = AccessionParsing::default();
        assert_eq!(
            accession_from_title("sp|Q1|x text\u{0001}UniRef50_A0.1", true, &mut s),
            vec!["Q1".to_string(), "A0".to_string()]
        );
        assert_eq!(
            accession_from_title("sp|Q1|x text", false, &mut s),
            vec!["sp|Q1|x".to_string()]
        );
    }

    #[test]
    fn test_format_fasta_and_fastq() {
        let traits = ValueTraits::new(
            AMINO_ACID_ALPHABET,
            MASK_LETTER,
            b"UO-",
            SequenceType::AminoAcid,
        );
        let mut out = TextBuffer::new();
        format(
            &[0, 1, 2, 3, 4],
            "id desc",
            None,
            &mut out,
            "fasta",
            &traits,
            2,
        )
        .unwrap();
        assert_eq!(
            std::str::from_utf8(out.data()).unwrap(),
            ">id desc\nAR\nND\nC\n"
        );

        out.clear();
        format(
            &[0, -127, 2],
            "id",
            Some(b"abc"),
            &mut out,
            "fastq",
            &traits,
            160,
        )
        .unwrap();
        assert_eq!(
            std::str::from_utf8(out.data()).unwrap(),
            "@id\nArN\n+\nabc\n"
        );
        assert_eq!(
            format(&[0], "id", None, &mut out, "bad", &traits, 160).unwrap_err(),
            "Invalid sequence file format"
        );
    }

    #[test]
    fn test_clip_title_and_masking_utilities() {
        assert_eq!(clip(&[1, 2, 31, 3, 4, 31, 5], 4), &[3, 4]);
        assert_eq!(clip(&[1, 2, 31, 3, 4], 1), &[1, 2]);

        let mut title = "\t abc\tdef".to_string();
        assert_eq!(fix_title(&mut title), Some(TAB_ERR));
        assert_eq!(title, "abc\\tdef");
        let mut blank = "\r\n".to_string();
        assert_eq!(fix_title(&mut blank), Some(BLANK_ERR));
        assert_eq!(blank, "N/A");
        assert_eq!(
            get_title_def("abc def ghi"),
            ("abc".to_string(), "def ghi".to_string())
        );
        assert!(is_fully_masked(&[MASK_LETTER, STOP_LETTER]));
        assert!(!is_fully_masked(&[0, MASK_LETTER]));
    }

    #[test]
    fn test_find_orfs_dna_translate_and_string_conversion() {
        let mut seq = vec![0, 1, STOP_LETTER, 2, 3, 4, STOP_LETTER, 5];
        assert_eq!(find_orfs(&mut seq, 3), 3);
        assert_eq!(
            seq,
            vec![
                MASK_LETTER,
                MASK_LETTER,
                STOP_LETTER,
                2,
                3,
                4,
                STOP_LETTER,
                MASK_LETTER
            ]
        );

        let nuc = ValueTraits::new(NUCLEOTIDE_ALPHABET, 4, b"", SequenceType::Nucleotide);
        let dna = from_string("ACGTN", &nuc, 1).unwrap();
        assert!(looks_like_dna(&dna, &nuc).unwrap());
        let frames = translate(&dna);
        assert_eq!(frames.len(), 6);

        let mut out = Vec::new();
        from_string_into("AC\nGT\rN", &mut out, &nuc, 2).unwrap();
        assert_eq!(out, dna);
        assert_eq!(remove_newlines("a\nb\rc"), "abc");
    }

    #[test]
    fn test_window_scores_cpp_shape() {
        let q = vec![0, 1, 2, 3];
        let s = vec![0, 2, 2, 4];
        let scores = window_scores(&q, &s, 2, |a, b| if a == b { 2 } else { -1 });
        assert_eq!(scores, vec![2, 1, 1, 1]);
    }
}
