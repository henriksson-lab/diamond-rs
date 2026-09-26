use crate::basic::value::ValueTraits;
use crate::data::block::{AlignModeBlock, Block};
use crate::data::seed_set::{HashedSeedSet, SeedSet};
use crate::util::sequence;
use crate::util::text_buffer::TextBuffer;
use std::io::Write;
use std::sync::{LazyLock, Mutex};

pub static QUERY_ALIGNED: LazyLock<Mutex<Vec<bool>>> = LazyLock::new(|| Mutex::new(Vec::new()));
pub static QUERY_SEEDS_HASHED: LazyLock<Mutex<Option<HashedSeedSet>>> =
    LazyLock::new(|| Mutex::new(None));
pub static QUERY_SEEDS_BITSET: LazyLock<Mutex<Option<SeedSet>>> =
    LazyLock::new(|| Mutex::new(None));

pub fn write_unaligned<W: Write>(
    query: &Block,
    file: &mut W,
    format: &str,
    value_traits: &ValueTraits,
    align_mode: AlignModeBlock,
) -> Result<(), String> {
    let aligned = QUERY_ALIGNED.lock().map_err(|e| e.to_string())?;
    let n = query.ids()?.size() as usize;
    let mut buf = TextBuffer::new();
    for i in 0..n {
        if !aligned.get(i).copied().unwrap_or(false) {
            let seq = if align_mode.query_translated {
                query.source_seqs().get(i)
            } else {
                query.seqs().get(i)
            };
            let id = std::str::from_utf8(query.ids()?.get(i)).map_err(|e| e.to_string())?;
            let qual = if query.qual().empty() {
                None
            } else {
                Some(query.qual().get(i))
            };
            sequence::format(seq, id, qual, &mut buf, format, value_traits, 160)?;
            file.write_all(buf.data()).map_err(|e| e.to_string())?;
            buf.clear();
        }
    }
    Ok(())
}

pub fn write_aligned<W: Write>(
    query: &Block,
    file: &mut W,
    format: &str,
    value_traits: &ValueTraits,
    align_mode: AlignModeBlock,
) -> Result<(), String> {
    let aligned = QUERY_ALIGNED.lock().map_err(|e| e.to_string())?;
    let n = query.ids()?.size() as usize;
    let mut buf = TextBuffer::new();
    for i in 0..n {
        if aligned.get(i).copied().unwrap_or(false) {
            let seq = if align_mode.query_translated {
                query.source_seqs().get(i)
            } else {
                query.seqs().get(i)
            };
            let id = std::str::from_utf8(query.ids()?.get(i)).map_err(|e| e.to_string())?;
            let qual = if query.qual().empty() {
                None
            } else {
                Some(query.qual().get(i))
            };
            sequence::format(seq, id, qual, &mut buf, format, value_traits, 160)?;
            file.write_all(buf.data()).map_err(|e| e.to_string())?;
            buf.clear();
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::value::{
        SequenceType, AMINO_ACID_ALPHABET, MASK_LETTER, NUCLEOTIDE_ALPHABET,
    };

    fn amino_traits() -> ValueTraits {
        ValueTraits::new(
            AMINO_ACID_ALPHABET,
            MASK_LETTER,
            b"UO-",
            SequenceType::AminoAcid,
        )
    }

    fn nucleotide_traits() -> ValueTraits {
        ValueTraits::new(NUCLEOTIDE_ALPHABET, 4, b"", SequenceType::Nucleotide)
    }

    #[test]
    fn test_write_unaligned_and_aligned() {
        let mut block = Block::new();
        block
            .push_back(
                &[0, 1, 2],
                Some("q1"),
                None,
                0,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        block
            .push_back(
                &[3, 4],
                Some("q2"),
                None,
                1,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .unwrap();
        *QUERY_ALIGNED.lock().unwrap() = vec![false, true];

        let traits = amino_traits();
        let mut unaligned = Vec::new();
        write_unaligned(
            &block,
            &mut unaligned,
            "fasta",
            &traits,
            AlignModeBlock::blastp(),
        )
        .unwrap();
        assert_eq!(std::str::from_utf8(&unaligned).unwrap(), ">q1\nARN\n");

        let mut aligned = Vec::new();
        write_aligned(
            &block,
            &mut aligned,
            "fasta",
            &traits,
            AlignModeBlock::blastp(),
        )
        .unwrap();
        assert_eq!(std::str::from_utf8(&aligned).unwrap(), ">q2\nDC\n");

        // In translated mode the emitted record is the original nucleotide
        // query, not any of its six translated protein contexts. Qualities
        // likewise belong to that source record and must survive verbatim.
        let mut translated = Block::new();
        translated
            .push_back(
                &[0, 1, 2, 3, 4],
                Some("dna query"),
                Some(b"!#$%&"),
                0,
                SequenceType::Nucleotide,
                0b11_1111,
                true,
            )
            .unwrap();
        *QUERY_ALIGNED.lock().unwrap() = vec![true];

        let mut fastq = Vec::new();
        write_aligned(
            &translated,
            &mut fastq,
            "fastq",
            &nucleotide_traits(),
            AlignModeBlock::blastx(),
        )
        .unwrap();
        assert_eq!(fastq, b"@dna query\nACGTN\n+\n!#$%&\n");

        let mut excluded = Vec::new();
        write_unaligned(
            &translated,
            &mut excluded,
            "fastq",
            &nucleotide_traits(),
            AlignModeBlock::blastx(),
        )
        .unwrap();
        assert!(excluded.is_empty());

        // `Util::Seq::format` is called with its C++ default wrap of 160.
        let mut long = Block::new();
        long.push_back(
            &vec![0; 161],
            Some("long"),
            None,
            0,
            SequenceType::AminoAcid,
            0,
            false,
        )
        .unwrap();
        *QUERY_ALIGNED.lock().unwrap() = vec![false];
        let mut wrapped = Vec::new();
        write_unaligned(
            &long,
            &mut wrapped,
            "fasta",
            &traits,
            AlignModeBlock::blastp(),
        )
        .unwrap();
        let lines: Vec<&[u8]> = wrapped.split(|&b| b == b'\n').collect();
        assert_eq!(lines[0], b">long");
        assert_eq!(lines[1].len(), 160);
        assert!(lines[1].iter().all(|&b| b == b'A'));
        assert_eq!(lines[2], b"A");
    }
}
