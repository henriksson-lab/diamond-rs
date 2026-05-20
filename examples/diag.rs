fn main() {
    use diamond::basic::value::{SequenceType, MASK_LETTER, SEED_MASK};
    use diamond::data::{dmnd_reader, fasta};
    use diamond::dp::{banded_cbs::banded_sw_cbs, smith_waterman::smith_waterman_cbs};
    use diamond::stats::{cbs::hauser_correction, score_matrix::ScoreMatrix};
    use std::path::Path;

    let qs = fasta::read_fasta_file(
        Path::new("/tmp/diamond-bench/test_q.faa"),
        SequenceType::AminoAcid,
    )
    .unwrap();
    let (_h, refs) = dmnd_reader::read_dmnd(Path::new("/tmp/diamond-bench/rust.dmnd")).unwrap();
    let mut q = qs[0].sequence.clone();
    let prg4 = refs.iter().find(|r| r.id.contains("PRG4")).unwrap();
    let mut t = prg4.sequence.clone();
    let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 11415371).unwrap();
    diamond::masking::mask_sequence(&mut q, diamond::masking::MaskingAlgo::Tantan);
    diamond::masking::mask_sequence(&mut t, diamond::masking::MaskingAlgo::Tantan);

    let q_hm: Vec<i8> = q
        .iter()
        .map(|&l| if l & SEED_MASK != 0 { MASK_LETTER } else { l })
        .collect();
    let cbs = hauser_correction(&q_hm, &sm);

    let f = smith_waterman_cbs(&q, &t, &sm, &cbs);
    eprintln!(
        "FULL    SW: score={} len={} id={} mm={}  q=[{},{}) s=[{},{})",
        f.score,
        f.length,
        f.identities,
        f.mismatches,
        f.query_begin,
        f.query_end,
        f.subject_begin,
        f.subject_end
    );

    for band in [60, 80, 100, 120, 200, 500, 9999] {
        let b = banded_sw_cbs(&q, &t, 812, 894, band, &sm, &cbs);
        eprintln!(
            "BAND={:5} SW: score={} len={} id={} mm={}  q=[{},{}) s=[{},{})",
            band,
            b.score,
            b.length,
            b.identities,
            b.mismatches,
            b.query_begin,
            b.query_end,
            b.subject_begin,
            b.subject_end
        );
    }
}
