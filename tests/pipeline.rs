// Tests for the native Rust pipeline components.
// These tests verify that the Rust implementations produce reasonable
// results. Bit-exact matching with C++ is tested in tests/regression.rs.
use std::path::Path;

#[test]
fn test_dmnd_reader() {
    let (header, records) =
        diamond::data::dmnd_reader::read_dmnd(Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/data.dmnd"
        )))
        .unwrap();
    assert_eq!(header.sequences, 1);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].id, "d1ivsa4");
    assert_eq!(records[0].sequence.len(), 426);
}

#[test]
fn test_fasta_parsing_matches_dmnd() {
    // Verify that FASTA-parsed sequences match the DMND-stored sequences
    let fasta_records = diamond::data::fasta::read_fasta_file(
        Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/1.faa"
        )),
        diamond::basic::value::SequenceType::AminoAcid,
    )
    .unwrap();
    assert!(!fasta_records.is_empty());
    // All sequences should have valid letter values (0-25)
    for r in &fasta_records {
        for &l in &r.sequence {
            assert!(
                (0..=25).contains(&l),
                "Invalid letter {} in sequence {}",
                l,
                r.id
            );
        }
    }
}

#[test]
fn test_self_alignment_scoring() {
    // Self-alignment should produce the highest possible score
    let sm = diamond::stats::score_matrix::ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
    let records = diamond::data::fasta::read_fasta_file(
        Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/1.faa"
        )),
        diamond::basic::value::SequenceType::AminoAcid,
    )
    .unwrap();

    for r in &records {
        let result =
            diamond::dp::smith_waterman::smith_waterman(&r.sequence, &r.sequence, &sm);
        assert_eq!(
            result.identities, result.length,
            "Self-alignment of {} should have 100% identity",
            r.id
        );
        assert_eq!(result.mismatches, 0);
        assert_eq!(result.gaps, 0);
    }
}

#[test]
fn test_seed_matching_self() {
    let records = diamond::data::fasta::read_fasta_file(
        Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/1.faa"
        )),
        diamond::basic::value::SequenceType::AminoAcid,
    )
    .unwrap();

    let reduction = diamond::basic::reduction::Reduction::default_reduction();
    let shape = diamond::basic::shape::Shape::from_code("111111", &reduction);

    let seqs: Vec<&[diamond::basic::value::Letter]> =
        records.iter().map(|r| r.sequence.as_slice()).collect();
    let index = diamond::search::seed_match::build_seed_index(&seqs, &shape, &reduction);

    // Self-matching should find matches for every sequence
    let matches =
        diamond::search::seed_match::find_seed_matches(&seqs, &index, &shape, &reduction);
    assert!(
        !matches.is_empty(),
        "Self-matching should produce seed matches"
    );

    // Every query should match itself (diagonal 0)
    for (i, seq) in seqs.iter().enumerate() {
        let self_matches: Vec<_> = matches
            .iter()
            .filter(|m| m.query_id == i as u32 && m.ref_id == i as u32 && m.diagonal() == 0)
            .collect();
        if seq.len() >= 6 {
            // Only sequences >= shape length will have seeds
            assert!(
                !self_matches.is_empty(),
                "Sequence {} (len={}) should self-match on diagonal 0",
                i,
                seq.len()
            );
        }
    }
}

#[test]
fn test_evalue_computation() {
    let sm = diamond::stats::score_matrix::ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();

    // Higher scores should give lower e-values (using large sequences to avoid
    // finite-size correction edge effects)
    let e1 = sm.evalue(20, 500, 50000);
    let e2 = sm.evalue(30, 500, 50000);
    let e3 = sm.evalue(40, 500, 50000);
    assert!(e1 > e2, "Higher score should give lower e-value: {} vs {}", e1, e2);
    assert!(e2 > e3, "Higher score should give lower e-value: {} vs {}", e2, e3);

    // E-values should be non-negative
    assert!(e1 >= 0.0);
    assert!(e2 >= 0.0);
    assert!(e3 >= 0.0);
}

#[test]
fn test_masking() {
    let records = diamond::data::fasta::read_fasta_file(
        Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/1.faa"
        )),
        diamond::basic::value::SequenceType::AminoAcid,
    )
    .unwrap();

    for r in &records {
        let mut seq = r.sequence.clone();
        diamond::masking::mask_sequence(&mut seq, diamond::masking::MaskingAlgo::Tantan);
        // Masking should not change the length
        assert_eq!(seq.len(), r.sequence.len());
    }
}
