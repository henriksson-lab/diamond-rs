//! Frequency-ratio matrix lookup from NCBI BLAST.
//!
//! The eight 28x28 numeric payloads are shared with [`crate::stats::matrix_data`]
//! so there is one authoritative Rust copy. This module preserves the vendor
//! aliases, scaled BLOSUM62 variants, bit-scale metadata, and owned mutable
//! result returned by `_PSIMatrixFrequencyRatiosNew`.

use crate::stats::matrix_data;

pub const BLASTAA_SIZE: usize = 28;
pub const BLOSUM62_20A_SCALE_MULTIPLIER: f64 = 0.9666;
pub const BLOSUM62_20B_SCALE_MULTIPLIER: f64 = 0.9344;

/// Owned counterpart of the vendor `SFreqRatios` allocation.
#[derive(Debug, Clone, PartialEq)]
pub struct SFreqRatios {
    pub data: Box<[[f64; BLASTAA_SIZE]; BLASTAA_SIZE]>,
    pub bit_scale_factor: i32,
}

type Matrix = [[f64; BLASTAA_SIZE]; BLASTAA_SIZE];

/// Allocate an owned frequency-ratio matrix for a vendor matrix name.
///
/// A Rust `&str` represents the C routine's asserted non-null name. Unknown
/// names return `None`, matching its null return.
pub fn psi_matrix_frequency_ratios_new(matrix_name: &str) -> Option<SFreqRatios> {
    let (source, bit_scale_factor, multiplier): (&Matrix, i32, f64) = if matrix_name
        .eq_ignore_ascii_case("BLOSUM62")
        || matrix_name.eq_ignore_ascii_case("BLOSUM62_20")
    {
        (&matrix_data::BLOSUM62_FREQ_RATIOS, 2, 1.0)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM62_20A") {
        (
            &matrix_data::BLOSUM62_FREQ_RATIOS,
            2,
            BLOSUM62_20A_SCALE_MULTIPLIER,
        )
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM62_20B") {
        (
            &matrix_data::BLOSUM62_FREQ_RATIOS,
            2,
            BLOSUM62_20B_SCALE_MULTIPLIER,
        )
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM45") {
        (&matrix_data::BLOSUM45_FREQ_RATIOS, 3, 1.0)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM80") {
        (&matrix_data::BLOSUM80_FREQ_RATIOS, 2, 1.0)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM50") {
        (&matrix_data::BLOSUM50_FREQ_RATIOS, 2, 1.0)
    } else if matrix_name.eq_ignore_ascii_case("BLOSUM90") {
        (&matrix_data::BLOSUM90_FREQ_RATIOS, 2, 1.0)
    } else if matrix_name.eq_ignore_ascii_case("PAM30") {
        (&matrix_data::PAM30_FREQ_RATIOS, 2, 1.0)
    } else if matrix_name.eq_ignore_ascii_case("PAM70") {
        (&matrix_data::PAM70_FREQ_RATIOS, 2, 1.0)
    } else if matrix_name.eq_ignore_ascii_case("PAM250") {
        (&matrix_data::PAM250_FREQ_RATIOS, 2, 1.0)
    } else {
        return None;
    };
    let mut data = Box::new(*source);

    if multiplier != 1.0 {
        for row in data.iter_mut() {
            for value in row {
                // Keep the operand order used by the C initializer loop.
                *value = multiplier * *value;
            }
        }
    }

    Some(SFreqRatios {
        data,
        bit_scale_factor,
    })
}

/// Consume a frequency-ratio allocation and return null, as the C free does.
pub fn psi_matrix_frequency_ratios_free(_freq_ratios: Option<SFreqRatios>) -> Option<SFreqRatios> {
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_upstream_table(name: &str) -> Vec<f64> {
        let source = include_str!("../../../diamond/src/lib/blast/matrix_freq_ratios.c");
        let declaration = format!("static const double {name}");
        let begin = source
            .find(&declaration)
            .expect("upstream table declaration");
        let initializer = &source[begin..];
        let initializer = &initializer[initializer.find('=').expect("table equals") + 1..];
        let initializer = &initializer[..initializer.find(';').expect("table semicolon")];

        initializer
            .split(|c: char| c == '{' || c == '}' || c == ',' || c.is_whitespace())
            .filter(|field| !field.is_empty())
            .map(|field| field.parse::<f64>().expect("floating-point literal"))
            .collect()
    }

    fn assert_full_table(name: &str, rust: &Matrix) {
        let upstream = parse_upstream_table(name);
        assert_eq!(upstream.len(), BLASTAA_SIZE * BLASTAA_SIZE);
        for (index, (actual, expected)) in rust.iter().flatten().zip(upstream).enumerate() {
            assert_eq!(
                actual.to_bits(),
                expected.to_bits(),
                "{name} differs at flat index {index}"
            );
        }
    }

    #[test]
    fn all_eight_payloads_are_bit_exact_to_c_source() {
        assert_full_table("BLOSUM45_FREQRATIOS", &matrix_data::BLOSUM45_FREQ_RATIOS);
        assert_full_table("BLOSUM50_FREQRATIOS", &matrix_data::BLOSUM50_FREQ_RATIOS);
        assert_full_table("BLOSUM62_FREQRATIOS", &matrix_data::BLOSUM62_FREQ_RATIOS);
        assert_full_table("BLOSUM80_FREQRATIOS", &matrix_data::BLOSUM80_FREQ_RATIOS);
        assert_full_table("BLOSUM90_FREQRATIOS", &matrix_data::BLOSUM90_FREQ_RATIOS);
        assert_full_table("PAM30_FREQRATIOS", &matrix_data::PAM30_FREQ_RATIOS);
        assert_full_table("PAM70_FREQRATIOS", &matrix_data::PAM70_FREQ_RATIOS);
        assert_full_table("PAM250_FREQRATIOS", &matrix_data::PAM250_FREQ_RATIOS);
    }

    #[test]
    fn lookup_preserves_alias_case_bit_scales_and_owned_copy() {
        let base = psi_matrix_frequency_ratios_new("blosum62").unwrap();
        let alias = psi_matrix_frequency_ratios_new("BLOSUM62_20").unwrap();
        assert_eq!(base.bit_scale_factor, 2);
        assert_eq!(base, alias);
        assert_eq!(*base.data, matrix_data::BLOSUM62_FREQ_RATIOS);

        let blosum45 = psi_matrix_frequency_ratios_new("bLoSuM45").unwrap();
        assert_eq!(blosum45.bit_scale_factor, 3);
        assert_eq!(*blosum45.data, matrix_data::BLOSUM45_FREQ_RATIOS);

        for name in [
            "BLOSUM50", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250",
        ] {
            assert_eq!(
                psi_matrix_frequency_ratios_new(name)
                    .unwrap()
                    .bit_scale_factor,
                2
            );
        }

        let mut first = psi_matrix_frequency_ratios_new("BLOSUM62").unwrap();
        first.data[1][1] = -1.0;
        let second = psi_matrix_frequency_ratios_new("BLOSUM62").unwrap();
        assert_ne!(first.data[1][1], second.data[1][1]);
    }

    #[test]
    fn scaled_variants_match_every_c_loop_operation() {
        for (name, multiplier) in [
            ("BLOSUM62_20A", BLOSUM62_20A_SCALE_MULTIPLIER),
            ("BLOSUM62_20B", BLOSUM62_20B_SCALE_MULTIPLIER),
        ] {
            let actual = psi_matrix_frequency_ratios_new(name).unwrap();
            assert_eq!(actual.bit_scale_factor, 2);
            for (value, source) in actual
                .data
                .iter()
                .flatten()
                .zip(matrix_data::BLOSUM62_FREQ_RATIOS.iter().flatten())
            {
                assert_eq!(value.to_bits(), (multiplier * *source).to_bits());
            }
        }
    }

    #[test]
    fn unknown_and_free_preserve_null_semantics() {
        assert!(psi_matrix_frequency_ratios_new("BLOSUM100").is_none());
        assert!(psi_matrix_frequency_ratios_free(None).is_none());
        assert!(
            psi_matrix_frequency_ratios_free(psi_matrix_frequency_ratios_new("PAM250")).is_none()
        );
    }
}
