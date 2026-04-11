use crate::config::Sensitivity;

/// Sensitivity-specific search parameters.
#[derive(Debug, Clone)]
pub struct SensitivityTraits {
    pub query_indexed: bool,
    pub motif_masking: bool,
    pub freq_sd: f64,
    pub min_identities: u32,
    pub ungapped_evalue: u64,
    pub ungapped_evalue_short: u64,
    pub gapped_filter_evalue: u64,
    pub index_chunks: u32,
    pub query_bins: u32,
    pub contiguous_seed: Option<&'static str>,
    pub seed_cut: f64,
    pub block_size: f64,
}

/// Get sensitivity traits for a given level.
pub fn get_traits(sens: Sensitivity) -> SensitivityTraits {
    match sens {
        Sensitivity::Faster => SensitivityTraits {
            query_indexed: true, motif_masking: true, freq_sd: 50.0, min_identities: 11,
            ungapped_evalue: 0, ungapped_evalue_short: 0, gapped_filter_evalue: 0,
            index_chunks: 4, query_bins: 16, contiguous_seed: None, seed_cut: 0.9, block_size: 2.0,
        },
        Sensitivity::Fast => SensitivityTraits {
            query_indexed: true, motif_masking: true, freq_sd: 50.0, min_identities: 11,
            ungapped_evalue: 0, ungapped_evalue_short: 0, gapped_filter_evalue: 0,
            index_chunks: 4, query_bins: 16, contiguous_seed: None, seed_cut: 0.9, block_size: 2.0,
        },
        Sensitivity::Default => SensitivityTraits {
            query_indexed: true, motif_masking: true, freq_sd: 50.0, min_identities: 11,
            ungapped_evalue: 10000, ungapped_evalue_short: 10000, gapped_filter_evalue: 0,
            index_chunks: 4, query_bins: 16, contiguous_seed: Some("111111"), seed_cut: 0.8, block_size: 2.0,
        },
        Sensitivity::MidSensitive => SensitivityTraits {
            query_indexed: true, motif_masking: true, freq_sd: 20.0, min_identities: 11,
            ungapped_evalue: 10000, ungapped_evalue_short: 10000, gapped_filter_evalue: 0,
            index_chunks: 4, query_bins: 16, contiguous_seed: None, seed_cut: 1.0, block_size: 2.0,
        },
        Sensitivity::Sensitive => SensitivityTraits {
            query_indexed: true, motif_masking: true, freq_sd: 20.0, min_identities: 11,
            ungapped_evalue: 10000, ungapped_evalue_short: 10000, gapped_filter_evalue: 1,
            index_chunks: 4, query_bins: 16, contiguous_seed: Some("11111"), seed_cut: 1.0, block_size: 2.0,
        },
        Sensitivity::MoreSensitive => SensitivityTraits {
            query_indexed: true, motif_masking: false, freq_sd: 200.0, min_identities: 11,
            ungapped_evalue: 10000, ungapped_evalue_short: 10000, gapped_filter_evalue: 1,
            index_chunks: 4, query_bins: 16, contiguous_seed: Some("11111"), seed_cut: 1.0, block_size: 2.0,
        },
        Sensitivity::VerySensitive => SensitivityTraits {
            query_indexed: true, motif_masking: false, freq_sd: 15.0, min_identities: 9,
            ungapped_evalue: 100000, ungapped_evalue_short: 30000, gapped_filter_evalue: 1,
            index_chunks: 1, query_bins: 16, contiguous_seed: None, seed_cut: 1.0, block_size: 0.4,
        },
        Sensitivity::UltraSensitive => SensitivityTraits {
            query_indexed: true, motif_masking: false, freq_sd: 20.0, min_identities: 9,
            ungapped_evalue: 300000, ungapped_evalue_short: 30000, gapped_filter_evalue: 1,
            index_chunks: 1, query_bins: 64, contiguous_seed: None, seed_cut: 1.0, block_size: 0.4,
        },
    }
}

/// Shape codes for each sensitivity level.
pub fn get_shape_codes(sens: Sensitivity) -> &'static [&'static str] {
    match sens {
        Sensitivity::Faster | Sensitivity::Fast => &["1101110101101111"],
        Sensitivity::Default => &["111101110111", "111011010010111"],
        Sensitivity::MidSensitive => &[
            "11110110111", "1101100111101", "1110010101111", "11010101100111",
            "11101110001011", "1110100100010111", "1101000011010111", "1110011000011011",
        ],
        Sensitivity::Sensitive | Sensitivity::MoreSensitive => &[
            "1011110111", "110100100010111", "11001011111", "101110001111",
            "11011101100001", "1111010010101", "111001001001011", "10101001101011",
            "111101010011", "1111000010000111", "1100011011011", "1101010000011011",
            "1110001010101001", "110011000110011", "11011010001101", "1101001100010011",
        ],
        Sensitivity::VerySensitive => &[
            "11101111", "110110111", "111111001", "1010111011",
            "11110001011", "110100101011", "110110001101", "1010101000111",
            "1100101001011", "1101010101001", "1110010010011", "110110000010011",
            "111001000100011", "1101000100010011",
        ],
        Sensitivity::UltraSensitive => &[
            "1111111", "11101111", "110011111", "110110111",
            "111111001", "1010111011", "1011110101", "1111000111",
            "10011110011", "10101101101", "10111010101", "11001010111",
            "11001100111", "11010101101", "11110001011", "100111010011",
            "101100110101", "101110000111", "110100101011", "110110001101",
            "111000110011", "1010001011011", "1010101000111", "1010110100011",
            "1100100110011", "1100101001011", "1101001100101", "1101010101001",
            "1110001010101", "1110010010011", "10100001101101", "11000100010111",
            "11010000100111", "11010100110001", "11101000011001", "11110000001101",
            "11110100000011", "101001000001111", "110000100101011", "110010010000111",
            "110101100001001", "110110000010011", "111001000100011", "111100000100101",
            "1000110010010101", "1001000100101101", "1001000110011001", "1010001001001011",
            "1010001010010011", "1010010001010101", "1010010100010011", "1010010101001001",
            "1010100000101011", "1010100011000101", "1011000010001011", "1100010000111001",
            "1100010010001011", "1100100001001011", "1100100100100011", "1100110000001101",
            "1101000100010011", "1101000110000101", "1110000001010011", "1110100000010101",
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sensitivity_traits() {
        let t = get_traits(Sensitivity::Default);
        assert_eq!(t.index_chunks, 4);
        assert_eq!(t.contiguous_seed, Some("111111"));
    }

    #[test]
    fn test_shape_codes() {
        let codes = get_shape_codes(Sensitivity::Default);
        assert_eq!(codes.len(), 2);
        assert_eq!(codes[0], "111101110111");

        let codes = get_shape_codes(Sensitivity::UltraSensitive);
        assert_eq!(codes.len(), 64);
    }

    #[test]
    fn test_fast_shape() {
        let codes = get_shape_codes(Sensitivity::Fast);
        assert_eq!(codes.len(), 1);
        assert_eq!(codes[0], "1101110101101111");
    }
}
