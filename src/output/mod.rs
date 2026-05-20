pub mod edge;
pub mod format;
pub mod intermediate;
pub mod join_blocks;
pub mod output_sink;
pub mod paf;
pub mod pairwise;
pub mod sam;
pub mod target_culling;
pub mod taxon;
pub mod xml;

use crate::align::hsp::{Hsp, HspContext};

pub const INTERMEDIATE_RECORD_SEED_ONLY: u8 = 1u8 << 7;

pub fn get_length_flag(x: u32) -> u32 {
    if x <= u8::MAX as u32 {
        0
    } else if x <= u16::MAX as u32 {
        1
    } else {
        2
    }
}

pub fn get_rev_flag(frame: u32) -> u32 {
    if frame > 2 {
        1
    } else {
        0
    }
}

pub fn get_segment_flag_hsp(m: &Hsp) -> u8 {
    let rev = get_rev_flag(m.frame as u32);
    (get_length_flag(m.score as u32)
        | (get_length_flag(m.oriented_range().begin as u32) << 2)
        | (get_length_flag(m.subject_range.begin as u32) << 4)
        | (rev << 6)
        | if m.seed_only {
            INTERMEDIATE_RECORD_SEED_ONLY as u32
        } else {
            0
        }) as u8
}

pub fn get_segment_flag_context(m: &HspContext) -> u8 {
    let rev = get_rev_flag(m.frame());
    (get_length_flag(m.score())
        | (get_length_flag(m.oriented_query_range().begin as u32) << 2)
        | (get_length_flag(m.subject_range().begin as u32) << 4)
        | (rev << 6)
        | if m.seed_only() {
            INTERMEDIATE_RECORD_SEED_ONLY as u32
        } else {
            0
        }) as u8
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::interval::Interval;

    #[test]
    fn test_output_length_and_reverse_flags() {
        assert_eq!(get_length_flag(255), 0);
        assert_eq!(get_length_flag(256), 1);
        assert_eq!(get_length_flag(65_535), 1);
        assert_eq!(get_length_flag(65_536), 2);
        assert_eq!(get_rev_flag(2), 0);
        assert_eq!(get_rev_flag(3), 1);
    }

    #[test]
    fn test_output_segment_flags_hsp_and_context() {
        let mut hsp = Hsp::new();
        hsp.score = 300;
        hsp.frame = 4;
        hsp.query_source_range = Interval::new(70_000, 70_010);
        hsp.subject_range = Interval::new(100, 110);
        hsp.seed_only = true;
        assert_eq!(
            get_segment_flag_hsp(&hsp),
            1 | (2 << 2) | (0 << 4) | (1 << 6) | INTERMEDIATE_RECORD_SEED_ONLY
        );

        let ctx = HspContext::new(
            hsp.clone(),
            0,
            0,
            Vec::new(),
            0,
            "",
            0,
            0,
            "",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );
        assert_eq!(get_segment_flag_context(&ctx), get_segment_flag_hsp(&hsp));
    }
}
