//! Output-facing home for `diamond/src/output/daa/daa_record.cpp`.
//!
//! The wire decoder predates this mirrored hierarchy and remains implemented
//! in [`crate::data::daa`]. Re-exporting it here keeps a single parser while
//! making the original source-file mapping explicit.

pub use crate::data::daa::{
    copy_match_record_raw, translate_query, DaaMatch, DaaMatchIterator, DaaQueryRecord, DaaRawMatch,
};

use crate::output::format::OutputFormatSpec;

/// C++ `DAAFormat::DAAFormat` with configuration globals made explicit.
pub fn daa_format(salltitles: bool, sallseqid: bool) -> OutputFormatSpec {
    OutputFormatSpec::daa_format(salltitles, sallseqid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_transcript::{EditOperation, PackedOperation};
    use crate::data::daa::compute_flag;
    use crate::output::format::{OutputFlags, OutputFormatKind};
    use crate::util::binary_buffer::BinaryBufferIterator;
    use std::collections::HashMap;

    #[test]
    fn daa_format_preserves_title_flags() {
        let default = daa_format(false, false);
        assert_eq!(default.code, OutputFormatKind::Daa);
        assert!(default.flags.any(OutputFlags::SSEQID));
        assert!(!default.flags.any(OutputFlags::FULL_TITLES));

        let all_titles = daa_format(true, false);
        assert!(all_titles.flags.any(OutputFlags::FULL_TITLES));
        assert!(all_titles.flags.any(OutputFlags::ALL_SEQIDS));

        let all_ids = daa_format(false, true);
        assert!(!all_ids.flags.any(OutputFlags::FULL_TITLES));
        assert!(all_ids.flags.any(OutputFlags::ALL_SEQIDS));
        assert_eq!(all_ids.hsp_values, crate::dp::swipe::HspValues::TRANSCRIPT);
    }

    #[test]
    fn translate_query_maps_all_six_frames() {
        let query = vec![0, 3, 2, 0, 1, 2, 3, 0, 0];
        let context = translate_query(&query);
        assert_eq!(context.len(), 6);
        assert_eq!(context[0].len(), 3);
        assert_eq!(context[1].len(), 2);
        assert_eq!(context[2].len(), 2);
        assert_eq!(context[3].len(), 3);
        assert_eq!(context[4].len(), 2);
        assert_eq!(context[5].len(), 2);
    }

    #[test]
    fn raw_copy_remaps_subject_and_preserves_canonical_wire_record() {
        let flag = compute_flag(65_536, 256, 65_535, true);
        let mut input = Vec::new();
        input.extend_from_slice(&7u32.to_ne_bytes());
        input.push(flag);
        input.extend_from_slice(&65_536u32.to_ne_bytes());
        input.extend_from_slice(&256u16.to_ne_bytes());
        input.extend_from_slice(&65_535u16.to_ne_bytes());
        input.push(PackedOperation::from_op_count(EditOperation::Match, 9).code);
        input.push(PackedOperation::terminator().code);

        let mut map = HashMap::new();
        map.insert(7, 19);
        let mut iterator = BinaryBufferIterator::new(&input);
        let mut output = Vec::new();
        copy_match_record_raw(&mut iterator, &mut output, &map).unwrap();

        let mut expected = 19u32.to_ne_bytes().to_vec();
        expected.extend_from_slice(&input[4..]);
        assert_eq!(output, expected);
        assert!(!iterator.good());
    }

    #[test]
    fn raw_match_rejects_unterminated_transcript() {
        let flag = compute_flag(1, 2, 3, false);
        let mut input = 0u32.to_ne_bytes().to_vec();
        input.extend_from_slice(&[flag, 1, 2, 3]);
        input.push(PackedOperation::from_op_count(EditOperation::Match, 1).code);
        let error = DaaRawMatch::read(&mut BinaryBufferIterator::new(&input)).unwrap_err();
        assert!(!error.is_empty());
    }
}
