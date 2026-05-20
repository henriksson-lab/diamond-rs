use crate::align::hsp::HspContext;
use crate::align::hsp::Match;
use crate::basic::statistics::{StatValue, Statistics};
use crate::basic::value::{BlockId, DictId, Letter, OId, SequenceType};
use crate::data::daa::{finish_daa_query_record, write_daa_query_record, write_daa_record_hsp};
use crate::output::format::FormatCode;
use crate::output::intermediate::{IntermediateRecord, OutputFormat};
use std::io;

/// Target/query access surface used by C++ `Search::Config` in `generate_output`.
pub trait GenerateOutputTarget {
    fn block_id2oid(&self, subject_id: BlockId) -> OId;
    fn subject_len(&self, subject_id: BlockId) -> usize;
    fn subject_seq(&self, subject_id: BlockId) -> Vec<Letter>;
    fn target_title(&self, subject_id: BlockId, all_seqids: bool) -> io::Result<String>;
    fn target_self_aln_score(&self, subject_id: BlockId) -> f64;
    fn dict_id(
        &mut self,
        current_ref_block: usize,
        subject_id: BlockId,
        output_format: &OutputFormat,
    ) -> io::Result<DictId>;
}

/// Output format virtual printing surface used by C++ `generate_output`.
pub trait GenerateOutputFormat {
    fn print_query_intro(
        &mut self,
        out: &mut Vec<u8>,
        query_block_id: BlockId,
        unaligned: bool,
    ) -> io::Result<()>;

    fn print_match(&mut self, out: &mut Vec<u8>, context: &HspContext) -> io::Result<()>;

    fn print_query_epilog(
        &mut self,
        out: &mut Vec<u8>,
        query_block_id: BlockId,
        unaligned: bool,
    ) -> io::Result<()>;
}

/// Matches C++ `generate_output(targets, query_block_id, stat, iterated)`.
#[allow(clippy::too_many_arguments)]
pub fn generate_output<T, F>(
    targets: &[Match],
    query_block_id: BlockId,
    stat: &mut Statistics,
    iterated: bool,
    report_unaligned: bool,
    current_ref_block: usize,
    output_format: &OutputFormat,
    target: &mut T,
    printer: Option<&mut F>,
    query: Vec<Vec<Letter>>,
    query_source: &[Letter],
    query_title: &str,
    query_oid: OId,
    query_self_aln_score: f64,
    input_sequence_type: SequenceType,
    all_seqids: bool,
) -> io::Result<Option<Vec<u8>>>
where
    T: GenerateOutputTarget,
    F: GenerateOutputFormat,
{
    let aligned = !targets.is_empty();
    if !aligned && !report_unaligned {
        return Ok(None);
    }

    let mut out = Vec::new();
    let mut printer = printer;
    let mut seek_pos = 0usize;
    let mut n_hsp = 0i64;

    if iterated {
        if aligned {
            seek_pos = IntermediateRecord::write_query_intro(&mut out, query_block_id);
        }
    } else if output_format.code == FormatCode::Daa {
        if aligned {
            seek_pos =
                write_daa_query_record(&mut out, query_title, query_source, input_sequence_type);
        }
    } else if aligned || report_unaligned {
        let printer = printer.as_deref_mut().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "generate_output: missing text output format",
            )
        })?;
        printer.print_query_intro(&mut out, query_block_id, !aligned)?;
    }

    for (i, target_match) in targets.iter().enumerate() {
        if target_match.hsps.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "generate_output: target with no hsps.",
            ));
        }

        let subject_id = target_match.target_block_id;
        let database_id = target.block_id2oid(subject_id);
        let subject_len = target.subject_len(subject_id);
        let target_self_aln_score = target.target_self_aln_score(subject_id);

        for (hit_hsps, hsp) in target_match.hsps.iter().enumerate() {
            if output_format.code == FormatCode::Daa {
                let dict = target.dict_id(current_ref_block, subject_id, output_format)? as u32;
                write_daa_record_hsp(&mut out, hsp, dict);
            } else if iterated {
                let dict = target.dict_id(current_ref_block, subject_id, output_format)?;
                IntermediateRecord::write(
                    &mut out,
                    hsp,
                    query_block_id,
                    dict,
                    database_id,
                    output_format,
                )?;
            } else {
                let target_title = target.target_title(subject_id, all_seqids)?;
                let context = HspContext::new(
                    hsp.clone(),
                    query_block_id,
                    query_oid,
                    query.clone(),
                    query_source.len() as i32,
                    query_title,
                    database_id,
                    subject_len as i32,
                    target_title,
                    i as u32,
                    hit_hsps as u32,
                    target.subject_seq(subject_id),
                    query_self_aln_score,
                    target_self_aln_score,
                );
                let printer = printer.as_deref_mut().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "generate_output: missing text output format",
                    )
                })?;
                printer.print_match(&mut out, &context)?;
            }

            n_hsp += 1;
        }
    }

    if iterated {
        if aligned {
            IntermediateRecord::finish_query(&mut out, seek_pos);
        }
    } else {
        stat.inc(StatValue::Matches, n_hsp);
        stat.inc(StatValue::Pairwise, targets.len() as i64);
        if aligned {
            stat.inc(StatValue::Aligned, 1);
        }
        if output_format.code == FormatCode::Daa {
            if aligned {
                finish_daa_query_record(&mut out, seek_pos);
            }
        } else if aligned || report_unaligned {
            let printer = printer.as_deref_mut().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "generate_output: missing text output format",
                )
            })?;
            printer.print_query_epilog(&mut out, query_block_id, !aligned)?;
        }
    }

    Ok(Some(out))
}

/// Matches C++ `generate_intermediate_output(targets, query_block_id, current_ref_block)`.
pub fn generate_intermediate_output<F, G>(
    targets: &[Match],
    query_block_id: BlockId,
    current_ref_block: usize,
    output_format: &OutputFormat,
    mut dict_id: F,
    mut block_id2oid: G,
) -> io::Result<Vec<u8>>
where
    F: FnMut(usize, BlockId, &OutputFormat) -> io::Result<DictId>,
    G: FnMut(BlockId) -> OId,
{
    let mut out = Vec::new();
    if targets.is_empty() {
        return Ok(out);
    }
    let seek_pos = IntermediateRecord::write_query_intro(&mut out, query_block_id);

    for target in targets {
        let block_id = target.target_block_id;
        let dict = dict_id(current_ref_block, block_id, output_format)?;
        let oid: OId = block_id2oid(block_id);
        for hsp in &target.hsps {
            IntermediateRecord::write(&mut out, hsp, query_block_id, dict, oid, output_format)?;
        }
    }

    IntermediateRecord::finish_query(&mut out, seek_pos);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::{Hsp, Match};
    use crate::basic::statistics::StatValue;
    use crate::basic::value::SequenceType;
    use crate::dp::swipe::HspValues;
    use crate::output::format::FormatCode;
    use crate::util::interval::Interval;
    use std::collections::HashMap;

    #[derive(Default)]
    struct TestTarget {
        dict_calls: Vec<(usize, BlockId)>,
        titles: HashMap<BlockId, String>,
        seqs: HashMap<BlockId, Vec<Letter>>,
    }

    impl GenerateOutputTarget for TestTarget {
        fn block_id2oid(&self, subject_id: BlockId) -> OId {
            1000 + subject_id as OId
        }

        fn subject_len(&self, subject_id: BlockId) -> usize {
            self.seqs.get(&subject_id).map(Vec::len).unwrap_or(0)
        }

        fn subject_seq(&self, subject_id: BlockId) -> Vec<Letter> {
            self.seqs.get(&subject_id).cloned().unwrap_or_default()
        }

        fn target_title(&self, subject_id: BlockId, _all_seqids: bool) -> io::Result<String> {
            Ok(self.titles.get(&subject_id).cloned().unwrap_or_default())
        }

        fn target_self_aln_score(&self, subject_id: BlockId) -> f64 {
            subject_id as f64 * 2.0
        }

        fn dict_id(
            &mut self,
            current_ref_block: usize,
            subject_id: BlockId,
            _output_format: &OutputFormat,
        ) -> io::Result<DictId> {
            self.dict_calls.push((current_ref_block, subject_id));
            Ok(700 + subject_id as DictId)
        }
    }

    #[derive(Default)]
    struct TestFormat {
        contexts: Vec<HspContext>,
    }

    impl GenerateOutputFormat for TestFormat {
        fn print_query_intro(
            &mut self,
            out: &mut Vec<u8>,
            query_block_id: BlockId,
            unaligned: bool,
        ) -> io::Result<()> {
            out.extend_from_slice(format!("intro:{query_block_id}:{unaligned}\n").as_bytes());
            Ok(())
        }

        fn print_match(&mut self, out: &mut Vec<u8>, context: &HspContext) -> io::Result<()> {
            out.extend_from_slice(
                format!(
                    "match:{}:{}:{}:{}\n",
                    context.query_id, context.subject_oid, context.hit_num, context.hsp_num
                )
                .as_bytes(),
            );
            self.contexts.push(context.clone());
            Ok(())
        }

        fn print_query_epilog(
            &mut self,
            out: &mut Vec<u8>,
            query_block_id: BlockId,
            unaligned: bool,
        ) -> io::Result<()> {
            out.extend_from_slice(format!("epilog:{query_block_id}:{unaligned}\n").as_bytes());
            Ok(())
        }
    }

    fn one_hsp_match() -> Match {
        let mut m = Match::new_extension(7, &[0, 1, 2, 3], None, 44, 0, f64::MAX);
        m.target_oid = 77;
        let mut hsp = Hsp::new();
        hsp.score = 80;
        hsp.evalue = 1.0e-9;
        hsp.frame = 0;
        hsp.query_range = Interval::new(1, 4);
        hsp.query_source_range = Interval::new(1, 4);
        hsp.subject_range = Interval::new(2, 5);
        hsp.identities = 3;
        hsp.length = 3;
        m.hsps.push(hsp);
        m
    }

    fn test_target() -> TestTarget {
        let mut target = TestTarget::default();
        target.titles.insert(7, "target7".to_string());
        target.seqs.insert(7, vec![0, 1, 2, 3, 4]);
        target
    }

    #[test]
    fn test_generate_intermediate_output_empty_is_empty_buffer() {
        let fmt = OutputFormat::new(FormatCode::Tabular, HspValues::COORDS);
        let out = generate_intermediate_output(&[], 3, 0, &fmt, |_, _, _| Ok(0), |_| 0).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn test_generate_intermediate_output_writes_query_and_records() {
        let fmt = OutputFormat::new(FormatCode::Tabular, HspValues::COORDS);
        let m = one_hsp_match();

        let out = generate_intermediate_output(
            &[m],
            5,
            2,
            &fmt,
            |block, block_id, _| {
                assert_eq!(block, 2);
                assert_eq!(block_id, 7);
                Ok(99)
            },
            |block_id| {
                assert_eq!(block_id, 7);
                7700
            },
        )
        .unwrap();

        assert_eq!(u32::from_ne_bytes(out[0..4].try_into().unwrap()), 5);
        assert!(u32::from_ne_bytes(out[4..8].try_into().unwrap()) > 0);
        let mut cursor = &out[8..];
        let record = IntermediateRecord::read(&mut cursor, &fmt).unwrap();
        assert_eq!(record.target_dict_id, 99);
        assert_eq!(record.target_oid, 0);
        assert_eq!(record.score, 80);
        assert_eq!(record.query_begin, 1);
        assert_eq!(record.query_end, 3);
        assert_eq!(record.subject_begin, 2);
        assert_eq!(record.subject_end, 5);
    }

    #[test]
    fn test_generate_intermediate_output_daa_uses_block_oid_mapping() {
        use crate::basic::packed_transcript::EditOperation;

        let fmt = OutputFormat::daa();
        let mut m = one_hsp_match();
        m.target_oid = 1;
        m.hsps[0]
            .transcript
            .push_with_count(EditOperation::Match, 3);
        m.hsps[0].transcript.push_terminator();

        let out =
            generate_intermediate_output(&[m], 5, 2, &fmt, |_, _, _| Ok(99), |_| 7700).unwrap();

        let mut cursor = &out[8..];
        let record = IntermediateRecord::read(&mut cursor, &fmt).unwrap();
        assert_eq!(record.target_dict_id, 99);
        assert_eq!(record.target_oid, 7700);
    }

    #[test]
    fn test_generate_output_empty_without_unaligned_returns_none() {
        let fmt = OutputFormat::new(FormatCode::Tabular, HspValues::COORDS);
        let mut stat = Statistics::new();
        let mut target = test_target();
        let out = generate_output::<_, TestFormat>(
            &[],
            5,
            &mut stat,
            false,
            false,
            2,
            &fmt,
            &mut target,
            None,
            vec![vec![0, 1, 2]],
            &[0, 1, 2],
            "query1",
            11,
            0.0,
            SequenceType::AminoAcid,
            false,
        )
        .unwrap();
        assert!(out.is_none());
        assert_eq!(stat.get(StatValue::Matches), 0);
    }

    #[test]
    fn test_generate_output_text_prints_context_and_stats() {
        let fmt = OutputFormat::new(FormatCode::Tabular, HspValues::COORDS);
        let mut stat = Statistics::new();
        let mut target = test_target();
        let mut printer = TestFormat::default();
        let out = generate_output(
            &[one_hsp_match()],
            5,
            &mut stat,
            false,
            false,
            2,
            &fmt,
            &mut target,
            Some(&mut printer),
            vec![vec![0, 1, 2]],
            &[0, 1, 2],
            "query1",
            11,
            99.0,
            SequenceType::AminoAcid,
            false,
        )
        .unwrap()
        .unwrap();

        assert_eq!(
            std::str::from_utf8(&out).unwrap(),
            "intro:5:false\nmatch:5:1007:0:0\nepilog:5:false\n"
        );
        assert_eq!(printer.contexts.len(), 1);
        assert_eq!(printer.contexts[0].query_title, "query1");
        assert_eq!(printer.contexts[0].target_title, "target7");
        assert_eq!(printer.contexts[0].query_self_aln_score, 99.0);
        assert_eq!(printer.contexts[0].target_self_aln_score, 14.0);
        assert_eq!(stat.get(StatValue::Matches), 1);
        assert_eq!(stat.get(StatValue::Pairwise), 1);
        assert_eq!(stat.get(StatValue::Aligned), 1);
    }

    #[test]
    fn test_generate_output_iterated_writes_intermediate_without_stats() {
        let fmt = OutputFormat::new(FormatCode::Tabular, HspValues::COORDS);
        let mut stat = Statistics::new();
        let mut target = test_target();
        let out = generate_output::<_, TestFormat>(
            &[one_hsp_match()],
            5,
            &mut stat,
            true,
            false,
            2,
            &fmt,
            &mut target,
            None,
            vec![vec![0, 1, 2]],
            &[0, 1, 2],
            "query1",
            11,
            0.0,
            SequenceType::AminoAcid,
            false,
        )
        .unwrap()
        .unwrap();

        assert_eq!(u32::from_ne_bytes(out[0..4].try_into().unwrap()), 5);
        let mut cursor = &out[8..];
        let record = IntermediateRecord::read(&mut cursor, &fmt).unwrap();
        assert_eq!(record.target_dict_id, 707);
        assert_eq!(record.score, 80);
        assert_eq!(stat.get(StatValue::Matches), 0);
        assert_eq!(target.dict_calls, vec![(2, 7)]);
    }

    #[test]
    fn test_generate_output_daa_writes_query_record_hsp_and_stats() {
        let fmt = OutputFormat::daa();
        let mut stat = Statistics::new();
        let mut target = test_target();
        let out = generate_output::<_, TestFormat>(
            &[one_hsp_match()],
            5,
            &mut stat,
            false,
            false,
            2,
            &fmt,
            &mut target,
            None,
            vec![vec![0, 1, 2]],
            &[0, 1, 2],
            "query1 title",
            11,
            0.0,
            SequenceType::AminoAcid,
            false,
        )
        .unwrap()
        .unwrap();

        let record_size = u32::from_ne_bytes(out[0..4].try_into().unwrap()) as usize;
        assert_eq!(record_size, out.len() - 4);
        assert_eq!(u32::from_ne_bytes(out[4..8].try_into().unwrap()), 3);
        assert!(out.windows(7).any(|w| w == b"query1\0"));
        assert_eq!(stat.get(StatValue::Matches), 1);
        assert_eq!(stat.get(StatValue::Pairwise), 1);
        assert_eq!(stat.get(StatValue::Aligned), 1);
        assert_eq!(target.dict_calls, vec![(2, 7)]);
    }
}
