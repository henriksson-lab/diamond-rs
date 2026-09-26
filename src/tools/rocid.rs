//! Family-stratified identity ROC tool from `diamond/src/tools/rocid.cpp`.
//!
//! The upstream implementation stores its state and streams in globals. This
//! module uses owned state and explicit readers/writers while preserving its
//! ordered-map look-ahead behavior.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

use crate::util::string::{CharDelimiter, Tokenizer};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Assoc {
    pub id: u8,
    pub fam_idx: u8,
}

#[derive(Debug, Clone, Default)]
pub struct RocIdState {
    pub query_aln: String,
    pub query_mapped: String,
    pub totals: Vec<[i32; 10]>,
    pub counts: Vec<[i32; 10]>,
    pub fam2idx: BTreeMap<String, u8>,
    pub acc2id: HashMap<String, Vec<Assoc>>,
    pub unmapped_query: usize,
    pub total_unmapped: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RocIdConfig {
    pub query_file: PathBuf,
    pub family_map: PathBuf,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RocIdStats {
    pub queries: usize,
    pub hits: usize,
    pub unmapped_hits: usize,
    pub unmapped_queries: usize,
}

/// Line reader with the one-line putback behavior used by `TextInputFile`.
pub struct LineInput<R> {
    reader: R,
    line: String,
    putback: bool,
    eof: bool,
}

impl<R: BufRead> LineInput<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line: String::new(),
            putback: false,
            eof: false,
        }
    }

    fn getline(&mut self) -> Result<bool, String> {
        if self.putback {
            self.putback = false;
            return Ok(!self.eof);
        }
        let mut bytes = Vec::new();
        let read = self
            .reader
            .read_until(b'\n', &mut bytes)
            .map_err(|error| error.to_string())?;
        let terminated = bytes.last() == Some(&b'\n');
        if terminated {
            bytes.pop();
        }
        if bytes.last() == Some(&b'\r') {
            bytes.pop();
        }
        self.eof = read == 0 || !terminated;
        self.line = String::from_utf8(bytes).map_err(|error| error.to_string())?;
        Ok(!self.eof)
    }

    fn putback_line(&mut self) {
        self.putback = true;
    }
}

fn parse_map_line(line: &str) -> Result<(String, String, f32, String), String> {
    let mut tokens = Tokenizer::new(line, CharDelimiter::new('\t'));
    Ok((
        tokens.read_string().map_err(|error| error.to_string())?,
        tokens.read_string().map_err(|error| error.to_string())?,
        tokens.read_f32().map_err(|error| error.to_string())?,
        tokens.read_string().map_err(|error| error.to_string())?,
    ))
}

fn parse_alignment_line(line: &str) -> Result<(String, String), String> {
    let mut tokens = Tokenizer::new(line, CharDelimiter::new('\t'));
    Ok((
        tokens.read_string().map_err(|error| error.to_string())?,
        tokens.read_string().map_err(|error| error.to_string())?,
    ))
}

fn format_default_float(value: f64) -> String {
    if value == 0.0 {
        return "0".to_owned();
    }
    let exponent = value.abs().log10().floor() as i32;
    if !(-4..6).contains(&exponent) {
        let raw = format!("{value:.5e}");
        let (mantissa, exponent) = raw.split_once('e').expect("scientific format has e");
        let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
        let exponent: i32 = exponent.parse().expect("scientific exponent is numeric");
        format!("{mantissa}e{exponent:+03}")
    } else {
        let decimals = (5 - exponent).max(0) as usize;
        let formatted = format!("{value:.decimals$}");
        if decimals == 0 {
            formatted
        } else {
            formatted
                .trim_end_matches('0')
                .trim_end_matches('.')
                .to_owned()
        }
    }
}

/// Load the next query group from the ordered family map.
pub fn fetch_map<R: BufRead>(
    map_in: &mut LineInput<R>,
    query: &str,
    state: &mut RocIdState,
) -> Result<bool, String> {
    state.acc2id.clear();
    state.fam2idx.clear();
    state.counts.clear();
    state.totals.clear();
    state.query_mapped.clear();
    let mut next_query = String::new();

    while map_in.getline()? {
        let (q, target, identity, family) = parse_map_line(&map_in.line)?;
        if next_query.is_empty() {
            next_query.clone_from(&q);
            state.query_mapped.clone_from(&q);
            if next_query.as_str() > query {
                return Ok(true);
            }
        }
        if q != next_query {
            map_in.putback_line();
            return Ok(next_query == query);
        }

        let fam_idx = if let Some(&idx) = state.fam2idx.get(&family) {
            idx
        } else {
            let idx = u8::try_from(state.fam2idx.len())
                .map_err(|_| "more than 256 families in one query".to_owned())?;
            state.fam2idx.insert(family, idx);
            state.totals.push([0; 10]);
            state.counts.push([0; 10]);
            idx
        };
        if !identity.is_finite() || identity < 0.0 {
            return Err("identity must be a finite non-negative value".to_owned());
        }
        let bin = (((identity * 100.0) as i32) / 10).min(9) as usize;
        state.acc2id.entry(target).or_default().push(Assoc {
            id: bin as u8,
            fam_idx,
        });
        state.totals[fam_idx as usize][bin] = state.totals[fam_idx as usize][bin]
            .checked_add(1)
            .ok_or_else(|| "family-map count overflow".to_owned())?;
    }
    Ok(next_query == query || (next_query.is_empty() && map_in.eof))
}

/// Emit the current query's ten averaged identity-bin recovery fractions.
pub fn print<W: Write, M: Write>(
    state: &mut RocIdState,
    out: &mut W,
    message: &mut M,
) -> Result<(), String> {
    if state.unmapped_query != 0
        || state.query_mapped > state.query_aln
        || state.query_mapped.is_empty()
    {
        writeln!(message, "Unmapped query: {}", state.query_aln)
            .map_err(|error| error.to_string())?;
        state.total_unmapped = state
            .total_unmapped
            .checked_add(1)
            .ok_or_else(|| "unmapped query count overflow".to_owned())?;
        return Ok(());
    }

    write!(out, "{}", state.query_mapped).map_err(|error| error.to_string())?;
    for bin in 0..10 {
        let mut sum = 0.0;
        let mut families = 0.0;
        for fam_idx in 0..state.fam2idx.len() {
            if state.totals[fam_idx][bin] > 0 {
                sum += state.counts[fam_idx][bin] as f64 / state.totals[fam_idx][bin] as f64;
                families += 1.0;
            }
        }
        let value = if families > 0.0 { sum / families } else { -1.0 };
        write!(out, "\t{}", format_default_float(value)).map_err(|error| error.to_string())?;
    }
    writeln!(out).map_err(|error| error.to_string())
}

/// Reader-based implementation used by the path-configured command and tests.
pub fn roc_id_with_io<A: BufRead, F: BufRead, W: Write, M: Write>(
    alignment: A,
    family_map: F,
    out: &mut W,
    message: &mut M,
) -> Result<RocIdStats, String> {
    let mut input = LineInput::new(alignment);
    let mut map_input = LineInput::new(family_map);
    let mut state = RocIdState::default();
    let mut stats = RocIdStats::default();

    while input.getline()? {
        if input.line.is_empty() {
            break;
        }
        let (query, target) = parse_alignment_line(&input.line)?;
        stats.hits = stats.hits.checked_add(1).ok_or("hit count overflow")?;
        if query != state.query_aln {
            print(&mut state, out, message)?;
            state.unmapped_query = 0;
            state.query_aln = query.clone();
            while !fetch_map(&mut map_input, &query, &mut state)? {
                print(&mut state, out, message)?;
            }
            stats.queries = stats.queries.checked_add(1).ok_or("query count overflow")?;
            if stats.queries % 1000 == 0 {
                writeln!(
                    message,
                    "{} {} {}",
                    stats.queries, stats.hits, stats.unmapped_hits
                )
                .map_err(|error| error.to_string())?;
            }
        }

        let Some(associations) = state.acc2id.get(&target) else {
            state.unmapped_query = state
                .unmapped_query
                .checked_add(1)
                .ok_or("query-unmapped count overflow")?;
            stats.unmapped_hits = stats
                .unmapped_hits
                .checked_add(1)
                .ok_or("unmapped hit count overflow")?;
            continue;
        };
        for association in associations {
            let count = &mut state.counts[association.fam_idx as usize][association.id as usize];
            *count = count.checked_add(1).ok_or("matched count overflow")?;
        }
    }
    print(&mut state, out, message)?;
    stats.unmapped_queries = state.total_unmapped;
    writeln!(message, "Queries = {}", stats.queries).map_err(|error| error.to_string())?;
    writeln!(message, "Unmapped = {}", stats.unmapped_queries)
        .map_err(|error| error.to_string())?;
    Ok(stats)
}

/// Run `roc-id` using explicit input paths and output destinations.
pub fn roc_id<W: Write, M: Write>(
    config: &RocIdConfig,
    out: &mut W,
    message: &mut M,
) -> Result<RocIdStats, String> {
    let alignment = File::open(&config.query_file).map_err(|error| error.to_string())?;
    let family_map = File::open(&config.family_map).map_err(|error| error.to_string())?;
    roc_id_with_io(
        BufReader::new(alignment),
        BufReader::new(family_map),
        out,
        message,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn fetch_map_groups_families_targets_and_identity_bins() {
        let data = b"q1\ta\t0.05\tf1\nq1\ta\t0.95\tf2\nq1\tb\t1.00\tf1\nq2\tc\t0.40\tf3\n";
        let mut input = LineInput::new(Cursor::new(data));
        let mut state = RocIdState::default();
        assert!(fetch_map(&mut input, "q1", &mut state).unwrap());
        assert_eq!(state.fam2idx.len(), 2);
        assert_eq!(state.acc2id["a"].len(), 2);
        assert_eq!(state.totals[0][0], 1);
        assert_eq!(state.totals[0][9], 1);
        assert!(fetch_map(&mut input, "q2", &mut state).unwrap());
        assert_eq!(state.totals[0][4], 1);
    }

    #[test]
    fn print_averages_only_families_present_in_each_bin() {
        let mut state = RocIdState {
            query_aln: "q".into(),
            query_mapped: "q".into(),
            totals: vec![[0; 10], [0; 10]],
            counts: vec![[0; 10], [0; 10]],
            fam2idx: BTreeMap::from([("a".into(), 0), ("b".into(), 1)]),
            ..RocIdState::default()
        };
        state.totals[0][3] = 4;
        state.counts[0][3] = 2;
        state.totals[1][3] = 2;
        state.counts[1][3] = 2;
        let (mut out, mut messages) = (Vec::new(), Vec::new());
        print(&mut state, &mut out, &mut messages).unwrap();
        let fields: Vec<_> = std::str::from_utf8(&out)
            .unwrap()
            .trim()
            .split('\t')
            .collect();
        assert_eq!(fields[4], "0.75");
        assert!(messages.is_empty());
    }

    #[test]
    fn print_uses_cpp_default_six_significant_digit_formatting() {
        assert_eq!(format_default_float(1.0 / 3.0), "0.333333");
        assert_eq!(format_default_float(1.0), "1");
        assert_eq!(format_default_float(0.00001), "1e-05");
    }

    #[test]
    fn roc_id_preserves_reporting_and_duplicate_associations() {
        let alignments = b"q1\ta\nq1\tb\nq2\tc\n";
        let map = b"q1\ta\t0.25\tf1\nq1\ta\t0.25\tf2\nq1\tb\t0.25\tf1\nq2\tc\t0.85\tf3\n";
        let (mut out, mut messages) = (Vec::new(), Vec::new());
        let stats = roc_id_with_io(
            Cursor::new(alignments),
            Cursor::new(map),
            &mut out,
            &mut messages,
        )
        .unwrap();
        assert_eq!((stats.queries, stats.hits, stats.unmapped_hits), (2, 3, 0));
        assert_eq!(stats.unmapped_queries, 1); // Initial source `print()` call.
        let rows: Vec<_> = std::str::from_utf8(&out).unwrap().lines().collect();
        assert_eq!(rows[0].split('\t').nth(3), Some("1"));
        assert_eq!(rows[1].split('\t').nth(9), Some("1"));
        assert_eq!(
            String::from_utf8(messages).unwrap(),
            "Unmapped query: \nQueries = 2\nUnmapped = 1\n"
        );
    }

    #[test]
    fn unknown_target_marks_whole_query_unmapped() {
        let (mut out, mut messages) = (Vec::new(), Vec::new());
        let stats = roc_id_with_io(
            Cursor::new(b"q1\ta\nq1\tmissing\n"),
            Cursor::new(b"q1\ta\t0.50\tf1\n"),
            &mut out,
            &mut messages,
        )
        .unwrap();
        assert!(out.is_empty());
        assert_eq!((stats.unmapped_hits, stats.unmapped_queries), (1, 2));
    }

    #[test]
    fn unterminated_final_lines_match_source_eof_semantics() {
        let (mut out, mut messages) = (Vec::new(), Vec::new());
        let stats = roc_id_with_io(
            Cursor::new(b"q1\ta"),
            Cursor::new(b"q1\ta\t0.5\tf\n"),
            &mut out,
            &mut messages,
        )
        .unwrap();
        assert_eq!(
            (stats.hits, stats.queries, stats.unmapped_queries),
            (0, 0, 1)
        );
    }

    #[test]
    fn source_function_inventory_remains_audited() {
        let source = String::from_utf8_lossy(include_bytes!("../../diamond/src/tools/rocid.cpp"));
        for signature in [
            "static bool fetch_map(",
            "static void print()",
            "void roc_id()",
        ] {
            assert!(
                source.contains(signature),
                "missing source function {signature}"
            );
        }
    }
}
