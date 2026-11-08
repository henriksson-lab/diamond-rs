//! Concrete native search boundary for the translated upstream regression driver.
#[cfg(not(feature = "extra"))]
use std::cell::RefCell;
use std::ffi::OsString;
use std::io;
use std::path::PathBuf;
#[cfg(not(feature = "extra"))]
use std::rc::Rc;

use crate::basic::statistics::STATISTICS;
use crate::basic::value::Letter;
use crate::commands::blastp;
use crate::data::fasta::FastaRecord;
use crate::upstream_test::test::{RegressionBackend, RegressionSequenceFile};
#[cfg(not(feature = "extra"))]
use crate::util::command_line_parser::CommandLineParser;
use crate::util::io::TempFileData;

/// The source regression database is an in-memory SequenceFile. Keep encoded
/// records in their original order, with an explicit readable cursor.
#[derive(Default)]
pub struct RegressionFile {
    records: Vec<FastaRecord>,
    cursor: usize,
}

impl RegressionSequenceFile for RegressionFile {
    fn init_write(&mut self) -> Result<(), String> {
        self.records.clear();
        self.cursor = 0;
        Ok(())
    }
    fn write_seq(&mut self, sequence: &[Letter], id: &str) -> Result<(), String> {
        self.records.push(FastaRecord {
            quality: None,
            id: id.into(),
            sequence: sequence.to_vec(),
        });
        Ok(())
    }
    fn set_seqinfo_ptr(&mut self, oid: u32) -> Result<(), String> {
        if oid as usize > self.records.len() {
            return Err("regression sequence cursor is out of bounds".into());
        }
        self.cursor = oid as usize;
        Ok(())
    }
    fn close(&mut self) -> Result<(), String> {
        self.records.clear();
        self.cursor = 0;
        Ok(())
    }
}

struct OutputPath(PathBuf, bool);
impl Drop for OutputPath {
    fn drop(&mut self) {
        if !self.1 {
            let _ = std::fs::remove_file(&self.0);
        }
    }
}

/// Concrete seed-search/extension adapter using the shared command configuration
/// and in-memory query/reference lifecycle.
#[derive(Default)]
pub struct NativeRegressionBackend {
    keep_output_files: bool,
    to_cout: bool,
    session: super::source_config::Session,
    errors: crate::util::io::IoErrorOwner,
}

impl RegressionBackend for NativeRegressionBackend {
    type SequenceFile = RegressionFile;
    fn create_sequence_file(&mut self, _name: &str) -> Result<RegressionFile, String> {
        Ok(RegressionFile::default())
    }
    fn reset_statistics(&mut self) {
        STATISTICS
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .reset();
    }
    fn run_search(
        &mut self,
        args: &[String],
        db: &mut RegressionFile,
        query: &mut RegressionFile,
    ) -> Result<Vec<u8>, String> {
        // Each Source testcase constructs a complete Config. Reuse the same
        // native factory as the production command so tuning and algorithms do
        // not silently diverge in this in-memory dependency adapter.
        // Test::run_testcase constructs Config with check_io=false because
        // both sequence-file owners are supplied directly to Search::run.
        let startup = self
            .errors
            .capture(super::source_config::prepare_without_io_in(
                args,
                &mut self.session,
            ))?;
        let invocation = args.join(" ");
        let mut config = self
            .errors
            .capture(super::dispatch::parse_native_blastp_config(
                args,
                &invocation,
            ))?;
        startup.apply_search(&mut config);
        let path = if self.to_cout {
            None
        } else {
            let temp = self
                .errors
                .capture(TempFileData::init(false).map_err(io::Error::other))?;
            let path = OutputPath(temp.name.clone(), self.keep_output_files);
            drop(temp);
            Some(path)
        };
        config.output = path.as_ref().map(|path| path.0.clone());
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads.max(0) as usize)
            .build()
            .map_err(|error| error.to_string())?;
        self.errors.capture(pool.install(|| {
            blastp::run_with_records(
                &config,
                db.records[db.cursor..].to_vec(),
                query.records[query.cursor..].to_vec(),
            )
        }))?;
        db.cursor = db.records.len();
        query.cursor = query.records.len();
        match path {
            Some(path) => self.errors.capture(std::fs::read(&path.0)),
            None => Ok(Vec::new()),
        }
    }
}

/// Standard source workflows permit the three advanced groups here, but not
/// general/search options such as --threads, --log, --out, or --query. Their
/// values are reset by the source's per-case Config construction; bootstrap
/// is copied before that reset and controls reference-hash/file generation.
pub fn parse_cli(args: &[String]) -> io::Result<bool> {
    parse_cli_os(&args.iter().map(OsString::from).collect::<Vec<_>>())
}

pub fn parse_cli_os(args: &[OsString]) -> io::Result<bool> {
    #[cfg(feature = "extra")]
    {
        let values = super::argument_registry::extra_values_os(args)?;
        return Ok(values.contains_key("bootstrap"));
    }
    #[cfg(not(feature = "extra"))]
    parse_baseline_cli_os(args)
}

#[cfg(not(feature = "extra"))]
fn parse_baseline_cli_os(args: &[OsString]) -> io::Result<bool> {
    const ALLOWED: &str = "file-buffer-size no-unlink ignore-warnings no-parse-seqids parallel-tmpdir bin ext-chunk-size no-ranking dbsize no-auto-append tantan-minMaskProb oid-output swipe-task-size anchored-swipe query-match-distance-threshold length-ratio-threshold cbs-angle linclust-banded-ext linclust-chunk-size hit-membuf algo min-orf min-query-len load-threads minichunk seed-cut freq-masking freq-sd sketch-size tile-size id2 linsearch lin-stage1 lin-combo xdrop ungapped-evalue ungapped-evalue-short short-query-ungapped-bitscore gapped-filter-evalue band shape-mask multiprocessing mp-init mp-recover mp-query-chunk culling-overlap taxon-k range-cover xml-blord-format sam-query-len stop-match-score target-indexed unaligned-targets cut-bar check-multi-target roc-file family-map family-map-query query-parallel-limit log-evalue-scale bootstrap heartbeat mp-self query-or-subject-cover";
    const BOOLEANS: &str = "no-unlink ignore-warnings no-parse-seqids no-ranking no-auto-append oid-output anchored-swipe linclust-banded-ext hit-membuf freq-masking linsearch lin-stage1 lin-combo multiprocessing mp-init mp-recover xml-blord-format sam-query-len target-indexed cut-bar check-multi-target bootstrap heartbeat mp-self";
    let allowed = ALLOWED
        .split_whitespace()
        .collect::<std::collections::HashSet<_>>();
    let booleans = BOOLEANS
        .split_whitespace()
        .collect::<std::collections::HashSet<_>>();
    let bootstrap = Rc::new(RefCell::new(false));
    let target_indexed = Rc::new(RefCell::new(false));
    let algorithm = Rc::new(RefCell::new(Vec::<u8>::new()));
    let mut parser = CommandLineParser::new();
    parser.add_command("test", "", 1);
    let permitted = parser.add_group("", vec![1], false);
    let other = parser.add_group("", vec![2], false);
    let hidden = parser.add_group("", vec![], true);
    let mut disabled = false;
    for name in crate::basic::config::registered_options(false) {
        disabled |= name == "match1";
        let group = if disabled {
            hidden
        } else if allowed.contains(name) {
            permitted
        } else {
            other
        };
        let short = super::deprecated::short_option(name);
        if group == permitted && booleans.contains(name) {
            let value = match name {
                "bootstrap" => bootstrap.clone(),
                "target-indexed" => target_indexed.clone(),
                _ => Rc::new(RefCell::new(false)),
            };
            parser.add_option(group, name, short, "", value, false, 0);
        } else if group == permitted && name == "shape-mask" {
            parser.add_option(
                group,
                name,
                short,
                "",
                Rc::new(RefCell::new(Vec::<Vec<u8>>::new())),
                vec![],
                1,
            );
        } else {
            let value = if name == "algo" {
                algorithm.clone()
            } else {
                Rc::new(RefCell::new(Vec::<u8>::new()))
            };
            parser.add_option(group, name, short, "", value, Vec::<u8>::new(), 1);
        }
    }
    let argv = args.iter().map(OsString::as_os_str).collect::<Vec<_>>();
    parser.store_os(&argv, &mut 0)?;
    if !matches!(algorithm.borrow().as_slice(), b"" | b"0" | b"1" | b"ctg") {
        let mut message = b"Invalid value for string field: ".to_vec();
        message.extend_from_slice(&algorithm.borrow());
        message.extend_from_slice(b". Permitted values: , 0, 1, ctg");
        return Err(crate::util::io::RawDiagnosticError(message).into_io());
    }
    if *target_indexed.borrow() {
        return Err(io::Error::other("--target-indexed requires -c1."));
    }
    let value = *bootstrap.borrow();
    Ok(value)
}

pub fn run(args: &[String]) -> io::Result<i32> {
    run_os(&args.iter().map(OsString::from).collect::<Vec<_>>())
}

pub fn run_os(args: &[OsString]) -> io::Result<i32> {
    let mut session = super::source_config::Session::default();
    let _startup = super::source_config::prepare_os_in(args, &mut session)?;
    run_with_session_os(args, &mut session)
}

pub(crate) fn run_with_session(
    args: &[String],
    session: &mut super::source_config::Session,
) -> io::Result<i32> {
    run_with_session_os(
        &args.iter().map(OsString::from).collect::<Vec<_>>(),
        session,
    )
}

pub(crate) fn run_with_session_os(
    args: &[OsString],
    session: &mut super::source_config::Session,
) -> io::Result<i32> {
    let bootstrap = parse_cli_os(args)?;

    #[cfg(feature = "extra")]
    let (debug_log, to_cout) = {
        let values = super::argument_registry::extra_values_os(args)?;
        (
            values.contains_key("log"),
            values
                .get("out")
                .and_then(|values| values.first())
                .is_some_and(|value| value == "stdout"),
        )
    };
    #[cfg(not(feature = "extra"))]
    let (debug_log, to_cout) = (false, false);
    let mut backend = NativeRegressionBackend {
        keep_output_files: bootstrap,
        to_cout,
        session: std::mem::take(session),
        errors: Default::default(),
    };
    // Source output always uses its color operations, even when redirected.
    let result = crate::upstream_test::test::run(
        &mut backend,
        crate::upstream_test::data::seqs(),
        crate::upstream_test::test_cases::TEST_CASES,
        crate::upstream_test::test_cases::REF_HASHES,
        crate::upstream_test::test::TestRunConfig {
            bootstrap,
            debug_log,
            to_cout,
            color: true,
        },
        &mut crate::util::io::console::stdout_endl(),
    );
    let result = backend.errors.finish(result);
    // Keep the current testcase's Source log alive for the enclosing main's
    // result reporting, then restore streams once for the complete command.
    *session = std::mem::take(&mut backend.session);
    result
}

#[cfg(test)]
mod cli_tests {
    use super::*;

    #[cfg(all(unix, not(feature = "extra")))]
    #[test]
    fn raw_baseline_test_cli_keeps_option_cells_and_diagnostic_bytes() {
        use std::os::unix::ffi::OsStringExt;
        let args = vec![
            OsString::from("diamond"),
            OsString::from("test"),
            OsString::from("--shape-mask"),
            OsString::from_vec(b"shape_\xff".to_vec()),
        ];
        assert!(!parse_cli_os(&args).unwrap());
        let args = vec![
            OsString::from("diamond"),
            OsString::from("test"),
            OsString::from("--algo"),
            OsString::from_vec(b"algo_\xff".to_vec()),
        ];
        let error = parse_cli_os(&args).unwrap_err();
        assert_eq!(
            crate::util::io::RawDiagnosticError::find(&error),
            Some(
                b"Invalid value for string field: algo_\xff. Permitted values: , 0, 1, ctg"
                    .as_slice()
            )
        );
    }

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }
    #[test]
    fn standard_test_workflow_and_bootstrap_are_native() {
        assert!(!parse_cli(&args(&["diamond", "test"])).unwrap());
        assert!(parse_cli(&args(&["diamond", "--test", "--bootstrap"])).unwrap());
        assert!(parse_cli(&args(&["diamond", "test", "--bootstrap", "--bootstrap"])).unwrap());
        assert!(!parse_cli(&args(&[
            "diamond",
            "test",
            "--algo",
            "ctg",
            "--shape-mask",
            "111",
            "101",
            "--load-threads",
            "4"
        ]))
        .unwrap());
    }
    #[test]
    fn test_parser_preserves_source_errors() {
        // EXTRA enables hidden options and removes workflow membership checks.
        // These stores have no bootstrap effect, and do not execute tests.
        for (options, expected) in [
            (
                vec!["--log"],
                "Option is not permitted for this workflow: log",
            ),
            (vec!["-p4"], "Option is not permitted for this workflow: p"),
            (
                vec!["--query", "x"],
                "Option is not permitted for this workflow: query",
            ),
            (vec!["--reverse"], "Invalid option: reverse"),
        ] {
            let mut command = args(&["diamond", "test"]);
            command.extend(options.into_iter().map(str::to_string));
            if cfg!(feature = "extra") {
                assert!(!parse_cli(&command).unwrap());
            } else {
                assert_eq!(parse_cli(&command).unwrap_err().to_string(), expected);
            }
        }
        // Cardinality remains strict in both profiles.
        for (options, expected) in [
            (
                vec!["--bootstrap", "1"],
                "Invalid parameter count for option '--bootstrap'",
            ),
            (
                vec!["--min-orf"],
                "Invalid parameter count for option '-l/--min-orf'",
            ),
        ] {
            let mut command = args(&["diamond", "test"]);
            command.extend(options.into_iter().map(str::to_string));
            assert_eq!(parse_cli(&command).unwrap_err().to_string(), expected);
        }
    }

    #[test]
    fn target_indexed_constraint_is_checked_by_source_startup() {
        let error = super::super::source_config::prepare_without_io(&args(&[
            "diamond",
            "test",
            "--target-indexed",
        ]))
        .err()
        .expect("Source Config must reject its default four chunks");
        assert_eq!(error.to_string(), "--target-indexed requires -c1.");
        #[cfg(feature = "extra")]
        {
            let startup = super::super::source_config::prepare_without_io(&args(&[
                "diamond",
                "test",
                "--target-indexed",
                "--index-chunks",
                "1",
            ]))
            .unwrap();
            assert_eq!(
                startup.config.command,
                crate::basic::config::Command::RegressionTest
            );
        }
    }
}
