//! Test-only regression driver translated from `diamond/src/test/test.cpp`.
//!
//! The source relies on process-global configuration, statistics, console
//! streams, temporary files, and the search engine.  Those dependencies are
//! explicit here so the orchestration can be tested without mutating global
//! state or recursively invoking the executable.

use std::io::Write;

use crate::basic::sequence::Sequence;
use crate::basic::value::{Letter, SequenceType, ValueTraits, AMINO_ACID_ALPHABET};
use crate::util::hash::file_hash;

use super::test_cases::TestCase;

/// Minimal writable sequence-file surface used by the upstream driver.
pub trait RegressionSequenceFile {
    fn init_write(&mut self) -> Result<(), String>;
    fn write_seq(&mut self, sequence: &[Letter], id: &str) -> Result<(), String>;
    fn set_seqinfo_ptr(&mut self, oid: u32) -> Result<(), String>;
    fn close(&mut self) -> Result<(), String>;
}

/// Safe adapter for the test-only side effects performed by the C++ driver.
pub trait RegressionBackend {
    type SequenceFile: RegressionSequenceFile;

    fn run_hit_buffer_stress_test(&mut self) -> Result<i32, String> {
        Ok(super::hit_buffer_stress::run_hit_buffer_stress_test())
    }
    fn run_queue_stress_test(&mut self) -> Result<i32, String> {
        Ok(super::queue::run_queue_stress_test())
    }
    fn create_sequence_file(&mut self, name: &str) -> Result<Self::SequenceFile, String>;
    fn reset_statistics(&mut self);
    fn run_search(
        &mut self,
        args: &[String],
        db: &mut Self::SequenceFile,
        query_file: &mut Self::SequenceFile,
    ) -> Result<Vec<u8>, String>;
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TestRunConfig {
    pub bootstrap: bool,
    pub debug_log: bool,
    pub to_cout: bool,
    /// Whether result labels should use the colors selected by C++
    /// `set_color`. Keeping this explicit makes captured output deterministic.
    pub color: bool,
}

fn passed_label(passed: bool, color: bool) -> String {
    let label = if passed { "Passed" } else { "Failed" };
    if color {
        let code = if passed { 32 } else { 31 };
        format!("\x1b[{code}m{label}\x1b[0m")
    } else {
        label.to_owned()
    }
}

/// Run one hash-verified regression case.
pub fn run_testcase<B: RegressionBackend, W: Write>(
    case: &TestCase,
    reference_hash: u64,
    db: &mut B::SequenceFile,
    query_file: &mut B::SequenceFile,
    max_width: usize,
    config: TestRunConfig,
    backend: &mut B,
    out: &mut W,
) -> Result<usize, String> {
    let mut args = Vec::new();
    args.push("diamond".to_owned());
    args.extend(
        case.command_line
            .split_ascii_whitespace()
            .map(str::to_owned),
    );
    if config.debug_log {
        args.push("--log".to_owned());
    }

    backend.reset_statistics();
    query_file.set_seqinfo_ptr(0)?;
    db.set_seqinfo_ptr(0)?;
    let output = backend.run_search(&args, db, query_file)?;

    if config.to_cout {
        out.write_all(&output).map_err(|error| error.to_string())?;
        return Ok(0);
    }

    let hash = file_hash(&output);
    if config.bootstrap {
        writeln!(out, "0x{hash:x},").map_err(|error| error.to_string())?;
        return Ok(0);
    }

    let passed = hash == reference_hash;
    writeln!(
        out,
        "{:<width$} [ {} ]",
        case.description,
        passed_label(passed, config.color),
        width = max_width
    )
    .map_err(|error| error.to_string())?;
    Ok(usize::from(passed))
}

/// Populate a writable protein sequence file from the shared test inventory.
pub fn load_seqs<F: RegressionSequenceFile>(
    file: &mut F,
    seqs: &[(String, String)],
) -> Result<(), String> {
    let traits = ValueTraits::new(AMINO_ACID_ALPHABET, 23, b"-U", SequenceType::AminoAcid);
    file.init_write()?;
    for (id, text) in seqs {
        let sequence = Sequence::from_string(text, &traits, 0)?;
        file.write_seq(&sequence, id)?;
    }
    Ok(())
}

/// Run the complete built-in regression inventory and return its process code.
pub fn run<B: RegressionBackend, W: Write>(
    backend: &mut B,
    seqs: &[(String, String)],
    test_cases: &[TestCase],
    ref_hashes: &[u64],
    config: TestRunConfig,
    out: &mut W,
) -> Result<i32, String> {
    if test_cases.len() != ref_hashes.len() {
        return Err("test case and reference hash counts differ".to_owned());
    }

    // The source intentionally ignores stress-test return codes; a transport
    // or panic-equivalent failure is still propagated by the safe adapter.
    let _ = backend.run_hit_buffer_stress_test()?;
    let _ = backend.run_queue_stress_test()?;

    // `proteins` is constructed by the source as part of dataset setup even
    // though this translation unit never writes to it.
    let _proteins = backend.create_sequence_file("test1")?;
    let mut query_file = backend.create_sequence_file("test2")?;
    let mut db = backend.create_sequence_file("test3")?;
    load_seqs(&mut query_file, seqs)?;
    load_seqs(&mut db, seqs)?;

    let max_width = test_cases
        .iter()
        .map(|case| case.description.len())
        .max()
        .unwrap_or(0);
    let mut passed = 0;
    for (case, &reference_hash) in test_cases.iter().zip(ref_hashes) {
        passed += run_testcase(
            case,
            reference_hash,
            &mut db,
            &mut query_file,
            max_width,
            config,
            backend,
            out,
        )?;
    }

    writeln!(out, "\n#Test cases passed: {passed}/{}", test_cases.len())
        .map_err(|error| error.to_string())?;
    query_file.close()?;
    db.close()?;
    Ok(i32::from(passed != test_cases.len()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Default)]
    struct File {
        initialized: bool,
        records: Vec<(String, Vec<Letter>)>,
        rewinds: usize,
        closed: bool,
    }

    impl RegressionSequenceFile for File {
        fn init_write(&mut self) -> Result<(), String> {
            self.initialized = true;
            Ok(())
        }

        fn write_seq(&mut self, sequence: &[Letter], id: &str) -> Result<(), String> {
            if !self.initialized {
                return Err("not initialized".to_owned());
            }
            self.records.push((id.to_owned(), sequence.to_vec()));
            Ok(())
        }

        fn set_seqinfo_ptr(&mut self, oid: u32) -> Result<(), String> {
            if oid != 0 {
                return Err("mock only supports rewind".to_owned());
            }
            self.rewinds += 1;
            Ok(())
        }

        fn close(&mut self) -> Result<(), String> {
            self.closed = true;
            Ok(())
        }
    }

    #[derive(Default)]
    struct Backend {
        created: Vec<String>,
        calls: Vec<Vec<String>>,
        outputs: Vec<Vec<u8>>,
        resets: usize,
        stress: Vec<&'static str>,
    }

    impl RegressionBackend for Backend {
        type SequenceFile = File;

        fn run_hit_buffer_stress_test(&mut self) -> Result<i32, String> {
            self.stress.push("hit_buffer");
            Ok(0)
        }

        fn run_queue_stress_test(&mut self) -> Result<i32, String> {
            self.stress.push("queue");
            Ok(0)
        }

        fn create_sequence_file(&mut self, name: &str) -> Result<File, String> {
            self.created.push(name.to_owned());
            Ok(File::default())
        }

        fn reset_statistics(&mut self) {
            self.resets += 1;
        }

        fn run_search(
            &mut self,
            args: &[String],
            _db: &mut File,
            _query_file: &mut File,
        ) -> Result<Vec<u8>, String> {
            self.calls.push(args.to_vec());
            Ok(self.outputs.remove(0))
        }
    }

    const CASE: TestCase = TestCase {
        description: "blastp (default)",
        command_line: "blastp -p1",
    };

    #[test]
    fn load_seqs_initializes_and_encodes_proteins() {
        let mut file = File::default();
        load_seqs(&mut file, &[("one".to_owned(), "ARu-".to_owned())]).unwrap();
        assert!(file.initialized);
        assert_eq!(file.records[0].0, "one");
        assert_eq!(file.records[0].1, vec![0, 1, 23, 23]);
    }

    #[test]
    fn run_testcase_rewinds_resets_hashes_and_builds_arguments() {
        let bytes = b"alignment\n".to_vec();
        let expected = file_hash(&bytes);
        let mut backend = Backend {
            outputs: vec![bytes],
            ..Backend::default()
        };
        let mut db = File::default();
        let mut query = File::default();
        let mut out = Vec::new();
        let passed = run_testcase(
            &CASE,
            expected,
            &mut db,
            &mut query,
            CASE.description.len(),
            TestRunConfig {
                debug_log: true,
                ..TestRunConfig::default()
            },
            &mut backend,
            &mut out,
        )
        .unwrap();
        assert_eq!(passed, 1);
        assert_eq!(backend.resets, 1);
        assert_eq!(db.rewinds, 1);
        assert_eq!(query.rewinds, 1);
        assert_eq!(backend.calls[0], ["diamond", "blastp", "-p1", "--log"]);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "blastp (default) [ Passed ]\n"
        );
    }

    #[test]
    fn bootstrap_and_stdout_modes_skip_comparison() {
        let bytes = b"result".to_vec();
        let mut backend = Backend {
            outputs: vec![bytes.clone(), bytes.clone()],
            ..Backend::default()
        };
        let mut db = File::default();
        let mut query = File::default();
        let mut out = Vec::new();
        assert_eq!(
            run_testcase(
                &CASE,
                0,
                &mut db,
                &mut query,
                0,
                TestRunConfig {
                    bootstrap: true,
                    ..TestRunConfig::default()
                },
                &mut backend,
                &mut out,
            )
            .unwrap(),
            0
        );
        assert_eq!(
            String::from_utf8(out).unwrap(),
            format!("0x{:x},\n", file_hash(&bytes))
        );

        let mut direct = Vec::new();
        assert_eq!(
            run_testcase(
                &CASE,
                0,
                &mut db,
                &mut query,
                0,
                TestRunConfig {
                    to_cout: true,
                    ..TestRunConfig::default()
                },
                &mut backend,
                &mut direct,
            )
            .unwrap(),
            0
        );
        assert_eq!(direct, bytes);
    }

    #[test]
    fn run_preserves_stress_setup_order_and_summary_status() {
        let good = b"good".to_vec();
        let bad = b"bad".to_vec();
        let cases = [
            CASE,
            TestCase {
                description: "blastp (other)",
                command_line: "blastp -p4",
            },
        ];
        let hashes = [file_hash(&good), 0];
        let mut backend = Backend {
            outputs: vec![good, bad],
            ..Backend::default()
        };
        let mut out = Vec::new();
        let status = run(
            &mut backend,
            &[("id".to_owned(), "ARND".to_owned())],
            &cases,
            &hashes,
            TestRunConfig::default(),
            &mut out,
        )
        .unwrap();
        assert_eq!(status, 1);
        assert_eq!(backend.stress, ["hit_buffer", "queue"]);
        assert_eq!(backend.created, ["test1", "test2", "test3"]);
        assert_eq!(backend.resets, 2);
        assert!(String::from_utf8(out)
            .unwrap()
            .ends_with("\n#Test cases passed: 1/2\n"));
    }

    #[test]
    fn source_function_inventory_remains_audited() {
        let source = String::from_utf8_lossy(include_bytes!("../../diamond/src/test/test.cpp"));
        for signature in ["run_testcase(", "load_seqs(", "int run()"] {
            assert!(
                source.contains(signature),
                "missing source function {signature}"
            );
        }
        assert_eq!(source.matches("static size_t run_testcase(").count(), 1);
        assert_eq!(source.matches("static void load_seqs(").count(), 1);
        assert_eq!(source.matches("int run()").count(), 1);
    }
}
