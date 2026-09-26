//! Rust translation facade for `diamond/src/ffi/ffi.cpp`.
//!
//! Upstream duplicates the complete `run/main.cpp` dispatcher solely to give
//! it an `extern "C"` entry point.  Rust keeps one audited dispatcher and
//! exposes that same behavior here; ABI conversion remains in the parent
//! module's optional C++ conformance adapter.

use std::io::Write;

use crate::run::main::{RunBackend, RunFeatures};

/// Safe Rust counterpart of C++ `diamond_main`.
pub fn diamond_main<B: RunBackend, O: Write, E: Write, L: Write>(
    backend: &mut B,
    args: &[String],
    features: RunFeatures,
    out: &mut O,
    err: &mut E,
    log: &mut L,
) -> i32 {
    crate::run::main::main(backend, args, features, out, err, log)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::run::main::{MainConfig, RunCommand, RunError};

    #[derive(Clone)]
    struct Backend {
        config: Result<MainConfig, RunError>,
        init_error: Option<RunError>,
        regression_status: i32,
        messages: Vec<String>,
    }

    impl Backend {
        fn command(command: RunCommand) -> Self {
            Self {
                config: Ok(MainConfig {
                    command,
                    program_name: "diamond".to_owned(),
                    version_string: "2.1.14".to_owned(),
                    daa_file: String::new(),
                    cluster_similarity: String::new(),
                    cluster_algorithm: None,
                }),
                init_error: None,
                regression_status: 0,
                messages: Vec::new(),
            }
        }
    }

    macro_rules! ok_backend_methods {
        ($($name:ident),* $(,)?) => {$(
            fn $name(&mut self) -> Result<(), RunError> { Ok(()) }
        )* };
    }

    impl RunBackend for Backend {
        fn init_motif_table(&mut self) -> Result<(), RunError> {
            match &self.init_error {
                Some(error) => Err(error.clone()),
                None => Ok(()),
            }
        }

        fn parse_config(&mut self, _args: &[String]) -> Result<MainConfig, RunError> {
            self.config.clone()
        }

        ok_backend_methods!(
            make_db,
            search,
            view_daa,
            get_seq,
            random_seqs,
            run_masker,
            fastq2fasta,
            db_info,
            info,
            pairwise,
            benchmark,
            split,
            reverse,
            makeindex,
            find_shapes,
            hash_seqs,
            set_warning_color,
            reset_warning_color,
            list_seeds,
            cluster_realign,
            greedy_vertex_cover,
            cluster_reassign,
            recluster,
            merge_daa,
            word_count,
            cut,
            profile_recluster,
        );

        fn evaluate_cluster_similarity(&mut self, _expression: &str) -> Result<(), RunError> {
            Ok(())
        }

        fn run_cluster(&mut self, _algorithm: &str) -> Result<(), RunError> {
            Ok(())
        }

        fn regression_test(&mut self) -> Result<i32, RunError> {
            Ok(self.regression_status)
        }

        fn message(&mut self, text: &str) {
            self.messages.push(text.to_owned());
        }
    }

    type ResultTuple = (i32, Vec<u8>, Vec<u8>, Vec<u8>);

    fn through_facade(mut backend: Backend, features: RunFeatures) -> ResultTuple {
        let mut out = Vec::new();
        let mut err = Vec::new();
        let mut log = Vec::new();
        let status = diamond_main(
            &mut backend,
            &["diamond".to_owned()],
            features,
            &mut out,
            &mut err,
            &mut log,
        );
        (status, out, err, log)
    }

    fn through_dispatcher(mut backend: Backend, features: RunFeatures) -> ResultTuple {
        let mut out = Vec::new();
        let mut err = Vec::new();
        let mut log = Vec::new();
        let status = crate::run::main::main(
            &mut backend,
            &["diamond".to_owned()],
            features,
            &mut out,
            &mut err,
            &mut log,
        );
        (status, out, err, log)
    }

    fn assert_facade_parity(backend: Backend, features: RunFeatures) -> ResultTuple {
        let expected = through_dispatcher(backend.clone(), features);
        let actual = through_facade(backend, features);
        assert_eq!(actual, expected);
        actual
    }

    #[test]
    fn facade_matches_success_output_and_nonzero_command_status() {
        let version = assert_facade_parity(
            Backend::command(RunCommand::Version),
            RunFeatures::default(),
        );
        assert_eq!(
            version,
            (0, b"diamond version 2.1.14\n".to_vec(), vec![], vec![])
        );

        let mut regression = Backend::command(RunCommand::RegressionTest);
        regression.regression_status = 17;
        let regression = assert_facade_parity(regression, RunFeatures::default());
        assert_eq!(regression, (17, vec![], vec![], vec![]));

        let prepdb =
            assert_facade_parity(Backend::command(RunCommand::PrepDb), RunFeatures::default());
        assert_eq!(prepdb.0, 0);
        assert_eq!(prepdb.1, Vec::<u8>::new());
        assert_eq!(prepdb.2, b"Warning: prepdb is deprecated since v2.1.14 and no longer needed to use BLAST databases. No action was taken.\n");
    }

    #[test]
    fn facade_matches_standard_allocation_file_and_unknown_error_paths() {
        let view = assert_facade_parity(Backend::command(RunCommand::View), RunFeatures::default());
        assert_eq!(view.0, 1);
        assert_eq!(
            view.2,
            b"Error: The view command requires a DAA (option -a) input file.\n"
        );
        assert_eq!(view.3, view.2);

        let mut allocation = Backend::command(RunCommand::Help);
        allocation.init_error = Some(RunError::BadAlloc("allocator exhausted".to_owned()));
        let allocation = assert_facade_parity(allocation, RunFeatures::default());
        assert_eq!(allocation.0, 1);
        assert!(String::from_utf8(allocation.2)
            .unwrap()
            .starts_with("Failed to allocate"));
        assert_eq!(allocation.3, b"Error: allocator exhausted\n");

        let mut file = Backend::command(RunCommand::Help);
        file.init_error = Some(RunError::FileOpen);
        assert_eq!(
            assert_facade_parity(file, RunFeatures::default()),
            (1, vec![], vec![], vec![])
        );

        let mut unknown = Backend::command(RunCommand::Help);
        unknown.init_error = Some(RunError::Unknown);
        let unknown = assert_facade_parity(unknown, RunFeatures::default());
        assert_eq!(unknown.0, 1);
        assert_eq!(unknown.2, b"Exception of unknown type!\n");
        assert!(unknown.3.is_empty());
    }
}
