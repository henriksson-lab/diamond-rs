//! Standard benchmark command, including the pinned source's failure boundary.
use std::ffi::OsString;
use std::io::{self, Write};
pub fn run(args: &[String]) -> io::Result<()> {
    run_os(&args.iter().map(OsString::from).collect::<Vec<_>>())
}

pub fn run_os(args: &[OsString]) -> io::Result<()> {
    let startup = super::source_config::prepare_os(args)?;
    run_with_startup_os(args, &startup)
}

pub(crate) fn run_with_startup(
    args: &[String],
    startup: &super::source_config::Startup,
) -> io::Result<()> {
    run_with_startup_os(
        &args.iter().map(OsString::from).collect::<Vec<_>>(),
        startup,
    )
}

pub(crate) fn run_with_startup_os(
    args: &[OsString],
    startup: &super::source_config::Startup,
) -> io::Result<()> {
    let values = super::argument_registry::values_os(args)?;
    let arch = crate::util::simd::arch();
    if values
        .get("type")
        .and_then(|items| items.first())
        .is_some_and(|kind| kind == "swipe")
    {
        #[cfg(feature = "extra")]
        if matches!(
            arch,
            crate::util::simd::Arch::Sse4_1
                | crate::util::simd::Arch::Avx2
                | crate::util::simd::Arch::Avx512
        ) {
            // Source's benchmark leaves vector cells/gap penalties undefined.
            // Execute the translated initialized recurrence; retain the real
            // rand() sequence and Source's one-million iteration count.
            std::hint::black_box(crate::tools::benchmark_swipe::swipe_cell_update_with_rng(
                1_000_000,
                &mut || crate::util::compat_rng::next_i31() as u32,
            ));
        }
        return Ok(());
    }
    let matrix = startup
        .matrix
        .as_deref()
        .ok_or_else(|| io::Error::other("Benchmark scoring matrix was not initialized."))?;
    let result = crate::tools::benchmark_native::run_with_matrix_and_controls(
        &mut crate::util::io::console::printf_stdout(),
        &mut BenchmarkMessages,
        arch,
        Default::default(),
        matrix,
        super::dispatch::parse_signed_search_option_os(args, "cbs-it-limit", 2000)?,
        super::dispatch::parse_float_search_option_os(args, "transcript-len-estimate", 1.0)?,
        super::dispatch::parse_signed_search_option_os(args, "cbs-window", 40)?,
    );
    result
}

// Source benchmark timings are MessageStream insertions, not unconditional
// cerr: preserve the active quiet/log routing and cout tie through that owner.
struct BenchmarkMessages;
impl Write for BenchmarkMessages {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        crate::util::log_stream::message_stream()
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .write_bytes(bytes)?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        crate::util::log_stream::message_stream()
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    #[test]
    #[cfg(not(feature = "extra"))]
    fn standard_build_rejects_options_before_running_workloads() {
        // Frozen diagnostics from DIAMOND 2.1.24 build 178, standard build.
        for (options, expected) in [
            (
                vec!["-p", "1"],
                "Option is not permitted for this workflow: p",
            ),
            (vec!["--type", "swipe"], "Invalid option: type"),
            (
                vec!["--matrix", "BLOSUM62"],
                "Option is not permitted for this workflow: matrix",
            ),
            (
                vec!["-q", "missing.faa"],
                "Option is not permitted for this workflow: q",
            ),
            (vec!["--help"], "Invalid option: help"),
            (vec!["--unknown"], "Invalid option: unknown"),
        ] {
            let args = ["diamond", "benchmark"]
                .into_iter()
                .chain(options)
                .map(str::to_owned)
                .collect::<Vec<_>>();
            assert_eq!(super::run(&args).unwrap_err().to_string(), expected);
        }
    }
}
