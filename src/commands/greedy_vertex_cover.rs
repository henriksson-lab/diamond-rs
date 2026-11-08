//! Native CLI adapter for the translated greedy-vertex-cover workflow.

use std::ffi::OsString;
use std::io;
use std::path::PathBuf;

use crate::tools::GreedyVertexCoverConfig;
use crate::util::command_line_parser::ReadOption;
use crate::util::log_stream::{message_stream, MessageStream};
use crate::util::options::OptionValue;

pub fn run(args: &[String]) -> io::Result<()> {
    let args: Vec<OsString> = args.iter().map(OsString::from).collect();
    run_os(&args)
}

fn option<T: ReadOption>(
    values: &super::argument_registry::OsValues,
    name: &str,
    mut default: T,
) -> io::Result<T> {
    if let Some(input) = values.get(name) {
        T::read_option_os(&mut default, input)?;
    }
    Ok(default)
}

/// The canonical registry owns every Source option, permission and cardinality.
/// Project only globals read by this workflow; unrelated accepted options retain
/// their shared startup semantics without a second, narrower option registry.
pub fn run_os(args: &[OsString]) -> io::Result<()> {
    let values = super::argument_registry::values_os(args)?;
    if args.len() == 2 {
        crate::util::io::console::source_print!(
            "{}",
            super::argument_registry::documentation("greedy-vertex-cover")
        );
    }
    if !values.contains_key("db") {
        return Err(io::Error::other("Missing parameter: --db/-d"));
    }
    let database = option(&values, "db", PathBuf::new())?;
    let edges = option(&values, "edges", PathBuf::new())?;
    let output = option(&values, "out", PathBuf::new())?;
    let centroids = option(&values, "centroid-out", PathBuf::new())?;
    let threads = option(&values, "threads", 0i32)?;
    let file_buffer = option(&values, "file-buffer-size", 67_108_864i64)?;
    let member_cover = option(&values, "member-cover", OptionValue::<f64>::new())?;
    let query_or_target_cover = option(&values, "query-or-subject-cover", 0.0f64)?;
    let edge_format = option(&values, "edge-format", Vec::<u8>::new())?;
    let depth = option(&values, "connected-component-depth", Vec::<Vec<u8>>::new())?;
    let log = values.contains_key("log");
    let quiet = values.contains_key("quiet");
    *message_stream().lock().unwrap() = MessageStream::new(log || !quiet, log);
    let optional_path = |value: PathBuf| (!value.as_os_str().is_empty()).then_some(value);
    let config = GreedyVertexCoverConfig {
        database,
        edges,
        query_or_target_cover,
        member_cover: member_cover.get_present().ok(),
        edge_format,
        symmetric: values.contains_key("symmetric"),
        strict_gvc: values.contains_key("strict-gvc"),
        no_gvc_reassign: values.contains_key("no-reassign"),
        connected_component_depth: depth,
        centroid_out: optional_path(centroids),
        output_file: optional_path(output),
    };
    let worker_threads = if threads == 0 {
        crate::util::system::hardware_concurrency() as usize
    } else {
        usize::try_from(threads).unwrap_or(0)
    };
    let result =
        crate::tools::greedy_vertex_cover::greedy_vertex_cover_from_source_files_with_threads(
            &config,
            file_buffer as usize,
            worker_threads,
        )?;
    message_stream()
        .lock()
        .unwrap_or_else(|error| error.into_inner())
        .write(format!("#Centroids: {}\n", result.centroid_count()))?;
    Ok(())
}
