//! Native CLI adapter for the translated greedy-vertex-cover workflow.

use std::cell::RefCell;
use std::io;
use std::path::PathBuf;
use std::rc::Rc;

use crate::tools::GreedyVertexCoverConfig;
use crate::util::command_line_parser::CommandLineParser;
use crate::util::log_stream::{message_stream, MessageStream};
use crate::util::options::OptionValue;

/// Parse with the translated upstream parser so scalar cardinality, repeated
/// vectors, short options, and numeric conversion follow the C++ rules.
pub fn run(args: &[String]) -> io::Result<()> {
    let mut parser = CommandLineParser::new();
    parser.add_command("greedy-vertex-cover", "Compute greedy vertex cover", 1);
    let group = parser.add_group("General options", vec![1], false);
    macro_rules! option {
        ($name:ident, $id:literal, $short:expr, $desc:literal, $default:expr, $min:expr) => {
            let $name = Rc::new(RefCell::new($default));
            parser.add_option(group, $id, $short, $desc, $name.clone(), $default, $min);
        };
    }
    option!(threads, "threads", 'p', "number of CPU threads", 0i32, 1);
    option!(verbose, "verbose", 'v', "verbose console output", false, 0);
    option!(log, "log", '\0', "enable debug log", false, 0);
    option!(quiet, "quiet", '\0', "disable console output", false, 0);
    option!(
        tmpdir,
        "tmpdir",
        't',
        "directory for temporary files",
        String::new(),
        1
    );
    option!(database, "db", 'd', "database file", String::new(), 1);
    option!(
        file_buffer,
        "file-buffer-size",
        '\0',
        "file buffer size in bytes",
        67_108_864i64,
        1
    );
    option!(output, "out", 'o', "output file", String::new(), 1);
    option!(
        header,
        "header",
        '\0',
        "Use header lines in tabular output format (0/simple/verbose).",
        OptionValue::<Vec<String>>::new(),
        0
    );
    option!(
        memory,
        "memory-limit",
        'M',
        "Memory limit in GB (default = 16G)",
        String::new(),
        1
    );
    option!(
        member_cover,
        "member-cover",
        '\0',
        "Minimum coverage% of the cluster member sequence (default=80.0)",
        OptionValue::<f64>::new(),
        1
    );
    option!(
        mutual_cover,
        "mutual-cover",
        '\0',
        "Minimum mutual coverage% of the cluster member and representative sequence",
        OptionValue::<f64>::new(),
        1
    );
    option!(
        depth,
        "connected-component-depth",
        '\0',
        "Depth to cluster connected components",
        Vec::<String>::new(),
        1
    );
    option!(
        no_reassign,
        "no-reassign",
        '\0',
        "Do not reassign to closest representative",
        false,
        0
    );
    option!(
        centroids,
        "centroid-out",
        '\0',
        "Output file for centroids",
        String::new(),
        1
    );
    option!(
        edges,
        "edges",
        '\0',
        "Input file for greedy vertex cover",
        String::new(),
        1
    );
    option!(
        edge_format,
        "edge-format",
        '\0',
        "Edge format for greedy vertex cover (default/triplet)",
        String::new(),
        1
    );
    option!(
        symmetric,
        "symmetric",
        '\0',
        "Edges are symmetric",
        false,
        0
    );

    let argv = args.iter().map(String::as_str).collect::<Vec<_>>();
    let mut command = 0;
    parser
        .store(&argv, &mut command)
        .map_err(io::Error::other)?;
    if args.len() == 2 {
        print!("{}", parser.print_documentation(1));
    }
    parser.require("db").map_err(io::Error::other)?;
    // These general options are parsed upstream but do not affect this
    // command's graph construction or its two plain-text output files.
    let _ = (threads, tmpdir, header, memory, mutual_cover, verbose);
    *message_stream().lock().unwrap() =
        MessageStream::new(*log.borrow() || !*quiet.borrow(), *log.borrow());
    let optional_path = |value: &str| (!value.is_empty()).then(|| PathBuf::from(value));
    let config = GreedyVertexCoverConfig {
        database: PathBuf::from(&*database.borrow()),
        edges: PathBuf::from(&*edges.borrow()),
        query_or_target_cover: 0.0,
        member_cover: member_cover.borrow().get_present().ok(),
        edge_format: edge_format.borrow().clone(),
        symmetric: *symmetric.borrow(),
        strict_gvc: false,
        no_gvc_reassign: *no_reassign.borrow(),
        connected_component_depth: depth.borrow().clone(),
        centroid_out: optional_path(&centroids.borrow()),
        output_file: optional_path(&output.borrow()),
    };
    let result = crate::tools::greedy_vertex_cover::greedy_vertex_cover_from_source_files(
        &config,
        *file_buffer.borrow() as usize,
    )?;
    message_stream()
        .lock()
        .unwrap_or_else(|error| error.into_inner())
        .write(format!("#Centroids: {}\n", result.centroid_count()))?;
    Ok(())
}
