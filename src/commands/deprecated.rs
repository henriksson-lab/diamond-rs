//! Standard upstream workflows with no permitted command-line options.
use std::cell::RefCell;
use std::io;
use std::rc::Rc;

use crate::util::command_line_parser::CommandLineParser;

/// Preserve registration/disabled/workflow errors before entering `roc` or
/// `prepdb`. Neither workflow is a member of any standard option group.
pub fn validate(args: &[String], command: &str) -> io::Result<()> {
    let mut parser = CommandLineParser::new();
    parser.add_command(command, "", 1);
    let other = parser.add_group("", vec![2], false);
    let hidden = parser.add_group("", vec![], true);
    let mut disabled = false;
    for name in crate::basic::config::registered_options(false) {
        disabled |= name == "match1";
        let short = short_option(name);
        parser.add_option(
            if disabled { hidden } else { other },
            name,
            short,
            "",
            Rc::new(RefCell::new(String::new())),
            String::new(),
            1,
        );
    }
    let argv = args.iter().map(String::as_str).collect::<Vec<_>>();
    parser.store(&argv, &mut 0).map_err(io::Error::other)
}

/// Source standard short-option aliases, shared by workflow validators.
pub(crate) fn short_option(name: &str) -> char {
    match name {
        "threads" => 'p',
        "verbose" => 'v',
        "tmpdir" => 't',
        "db" => 'd',
        "out" => 'o',
        "evalue" => 'e',
        "max-target-seqs" => 'k',
        "shapes" => 's',
        "query" => 'q',
        "global-ranking" => 'g',
        "block-size" => 'b',
        "index-chunks" => 'c',
        "frameshift" => 'F',
        "outfmt" => 'f',
        "memory-limit" => 'M',
        "min-orf" => 'l',
        "xdrop" => 'x',
        "daa" => 'a',
        "window" | "seed-weight" => 'w',
        "gapped-xdrop" => 'X',
        _ => '\0',
    }
}
