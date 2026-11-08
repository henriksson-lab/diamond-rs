//! Source makeidx argv ownership and the concrete persistent-index backend.
use super::source_config::Startup;
use crate::util::command_line_parser::ReadOption;
use std::ffi::OsString;
use std::io;
use std::sync::Arc;

pub(crate) fn run(args: &[String]) -> io::Result<()> {
    run_os(&args.iter().map(OsString::from).collect::<Vec<_>>())
}
pub(crate) fn run_os(args: &[OsString]) -> io::Result<()> {
    run_with_startup_os(args, None)
}
pub(crate) fn run_with_startup_os(
    args: &[OsString],
    startup: Option<Arc<Startup>>,
) -> io::Result<()> {
    let startup = match startup {
        Some(startup) => startup,
        None => super::source_config::prepare_os(args)?,
    };
    let values = super::argument_registry::values_os(args)?;
    let mut config = crate::data::index::MakeIndexConfig::new(
        startup.config.database.clone(),
        startup.config.sensitivity,
    );
    config.threads = startup.config.threads.max(0) as usize;
    if let Some(items) = values.get("shape-mask") {
        config.shape_mask = items
            .iter()
            .map(|s| {
                crate::util::io::source_os_bytes(s)
                    .map(|s| s.into_owned())
                    .map_err(io::Error::other)
            })
            .collect::<io::Result<_>>()?;
    }
    if let Some(items) = values.get("shapes") {
        u32::read_option_os(&mut config.shapes, items)?;
    }
    if let Some(items) = values.get("file-buffer-size") {
        let mut size = 67_108_864i64;
        i64::read_option_os(&mut size, items)?;
        config.file_buffer_size = size as usize;
    }
    crate::data::index::make_index(&config).map(|_| ())
}
