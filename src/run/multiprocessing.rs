//! Concrete file-stack lifecycle for `run/double_indexed.cpp` multiprocessing.
//! Work records and named block files remain on disk until their actual join
//! succeeds, permitting a later `--mp-recover` process to restore in-flight work.

use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use super::double_indexed::{
    get_ref_part_file_name, STACK_ALIGN_DONE, STACK_ALIGN_TODO, STACK_ALIGN_WIP, STACK_JOIN_DONE,
    STACK_JOIN_TODO, STACK_JOIN_WIP,
};
use crate::data::dmnd::{to_chunk, Chunk};
use crate::util::parallel::{filestack::FileStack, parallelizer::Parallelizer};

pub struct MultiprocessingQueues {
    parallelizer: Parallelizer,
    directory: PathBuf,
}
pub struct QueryQueues {
    query: usize,
    todo: Arc<FileStack>,
    wip: Arc<FileStack>,
    done: Arc<FileStack>,
    join_todo: Arc<FileStack>,
    join_wip: Arc<FileStack>,
    join_done: Arc<FileStack>,
}

impl MultiprocessingQueues {
    pub fn new(directory: &Path) -> io::Result<Self> {
        let mut parallelizer = Parallelizer::new();
        parallelizer.init(directory).map_err(io::Error::other)?;
        Ok(Self {
            parallelizer,
            directory: directory.into(),
        })
    }

    pub fn query(&mut self, query: usize) -> io::Result<QueryQueues> {
        let mut open = |tag| -> io::Result<Arc<FileStack>> {
            self.parallelizer
                .create_stack_from_file(
                    tag,
                    get_ref_part_file_name(&self.directory, tag, query, ""),
                )
                .map_err(io::Error::other)?;
            self.parallelizer.get_stack(tag).map_err(io::Error::other)
        };
        Ok(QueryQueues {
            query,
            todo: open(STACK_ALIGN_TODO)?,
            wip: open(STACK_ALIGN_WIP)?,
            done: open(STACK_ALIGN_DONE)?,
            join_todo: open(STACK_JOIN_TODO)?,
            join_wip: open(STACK_JOIN_WIP)?,
            join_done: open(STACK_JOIN_DONE)?,
        })
    }

    pub fn claim_reference(&self, queues: &QueryQueues) -> io::Result<Option<(String, Chunk)>> {
        if stop_exists() {
            return Ok(None);
        }
        let Some(record) = queues.todo.pop_string().map_err(io::Error::other)? else {
            return Ok(None);
        };
        queues.wip.push_string(&record).map_err(io::Error::other)?;
        let chunk = to_chunk(&record).map_err(io::Error::other)?;
        self.parallelizer
            .log(&format!("SEARCH BEGIN {} {}", queues.query, chunk.i))
            .map_err(io::Error::other)?;
        Ok(Some((record, chunk)))
    }

    pub fn finish_reference(
        &self,
        queues: &QueryQueues,
        record: &str,
        chunk: Chunk,
        reference_count: usize,
    ) -> io::Result<()> {
        let (_, done_count) = queues
            .done
            .push_string_size(record)
            .map_err(io::Error::other)?;
        if done_count == reference_count {
            queues
                .join_todo
                .push_string("TOKEN")
                .map_err(io::Error::other)?;
        }
        queues.wip.remove(record).map_err(io::Error::other)?;
        self.parallelizer
            .log(&format!("SEARCH END {} {}", queues.query, chunk.i))
            .map_err(io::Error::other)
    }

    pub fn claim_join(&self, queues: &QueryQueues) -> io::Result<Option<String>> {
        if stop_exists() {
            return Ok(None);
        }
        let Some(token) = queues.join_todo.pop_string().map_err(io::Error::other)? else {
            return Ok(None);
        };
        self.parallelizer
            .log(&format!("JOIN BEGIN {}", queues.query))
            .map_err(io::Error::other)?;
        queues.join_wip.clear().map_err(io::Error::other)?;
        queues.join_done.clear().map_err(io::Error::other)?;
        queues
            .join_wip
            .push_string(&token)
            .map_err(io::Error::other)?;
        queues.join_todo.clear().map_err(io::Error::other)?;
        Ok(Some(token))
    }

    pub fn finish_join(&self, queues: &QueryQueues, token: &str) -> io::Result<()> {
        queues
            .join_done
            .push_string(token)
            .map_err(io::Error::other)?;
        queues.join_wip.pop_string().map_err(io::Error::other)?;
        self.parallelizer
            .log(&format!("JOIN END {}", queues.query))
            .map_err(io::Error::other)
    }
}

pub fn stop_exists() -> bool {
    std::fs::File::open("stop").is_ok()
}

/// Restore source-compatible alignment and join work lists without replaying
/// any completed record or deleting its named intermediate/dictionary files.
pub fn recover(directory: &Path) -> io::Result<()> {
    for query in 0..65_536 {
        let todo_path = get_ref_part_file_name(directory, STACK_ALIGN_TODO, query, "");
        if !todo_path.exists() {
            break;
        }
        let todo = FileStack::new(todo_path);
        let wip = FileStack::new(get_ref_part_file_name(
            directory,
            STACK_ALIGN_WIP,
            query,
            "",
        ));
        let mut count = 0;
        while let Some(record) = wip.pop_string().map_err(io::Error::other)? {
            todo.push_string(&record).map_err(io::Error::other)?;
            count += 1;
        }
        if count != 0 {
            crate::util::log_stream::message_stream()
                .lock()
                .unwrap_or_else(|error| error.into_inner())
                .write(format!("Restored {count} align chunks for query {query}\n"))?;
        }
        let wip = FileStack::new(get_ref_part_file_name(directory, STACK_JOIN_WIP, query, ""));
        if wip.size().map_err(io::Error::other)? != 0 {
            let todo = FileStack::new(get_ref_part_file_name(
                directory,
                STACK_JOIN_TODO,
                query,
                "",
            ));
            let mut count = 0;
            while let Some(record) = wip.pop_string().map_err(io::Error::other)? {
                todo.push_string(&record).map_err(io::Error::other)?;
                count += 1;
            }
            if count != 0 {
                crate::util::log_stream::message_stream()
                    .lock()
                    .unwrap_or_else(|error| error.into_inner())
                    .write(format!("Restored join of query {query}\n"))?;
            }
        }
    }
    if stop_exists() {
        std::fs::remove_file("stop")?;
        crate::util::log_stream::message_stream()
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .write("Removed 'stop' file\n")?;
    }
    Ok(())
}
