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
    errors: crate::util::io::IoErrorOwner,
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
        Self::new_owned(directory, Default::default())
    }

    pub fn new_owned(directory: &Path, errors: crate::util::io::IoErrorOwner) -> io::Result<Self> {
        let mut parallelizer = Parallelizer::new_owned(errors.clone());
        errors.finish(parallelizer.init(directory))?;
        Ok(Self {
            parallelizer,
            directory: directory.into(),
            errors,
        })
    }

    pub fn query(&mut self, query: usize) -> io::Result<QueryQueues> {
        let mut open = |tag| -> io::Result<Arc<FileStack>> {
            self.parallelizer
                .create_stack_from_file(
                    tag,
                    get_ref_part_file_name(&self.directory, tag, query, ""),
                )
                .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
            self.parallelizer
                .get_stack(tag)
                .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())
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
        let Some(record) = queues
            .todo
            .pop_string()
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?
        else {
            return Ok(None);
        };
        queues
            .wip
            .push_string(&record)
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        let chunk =
            to_chunk(&record).map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        self.parallelizer
            .log(&format!("SEARCH BEGIN {} {}", queues.query, chunk.i))
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
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
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        if done_count == reference_count {
            queues
                .join_todo
                .push_string("TOKEN")
                .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        }
        queues
            .wip
            .remove(record)
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        self.parallelizer
            .log(&format!("SEARCH END {} {}", queues.query, chunk.i))
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())
    }

    pub fn claim_join(&self, queues: &QueryQueues) -> io::Result<Option<String>> {
        if stop_exists() {
            return Ok(None);
        }
        let Some(token) = queues
            .join_todo
            .pop_string()
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?
        else {
            return Ok(None);
        };
        self.parallelizer
            .log(&format!("JOIN BEGIN {}", queues.query))
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        queues
            .join_wip
            .clear()
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        queues
            .join_done
            .clear()
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        queues
            .join_wip
            .push_string(&token)
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        queues
            .join_todo
            .clear()
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        Ok(Some(token))
    }

    pub fn finish_join(&self, queues: &QueryQueues, token: &str) -> io::Result<()> {
        queues
            .join_done
            .push_string(token)
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        queues
            .join_wip
            .pop_string()
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())?;
        self.parallelizer
            .log(&format!("JOIN END {}", queues.query))
            .map_err(|error| self.errors.finish::<()>(Err(error)).unwrap_err())
    }
}

fn source_file_exists(path: &Path) -> bool {
    // Source file_exists constructs ifstream and tests good(), not stat.
    crate::util::io::source_text::SourceTextInput::open_unchecked(path).good()
}

pub fn stop_exists() -> bool {
    source_file_exists(Path::new("stop"))
}

/// Restore source-compatible alignment and join work lists without replaying
/// any completed record or deleting its named intermediate/dictionary files.
pub fn recover(directory: &Path) -> io::Result<()> {
    recover_owned(directory, &Default::default())
}

pub fn recover_owned(directory: &Path, errors: &crate::util::io::IoErrorOwner) -> io::Result<()> {
    recover_with_stop(directory, errors, Path::new("stop"))
}

fn recover_with_stop(
    directory: &Path,
    errors: &crate::util::io::IoErrorOwner,
    stop: &Path,
) -> io::Result<()> {
    for query in 0..65_536 {
        let todo_path = get_ref_part_file_name(directory, STACK_ALIGN_TODO, query, "");
        if !source_file_exists(&todo_path) {
            break;
        }
        let todo = errors.finish(FileStack::try_new_owned(todo_path, errors.clone()))?;
        let wip = errors.finish(FileStack::try_new_owned(
            get_ref_part_file_name(directory, STACK_ALIGN_WIP, query, ""),
            errors.clone(),
        ))?;
        let mut count = 0;
        while let Some(record) = wip
            .pop_string()
            .map_err(|error| errors.finish::<()>(Err(error)).unwrap_err())?
        {
            todo.push_string(&record)
                .map_err(|error| errors.finish::<()>(Err(error)).unwrap_err())?;
            count += 1;
        }
        if count != 0 {
            crate::util::log_stream::message_stream()
                .lock()
                .unwrap_or_else(|error| error.into_inner())
                .write(format!("Restored {count} align chunks for query {query}\n"))?;
        }
        let wip = errors.finish(FileStack::try_new_owned(
            get_ref_part_file_name(directory, STACK_JOIN_WIP, query, ""),
            errors.clone(),
        ))?;
        if wip
            .size()
            .map_err(|error| errors.finish::<()>(Err(error)).unwrap_err())?
            != 0
        {
            let todo = errors.finish(FileStack::try_new_owned(
                get_ref_part_file_name(directory, STACK_JOIN_TODO, query, ""),
                errors.clone(),
            ))?;
            let mut count = 0;
            while let Some(record) = wip
                .pop_string()
                .map_err(|error| errors.finish::<()>(Err(error)).unwrap_err())?
            {
                todo.push_string(&record)
                    .map_err(|error| errors.finish::<()>(Err(error)).unwrap_err())?;
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
    if source_file_exists(stop) {
        let _ = std::fs::remove_file(crate::util::io::source_c_path(stop).as_ref());
        crate::util::log_stream::message_stream()
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .write("Removed 'stop' file\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn root(tag: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-mpp-{tag}-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ))
    }

    #[test]
    fn recovery_of_absent_work_directory_creates_nothing() {
        let root = root("absent");
        recover_with_stop(
            &root.join("nested"),
            &Default::default(),
            &root.join("stop"),
        )
        .unwrap();
        assert!(!root.exists());
    }

    #[test]
    fn recovery_restores_only_inflight_reference_and_join_records() {
        let root = root("restore");
        std::fs::create_dir(&root).unwrap();
        let path = |tag| get_ref_part_file_name(&root, tag, 0, "");
        let todo = FileStack::try_new(path(STACK_ALIGN_TODO)).unwrap();
        todo.push_string("0 0 1").unwrap();
        let wip = FileStack::try_new(path(STACK_ALIGN_WIP)).unwrap();
        wip.push_string("1 1 1").unwrap();
        wip.push_string("2 2 1").unwrap();
        let join = FileStack::try_new(path(STACK_JOIN_WIP)).unwrap();
        join.push_string("TOKEN").unwrap();
        let done = path(STACK_ALIGN_DONE);
        std::fs::write(&done, b"already completed\n").unwrap();
        let stop = root.join("stop");
        std::fs::write(&stop, []).unwrap();
        recover_with_stop(&root, &Default::default(), &stop).unwrap();
        assert_eq!(todo.pop_string().unwrap().as_deref(), Some("1 1 1"));
        assert_eq!(todo.pop_string().unwrap().as_deref(), Some("2 2 1"));
        assert_eq!(todo.pop_string().unwrap().as_deref(), Some("0 0 1"));
        assert_eq!(todo.pop_string().unwrap(), None);
        assert_eq!(wip.size().unwrap(), 0);
        assert_eq!(join.size().unwrap(), 0);
        let restored_join = FileStack::try_new(path(STACK_JOIN_TODO)).unwrap();
        assert_eq!(
            restored_join.pop_string().unwrap().as_deref(),
            Some("TOKEN")
        );
        assert_eq!(std::fs::read(done).unwrap(), b"already completed\n");
        assert!(!stop.exists());
        assert!(!get_ref_part_file_name(&root, STACK_ALIGN_TODO, 1, "").exists());
        drop((todo, wip, join, restored_join));
        std::fs::remove_dir_all(root).unwrap();
    }
}
