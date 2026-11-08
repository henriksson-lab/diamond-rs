//! Forward-only Source TSV input with a bounded producer/consumer queue.
//!
//! Source-profile reader/consumer errors terminate at their actual worker.
//! The reusable library profile returns those errors after joining workers.
use super::read_text_mt::READ_TEXT_MT_SIZE;
use super::{Schema, Table};
use crate::util::io::InputFile;
use std::collections::BTreeMap;
use std::sync::Mutex;

pub struct StreamingFile {
    input: InputFile,
    schema: Schema,
    size: i64,
    source_workers: bool,
}

impl StreamingFile {
    /// Read an existing Source stream (including an unlinked temporary TSV)
    /// without reopening its filename or changing its codec/profile.
    pub fn from_input(schema: Schema, mut input: InputFile) -> Result<Self, String> {
        let size = input.file_size().map_err(|error| error.to_string())?;
        Ok(Self {
            input,
            schema,
            size,
            source_workers: false,
        })
    }

    pub fn open(schema: Schema, path: &str, buffer_size: usize) -> Result<Self, String> {
        let mut input =
            InputFile::with_buffer_size(path, 0, buffer_size).map_err(|error| error.to_string())?;
        let size = input.file_size().map_err(|error| error.to_string())?;
        Ok(Self {
            input,
            schema,
            size,
            source_workers: false,
        })
    }

    pub fn open_source(schema: Schema, path: &str, buffer_size: usize) -> Result<Self, String> {
        let mut input = InputFile::with_source_profile(path, 0, buffer_size)
            .map_err(|error| error.to_string())?;
        let size = input.file_size().map_err(|error| error.to_string())?;
        Ok(Self {
            input,
            schema,
            size,
            source_workers: true,
        })
    }

    /// Concrete command opener retaining FileOpenException metadata across
    /// the legacy string-returning table API boundary.
    pub fn open_source_io(schema: Schema, path: &str, buffer_size: usize) -> std::io::Result<Self> {
        Self::open_source_path_io(schema, std::path::Path::new(path), buffer_size)
    }
    pub fn open_source_path_io(
        schema: Schema,
        path: &std::path::Path,
        buffer_size: usize,
    ) -> std::io::Result<Self> {
        let mut input = InputFile::with_source_profile_path(path, 0, buffer_size)
            .map_err(std::io::Error::other)?;
        let size = input.file_size().map_err(std::io::Error::other)?;
        Ok(Self {
            input,
            schema,
            size,
            source_workers: true,
        })
    }

    pub fn size(&self) -> i64 {
        self.size
    }

    /// Source read_text_mt performs at least one 1 MiB read, extends full
    /// buffers to the next LF, and tests max_size between complete buffers.
    /// Parsed chunks are reordered before building the returned table.
    pub fn read(&mut self, max_size: i64, threads: usize) -> Result<Table, String> {
        Self::read_input(
            &mut self.input,
            self.schema.clone(),
            max_size,
            threads,
            &crate::util::io::IoErrorOwner::default(),
            self.source_workers,
        )
    }

    pub fn read_io(&mut self, max_size: i64, threads: usize) -> std::io::Result<Table> {
        let errors = crate::util::io::IoErrorOwner::default();
        let result = Self::read_input(
            &mut self.input,
            self.schema.clone(),
            max_size,
            threads,
            &errors,
            self.source_workers,
        );
        errors.finish(result)
    }
    /// Borrow the retained Source TSV input owner. No second constructor,
    /// codec detection, descriptor or closing boundary is introduced.
    pub(crate) fn read_input_io(
        input: &mut InputFile,
        schema: Schema,
        max_size: i64,
        threads: usize,
    ) -> std::io::Result<Table> {
        let errors = crate::util::io::IoErrorOwner::default();
        let result = Self::read_input(input, schema, max_size, threads, &errors, true);
        errors.finish(result)
    }
    fn read_input(
        input: &mut InputFile,
        schema: Schema,
        max_size: i64,
        threads: usize,
        errors: &crate::util::io::IoErrorOwner,
        source_workers: bool,
    ) -> Result<Table, String> {
        // Source BuildHelper(file_size/min(max_size)) is constructed on the
        // caller before any reader/consumer thread or physical data read.
        let file_size = errors.capture(input.file_size().map_err(std::io::Error::other))?;
        if max_size < 0 {
            // Source BuildHelper::reserve converts its negative size to an
            // impossible size_t before the producer is started.
            return Err("vector::reserve".into());
        }
        struct OrderedTables {
            next: i64,
            pending: BTreeMap<i64, Table>,
            table: Table,
        }
        let capacity = (max_size.min(file_size) as f64 * 1.1) as usize;
        if capacity > isize::MAX as usize {
            return Err("vector::reserve".into());
        }
        let mut table = Table::new(schema.clone());
        errors.capture(
            table
                .data
                .try_reserve_exact(capacity)
                .map_err(|_| crate::util::io::BadAllocError::default().into_io()),
        )?;
        let ordered = Mutex::new(OrderedTables {
            next: 0,
            pending: BTreeMap::new(),
            table,
        });
        let result =
            read_input_parallel_io(input, max_size, threads, source_workers, |index, bytes| {
                let chunk = parse_chunk_io(schema.clone(), bytes)?;
                let mut ordered = ordered
                    .lock()
                    .map_err(|_| std::io::Error::other("TSV reorder queue poisoned"))?;
                ordered.pending.insert(index, chunk);
                loop {
                    let next = ordered.next;
                    let Some(chunk) = ordered.pending.remove(&next) else {
                        break;
                    };
                    ordered
                        .table
                        .append(&chunk)
                        .map_err(std::io::Error::other)?;
                    ordered.next += 1;
                }
                Ok(())
            });
        errors.capture(result)?;
        let table = ordered
            .into_inner()
            .map_err(|_| "TSV reorder queue poisoned".to_owned())?
            .table;
        Ok(table)
    }
}

/// Shared actual read_text_mt engine: one reader and N consumers, with at
/// most N complete buffers waiting. Only the reader borrows the input owner.
struct ChunkQueue {
    state: Mutex<ChunkState>,
    ready: std::sync::Condvar,
    space: std::sync::Condvar,
    capacity: usize,
}
#[derive(Default)]
struct ChunkState {
    buffers: std::collections::VecDeque<Vec<u8>>,
    next_chunk: i64,
    done: bool,
    cancelled: bool,
    error: Option<std::io::Error>,
}
impl ChunkQueue {
    fn fail(&self, error: std::io::Error, source_workers: bool) {
        if source_workers {
            source_worker_failure(&error);
        }
        let mut state = self.state.lock().unwrap_or_else(|error| error.into_inner());
        if state.error.is_none() {
            state.error = Some(error);
        }
        state.cancelled = true;
        self.ready.notify_all();
        self.space.notify_all();
    }
    fn finish(&self) {
        self.state
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .done = true;
        self.ready.notify_all();
    }
}
struct WorkerLifetime<'a> {
    queue: &'a ChunkQueue,
    producer: bool,
    source_workers: bool,
}
impl Drop for WorkerLifetime<'_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            if self.source_workers {
                crate::util::io::abort_source_process();
            }
            self.queue
                .fail(std::io::Error::other("TSV worker panicked"), false);
        }
        if self.producer {
            self.queue.finish();
        }
    }
}
fn source_worker_failure(error: &std::io::Error) -> ! {
    // Source fread/read_to emits perror at the failed physical operation,
    // then an exception escaping std::thread reaches main's abort handler.
    if let Err(error) = crate::util::io::emit_source_perror(error) {
        crate::util::io::abort_source_diagnostic_error(&error);
    }
    crate::util::io::abort_source_process()
}
fn input_chunks(
    input: &mut InputFile,
    max_size: i64,
) -> impl FnMut() -> std::io::Result<Option<Vec<u8>>> + Send + '_ {
    let mut finished = false;
    let mut total = 0i64;
    move || {
        if finished {
            return Ok(None);
        }
        let mut bytes = vec![0; READ_TEXT_MT_SIZE];
        let count = input.read_raw(&mut bytes).map_err(std::io::Error::other)?;
        bytes.truncate(count);
        if count == READ_TEXT_MT_SIZE {
            input
                .read_to(&mut bytes, b'\n')
                .map_err(std::io::Error::other)?;
        }
        total = total.saturating_add(bytes.len() as i64);
        finished = bytes.len() < READ_TEXT_MT_SIZE
            || total.saturating_add(READ_TEXT_MT_SIZE as i64) > max_size;
        Ok((!bytes.is_empty()).then_some(bytes))
    }
}
fn parallel_chunks_io<R, F>(
    mut reader: R,
    threads: usize,
    source_workers: bool,
    callback: F,
) -> std::io::Result<()>
where
    R: FnMut() -> std::io::Result<Option<Vec<u8>>> + Send,
    F: Fn(i64, &[u8]) -> std::io::Result<()> + Sync,
{
    let consumers = threads.saturating_sub(1).max(1);
    let queue = ChunkQueue {
        state: Mutex::new(ChunkState::default()),
        ready: Default::default(),
        space: Default::default(),
        capacity: consumers,
    };
    std::thread::scope(|scope| -> std::io::Result<()> {
        let q = &queue;
        let producer = std::thread::Builder::new().spawn_scoped(scope, move || {
            let _lifetime = WorkerLifetime {
                queue: q,
                producer: true,
                source_workers,
            };
            loop {
                if q.state
                    .lock()
                    .unwrap_or_else(|error| error.into_inner())
                    .cancelled
                {
                    break;
                }
                let bytes = match reader() {
                    Ok(Some(bytes)) => bytes,
                    Ok(None) => break,
                    Err(error) => {
                        q.fail(error, source_workers);
                        break;
                    }
                };
                let mut state = q.state.lock().unwrap_or_else(|error| error.into_inner());
                while state.buffers.len() >= q.capacity && !state.cancelled {
                    state = q
                        .space
                        .wait(state)
                        .unwrap_or_else(|error| error.into_inner());
                }
                if state.cancelled {
                    break;
                }
                state.buffers.push_back(bytes);
                q.ready.notify_one();
            }
        })?;
        let mut workers = Vec::with_capacity(consumers);
        for _ in 0..consumers {
            let q = &queue;
            let callback = &callback;
            match std::thread::Builder::new().spawn_scoped(scope, move || {
                let _lifetime = WorkerLifetime {
                    queue: q,
                    producer: false,
                    source_workers,
                };
                loop {
                    let (index, bytes) = {
                        let mut state = q.state.lock().unwrap_or_else(|error| error.into_inner());
                        while state.buffers.is_empty() && !state.done && !state.cancelled {
                            state = q
                                .ready
                                .wait(state)
                                .unwrap_or_else(|error| error.into_inner());
                        }
                        if state.cancelled {
                            return;
                        }
                        let Some(bytes) = state.buffers.pop_front() else {
                            return;
                        };
                        let index = state.next_chunk;
                        state.next_chunk += 1;
                        q.space.notify_one();
                        (index, bytes)
                    };
                    if let Err(error) = callback(index, &bytes) {
                        q.fail(error, source_workers);
                        return;
                    }
                }
            }) {
                Ok(worker) => workers.push(worker),
                Err(error) => {
                    queue.fail(error, false);
                    break;
                }
            }
        }
        if producer.join().is_err() {
            queue.fail(std::io::Error::other("TSV reader worker panicked"), false);
        }
        for worker in workers {
            if worker.join().is_err() {
                queue.fail(std::io::Error::other("TSV callback worker panicked"), false);
            }
        }
        let error = queue
            .state
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .error
            .take();
        error.map_or(Ok(()), Err)
    })
}

/// Canonical Source raw callback service. Callback invocations can overlap;
/// the chunk index is assigned under the FIFO queue lock when dequeued.
pub(crate) fn read_input_parallel_io<F>(
    input: &mut InputFile,
    max_size: i64,
    threads: usize,
    source_workers: bool,
    callback: F,
) -> std::io::Result<()>
where
    F: Fn(i64, &[u8]) -> std::io::Result<()> + Sync,
{
    parallel_chunks_io(
        input_chunks(input, max_size),
        threads,
        source_workers,
        callback,
    )
}
pub(crate) fn read_bytes_parallel_io<F>(
    data: &[u8],
    max_size: i64,
    threads: usize,
    callback: F,
) -> std::io::Result<()>
where
    F: Fn(i64, &[u8]) -> std::io::Result<()> + Sync,
{
    let mut position = 0;
    let mut total = 0i64;
    let mut finished = false;
    parallel_chunks_io(
        move || {
            if finished || position >= data.len() {
                return Ok(None);
            }
            let raw_end = position.saturating_add(READ_TEXT_MT_SIZE).min(data.len());
            let mut end = raw_end;
            if raw_end - position == READ_TEXT_MT_SIZE {
                while end < data.len() && data[end] != b'\n' {
                    end += 1;
                }
            }
            let bytes = data[position..end].to_vec();
            total = total.saturating_add(bytes.len() as i64);
            finished = bytes.len() < READ_TEXT_MT_SIZE
                || total.saturating_add(READ_TEXT_MT_SIZE as i64) > max_size;
            position = if end < data.len() && data[end] == b'\n' {
                end + 1
            } else {
                end
            };
            Ok(Some(bytes))
        },
        threads,
        false,
        callback,
    )
}
/// Compatibility FnMut service: the callback stays on its caller's thread,
/// but physical Source reads still belong to a spawned, bounded producer.
pub(crate) fn read_input_serial_io<F>(
    input: &mut InputFile,
    max_size: i64,
    threads: usize,
    source_workers: bool,
    mut callback: F,
) -> std::io::Result<()>
where
    F: FnMut(i64, &[u8]),
{
    let (send, receive) = std::sync::mpsc::sync_channel(threads.saturating_sub(1).max(1));
    let mut reader = input_chunks(input, max_size);
    std::thread::scope(|scope| -> std::io::Result<()> {
        let producer =
            std::thread::Builder::new().spawn_scoped(scope, move || -> std::io::Result<()> {
                loop {
                    let bytes = match reader() {
                        Ok(Some(bytes)) => bytes,
                        Ok(None) => return Ok(()),
                        Err(error) => {
                            if source_workers {
                                source_worker_failure(&error);
                            }
                            return Err(error);
                        }
                    };
                    if send.send(bytes).is_err() {
                        return Ok(());
                    }
                }
            })?;
        for (index, bytes) in receive.into_iter().enumerate() {
            callback(index as i64, &bytes);
        }
        producer
            .join()
            .map_err(|_| std::io::Error::other("TSV reader worker panicked"))?
    })
}

#[cfg(test)]
fn parse_chunk(schema: Schema, bytes: &[u8]) -> Result<Table, String> {
    parse_chunk_io(schema, bytes).map_err(|error| error.to_string())
}
pub(crate) fn parse_chunk_io(schema: Schema, bytes: &[u8]) -> std::io::Result<Table> {
    let mut table = Table::new(schema);
    let mut begin = 0usize;
    while begin < bytes.len() {
        let end = bytes[begin..]
            .iter()
            .position(|&byte| byte == b'\n')
            .map_or(bytes.len(), |length| begin + length);
        let mut line = &bytes[begin..end];
        // Source CharTokenizer removes CR before LF, and also in the final
        // unterminated record through its end-pointer branch.
        if line.last() == Some(&b'\r') {
            line = &line[..line.len() - 1];
        }
        table.push_back_raw_line_io(line, None)?;
        begin = end.saturating_add(1);
    }
    Ok(table)
}

#[cfg(test)]
mod profile_tests {
    use super::*;
    #[test]
    fn reusable_library_parser_preserves_raw_error_without_process_abort() {
        let error = parse_chunk_io(vec![super::super::Type::Int64], b"7\xff\n").unwrap_err();
        assert_eq!(
            crate::util::io::RawDiagnosticError::find(&error).unwrap(),
            b"Error converting integer value: 7\xff"
        );
    }
}
