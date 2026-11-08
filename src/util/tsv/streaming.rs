//! Forward-only Source TSV input with a bounded producer/consumer queue.
//!
//! Source lets a malformed-row exception escape its consumer std::thread,
//! terminating the process. Native consumers retain the schema diagnostic and
//! return it after joining workers and releasing their input buffers instead.
use super::read_text_mt::READ_TEXT_MT_SIZE;
use super::{Schema, Table};
use crate::util::io::InputFile;
use std::collections::BTreeMap;
use std::sync::{mpsc, Mutex};

pub struct StreamingFile {
    input: InputFile,
    schema: Schema,
    size: i64,
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
        input.file_size().map_err(std::io::Error::other)?;
        let errors = crate::util::io::IoErrorOwner::default();
        let result = Self::read_input(input, schema, max_size, threads, &errors);
        errors.finish(result)
    }
    fn read_input(
        input: &mut InputFile,
        schema: Schema,
        max_size: i64,
        threads: usize,
        errors: &crate::util::io::IoErrorOwner,
    ) -> Result<Table, String> {
        if max_size < 0 {
            // Source BuildHelper::reserve converts its negative size to an
            // impossible size_t before the producer is started.
            return Err("vector::reserve".into());
        }
        let consumers = threads.saturating_sub(1).max(1);
        let (jobs, receiver) = mpsc::sync_channel::<(usize, Vec<u8>)>(consumers);
        let receiver = Mutex::new(receiver);
        let (results, completed) = mpsc::channel();
        let schema = &schema;
        let chunks = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(consumers);
            for _ in 0..consumers {
                let receiver = &receiver;
                let results = results.clone();
                handles.push(scope.spawn(move || loop {
                    let job = receiver
                        .lock()
                        .map_err(|_| "TSV input queue lock poisoned".to_owned())?
                        .recv();
                    let Ok((index, bytes)) = job else {
                        return Ok::<_, String>(());
                    };
                    let table = parse_chunk(schema.clone(), &bytes);
                    if results.send((index, table)).is_err() {
                        return Ok(());
                    }
                }));
            }
            drop(results);
            let loaded = (|| {
                let mut index = 0usize;
                let mut total = 0i64;
                loop {
                    let mut bytes = vec![0; READ_TEXT_MT_SIZE];
                    let n = errors
                        .capture(input.read_raw(&mut bytes).map_err(std::io::Error::other))?;
                    bytes.truncate(n);
                    if n == READ_TEXT_MT_SIZE {
                        errors.capture(
                            input
                                .read_to(&mut bytes, b'\n')
                                .map_err(std::io::Error::other),
                        )?;
                    }
                    let n = bytes.len();
                    total = total.saturating_add(n as i64);
                    if n > 0 {
                        jobs.send((index, bytes))
                            .map_err(|_| "TSV input consumers stopped".to_owned())?;
                        index += 1;
                    }
                    if n < READ_TEXT_MT_SIZE
                        || total.saturating_add(READ_TEXT_MT_SIZE as i64) > max_size
                    {
                        break;
                    }
                }
                Ok::<_, String>(())
            })();
            drop(jobs);
            let chunks = completed.into_iter().collect::<BTreeMap<_, _>>();
            for handle in handles {
                handle
                    .join()
                    .map_err(|_| "TSV parser worker panicked".to_owned())??;
            }
            loaded?;
            Ok::<_, String>(chunks)
        })?;
        let mut table = Table::new(schema.clone());
        for (_, chunk) in chunks {
            table.append(&chunk?).map_err(|error| error.to_string())?;
        }
        Ok(table)
    }
}

fn parse_chunk(schema: Schema, bytes: &[u8]) -> Result<Table, String> {
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
        table.push_back_raw_line(line, None)?;
        begin = end.saturating_add(1);
    }
    Ok(table)
}
