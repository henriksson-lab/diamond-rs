use std::collections::HashMap;
use std::collections::VecDeque;
use std::fs::File as StdFile;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicI64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex as StdMutex, OnceLock};
use std::thread::{JoinHandle, ThreadId};
use std::time::Duration;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Sync;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Async;

pub trait MutexTag {
    type State: Default;
    fn lock(state: &Self::State);
    fn unlock(state: &Self::State);
}

impl MutexTag for Sync {
    type State = ();

    fn lock(_: &Self::State) {}
    fn unlock(_: &Self::State) {}
}

impl MutexTag for Async {
    type State = AtomicBool;

    fn lock(state: &Self::State) {
        while state
            .compare_exchange(false, true, Ordering::Acquire, Ordering::Relaxed)
            .is_err()
        {
            std::hint::spin_loop();
        }
    }

    fn unlock(state: &Self::State) {
        state.store(false, Ordering::Release);
    }
}

#[derive(Debug, Default)]
pub struct Mutex<Tag: MutexTag> {
    state: Tag::State,
}

impl<Tag: MutexTag> Mutex<Tag> {
    pub fn new() -> Self {
        Self {
            state: Tag::State::default(),
        }
    }

    pub fn lock(&self) {
        Tag::lock(&self.state);
    }

    pub fn unlock(&self) {
        Tag::unlock(&self.state);
    }
}

#[derive(Debug)]
pub struct CountingSemaphore<const LEAST_MAX_VALUE: isize = 2_147_483_647> {
    counter: StdMutex<isize>,
    cv: Condvar,
}

impl<const LEAST_MAX_VALUE: isize> CountingSemaphore<LEAST_MAX_VALUE> {
    pub const fn max() -> isize {
        LEAST_MAX_VALUE
    }

    pub fn new(desired: isize) -> Self {
        Self {
            counter: StdMutex::new(desired),
            cv: Condvar::new(),
        }
    }

    pub fn release(&self, update: isize) {
        let mut counter = self.counter.lock().unwrap();
        *counter += update;
        if update > 1 {
            self.cv.notify_all();
        } else {
            self.cv.notify_one();
        }
    }

    pub fn acquire(&self) {
        let mut counter = self.counter.lock().unwrap();
        while *counter <= 0 {
            counter = self.cv.wait(counter).unwrap();
        }
        *counter -= 1;
    }
}

pub fn pool_worker<F>(partition: &AtomicUsize, thread_id: usize, partition_count: usize, f: &F)
where
    F: Fn(usize, usize) + std::marker::Sync,
{
    loop {
        let p = partition.fetch_add(1, Ordering::Relaxed);
        if p >= partition_count {
            return;
        }
        f(p, thread_id);
    }
}

pub fn scheduled_thread_pool<F>(thread_count: usize, f: &F)
where
    F: Fn(&AtomicUsize, usize) + std::marker::Sync,
{
    let partition = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for i in 0..thread_count {
            let partition = &partition;
            scope.spawn(move || f(partition, i));
        }
    });
}

pub fn scheduled_thread_pool_auto<F>(thread_count: usize, partition_count: usize, f: &F)
where
    F: Fn(usize, usize) + std::marker::Sync,
{
    let partition = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for i in 0..thread_count {
            let partition = &partition;
            scope.spawn(move || pool_worker(partition, i, partition_count, f));
        }
    });
}

pub fn launch_threads<F>(thread_count: i32, f: &F)
where
    F: Fn() + std::marker::Sync,
{
    std::thread::scope(|scope| {
        for _ in 0..thread_count {
            scope.spawn(f);
        }
    });
}

pub struct SimpleThreadPool {
    threads: HashMap<ThreadId, JoinHandle<()>>,
    stop_flag: Arc<AtomicBool>,
    first_exception: Arc<StdMutex<Option<String>>>,
}

impl SimpleThreadPool {
    pub fn new() -> Self {
        Self {
            threads: HashMap::new(),
            stop_flag: Arc::new(AtomicBool::new(false)),
            first_exception: Arc::new(StdMutex::new(None)),
        }
    }

    pub fn stop(&self) -> &AtomicBool {
        &self.stop_flag
    }

    pub fn request_stop(&self) {
        self.stop_flag.store(true, Ordering::Relaxed);
    }

    pub fn spawn<F>(&mut self, func: F) -> ThreadId
    where
        F: FnOnce(&AtomicBool) -> Result<(), String> + Send + 'static,
    {
        let stop = Arc::clone(&self.stop_flag);
        let first_exception = Arc::clone(&self.first_exception);
        let handle = std::thread::spawn(move || {
            let result =
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| func(stop.as_ref())));
            let message = match result {
                Ok(Ok(())) => return,
                Ok(Err(e)) => e,
                Err(payload) => {
                    if let Some(s) = payload.downcast_ref::<&str>() {
                        (*s).to_string()
                    } else if let Some(s) = payload.downcast_ref::<String>() {
                        s.clone()
                    } else {
                        "Thread panicked.".to_string()
                    }
                }
            };
            let mut exception = first_exception.lock().unwrap();
            if exception.is_none() {
                *exception = Some(message);
                stop.store(true, Ordering::Relaxed);
            }
        });
        let id = handle.thread().id();
        self.threads.insert(id, handle);
        id
    }

    pub fn spawn_method<T, F>(&mut self, obj: Arc<StdMutex<T>>, method: F) -> ThreadId
    where
        T: Send + 'static,
        F: FnOnce(&mut T, &AtomicBool) -> Result<(), String> + Send + 'static,
    {
        self.spawn(move |stop| {
            let mut guard = obj.lock().unwrap();
            method(&mut *guard, stop)
        })
    }

    pub fn join_all(&mut self) -> Result<(), String> {
        for (_, handle) in self.threads.drain() {
            if let Err(payload) = handle.join() {
                let message = if let Some(s) = payload.downcast_ref::<&str>() {
                    (*s).to_string()
                } else if let Some(s) = payload.downcast_ref::<String>() {
                    s.clone()
                } else {
                    "Thread panicked.".to_string()
                };
                let mut exception = self.first_exception.lock().unwrap();
                if exception.is_none() {
                    *exception = Some(message);
                    self.stop_flag.store(true, Ordering::Relaxed);
                }
            }
        }
        if let Some(message) = self.first_exception.lock().unwrap().take() {
            Err(message)
        } else {
            Ok(())
        }
    }

    pub fn join_ids<I>(&mut self, ids: I) -> Result<(), String>
    where
        I: IntoIterator<Item = ThreadId>,
    {
        for id in ids {
            self.join(id)?;
        }
        if let Some(message) = self.first_exception.lock().unwrap().take() {
            Err(message)
        } else {
            Ok(())
        }
    }

    pub fn join(&mut self, thread_id: ThreadId) -> Result<(), String> {
        let handle = self
            .threads
            .remove(&thread_id)
            .ok_or_else(|| "Thread ID not found in thread pool.".to_string())?;
        if let Err(payload) = handle.join() {
            let message = if let Some(s) = payload.downcast_ref::<&str>() {
                (*s).to_string()
            } else if let Some(s) = payload.downcast_ref::<String>() {
                s.clone()
            } else {
                "Thread panicked.".to_string()
            };
            let mut exception = self.first_exception.lock().unwrap();
            if exception.is_none() {
                *exception = Some(message);
                self.stop_flag.store(true, Ordering::Relaxed);
            }
        }
        Ok(())
    }
}

impl Default for SimpleThreadPool {
    fn default() -> Self {
        Self::new()
    }
}

impl Drop for SimpleThreadPool {
    fn drop(&mut self) {
        self.stop_flag.store(true, Ordering::Relaxed);
        for (_, handle) in self.threads.drain() {
            let _ = handle.join();
        }
    }
}

#[derive(Debug)]
pub struct FileStack {
    file_name: PathBuf,
    file: StdMutex<StdFile>,
    max_line_length: usize,
    mtx: StdMutex<()>,
}

impl Default for FileStack {
    fn default() -> Self {
        Self::new("default_stack.idx")
    }
}

impl FileStack {
    pub const DEFAULT_MAX_LINE_LENGTH: usize = 4096;

    pub fn new(file_name: impl AsRef<Path>) -> Self {
        Self::with_max_line_length(file_name, Self::DEFAULT_MAX_LINE_LENGTH)
    }

    pub fn with_max_line_length(file_name: impl AsRef<Path>, maximum_line_length: usize) -> Self {
        let file_name = file_name.as_ref().to_path_buf();
        let file = std::fs::OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .open(&file_name)
            .unwrap_or_else(|_| panic!("could not open file {}", file_name.display()));
        let mut out = Self {
            file_name,
            file: StdMutex::new(file),
            max_line_length: Self::DEFAULT_MAX_LINE_LENGTH,
            mtx: StdMutex::new(()),
        };
        out.set_max_line_length(maximum_line_length);
        out
    }

    pub fn size(&self) -> Result<usize, String> {
        let _lock = self.mtx.lock().unwrap();
        let data = std::fs::read(&self.file_name).map_err(|e| e.to_string())?;
        Ok(data.iter().filter(|&&c| c == b'\n').count())
    }

    pub fn pop_i64(&self) -> Result<i64, String> {
        let mut i = 0i64;
        self.pop_i64_into(&mut i)?;
        Ok(i)
    }

    pub fn pop_i64_into(&self, i: &mut i64) -> Result<i64, String> {
        let _lock = self.mtx.lock().unwrap();
        match self.pop_non_locked_string(false)? {
            Some(buf) => {
                *i = buf.parse::<i64>().map_err(|e| e.to_string())?;
                Ok(*i)
            }
            None => {
                *i = -1;
                Ok(-1)
            }
        }
    }

    pub fn pop_string(&self) -> Result<Option<String>, String> {
        let _lock = self.mtx.lock().unwrap();
        self.pop_non_locked_string(false)
    }

    pub fn pop_string_size(&self) -> Result<(Option<String>, usize), String> {
        let _lock = self.mtx.lock().unwrap();
        let out = self.pop_non_locked_string(false)?;
        let size_after_pop = std::fs::read(&self.file_name)
            .map_err(|e| e.to_string())?
            .iter()
            .filter(|&&c| c == b'\n')
            .count();
        Ok((out, size_after_pop))
    }

    pub fn top_i64(&self) -> Result<i64, String> {
        let _lock = self.mtx.lock().unwrap();
        match self.pop_non_locked_string(true)? {
            Some(buf) => buf.parse::<i64>().map_err(|e| e.to_string()),
            None => Ok(-1),
        }
    }

    pub fn top_i64_into(&self, i: &mut i64) -> Result<i64, String> {
        *i = self.top_i64()?;
        Ok(*i)
    }

    pub fn top_string(&self) -> Result<Option<String>, String> {
        let _lock = self.mtx.lock().unwrap();
        self.pop_non_locked_string(true)
    }

    pub fn remove(&self, line: &str) -> Result<(), String> {
        let _lock = self.mtx.lock().unwrap();
        let data = std::fs::read_to_string(&self.file_name).map_err(|e| e.to_string())?;
        let mut out = String::new();
        for token in data.split('\n') {
            if !token.is_empty() && token != line {
                out.push_str(token);
                out.push('\n');
            }
        }
        std::fs::write(&self.file_name, out).map_err(|e| e.to_string())
    }

    pub fn push_i64(&self, i: i64) -> Result<i64, String> {
        let _lock = self.mtx.lock().unwrap();
        self.push_non_locked_string(&i.to_string())
    }

    pub fn push_string(&self, buf: &str) -> Result<i64, String> {
        let _lock = self.mtx.lock().unwrap();
        self.push_non_locked_string(buf)
    }

    pub fn push_string_size(&self, buf: &str) -> Result<(i64, usize), String> {
        let _lock = self.mtx.lock().unwrap();
        let written = self.push_non_locked_string(buf)?;
        let size_after_push = std::fs::read(&self.file_name)
            .map_err(|e| e.to_string())?
            .iter()
            .filter(|&&c| c == b'\n')
            .count();
        Ok((written, size_after_push))
    }

    pub fn fetch_add(&self, n: i64) -> Result<i64, String> {
        let _lock = self.mtx.lock().unwrap();
        let mut i = match self.pop_non_locked_string(false)? {
            Some(buf) => buf.parse::<i64>().map_err(|e| e.to_string())?,
            None => 0,
        };
        let out = i;
        i += n;
        self.push_non_locked_string(&i.to_string())?;
        Ok(out)
    }

    pub fn fetch_add_one(&self) -> Result<i64, String> {
        self.fetch_add(1)
    }

    pub fn get_max_line_length(&self) -> usize {
        self.max_line_length
    }

    pub fn set_max_line_length(&mut self, n: usize) -> usize {
        self.max_line_length = n.max(8);
        self.max_line_length
    }

    pub fn clear(&self) -> Result<(), String> {
        let _lock = self.mtx.lock().unwrap();
        std::fs::write(&self.file_name, []).map_err(|e| e.to_string())
    }

    pub fn seek(&self, offset: i64, mode: SeekFrom) -> Result<i64, String> {
        let target = match mode {
            SeekFrom::Start(_) => SeekFrom::Start(offset as u64),
            SeekFrom::End(_) => SeekFrom::End(offset),
            SeekFrom::Current(_) => SeekFrom::Current(offset),
        };
        self.file
            .lock()
            .unwrap()
            .seek(target)
            .map(|p| p as i64)
            .map_err(|e| e.to_string())
    }

    pub fn read(&self, buf: &mut [u8]) -> Result<usize, String> {
        self.file
            .lock()
            .unwrap()
            .read(buf)
            .map_err(|e| format!("Error reading file {}: {e}", self.file_name.display()))
    }

    pub fn write(&self, buf: &[u8]) -> Result<i64, String> {
        self.file
            .lock()
            .unwrap()
            .write(buf)
            .map(|n| n as i64)
            .map_err(|e| format!("Error writing file {}: {e}", self.file_name.display()))
    }

    pub fn truncate(&self, size: usize) -> Result<(), String> {
        self.file
            .lock()
            .unwrap()
            .set_len(size as u64)
            .map_err(|e| e.to_string())
    }

    pub fn poll_query(&self, query: &str, sleep_s: f64, max_iter: usize) -> Result<bool, String> {
        for _ in 0..max_iter {
            if let Some(buf) = self.top_string()? {
                if buf.contains(query) {
                    return Ok(true);
                }
                if buf.contains("STOP") {
                    return Err(format!("STOP on FileStack {}", self.file_name.display()));
                }
            }
            std::thread::sleep(Duration::from_secs_f64(sleep_s));
        }
        Err(format!(
            "Could not discover keyword {} on FileStack {} within {} seconds.",
            query,
            self.file_name.display(),
            max_iter as f64 * sleep_s
        ))
    }

    pub fn poll_size(&self, size: usize, sleep_s: f64, max_iter: usize) -> Result<bool, String> {
        for _ in 0..max_iter {
            if self.size()? == size {
                return Ok(true);
            }
            std::thread::sleep(Duration::from_secs_f64(sleep_s));
        }
        Err(format!(
            "Could not detect size {} of FileStack {} within {} seconds.",
            size,
            self.file_name.display(),
            max_iter as f64 * sleep_s
        ))
    }

    pub fn file_name(&self) -> &Path {
        &self.file_name
    }

    pub fn pop_exclusive(&self) -> Result<Option<String>, String> {
        self.pop_non_locked_string(false)
    }

    pub fn push_exclusive(&self, buf: &str) -> Result<i64, String> {
        self.push_non_locked_string(buf)
    }

    fn pop_non_locked_string(&self, keep_flag: bool) -> Result<Option<String>, String> {
        let mut data = std::fs::read(&self.file_name).map_err(|e| e.to_string())?;
        if data.is_empty() {
            return Ok(None);
        }
        if data.last() == Some(&b'\n') {
            data.pop();
        }
        let begin = data.iter().rposition(|&c| c == b'\n').map_or(0, |i| i + 1);
        let line = String::from_utf8(data[begin..].to_vec()).map_err(|e| e.to_string())?;
        if !keep_flag {
            data.truncate(begin);
            std::fs::write(&self.file_name, data).map_err(|e| e.to_string())?;
        }
        Ok(Some(line))
    }

    fn push_non_locked_string(&self, buf: &str) -> Result<i64, String> {
        use std::io::Write;

        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.file_name)
            .map_err(|e| e.to_string())?;
        file.write_all(buf.as_bytes()).map_err(|e| e.to_string())?;
        let mut n = buf.len() as i64;
        if !buf.ends_with('\n') {
            file.write_all(b"\n").map_err(|e| e.to_string())?;
            n += 1;
        }
        Ok(n)
    }
}

#[derive(Debug)]
pub struct Atomic {
    stack: FileStack,
}

impl Atomic {
    pub fn new(file_name: impl AsRef<Path>) -> Self {
        Self {
            stack: FileStack::new(file_name),
        }
    }

    pub fn get(&self) -> Result<i64, String> {
        let i = self.stack.top_i64()?;
        Ok(if i >= 0 { i } else { 0 })
    }

    pub fn fetch_add(&self, n: i64) -> Result<i64, String> {
        self.stack.fetch_add(n)
    }

    pub fn fetch_add_one(&self) -> Result<i64, String> {
        self.fetch_add(1)
    }

    pub fn await_value(&self, n: i64) -> Result<(), String> {
        loop {
            if self.stack.top_i64()? >= n {
                return Ok(());
            }
            std::thread::sleep(Duration::from_secs(1));
        }
    }
}

#[derive(Debug)]
pub struct Parallelizer {
    work_directory: PathBuf,
    barrier_file: PathBuf,
    rank: i32,
    id: String,
    n_registered: i32,
    master_flag: bool,
    i_barrier: i32,
    initialized: bool,
    continuous_cleanup_list: Vec<PathBuf>,
    final_cleanup_list: Vec<PathBuf>,
    fs_map: HashMap<String, Arc<FileStack>>,
}

static PARALLELIZER_INSTANCE: OnceLock<Arc<StdMutex<Parallelizer>>> = OnceLock::new();

impl Default for Parallelizer {
    fn default() -> Self {
        Self::new()
    }
}

impl Parallelizer {
    pub const LOG: &'static str = "log";
    pub const COMMAND: &'static str = "command";
    pub const WORKERS: &'static str = "workers";
    pub const REGISTER: &'static str = "register";
    const BARRIER: &'static str = "barrier";

    pub fn get() -> Arc<StdMutex<Parallelizer>> {
        Arc::clone(PARALLELIZER_INSTANCE.get_or_init(|| Arc::new(StdMutex::new(Self::new()))))
    }

    pub fn new() -> Self {
        Self {
            work_directory: PathBuf::from("parallelizer"),
            barrier_file: PathBuf::new(),
            rank: 0,
            id: String::new(),
            n_registered: 0,
            master_flag: true,
            i_barrier: 0,
            initialized: false,
            continuous_cleanup_list: Vec::new(),
            final_cleanup_list: Vec::new(),
            fs_map: HashMap::new(),
        }
    }

    pub fn init(&mut self, tempdir: impl AsRef<Path>) -> Result<(), String> {
        let tempdir = tempdir.as_ref();
        if !tempdir.as_os_str().is_empty() {
            self.work_directory = tempdir.join(&self.work_directory);
        }
        std::fs::create_dir_all(&self.work_directory).map_err(|e| e.to_string())?;

        let hostname = std::env::var("HOSTNAME")
            .or_else(|_| std::env::var("COMPUTERNAME"))
            .unwrap_or_else(|_| "localhost".to_string());
        self.id = format!("{}_{}", hostname, std::process::id());

        self.create_stack(Self::LOG, &self.id.clone())?;
        self.create_stack(Self::COMMAND, "")?;
        self.create_stack(Self::WORKERS, "")?;
        self.create_stack(Self::REGISTER, "")?;

        self.barrier_file = self.work_directory.join(Self::BARRIER);
        self.log("PARALLELIZER BEGIN")?;
        self.initialized = true;
        Ok(())
    }

    pub fn clear(&mut self) {}

    pub fn get_rank(&self) -> i32 {
        self.rank
    }

    pub fn get_id(&self) -> &str {
        &self.id
    }

    pub fn get_work_directory(&self) -> &Path {
        &self.work_directory
    }

    pub fn get_n_registered(&self) -> i32 {
        self.n_registered
    }

    pub fn is_master(&self) -> bool {
        self.master_flag
    }

    pub fn register_workers(&mut self, sleep_s: f64) -> Result<bool, String> {
        self.get_stack(Self::REGISTER)?
            .push_string(&self.id)
            .map(|_| ())?;
        Self::sleep(sleep_s);
        if self.is_master() {
            while let Some(line) = self.get_stack(Self::REGISTER)?.pop_string()? {
                self.get_stack(Self::WORKERS)?.push_string(&line)?;
                self.n_registered += 1;
            }
        }
        Ok(true)
    }

    pub fn barrier(&mut self, tag: &str) -> Result<bool, String> {
        if !self.initialized {
            return Ok(false);
        }

        let cmd_file_name = self.get_barrier_file_name("cmd", tag, self.i_barrier);
        let cmd_fs = FileStack::new(&cmd_file_name);
        let ack_file_name = self.get_barrier_file_name("ack", tag, self.i_barrier);
        let ack_fs = FileStack::new(&ack_file_name);

        let msg = "WAIT";
        if self.is_master() {
            ack_fs.clear()?;
            cmd_fs.push_string(msg)?;
        }
        cmd_fs.poll_query(msg, 0.0, 1)?;
        ack_fs.push_string(&self.id)?;

        let msg_ok = "GOON";
        if self.is_master() {
            let n_workers = self.get_stack(Self::WORKERS)?.size()?;
            ack_fs.poll_size(n_workers, 0.0, 1)?;
            cmd_fs.push_string(msg_ok)?;
        }
        cmd_fs.poll_query(msg_ok, 0.0, 1)?;

        if self.is_master() {
            Self::clean(&mut self.continuous_cleanup_list);
            self.continuous_cleanup_list.push(cmd_file_name);
            self.continuous_cleanup_list.push(ack_file_name);
        }

        self.i_barrier += 1;
        Ok(true)
    }

    pub fn create_stack(&mut self, tag: &str, sfx: &str) -> Result<bool, String> {
        if self.fs_map.contains_key(tag) {
            return Ok(false);
        }
        let suffix = if sfx.is_empty() {
            String::new()
        } else {
            format!("_{sfx}")
        };
        let file_name = self.work_directory.join(format!("{tag}{suffix}"));
        self.create_stack_from_file(tag, &file_name)
    }

    pub fn create_stack_from_file(
        &mut self,
        tag: &str,
        file_name: impl AsRef<Path>,
    ) -> Result<bool, String> {
        self.delete_stack(tag);
        let file_stack = Arc::new(FileStack::new(file_name.as_ref()));
        self.fs_map.insert(tag.to_string(), file_stack);
        Ok(true)
    }

    pub fn delete_stack(&mut self, tag: &str) -> bool {
        self.fs_map.remove(tag).is_some()
    }

    pub fn get_stack(&self, tag: &str) -> Result<Arc<FileStack>, String> {
        self.fs_map
            .get(tag)
            .cloned()
            .ok_or_else(|| format!("FileStack tag not found: {tag}"))
    }

    pub fn sleep(sleep_s: f64) {
        std::thread::sleep(Duration::from_secs_f64(sleep_s));
    }

    pub fn log(&self, buf: &str) -> Result<(), String> {
        let log_stack = self.get_stack(Self::LOG)?;
        let ms = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_err(|e| e.to_string())?
            .as_millis();
        let tagged_buf = format!("{ms} {buf}\n");
        log_stack.push_string(&tagged_buf).map(|_| ())
    }

    pub fn list_filestacks(&self) -> Vec<(String, PathBuf)> {
        let mut out = self
            .fs_map
            .iter()
            .map(|(tag, stack)| (tag.clone(), stack.file_name().to_path_buf()))
            .collect::<Vec<_>>();
        out.sort();
        out
    }

    fn get_barrier_file_name(&self, step: &str, tag: &str, i: i32) -> PathBuf {
        PathBuf::from(format!(
            "{}_{}_{}_{}",
            self.barrier_file.display(),
            step,
            tag,
            i
        ))
    }

    fn clean(file_list: &mut Vec<PathBuf>) -> bool {
        for path in file_list.iter() {
            let _ = std::fs::remove_file(path);
        }
        file_list.clear();
        true
    }
}

impl Drop for Parallelizer {
    fn drop(&mut self) {
        if self.initialized {
            let _ = self.log("PARALLELIZER END");
            Self::clean(&mut self.continuous_cleanup_list);
            Self::clean(&mut self.final_cleanup_list);
        }
    }
}

pub type ThreadPoolDefaultTask = Arc<dyn Fn(ThreadPool, i64) + Send + std::marker::Sync + 'static>;

pub struct TaskSet {
    inner: Arc<TaskSetInner>,
    thread_pool: ThreadPool,
}

struct TaskSetInner {
    priority: usize,
    total: AtomicI64,
    finished: AtomicI64,
    cv: Condvar,
}

impl TaskSet {
    pub fn new(thread_pool: &ThreadPool, priority: i32) -> Self {
        Self {
            inner: Arc::new(TaskSetInner {
                priority: priority as usize,
                total: AtomicI64::new(0),
                finished: AtomicI64::new(0),
                cv: Condvar::new(),
            }),
            thread_pool: thread_pool.clone(),
        }
    }

    pub fn finish(&self) {
        self.inner.finish();
    }

    pub fn finished(&self) -> bool {
        self.inner.finished()
    }

    pub fn total(&self) -> i64 {
        self.inner.total.load(Ordering::Relaxed)
    }

    pub fn run(&self) {
        {
            let _lock = self.thread_pool.inner.state.lock().unwrap();
            if self.finished() {
                return;
            }
        }
        self.thread_pool
            .run_set_inner(Some(Arc::clone(&self.inner)));
    }

    pub fn enqueue<F>(&self, f: F)
    where
        F: FnOnce() + Send + 'static,
    {
        self.thread_pool.enqueue(self, f);
    }
}

impl TaskSetInner {
    fn finish(&self) {
        let finished = self.finished.fetch_add(1, Ordering::Relaxed) + 1;
        if finished == self.total.load(Ordering::Relaxed) {
            self.cv.notify_all();
        }
    }

    fn finished(&self) -> bool {
        self.total.load(Ordering::Relaxed) == self.finished.load(Ordering::Relaxed)
    }
}

struct ThreadPoolTask {
    f: Option<Box<dyn FnOnce() + Send + 'static>>,
    task_set: Option<Arc<TaskSetInner>>,
}

impl ThreadPoolTask {
    fn with_task_set<F>(f: F, task_set: Arc<TaskSetInner>) -> Self
    where
        F: FnOnce() + Send + 'static,
    {
        Self {
            f: Some(Box::new(f)),
            task_set: Some(task_set),
        }
    }

    fn run(mut self) {
        if let Some(f) = self.f.take() {
            f();
        }
        if let Some(task_set) = self.task_set {
            task_set.finish();
        }
    }
}

struct ThreadPoolState {
    tasks: [VecDeque<ThreadPoolTask>; ThreadPool::PRIORITY_COUNT],
}

struct ThreadPoolInner {
    pop_before_enqueue: bool,
    default_end: i64,
    default_count: i64,
    default_task: Option<ThreadPoolDefaultTask>,
    state: StdMutex<ThreadPoolState>,
    workers: StdMutex<Vec<JoinHandle<()>>>,
    heartbeat: StdMutex<Option<JoinHandle<()>>>,
    default_begin: AtomicI64,
    default_finished: AtomicI64,
    threads_finished: AtomicI64,
}

pub struct ThreadPool {
    inner: Arc<ThreadPoolInner>,
    join_on_drop: bool,
}

impl Clone for ThreadPool {
    fn clone(&self) -> Self {
        Self {
            inner: Arc::clone(&self.inner),
            join_on_drop: false,
        }
    }
}

impl Default for ThreadPool {
    fn default() -> Self {
        Self::new(None, 0, 0, false)
    }
}

impl ThreadPool {
    pub const PRIORITY_COUNT: usize = 2;

    pub fn new(
        default_task: Option<ThreadPoolDefaultTask>,
        default_begin: i64,
        default_end: i64,
        pop_before_enqueue: bool,
    ) -> Self {
        Self {
            inner: Arc::new(ThreadPoolInner {
                pop_before_enqueue,
                default_end,
                default_count: default_end - default_begin,
                default_task,
                state: StdMutex::new(ThreadPoolState {
                    tasks: [VecDeque::new(), VecDeque::new()],
                }),
                workers: StdMutex::new(Vec::new()),
                heartbeat: StdMutex::new(None),
                default_begin: AtomicI64::new(default_begin),
                default_finished: AtomicI64::new(0),
                threads_finished: AtomicI64::new(0),
            }),
            join_on_drop: true,
        }
    }

    pub fn enqueue<F>(&self, task_set: &TaskSet, f: F)
    where
        F: FnOnce() + Send + 'static,
    {
        if self.inner.pop_before_enqueue {
            loop {
                let task = {
                    let mut state = self.inner.state.lock().unwrap();
                    self.pop_task_locked(&mut state, Self::PRIORITY_COUNT - 1)
                };
                if let Some(task) = task {
                    task.run();
                } else {
                    break;
                }
            }
        }

        {
            let mut state = self.inner.state.lock().unwrap();
            task_set.inner.total.fetch_add(1, Ordering::Relaxed);
            state.tasks[task_set.inner.priority].push_back(ThreadPoolTask::with_task_set(
                f,
                Arc::clone(&task_set.inner),
            ));
            task_set.inner.cv.notify_one();
        }
    }

    pub fn run_set(&self, task_set: Option<&TaskSet>) {
        self.run_set_inner(task_set.map(|task_set| Arc::clone(&task_set.inner)));
    }

    fn run_set_inner(&self, task_set: Option<Arc<TaskSetInner>>) {
        loop {
            let task;

            if task_set.is_none() {
                if self.inner.default_finished.load(Ordering::Relaxed) >= self.inner.default_count {
                    self.inner.threads_finished.fetch_add(1, Ordering::Relaxed);
                    return;
                }
                task = {
                    let mut state = self.inner.state.lock().unwrap();
                    self.pop_task_locked(&mut state, Self::PRIORITY_COUNT - 1)
                };
                if task.is_none() {
                    let next = self.inner.default_begin.fetch_add(1, Ordering::Relaxed);
                    if next < self.inner.default_end {
                        if let Some(default_task) = &self.inner.default_task {
                            default_task(self.clone(), next);
                        }
                        self.inner.default_finished.fetch_add(1, Ordering::Relaxed);
                    }
                    continue;
                }
            } else {
                let task_set = task_set.as_ref().unwrap();
                let mut state = self.inner.state.lock().unwrap();
                while self.queue_empty_locked(&state, task_set.priority) && !task_set.finished() {
                    state = task_set.cv.wait(state).unwrap();
                }
                if task_set.finished() {
                    return;
                }
                task = self.pop_task_locked(&mut state, task_set.priority);
            }

            if let Some(task) = task {
                task.run();
            }
        }
    }

    pub fn run(&self, threads: i32, heartbeat: bool, task_set: Option<&TaskSet>) {
        let task_set_inner = task_set.map(|task_set| Arc::clone(&task_set.inner));
        {
            let mut workers = self.inner.workers.lock().unwrap();
            for _ in 0..threads {
                let pool = self.clone();
                let task_set_inner = task_set_inner.clone();
                workers.push(std::thread::spawn(move || {
                    pool.run_set_inner(task_set_inner);
                }));
            }
        }
        if heartbeat {
            let pool = self.clone();
            *self.inner.heartbeat.lock().unwrap() = Some(std::thread::spawn(move || {
                while pool.inner.default_finished.load(Ordering::Relaxed) < pool.inner.default_count
                {
                    eprintln!(
                        "Workers={}/{} begin = {} finished = {} queue={}/{}",
                        pool.inner.workers.lock().unwrap().len(),
                        pool.inner.threads_finished.load(Ordering::Relaxed),
                        pool.inner.default_begin.load(Ordering::Relaxed),
                        pool.inner.default_finished.load(Ordering::Relaxed),
                        pool.queue_len(0),
                        pool.queue_len(1)
                    );
                    std::thread::sleep(Duration::from_secs(1));
                }
            }));
        }
    }

    pub fn join(&self) {
        let workers = {
            let mut workers = self.inner.workers.lock().unwrap();
            std::mem::take(&mut *workers)
        };
        for worker in workers {
            let _ = worker.join();
        }
        if let Some(heartbeat) = self.inner.heartbeat.lock().unwrap().take() {
            let _ = heartbeat.join();
        }
    }

    pub fn queue_len(&self, priority: usize) -> i64 {
        let state = self.inner.state.lock().unwrap();
        state.tasks[priority].len() as i64
    }

    fn queue_empty_locked(&self, state: &ThreadPoolState, priority: usize) -> bool {
        (0..=priority).all(|i| state.tasks[i].is_empty())
    }

    fn pop_task_locked(
        &self,
        state: &mut ThreadPoolState,
        priority: usize,
    ) -> Option<ThreadPoolTask> {
        for i in 0..=priority {
            if let Some(task) = state.tasks[i].pop_front() {
                return Some(task);
            }
        }
        None
    }
}

impl Drop for ThreadPool {
    fn drop(&mut self) {
        if self.join_on_drop {
            self.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicI32, AtomicUsize};
    use std::sync::Arc;

    #[test]
    fn test_mutex_tags() {
        let sync = Mutex::<Sync>::new();
        sync.lock();
        sync.unlock();

        let async_mutex = Mutex::<Async>::new();
        async_mutex.lock();
        async_mutex.unlock();
    }

    #[test]
    fn test_counting_semaphore() {
        let sem = Arc::new(CountingSemaphore::<10>::new(0));
        let hit = Arc::new(AtomicI32::new(0));
        std::thread::scope(|scope| {
            let sem_wait = Arc::clone(&sem);
            let hit = Arc::clone(&hit);
            scope.spawn(move || {
                sem_wait.acquire();
                hit.fetch_add(1, Ordering::Relaxed);
            });
            sem.release(1);
        });
        assert_eq!(hit.load(Ordering::Relaxed), 1);
        assert_eq!(CountingSemaphore::<10>::max(), 10);
    }

    #[test]
    fn test_scheduled_thread_pool_auto() {
        let seen = AtomicUsize::new(0);
        scheduled_thread_pool_auto(4, 17, &|_, _| {
            seen.fetch_add(1, Ordering::Relaxed);
        });
        assert_eq!(seen.load(Ordering::Relaxed), 17);
    }

    #[test]
    fn test_launch_threads() {
        let seen = AtomicI32::new(0);
        launch_threads(3, &|| {
            seen.fetch_add(1, Ordering::Relaxed);
        });
        assert_eq!(seen.load(Ordering::Relaxed), 3);
    }

    #[test]
    fn test_simple_thread_pool_join_all() {
        let seen = Arc::new(AtomicUsize::new(0));
        let mut pool = SimpleThreadPool::new();
        for _ in 0..8 {
            let seen = Arc::clone(&seen);
            pool.spawn(move |_| {
                seen.fetch_add(1, Ordering::Relaxed);
                Ok(())
            });
        }
        pool.join_all().unwrap();
        assert_eq!(seen.load(Ordering::Relaxed), 8);
    }

    #[test]
    fn test_simple_thread_pool_join_ids_reports_first_exception() {
        let mut pool = SimpleThreadPool::new();
        let ok = pool.spawn(|_| Ok(()));
        let bad = pool.spawn(|_| Err("boom".to_string()));

        let err = pool.join_ids([ok, bad]).unwrap_err();
        assert_eq!(err, "boom");
        assert!(pool.stop().load(Ordering::Relaxed));
    }

    #[test]
    fn test_simple_thread_pool_join_one_and_missing_id() {
        let mut pool = SimpleThreadPool::new();
        let id = pool.spawn(|stop| {
            assert!(!stop.load(Ordering::Relaxed));
            Ok(())
        });
        pool.join(id).unwrap();

        let err = pool.join(id).unwrap_err();
        assert_eq!(err, "Thread ID not found in thread pool.");
        pool.request_stop();
        assert!(pool.stop().load(Ordering::Relaxed));
    }

    #[test]
    fn test_simple_thread_pool_spawn_method() {
        struct Worker {
            value: usize,
        }

        let worker = Arc::new(StdMutex::new(Worker { value: 0 }));
        let mut pool = SimpleThreadPool::new();
        pool.spawn_method(Arc::clone(&worker), |worker, stop| {
            assert!(!stop.load(Ordering::Relaxed));
            worker.value = 11;
            Ok(())
        });
        pool.join_all().unwrap();
        assert_eq!(worker.lock().unwrap().value, 11);
    }

    #[test]
    fn test_thread_pool_default_range() {
        let sum = Arc::new(AtomicI64::new(0));
        let seen = Arc::new(AtomicI64::new(0));
        let task: ThreadPoolDefaultTask = {
            let sum = Arc::clone(&sum);
            let seen = Arc::clone(&seen);
            Arc::new(move |_, i| {
                sum.fetch_add(i, Ordering::Relaxed);
                seen.fetch_add(1, Ordering::Relaxed);
            })
        };
        let pool = ThreadPool::new(Some(task), 0, 11, false);
        pool.run(3, false, None);
        pool.join();
        assert_eq!(seen.load(Ordering::Relaxed), 11);
        assert_eq!(sum.load(Ordering::Relaxed), 55);
    }

    #[test]
    fn test_thread_pool_task_set_priority_and_run() {
        let pool = ThreadPool::default();
        let low = TaskSet::new(&pool, 1);
        let high = TaskSet::new(&pool, 0);
        let order = Arc::new(StdMutex::new(Vec::new()));

        {
            let order = Arc::clone(&order);
            low.enqueue(move || order.lock().unwrap().push(1));
        }
        {
            let order = Arc::clone(&order);
            high.enqueue(move || order.lock().unwrap().push(0));
        }

        assert_eq!(low.total(), 1);
        assert_eq!(high.total(), 1);
        low.run();
        assert!(low.finished());
        assert!(high.finished());
        assert_eq!(*order.lock().unwrap(), vec![0, 1]);
    }

    #[test]
    fn test_thread_pool_run_task_set_with_workers() {
        let pool = ThreadPool::default();
        let task_set = TaskSet::new(&pool, 0);
        let seen = Arc::new(AtomicI64::new(0));

        for _ in 0..8 {
            let seen = Arc::clone(&seen);
            task_set.enqueue(move || {
                seen.fetch_add(1, Ordering::Relaxed);
            });
        }

        pool.run(2, false, Some(&task_set));
        pool.join();
        assert!(task_set.finished());
        assert_eq!(seen.load(Ordering::Relaxed), 8);
        assert_eq!(pool.queue_len(0), 0);
    }

    #[test]
    fn test_file_stack_string_operations() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-filestack-{}-{}.idx",
            std::process::id(),
            1
        ));
        let _ = std::fs::remove_file(&path);
        let stack = FileStack::with_max_line_length(&path, 2);
        assert_eq!(stack.get_max_line_length(), 8);
        assert_eq!(stack.size().unwrap(), 0);

        assert_eq!(stack.push_string("alpha").unwrap(), 6);
        assert_eq!(stack.push_string("beta\n").unwrap(), 5);
        let (written, size_after_push) = stack.push_string_size("gamma").unwrap();
        assert_eq!(written, 6);
        assert_eq!(size_after_push, 3);
        assert_eq!(stack.size().unwrap(), 3);
        assert_eq!(stack.top_string().unwrap().as_deref(), Some("gamma"));
        assert_eq!(
            stack.pop_string_size().unwrap(),
            (Some("gamma".to_string()), 2)
        );
        stack.remove("alpha").unwrap();
        assert_eq!(stack.size().unwrap(), 1);
        assert_eq!(stack.pop_string().unwrap().as_deref(), Some("beta"));
        assert_eq!(stack.pop_string().unwrap(), None);

        stack.push_string("ready").unwrap();
        assert!(stack.poll_query("ready", 0.0, 1).unwrap());
        assert!(stack.poll_size(1, 0.0, 1).unwrap());
        stack.clear().unwrap();
        assert_eq!(stack.size().unwrap(), 0);
        assert_eq!(stack.file_name(), path.as_path());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_file_stack_raw_file_operations() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-filestack-{}-{}.idx",
            std::process::id(),
            3
        ));
        let _ = std::fs::remove_file(&path);
        let stack = FileStack::new(&path);

        assert_eq!(stack.write(b"abcdef").unwrap(), 6);
        assert_eq!(stack.seek(2, SeekFrom::Start(0)).unwrap(), 2);
        let mut buf = [0u8; 3];
        assert_eq!(stack.read(&mut buf).unwrap(), 3);
        assert_eq!(&buf, b"cde");
        assert_eq!(stack.seek(-2, SeekFrom::End(0)).unwrap(), 4);
        assert_eq!(stack.write(b"XY").unwrap(), 2);
        stack.truncate(5).unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"abcdX");

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_file_stack_i64_and_atomic() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-filestack-{}-{}.idx",
            std::process::id(),
            2
        ));
        let _ = std::fs::remove_file(&path);
        let stack = FileStack::new(&path);
        assert_eq!(stack.pop_i64().unwrap(), -1);
        assert_eq!(stack.push_i64(3).unwrap(), 2);
        assert_eq!(stack.fetch_add(4).unwrap(), 3);
        assert_eq!(stack.top_i64().unwrap(), 7);
        let mut value = 0;
        assert_eq!(stack.top_i64_into(&mut value).unwrap(), 7);
        assert_eq!(value, 7);
        assert_eq!(stack.pop_i64_into(&mut value).unwrap(), 7);
        assert_eq!(value, 7);

        let atomic = Atomic::new(&path);
        assert_eq!(atomic.get().unwrap(), 0);
        assert_eq!(atomic.fetch_add_one().unwrap(), 0);
        assert_eq!(atomic.fetch_add(5).unwrap(), 1);
        assert_eq!(atomic.get().unwrap(), 6);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_parallelizer_init_register_barrier_and_log() {
        let root = std::env::temp_dir().join(format!(
            "diamond-rs-parallelizer-{}-{}",
            std::process::id(),
            1
        ));
        let _ = std::fs::remove_dir_all(&root);

        {
            let mut p = Parallelizer::new();
            p.init(&root).unwrap();
            assert!(p.is_master());
            assert_eq!(p.get_rank(), 0);
            assert!(p.get_id().contains('_'));
            assert!(p.get_work_directory().ends_with("parallelizer"));

            let stacks = p.list_filestacks();
            assert!(stacks.iter().any(|(tag, _)| tag == Parallelizer::LOG));
            assert!(stacks.iter().any(|(tag, _)| tag == Parallelizer::COMMAND));
            assert!(stacks.iter().any(|(tag, _)| tag == Parallelizer::WORKERS));
            assert!(stacks.iter().any(|(tag, _)| tag == Parallelizer::REGISTER));

            assert!(p.register_workers(0.0).unwrap());
            assert_eq!(p.get_n_registered(), 1);
            assert_eq!(
                p.get_stack(Parallelizer::WORKERS).unwrap().size().unwrap(),
                1
            );
            assert!(p.barrier("unit").unwrap());
            assert_eq!(p.i_barrier, 1);
            p.log("custom message").unwrap();

            let log_top = p
                .get_stack(Parallelizer::LOG)
                .unwrap()
                .top_string()
                .unwrap();
            assert!(log_top.unwrap().contains("custom message"));

            let custom_path = root.join("custom-stack");
            assert!(p.create_stack_from_file("custom", &custom_path).unwrap());
            assert!(p.get_stack("custom").is_ok());
            assert!(!p.create_stack("custom", "").unwrap());
            assert!(p.delete_stack("custom"));
            assert!(p.get_stack("custom").is_err());
        }

        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn test_parallelizer_get_singleton_and_uninitialized_barrier() {
        let singleton = Parallelizer::get();
        {
            let p = singleton.lock().unwrap();
            assert!(p.get_work_directory().ends_with("parallelizer"));
        }

        let mut p = Parallelizer::new();
        assert!(!p.barrier("before-init").unwrap());
        p.clear();
    }
}
