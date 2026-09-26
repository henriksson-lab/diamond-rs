use std::collections::HashMap;
use std::collections::VecDeque;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicI64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex as StdMutex};
use std::thread::{JoinHandle, ThreadId};
use std::time::Duration;

pub mod filestack;
/// Serialization and path helpers from `util/parallel/multiprocessing.{h,cpp}`.
pub mod multiprocessing;
pub mod parallelizer;

pub use filestack::FileStack;
pub use parallelizer::Parallelizer;

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
    use std::io::SeekFrom;
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
