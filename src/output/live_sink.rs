//! Actual ordered output storage and Source heartbeat over that live storage.
use crate::util::data_structures::{AllocSize, ReorderQueue};
use std::cell::RefCell;
use std::io::{self, Write};
use std::sync::{
    atomic::{AtomicBool, AtomicUsize, Ordering},
    Arc, Mutex,
};

#[derive(Default)]
pub(crate) struct TaskMonitor {
    dp_queued: AtomicUsize,
}
thread_local! { static CURRENT: RefCell<Option<Arc<TaskMonitor>>> = const { RefCell::new(None) }; }
pub(crate) struct TaskScope(Option<Arc<TaskMonitor>>);
struct StopHeartbeat<'a>(&'a AtomicBool);
impl Drop for StopHeartbeat<'_> {
    fn drop(&mut self) {
        self.0.store(true, Ordering::Release);
    }
}
impl Drop for TaskScope {
    fn drop(&mut self) {
        CURRENT.with(|slot| *slot.borrow_mut() = self.0.take());
    }
}
pub(crate) fn current_monitor() -> Option<Arc<TaskMonitor>> {
    CURRENT.with(|slot| slot.borrow().clone())
}
pub(crate) fn enter(monitor: Option<Arc<TaskMonitor>>) -> TaskScope {
    TaskScope(CURRENT.with(|slot| slot.replace(monitor)))
}
impl TaskMonitor {
    pub(crate) fn enqueue_dp(&self, count: usize) {
        self.dp_queued.fetch_add(count, Ordering::Relaxed);
    }
    pub(crate) fn start_dp(&self) {
        self.dp_queued.fetch_sub(1, Ordering::Relaxed);
    }
}

struct Buffer<T> {
    value: T,
    bytes: usize,
}
impl<T> AllocSize for Buffer<T> {
    fn alloc_size(&self) -> usize {
        self.bytes
    }
}
struct State<T, F> {
    queue: ReorderQueue<Buffer<T>>,
    consume: F,
    error: Option<io::Error>,
}
pub(crate) struct LiveSink<T, F> {
    state: Mutex<State<T, F>>,
    monitor: Option<Arc<TaskMonitor>>,
    stopped: AtomicBool,
}
impl<T: Send, F: FnMut(T) -> io::Result<()> + Send> LiveSink<T, F> {
    pub(crate) fn new(begin: usize, heartbeat: bool, consume: F) -> Self {
        Self {
            state: Mutex::new(State {
                queue: ReorderQueue::new(begin),
                consume,
                error: None,
            }),
            monitor: heartbeat.then(|| Arc::new(TaskMonitor::default())),
            stopped: AtomicBool::new(false),
        }
    }
    pub(crate) fn task_scope(&self) -> TaskScope {
        enter(self.monitor.clone())
    }
    pub(crate) fn push(&self, query: usize, value: T, bytes: usize) -> io::Result<()> {
        let mut state = self
            .state
            .lock()
            .map_err(|_| io::Error::other("Output sink lock poisoned."))?;
        let State {
            queue,
            consume,
            error,
        } = &mut *state;
        queue.push(query, Some(Buffer { value, bytes }), |buffer| {
            if error.is_none() {
                if let Err(failure) = consume(buffer.value) {
                    *error = Some(failure);
                }
            }
        });
        Ok(())
    }
    pub(crate) fn execute<P>(&self, begin: usize, titles: &[&str], process: P) -> io::Result<()>
    where
        P: FnOnce() -> io::Result<()>,
    {
        let titles: Vec<&[u8]> = titles.iter().map(|title| title.as_bytes()).collect();
        self.execute_bytes(begin, &titles, process)
    }
    pub(crate) fn execute_bytes<P>(
        &self,
        begin: usize,
        titles: &[&[u8]],
        process: P,
    ) -> io::Result<()>
    where
        P: FnOnce() -> io::Result<()>,
    {
        std::thread::scope(|scope| {
            let heartbeat = self.monitor.as_ref().map(|monitor| {
                scope.spawn(move || {
                    let mut sink = HeartbeatSink {
                        sink: self,
                        end: begin + titles.len(),
                    };
                    let cfg = HeartbeatConfig {
                        begin,
                        titles,
                        monitor,
                    };
                    crate::output::output_sink::heartbeat_worker_bytes(
                        begin + titles.len(),
                        &mut sink,
                        &cfg,
                        &mut VerboseOutput,
                        std::thread::sleep,
                    )
                })
            });
            let stop = StopHeartbeat(&self.stopped);
            let result = process();
            drop(stop);
            let logged = match heartbeat {
                Some(thread) => thread
                    .join()
                    .map_err(|_| io::Error::other("Heartbeat thread panicked."))?,
                None => Ok(()),
            };
            result?;
            logged
        })
    }
    pub(crate) fn finish(self, end: usize) -> io::Result<()> {
        let state = self
            .state
            .into_inner()
            .map_err(|_| io::Error::other("Output sink lock poisoned."))?;
        if let Some(error) = state.error {
            return Err(error);
        }
        if state.queue.next() != end {
            return Err(io::Error::other("Missing ordered query output."));
        }
        Ok(())
    }
}
struct HeartbeatSink<'a, T, F> {
    sink: &'a LiveSink<T, F>,
    end: usize,
}
impl<T, F> crate::output::output_sink::HeartbeatOutputSink for HeartbeatSink<'_, T, F> {
    fn next(&mut self) -> usize {
        if self.sink.stopped.load(Ordering::Acquire) {
            self.end
        } else {
            self.sink
                .state
                .lock()
                .map_or(self.end, |state| state.queue.next())
        }
    }
    fn size(&self) -> usize {
        self.sink.state.lock().map_or(0, |state| state.queue.size())
    }
    fn max_size(&self) -> usize {
        self.sink
            .state
            .lock()
            .map_or(0, |state| state.queue.max_size())
    }
}
struct HeartbeatConfig<'a> {
    begin: usize,
    titles: &'a [&'a [u8]],
    monitor: &'a TaskMonitor,
}
impl crate::output::output_sink::HeartbeatByteConfig for HeartbeatConfig<'_> {
    fn query_title_bytes(&self, query: usize) -> &[u8] {
        self.titles[query - self.begin]
    }
    fn queue_len(&self, priority: usize) -> usize {
        // Default query jobs are atomically claimed, not put into priority0.
        if priority == 1 {
            self.monitor.dp_queued.load(Ordering::Relaxed)
        } else {
            0
        }
    }
}
struct VerboseOutput;
impl Write for VerboseOutput {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        crate::util::log_stream::verbose_stream()
            .lock()
            .map_err(|_| io::Error::other("Verbose log lock poisoned."))?
            .write_bytes(bytes)?;
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
