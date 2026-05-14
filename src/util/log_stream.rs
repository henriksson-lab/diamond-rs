use std::fmt::Display;
use std::fs::OpenOptions;
use std::io::{self, Write};
use std::sync::{LazyLock, Mutex};
use std::time::Instant;

pub static MESSAGE_STREAM: LazyLock<Mutex<MessageStream>> =
    LazyLock::new(|| Mutex::new(MessageStream::new(true, false)));
pub static VERBOSE_STREAM: LazyLock<Mutex<MessageStream>> =
    LazyLock::new(|| Mutex::new(MessageStream::new(true, false)));
pub static LOG_STREAM: LazyLock<Mutex<MessageStream>> =
    LazyLock::new(|| Mutex::new(MessageStream::new(true, false)));

#[derive(Debug, Clone)]
pub struct MessageStream {
    to_cout: bool,
    to_file: bool,
}

impl Default for MessageStream {
    fn default() -> Self {
        Self::new(true, false)
    }
}

impl MessageStream {
    pub fn new(to_cout: bool, to_file: bool) -> Self {
        Self { to_cout, to_file }
    }

    pub fn write<T: Display>(&mut self, x: T) -> io::Result<&mut Self> {
        let s = x.to_string();
        if self.to_cout {
            eprint!("{s}");
            io::stderr().flush()?;
        }
        if self.to_file {
            let mut f = OpenOptions::new()
                .create(true)
                .append(true)
                .open("diamond.log")?;
            f.write_all(s.as_bytes())?;
        }
        Ok(self)
    }

    pub fn endl(&mut self) -> io::Result<&mut Self> {
        self.write('\n')
    }

    pub fn flush(&mut self) -> io::Result<&mut Self> {
        if self.to_cout {
            io::stderr().flush()?;
        }
        Ok(self)
    }

    pub fn to_cout(&self) -> bool {
        self.to_cout
    }

    pub fn to_file(&self) -> bool {
        self.to_file
    }
}

pub fn message_stream() -> &'static Mutex<MessageStream> {
    &MESSAGE_STREAM
}

pub fn verbose_stream() -> &'static Mutex<MessageStream> {
    &VERBOSE_STREAM
}

pub fn log_stream() -> &'static Mutex<MessageStream> {
    &LOG_STREAM
}

#[derive(Debug)]
pub struct TaskTimer {
    level: u32,
    msg: Option<String>,
    stream: &'static Mutex<MessageStream>,
    t: Instant,
}

impl Default for TaskTimer {
    fn default() -> Self {
        Self::new(1)
    }
}

impl TaskTimer {
    pub fn new(level: u32) -> Self {
        let stream = Self::get_stream(level);
        let mut out = Self {
            level,
            msg: None,
            stream,
            t: Instant::now(),
        };
        out.start(None);
        out
    }

    pub fn with_stream(stream: &'static Mutex<MessageStream>, level: u32) -> Self {
        let mut out = Self {
            level,
            msg: None,
            stream,
            t: Instant::now(),
        };
        out.start(None);
        out
    }

    pub fn with_msg(msg: &str, level: u32) -> Self {
        Self::with_msg_stream(msg, Self::get_stream(level), level)
    }

    pub fn with_msg_stream(msg: &str, stream: &'static Mutex<MessageStream>, level: u32) -> Self {
        let mut out = Self {
            level,
            msg: None,
            stream,
            t: Instant::now(),
        };
        out.start(Some(msg));
        out
    }

    pub fn go(&mut self, msg: Option<&str>) {
        self.finish();
        self.start(msg);
    }

    pub fn go_string(&mut self, s: &str) {
        self.go(Some(s));
    }

    pub fn finish(&mut self) {
        if self.msg.is_none() || self.level == u32::MAX {
            return;
        }
        if let Ok(mut stream) = self.stream.lock() {
            let _ = stream
                .write(" [")
                .and_then(|s| s.write(self.get()))
                .and_then(|s| s.write("s]"))
                .and_then(|s| s.endl());
        }
        self.msg = None;
    }

    pub fn get(&self) -> f64 {
        self.t.elapsed().as_millis() as f64 / 1000.0
    }

    pub fn seconds(&self) -> i64 {
        self.t.elapsed().as_secs() as i64
    }

    pub fn milliseconds(&self) -> i64 {
        self.t.elapsed().as_millis() as i64
    }

    pub fn microseconds(&self) -> i64 {
        self.t.elapsed().as_micros() as i64
    }

    pub fn nanoseconds(&self) -> i64 {
        self.t.elapsed().as_nanos() as i64
    }

    fn start(&mut self, msg: Option<&str>) {
        self.t = Instant::now();
        if self.level == u32::MAX {
            return;
        }
        let Some(msg) = msg else {
            return;
        };
        if let Ok(mut stream) = self.stream.lock() {
            let _ = stream
                .write(msg)
                .and_then(|s| s.write("... "))
                .and_then(|s| s.flush());
        }
        self.msg = Some(msg.to_string());
    }

    fn get_stream(level: u32) -> &'static Mutex<MessageStream> {
        match level {
            1 => &MESSAGE_STREAM,
            2 => &VERBOSE_STREAM,
            3 => &LOG_STREAM,
            _ => &MESSAGE_STREAM,
        }
    }
}

impl Drop for TaskTimer {
    fn drop(&mut self) {
        self.finish();
    }
}

pub fn exit_with_error<E: Display>(e: E) -> ! {
    if let Ok(mut stream) = LOG_STREAM.lock() {
        let _ = stream
            .write("Error: ")
            .and_then(|s| s.write(e))
            .and_then(|s| s.endl());
    }
    std::process::exit(1);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_message_stream_write_and_flags() {
        let mut stream = MessageStream::new(false, false);
        assert!(!stream.to_cout());
        assert!(!stream.to_file());
        stream.write("abc").unwrap().endl().unwrap();
        stream.flush().unwrap();
    }

    #[test]
    fn test_task_timer_methods() {
        let mut timer = TaskTimer::new(u32::MAX);
        assert!(timer.get() >= 0.0);
        assert!(timer.seconds() >= 0);
        assert!(timer.milliseconds() >= 0);
        assert!(timer.microseconds() >= 0);
        assert!(timer.nanoseconds() >= 0);
        timer.go(Some("suppressed"));
        timer.finish();
    }
}
