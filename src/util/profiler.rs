use crate::util::log_stream::{message_stream, TaskTimer};
use std::collections::HashMap;
use std::sync::{LazyLock, Mutex};

pub static TIMES: LazyLock<Mutex<HashMap<String, u64>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

#[derive(Debug)]
pub struct Profiler {
    pub timer: TaskTimer,
    key: Option<String>,
}

impl Profiler {
    pub fn new(key: &str) -> Self {
        Self {
            timer: TaskTimer::new(u32::MAX),
            key: Some(key.to_string()),
        }
    }

    pub fn finish(&mut self) {
        let Some(key) = self.key.take() else {
            return;
        };
        let elapsed = self.timer.nanoseconds().max(0) as u64;
        let mut times = TIMES.lock().unwrap();
        *times.entry(key).or_insert(0) += elapsed;
    }

    pub fn print(n: usize) {
        let times = TIMES.lock().unwrap();
        if n == 0 {
            return;
        }
        if let Ok(mut stream) = message_stream().lock() {
            for (key, value) in times.iter() {
                let micros = *value as f64 / n as f64 / 1e3;
                let _ = stream
                    .write(key)
                    .and_then(|s| s.write(": "))
                    .and_then(|s| s.write(micros))
                    .and_then(|s| s.write(" micros"))
                    .and_then(|s| s.endl());
            }
        }
    }

    pub fn times() -> &'static Mutex<HashMap<String, u64>> {
        &TIMES
    }
}

impl Drop for Profiler {
    fn drop(&mut self) {
        self.finish();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_profiler_finish_accumulates_once() {
        TIMES.lock().unwrap().clear();
        let mut profiler = Profiler::new("x");
        profiler.finish();
        profiler.finish();
        assert!(TIMES.lock().unwrap().get("x").copied().unwrap() > 0);
    }

    #[test]
    fn test_profiler_drop_and_print() {
        TIMES.lock().unwrap().clear();
        {
            let _profiler = Profiler::new("drop");
        }
        assert!(TIMES.lock().unwrap().contains_key("drop"));
        Profiler::print(0);
    }
}
