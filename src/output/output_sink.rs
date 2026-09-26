use crate::util::data_structures::ReorderQueue;
use crate::util::misc::megabytes;
use std::io::{self, Write};
use std::time::Duration;

/// Rust translation of C++ `OutputWriter`.
#[derive(Debug)]
pub struct OutputWriter<W> {
    pub file: W,
    pub first: bool,
    pub sep: u8,
}

impl<W> OutputWriter<W> {
    /// Matches C++ `OutputWriter::OutputWriter(file, sep, first)`.
    pub fn new(file: W, sep: u8, first: bool) -> Self {
        Self { file, first, sep }
    }

    /// C++ constructor defaults: `sep = '\0'`, `first = true`.
    pub fn from_writer(file: W) -> Self {
        Self::new(file, 0, true)
    }

    pub fn into_inner(self) -> W {
        self.file
    }
}

impl<W: Write> OutputWriter<W> {
    /// Matches C++ `OutputWriter::operator()(buf)`.
    pub fn consume(&mut self, buf: &[u8]) -> io::Result<()> {
        if !self.first && self.sep != 0 {
            self.file.write_all(&[self.sep])?;
        }
        self.file.write_all(buf)?;
        self.first = false;
        Ok(())
    }
}

/// `ReorderQueue` surface used by C++ `heartbeat_worker`.
pub trait HeartbeatOutputSink {
    fn next(&mut self) -> usize;
    fn size(&self) -> usize;
    fn max_size(&self) -> usize;
}

impl<T> HeartbeatOutputSink for ReorderQueue<T> {
    fn next(&mut self) -> usize {
        ReorderQueue::next(self)
    }

    fn size(&self) -> usize {
        ReorderQueue::size(self)
    }

    fn max_size(&self) -> usize {
        ReorderQueue::max_size(self)
    }
}

/// `Search::Config` surface used by C++ `heartbeat_worker`.
pub trait HeartbeatConfig {
    fn query_title(&self, query: usize) -> &str;
    fn queue_len(&self, queue: usize) -> usize;
}

/// Matches C++ `heartbeat_worker(qend, output_sink, cfg, verbose_stream)`.
pub fn heartbeat_worker<S, C, W, F>(
    qend: usize,
    output_sink: &mut S,
    cfg: &C,
    verbose_stream: &mut W,
    mut sleep_for: F,
) -> io::Result<()>
where
    S: HeartbeatOutputSink,
    C: HeartbeatConfig,
    W: Write,
    F: FnMut(Duration),
{
    const INTERVAL: i32 = 100;
    let mut n = 0;
    loop {
        let next = output_sink.next();
        if next >= qend {
            break;
        }
        if n == INTERVAL {
            let title = cfg.query_title(next);
            let title = title.split(' ').next().unwrap_or(title);
            writeln!(
                verbose_stream,
                "Queries={} size={} max_size={} next={} queue={}/{}",
                next,
                megabytes(output_sink.size()),
                megabytes(output_sink.max_size()),
                title,
                cfg.queue_len(0),
                cfg.queue_len(1)
            )?;
            n = 0;
        } else {
            n += 1;
        }
        sleep_for(Duration::from_millis(10));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::data_structures::AllocSize;

    #[derive(Debug)]
    struct TestSink {
        next_values: Vec<usize>,
        size: usize,
        max_size: usize,
    }

    impl HeartbeatOutputSink for TestSink {
        fn next(&mut self) -> usize {
            self.next_values.remove(0)
        }

        fn size(&self) -> usize {
            self.size
        }

        fn max_size(&self) -> usize {
            self.max_size
        }
    }

    #[derive(Debug)]
    struct TestConfig {
        titles: Vec<String>,
        queues: [usize; 2],
    }

    impl HeartbeatConfig for TestConfig {
        fn query_title(&self, query: usize) -> &str {
            &self.titles[query]
        }

        fn queue_len(&self, queue: usize) -> usize {
            self.queues[queue]
        }
    }

    #[test]
    fn test_output_writer_consumes_separator_after_first_buffer() {
        let mut writer = OutputWriter::new(Vec::new(), b'\n', true);
        writer.consume(b"one").unwrap();
        writer.consume(b"two").unwrap();
        writer.consume(b"three").unwrap();
        assert_eq!(writer.file, b"one\ntwo\nthree");
        assert!(!writer.first);
    }

    #[test]
    fn test_output_writer_omits_zero_separator() {
        let mut writer = OutputWriter::new(Vec::new(), 0, false);
        writer.consume(b"one").unwrap();
        writer.consume(b"two").unwrap();
        assert_eq!(writer.file, b"onetwo");
    }

    #[derive(Debug)]
    struct Buffer(Vec<u8>);

    impl AllocSize for Buffer {
        fn alloc_size(&self) -> usize {
            self.0.capacity()
        }
    }

    #[test]
    fn test_reorder_queue_and_output_writer_preserve_cpp_binary_order() {
        let mut queue = ReorderQueue::new(1);
        let mut writer = OutputWriter::new(Vec::new(), b'|', true);

        queue.push(2, Some(Buffer(b"two".to_vec())), |buf| {
            writer.consume(&buf.0).unwrap()
        });
        assert_eq!(queue.next(), 1);
        assert!(queue.size() >= 3);
        assert!(queue.max_size() >= 3);

        queue.push(1, Some(Buffer(b"one".to_vec())), |buf| {
            writer.consume(&buf.0).unwrap()
        });
        assert_eq!(queue.next(), 3);
        assert_eq!(writer.file, b"one|two");

        queue.push(4, Some(Buffer(b"four".to_vec())), |buf| {
            writer.consume(&buf.0).unwrap()
        });
        queue.push(3, None, |buf: Buffer| writer.consume(&buf.0).unwrap());
        assert_eq!(queue.next(), 5);
        assert_eq!(writer.file, b"one|two|four");
    }

    #[test]
    fn test_output_writer_cpp_constructor_defaults() {
        let mut writer = OutputWriter::from_writer(Vec::new());
        writer.consume(b"a").unwrap();
        writer.consume(b"b").unwrap();
        assert_eq!(writer.into_inner(), b"ab");
    }

    #[test]
    fn test_heartbeat_worker_logs_every_cpp_interval_and_stops_at_qend() {
        let mut next_values = vec![0usize; 101];
        next_values.push(3);
        let mut sink = TestSink {
            next_values,
            size: 2 << 20,
            max_size: 3 << 20,
        };
        let cfg = TestConfig {
            titles: vec!["query zero".to_string()],
            queues: [4, 5],
        };
        let mut out = Vec::new();
        let mut sleeps = Vec::new();
        heartbeat_worker(3, &mut sink, &cfg, &mut out, |d| sleeps.push(d)).unwrap();

        assert_eq!(
            std::str::from_utf8(&out).unwrap(),
            "Queries=0 size=2 max_size=3 next=query queue=4/5\n"
        );
        assert_eq!(sleeps.len(), 101);
        assert!(sleeps.iter().all(|d| *d == Duration::from_millis(10)));
    }
}
