//! Rust counterpart of the upstream queue stress test in
//! `diamond/src/test/queue.cpp`.

use crate::util::data_structures::Queue;
use crate::util::parallel::filestack::FileStack;
use std::path::Path;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Arc;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QueueStressTestResult {
    pub passed: bool,
    pub items_sent: usize,
    pub items_received: usize,
    pub expected_checksum: u64,
    pub received_checksum: u64,
}

fn test_many_producers_one_consumer(
    thread_count: usize,
    items_per_producer: usize,
) -> QueueStressTestResult {
    let producer_count = thread_count.saturating_sub(1);
    if producer_count < 1 {
        return QueueStressTestResult {
            passed: false,
            items_sent: 0,
            items_received: 0,
            expected_checksum: 0,
            received_checksum: 0,
        };
    }
    let queue = Arc::new(Queue::new(1024, producer_count as i32, 1, -1_i64));
    let total_sent = Arc::new(AtomicUsize::new(0));
    let total_received = Arc::new(AtomicUsize::new(0));
    let sent_checksum = Arc::new(AtomicU64::new(0));
    let received_checksum = Arc::new(AtomicU64::new(0));
    let mut producers = Vec::with_capacity(producer_count);
    for producer in 0..producer_count {
        let queue = Arc::clone(&queue);
        let total_sent = Arc::clone(&total_sent);
        let sent_checksum = Arc::clone(&sent_checksum);
        producers.push(std::thread::spawn(move || {
            let mut checksum = 0_u64;
            for item in 0..items_per_producer {
                let value = (producer * items_per_producer + item) as i64;
                queue.enqueue(value);
                checksum = checksum.wrapping_add(value as u64);
            }
            total_sent.fetch_add(items_per_producer, Ordering::Relaxed);
            sent_checksum.fetch_add(checksum, Ordering::Relaxed);
            queue.enqueue(-1);
        }));
    }
    let consumer_queue = Arc::clone(&queue);
    let consumer_count = Arc::clone(&total_received);
    let consumer_checksum = Arc::clone(&received_checksum);
    let consumer = std::thread::spawn(move || {
        let mut count = 0;
        let mut checksum = 0_u64;
        while let Some(value) = consumer_queue.wait_and_dequeue() {
            count += 1;
            checksum = checksum.wrapping_add(value as u64);
        }
        consumer_count.fetch_add(count, Ordering::Relaxed);
        consumer_checksum.fetch_add(checksum, Ordering::Relaxed);
    });
    for producer in producers {
        producer.join().unwrap();
    }
    consumer.join().unwrap();
    let expected_count = producer_count * items_per_producer;
    let sent = total_sent.load(Ordering::Relaxed);
    let received = total_received.load(Ordering::Relaxed);
    let sent_sum = sent_checksum.load(Ordering::Relaxed);
    let received_sum = received_checksum.load(Ordering::Relaxed);
    QueueStressTestResult {
        passed: sent == expected_count && received == expected_count && sent_sum == received_sum,
        items_sent: sent,
        items_received: received,
        expected_checksum: sent_sum,
        received_checksum: received_sum,
    }
}

fn test_one_producer_many_consumers(
    thread_count: usize,
    total_items: usize,
) -> QueueStressTestResult {
    let consumer_count = thread_count.saturating_sub(1);
    if consumer_count < 1 {
        return QueueStressTestResult {
            passed: false,
            items_sent: 0,
            items_received: 0,
            expected_checksum: 0,
            received_checksum: 0,
        };
    }
    let queue = Arc::new(Queue::new(1024, 1, consumer_count as i32, -1_i64));
    let total_received = Arc::new(AtomicUsize::new(0));
    let received_checksum = Arc::new(AtomicU64::new(0));
    let producer_queue = Arc::clone(&queue);
    let producer = std::thread::spawn(move || {
        let mut checksum = 0_u64;
        for item in 0..total_items {
            producer_queue.enqueue(item as i64);
            checksum = checksum.wrapping_add(item as u64);
        }
        producer_queue.close();
        checksum
    });
    let mut consumers = Vec::with_capacity(consumer_count);
    for _ in 0..consumer_count {
        let queue = Arc::clone(&queue);
        let count = Arc::clone(&total_received);
        let checksum = Arc::clone(&received_checksum);
        consumers.push(std::thread::spawn(move || {
            let mut local_count = 0;
            let mut local_checksum = 0_u64;
            while let Some(value) = queue.wait_and_dequeue() {
                local_count += 1;
                local_checksum = local_checksum.wrapping_add(value as u64);
            }
            count.fetch_add(local_count, Ordering::Relaxed);
            checksum.fetch_add(local_checksum, Ordering::Relaxed);
        }));
    }
    let expected_checksum = producer.join().unwrap();
    for consumer in consumers {
        consumer.join().unwrap();
    }
    let received = total_received.load(Ordering::Relaxed);
    let received_sum = received_checksum.load(Ordering::Relaxed);
    QueueStressTestResult {
        passed: received == total_items && received_sum == expected_checksum,
        items_sent: total_items,
        items_received: received,
        expected_checksum,
        received_checksum: received_sum,
    }
}

/// Runs both upstream queue stress scenarios. A thread count is explicit so
/// tests remain deterministic on single-core CI hosts.
pub fn run_queue_stress_test_with_threads(thread_count: usize) -> i32 {
    let first = test_many_producers_one_consumer(thread_count, 300);
    let second = test_one_producer_many_consumers(thread_count, 10_000);
    (!first.passed) as i32 + (!second.passed) as i32
}

/// Translation of the upstream no-argument entry point.
pub fn run_queue_stress_test() -> i32 {
    let threads = std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1);
    println!("Queue Stress Test");
    println!("=================");
    println!("Hardware threads: {threads}");
    println!();
    if threads < 2 {
        println!("Error: Need at least 2 threads for stress test");
        return 1;
    }

    println!("Test 1: Many producers ({}), one consumer", threads - 1);
    println!("  Items per producer: 300");
    println!("  Total items: {}", (threads - 1) * 300);
    let first = test_many_producers_one_consumer(threads, 300);
    println!("  Items sent: {}", first.items_sent);
    println!("  Items received: {}", first.items_received);
    println!("  Expected checksum: {}", first.expected_checksum);
    println!("  Received checksum: {}", first.received_checksum);
    println!(
        "  Result: {}",
        if first.passed { "PASSED" } else { "FAILED" }
    );
    println!();

    println!("Test 2: One producer, many consumers ({})", threads - 1);
    println!("  Total items: 10000");
    let second = test_one_producer_many_consumers(threads, 10_000);
    println!("  Items sent: {}", second.items_sent);
    println!("  Items received: {}", second.items_received);
    println!("  Expected checksum: {}", second.expected_checksum);
    println!("  Received checksum: {}", second.received_checksum);
    println!(
        "  Result: {}",
        if second.passed { "PASSED" } else { "FAILED" }
    );
    println!();

    let failures = (!first.passed) as i32 + (!second.passed) as i32;
    println!("=================");
    println!("Tests passed: {}/2", 2 - failures);
    failures
}

/// File-stack concurrency exercise with the source process globals made
/// explicit. Values are deterministic; the upstream random values are not an
/// assertion target and only serve to vary line contents.
pub fn filestack(path: &Path, thread_count: usize) -> Result<(), String> {
    let stack = Arc::new(FileStack::new(path));
    stack.clear()?;
    let mut workers = Vec::with_capacity(thread_count);
    for thread_id in 0..thread_count {
        let stack = Arc::clone(&stack);
        workers.push(std::thread::spawn(move || -> Result<(), String> {
            for item in 0..100 {
                stack.push_string(&format!(
                    "{thread_id}\t{item}\t{}\t{}\n",
                    thread_id * 100 + item,
                    item * 17
                ))?;
            }
            Ok(())
        }));
    }
    for worker in workers {
        worker
            .join()
            .map_err(|_| "filestack worker panicked".to_string())??;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_queue_stress_directions_preserve_counts_and_checksums() {
        let many_to_one = test_many_producers_one_consumer(4, 2_000);
        assert!(many_to_one.passed, "{many_to_one:?}");
        let one_to_many = test_one_producer_many_consumers(4, 10_000);
        assert!(one_to_many.passed, "{one_to_many:?}");
        assert_eq!(run_queue_stress_test_with_threads(4), 0);
    }

    #[test]
    fn filestack_workers_write_every_record() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-queue-filestack-{}-{}.tsv",
            std::process::id(),
            std::thread::current().name().unwrap_or("test")
        ));
        filestack(&path, 4).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(text.lines().count(), 400);
        std::fs::remove_file(path).unwrap();
    }
}
