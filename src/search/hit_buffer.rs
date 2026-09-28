//! Buffered seed-hit storage translated from `search/hit_buffer.cpp`.
//!
//! The C++ implementation overlaps disk I/O with search workers. Rust keeps
//! the same packet format and load/retrieve state machine, including loading
//! the next disk bin on a background thread while the previous bin is used.

use super::hit::Hit;
use crate::basic::seed::SeedOffset;
use crate::basic::value::BlockId;
use std::fs::{self, File};
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread::{self, JoinHandle};

static TEMP_ID: AtomicU64 = AtomicU64::new(0);
const DISK_BUFFER_SIZE: usize = 65_535;
const MEMORY_BUFFER_SIZE: usize = 8_192;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HitBufferMode {
    Disk,
    Memory,
    SwipeAll,
}

/// Query-local hit representation used by the native protein alignment path.
/// The query id is implicit in the containing bin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CompactHit {
    pub subject: u64,
    pub seed_offset: SeedOffset,
    pub score: u16,
}

/// Rust counterpart of C++ `Search::HitBuffer`.
pub struct HitBuffer {
    bins: Vec<Vec<Hit>>,
    count: Vec<usize>,
    key_partition: Vec<u32>,
    long_subject_offsets: bool,
    query_contexts: u32,
    max_query: u32,
    max_target: u64,
    mode: HitBufferMode,
    bins_processed: usize,
    input_range_next: (u32, u32),
    pending_hits: Vec<Hit>,
    spare_hits: Vec<Hit>,
    load_worker: Option<JoinHandle<Result<(Vec<Hit>, u64), String>>>,
    pending_grouped_hits: Vec<Vec<CompactHit>>,
    spare_grouped_hits: Vec<Vec<CompactHit>>,
    grouped_load_worker: Option<JoinHandle<Result<(Vec<Vec<CompactHit>>, u64), String>>>,
    pending: bool,
    writing_finished: bool,
    total_disk_size: u64,
    temp_files: Vec<PathBuf>,
    temp_writers: Vec<Option<BufWriter<File>>>,
    stream_text_buffers: Vec<Vec<u8>>,
    stream_buf_count: Vec<u32>,
    stream_last_key: Vec<Option<(u32, SeedOffset)>>,
    error: Option<String>,
    allocated: bool,
}

impl HitBuffer {
    /// Compatibility constructor for the historical in-memory Rust API.
    pub fn new(num_bins: usize) -> Self {
        Self::with_key_partition((1..=num_bins as u32).collect())
    }

    pub fn with_key_partition(key_partition: Vec<u32>) -> Self {
        Self::with_key_partition_and_contexts(key_partition, 1)
    }

    pub fn with_key_partition_and_contexts(key_partition: Vec<u32>, query_contexts: u32) -> Self {
        Self::with_limits(
            key_partition,
            "",
            true,
            query_contexts,
            u32::MAX,
            u64::MAX,
            HitBufferMode::Memory,
        )
        .expect("in-memory HitBuffer construction cannot fail")
    }

    /// Counterpart of the explicit C++ constructor. `thread_count` and the
    /// search pool only control overlap upstream; they do not alter the stored
    /// representation and are therefore absent here.
    #[allow(clippy::too_many_arguments)]
    pub fn with_limits(
        key_partition: Vec<u32>,
        tmpdir: impl AsRef<Path>,
        long_subject_offsets: bool,
        query_contexts: u32,
        max_query: u32,
        max_target: u64,
        mode: HitBufferMode,
    ) -> Result<Self, String> {
        if key_partition.is_empty() {
            return Err("HitBuffer requires at least one key partition".to_string());
        }
        if query_contexts == 0 {
            return Err("HitBuffer query_contexts must be positive".to_string());
        }
        if key_partition.windows(2).any(|w| w[0] >= w[1]) {
            return Err("HitBuffer key partitions must be strictly increasing".to_string());
        }

        let n = key_partition.len();
        let mut temp_files = Vec::new();
        let mut temp_writers = (0..n).map(|_| None).collect::<Vec<_>>();
        if mode == HitBufferMode::Disk {
            let dir = if tmpdir.as_ref().as_os_str().is_empty() {
                std::env::temp_dir()
            } else {
                tmpdir.as_ref().to_path_buf()
            };
            fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
            let id = TEMP_ID.fetch_add(1, Ordering::Relaxed);
            for bin in 0..n {
                let path = dir.join(format!(
                    "diamond-rs-hit-buffer-{}-{id}-{bin}.tmp",
                    std::process::id()
                ));
                let file = File::create(&path).map_err(|e| e.to_string())?;
                temp_files.push(path);
                temp_writers[bin] = Some(BufWriter::with_capacity(DISK_BUFFER_SIZE, file));
            }
        }

        Ok(Self {
            bins: vec![Vec::new(); n],
            count: vec![0; n],
            key_partition,
            long_subject_offsets,
            query_contexts,
            max_query,
            max_target,
            mode,
            bins_processed: 0,
            input_range_next: (0, 0),
            pending_hits: Vec::new(),
            spare_hits: Vec::new(),
            load_worker: None,
            pending_grouped_hits: Vec::new(),
            spare_grouped_hits: Vec::new(),
            grouped_load_worker: None,
            pending: false,
            writing_finished: false,
            total_disk_size: 0,
            temp_files,
            temp_writers,
            stream_text_buffers: vec![Vec::new(); n],
            stream_buf_count: vec![0; n],
            stream_last_key: vec![None; n],
            error: None,
            allocated: false,
        })
    }

    pub fn writer(&mut self) -> Writer<'_> {
        Writer::new(self)
    }

    /// Append directly to the persistent disk encoders. Upstream keeps one
    /// `TextBuffer` per bin for the lifetime of each search worker; retaining
    /// these buffers across partition batches avoids thousands of short
    /// packets and repeated allocation in forced-disk mode.
    pub fn append_disk_hit(&mut self, hit: Hit) -> Result<(), String> {
        if self.mode != HitBufferMode::Disk {
            return Err("HitBuffer::append_disk_hit(): disk mode required".to_string());
        }
        if hit.score == 0 {
            return Err("HitBuffer::append_disk_hit(): score must be positive".to_string());
        }
        if !self.long_subject_offsets && hit.subject > u32::MAX as u64 {
            return Err(
                "HitBuffer::append_disk_hit(): subject offset does not fit u32".to_string(),
            );
        }
        if self.long_subject_offsets && hit.subject >= (1u64 << 40) {
            return Err(
                "HitBuffer::append_disk_hit(): subject offset does not fit PackedLoc".to_string(),
            );
        }

        let bin = self.bin(hit.query / self.query_contexts)?;
        let key = (hit.query, hit.seed_offset);
        let subject_bytes = if self.long_subject_offsets { 5 } else { 4 };
        if self.stream_last_key[bin] != Some(key) {
            if self.stream_text_buffers[bin].len() + 10 + 2 + subject_bytes >= DISK_BUFFER_SIZE {
                self.flush_stream_bin(bin)?;
            }
            Self::start_stream_query(&mut self.stream_text_buffers[bin], key);
            self.stream_last_key[bin] = Some(key);
        }
        if self.stream_text_buffers[bin].len() + 2 + subject_bytes >= DISK_BUFFER_SIZE {
            self.flush_stream_bin(bin)?;
            Self::start_stream_query(&mut self.stream_text_buffers[bin], key);
        }

        self.stream_text_buffers[bin].extend_from_slice(&hit.score.to_ne_bytes());
        if self.long_subject_offsets {
            self.stream_text_buffers[bin].extend_from_slice(&hit.subject.to_le_bytes()[..5]);
        } else {
            self.stream_text_buffers[bin].extend_from_slice(&(hit.subject as u32).to_ne_bytes());
        }
        self.stream_buf_count[bin] += 1;
        self.count[bin] += 1;
        Ok(())
    }

    fn start_stream_query(buffer: &mut Vec<u8>, key: (u32, SeedOffset)) {
        buffer.extend_from_slice(&0u16.to_ne_bytes());
        buffer.extend_from_slice(&key.0.to_ne_bytes());
        buffer.extend_from_slice(&key.1.to_ne_bytes());
    }

    fn flush_stream_bin(&mut self, bin: usize) -> Result<(), String> {
        if self.stream_text_buffers[bin].is_empty() {
            return Ok(());
        }
        self.stream_text_buffers[bin].extend_from_slice(&0u16.to_ne_bytes());
        let count = self.stream_buf_count[bin];
        let payload = &self.stream_text_buffers[bin];
        let writer = self.temp_writers[bin]
            .as_mut()
            .ok_or_else(|| "HitBuffer::append_disk_hit(): writing already finished".to_string())?;
        writer
            .write_all(&payload.len().to_ne_bytes())
            .and_then(|_| writer.write_all(&count.to_ne_bytes()))
            .and_then(|_| writer.write_all(payload))
            .map_err(|error| error.to_string())?;
        self.stream_text_buffers[bin].clear();
        self.stream_buf_count[bin] = 0;
        Ok(())
    }

    pub fn mode(&self) -> HitBufferMode {
        self.mode
    }

    pub fn push(&mut self, bin: usize, hit: Hit) {
        self.bins[bin].push(hit);
        self.count[bin] += 1;
    }

    pub fn get_bin(&self, bin: usize) -> &[Hit] {
        &self.bins[bin]
    }

    pub fn get_bin_mut(&mut self, bin: usize) -> &mut Vec<Hit> {
        &mut self.bins[bin]
    }

    pub fn total_hits(&self) -> u64 {
        self.count.iter().map(|&n| n as u64).sum()
    }

    pub fn query_contexts(&self) -> u32 {
        self.query_contexts
    }

    pub fn num_bins(&self) -> usize {
        self.key_partition.len()
    }

    pub fn begin(&self, bin: usize) -> u32 {
        if bin == 0 {
            0
        } else {
            self.key_partition[bin - 1]
        }
    }

    pub fn end(&self, bin: usize) -> u32 {
        self.key_partition[bin]
    }

    pub fn bins(&self) -> i32 {
        self.key_partition.len() as i32
    }

    pub fn bin(&self, key: u32) -> Result<usize, String> {
        self.key_partition
            .iter()
            .position(|&end| key < end)
            .ok_or_else(|| "key_partition error".to_string())
    }

    pub fn bin_size(&self, bin: usize) -> i64 {
        self.count[bin] as i64
    }

    pub fn next_bin_size(&self) -> u64 {
        self.count.get(self.bins_processed).copied().unwrap_or(0) as u64
    }

    pub fn set_bins_processed(&mut self, bins_processed: usize) {
        self.bins_processed = bins_processed;
    }

    /// Synchronous counterpart of joining all C++ writer tasks.
    pub fn finish_writing(&mut self) {
        for bin in 0..self.num_bins() {
            if let Err(error) = self.flush_stream_bin(bin) {
                self.error = Some(error);
            }
        }
        for writer in &mut self.temp_writers {
            if let Some(mut writer) = writer.take() {
                if let Err(error) = writer.flush() {
                    self.error = Some(error.to_string());
                }
            }
        }
        self.writing_finished = true;
    }

    pub fn try_finish_writing(&mut self) -> Result<(), String> {
        self.finish_writing();
        self.take_error()
    }

    /// Prepare the next bin. C++ currently loads exactly one bin despite the
    /// `max_size` argument (`end - bins_processed_ == 0` is false immediately);
    /// this mirrors that observable behavior.
    pub fn load(&mut self, _max_size: usize) -> bool {
        if self.pending {
            self.error = Some("HitBuffer::load(): previous load still in progress".to_string());
            return false;
        }
        if self.bins_processed == self.num_bins() {
            return false;
        }
        let bin = self.bins_processed;
        self.input_range_next = (self.begin(bin), self.end(bin));
        match self.mode {
            HitBufferMode::Memory => {}
            HitBufferMode::SwipeAll => {}
            HitBufferMode::Disk => {
                let path = self.temp_files[bin].clone();
                let expected_count = self.count[bin];
                let long_subject_offsets = self.long_subject_offsets;
                let max_query = self.max_query;
                let max_target = self.max_target;
                let mut reuse = if self.spare_hits.capacity() >= self.pending_hits.capacity() {
                    std::mem::take(&mut self.spare_hits)
                } else {
                    std::mem::take(&mut self.pending_hits)
                };
                reuse.clear();
                match thread::Builder::new()
                    .name("diamond-rs-hit-loader".to_string())
                    .spawn(move || {
                        Self::load_bin_from(
                            path,
                            expected_count,
                            long_subject_offsets,
                            max_query,
                            max_target,
                            reuse,
                        )
                    }) {
                    Ok(worker) => self.load_worker = Some(worker),
                    Err(error) => self.error = Some(error.to_string()),
                }
                self.bins_processed += 1;
            }
        }
        self.pending = true;
        true
    }

    /// Native blastp loader: decode directly into query-local compact hit
    /// vectors instead of materializing a flat `Hit` array which the caller
    /// would immediately copy and regroup.
    pub(crate) fn load_grouped(&mut self) -> bool {
        if self.pending || self.load_worker.is_some() || self.grouped_load_worker.is_some() {
            self.error = Some("HitBuffer::load_grouped(): previous load still in progress".into());
            return false;
        }
        if self.mode != HitBufferMode::Disk {
            self.error = Some("HitBuffer::load_grouped(): disk mode required".into());
            return false;
        }
        if self.query_contexts != 1 {
            self.error =
                Some("HitBuffer::load_grouped(): only one query context is supported".into());
            return false;
        }
        if self.bins_processed == self.num_bins() {
            return false;
        }
        let bin = self.bins_processed;
        let begin = self.begin(bin);
        let end = self.end(bin);
        self.input_range_next = (begin, end);
        let path = self.temp_files[bin].clone();
        let expected_count = self.count[bin];
        let long_subject_offsets = self.long_subject_offsets;
        let max_query = self.max_query;
        let max_target = self.max_target;
        let reuse = std::mem::take(&mut self.spare_grouped_hits);
        match thread::Builder::new()
            .name("diamond-rs-grouped-hit-loader".to_string())
            .spawn(move || {
                Self::load_grouped_bin_from(
                    path,
                    expected_count,
                    long_subject_offsets,
                    max_query,
                    max_target,
                    begin,
                    end,
                    reuse,
                )
            }) {
            Ok(worker) => self.grouped_load_worker = Some(worker),
            Err(error) => self.error = Some(error.to_string()),
        }
        self.bins_processed += 1;
        self.pending = true;
        true
    }

    /// Fallible form of C++ `retrieve()`, including propagation of an
    /// exception raised by the load worker.
    pub fn try_retrieve(&mut self) -> Result<Option<(&[Hit], u32, u32)>, String> {
        self.finish_load_worker()?;
        if let Some(error) = self.error.take() {
            self.pending = false;
            return Err(error);
        }
        if !self.pending {
            if self.bins_processed >= self.num_bins() {
                return Ok(None);
            }
            return Err("HitBuffer retrieve w/o load".to_string());
        }
        self.pending = false;
        let (begin, end) = self.input_range_next;
        match self.mode {
            HitBufferMode::Memory => {
                let bin = self.bins_processed;
                self.bins_processed += 1;
                Ok(Some((&self.bins[bin], begin, end)))
            }
            HitBufferMode::SwipeAll => {
                self.bins_processed += 1;
                Ok(Some((&[], begin, end)))
            }
            HitBufferMode::Disk => Ok(Some((&self.pending_hits, begin, end))),
        }
    }

    /// Owned disk retrieval permits the caller to start loading the next bin
    /// before processing this one, matching upstream's double-buffered path.
    pub fn try_retrieve_owned(&mut self) -> Result<Option<(Vec<Hit>, u32, u32)>, String> {
        if self.mode != HitBufferMode::Disk {
            return Err("HitBuffer::try_retrieve_owned(): disk mode required".to_string());
        }
        self.finish_load_worker()?;
        if let Some(error) = self.error.take() {
            self.pending = false;
            return Err(error);
        }
        if !self.pending {
            if self.bins_processed >= self.num_bins() {
                return Ok(None);
            }
            return Err("HitBuffer retrieve w/o load".to_string());
        }
        self.pending = false;
        let (begin, end) = self.input_range_next;
        Ok(Some((std::mem::take(&mut self.pending_hits), begin, end)))
    }

    /// Return a consumed owned result buffer for reuse by a later load.
    pub fn recycle_hits(&mut self, mut hits: Vec<Hit>) {
        hits.clear();
        if hits.capacity() > self.spare_hits.capacity() {
            self.spare_hits = hits;
        }
    }

    pub(crate) fn try_retrieve_grouped_owned(
        &mut self,
    ) -> Result<Option<(Vec<Vec<CompactHit>>, u32, u32)>, String> {
        self.finish_grouped_load_worker()?;
        if let Some(error) = self.error.take() {
            self.pending = false;
            return Err(error);
        }
        if !self.pending {
            if self.bins_processed >= self.num_bins() {
                return Ok(None);
            }
            return Err("HitBuffer grouped retrieve w/o load".to_string());
        }
        self.pending = false;
        let (begin, end) = self.input_range_next;
        Ok(Some((
            std::mem::take(&mut self.pending_grouped_hits),
            begin,
            end,
        )))
    }

    pub(crate) fn recycle_grouped_hits(&mut self, mut hits: Vec<Vec<CompactHit>>) {
        for query_hits in &mut hits {
            query_hits.clear();
        }
        self.spare_grouped_hits = hits;
    }

    /// Compatibility wrapper for callers of the earlier in-memory API. New
    /// disk-backed callers should use `try_retrieve` so corruption is visible.
    pub fn retrieve(&mut self) -> Option<(&[Hit], u32, u32)> {
        self.try_retrieve().ok().flatten()
    }

    pub fn take_error(&mut self) -> Result<(), String> {
        match self.error.take() {
            Some(e) => Err(e),
            None => Ok(()),
        }
    }

    pub fn total_disk_size(&self) -> i64 {
        self.total_disk_size as i64
    }

    /// C++ reserves two maximum-bin buffers, using mmap where available. Vec
    /// owns the equivalent storage safely and grows only to the loaded count.
    pub fn alloc_buffer(&mut self) {
        self.allocated = true;
        if self.mode == HitBufferMode::Disk {
            let max_size = self.count.iter().copied().max().unwrap_or(0);
            self.pending_hits.reserve(max_size);
        }
    }

    pub fn free_buffer(&mut self) {
        let _ = self.finish_load_worker();
        let _ = self.finish_grouped_load_worker();
        self.pending_hits = Vec::new();
        self.spare_hits = Vec::new();
        self.pending_grouped_hits = Vec::new();
        self.spare_grouped_hits = Vec::new();
        self.allocated = false;
    }

    pub fn clear(&mut self) {
        let load_error = self.finish_load_worker().err();
        let grouped_load_error = self.finish_grouped_load_worker().err();
        for bin in &mut self.bins {
            bin.clear();
        }
        self.count.fill(0);
        self.bins_processed = 0;
        self.input_range_next = (0, 0);
        self.pending_hits.clear();
        self.spare_hits.clear();
        self.pending_grouped_hits.clear();
        self.spare_grouped_hits.clear();
        self.pending = false;
        self.writing_finished = false;
        self.total_disk_size = 0;
        for buffer in &mut self.stream_text_buffers {
            buffer.clear();
        }
        self.stream_buf_count.fill(0);
        self.stream_last_key.fill(None);
        self.error = load_error.or(grouped_load_error);
        if self.mode == HitBufferMode::Disk {
            for (bin, path) in self.temp_files.iter().enumerate() {
                match File::create(path) {
                    Ok(file) => {
                        self.temp_writers[bin] =
                            Some(BufWriter::with_capacity(DISK_BUFFER_SIZE, file));
                    }
                    Err(e) => {
                        self.error = Some(e.to_string());
                        break;
                    }
                }
            }
        }
    }

    pub fn sort_by_subject(&mut self) {
        for bin in &mut self.bins {
            bin.sort_by_key(|h| h.subject);
        }
    }

    /// Synchronous equivalent of `HitBuffer::write_worker`: consumes one
    /// queued packet for the selected bin.
    fn write_worker(&mut self, bin: usize, payload: &[u8], count: u32) -> Result<(), String> {
        let file = self.temp_writers[bin]
            .as_mut()
            .ok_or_else(|| "HitBuffer::write_worker(): writing already finished".to_string())?;
        file.write_all(&payload.len().to_ne_bytes())
            .and_then(|_| file.write_all(&count.to_ne_bytes()))
            .and_then(|_| file.write_all(payload))
            .map_err(|e| e.to_string())
    }

    fn finish_load_worker(&mut self) -> Result<(), String> {
        let Some(worker) = self.load_worker.take() else {
            return Ok(());
        };
        match worker.join() {
            Ok(Ok((hits, disk_size))) => {
                self.pending_hits = hits;
                self.total_disk_size += disk_size;
                Ok(())
            }
            Ok(Err(error)) => Err(error),
            Err(_) => Err("HitBuffer::load_bin(): background loader panicked".to_string()),
        }
    }

    fn finish_grouped_load_worker(&mut self) -> Result<(), String> {
        let Some(worker) = self.grouped_load_worker.take() else {
            return Ok(());
        };
        match worker.join() {
            Ok(Ok((hits, disk_size))) => {
                self.pending_grouped_hits = hits;
                self.total_disk_size += disk_size;
                Ok(())
            }
            Ok(Err(error)) => Err(error),
            Err(_) => Err("HitBuffer::load_grouped_bin(): background loader panicked".to_string()),
        }
    }

    fn load_bin_from(
        path: PathBuf,
        expected_count: usize,
        long_subject_offsets: bool,
        max_query: u32,
        max_target: u64,
        mut hits: Vec<Hit>,
    ) -> Result<(Vec<Hit>, u64), String> {
        let disk_size = fs::metadata(&path).map_err(|e| e.to_string())?.len();
        hits.clear();
        if expected_count == 0 {
            return Ok((hits, disk_size));
        }
        hits.reserve(expected_count.saturating_sub(hits.capacity()));
        let mut file = File::open(path).map_err(|e| e.to_string())?;
        file.seek(SeekFrom::Start(0)).map_err(|e| e.to_string())?;
        loop {
            let mut size_bytes = [0u8; std::mem::size_of::<usize>()];
            let n = file.read(&mut size_bytes).map_err(|e| e.to_string())?;
            if n == 0 {
                break;
            }
            if n != size_bytes.len() {
                return Err("HitBuffer::load_bin(): truncated packet size".to_string());
            }
            let len = usize::from_ne_bytes(size_bytes);
            let expected = read_u32_from(&mut file)? as usize;
            let mut payload = vec![0; len];
            file.read_exact(&mut payload).map_err(|_| {
                "HitBuffer::load_bin(): truncated packet / possibly corrupted temporary file"
                    .to_string()
            })?;
            Self::decode_packet(
                &payload,
                expected,
                long_subject_offsets,
                max_query,
                max_target,
                &mut hits,
            )?;
        }
        if hits.len() != expected_count {
            return Err("Mismatching hit count / possibly corrupted temporary file".to_string());
        }
        Ok((hits, disk_size))
    }

    fn decode_packet(
        payload: &[u8],
        expected: usize,
        long_subject_offsets: bool,
        max_query: u32,
        max_target: u64,
        out: &mut Vec<Hit>,
    ) -> Result<(), String> {
        let mut pos = 0usize;
        let nullscore = take_u16(payload, &mut pos)?;
        if nullscore != 0 {
            return Err("HitBuffer::load_bin(): invalid packet header".to_string());
        }
        let before = out.len();
        while pos < payload.len() {
            let query = take_u32(payload, &mut pos)?;
            if query >= max_query {
                return Err(
                    "HitBuffer::load_bin(): invalid query id / possibly corrupted temporary file"
                        .to_string(),
                );
            }
            let seed_offset = take_u32(payload, &mut pos)?;
            loop {
                let score = take_u16(payload, &mut pos)?;
                if score == 0 {
                    break;
                }
                let subject = if long_subject_offsets {
                    take_u40(payload, &mut pos)?
                } else {
                    take_u32(payload, &mut pos)? as u64
                };
                if subject >= max_target {
                    return Err("HitBuffer::load_bin(): invalid subject location / possibly corrupted temporary file".to_string());
                }
                if out.len() - before >= expected {
                    return Err("HitBuffer::load_bin(): buffer overflow / possibly corrupted temporary file".to_string());
                }
                out.push(Hit::with_score(query, subject, seed_offset, score));
            }
        }
        if out.len() - before != expected {
            return Err("Mismatching hit count / possibly corrupted temporary file".to_string());
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn load_grouped_bin_from(
        path: PathBuf,
        expected_count: usize,
        long_subject_offsets: bool,
        max_query: u32,
        max_target: u64,
        begin: u32,
        end: u32,
        mut hits: Vec<Vec<CompactHit>>,
    ) -> Result<(Vec<Vec<CompactHit>>, u64), String> {
        let disk_size = fs::metadata(&path).map_err(|e| e.to_string())?.len();
        let query_count = end.saturating_sub(begin) as usize;
        hits.resize_with(query_count, Vec::new);
        hits.truncate(query_count);
        for query_hits in &mut hits {
            query_hits.clear();
        }
        if expected_count == 0 {
            return Ok((hits, disk_size));
        }

        let mut file = File::open(path).map_err(|e| e.to_string())?;
        let mut decoded = 0usize;
        loop {
            let mut size_bytes = [0u8; std::mem::size_of::<usize>()];
            let n = file.read(&mut size_bytes).map_err(|e| e.to_string())?;
            if n == 0 {
                break;
            }
            if n != size_bytes.len() {
                return Err("HitBuffer::load_grouped_bin(): truncated packet size".to_string());
            }
            let len = usize::from_ne_bytes(size_bytes);
            let packet_count = read_u32_from(&mut file)? as usize;
            let mut payload = vec![0; len];
            file.read_exact(&mut payload).map_err(|_| {
                "HitBuffer::load_grouped_bin(): truncated packet / possibly corrupted temporary file"
                    .to_string()
            })?;
            let mut pos = 0usize;
            if take_u16(&payload, &mut pos)? != 0 {
                return Err("HitBuffer::load_grouped_bin(): invalid packet header".to_string());
            }
            let packet_begin = decoded;
            while pos < payload.len() {
                let query = take_u32(&payload, &mut pos)?;
                if query >= max_query || query < begin || query >= end {
                    return Err("HitBuffer::load_grouped_bin(): invalid query id / possibly corrupted temporary file".to_string());
                }
                let seed_offset = take_u32(&payload, &mut pos)?;
                let query_hits = &mut hits[(query - begin) as usize];
                loop {
                    let score = take_u16(&payload, &mut pos)?;
                    if score == 0 {
                        break;
                    }
                    let subject = if long_subject_offsets {
                        take_u40(&payload, &mut pos)?
                    } else {
                        take_u32(&payload, &mut pos)? as u64
                    };
                    if subject >= max_target {
                        return Err("HitBuffer::load_grouped_bin(): invalid subject location / possibly corrupted temporary file".to_string());
                    }
                    if decoded - packet_begin >= packet_count || decoded >= expected_count {
                        return Err("HitBuffer::load_grouped_bin(): buffer overflow / possibly corrupted temporary file".to_string());
                    }
                    query_hits.push(CompactHit {
                        subject,
                        seed_offset,
                        score,
                    });
                    decoded += 1;
                }
            }
            if decoded - packet_begin != packet_count {
                return Err("Mismatching hit count / possibly corrupted temporary file".to_string());
            }
        }
        if decoded != expected_count {
            return Err("Mismatching hit count / possibly corrupted temporary file".to_string());
        }
        Ok((hits, disk_size))
    }
}

impl Drop for HitBuffer {
    fn drop(&mut self) {
        let _ = self.finish_load_worker();
        let _ = self.finish_grouped_load_worker();
        for path in &self.temp_files {
            let _ = fs::remove_file(path);
        }
    }
}

/// C++ `HitBuffer::Writer`; buffers each bin independently and flushes all
/// remaining packets when dropped.
pub struct Writer<'a> {
    buffer_size: usize,
    last_bin: usize,
    seed_offset: SeedOffset,
    query: BlockId,
    buffers: Vec<Vec<Hit>>,
    text_buffers: Vec<Vec<u8>>,
    count: Vec<usize>,
    buf_count: Vec<u32>,
    parent: &'a mut HitBuffer,
}

impl<'a> Writer<'a> {
    pub fn new(parent: &'a mut HitBuffer) -> Self {
        let n = parent.num_bins();
        let memory = parent.mode == HitBufferMode::Memory;
        Self {
            buffer_size: if memory {
                MEMORY_BUFFER_SIZE
            } else {
                DISK_BUFFER_SIZE
            },
            last_bin: 0,
            seed_offset: 0,
            query: 0,
            buffers: vec![Vec::new(); n],
            text_buffers: vec![Vec::new(); n],
            count: vec![0; n],
            buf_count: vec![0; n],
            parent,
        }
    }

    pub fn new_query(&mut self, query: BlockId, seed_offset: SeedOffset) -> Result<(), String> {
        self.last_bin = self.parent.bin(query / self.parent.query_contexts)?;
        self.seed_offset = seed_offset;
        self.query = query;
        if self.parent.mode == HitBufferMode::Disk {
            self.start_query(self.last_bin, query, seed_offset);
        }
        Ok(())
    }

    pub fn write(
        &mut self,
        query: BlockId,
        subject: u64,
        score: u16,
        _target_block_id: u32,
    ) -> Result<(), String> {
        if score == 0 {
            return Err("HitBuffer::Writer::write(): score must be positive".to_string());
        }
        if self.last_bin >= self.parent.num_bins() {
            return Err("HitBuffer::Writer::write(): invalid bin".to_string());
        }
        if !self.parent.long_subject_offsets && subject > u32::MAX as u64 {
            return Err("HitBuffer::Writer::write(): subject offset does not fit u32".to_string());
        }
        if self.parent.long_subject_offsets && subject >= (1u64 << 40) {
            return Err(
                "HitBuffer::Writer::write(): subject offset does not fit PackedLoc".to_string(),
            );
        }

        let bin = self.last_bin;
        match self.parent.mode {
            HitBufferMode::Memory => {
                if self.buffers[bin].len() >= self.buffer_size {
                    self.flush(bin, false)?;
                }
                self.buffers[bin].push(Hit::with_score(query, subject, self.seed_offset, score));
            }
            HitBufferMode::Disk => {
                let subject_bytes = if self.parent.long_subject_offsets {
                    5
                } else {
                    4
                };
                if self.text_buffers[bin].len() + 2 + subject_bytes >= self.buffer_size {
                    self.flush(bin, false)?;
                    self.start_query(bin, self.query, self.seed_offset);
                }
                self.text_buffers[bin].extend_from_slice(&score.to_ne_bytes());
                if self.parent.long_subject_offsets {
                    self.text_buffers[bin].extend_from_slice(&subject.to_le_bytes()[..5]);
                } else {
                    self.text_buffers[bin].extend_from_slice(&(subject as u32).to_ne_bytes());
                }
            }
            HitBufferMode::SwipeAll => {}
        }
        self.count[bin] += 1;
        self.buf_count[bin] += 1;
        Ok(())
    }

    pub fn count(&self, bin: usize) -> usize {
        self.count[bin]
    }

    pub fn flush(&mut self, bin: usize, done: bool) -> Result<(), String> {
        match self.parent.mode {
            HitBufferMode::Memory => {
                if !self.buffers[bin].is_empty() {
                    self.parent.bins[bin].append(&mut self.buffers[bin]);
                }
                if !done {
                    self.buffers[bin].reserve(self.buffer_size);
                }
            }
            HitBufferMode::Disk => {
                if !self.text_buffers[bin].is_empty() {
                    self.text_buffers[bin].extend_from_slice(&0u16.to_ne_bytes());
                    let payload = std::mem::take(&mut self.text_buffers[bin]);
                    let n = self.buf_count[bin];
                    self.parent.write_worker(bin, &payload, n)?;
                    self.buf_count[bin] = 0;
                }
            }
            HitBufferMode::SwipeAll => {}
        }
        Ok(())
    }

    fn start_query(&mut self, bin: usize, query: u32, seed_offset: u32) {
        self.text_buffers[bin].extend_from_slice(&0u16.to_ne_bytes());
        self.text_buffers[bin].extend_from_slice(&query.to_ne_bytes());
        self.text_buffers[bin].extend_from_slice(&seed_offset.to_ne_bytes());
    }
}

impl Drop for Writer<'_> {
    fn drop(&mut self) {
        for bin in 0..self.parent.num_bins() {
            if let Err(e) = self.flush(bin, true) {
                self.parent.error = Some(e);
            }
            self.parent.count[bin] += self.count[bin];
        }
    }
}

fn read_u32_from(file: &mut File) -> Result<u32, String> {
    let mut bytes = [0; 4];
    file.read_exact(&mut bytes)
        .map_err(|_| "HitBuffer::load_bin(): truncated packet count".to_string())?;
    Ok(u32::from_ne_bytes(bytes))
}

fn take_u16(input: &[u8], pos: &mut usize) -> Result<u16, String> {
    let bytes = take(input, pos, 2)?;
    Ok(u16::from_ne_bytes([bytes[0], bytes[1]]))
}

fn take_u32(input: &[u8], pos: &mut usize) -> Result<u32, String> {
    let bytes = take(input, pos, 4)?;
    Ok(u32::from_ne_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
}

fn take_u40(input: &[u8], pos: &mut usize) -> Result<u64, String> {
    let bytes = take(input, pos, 5)?;
    Ok(u64::from_le_bytes([
        bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], 0, 0, 0,
    ]))
}

fn take<'a>(input: &'a [u8], pos: &mut usize, n: usize) -> Result<&'a [u8], String> {
    let end = pos
        .checked_add(n)
        .filter(|&end| end <= input.len())
        .ok_or_else(|| {
            "HitBuffer::load_bin(): truncated payload / possibly corrupted temporary file"
                .to_string()
        })?;
    let out = &input[*pos..end];
    *pos = end;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_dir() -> PathBuf {
        std::env::temp_dir()
    }

    #[test]
    fn disk_packets_roundtrip_short_and_long_subject_offsets() {
        for long in [false, true] {
            let max_target = if long { 1u64 << 39 } else { u32::MAX as u64 };
            let subject = if long { (1u64 << 35) + 17 } else { 123_456 };
            let mut buffer = HitBuffer::with_limits(
                vec![5, 10],
                temp_dir(),
                long,
                2,
                20,
                max_target,
                HitBufferMode::Disk,
            )
            .unwrap();
            {
                let mut writer = buffer.writer();
                writer.new_query(2, 7).unwrap();
                writer.write(2, subject, 11, 0).unwrap();
                writer.write(2, subject + 1, 12, 0).unwrap();
                writer.new_query(12, 9).unwrap();
                writer.write(12, 42, 13, 0).unwrap();
            }
            buffer.try_finish_writing().unwrap();
            buffer.alloc_buffer();
            assert!(buffer.load(1));
            let (hits, begin, end) = buffer.retrieve().unwrap();
            assert_eq!((begin, end), (0, 5));
            assert_eq!(hits.len(), 2);
            assert_eq!(hits[0], Hit::with_score(2, subject, 7, 11));
            assert_eq!(hits[1], Hit::with_score(2, subject + 1, 7, 12));
            assert!(buffer.load(usize::MAX));
            let (hits, begin, end) = buffer.retrieve().unwrap();
            assert_eq!((begin, end), (5, 10));
            assert_eq!(hits, &[Hit::with_score(12, 42, 9, 13)]);
            assert!(buffer.total_disk_size() > 0);
            assert!(!buffer.load(1));
            buffer.free_buffer();
        }
    }

    #[test]
    fn persistent_disk_encoder_roundtrips_flush_boundary_and_overlapped_bins() {
        let mut buffer = HitBuffer::with_limits(
            vec![10, 20],
            temp_dir(),
            false,
            1,
            20,
            50_000,
            HitBufferMode::Disk,
        )
        .unwrap();
        for subject in 0..11_000u64 {
            buffer
                .append_disk_hit(Hit::with_score(2, subject, 7, 11))
                .unwrap();
        }
        buffer
            .append_disk_hit(Hit::with_score(12, 42, 9, 13))
            .unwrap();
        buffer
            .append_disk_hit(Hit::with_score(12, 43, 9, 14))
            .unwrap();
        buffer.try_finish_writing().unwrap();
        buffer.alloc_buffer();

        assert!(buffer.load(usize::MAX));
        let (first, begin, end) = buffer.try_retrieve_owned().unwrap().unwrap();
        assert_eq!((begin, end), (0, 10));
        assert!(buffer.load(usize::MAX));
        assert_eq!(first.len(), 11_000);
        assert_eq!(first[0], Hit::with_score(2, 0, 7, 11));
        assert_eq!(first[10_999], Hit::with_score(2, 10_999, 7, 11));
        buffer.recycle_hits(first);

        let (second, begin, end) = buffer.try_retrieve_owned().unwrap().unwrap();
        assert_eq!((begin, end), (10, 20));
        assert_eq!(
            second,
            vec![
                Hit::with_score(12, 42, 9, 13),
                Hit::with_score(12, 43, 9, 14)
            ]
        );
        assert!(!buffer.load(usize::MAX));
        assert!(buffer.total_disk_size() > 0);
    }

    #[test]
    fn grouped_disk_decoder_preserves_query_order_and_reuses_bins() {
        let mut buffer = HitBuffer::with_limits(
            vec![4, 8],
            temp_dir(),
            false,
            1,
            8,
            10_000,
            HitBufferMode::Disk,
        )
        .unwrap();
        for hit in [
            Hit::with_score(2, 20, 5, 11),
            Hit::with_score(1, 10, 3, 7),
            Hit::with_score(2, 21, 6, 12),
            Hit::with_score(6, 60, 9, 13),
        ] {
            buffer.append_disk_hit(hit).unwrap();
        }
        buffer.try_finish_writing().unwrap();

        assert!(buffer.load_grouped());
        let (first, begin, end) = buffer.try_retrieve_grouped_owned().unwrap().unwrap();
        assert_eq!((begin, end), (0, 4));
        assert_eq!(
            first[1],
            vec![CompactHit {
                subject: 10,
                seed_offset: 3,
                score: 7
            }]
        );
        assert_eq!(
            first[2].iter().map(|hit| hit.subject).collect::<Vec<_>>(),
            vec![20, 21]
        );
        buffer.recycle_grouped_hits(first);

        assert!(buffer.load_grouped());
        let (second, begin, end) = buffer.try_retrieve_grouped_owned().unwrap().unwrap();
        assert_eq!((begin, end), (4, 8));
        assert_eq!(
            second[2],
            vec![CompactHit {
                subject: 60,
                seed_offset: 9,
                score: 13
            }]
        );
        assert!(!buffer.load_grouped());
    }

    #[test]
    fn validates_partitions_bounds_and_load_state() {
        assert!(
            HitBuffer::with_limits(vec![], temp_dir(), false, 1, 1, 1, HitBufferMode::Disk)
                .is_err()
        );
        assert!(HitBuffer::with_limits(
            vec![2, 2],
            temp_dir(),
            false,
            1,
            1,
            1,
            HitBufferMode::Disk
        )
        .is_err());

        let mut buffer =
            HitBuffer::with_limits(vec![10], temp_dir(), false, 1, 2, 100, HitBufferMode::Disk)
                .unwrap();
        {
            let mut writer = buffer.writer();
            writer.new_query(2, 0).unwrap();
            writer.write(2, 1, 1, 0).unwrap();
            assert!(writer.write(2, 1, 0, 0).is_err());
        }
        assert!(buffer.load(1));
        assert!(buffer.try_retrieve().is_err());
    }

    #[test]
    fn swipe_all_yields_ranges_without_hits() {
        let mut buffer = HitBuffer::with_limits(
            vec![3, 7],
            temp_dir(),
            false,
            1,
            8,
            8,
            HitBufferMode::SwipeAll,
        )
        .unwrap();
        assert!(buffer.load(1));
        assert_eq!(buffer.retrieve().unwrap(), (&[][..], 0, 3));
        assert!(buffer.load(1));
        assert_eq!(buffer.retrieve().unwrap(), (&[][..], 3, 7));
        assert!(!buffer.load(1));
    }
}
