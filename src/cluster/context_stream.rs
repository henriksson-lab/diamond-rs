//! Source basic/match.h serialized HspContext streams used by cluster/output.cpp.
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Seek, SeekFrom, Write};
use std::path::PathBuf;

use crate::align::hsp::{ByteHspContext as HspContext, Hsp};
use crate::util::interval::Interval;

#[derive(Clone, Debug)]
pub struct ClusterOutputIo {
    pub tmpdir: PathBuf,
    pub file_buffer_size: usize,
    pub no_unlink: bool,
}
impl Default for ClusterOutputIo {
    fn default() -> Self {
        Self {
            tmpdir: PathBuf::new(),
            file_buffer_size: 67_108_864,
            no_unlink: false,
        }
    }
}

pub(super) struct ContextFile {
    file: Option<File>,
    path: PathBuf,
    unlinked: bool,
    no_unlink: bool,
    capacity: usize,
}
impl ContextFile {
    pub fn create(config: &ClusterOutputIo) -> io::Result<(Self, BufWriter<File>)> {
        let mut reserved =
            crate::util::io::TempFileData::init_in(false, &config.tmpdir, config.no_unlink)
                .map_err(io::Error::other)?;
        let path = PathBuf::from(&reserved.name);
        // mkstemp's actual descriptor is the stream owner. Closing it and
        // reopening its pathname introduces a race and changes failure order.
        #[cfg(unix)]
        let file = reserved.take_file().map_err(io::Error::other)?;
        // Source's non-POSIX reservation is a name-only temporary; its wb+
        // stream constructor creates the file at this boundary.
        #[cfg(not(unix))]
        let file = File::options()
            .read(true)
            .write(true)
            .create(true)
            .truncate(true)
            .open(&path)
            .map_err(|error| crate::util::io::FileOpenError::into_io(&path, error))?;
        drop(reserved);
        let mut output = Self {
            file: Some(file),
            path,
            unlinked: false,
            no_unlink: config.no_unlink,
            capacity: config.file_buffer_size,
        };
        #[cfg(unix)]
        if !config.no_unlink {
            std::fs::remove_file(&output.path)?;
            output.unlinked = true;
        }
        let writer = BufWriter::with_capacity(
            config.file_buffer_size,
            output.file.as_ref().unwrap().try_clone()?,
        );
        Ok((output, writer))
    }
    fn reader(&self) -> io::Result<BufReader<File>> {
        let mut file = self.file.as_ref().unwrap().try_clone()?;
        file.seek(SeekFrom::Start(0))?;
        Ok(BufReader::with_capacity(self.capacity, file))
    }
    fn close_and_delete(mut self) {
        // Source explicitly deletes successful streams even with --no-unlink.
        drop(self.file.take());
        if !self.unlinked {
            if std::fs::remove_file(&self.path).is_err() {
                let path = crate::util::io::path_diagnostic_bytes(&self.path)
                    .unwrap_or_else(|error| crate::util::io::abort_source_diagnostic_error(&error));
                let mut warning = b"Warning: Failed to delete temporary file ".to_vec();
                warning.extend_from_slice(&path);
                warning.push(b'\n');
                let _ = crate::util::io::console::stderr().write_all(&warning);
            } else {
                self.unlinked = true;
            }
        }
    }
}
impl Drop for ContextFile {
    fn drop(&mut self) {
        drop(self.file.take());
        if !self.unlinked && !self.no_unlink {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

pub(super) fn capture_io<T>(
    result: io::Result<T>,
    errors: Option<&crate::util::io::IoErrorOwner>,
) -> Result<T, String> {
    match errors {
        Some(errors) => errors.capture(result),
        None => result.map_err(|error| error.to_string()),
    }
}

pub(super) fn write_context<T: AsRef<[u8]>>(
    out: &mut impl Write,
    h: &crate::align::hsp::HspContext<T>,
) -> io::Result<()> {
    out.write_all(&h.query_id.to_ne_bytes())?;
    out.write_all(&h.query_oid.to_ne_bytes())?;
    out.write_all(&h.subject_oid.to_ne_bytes())?;
    for title in [&h.query_title, &h.target_title] {
        out.write_all(crate::util::sequence::c_string_bytes(title.as_ref()))?;
        out.write_all(&[0])?;
    }
    for value in [
        h.query_len,
        h.subject_len,
        h.hsp.identities,
        h.hsp.mismatches,
        h.hsp.positives,
        h.hsp.gaps,
        h.hsp.length,
        h.hsp.gap_openings,
        h.hsp.query_range.begin,
        h.hsp.query_range.end,
        h.hsp.subject_range.begin,
        h.hsp.subject_range.end,
    ] {
        out.write_all(&value.to_ne_bytes())?;
    }
    out.write_all(&h.hsp.bit_score.to_ne_bytes())?;
    out.write_all(&h.hsp.evalue.to_ne_bytes())?;
    out.write_all(&h.hsp.score.to_ne_bytes())?;
    out.write_all(&h.hsp.approx_id.to_ne_bytes())?;
    #[cfg(feature = "dna")]
    for value in [u32::from(h.hsp.mapping_quality), h.hsp.n_anchors as u32] {
        out.write_all(&value.to_ne_bytes())?;
    }
    Ok(())
}

fn read_context(input: &mut impl BufRead) -> io::Result<Option<HspContext>> {
    let mut query_id = [0u8; 4];
    if input.read(&mut query_id[..1])? == 0 {
        return Ok(None);
    }
    input.read_exact(&mut query_id[1..])?;
    macro_rules! read {
        ($ty:ty,$n:expr) => {{
            let mut bytes = [0u8; $n];
            input.read_exact(&mut bytes)?;
            <$ty>::from_ne_bytes(bytes)
        }};
    }
    let mut out = HspContext::default();
    out.query_id = u32::from_ne_bytes(query_id);
    out.query_oid = read!(u64, 8);
    out.subject_oid = read!(u64, 8);
    for title in [&mut out.query_title, &mut out.target_title] {
        let mut bytes = Vec::new();
        if input.read_until(0, &mut bytes)? == 0 || bytes.pop() != Some(0) {
            return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
        }
        *title = bytes;
    }
    out.query_len = read!(i32, 4);
    out.subject_len = read!(i32, 4);
    let mut hsp = Hsp::default();
    hsp.identities = read!(i32, 4);
    hsp.mismatches = read!(i32, 4);
    hsp.positives = read!(i32, 4);
    hsp.gaps = read!(i32, 4);
    hsp.length = read!(i32, 4);
    hsp.gap_openings = read!(i32, 4);
    hsp.query_range = Interval::new(read!(i32, 4), read!(i32, 4));
    hsp.subject_range = Interval::new(read!(i32, 4), read!(i32, 4));
    hsp.bit_score = read!(f64, 8);
    hsp.evalue = read!(f64, 8);
    hsp.score = read!(i32, 4);
    hsp.approx_id = read!(f64, 8);
    #[cfg(feature = "dna")]
    {
        hsp.mapping_quality = read!(u32, 4) as u8;
        hsp.n_anchors = read!(i32, 4);
    }
    hsp.query_source_range = hsp.query_range;
    hsp.subject_source_range = hsp.subject_range;
    out.hsp = hsp;
    // Source's default TranslatedSequence contains six empty Sequence objects.
    out.query = vec![Vec::new(); 6];
    Ok(Some(out))
}

pub(super) fn merge_files<F>(
    files: Vec<ContextFile>,
    mut callback: F,
    errors: Option<&crate::util::io::IoErrorOwner>,
) -> Result<(), String>
where
    F: FnMut(&HspContext) -> Result<(), String>,
{
    let mut readers = capture_io(
        files
            .iter()
            .map(ContextFile::reader)
            .collect::<io::Result<Vec<_>>>(),
        errors,
    )?;
    // Same push/adjust operations as Source merge_sorted_files' priority_queue:
    // equal query OIDs rise, without an invented stream-index tie breaker.
    let mut heap: Vec<(HspContext, usize)> = Vec::new();
    fn push(heap: &mut Vec<(HspContext, usize)>, value: (HspContext, usize)) {
        heap.push(value);
        let mut hole = heap.len() - 1;
        while hole > 0 {
            let parent = (hole - 1) / 2;
            if heap[parent].0 < heap[hole].0 {
                break;
            }
            heap.swap(parent, hole);
            hole = parent;
        }
    }
    fn pop(heap: &mut Vec<(HspContext, usize)>) -> Option<(HspContext, usize)> {
        let last = heap.pop()?;
        if heap.is_empty() {
            return Some(last);
        }
        let out = std::mem::replace(&mut heap[0], last);
        let mut hole = 0;
        while 2 * hole + 1 < heap.len() {
            let left = 2 * hole + 1;
            let right = left + 1;
            let child = if right < heap.len() && heap[right].0 < heap[left].0 {
                right
            } else {
                left
            };
            heap.swap(hole, child);
            hole = child;
        }
        while hole > 0 {
            let parent = (hole - 1) / 2;
            if heap[parent].0 < heap[hole].0 {
                break;
            }
            heap.swap(parent, hole);
            hole = parent;
        }
        Some(out)
    }
    for (index, reader) in readers.iter_mut().enumerate() {
        if let Some(value) = capture_io(read_context(reader), errors)? {
            push(&mut heap, (value, index));
        }
    }
    while let Some((value, index)) = pop(&mut heap) {
        callback(&value)?;
        if let Some(next) = capture_io(read_context(&mut readers[index]), errors)? {
            push(&mut heap, (next, index));
        }
    }
    drop(readers);
    for file in files {
        file.close_and_delete();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn context_wire_preserves_title_bytes_and_source_c_string_boundary() {
        let mut context = HspContext::default();
        context.query_title = b"q\xff\0ignored".to_vec();
        context.target_title = b"s\xfe".to_vec();
        context.query_oid = 7;
        context.subject_oid = 9;
        context.hsp.bit_score = 80.0;
        let mut wire = Vec::new();
        write_context(&mut wire, &context).unwrap();
        let decoded = read_context(&mut std::io::Cursor::new(wire))
            .unwrap()
            .unwrap();
        assert_eq!(decoded.query_title, b"q\xff");
        assert_eq!(decoded.target_title, b"s\xfe");
        assert_eq!((decoded.query_oid, decoded.subject_oid), (7, 9));
        assert_eq!(decoded.hsp.bit_score.to_bits(), 80.0_f64.to_bits());
    }
}
