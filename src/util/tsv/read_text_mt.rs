//! Record-boundary-aware chunking for TSV text.
//!
//! This mirrors `diamond/src/util/tsv/read_text_mt.cpp`. The upstream reader
//! uses a reader thread and one or more consumer threads; this slice-backed
//! facade invokes the mutable callback synchronously in assigned chunk order.

pub const READ_TEXT_MT_SIZE: usize = 1 << 20;

/// Read text in 1 MiB blocks, extending full blocks through the next newline.
///
/// A newline used to extend a full block is consumed but excluded from the
/// callback slice, matching `Deserializer::read_to`. Newlines already inside a
/// short final block remain part of that slice.
pub fn read_text_mt<F>(data: &[u8], max_size: i64, _threads: usize, mut callback: F)
where
    F: FnMut(i64, &[u8]),
{
    let mut next_chunk = 0i64;
    let mut total = 0i64;
    let mut position = 0usize;

    loop {
        if position >= data.len() {
            break;
        }

        let raw_end = (position + READ_TEXT_MT_SIZE).min(data.len());
        let mut end = raw_end;
        if raw_end - position == READ_TEXT_MT_SIZE {
            while end < data.len() && data[end] != b'\n' {
                end += 1;
            }
        }

        let size = end - position;
        total += size as i64;
        if size > 0 {
            callback(next_chunk, &data[position..end]);
            next_chunk += 1;
        }

        if size < READ_TEXT_MT_SIZE || total + READ_TEXT_MT_SIZE as i64 > max_size {
            break;
        }
        position = if end < data.len() && data[end] == b'\n' {
            end + 1
        } else {
            end
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_input_never_calls_callback() {
        let mut calls = 0;
        read_text_mt(b"", i64::MAX, 8, |_, _| calls += 1);
        assert_eq!(calls, 0);
    }

    #[test]
    fn short_final_read_preserves_all_bytes_and_newlines() {
        let mut chunks = Vec::new();
        read_text_mt(b"a\nb\n", i64::MAX, 3, |chunk, bytes| {
            chunks.push((chunk, bytes.to_vec()));
        });
        assert_eq!(chunks, vec![(0, b"a\nb\n".to_vec())]);
    }

    #[test]
    fn full_block_extends_to_and_consumes_boundary_newline() {
        let mut data = vec![b'a'; READ_TEXT_MT_SIZE];
        data.extend_from_slice(b"tail\nnext\n");
        let mut chunks = Vec::new();

        read_text_mt(&data, i64::MAX, 4, |chunk, bytes| {
            chunks.push((chunk, bytes.len(), bytes.to_vec()));
        });

        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].0, 0);
        assert_eq!(chunks[0].1, READ_TEXT_MT_SIZE + 4);
        assert!(chunks[0].2.ends_with(b"tail"));
        assert_eq!(chunks[1], (1, 5, b"next\n".to_vec()));
    }

    #[test]
    fn exact_block_at_eof_is_emitted_once() {
        let data = vec![b'x'; READ_TEXT_MT_SIZE];
        let mut sizes = Vec::new();
        read_text_mt(&data, i64::MAX, 0, |chunk, bytes| {
            sizes.push((chunk, bytes.len()));
        });
        assert_eq!(sizes, vec![(0, READ_TEXT_MT_SIZE)]);
    }

    #[test]
    fn max_size_is_checked_between_full_extended_blocks() {
        let mut data = Vec::new();
        for byte in [b'a', b'b', b'c'] {
            data.extend(std::iter::repeat_n(byte, READ_TEXT_MT_SIZE));
            data.push(b'\n');
        }
        let mut chunks = Vec::new();

        read_text_mt(&data, (2 * READ_TEXT_MT_SIZE) as i64, 2, |chunk, bytes| {
            chunks.push((chunk, bytes[0], bytes.len()))
        });

        assert_eq!(
            chunks,
            vec![(0, b'a', READ_TEXT_MT_SIZE), (1, b'b', READ_TEXT_MT_SIZE)]
        );
    }

    #[test]
    fn tiny_maximum_still_delivers_first_complete_record_block() {
        let mut data = vec![b'z'; READ_TEXT_MT_SIZE + 7];
        data.push(b'\n');
        data.extend_from_slice(b"later\n");
        let mut chunks = Vec::new();

        read_text_mt(&data, 0, 1, |chunk, bytes| {
            chunks.push((chunk, bytes.len(), bytes.last().copied()));
        });

        assert_eq!(chunks, vec![(0, READ_TEXT_MT_SIZE + 7, Some(b'z'))]);
    }
}
