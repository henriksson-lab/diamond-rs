use crate::basic::value::TaxId;
use crate::util::io::source_c_path;
use std::fs::File;
use std::io::{BufRead, BufReader, Error, ErrorKind};

pub fn read_nodes_dmp_path<F>(file_name: &std::path::Path, mut f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &str),
{
    read_nodes_bytes_reader(
        BufReader::new(File::open(source_c_path(&(file_name)).as_ref())?),
        |id, parent, rank| {
            f(
                id,
                parent,
                std::str::from_utf8(rank)
                    .map_err(|error| Error::new(ErrorKind::InvalidData, error))?,
            );
            Ok(())
        },
    )
}
pub fn read_nodes_dmp_bytes_path<F>(file_name: &std::path::Path, mut f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &[u8]),
{
    read_nodes_bytes_reader(
        BufReader::new(File::open(source_c_path(&(file_name)).as_ref())?),
        |id, parent, rank| {
            f(id, parent, rank);
            Ok(())
        },
    )
}
pub fn read_nodes_dmp_source_path<F>(file_name: &std::path::Path, f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &str),
{
    read_nodes_dmp_with_source_buffer_path(file_name, 1 << 20, f)
}
pub fn read_nodes_dmp_with_source_buffer_path<F>(
    file_name: &std::path::Path,
    buffer_size: usize,
    mut f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &str),
{
    read_nodes_source_impl(file_name, buffer_size, |id, parent, rank| {
        f(
            id,
            parent,
            std::str::from_utf8(rank).map_err(|error| Error::new(ErrorKind::InvalidData, error))?,
        );
        Ok(())
    })
}
pub fn read_nodes_dmp_bytes_with_source_buffer_path<F>(
    file_name: &std::path::Path,
    buffer_size: usize,
    mut f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &[u8]),
{
    read_nodes_source_impl(file_name, buffer_size, |id, parent, rank| {
        f(id, parent, rank);
        Ok(())
    })
}
pub(crate) fn read_nodes_dmp_bytes_io_with_source_buffer_path<F>(
    file_name: &std::path::Path,
    buffer_size: usize,
    f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &[u8]) -> std::io::Result<()>,
{
    read_nodes_source_impl(file_name, buffer_size, f)
}
pub(crate) fn read_nodes_dmp_bytes_io_path<F>(
    file_name: &std::path::Path,
    f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &[u8]) -> std::io::Result<()>,
{
    read_nodes_bytes_reader(
        BufReader::new(File::open(source_c_path(&(file_name)).as_ref())?),
        f,
    )
}
fn read_nodes_source_impl<F>(
    file_name: &std::path::Path,
    buffer_size: usize,
    mut f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &[u8]) -> std::io::Result<()>,
{
    let input = crate::util::io::InputFile::with_source_profile_path(file_name, 0, buffer_size)
        .map_err(Error::other)?;
    let mut input =
        crate::util::io::ByteTextInputFile::with_default_separator_bytes(input.into_deserializer());
    while !input.eof() {
        input.getline().map_err(Error::other)?;
        if input.line.is_empty() {
            break;
        }
        nodes_row_bytes(&input.line, &mut f)?;
    }
    input.into_inner().close().map_err(Error::other)
}
fn read_nodes_bytes_reader<R: BufRead, F>(mut reader: R, mut f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &[u8]) -> std::io::Result<()>,
{
    loop {
        let mut line = Vec::new();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        if line.is_empty() {
            break;
        }
        nodes_row_bytes(&line, &mut f)?;
    }
    Ok(())
}
fn string_fields(line: &[u8]) -> Vec<&[u8]> {
    // Source StringDelimiter scans a C string, then returns std::string fields.
    let line = &line[..line
        .iter()
        .position(|&byte| byte == 0)
        .unwrap_or(line.len())];
    let mut fields = Vec::new();
    let mut begin = 0;
    while let Some(offset) = line[begin..].windows(3).position(|part| part == b"\t|\t") {
        let end = begin + offset;
        fields.push(&line[begin..end]);
        begin = end + 3;
    }
    fields.push(&line[begin..]);
    fields
}
fn field<'a>(fields: &[&'a [u8]], i: usize) -> std::io::Result<&'a [u8]> {
    fields.get(i).copied().ok_or_else(|| {
        Error::new(
            ErrorKind::InvalidData,
            crate::util::string::TokenizerException::new(),
        )
    })
}
fn field_taxid(bytes: &[u8]) -> std::io::Result<TaxId> {
    crate::util::string::ByteTokenizer::new(bytes, 0)
        .read_i32()
        .map_err(|error| Error::new(ErrorKind::InvalidData, error))
}
fn nodes_row_bytes<F>(line: &[u8], f: &mut F) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &[u8]) -> std::io::Result<()>,
{
    let fields = string_fields(line);
    let id = field_taxid(field(&fields, 0)?)?;
    let parent = field_taxid(field(&fields, 1)?)?;
    f(id, parent, field(&fields, 2)?)
}

pub fn read_names_dmp_path<F>(file_name: &std::path::Path, f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, &str),
{
    read_names_reader(
        BufReader::new(File::open(source_c_path(&(file_name)).as_ref())?),
        f,
    )
}
pub fn read_names_dmp_source_path<F>(file_name: &std::path::Path, f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, &str),
{
    read_names_dmp_with_source_buffer_path(file_name, 1 << 20, f)
}
pub fn read_names_dmp_with_source_buffer_path<F>(
    file_name: &std::path::Path,
    buffer_size: usize,
    mut f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, &str),
{
    read_names_source_impl(file_name, buffer_size, |id, bytes| {
        let name = std::str::from_utf8(bytes)
            .map_err(|error| Error::new(ErrorKind::InvalidData, error))?;
        f(id, name);
        Ok(())
    })
}
/// Actual Source TextInputFile/Tokenizer byte path. A final unterminated line
/// sets eof before the loop body and is therefore not consumed by Source.
pub fn read_names_dmp_bytes_with_source_buffer_path<F>(
    file_name: &std::path::Path,
    buffer_size: usize,
    mut f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, &[u8]),
{
    read_names_source_impl(file_name, buffer_size, |id, bytes| {
        f(id, bytes);
        Ok(())
    })
}
fn read_names_source_impl<F>(
    file_name: &std::path::Path,
    buffer_size: usize,
    mut f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, &[u8]) -> std::io::Result<()>,
{
    let input = crate::util::io::InputFile::with_source_profile_path(file_name, 0, buffer_size)
        .map_err(Error::other)?;
    let mut input =
        crate::util::io::ByteTextInputFile::with_default_separator_bytes(input.into_deserializer());
    loop {
        input.getline().map_err(Error::other)?;
        if input.eof() {
            break;
        }
        if !input.line.is_empty() {
            names_row_bytes(&input.line, &mut f)?;
        }
    }
    input.into_inner().close().map_err(Error::other)
}
fn names_row_bytes<F>(line: &[u8], f: &mut F) -> std::io::Result<()>
where
    F: FnMut(TaxId, &[u8]) -> std::io::Result<()>,
{
    let fields = string_fields(line);
    let id = field_taxid(field(&fields, 0)?)?;
    let name = field(&fields, 1)?;
    let _unique = field(&fields, 2)?;
    let mut kind = field(&fields, 3)?;
    while kind.last().is_some_and(|byte| matches!(byte, b'\t' | b'|')) {
        kind = &kind[..kind.len() - 1];
    }
    if kind == b"scientific name" {
        f(id, name)?;
    }
    Ok(())
}
fn read_names_reader<R: BufRead, F>(mut reader: R, mut f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, &str),
{
    loop {
        let mut line = Vec::new();
        let len = reader.read_until(b'\n', &mut line)?;
        if len == 0 {
            break;
        }
        if line.last() != Some(&b'\n') {
            break;
        }
        line.pop();
        if line.last() == Some(&b'\r') {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }
        names_row_bytes(&line, &mut |id, bytes| {
            f(
                id,
                std::str::from_utf8(bytes)
                    .map_err(|error| Error::new(ErrorKind::InvalidData, error))?,
            );
            Ok(())
        })?;
    }
    Ok(())
}

pub fn read_nodes_dmp<F>(file_name: &str, f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &str),
{
    read_nodes_dmp_path(std::path::Path::new(file_name), f)
}
pub fn read_nodes_dmp_bytes<F>(file_name: &str, f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &[u8]),
{
    read_nodes_dmp_bytes_path(std::path::Path::new(file_name), f)
}
pub fn read_nodes_dmp_source<F>(file_name: &str, f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &str),
{
    read_nodes_dmp_source_path(std::path::Path::new(file_name), f)
}
pub fn read_nodes_dmp_with_source_buffer<F>(
    file_name: &str,
    buffer_size: usize,
    f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &str),
{
    read_nodes_dmp_with_source_buffer_path(std::path::Path::new(file_name), buffer_size, f)
}
pub fn read_nodes_dmp_bytes_with_source_buffer<F>(
    file_name: &str,
    buffer_size: usize,
    f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &[u8]),
{
    read_nodes_dmp_bytes_with_source_buffer_path(std::path::Path::new(file_name), buffer_size, f)
}
pub(crate) fn read_nodes_dmp_bytes_io_with_source_buffer<F>(
    file_name: &str,
    buffer_size: usize,
    f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &[u8]) -> std::io::Result<()>,
{
    read_nodes_dmp_bytes_io_with_source_buffer_path(std::path::Path::new(file_name), buffer_size, f)
}
pub(crate) fn read_nodes_dmp_bytes_io<F>(file_name: &str, f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, TaxId, &[u8]) -> std::io::Result<()>,
{
    read_nodes_dmp_bytes_io_path(std::path::Path::new(file_name), f)
}
pub fn read_names_dmp<F>(file_name: &str, f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, &str),
{
    read_names_dmp_path(std::path::Path::new(file_name), f)
}
pub fn read_names_dmp_source<F>(file_name: &str, f: F) -> std::io::Result<()>
where
    F: FnMut(TaxId, &str),
{
    read_names_dmp_source_path(std::path::Path::new(file_name), f)
}
pub fn read_names_dmp_with_source_buffer<F>(
    file_name: &str,
    buffer_size: usize,
    f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, &str),
{
    read_names_dmp_with_source_buffer_path(std::path::Path::new(file_name), buffer_size, f)
}
pub fn read_names_dmp_bytes_with_source_buffer<F>(
    file_name: &str,
    buffer_size: usize,
    f: F,
) -> std::io::Result<()>
where
    F: FnMut(TaxId, &[u8]),
{
    read_names_dmp_bytes_with_source_buffer_path(std::path::Path::new(file_name), buffer_size, f)
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_read_nodes_dmp_template_callback() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-nodes-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "1\t|\t1\t|\tno rank\t|\t").unwrap();
            writeln!(f, "2\t|\t1\t|\tsuperkingdom\t|\t").unwrap();
        }

        let mut rows = Vec::new();
        read_nodes_dmp(path.to_str().unwrap(), |taxid, parent, rank| {
            rows.push((taxid, parent, rank.to_string()));
        })
        .unwrap();
        assert_eq!(
            rows,
            vec![
                (1, 1, "no rank".to_string()),
                (2, 1, "superkingdom".to_string())
            ]
        );
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_read_nodes_dmp_stops_at_empty_line() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-nodes-empty-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "1\t|\t1\t|\tno rank\t|\t").unwrap();
            writeln!(f).unwrap();
            writeln!(f, "2\t|\t1\t|\tsuperkingdom\t|\t").unwrap();
        }

        let mut rows = Vec::new();
        read_nodes_dmp(path.to_str().unwrap(), |taxid, parent, rank| {
            rows.push((taxid, parent, rank.to_string()));
        })
        .unwrap();
        assert_eq!(rows, vec![(1, 1, "no rank".to_string())]);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_read_nodes_dmp_rejects_malformed_taxid() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-nodes-bad-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "not-a-taxid\t|\t1\t|\tno rank\t|\t").unwrap();
        }

        let err = read_nodes_dmp(path.to_str().unwrap(), |_, _, _| {}).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn source_names_keep_raw_bytes_and_skip_unterminated_final_line() {
        let path = std::env::temp_dir().join(format!(
            "diamond-raw-taxnames-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, b"2\t|\tName\xff\t|\t\t|\tscientific name\t|\r\n3\t|\tskipped\t|\t\t|\tscientific name\t|").unwrap();
        let mut rows = Vec::new();
        read_names_dmp_bytes_with_source_buffer(path.to_str().unwrap(), 32, |id, name| {
            rows.push((id, name.to_vec()))
        })
        .unwrap();
        assert_eq!(rows, vec![(2, b"Name\xff".to_vec())]);
        assert!(read_names_dmp_with_source_buffer(path.to_str().unwrap(), 32, |_, _| {}).is_err());
        let mut names = crate::data::taxonomy::ByteTaxonomy::default();
        assert_eq!(
            names
                .init_source_bytes_io(path.to_str().unwrap(), 32)
                .unwrap(),
            1
        );
        assert_eq!(names.name[2], b"Name\xff");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_read_names_dmp_filters_scientific_name() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-names-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "2\t|\tBacteria\t|\t\t|\tscientific name\t|").unwrap();
            writeln!(f, "2\t|\tbacteria\t|\t\t|\tcommon name\t|").unwrap();
            writeln!(f, "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|").unwrap();
        }

        let mut rows = Vec::new();
        read_names_dmp(path.to_str().unwrap(), |id, name| {
            rows.push((id, name.to_string()));
        })
        .unwrap();
        assert_eq!(
            rows,
            vec![
                (2, "Bacteria".to_string()),
                (9606, "Homo sapiens".to_string())
            ]
        );
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_read_names_dmp_rejects_malformed_taxid() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-names-bad-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "bad\t|\tBacteria\t|\t\t|\tscientific name\t|").unwrap();
        }

        let err = read_names_dmp(path.to_str().unwrap(), |_, _| {}).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::InvalidData);
        std::fs::remove_file(path).unwrap();
    }
}
