//! Source-compatible MCL graph checkpoints and component-at-a-time access.
//!
//! The memory setting bounds the pending edge table, not an individual connected
//! component. As in Source, a component must fit in memory for matrix expansion.

use super::mcl::Triplet;
use std::collections::BTreeMap;
use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

pub struct SparseMatrixStream {
    path: PathBuf,
    temporary: bool,
    symmetric: bool,
    parents: Vec<usize>,
    pending: BTreeMap<(usize, usize), Triplet>,
    writer: Option<BufWriter<File>>,
    max_bytes: f64,
    stored_edges: u64,
    warned: bool,
}

impl SparseMatrixStream {
    pub fn create(
        sequence_count: usize,
        symmetric: bool,
        checkpoint: Option<&Path>,
        tmpdir: &Path,
        memory_gib: f64,
    ) -> io::Result<Self> {
        if sequence_count > u32::MAX as usize {
            return Err(invalid("MCL graph exceeds the 32-bit sequence index range"));
        }
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let (path, file) = if let Some(path) = checkpoint {
            (path.to_path_buf(), File::create(path)?)
        } else {
            loop {
                let path = tmpdir.join(format!(
                    "diamond-mcl-{}-{}.bin",
                    std::process::id(),
                    NEXT.fetch_add(1, Ordering::Relaxed)
                ));
                match OpenOptions::new().write(true).create_new(true).open(&path) {
                    Ok(file) => break (path, file),
                    Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                    Err(error) => return Err(error),
                }
            }
        };
        let mut writer = BufWriter::new(file);
        writer.write_all(&sequence_count.to_ne_bytes())?;
        writer.write_all(&0u32.to_ne_bytes())?;
        Ok(Self {
            path,
            temporary: checkpoint.is_none(),
            symmetric,
            parents: (0..sequence_count).collect(),
            pending: BTreeMap::new(),
            writer: Some(writer),
            stored_edges: 0,
            warned: false,
            max_bytes: (if memory_gib > 0.0 {
                memory_gib as f32
            } else {
                2.0
            }) as f64
                * 1024.0
                * 1024.0
                * 1024.0,
        })
    }

    pub fn open(path: &Path, symmetric: bool) -> io::Result<Self> {
        let mut reader = BufReader::new(File::open(path)?);
        let mut count = [0; std::mem::size_of::<usize>()];
        reader.read_exact(&mut count)?;
        let count = usize::from_ne_bytes(count);
        let mut version = [0; 4];
        reader.read_exact(&mut version)?;
        if u32::from_ne_bytes(version) != 0 {
            return Err(invalid("Unsupported MCL graph index version"));
        }
        if count > u32::MAX as usize {
            return Err(invalid("MCL graph exceeds the 32-bit sequence index range"));
        }
        let mut graph = Self {
            path: path.to_path_buf(),
            temporary: false,
            symmetric,
            parents: (0..count).collect(),
            pending: BTreeMap::new(),
            writer: None,
            max_bytes: 0.0,
            stored_edges: 0,
            warned: false,
        };
        visit_blocks(&mut reader, |edge| {
            graph.merge(edge.row, edge.col)?;
            graph.stored_edges += 1;
            Ok(())
        })?;
        Ok(graph)
    }

    pub fn sequence_count(&self) -> usize {
        self.parents.len()
    }

    pub fn edge_count(&self) -> u64 {
        self.stored_edges
    }

    /// Accept the exact 16-byte records emitted by `ClusteringFormat`.
    pub fn consume(&mut self, records: &[u8]) -> io::Result<()> {
        if self.writer.is_none() {
            return Err(invalid("MCL graph is already finalized"));
        }
        if records.len() % 16 != 0 {
            return Err(invalid("Truncated MCL graph record"));
        }
        for record in records.chunks_exact(16) {
            let row = u32::from_ne_bytes(record[..4].try_into().unwrap()) as usize;
            let col = u32::from_ne_bytes(record[4..8].try_into().unwrap()) as usize;
            let original = f64::from_ne_bytes(record[8..16].try_into().unwrap());
            if !self.warned && (original > f32::MAX as f64 || original < f32::MIN_POSITIVE as f64) {
                eprintln!("\nWARNING: The clustering similarity measure cannot be stored in a float, results may become unreliable\n         Please modify --clustering-similarity accordingly.\n");
                self.warned = true;
            }
            let value = original as f32;
            self.merge(row, col)?;
            let key = if self.symmetric {
                (row.max(col), row.min(col))
            } else {
                (row, col)
            };
            let edge = Triplet::new(row, col, value);
            match self.pending.get_mut(&key) {
                Some(old) if value > old.value => *old = edge,
                None => {
                    self.pending.insert(key, edge);
                }
                _ => (),
            }
            if self.pending.len() as f64 * 16.0 > self.max_bytes {
                self.dump()?;
            }
        }
        Ok(())
    }

    pub fn finish(&mut self) -> io::Result<()> {
        self.dump()?;
        if let Some(mut writer) = self.writer.take() {
            writer.flush()?;
        }
        Ok(())
    }

    pub fn components(&mut self) -> Vec<Vec<usize>> {
        let mut components: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for index in 0..self.parents.len() {
            let root = self.root(index);
            components.entry(root).or_default().push(index);
        }
        let components: Vec<_> = components.into_values().collect();
        super::source_order::components(components)
    }

    /// Read only one component into RAM. Source rereads spills into a fresh
    /// max-deduplicating set before constructing the component matrix.
    pub fn component(&self, indices: &[usize]) -> io::Result<Vec<Triplet>> {
        if self.writer.is_some() {
            return Err(invalid("MCL graph is not finalized"));
        }
        let remap: BTreeMap<_, _> = indices.iter().enumerate().map(|(i, &id)| (id, i)).collect();
        let mut reader = BufReader::new(File::open(&self.path)?);
        let mut header = vec![0; std::mem::size_of::<usize>() + 4];
        reader.read_exact(&mut header)?;
        let mut edges: BTreeMap<(usize, usize), Triplet> = BTreeMap::new();
        visit_blocks(&mut reader, |edge| {
            if let Some(&row) = remap.get(&edge.row) {
                let col = *remap
                    .get(&edge.col)
                    .ok_or_else(|| invalid("MCL graph component is inconsistent"))?;
                let key = if self.symmetric {
                    (edge.row.max(edge.col), edge.row.min(edge.col))
                } else {
                    (edge.row, edge.col)
                };
                match edges.get_mut(&key) {
                    Some(previous) if edge.value > previous.value => {
                        *previous = Triplet::new(row, col, edge.value)
                    }
                    None => {
                        edges.insert(key, Triplet::new(row, col, edge.value));
                    }
                    _ => (),
                }
            }
            Ok(())
        })?;
        Ok(edges.into_values().collect())
    }

    fn root(&mut self, mut node: usize) -> usize {
        while self.parents[node] != node {
            self.parents[node] = self.parents[self.parents[node]];
            node = self.parents[node];
        }
        node
    }

    fn merge(&mut self, row: usize, col: usize) -> io::Result<()> {
        if row >= self.parents.len() || col >= self.parents.len() {
            return Err(invalid("MCL graph edge is outside the sequence range"));
        }
        let a = self.root(row);
        let b = self.root(col);
        self.parents[b] = a;
        Ok(())
    }

    fn dump(&mut self) -> io::Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        // Source dump() calls get_indices() before splitting its ordered edge
        // set. Component IDs follow first physical members, while each block's
        // header is the first GNU unordered member, not the smallest OID or the
        // union-find representative. The later MCL work-size sort is separate.
        let mut by_root: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for oid in 0..self.parents.len() {
            let root = self.root(oid);
            by_root.entry(root).or_default().push(oid);
        }
        let components = super::source_order::indexed_components(by_root.into_values().collect());
        let mut index_to_component = vec![0; self.parents.len()];
        for (index, members) in components.iter().enumerate() {
            for &oid in members {
                index_to_component[oid] = index;
            }
        }
        let mut groups = vec![Vec::new(); components.len()];
        for edge in std::mem::take(&mut self.pending).into_values() {
            groups[index_to_component[edge.row]].push(edge);
        }
        let writer = self
            .writer
            .as_mut()
            .ok_or_else(|| invalid("MCL graph is already finalized"))?;
        for (members, edges) in components.into_iter().zip(groups) {
            if edges.is_empty() {
                continue;
            }
            let first = members[0];
            writer.write_all(&(first as i64).to_ne_bytes())?;
            let count =
                i32::try_from(edges.len()).map_err(|_| invalid("MCL graph block is too large"))?;
            writer.write_all(&count.to_ne_bytes())?;
            self.stored_edges += edges.len() as u64;
            for edge in edges {
                writer.write_all(&(edge.row as u32).to_ne_bytes())?;
                writer.write_all(&(edge.col as u32).to_ne_bytes())?;
                writer.write_all(&(edge.value as f64).to_ne_bytes())?;
            }
        }
        writer.flush()
    }
}

impl Drop for SparseMatrixStream {
    fn drop(&mut self) {
        self.writer.take();
        if self.temporary {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

fn visit_blocks(
    reader: &mut impl Read,
    mut visit: impl FnMut(Triplet) -> io::Result<()>,
) -> io::Result<()> {
    loop {
        let mut first = [0; 8];
        match reader.read(&mut first[..1])? {
            0 => return Ok(()),
            _ => reader.read_exact(&mut first[1..])?,
        }
        let mut count = [0; 4];
        reader.read_exact(&mut count)?;
        let count = i32::from_ne_bytes(count);
        if count < 0 {
            return Err(invalid("Negative MCL graph block size"));
        }
        for _ in 0..count {
            let mut record = [0; 16];
            reader.read_exact(&mut record)?;
            visit(Triplet::new(
                u32::from_ne_bytes(record[..4].try_into().unwrap()) as usize,
                u32::from_ne_bytes(record[4..8].try_into().unwrap()) as usize,
                f64::from_ne_bytes(record[8..].try_into().unwrap()) as f32,
            ))?;
        }
    }
}
