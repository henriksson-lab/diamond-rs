use std::collections::{BTreeMap, HashMap};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BitVector {
    data: Vec<u64>,
    size: u64,
}

pub trait ModuloPolicy {
    fn modulo(x: u64, y: u64) -> u64;
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Modulo2;

impl ModuloPolicy for Modulo2 {
    fn modulo(x: u64, y: u64) -> u64 {
        x & (y - 1)
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Modulo;

impl ModuloPolicy for Modulo {
    fn modulo(x: u64, y: u64) -> u64 {
        x % y
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NoModulo;

impl ModuloPolicy for NoModulo {
    fn modulo(x: u64, _y: u64) -> u64 {
        x
    }
}

pub trait U64Hash {
    fn hash(x: u64) -> u64;
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Identity;

impl U64Hash for Identity {
    fn hash(x: u64) -> u64 {
        x
    }
}

impl U64Hash for crate::util::hash::MurmurHash {
    fn hash(x: u64) -> u64 {
        crate::util::hash::murmur_hash_u64(x)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HashSet<M = Modulo, H = Identity> {
    table: Vec<u8>,
    size: usize,
    _marker: std::marker::PhantomData<(M, H)>,
}

impl<M, H> Default for HashSet<M, H> {
    fn default() -> Self {
        Self {
            table: Vec::new(),
            size: 0,
            _marker: std::marker::PhantomData,
        }
    }
}

impl<M, H> HashSet<M, H>
where
    M: ModuloPolicy,
    H: U64Hash,
{
    pub const PADDING: usize = 16;

    pub fn new(size: usize) -> Self {
        Self {
            table: vec![0; size + Self::PADDING],
            size,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn from_data(data: Vec<u8>, size: usize) -> Self {
        Self {
            table: data,
            size,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn contains(&self, key: u64) -> bool {
        let (found, _) = self.get_entry(key).expect("hash table overflow");
        found
    }

    pub fn insert(&mut self, key: u64) {
        let hash = H::hash(key);
        let f = Self::finger_print(hash);
        let (_, entry) = self.get_entry(key).expect("hash table overflow");
        if self.table[entry] == 0 {
            self.table[entry] = f;
        }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn load(&self) -> usize {
        self.table[..self.size].iter().filter(|&&x| x != 0).count()
    }

    pub fn data(&self) -> &[u8] {
        &self.table
    }

    pub fn table(&self) -> &[u8] {
        &self.table
    }

    pub fn table_mut(&mut self) -> &mut [u8] {
        &mut self.table
    }

    pub fn finish(&mut self) {
        for i in 0..Self::PADDING {
            self.table[self.size + i] = self.table[i];
        }
    }

    fn finger_print(hash: u64) -> u8 {
        let x = (hash & ((1u64 << 8) - 1)) as u8;
        x.max(1)
    }

    fn get_entry(&self, key: u64) -> Result<(bool, usize), &'static str> {
        let hash = H::hash(key);
        let f = Self::finger_print(hash);
        let mut p = M::modulo(hash >> 8, self.size as u64) as usize;
        let mut wrapped = false;
        loop {
            if self.table[p] == f {
                return Ok((true, p));
            }
            if self.table[p] == 0 {
                return Ok((false, p));
            }
            p += 1;
            if p == self.size {
                if wrapped {
                    return Err("Hash table overflow");
                }
                p = 0;
                wrapped = true;
            }
        }
    }
}

pub trait HashTableHash<Key> {
    fn hash(key: Key) -> usize;
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct IdentityHash;

impl HashTableHash<usize> for IdentityHash {
    fn hash(key: usize) -> usize {
        key
    }
}

impl HashTableHash<u64> for IdentityHash {
    fn hash(key: u64) -> usize {
        key as usize
    }
}

impl HashTableHash<u32> for IdentityHash {
    fn hash(key: u32) -> usize {
        key as usize
    }
}

impl HashTableHash<u64> for crate::util::hash::MurmurHash {
    fn hash(key: u64) -> usize {
        crate::util::hash::murmur_hash_u64(key) as usize
    }
}

impl HashTableHash<u32> for crate::util::hash::MurmurHash {
    fn hash(key: u32) -> usize {
        crate::util::hash::murmur_hash_u64(key as u64) as usize
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HashTableEntry<Key, Value> {
    pub key: Key,
    pub value: Value,
}

impl<Key, Value> HashTableEntry<Key, Value>
where
    Value: Default + PartialEq,
{
    pub fn blank(&self) -> bool {
        self.value == Value::default()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HashTable<Key, Value, H = IdentityHash, M = Modulo> {
    table: Vec<HashTableEntry<Key, Value>>,
    size: usize,
    _marker: std::marker::PhantomData<(H, M)>,
}

impl<Key, Value, H, M> HashTable<Key, Value, H, M>
where
    Key: Copy + Default + PartialEq,
    Value: Copy + Default + PartialEq,
    H: HashTableHash<Key>,
    M: ModuloPolicy,
{
    pub fn new(size: usize) -> Self {
        Self {
            table: vec![
                HashTableEntry {
                    key: Key::default(),
                    value: Value::default(),
                };
                size
            ],
            size,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn get_mut(&mut self, key: Key) -> Result<&mut Value, &'static str> {
        let idx = self.get_entry(key)?;
        if self.table[idx].blank() {
            self.table[idx].key = key;
        }
        Ok(&mut self.table[idx].value)
    }

    pub fn find(&mut self, key: Key) -> Result<Option<&mut Value>, &'static str> {
        let idx = self.get_present_entry(key)?;
        Ok(idx.map(|i| &mut self.table[i].value))
    }

    pub fn find_entry(
        &mut self,
        key: Key,
    ) -> Result<Option<&mut HashTableEntry<Key, Value>>, &'static str> {
        let idx = self.get_present_entry(key)?;
        Ok(idx.map(|i| &mut self.table[i]))
    }

    pub fn insert(&mut self, key: Key) -> Result<&mut HashTableEntry<Key, Value>, &'static str> {
        let idx = self.get_or_insert_entry(key)?;
        Ok(&mut self.table[idx])
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn count(&self) -> usize {
        self.table.iter().filter(|entry| !entry.blank()).count()
    }

    pub fn data(&self) -> &[HashTableEntry<Key, Value>] {
        &self.table
    }

    pub fn data_mut(&mut self) -> &mut [HashTableEntry<Key, Value>] {
        &mut self.table
    }

    fn start(&self, key: Key) -> usize {
        M::modulo(H::hash(key) as u64, self.size as u64) as usize
    }

    fn get_entry(&mut self, key: Key) -> Result<usize, &'static str> {
        let mut p = self.start(key);
        let mut wrapped = false;
        while self.table[p].key != key && !self.table[p].blank() {
            p += 1;
            if p == self.size {
                if wrapped {
                    return Err("Hash table overflow.");
                }
                p = 0;
                wrapped = true;
            }
        }
        Ok(p)
    }

    fn get_present_entry(&mut self, key: Key) -> Result<Option<usize>, &'static str> {
        let mut p = self.start(key);
        let mut wrapped = false;
        loop {
            if self.table[p].blank() {
                return Ok(None);
            }
            if self.table[p].key == key {
                return Ok(Some(p));
            }
            p += 1;
            if p == self.size {
                if wrapped {
                    return Err("Hash table overflow.");
                }
                p = 0;
                wrapped = true;
            }
        }
    }

    fn get_or_insert_entry(&mut self, key: Key) -> Result<usize, &'static str> {
        let mut p = self.start(key);
        let mut wrapped = false;
        loop {
            if self.table[p].key == key {
                return Ok(p);
            }
            if self.table[p].blank() {
                self.table[p].key = key;
                return Ok(p);
            }
            p += 1;
            if p == self.size {
                if wrapped {
                    return Err("Hash table overflow.");
                }
                p = 0;
                wrapped = true;
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Queue<T> {
    capacity: usize,
    producer_count: i32,
    consumer_count: i32,
    poison_pill: T,
    data: std::collections::VecDeque<T>,
    pills_received: usize,
    closed: bool,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PtrVector<T> {
    data: Vec<Box<T>>,
}

impl<T> PtrVector<T> {
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    pub fn with_size(n: usize) -> Self
    where
        T: Default,
    {
        let mut data = Vec::with_capacity(n);
        for _ in 0..n {
            data.push(Box::<T>::default());
        }
        Self { data }
    }

    pub fn push(&mut self, value: T) {
        self.data.push(Box::new(value));
    }

    pub fn get_ptr(&mut self, n: usize) -> &mut Box<T> {
        &mut self.data[n]
    }

    pub fn get(&mut self, n: usize) -> &mut Box<T> {
        self.get_ptr(n)
    }

    pub fn back(&mut self) -> &mut T {
        self.data.last_mut().unwrap()
    }

    pub fn erase(&mut self, first: usize, last: usize) {
        self.data.drain(first..last);
    }

    pub fn clear(&mut self) {
        self.data.clear();
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

impl<T> std::ops::Index<usize> for PtrVector<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<T> std::ops::IndexMut<usize> for PtrVector<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl<T: Clone + PartialEq> Queue<T> {
    pub fn new(capacity: usize, producer_count: i32, consumer_count: i32, poison_pill: T) -> Self {
        assert!(producer_count == 1 || consumer_count == 1);
        Self {
            capacity: Self::round_up_to_power_of_two(if capacity != 0 { capacity } else { 1 }),
            producer_count,
            consumer_count,
            poison_pill,
            data: std::collections::VecDeque::new(),
            pills_received: 0,
            closed: false,
        }
    }

    pub const fn is_power_of_two(x: usize) -> bool {
        (x & (x - 1)) == 0
    }

    pub fn round_up_to_power_of_two(mut x: usize) -> usize {
        if Self::is_power_of_two(x) {
            return x;
        }
        x -= 1;
        let mut i = 1usize;
        while i < std::mem::size_of::<usize>() * 8 {
            x |= x >> i;
            i <<= 1;
        }
        x + 1
    }

    pub fn enqueue(&mut self, v: T) {
        assert!(self.data.len() < self.capacity);
        self.data.push_back(v);
    }

    pub fn wait_and_dequeue(&mut self) -> Option<T> {
        loop {
            let out = self.data.pop_front()?;
            if out == self.poison_pill {
                if self.producer_count > 1 {
                    self.pills_received += 1;
                    if self.pills_received == self.producer_count as usize {
                        return None;
                    }
                    continue;
                }
                return None;
            }
            return Some(out);
        }
    }

    pub fn empty(&self) -> bool {
        self.approx_size() == 0
    }

    pub fn capacity(&self) -> usize {
        self.capacity
    }

    pub fn approx_size(&self) -> usize {
        self.data.len()
    }

    pub fn close(&mut self) {
        for _ in 0..self.consumer_count {
            self.enqueue(self.poison_pill.clone());
        }
        self.closed = true;
    }

    pub fn producer_count(&self) -> i32 {
        self.producer_count
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Deque<T, const E: usize> {
    buckets: Vec<Vec<T>>,
    total: usize,
}

impl<T, const E: usize> Default for Deque<T, E> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T, const E: usize> Deque<T, E> {
    pub const EXPONENT: usize = E;
    pub const N: usize = 1usize << E;
    pub const SHIFT: isize = E as isize;
    pub const MASK: isize = ((1isize) << E) - 1;

    pub fn new() -> Self {
        let mut out = Self {
            buckets: Vec::new(),
            total: 0,
        };
        out.new_bucket();
        out
    }

    pub fn reserve(&mut self, _n: usize) {}

    pub fn push_back(&mut self, x: T) {
        if self.buckets.last().unwrap().len() >= Self::N {
            self.new_bucket();
        }
        self.buckets.last_mut().unwrap().push(x);
    }

    pub fn push_back_slice(&mut self, ptr: &[T])
    where
        T: Clone,
    {
        if self.buckets.last().unwrap().len() + ptr.len() > Self::N {
            self.new_bucket();
        }
        self.buckets.last_mut().unwrap().extend_from_slice(ptr);
    }

    pub fn push_back_iter(&mut self, mut begin: &[T])
    where
        T: Clone,
    {
        while !begin.is_empty() {
            if self.buckets.last().unwrap().len() == Self::N {
                self.new_bucket();
            }
            let n = begin
                .len()
                .min(Self::N - self.buckets.last().unwrap().len());
            self.buckets
                .last_mut()
                .unwrap()
                .extend_from_slice(&begin[..n]);
            begin = &begin[n..];
        }
    }

    pub fn size(&self) -> usize {
        let mut n = 0;
        for b in &self.buckets {
            n += b.len();
        }
        n
    }

    pub fn move_to(&mut self, dst: &mut Vec<T>)
    where
        T: Clone,
    {
        if self.buckets.len() == 1 && dst.is_empty() {
            *dst = std::mem::take(&mut self.buckets[0]);
        } else {
            for b in &self.buckets {
                dst.extend_from_slice(b);
            }
        }
        self.buckets.clear();
    }

    pub fn begin(&mut self) -> DequeIterator<T, E>
    where
        T: Clone,
    {
        self.init();
        DequeIterator::new(0, self.flattened())
    }

    pub fn end(&mut self) -> DequeIterator<T, E>
    where
        T: Clone,
    {
        self.init();
        DequeIterator::new(self.total, self.flattened())
    }

    fn init(&mut self) {
        self.total = self.size();
    }

    fn flattened(&self) -> Vec<T>
    where
        T: Clone,
    {
        let mut out = Vec::with_capacity(self.size());
        for b in &self.buckets {
            out.extend_from_slice(b);
        }
        out
    }

    fn new_bucket(&mut self) {
        self.buckets.push(Vec::with_capacity(Self::N));
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DequeIterator<T, const E: usize> {
    i: isize,
    data: Vec<T>,
}

impl<T, const E: usize> DequeIterator<T, E> {
    fn new(i: usize, data: Vec<T>) -> Self {
        Self {
            i: i as isize,
            data,
        }
    }

    pub fn get(&self) -> &T {
        &self.data[self.i as usize]
    }

    pub fn get_offset(&self, i: isize) -> &T {
        &self.data[(self.i + i) as usize]
    }

    pub fn increment(&mut self) -> &mut Self {
        self.i += 1;
        self
    }

    pub fn decrement(&mut self) -> &mut Self {
        self.i -= 1;
        self
    }

    pub fn add(&self, i: isize) -> Self
    where
        T: Clone,
    {
        Self {
            i: self.i + i,
            data: self.data.clone(),
        }
    }

    pub fn sub(&self, i: isize) -> Self
    where
        T: Clone,
    {
        Self {
            i: self.i - i,
            data: self.data.clone(),
        }
    }

    pub fn distance(&self, it: &Self) -> isize {
        self.i - it.i
    }
}

pub struct AsyncWriter<'a, T: Clone, const E: usize> {
    dst: &'a mut Deque<T, E>,
    buf: Vec<T>,
}

impl<'a, T: Clone, const E: usize> AsyncWriter<'a, T, E> {
    const BUF_SIZE: usize = 4096;

    pub fn new(dst: &'a mut Deque<T, E>) -> Self {
        Self {
            dst,
            buf: Vec::new(),
        }
    }

    pub fn write(&mut self, v: T) {
        self.buf.push(v);
        if self.buf.len() >= Self::BUF_SIZE {
            self.dst.push_back_iter(&self.buf);
            self.buf.clear();
        }
    }
}

impl<'a, T: Clone, const E: usize> Drop for AsyncWriter<'a, T, E> {
    fn drop(&mut self) {
        self.dst.push_back_iter(&self.buf);
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MemBuffer<T> {
    data: Vec<T>,
    alloc_size: usize,
}

impl<T> Default for MemBuffer<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> MemBuffer<T> {
    pub const ALIGN: usize = 32;

    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            alloc_size: 0,
        }
    }

    pub fn with_size(n: usize) -> Self
    where
        T: Default + Clone,
    {
        Self {
            data: vec![T::default(); n],
            alloc_size: n,
        }
    }

    pub fn resize(&mut self, n: usize)
    where
        T: Default + Clone,
    {
        if self.alloc_size < n {
            self.data = vec![T::default(); n];
            self.alloc_size = n;
        } else {
            self.data.resize(n, T::default());
        }
    }

    pub fn size(&self) -> usize {
        self.data.len()
    }

    pub fn begin(&self) -> &[T] {
        &self.data
    }

    pub fn begin_mut(&mut self) -> &mut [T] {
        &mut self.data
    }

    pub fn end(&self) -> &[T] {
        &self.data[self.data.len()..]
    }

    pub fn end_mut(&mut self) -> &mut [T] {
        let len = self.data.len();
        &mut self.data[len..]
    }

    pub fn data(&self) -> &[T] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut [T] {
        &mut self.data
    }

    pub fn get(&self, i: usize) -> &T {
        &self.data[i]
    }

    pub fn get_mut(&mut self, i: usize) -> &mut T {
        &mut self.data[i]
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleArray {
    data: Vec<u8>,
    size: usize,
}

impl Default for DoubleArray {
    fn default() -> Self {
        Self::new()
    }
}

impl DoubleArray {
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            size: 0,
        }
    }

    pub fn from_bytes(data: Vec<u8>, size: usize) -> Self {
        Self { data, size }
    }

    pub fn begin(&self) -> DoubleArrayIterator<'_> {
        self.begin_with_elem_size(4)
    }

    pub fn begin_with_elem_size(&self, elem_size: usize) -> DoubleArrayIterator<'_> {
        DoubleArrayIterator::new(&self.data, 0, self.size, elem_size)
    }

    pub fn begin_mut(&mut self) -> DoubleArrayIteratorMut<'_> {
        self.begin_mut_with_elem_size(4)
    }

    pub fn begin_mut_with_elem_size(&mut self, elem_size: usize) -> DoubleArrayIteratorMut<'_> {
        DoubleArrayIteratorMut::new(&mut self.data, 0, self.size, elem_size)
    }

    pub fn set_end(&mut self, it: &DoubleArrayIterator<'_>) {
        self.size = it.ptr;
    }

    pub fn append(&mut self, d: &DoubleArray) {
        if self.data.len() < self.size + d.size {
            self.data.resize(self.size + d.size, 0);
        }
        self.data[self.size..self.size + d.size].copy_from_slice(&d.data[..d.size]);
        self.size += d.size;
    }

    pub fn offset(&self, it: &DoubleArrayIterator<'_>) -> u32 {
        it.ptr as u32
    }

    pub fn get_u32_at(&self, i: u32) -> u32 {
        let i = i as usize;
        u32::from_ne_bytes(self.data[i..i + 4].try_into().unwrap())
    }

    pub fn data(&self) -> &[u8] {
        &self.data[..self.size]
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleArrayIterator<'a> {
    data: &'a [u8],
    ptr: usize,
    end: usize,
    elem_size: usize,
}

impl<'a> DoubleArrayIterator<'a> {
    fn new(data: &'a [u8], ptr: usize, end: usize, elem_size: usize) -> Self {
        let mut out = Self {
            data,
            ptr,
            end,
            elem_size,
        };
        out.skip_del();
        out
    }

    pub fn count(&self) -> u32 {
        u32::from_ne_bytes(self.data[self.ptr..self.ptr + 4].try_into().unwrap())
    }

    pub fn range_bytes(&self, elem_size: usize) -> &'a [u8] {
        let begin = self.ptr + 4;
        let end = begin + self.count() as usize * elem_size;
        &self.data[begin..end]
    }

    pub fn increment(&mut self, elem_size: usize) -> &mut Self {
        self.elem_size = elem_size;
        self.next(elem_size);
        self.skip_del();
        self
    }

    pub fn next(&mut self, elem_size: usize) {
        self.ptr += self.count() as usize * elem_size + 4;
    }

    pub fn good(&self) -> bool {
        self.ptr < self.end
    }

    pub fn distance(&self, x: &Self) -> isize {
        self.ptr as isize - x.ptr as isize
    }

    fn skip_del(&mut self) {
        while self.ptr < self.end && self.count() == 0 {
            let n = u32::from_ne_bytes(self.data[self.ptr + 4..self.ptr + 8].try_into().unwrap());
            self.ptr += n as usize * self.elem_size + 4;
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct DoubleArrayIteratorMut<'a> {
    data: &'a mut [u8],
    ptr: usize,
    end: usize,
    elem_size: usize,
}

impl<'a> DoubleArrayIteratorMut<'a> {
    fn new(data: &'a mut [u8], ptr: usize, end: usize, elem_size: usize) -> Self {
        let mut out = Self {
            data,
            ptr,
            end,
            elem_size,
        };
        out.skip_del();
        out
    }

    pub fn count(&self) -> u32 {
        u32::from_ne_bytes(self.data[self.ptr..self.ptr + 4].try_into().unwrap())
    }

    pub fn range_bytes(&self, elem_size: usize) -> &[u8] {
        let begin = self.ptr + 4;
        let end = begin + self.count() as usize * elem_size;
        &self.data[begin..end]
    }

    pub fn increment(&mut self, elem_size: usize) -> &mut Self {
        self.elem_size = elem_size;
        self.next(elem_size);
        self.skip_del();
        self
    }

    pub fn next(&mut self, elem_size: usize) {
        self.ptr += self.count() as usize * elem_size + 4;
    }

    pub fn good(&self) -> bool {
        self.ptr < self.end
    }

    pub fn erase(&mut self, elem_size: usize) {
        self.elem_size = elem_size;
        let n = self.count();
        let bytes = n.to_ne_bytes();
        self.data[self.ptr + 4..self.ptr + 8].copy_from_slice(&bytes);
        let zero = 0u32.to_ne_bytes();
        self.data[self.ptr..self.ptr + 4].copy_from_slice(&zero);
        self.ptr += n as usize * elem_size + 4;
    }

    pub fn distance(&self, x: &Self) -> isize {
        self.ptr as isize - x.ptr as isize
    }

    fn skip_del(&mut self) {
        while self.ptr < self.end && self.count() == 0 {
            let n = u32::from_ne_bytes(self.data[self.ptr + 4..self.ptr + 8].try_into().unwrap());
            self.ptr += n as usize * self.elem_size + 4;
        }
    }
}

impl<T> std::ops::Index<usize> for MemBuffer<T> {
    type Output = T;

    fn index(&self, i: usize) -> &Self::Output {
        &self.data[i]
    }
}

impl<T> std::ops::IndexMut<usize> for MemBuffer<T> {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.data[i]
    }
}

impl Default for BitVector {
    fn default() -> Self {
        Self::new()
    }
}

impl BitVector {
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            size: 0,
        }
    }

    pub fn with_size(n: usize) -> Self {
        Self {
            data: vec![0; (n + 63) / 64],
            size: n as u64,
        }
    }

    pub fn set(&mut self, i: usize) {
        self.data[i >> 6] |= 1u64 << (i & 63);
    }

    pub fn get(&self, i: usize) -> bool {
        (self.data[i >> 6] & (1u64 << (i & 63))) != 0
    }

    pub fn or_assign(&mut self, v: &BitVector) {
        for i in 0..self.data.len() {
            self.data[i] |= v.data[i];
        }
    }

    pub fn reset(&mut self) {
        self.data.fill(0);
    }

    pub fn one_count(&self) -> usize {
        self.data.iter().map(|x| x.count_ones() as usize).sum()
    }

    pub fn empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn size(&self) -> u64 {
        self.size
    }

    pub fn negative_list(&self) -> Vec<u64> {
        let mut v = Vec::new();
        for i in 0..self.size {
            if !self.get(i as usize) {
                v.push(i);
            }
        }
        v
    }
}

impl std::ops::BitOrAssign<&BitVector> for BitVector {
    fn bitor_assign(&mut self, rhs: &BitVector) {
        self.or_assign(rhs);
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Array<T> {
    data: Vec<T>,
    alloc_size: usize,
}

impl<T> Default for Array<T> {
    fn default() -> Self {
        Self {
            data: Vec::new(),
            alloc_size: 0,
        }
    }
}

impl<T: Clone> Array<T> {
    pub fn with_alloc_size(alloc_size: i64) -> Self {
        Self {
            data: Vec::with_capacity(alloc_size as usize),
            alloc_size: alloc_size as usize,
        }
    }

    pub fn data(&self) -> &[T] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut [T] {
        &mut self.data
    }

    pub fn begin(&self) -> &[T] {
        &self.data
    }

    pub fn end(&self) -> &[T] {
        &self.data[self.data.len()..]
    }

    pub fn begin_mut(&mut self) -> &mut [T] {
        &mut self.data
    }

    pub fn end_mut(&mut self) -> &mut [T] {
        let len = self.data.len();
        &mut self.data[len..]
    }

    pub fn size(&self) -> i64 {
        self.data.len() as i64
    }

    pub fn assign_value(&mut self, x: &T) {
        self.data.clear();
        self.data.push(x.clone());
    }

    pub fn assign(&mut self, begin: &[T]) {
        self.data.clear();
        self.data.extend_from_slice(begin);
        assert!(self.data.len() <= self.alloc_size || self.alloc_size == 0);
    }

    pub fn assign_reversed(&mut self, begin: &[T]) {
        self.data.clear();
        self.data.extend(begin.iter().rev().cloned());
        assert!(self.data.len() <= self.alloc_size || self.alloc_size == 0);
    }

    pub fn push_back(&mut self, begin: &[T]) {
        self.data.extend_from_slice(begin);
        assert!(self.data.len() <= self.alloc_size || self.alloc_size == 0);
    }

    pub fn push_back_reversed(&mut self, begin: &[T]) {
        self.data.extend(begin.iter().rev().cloned());
        assert!(self.data.len() <= self.alloc_size || self.alloc_size == 0);
    }

    pub fn push_back_value(&mut self, n: i64, value: &T) {
        self.data
            .extend(std::iter::repeat(value.clone()).take(n as usize));
        assert!(self.data.len() <= self.alloc_size || self.alloc_size == 0);
    }

    pub fn get(&self, i: i32) -> T
    where
        T: Copy,
    {
        self.data[i as usize]
    }
}

impl<T> std::ops::Index<usize> for Array<T> {
    type Output = T;

    fn index(&self, i: usize) -> &Self::Output {
        &self.data[i]
    }
}

impl<T> std::ops::Index<i32> for Array<T> {
    type Output = T;

    fn index(&self, i: i32) -> &Self::Output {
        &self.data[i as usize]
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DisjointSet<Int = usize> {
    nodes: Vec<DisjointSetNode<Int>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct DisjointSetNode<Int> {
    parent: Int,
    size: Int,
}

impl<Int> DisjointSet<Int>
where
    Int: Copy + Ord + From<u8> + TryFrom<usize> + TryInto<usize> + std::hash::Hash + Eq,
    <Int as TryFrom<usize>>::Error: std::fmt::Debug,
    <Int as TryInto<usize>>::Error: std::fmt::Debug,
{
    pub fn new(size: Int) -> Self {
        let size_usize: usize = size.try_into().unwrap();
        let mut nodes = Vec::with_capacity(size_usize);
        for i in 0..size_usize {
            nodes.push(DisjointSetNode {
                parent: Int::try_from(i).unwrap(),
                size: Int::from(1),
            });
        }
        Self { nodes }
    }

    pub fn find(&mut self, i: Int) -> Int {
        let idx: usize = i.try_into().unwrap();
        if self.nodes[idx].parent != i {
            self.nodes[idx].parent = self.find(self.nodes[idx].parent);
            self.nodes[idx].parent
        } else {
            i
        }
    }

    pub fn merge(&mut self, i: Int, j: Int) {
        let mut i = self.find(i);
        let mut j = self.find(j);
        if i == j {
            return;
        }
        let mut i_idx: usize = i.try_into().unwrap();
        let mut j_idx: usize = j.try_into().unwrap();
        if self.nodes[i_idx].size < self.nodes[j_idx].size {
            std::mem::swap(&mut i, &mut j);
            std::mem::swap(&mut i_idx, &mut j_idx);
        }
        self.nodes[j_idx].parent = i;
        let size_i: usize = self.nodes[i_idx].size.try_into().unwrap();
        let size_j: usize = self.nodes[j_idx].size.try_into().unwrap();
        self.nodes[i_idx].size = Int::try_from(size_i + size_j).unwrap();
    }

    pub fn sets(&mut self, _threads: i32) -> FlatArray<Int> {
        let mut v = Vec::with_capacity(self.nodes.len());
        for i in 0..self.nodes.len() {
            let item = Int::try_from(i).unwrap();
            v.push((self.find(item), item));
        }
        make_flat_array(&mut v).0
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RangePartition<const N: usize, T> {
    begin: [i32; N_PLUS_ONE_PLACEHOLDER],
    mask: [[T; N]; N],
    count: i32,
}

const N_PLUS_ONE_PLACEHOLDER: usize = 65;

impl<const N: usize, T> RangePartition<N, T> {
    pub fn begin(&self, i: i32) -> i32 {
        self.begin[i as usize]
    }

    pub fn end(&self, i: i32) -> i32 {
        self.begin[i as usize + 1]
    }

    pub fn count(&self) -> i32 {
        self.count
    }

    pub fn mask(&self, i: i32) -> &[T; N] {
        &self.mask[i as usize]
    }
}

pub trait RangePartitionValue: Copy + Ord {
    fn min_value() -> Self;
    fn zero() -> Self;
}

macro_rules! impl_range_partition_value {
    ($($t:ty),* $(,)?) => {
        $(
            impl RangePartitionValue for $t {
                fn min_value() -> Self {
                    <$t>::MIN
                }

                fn zero() -> Self {
                    0
                }
            }
        )*
    };
}

impl_range_partition_value!(i8, i16, i32, i64, isize);

impl<const N: usize, T> RangePartition<N, T>
where
    T: RangePartitionValue,
{
    pub fn new(begin: &[i32], count: i32, end: i32) -> Self {
        assert!(N + 1 <= N_PLUS_ONE_PLACEHOLDER);
        assert!(count > 0);
        let mut b = [(0i32, 0usize); N];
        for i in 0..count as usize {
            b[i] = (begin[i], i);
        }
        b[..count as usize].sort_unstable();

        let mut mask = [[T::min_value(); N]; N];
        let mut begin_out = [0i32; N_PLUS_ONE_PLACEHOLDER];
        let mut count_out = 1i32;
        mask[0][b[0].1] = T::zero();
        begin_out[0] = b[0].0;
        for item in b.iter().take(count as usize).skip(1) {
            if begin_out[count_out as usize - 1] < item.0 {
                begin_out[count_out as usize] = item.0;
                mask[count_out as usize] = mask[count_out as usize - 1];
                mask[count_out as usize][item.1] = T::zero();
                count_out += 1;
            } else {
                mask[count_out as usize - 1][item.1] = T::zero();
            }
        }
        begin_out[count_out as usize] = end;
        Self {
            begin: begin_out,
            mask,
            count: count_out,
        }
    }

    pub fn from_cpp(begin: &[i32], count: i32, end: i32) -> Self {
        Self::new(begin, count, end)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FlatArray<T, I = u64> {
    data: Vec<T>,
    limits: Vec<I>,
}

impl<T> Default for FlatArray<T, u64> {
    fn default() -> Self {
        Self::new()
    }
}

pub trait FlatIndex:
    Copy + Ord + Default + From<u8> + TryFrom<usize> + TryInto<usize> + std::ops::Add<Output = Self>
{
}

impl FlatIndex for usize {}
impl FlatIndex for u64 {}
impl FlatIndex for u32 {}

impl<T, I> FlatArray<T, I>
where
    I: FlatIndex,
    <I as TryInto<usize>>::Error: std::fmt::Debug,
    <I as TryFrom<usize>>::Error: std::fmt::Debug,
{
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            limits: vec![I::from(0)],
        }
    }

    pub fn from_sizes(mut sizes: Vec<I>) -> Self
    where
        T: Default + Clone,
    {
        let data_size = Self::sum(&mut sizes);
        Self {
            data: vec![T::default(); data_size.try_into().unwrap()],
            limits: sizes,
        }
    }

    pub fn from_limits_data(limits: Vec<I>, data: Vec<T>) -> Self {
        Self { data, limits }
    }

    pub fn push_item(&mut self, x: T) {
        self.data.push(x);
        let last = self.limits.last_mut().unwrap();
        let n: usize = (*last).try_into().unwrap();
        *last = I::try_from(n + 1).unwrap();
    }

    pub fn push_back(&mut self, values: &[T])
    where
        T: Clone,
    {
        self.data.extend_from_slice(values);
        let last: usize = self.limits.last().copied().unwrap().try_into().unwrap();
        self.limits.push(I::try_from(last + values.len()).unwrap());
    }

    pub fn next(&mut self) {
        self.limits.push(*self.limits.last().unwrap());
    }

    pub fn pop_back(&mut self) {
        self.limits.pop();
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.limits.clear();
        self.limits.push(I::from(0));
    }

    pub fn shrink_to_fit(&mut self) {
        self.data.shrink_to_fit();
        self.limits.shrink_to_fit();
    }

    pub fn size(&self) -> I {
        I::try_from(self.limits.len() - 1).unwrap()
    }

    pub fn data_size(&self) -> I {
        I::try_from(self.data.len()).unwrap()
    }

    pub fn begin(&self, i: I) -> usize {
        let i: usize = i.try_into().unwrap();
        self.limits[i].try_into().unwrap()
    }

    pub fn end(&self, i: I) -> usize {
        let i: usize = i.try_into().unwrap();
        self.limits[i + 1].try_into().unwrap()
    }

    pub fn cbegin(&self, i: I) -> usize {
        self.begin(i)
    }

    pub fn cend(&self, i: I) -> usize {
        self.end(i)
    }

    pub fn range(&self, i: I) -> &[T] {
        &self.data[self.begin(i)..self.end(i)]
    }

    pub fn range_mut(&mut self, i: I) -> &mut [T] {
        let begin = self.begin(i);
        let end = self.end(i);
        &mut self.data[begin..end]
    }

    pub fn count(&self, i: I) -> I {
        let i: usize = i.try_into().unwrap();
        let begin: usize = self.limits[i].try_into().unwrap();
        let end: usize = self.limits[i + 1].try_into().unwrap();
        I::try_from(end - begin).unwrap()
    }

    pub fn reserve(&mut self, size: I, data_size: I) {
        self.data.reserve(data_size.try_into().unwrap());
        self.limits.reserve(size.try_into().unwrap() + 1);
    }

    pub fn max_count(&self) -> I {
        let mut m = 0usize;
        let n: usize = self.size().try_into().unwrap();
        for i in 0..n {
            let c: usize = self.count(I::try_from(i).unwrap()).try_into().unwrap();
            m = m.max(c);
        }
        I::try_from(m).unwrap()
    }

    pub fn global_cbegin(&self) -> &[T] {
        &self.data
    }

    pub fn global_cend(&self) -> &[T] {
        &self.data[self.data.len()..]
    }

    pub fn data(&self) -> &[T] {
        &self.data
    }

    pub fn data_mut(&mut self) -> &mut [T] {
        &mut self.data
    }

    pub fn limits(&self) -> &[I] {
        &self.limits
    }

    pub fn limits_mut(&mut self) -> &mut [I] {
        &mut self.limits
    }

    fn sum(sizes: &mut [I]) -> I {
        assert!(sizes.last().copied().unwrap_or_else(|| I::from(0)) == I::from(0));
        let mut p = 0usize;
        for size in sizes {
            let n: usize = (*size).try_into().unwrap();
            *size = I::try_from(p).unwrap();
            p += n;
        }
        I::try_from(p).unwrap()
    }
}

impl<T, I> std::ops::Index<I> for FlatArray<T, I>
where
    I: FlatIndex,
    <I as TryInto<usize>>::Error: std::fmt::Debug,
    <I as TryFrom<usize>>::Error: std::fmt::Debug,
{
    type Output = [T];

    fn index(&self, i: I) -> &Self::Output {
        self.range(i)
    }
}

pub fn make_flat_array<K, V>(items: &mut [(K, V)]) -> (FlatArray<V>, Vec<K>)
where
    K: Copy + Ord + std::hash::Hash + Eq,
    V: Copy + Ord,
{
    items.sort_unstable_by(|a, b| a.cmp(b));
    let mut r = (FlatArray::new(), Vec::new());
    let mut i = 0usize;
    while i < items.len() {
        let key = items[i].0;
        r.1.push(key);
        let begin = i;
        while i < items.len() && items[i].0 == key {
            i += 1;
        }
        let values: Vec<V> = items[begin..i].iter().map(|x| x.1).collect();
        r.0.push_back(&values);
    }
    r
}

pub fn make_flat_array_dense<K, V>(items: &mut [(K, V)]) -> FlatArray<V>
where
    K: Copy + Ord + TryInto<usize>,
    <K as TryInto<usize>>::Error: std::fmt::Debug,
    V: Copy + Ord,
{
    items.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    let mut r = FlatArray::new();
    let mut k = 0usize;
    let mut i = 0usize;
    while i < items.len() {
        let key: usize = items[i].0.try_into().unwrap();
        while k != key {
            r.next();
            k += 1;
        }
        let begin = i;
        while i < items.len() && items[i].0.try_into().unwrap() == key {
            i += 1;
        }
        let values: Vec<V> = items[begin..i].iter().map(|x| x.1).collect();
        r.push_back(&values);
        k += 1;
    }
    r
}

pub fn make_flat_array_dense_data<T, K, GetKey>(
    mut data: Vec<T>,
    key_end: K,
    mut get_key: GetKey,
) -> FlatArray<T>
where
    T: Ord,
    K: Copy + Ord + TryInto<usize>,
    <K as TryInto<usize>>::Error: std::fmt::Debug,
    GetKey: FnMut(&T) -> K,
{
    data.sort_unstable();
    let key_end: usize = key_end.try_into().unwrap();
    let mut limits = Vec::with_capacity(key_end + 1);
    limits.push(0u64);
    let mut k = 0usize;
    let mut i = 0usize;
    while i < data.len() {
        let key: usize = get_key(&data[i]).try_into().unwrap();
        while k != key {
            limits.push(*limits.last().unwrap());
            k += 1;
        }
        while i < data.len() && get_key(&data[i]).try_into().unwrap() == key {
            i += 1;
        }
        limits.push(i as u64);
        k += 1;
    }
    while k < key_end {
        limits.push(*limits.last().unwrap());
        k += 1;
    }
    FlatArray::from_limits_data(limits, data)
}

#[derive(Debug, Clone)]
pub struct SparseFlatArray<Key, T> {
    map: HashMap<Key, i64>,
    data: FlatArray<T>,
}

impl<Key, T> SparseFlatArray<Key, T>
where
    Key: Copy + Ord + std::hash::Hash + Eq,
    T: Copy + Ord,
{
    pub fn new(items: &mut [(Key, T)]) -> Self {
        items.sort_unstable();
        let mut map = HashMap::new();
        let mut data = FlatArray::new();
        let mut i = 0usize;
        while i < items.len() {
            let key = items[i].0;
            map.insert(key, data.size() as i64);
            let begin = i;
            while i < items.len() && items[i].0 == key {
                i += 1;
            }
            let values: Vec<T> = items[begin..i].iter().map(|x| x.1).collect();
            data.push_back(&values);
        }
        Self { map, data }
    }

    pub fn empty(&self) -> bool {
        self.map.is_empty()
    }

    pub fn size(&self) -> i64 {
        self.map.len() as i64
    }

    pub fn data_size(&self) -> i64 {
        self.data.data_size() as i64
    }

    pub fn range(&self, key: Key) -> &[T] {
        self.data.range(*self.map.get(&key).unwrap() as u64)
    }

    pub fn cbegin(&self, key: Key) -> usize {
        self.data.cbegin(*self.map.get(&key).unwrap() as u64)
    }

    pub fn cend(&self, key: Key) -> usize {
        self.data.cend(*self.map.get(&key).unwrap() as u64)
    }

    pub fn iter(&self) -> impl Iterator<Item = (Key, &[T])> {
        self.map
            .iter()
            .map(|(&key, &idx)| (key, self.data.range(idx as u64)))
    }

    pub fn flat_iter(&self) -> impl Iterator<Item = (Key, T)> + '_ {
        self.map.iter().flat_map(|(&key, &idx)| {
            self.data
                .range(idx as u64)
                .iter()
                .copied()
                .map(move |value| (key, value))
        })
    }
}

pub trait AllocSize {
    fn alloc_size(&self) -> usize;
}

#[derive(Debug, Clone)]
pub struct ReorderQueue<T> {
    begin: usize,
    next: usize,
    size: usize,
    max_size: usize,
    backlog: BTreeMap<usize, T>,
}

impl<T> ReorderQueue<T> {
    pub fn new(begin: usize) -> Self {
        Self {
            begin,
            next: begin,
            size: 0,
            max_size: 0,
            backlog: BTreeMap::new(),
        }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn max_size(&self) -> usize {
        self.max_size
    }

    pub fn next(&self) -> usize {
        self.next
    }

    pub fn begin(&self) -> usize {
        self.begin
    }
}

impl<T> ReorderQueue<T>
where
    T: AllocSize,
{
    pub fn push<F>(&mut self, n: usize, value: Option<T>, mut f: F)
    where
        F: FnMut(T),
    {
        if n != self.next {
            if let Some(value) = value {
                self.size += value.alloc_size();
                self.max_size = self.max_size.max(self.size);
                self.backlog.insert(n, value);
            }
        } else {
            self.flush(value, &mut f);
        }
    }

    fn flush<F>(&mut self, value: Option<T>, f: &mut F)
    where
        F: FnMut(T),
    {
        let mut n = self.next + 1;
        if let Some(value) = value {
            f(value);
        }
        while let Some(value) = self.backlog.remove(&n) {
            let size = value.alloc_size();
            f(value);
            self.size -= size;
            n += 1;
        }
        self.next = n;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_vector() {
        let mut v = BitVector::with_size(70);
        assert!(!v.empty());
        assert_eq!(v.size(), 70);
        v.set(0);
        v.set(65);
        assert!(v.get(0));
        assert!(v.get(65));
        assert_eq!(v.one_count(), 2);
        assert_eq!(v.negative_list()[0], 1);

        let mut w = BitVector::with_size(70);
        w.set(1);
        v.or_assign(&w);
        assert!(v.get(1));
        v |= &w;
        assert!(v.get(1));
        v.reset();
        assert_eq!(v.one_count(), 0);
        assert!(BitVector::new().empty());
    }

    #[test]
    fn test_hash_set() {
        let mut set = HashSet::<Modulo, Identity>::new(8);
        assert_eq!(set.size(), 8);
        assert_eq!(set.load(), 0);
        assert!(!set.contains(0x1234));
        set.insert(0x1234);
        assert!(set.contains(0x1234));
        assert_eq!(set.load(), 1);
        set.insert(0x1234);
        assert_eq!(set.load(), 1);
        assert_eq!(set.table(), set.data());
        set.table_mut()[0] = set.table()[0];
        set.finish();
        assert_eq!(set.data().len(), HashSet::<Modulo, Identity>::PADDING + 8);

        let mut set2 = HashSet::<Modulo2, Identity>::new(8);
        set2.insert(0xff00);
        assert!(set2.contains(0xff00));

        let mut set3 = HashSet::<Modulo, crate::util::hash::MurmurHash>::new(8);
        set3.insert(0x1234);
        assert!(set3.contains(0x1234));
    }

    #[test]
    fn test_hash_table() {
        let mut table = HashTable::<usize, usize, IdentityHash, Modulo>::new(4);
        *table.get_mut(1).unwrap() = 10;
        *table.get_mut(5).unwrap() = 50;
        assert_eq!(table.count(), 2);
        assert_eq!(*table.find(1).unwrap().unwrap(), 10);
        assert_eq!(*table.find(5).unwrap().unwrap(), 50);
        assert!(table.find(9).unwrap().is_none());
        table.insert(9).unwrap().value = 90;
        assert_eq!(table.find_entry(9).unwrap().unwrap().value, 90);
        assert_eq!(table.size(), 4);

        let mut table = HashTable::<u64, u32, crate::util::hash::MurmurHash, Modulo>::new(8);
        table.insert(42).unwrap().value = 7;
        assert_eq!(table.find(42).unwrap().copied(), Some(7));
    }

    #[test]
    fn test_queue() {
        let mut q = Queue::new(3, 1, 2, -1);
        assert_eq!(q.capacity(), 4);
        assert_eq!(q.producer_count(), 1);
        assert!(q.empty());
        q.enqueue(7);
        q.enqueue(8);
        assert_eq!(q.approx_size(), 2);
        assert_eq!(q.wait_and_dequeue(), Some(7));
        q.close();
        assert_eq!(q.wait_and_dequeue(), Some(8));
        assert_eq!(q.wait_and_dequeue(), None);

        let mut q = Queue::new(2, 2, 1, 0);
        q.enqueue(1);
        q.enqueue(0);
        assert_eq!(q.wait_and_dequeue(), Some(1));
        assert_eq!(q.wait_and_dequeue(), None);
    }

    #[test]
    fn test_ptr_vector() {
        let mut v = PtrVector::new();
        v.push(1);
        v.push(2);
        v.push(3);
        assert_eq!(v[1], 2);
        *v.get_ptr(1) = Box::new(9);
        assert_eq!(v[1], 9);
        *v.get(1) = Box::new(8);
        assert_eq!(v[1], 8);
        *v.back() = 4;
        assert_eq!(v[2], 4);
        v.erase(0, 1);
        assert_eq!(v.len(), 2);
        assert_eq!(v[0], 8);
        v.clear();
        assert!(v.is_empty());

        let sized = PtrVector::<i32>::with_size(2);
        assert_eq!(sized.len(), 2);
        assert_eq!(sized[0], 0);
    }

    #[test]
    fn test_deque_and_async_writer() {
        let mut d = Deque::<i32, 2>::new();
        d.push_back(1);
        d.push_back_slice(&[2, 3]);
        d.push_back_iter(&[4, 5, 6, 7]);
        assert_eq!(d.size(), 7);
        let mut it = d.begin();
        assert_eq!(*it.get(), 1);
        it.increment();
        assert_eq!(*it.get(), 2);
        assert_eq!(*it.get_offset(2), 4);
        let end = d.end();
        assert_eq!(end.distance(&it), 6);
        let mut moved = Vec::new();
        d.move_to(&mut moved);
        assert_eq!(moved, vec![1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(d.size(), 0);

        let mut d = Deque::<i32, 3>::new();
        {
            let mut w = AsyncWriter::new(&mut d);
            w.write(10);
            w.write(11);
        }
        let mut moved = Vec::new();
        d.move_to(&mut moved);
        assert_eq!(moved, vec![10, 11]);
    }

    #[test]
    fn test_array() {
        let mut a = Array::with_alloc_size(10);
        a.assign_value(&3);
        assert_eq!(a.data(), &[3]);
        a.assign(&[1, 2, 3]);
        assert_eq!(a.size(), 3);
        assert_eq!(a.get(1), 2);
        assert_eq!(a[1usize], 2);
        assert_eq!(a[2i32], 3);
        a.begin_mut()[0] = 9;
        assert_eq!(a.data(), &[9, 2, 3]);
        assert_eq!(a.end_mut(), &mut []);
        a.assign_reversed(&[1, 2, 3]);
        assert_eq!(a.data(), &[3, 2, 1]);
        a.push_back(&[4, 5]);
        assert_eq!(a.data(), &[3, 2, 1, 4, 5]);
        a.push_back_reversed(&[6, 7]);
        assert_eq!(a.data(), &[3, 2, 1, 4, 5, 7, 6]);
        a.push_back_value(2, &9);
        assert_eq!(a.data(), &[3, 2, 1, 4, 5, 7, 6, 9, 9]);
    }

    #[test]
    fn test_disjoint_set() {
        let mut d = DisjointSet::new(5);
        d.merge(0, 1);
        d.merge(3, 4);
        assert_eq!(d.find(0), d.find(1));
        assert_ne!(d.find(0), d.find(2));
        d.merge(1, 4);
        let sets = d.sets(1);
        assert_eq!(sets.size(), 2);
        assert_eq!(sets.range(0), &[0, 1, 3, 4]);
        assert_eq!(sets.range(1), &[2]);

        let mut d = DisjointSet::<u32>::new(4);
        d.merge(0, 2);
        assert_eq!(d.find(0), d.find(2));
        let sets = d.sets(1);
        assert_eq!(sets.global_cbegin(), &[0, 2, 1, 3]);
    }

    #[test]
    fn test_mem_buffer() {
        let mut b = MemBuffer::<i32>::with_size(3);
        assert_eq!(b.size(), 3);
        *b.get_mut(1) = 9;
        assert_eq!(*b.get(1), 9);
        assert_eq!(b.begin(), &[0, 9, 0]);
        assert_eq!(b.data(), &[0, 9, 0]);
        b[2] = 8;
        assert_eq!(b[2], 8);
        b.begin_mut()[0] = 7;
        assert_eq!(b.data_mut(), &mut [7, 9, 8]);
        assert_eq!(b.end(), &[]);
        assert_eq!(b.end_mut(), &mut []);
        b.resize(5);
        assert_eq!(b.size(), 5);
        b.resize(2);
        assert_eq!(b.begin(), &[0, 0]);
    }

    #[test]
    fn test_double_array() {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&2u32.to_ne_bytes());
        bytes.extend_from_slice(&10u32.to_ne_bytes());
        bytes.extend_from_slice(&11u32.to_ne_bytes());
        bytes.extend_from_slice(&1u32.to_ne_bytes());
        bytes.extend_from_slice(&20u32.to_ne_bytes());
        let mut d = DoubleArray::from_bytes(bytes, 20);

        let mut it = d.begin();
        assert!(it.good());
        assert_eq!(it.count(), 2);
        assert_eq!(
            it.range_bytes(4),
            &[10u32.to_ne_bytes(), 11u32.to_ne_bytes()].concat()
        );
        it.increment(4);
        assert_eq!(it.count(), 1);
        assert_eq!(d.offset(&it), 12);
        assert_eq!(d.get_u32_at(16), 20);

        let mut it = d.begin_mut();
        it.erase(4);
        assert_eq!(it.count(), 1);
        drop(it);
        let it = d.begin();
        assert_eq!(it.count(), 1);
        assert_eq!(d.offset(&it), 12);

        let mut dst = DoubleArray::from_bytes(vec![0; 40], 0);
        dst.append(&d);
        assert_eq!(dst.data(), d.data());

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&1u32.to_ne_bytes());
        bytes.extend_from_slice(&10u64.to_ne_bytes());
        bytes.extend_from_slice(&1u32.to_ne_bytes());
        bytes.extend_from_slice(&20u64.to_ne_bytes());
        let mut d = DoubleArray::from_bytes(bytes, 24);
        let mut it = d.begin_mut_with_elem_size(8);
        it.erase(8);
        drop(it);
        let it = d.begin_with_elem_size(8);
        assert_eq!(it.count(), 1);
        assert_eq!(d.offset(&it), 12);
        assert_eq!(it.range_bytes(8), &20u64.to_ne_bytes());
    }

    #[test]
    fn test_range_partition() {
        let p = RangePartition::<4, i32>::from_cpp(&[8, 2, 2, 5], 4, 10);
        assert_eq!(p.count(), 3);
        assert_eq!(p.begin(0), 2);
        assert_eq!(p.end(0), 5);
        assert_eq!(p.begin(1), 5);
        assert_eq!(p.end(1), 8);
        assert_eq!(p.begin(2), 8);
        assert_eq!(p.end(2), 10);
        assert_eq!(p.mask(0), &[i32::MIN, 0, 0, i32::MIN]);
        assert_eq!(p.mask(1), &[i32::MIN, 0, 0, 0]);
        assert_eq!(p.mask(2), &[0, 0, 0, 0]);
    }

    #[test]
    fn test_flat_array() {
        let mut a = FlatArray::<i32>::new();
        assert_eq!(a.size(), 0);
        a.push_back(&[1, 2]);
        a.push_back(&[3]);
        a.next();
        a.push_item(4);
        assert_eq!(a.size(), 3);
        assert_eq!(a.data_size(), 4);
        assert_eq!(a.range(0), &[1, 2]);
        assert_eq!(a.range(1), &[3]);
        assert_eq!(a.range(2), &[4]);
        assert_eq!(&a[2u64], &[4]);
        assert_eq!(a.count(0), 2);
        assert_eq!(a.max_count(), 2);
        assert_eq!(a.global_cbegin(), &[1, 2, 3, 4]);
        assert_eq!(a.data(), &[1, 2, 3, 4]);
        assert_eq!(a.limits(), &[0, 2, 3, 4]);
        a.data_mut()[0] = 9;
        assert_eq!(a.range(0), &[9, 2]);
        a.limits_mut()[1] = 1;
        assert_eq!(a.range(0), &[9]);
        assert_eq!(a.range(1), &[2, 3]);

        let b = FlatArray::<i32>::from_sizes(vec![2u64, 1, 0]);
        assert_eq!(b.size(), 2);
        assert_eq!(b.data_size(), 3);

        let mut pairs = vec![(2u64, 20), (0, 1), (2, 21), (1, 10)];
        let (flat, keys) = make_flat_array(&mut pairs);
        assert_eq!(keys, vec![0, 1, 2]);
        assert_eq!(flat.range(0), &[1]);
        assert_eq!(flat.range(2), &[20, 21]);

        let mut dense_pairs = vec![(2u64, 20), (2, 21), (0, 1)];
        let dense = make_flat_array_dense(&mut dense_pairs);
        assert_eq!(dense.size(), 3);
        assert_eq!(dense.range(0), &[1]);
        assert_eq!(dense.range(1), &[] as &[i32]);
        assert_eq!(dense.range(2), &[20, 21]);

        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
        struct Pair {
            key: u64,
            value: i32,
        }
        let dense_data = make_flat_array_dense_data(
            vec![
                Pair { key: 2, value: 20 },
                Pair { key: 0, value: 1 },
                Pair { key: 2, value: 21 },
            ],
            4u64,
            |x| x.key,
        );
        assert_eq!(dense_data.size(), 4);
        assert_eq!(
            dense_data
                .range(0)
                .iter()
                .map(|x| x.value)
                .collect::<Vec<_>>(),
            vec![1]
        );
        assert_eq!(dense_data.range(1), &[] as &[Pair]);
        assert_eq!(
            dense_data
                .range(2)
                .iter()
                .map(|x| x.value)
                .collect::<Vec<_>>(),
            vec![20, 21]
        );
        assert_eq!(dense_data.range(3), &[] as &[Pair]);
    }

    #[test]
    fn test_sparse_flat_array() {
        let mut items = vec![("b", 20), ("a", 1), ("b", 21)];
        let s = SparseFlatArray::new(&mut items);
        assert!(!s.empty());
        assert_eq!(s.size(), 2);
        assert_eq!(s.data_size(), 3);
        assert_eq!(s.range("a"), &[1]);
        assert_eq!(s.range("b"), &[20, 21]);
        assert_eq!(s.cbegin("b"), 1);
        assert_eq!(s.cend("b"), 3);
        let mut flat = s.flat_iter().collect::<Vec<_>>();
        flat.sort();
        assert_eq!(flat, vec![("a", 1), ("b", 20), ("b", 21)]);
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    struct SizedItem {
        value: i32,
        bytes: usize,
    }

    impl AllocSize for SizedItem {
        fn alloc_size(&self) -> usize {
            self.bytes
        }
    }

    #[test]
    fn test_reorder_queue() {
        let mut q = ReorderQueue::new(10);
        let mut out = Vec::new();
        q.push(
            12,
            Some(SizedItem {
                value: 12,
                bytes: 5,
            }),
            |x| out.push(x.value),
        );
        assert_eq!(out, Vec::<i32>::new());
        assert_eq!(q.size(), 5);
        assert_eq!(q.max_size(), 5);
        assert_eq!(q.next(), 10);
        assert_eq!(q.begin(), 10);

        q.push(
            10,
            Some(SizedItem {
                value: 10,
                bytes: 1,
            }),
            |x| out.push(x.value),
        );
        assert_eq!(out, vec![10]);
        assert_eq!(q.next(), 11);

        q.push(
            11,
            Some(SizedItem {
                value: 11,
                bytes: 2,
            }),
            |x| out.push(x.value),
        );
        assert_eq!(out, vec![10, 11, 12]);
        assert_eq!(q.next(), 13);
        assert_eq!(q.size(), 0);
    }
}
