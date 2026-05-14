use std::alloc::{alloc, dealloc, handle_alloc_error, Layout};
use std::marker::PhantomData;
use std::mem;
use std::ptr::NonNull;

#[derive(Debug)]
pub struct AlignedAllocation {
    ptr: NonNull<u8>,
    layout: Layout,
    len: usize,
}

impl AlignedAllocation {
    pub fn data(&self) -> *mut u8 {
        self.ptr.as_ptr()
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn align(&self) -> usize {
        self.layout.align()
    }

    pub fn as_slice(&self) -> &[u8] {
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }

    pub fn as_mut_slice(&mut self) -> &mut [u8] {
        unsafe { std::slice::from_raw_parts_mut(self.ptr.as_ptr(), self.len) }
    }
}

impl Drop for AlignedAllocation {
    fn drop(&mut self) {
        unsafe {
            dealloc(self.ptr.as_ptr(), self.layout);
        }
    }
}

pub fn aligned_malloc(n: usize, align: usize) -> AlignedAllocation {
    let layout = Layout::from_size_align(n.max(1), align).expect("invalid alignment");
    let ptr = unsafe { alloc(layout) };
    let ptr = NonNull::new(ptr).unwrap_or_else(|| handle_alloc_error(layout));
    AlignedAllocation {
        ptr,
        layout,
        len: n,
    }
}

pub fn aligned_free(p: AlignedAllocation) {
    drop(p);
}

#[derive(Debug, Clone, Copy, Default)]
pub struct AlignmentAllocator<T, const N: usize = 16> {
    _marker: PhantomData<T>,
}

impl<T, const N: usize> AlignmentAllocator<T, N> {
    pub fn new() -> Self {
        Self {
            _marker: PhantomData,
        }
    }

    pub fn adress(&self, r: &T) -> *const T {
        r as *const T
    }

    pub fn adress_mut(&self, r: &mut T) -> *mut T {
        r as *mut T
    }

    pub fn address(&self, r: &T) -> *const T {
        self.adress(r)
    }

    pub fn address_mut(&self, r: &mut T) -> *mut T {
        self.adress_mut(r)
    }

    pub fn allocate(&self, n: usize) -> *mut T {
        let layout = Layout::from_size_align((n * mem::size_of::<T>()).max(1), N)
            .expect("invalid alignment");
        let ptr = unsafe { alloc(layout) };
        NonNull::new(ptr)
            .unwrap_or_else(|| handle_alloc_error(layout))
            .as_ptr() as *mut T
    }

    pub unsafe fn deallocate(&self, p: *mut T, n: usize) {
        let layout = Layout::from_size_align((n * mem::size_of::<T>()).max(1), N)
            .expect("invalid alignment");
        unsafe {
            dealloc(p as *mut u8, layout);
        }
    }

    pub unsafe fn construct(&self, p: *mut T, wert: T) {
        unsafe {
            p.write(wert);
        }
    }

    pub unsafe fn destroy(&self, p: *mut T) {
        unsafe {
            p.drop_in_place();
        }
    }

    pub fn max_size(&self) -> usize {
        usize::MAX / mem::size_of::<T>()
    }
}

impl<T, const N: usize> PartialEq for AlignmentAllocator<T, N> {
    fn eq(&self, _other: &Self) -> bool {
        true
    }
}

impl<T, const N: usize> Eq for AlignmentAllocator<T, N> {}

pub mod pmr {
    use std::collections::LinkedList;

    #[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
    pub struct MonotonicBufferResource;

    pub type List<T> = LinkedList<T>;
    pub type Vector<T> = Vec<T>;
    pub type String = std::string::String;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aligned_malloc_free() {
        let mut p = aligned_malloc(64, 32);
        assert_eq!((p.data() as usize) % 32, 0);
        assert_eq!(p.len(), 64);
        p.as_mut_slice()[0] = 7;
        assert_eq!(p.as_slice()[0], 7);
        aligned_free(p);
    }

    #[test]
    fn test_alignment_allocator() {
        let allocator = AlignmentAllocator::<u64, 32>::new();
        let p = allocator.allocate(2);
        assert_eq!((p as usize) % 32, 0);
        unsafe {
            allocator.construct(p, 11);
            allocator.construct(p.add(1), 13);
            assert_eq!(*p, 11);
            assert_eq!(*p.add(1), 13);
            allocator.destroy(p);
            allocator.destroy(p.add(1));
            allocator.deallocate(p, 2);
        }
        assert_eq!(allocator.max_size(), usize::MAX / mem::size_of::<u64>());
        assert_eq!(allocator, AlignmentAllocator::<u64, 32>::new());
    }

    #[test]
    fn test_pmr_aliases() {
        let _resource = pmr::MonotonicBufferResource;
        let mut vector: pmr::Vector<i32> = Vec::new();
        vector.push(3);
        let mut list: pmr::List<i32> = std::collections::LinkedList::new();
        list.push_back(vector[0]);
        let string: pmr::String = "x".to_string();
        assert_eq!(list.front(), Some(&3));
        assert_eq!(string, "x");
    }
}
