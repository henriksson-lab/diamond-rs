//! Pinned GNU libstdc++ `std::sort` comparison and tie movement.
//!
//! Extracted from the verified SketchIterator translation. C++ does not
//! standardize comparator-equivalent order. We use the pinned GNU reference
//! algorithm on every Rust target, including Apple; this does not claim the
//! permutation of an independently built Apple libc++ Source executable.

use std::cmp::Ordering;

const THRESHOLD: usize = 16;

pub fn sort_by<T: Copy, F: Fn(&T, &T) -> Ordering>(values: &mut [T], compare: F) {
    if values.len() < 2 {
        return;
    }
    let depth = 2 * (usize::BITS as usize - 1 - values.len().leading_zeros() as usize);
    sort_by_with_depth(values, depth, compare);
}

pub(crate) fn sort_by_with_depth<T: Copy, F: Fn(&T, &T) -> Ordering>(
    values: &mut [T],
    depth: usize,
    compare: F,
) {
    let len = values.len();
    if len < 2 {
        return;
    }
    let mut sorter = Sorter {
        values,
        less: |a: T, b: T| compare(&a, &b) == Ordering::Less,
    };
    sorter.introsort(0, len, depth);
    if len > THRESHOLD {
        sorter.insertion_sort(0, THRESHOLD);
        sorter.unguarded_insertion_sort(THRESHOLD, len);
    } else {
        sorter.insertion_sort(0, len);
    }
}

struct Sorter<'a, T, F> {
    values: &'a mut [T],
    less: F,
}
impl<T: Copy, F: Fn(T, T) -> bool> Sorter<'_, T, F> {
    fn move_median_to_first(&mut self, result: usize, a: usize, b: usize, c: usize) {
        if (self.less)(self.values[a], self.values[b]) {
            if (self.less)(self.values[b], self.values[c]) {
                self.values.swap(result, b);
            } else if (self.less)(self.values[a], self.values[c]) {
                self.values.swap(result, c);
            } else {
                self.values.swap(result, a);
            }
        } else if (self.less)(self.values[a], self.values[c]) {
            self.values.swap(result, a);
        } else if (self.less)(self.values[b], self.values[c]) {
            self.values.swap(result, c);
        } else {
            self.values.swap(result, b);
        }
    }

    fn partition(&mut self, mut first: usize, mut last: usize, pivot: usize) -> usize {
        loop {
            while (self.less)(self.values[first], self.values[pivot]) {
                first += 1;
            }
            last -= 1;
            while (self.less)(self.values[pivot], self.values[last]) {
                last -= 1;
            }
            if first >= last {
                return first;
            }
            self.values.swap(first, last);
            first += 1;
        }
    }

    fn insertion_sort(&mut self, first: usize, last: usize) {
        if first == last {
            return;
        }
        for i in first + 1..last {
            let value = self.values[i];
            if (self.less)(value, self.values[first]) {
                self.values.copy_within(first..i, first + 1);
                self.values[first] = value;
            } else {
                let mut hole = i;
                while (self.less)(value, self.values[hole - 1]) {
                    self.values[hole] = self.values[hole - 1];
                    hole -= 1;
                }
                self.values[hole] = value;
            }
        }
    }

    fn unguarded_insertion_sort(&mut self, first: usize, last: usize) {
        for i in first..last {
            let value = self.values[i];
            let mut hole = i;
            while (self.less)(value, self.values[hole - 1]) {
                self.values[hole] = self.values[hole - 1];
                hole -= 1;
            }
            self.values[hole] = value;
        }
    }

    fn push_heap(&mut self, first: usize, mut hole: usize, top: usize, value: T) {
        let mut parent = hole.saturating_sub(1) / 2;
        while hole > top && (self.less)(self.values[first + parent], value) {
            self.values[first + hole] = self.values[first + parent];
            hole = parent;
            parent = hole.saturating_sub(1) / 2;
        }
        self.values[first + hole] = value;
    }

    fn adjust_heap(&mut self, first: usize, mut hole: usize, len: usize, value: T) {
        let top = hole;
        let mut second_child = hole;
        while second_child < (len - 1) / 2 {
            second_child = 2 * (second_child + 1);
            if (self.less)(
                self.values[first + second_child],
                self.values[first + second_child - 1],
            ) {
                second_child -= 1;
            }
            self.values[first + hole] = self.values[first + second_child];
            hole = second_child;
        }
        if len & 1 == 0 && second_child == (len - 2) / 2 {
            second_child = 2 * (second_child + 1);
            self.values[first + hole] = self.values[first + second_child - 1];
            hole = second_child - 1;
        }
        self.push_heap(first, hole, top, value);
    }

    fn pop_heap(&mut self, first: usize, last: usize, result: usize) {
        let value = self.values[result];
        self.values[result] = self.values[first];
        self.adjust_heap(first, 0, last - first, value);
    }

    fn make_heap(&mut self, first: usize, last: usize) {
        if last - first < 2 {
            return;
        }
        let len = last - first;
        let mut parent = (len - 2) / 2;
        loop {
            let value = self.values[first + parent];
            self.adjust_heap(first, parent, len, value);
            if parent == 0 {
                return;
            }
            parent -= 1;
        }
    }

    fn partial_sort(&mut self, first: usize, middle: usize, last: usize) {
        self.make_heap(first, middle);
        for i in middle..last {
            if (self.less)(self.values[i], self.values[first]) {
                self.pop_heap(first, middle, i);
            }
        }
        let mut heap_last = middle;
        while heap_last - first > 1 {
            heap_last -= 1;
            self.pop_heap(first, heap_last, heap_last);
        }
    }

    fn introsort(&mut self, first: usize, mut last: usize, mut depth: usize) {
        while last - first > THRESHOLD {
            if depth == 0 {
                self.partial_sort(first, last, last);
                return;
            }
            depth -= 1;
            let mid = first + (last - first) / 2;
            self.move_median_to_first(first, first + 1, mid, last - 1);
            let cut = self.partition(first + 1, last, first);
            self.introsort(cut, last, depth);
            last = cut;
        }
    }
}
