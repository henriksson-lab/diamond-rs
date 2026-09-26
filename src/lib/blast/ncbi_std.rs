//! Safe Rust translation of `diamond/src/lib/blast/ncbi_std.cpp` and its
//! numeric compatibility definitions from `ncbi_std.h`.

pub const TRUE: i32 = 1;
pub const FALSE: i32 = 0;
pub const UINT4_MAX: u32 = u32::MAX;
pub const INT4_MAX: i32 = i32::MAX;
pub const INT4_MIN: i32 = i32::MIN;
pub const INT2_MAX: i16 = i16::MAX;
pub const INT2_MIN: i16 = i16::MIN;
pub const INT1_MAX: i8 = i8::MAX;
pub const INT1_MIN: i8 = i8::MIN;
pub const NCBIMATH_LN2: f64 = 0.693_147_180_559_945_3;
pub const NULLB: u8 = 0;

#[inline]
pub fn ncbi_min<T: Ord>(a: T, b: T) -> T {
    if a > b {
        b
    } else {
        a
    }
}

#[inline]
pub fn ncbi_max<T: Ord>(a: T, b: T) -> T {
    if a >= b {
        a
    } else {
        b
    }
}

#[inline]
pub fn ncbi_abs<T>(a: T) -> T
where
    T: Copy + Default + Ord + std::ops::Neg<Output = T>,
{
    if a >= T::default() {
        a
    } else {
        -a
    }
}

#[inline]
pub fn ncbi_sign<T: Ord + Default>(a: T) -> i32 {
    match a.cmp(&T::default()) {
        std::cmp::Ordering::Greater => 1,
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
    }
}

#[inline]
pub const fn dim<T, const N: usize>(_: &[T; N]) -> usize {
    N
}

/// C `BlastMemDup`: null input and zero length map to `None`; allocation
/// failure also returns `None` through `try_reserve_exact`.
pub fn blast_mem_dup(orig: Option<&[u8]>, size: usize) -> Option<Vec<u8>> {
    let orig = orig?;
    if size == 0 || size > orig.len() {
        return None;
    }
    let mut copy = Vec::new();
    copy.try_reserve_exact(size).ok()?;
    copy.extend_from_slice(&orig[..size]);
    Some(copy)
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ListNode<T> {
    pub choice: i32,
    pub ptr: Option<T>,
    pub next: Option<Box<ListNode<T>>>,
}

impl<T> Default for ListNode<T> {
    fn default() -> Self {
        Self {
            choice: 0,
            ptr: None,
            next: None,
        }
    }
}

/// Append a zero-initialized node after the last node, matching `ListNodeNew`.
pub fn list_node_new<T>(head: &mut Option<Box<ListNode<T>>>) -> &mut ListNode<T> {
    let mut slot = head;
    while slot.is_some() {
        slot = &mut slot.as_mut().expect("checked").next;
    }
    slot.insert(Box::new(ListNode::default())).as_mut()
}

/// `ListNodeAdd` has the same list mutation after Rust makes the nullable
/// pointer-to-head an explicit mutable `Option`.
pub fn list_node_add<T>(head: &mut Option<Box<ListNode<T>>>) -> &mut ListNode<T> {
    list_node_new(head)
}

pub fn list_node_add_pointer<T>(
    head: &mut Option<Box<ListNode<T>>>,
    choice: u32,
    value: T,
) -> &mut ListNode<T> {
    let node = list_node_add(head);
    node.choice = choice as i32;
    node.ptr = Some(value);
    node
}

pub fn list_node_copy_str<'a>(
    head: &'a mut Option<Box<ListNode<String>>>,
    choice: i32,
    value: Option<&str>,
) -> Option<&'a mut ListNode<String>> {
    let value = value?;
    let node = list_node_add(head);
    node.choice = choice;
    node.ptr = Some(value.to_owned());
    Some(node)
}

/// Free only nodes. Rust returns attached values to the caller instead of
/// leaking the C-owned pointers; the resulting list is the C return value NULL.
pub fn list_node_free<T>(head: Option<Box<ListNode<T>>>) -> (Option<Box<ListNode<T>>>, Vec<T>) {
    let mut values = Vec::new();
    let mut current = head;
    while let Some(mut node) = current {
        if let Some(value) = node.ptr.take() {
            values.push(value);
        }
        current = node.next.take();
    }
    (None, values)
}

/// Free nodes and attached data. Ownership makes both operations deterministic.
pub fn list_node_free_data<T>(head: Option<Box<ListNode<T>>>) -> Option<Box<ListNode<T>>> {
    drop(head);
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn memory_dup_preserves_null_zero_prefix_and_independence() {
        assert_eq!(blast_mem_dup(None, 3), None);
        assert_eq!(blast_mem_dup(Some(&[1, 2]), 0), None);
        assert_eq!(blast_mem_dup(Some(&[1, 2]), 3), None);
        let source = [1, 2, 3];
        let mut copy = blast_mem_dup(Some(&source), 2).unwrap();
        copy[0] = 9;
        assert_eq!(source, [1, 2, 3]);
        assert_eq!(copy, [9, 2]);
    }

    #[test]
    fn node_add_copy_and_free_variants_preserve_order_and_data_policy() {
        let mut head = None;
        list_node_add_pointer(&mut head, 7, "borrowed".to_owned());
        assert!(list_node_copy_str(&mut head, -1, None).is_none());
        list_node_copy_str(&mut head, 2, Some("copied")).unwrap();
        assert_eq!(head.as_ref().unwrap().choice, 7);
        assert_eq!(head.as_ref().unwrap().next.as_ref().unwrap().choice, 2);
        let (null, values) = list_node_free(head);
        assert!(null.is_none());
        assert_eq!(values, ["borrowed", "copied"]);

        let mut owned = None;
        list_node_add_pointer(&mut owned, 1, String::from("drop me"));
        assert!(list_node_free_data(owned).is_none());
    }

    #[test]
    fn header_numeric_macros_retain_boundary_and_tie_semantics() {
        assert_eq!((TRUE, FALSE), (1, 0));
        assert_eq!(
            (UINT4_MAX, INT4_MAX, INT4_MIN),
            (u32::MAX, i32::MAX, i32::MIN)
        );
        assert_eq!(
            (INT2_MAX, INT2_MIN, INT1_MAX, INT1_MIN),
            (32767, -32768, 127, -128)
        );
        assert_eq!(ncbi_min(4, 4), 4);
        assert_eq!(ncbi_max(4, 4), 4);
        assert_eq!(ncbi_abs(-9), 9);
        assert_eq!((ncbi_sign(-4), ncbi_sign(0), ncbi_sign(4)), (-1, 0, 1));
        assert_eq!(dim(&[0_u8; 11]), 11);
        assert_eq!(NULLB, b'\0');
        assert_eq!(NCBIMATH_LN2, std::f64::consts::LN_2);
    }
}
