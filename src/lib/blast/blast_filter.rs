//! Safe ownership translation of NCBI `lib/blast/blast_filter.cpp`.
//!
//! C++ represents mask intervals as individually allocated linked nodes.
//! [`BlastSeqLocList`] stores the same ordered nodes contiguously; append,
//! reversal, duplication, sorting, merging and inclusive coordinates retain
//! the source behavior without raw allocation/free operations.

use std::cmp::Ordering;

pub const K_NUCL_MASK: u8 = 14;
pub const K_PROT_MASK: u8 = 21;
pub const BLASTOPTIONS_BUFFER_SIZE: usize = 128;
pub const REPEATS_SEARCH_EVALUE: f64 = 0.1;
pub const REPEATS_SEARCH_MINSCORE: i32 = 26;
pub const REPEATS_SEARCH_PENALTY: i32 = -1;
pub const REPEATS_SEARCH_REWARD: i32 = 1;
pub const REPEATS_SEARCH_GAP_OPEN: i32 = 2;
pub const REPEATS_SEARCH_GAP_EXTEND: i32 = 1;
pub const REPEATS_SEARCH_WORD_SIZE: i32 = 11;
pub const REPEATS_SEARCH_XDROP_UNGAPPED: i32 = 40;
pub const REPEATS_SEARCH_XDROP_FINAL: i32 = 90;
pub const REPEATS_SEARCH_FILTER_STRING: &str = "F";
pub const REPEAT_MASK_LINK_VALUE: i32 = 5;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SeqRange {
    pub left: i32,
    pub right: i32,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BlastSeqLoc {
    pub ssr: SeqRange,
}

pub type BlastSeqLocList = Vec<BlastSeqLoc>;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BlastMaskLoc {
    pub total_size: i32,
    pub seqloc_array: Vec<BlastSeqLocList>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BlastSequenceBlk {
    pub sequence: Vec<u8>,
    pub length: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BlastFilterError {
    NegativeLength,
    LengthExceedsBuffer { length: usize, buffer_len: usize },
    MaskRangeOutOfBounds { start: i32, stop: i32, length: i32 },
}

/// C++ `BlastSeqLocNew`. The returned value is the newly-created tail value;
/// when `head` is supplied, an identical owned node is appended to that list.
pub fn blast_seq_loc_new(head: Option<&mut BlastSeqLocList>, from: i32, to: i32) -> BlastSeqLoc {
    let node = BlastSeqLoc {
        ssr: SeqRange {
            left: from,
            right: to,
        },
    };
    blast_seq_loc_append(head, Some(node)).expect("new node is present")
}

/// C++ `BlastSeqLocAppend`. `None` corresponds to a null node.
pub fn blast_seq_loc_append(
    head: Option<&mut BlastSeqLocList>,
    node: Option<BlastSeqLoc>,
) -> Option<BlastSeqLoc> {
    let node = node?;
    if let Some(head) = head {
        head.push(node);
    }
    Some(node)
}

fn s_blast_seq_loc_node_dup(source: Option<&BlastSeqLoc>) -> Option<BlastSeqLoc> {
    source.copied()
}

fn s_blast_seq_loc_len(list: &[BlastSeqLoc]) -> i32 {
    list.len() as i32
}

/// The C++ helper returns pointers into the linked list plus a null sentinel.
/// Safe Rust uses indices into the unchanged list; `None` is the sentinel.
fn s_blast_seq_loc_list_to_array_of_pointers(list: &[BlastSeqLoc]) -> Option<Vec<Option<usize>>> {
    if list.is_empty() {
        return None;
    }
    let count = s_blast_seq_loc_len(list) as usize;
    let mut pointers = (0..count).map(Some).collect::<Vec<_>>();
    pointers.push(None);
    Some(pointers)
}

pub fn blast_seq_loc_list_reverse(head: Option<&mut BlastSeqLocList>) {
    let Some(head) = head else {
        return;
    };
    if s_blast_seq_loc_list_to_array_of_pointers(head).is_none() {
        return;
    }
    head.reverse();
}

/// Ownership replacement for `BlastSeqLocNodeFree`; dropping consumes the
/// node and the C++ null return is represented by `None`.
pub fn blast_seq_loc_node_free(_loc: Option<BlastSeqLoc>) -> Option<BlastSeqLoc> {
    None
}

/// Ownership replacement for `BlastSeqLocFree`.
pub fn blast_seq_loc_free(_loc: Option<BlastSeqLocList>) -> Option<BlastSeqLocList> {
    None
}

pub fn blast_seq_loc_list_dup(head: Option<&BlastSeqLocList>) -> Option<BlastSeqLocList> {
    head.filter(|head| !head.is_empty()).map(|head| {
        head.iter()
            .filter_map(|node| s_blast_seq_loc_node_dup(Some(node)))
            .collect()
    })
}

pub fn blast_mask_loc_new(total: i32) -> BlastMaskLoc {
    BlastMaskLoc {
        total_size: total,
        seqloc_array: if total > 0 {
            vec![Vec::new(); total as usize]
        } else {
            Vec::new()
        },
    }
}

pub fn blast_mask_loc_dup(mask_loc: Option<&BlastMaskLoc>) -> Option<BlastMaskLoc> {
    mask_loc.map(|mask_loc| BlastMaskLoc {
        total_size: mask_loc.total_size,
        seqloc_array: mask_loc
            .seqloc_array
            .iter()
            .map(|list| blast_seq_loc_list_dup(Some(list)).unwrap_or_default())
            .collect(),
    })
}

/// Ownership replacement for `BlastMaskLocFree`.
pub fn blast_mask_loc_free(_mask_loc: Option<BlastMaskLoc>) -> Option<BlastMaskLoc> {
    None
}

fn s_seq_range_sort_by_start_position(left: &BlastSeqLoc, right: &BlastSeqLoc) -> Ordering {
    left.ssr.left.cmp(&right.ssr.left)
}

pub fn blast_seq_loc_combine(mask_loc: &mut BlastSeqLocList, link_value: i32) {
    if s_blast_seq_loc_list_to_array_of_pointers(mask_loc).is_none() {
        return;
    }
    mask_loc.sort_unstable_by(s_seq_range_sort_by_start_position);
    let mut merged: Vec<BlastSeqLoc> = Vec::with_capacity(mask_loc.len());
    for node in mask_loc.drain(..) {
        if let Some(tail) = merged.last_mut() {
            let stop = tail.ssr.right;
            if stop.wrapping_add(link_value) > node.ssr.left {
                tail.ssr.right = stop.max(node.ssr.right);
                continue;
            }
        }
        merged.push(node);
    }
    *mask_loc = merged;
}

/// C++ `BlastSeqLocReverse`, including its sequential in-place assignments.
/// The second statement observes the newly-written `left` value; this is not
/// the usual coordinate swap, but is the literal vendor behavior.
pub fn blast_seq_loc_reverse(masks: &mut BlastSeqLocList, query_length: i32) {
    for mask in masks {
        mask.ssr.left = query_length - 1 - mask.ssr.right;
        mask.ssr.right = query_length - 1 - mask.ssr.left;
    }
}

pub fn blast_mask_the_residues(
    buffer: &mut [u8],
    length: i32,
    is_na: bool,
    mask_loc: &[BlastSeqLoc],
    reverse: bool,
    offset: i32,
) -> Result<(), BlastFilterError> {
    if length < 0 {
        return Err(BlastFilterError::NegativeLength);
    }
    if length as usize > buffer.len() {
        return Err(BlastFilterError::LengthExceedsBuffer {
            length: length as usize,
            buffer_len: buffer.len(),
        });
    }
    let masking_letter = if is_na { K_NUCL_MASK } else { K_PROT_MASK };
    for mask in mask_loc {
        let (mut start, mut stop) = if reverse {
            (length - 1 - mask.ssr.right, length - 1 - mask.ssr.left)
        } else {
            (mask.ssr.left, mask.ssr.right)
        };
        start -= offset;
        stop -= offset;
        if start < 0 || stop < start || stop >= length {
            return Err(BlastFilterError::MaskRangeOutOfBounds {
                start,
                stop,
                length,
            });
        }
        buffer[start as usize..=stop as usize].fill(masking_letter);
    }
    Ok(())
}

pub fn blast_mask_unsupported_aa(
    sequence: &mut BlastSequenceBlk,
    min_invalid: u8,
) -> Result<(), BlastFilterError> {
    if sequence.length < 0 {
        return Err(BlastFilterError::NegativeLength);
    }
    if sequence.length as usize > sequence.sequence.len() {
        return Err(BlastFilterError::LengthExceedsBuffer {
            length: sequence.length as usize,
            buffer_len: sequence.sequence.len(),
        });
    }
    for letter in &mut sequence.sequence[..sequence.length as usize] {
        if *letter >= min_invalid {
            *letter = K_PROT_MASK;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn loc(left: i32, right: i32) -> BlastSeqLoc {
        BlastSeqLoc {
            ssr: SeqRange { left, right },
        }
    }

    #[test]
    fn list_construction_append_duplication_and_order_reversal_preserve_values() {
        let mut list = Vec::new();
        assert_eq!(blast_seq_loc_new(Some(&mut list), 1, 3), loc(1, 3));
        assert_eq!(
            blast_seq_loc_append(Some(&mut list), Some(loc(5, 8))),
            Some(loc(5, 8))
        );
        assert_eq!(blast_seq_loc_append(Some(&mut list), None), None);
        assert_eq!(s_blast_seq_loc_len(&list), 2);
        assert_eq!(
            s_blast_seq_loc_list_to_array_of_pointers(&list),
            Some(vec![Some(0), Some(1), None])
        );
        let duplicate = blast_seq_loc_list_dup(Some(&list)).unwrap();
        assert_eq!(blast_seq_loc_list_dup(Some(&Vec::new())), None);
        blast_seq_loc_list_reverse(Some(&mut list));
        assert_eq!(list, vec![loc(5, 8), loc(1, 3)]);
        assert_eq!(duplicate, vec![loc(1, 3), loc(5, 8)]);
        assert_eq!(blast_seq_loc_node_free(Some(loc(0, 0))), None);
        assert_eq!(blast_seq_loc_free(Some(list)), None);
    }

    #[test]
    fn combine_sorts_and_uses_the_vendor_strict_link_boundary() {
        let mut list = vec![loc(10, 12), loc(0, 4), loc(9, 9), loc(11, 20)];
        blast_seq_loc_combine(&mut list, 5);
        // 4 + 5 == 9 is not linked because the source comparison is strict
        // `>`, while the remaining nearby/nested ranges merge.
        assert_eq!(list, vec![loc(0, 4), loc(9, 20)]);

        let mut touching = vec![loc(0, 5), loc(5, 7)];
        blast_seq_loc_combine(&mut touching, 0);
        assert_eq!(touching, vec![loc(0, 5), loc(5, 7)]);
    }

    #[test]
    fn coordinate_reverse_preserves_sequential_assignment_quirk() {
        let mut list = vec![loc(2, 5), loc(0, 9)];
        blast_seq_loc_reverse(&mut list, 10);
        assert_eq!(list, vec![loc(4, 5), loc(0, 9)]);
    }

    #[test]
    fn mask_residues_uses_inclusive_ranges_reverse_coordinates_and_offset() {
        let masks = vec![loc(2, 4)];
        let mut protein = vec![1; 8];
        blast_mask_the_residues(&mut protein, 8, false, &masks, false, 1).unwrap();
        assert_eq!(
            protein,
            vec![1, K_PROT_MASK, K_PROT_MASK, K_PROT_MASK, 1, 1, 1, 1]
        );

        let mut nucleotide = vec![2; 8];
        blast_mask_the_residues(&mut nucleotide, 8, true, &masks, true, 1).unwrap();
        assert_eq!(
            nucleotide,
            vec![2, 2, K_NUCL_MASK, K_NUCL_MASK, K_NUCL_MASK, 2, 2, 2]
        );
    }

    #[test]
    fn safe_masking_rejects_the_cpp_out_of_bounds_cases() {
        let mut buffer = vec![0; 4];
        assert_eq!(
            blast_mask_the_residues(&mut buffer, 4, false, &[loc(3, 4)], false, 0),
            Err(BlastFilterError::MaskRangeOutOfBounds {
                start: 3,
                stop: 4,
                length: 4
            })
        );
        assert_eq!(
            blast_mask_the_residues(&mut buffer, 5, false, &[], false, 0),
            Err(BlastFilterError::LengthExceedsBuffer {
                length: 5,
                buffer_len: 4
            })
        );
    }

    #[test]
    fn mask_location_copy_is_deep_and_unsupported_aa_honors_declared_length() {
        let mut masks = blast_mask_loc_new(2);
        masks.seqloc_array[0].push(loc(1, 2));
        let duplicate = blast_mask_loc_dup(Some(&masks)).unwrap();
        masks.seqloc_array[0][0].ssr.left = 99;
        assert_eq!(duplicate.seqloc_array[0], vec![loc(1, 2)]);
        assert_eq!(blast_mask_loc_free(Some(masks)), None);

        let mut sequence = BlastSequenceBlk {
            sequence: vec![0, 20, 21, 25, 30],
            length: 4,
        };
        blast_mask_unsupported_aa(&mut sequence, 21).unwrap();
        assert_eq!(sequence.sequence, vec![0, 20, K_PROT_MASK, K_PROT_MASK, 30]);
    }
}
