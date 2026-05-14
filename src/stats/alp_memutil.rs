#![allow(non_snake_case)]

pub fn memNew(size: usize) -> Option<Vec<u8>> {
    if size == 0 {
        None
    } else {
        Some(vec![0u8; size])
    }
}

pub fn memMore(p: Option<Vec<u8>>, size: usize) -> Vec<u8> {
    let mut p = p.unwrap_or_default();
    p.resize(size, 0);
    p
}

pub fn memCpy(to: &mut [u8], from: &[u8], size: usize) {
    if size != 0 {
        to[..size].copy_from_slice(&from[..size]);
    }
}

pub fn memSet(p: &mut [u8], c: i32, size: usize) {
    for byte in p.iter_mut().take(size) {
        *byte = c as u8;
    }
}

pub fn memZero(p: &mut [u8], size: usize) {
    for byte in p.iter_mut().take(size) {
        *byte = 0;
    }
}

pub fn newMatrix<T: Default + Clone>(m_: usize, n_: usize) -> Vec<Vec<T>> {
    vec![vec![T::default(); n_]; m_]
}

pub fn deleteMatrix<T>(matrix: &mut Vec<Vec<T>>, _m_: usize, _n_: usize) {
    matrix.clear();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mem_new_more_copy_set_zero() {
        assert_eq!(memNew(0), None);
        assert_eq!(memNew(3), Some(vec![0, 0, 0]));

        let mut buf = memMore(Some(vec![1, 2]), 4);
        assert_eq!(buf, vec![1, 2, 0, 0]);
        buf = memMore(Some(buf), 1);
        assert_eq!(buf, vec![1]);

        let mut dst = [0u8; 4];
        memCpy(&mut dst, &[9, 8, 7, 6], 3);
        assert_eq!(dst, [9, 8, 7, 0]);
        memSet(&mut dst, 0xff, 2);
        assert_eq!(dst, [255, 255, 7, 0]);
        memZero(&mut dst, 3);
        assert_eq!(dst, [0, 0, 0, 0]);
    }

    #[test]
    fn test_matrix_allocation_and_delete() {
        let mut matrix = newMatrix::<i32>(2, 3);
        assert_eq!(matrix, vec![vec![0, 0, 0], vec![0, 0, 0]]);
        matrix[1][2] = 7;
        deleteMatrix(&mut matrix, 2, 3);
        assert!(matrix.is_empty());
    }
}
