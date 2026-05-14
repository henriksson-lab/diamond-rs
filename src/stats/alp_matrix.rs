#![allow(non_snake_case)]

use std::fmt::Display;
use std::str::FromStr;
use std::sync::Mutex;

use crate::stats::alp_approx;
use crate::stats::alp_function::AlpFloat;
use crate::stats::alp_ioutil;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    GENERAL,
    SYMMETRIC,
    FORMAT,
}

const FORMAT_DEFAULT: Format = Format::GENERAL;
static MATRIX_FORMAT_STATE: Mutex<Format> = Mutex::new(FORMAT_DEFAULT);

pub fn getFormat() -> Format {
    *MATRIX_FORMAT_STATE.lock().unwrap()
}

pub fn setFormat(format_: Format) {
    *MATRIX_FORMAT_STATE.lock().unwrap() = format_;
}

pub fn clearFormat() -> Format {
    let mut format = MATRIX_FORMAT_STATE.lock().unwrap();
    *format = FORMAT_DEFAULT;
    *format
}

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<T> {
    d_m: usize,
    d_n: usize,
    d_matrix_p: Vec<Vec<T>>,
    d_value: T,
}

impl<T: Default> Default for Matrix<T> {
    fn default() -> Self {
        Self {
            d_m: 0,
            d_n: 0,
            d_matrix_p: Vec::new(),
            d_value: T::default(),
        }
    }
}

impl<T: Clone + Default> Matrix<T> {
    pub fn matrix(k_: usize, m_: usize, n_: usize, a_: T) -> Vec<Self> {
        let mut matrix = vec![Self::default(); k_];
        for item in &mut matrix {
            item.copy_value(m_, n_, a_.clone());
        }
        matrix
    }

    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_vector(m_: usize, n_: usize, vector_: &[T]) -> Self {
        let mut matrix = Self::default();
        matrix.copy_vector(m_, n_, vector_);
        matrix
    }

    pub fn with_value(m_: usize, n_: usize, a_: T) -> Self {
        let mut matrix = Self::default();
        matrix.copy_value(m_, n_, a_);
        matrix
    }

    pub fn create(&self, isCopy_: bool) -> Box<Self> {
        if isCopy_ {
            Box::new(self.clone())
        } else {
            Box::new(Self::default())
        }
    }

    pub fn copy_matrix(&mut self, m_: usize, n_: usize, matrix_: &[Vec<T>]) {
        if m_ != self.getM() || n_ != self.getN() {
            self.free2();
            self.init(m_, n_);
        }
        for i in 0..self.getM() {
            for j in 0..self.getN() {
                self.d_matrix_p[i][j] = matrix_[i][j].clone();
            }
        }
    }

    pub fn copy_vector(&mut self, m_: usize, n_: usize, vector_: &[T]) {
        if m_ != self.getM() || n_ != self.getN() {
            self.free2();
            self.init(m_, n_);
        }
        let mut k = 0;
        for i in 0..self.getM() {
            for j in 0..self.getN() {
                self.d_matrix_p[i][j] = vector_[k].clone();
                k += 1;
            }
        }
    }

    pub fn copy_value(&mut self, m_: usize, n_: usize, a_: T) {
        if m_ != self.getM() || n_ != self.getN() {
            self.free2();
            self.init(m_, n_);
        }
        for i in 0..self.getM() {
            for j in 0..self.getN() {
                self.d_matrix_p[i][j] = a_.clone();
            }
        }
    }

    pub fn bool_(&self) -> bool {
        self.getM() > 0 && self.getN() > 0
    }

    pub fn assign_value(&mut self, a_: T) -> &mut Self {
        self.copy_value(self.getM(), self.getN(), a_);
        self
    }

    pub fn setValue(&mut self) -> &mut T {
        &mut self.d_value
    }

    pub fn getM(&self) -> usize {
        self.d_m
    }

    pub fn getN(&self) -> usize {
        self.d_n
    }

    pub fn getMatrix(&self) -> &[Vec<T>] {
        &self.d_matrix_p
    }

    pub fn getValue(&self) -> &T {
        &self.d_value
    }

    fn init(&mut self, m_: usize, n_: usize) {
        self.d_matrix_p = vec![vec![T::default(); n_]; m_];
        self.d_m = m_;
        self.d_n = n_;
    }

    fn free2(&mut self) {
        self.d_matrix_p.clear();
        self.d_m = 0;
        self.d_n = 0;
    }
}

impl<T: AlpFloat + Default> Matrix<T> {
    pub fn approx(x_: &Self, y_: &Self, eps_: T) -> bool {
        assert_eq!(x_.getM(), y_.getM());
        assert_eq!(x_.getN(), y_.getN());
        for i in 0..x_.getM() {
            for j in 0..x_.getN() {
                if !alp_approx::approx(x_[i][j], y_[i][j], eps_) {
                    return false;
                }
            }
        }
        true
    }

    pub fn relApprox(x_: &Self, y_: &Self, eps_: T) -> bool {
        assert_eq!(x_.getM(), y_.getM());
        assert_eq!(x_.getN(), y_.getN());
        for i in 0..x_.getM() {
            for j in 0..x_.getN() {
                if !alp_approx::relApprox(x_[i][j], y_[i][j], eps_) {
                    return false;
                }
            }
        }
        true
    }

    pub fn absRelApprox(x_: &Self, y_: &Self, tol_: T, rtol_: T) -> bool {
        assert_eq!(x_.getM(), y_.getM());
        assert_eq!(x_.getN(), y_.getN());
        for i in 0..x_.getM() {
            for j in 0..x_.getN() {
                if !alp_approx::absRelApprox(x_[i][j], y_[i][j], tol_, rtol_) {
                    return false;
                }
            }
        }
        true
    }
}

impl<T: Clone + Default + PartialEq> Matrix<T> {
    pub fn isSymmetric(x_: &Self) -> bool {
        if x_.getM() != x_.getN() {
            return false;
        }
        for i in 0..x_.getM() {
            for j in 0..i {
                if x_[i][j] != x_[j][i] {
                    return false;
                }
            }
        }
        true
    }
}

impl<T> std::ops::Index<usize> for Matrix<T> {
    type Output = [T];

    fn index(&self, i_: usize) -> &Self::Output {
        &self.d_matrix_p[i_]
    }
}

impl<T> std::ops::IndexMut<usize> for Matrix<T> {
    fn index_mut(&mut self, i_: usize) -> &mut Self::Output {
        &mut self.d_matrix_p[i_]
    }
}

impl<T> Matrix<T>
where
    T: Clone + Default + Display + FromStr + PartialEq,
{
    pub fn out(&self) -> String {
        let mut out = String::new();
        match getFormat() {
            Format::GENERAL => match alp_ioutil::getFormat() {
                alp_ioutil::Format::MACHINE => {
                    for i in 0..self.getM() {
                        if i != 0 {
                            out.push('\n');
                        }
                        for j in 0..self.getN() {
                            if j != 0 {
                                out.push('\t');
                            }
                            out.push_str(&self.getMatrix()[i][j].to_string());
                        }
                    }
                }
                alp_ioutil::Format::HUMAN | alp_ioutil::Format::FORMAT => {
                    out.push_str(&format!("{}\t! row dimension of matrix\n", self.getM()));
                    out.push_str(&format!("{}\t! column dimension of matrix\n", self.getN()));
                    for i in 0..self.getM() {
                        if i != 0 {
                            out.push_str("\t! matrix elements\n");
                        }
                        for j in 0..self.getN() {
                            if j != 0 {
                                out.push('\t');
                            }
                            out.push_str(&self.getMatrix()[i][j].to_string());
                        }
                    }
                    out.push_str("\t! matrix elements");
                }
            },
            Format::SYMMETRIC => {
                if !Self::isSymmetric(self) {
                    alp_ioutil::abort_msg("Matrix::out : matrix is not symmetric");
                }
                for i in 0..self.getM() {
                    if i != 0 {
                        out.push('\n');
                    }
                    for j in i..self.getN() {
                        if j != i {
                            out.push('\t');
                        }
                        out.push_str(&self.getMatrix()[i][j].to_string());
                    }
                }
            }
            Format::FORMAT => {
                alp_ioutil::abort_msg("Matrix::out : impossible MatrixIO::getFormat ()")
            }
        }
        alp_ioutil::clearFormat();
        clearFormat();
        out
    }

    pub fn in_(&mut self, input: &str) {
        match getFormat() {
            Format::GENERAL => match alp_ioutil::getFormat() {
                alp_ioutil::Format::MACHINE => self.in_machine(input),
                alp_ioutil::Format::HUMAN | alp_ioutil::Format::FORMAT => self.in_human(input),
            },
            Format::SYMMETRIC => self.in_symmetric(input),
            Format::FORMAT => {
                alp_ioutil::abort_msg("Matrix::in : impossible MatrixIO::getFormat ()")
            }
        }
        alp_ioutil::clearFormat();
        clearFormat();
    }

    fn in_machine(&mut self, input: &str) {
        let mut lines = input.lines();
        let Some(first_line) = lines.next() else {
            alp_ioutil::abort_msg("Njn::Matrix::in : bad n for the MACHINE");
        };
        let mut v: Vec<T> = first_line
            .split_whitespace()
            .map(|x| {
                x.parse::<T>().unwrap_or_else(|_| {
                    alp_ioutil::abort_msg("Njn::Matrix::in : bad value for the MACHINE")
                })
            })
            .collect();
        let n = v.len();
        if n == 0 {
            alp_ioutil::abort_msg("Njn::Matrix::in : bad n for the MACHINE");
        }
        for token in lines.flat_map(|line| line.split_whitespace()) {
            v.push(token.parse::<T>().unwrap_or_else(|_| {
                alp_ioutil::abort_msg("Njn::Matrix::in : bad value for the MACHINE")
            }));
        }
        let m = v.len() / n;
        if m * n != v.len() {
            alp_ioutil::abort_msg("Njn::Matrix::in : rows for the MACHINE have different lengths.");
        }
        self.copy_vector(m, n, &v);
    }

    fn in_human(&mut self, input: &str) {
        let mut offset = 0;
        let mut s = String::new();
        if !alp_ioutil::getLine(input, &mut offset, &mut s, alp_ioutil::getTerminator()) {
            alp_ioutil::abort_msg("Njn::Matrix::in : bad m");
        }
        let m: usize = s
            .split_whitespace()
            .next()
            .and_then(|x| x.parse().ok())
            .unwrap_or_else(|| alp_ioutil::abort_msg("Njn::Matrix::in : bad m"));
        if !alp_ioutil::getLine(input, &mut offset, &mut s, alp_ioutil::getTerminator()) {
            alp_ioutil::abort_msg("Njn::Matrix::in : bad n");
        }
        let n: usize = s
            .split_whitespace()
            .next()
            .and_then(|x| x.parse().ok())
            .unwrap_or_else(|| alp_ioutil::abort_msg("Njn::Matrix::in : bad n"));

        if m != self.getM() || n != self.getN() {
            self.free2();
            self.init(m, n);
        }

        for i in 0..self.getM() {
            if !alp_ioutil::getLine(input, &mut offset, &mut s, alp_ioutil::getTerminator()) {
                alp_ioutil::abort_msg(&format!("Njn::Matrix::in : bad d_matrix_p [{}][0]", i));
            }
            let mut tokens = s.split_whitespace();
            for j in 0..self.getN() {
                self.d_matrix_p[i][j] = tokens
                    .next()
                    .and_then(|x| x.parse::<T>().ok())
                    .unwrap_or_else(|| {
                        alp_ioutil::abort_msg(&format!(
                            "Njn::Matrix::in : bad d_matrix_p [{}][{}]",
                            i, j
                        ))
                    });
            }
        }
    }

    fn in_symmetric(&mut self, input: &str) {
        let mut lines = input.lines();
        let Some(first_line) = lines.next() else {
            alp_ioutil::abort_msg("Njn::Matrix::in : bad n for the MatrixIO::SYMMETRIC");
        };
        let mut v: Vec<T> = first_line
            .split_whitespace()
            .map(|x| {
                x.parse::<T>().unwrap_or_else(|_| {
                    alp_ioutil::abort_msg("Njn::Matrix::in : bad value for the MatrixIO::SYMMETRIC")
                })
            })
            .collect();
        let m = v.len();
        let n = m;
        if n == 0 {
            alp_ioutil::abort_msg("Njn::Matrix::in : bad n for the MatrixIO::SYMMETRIC");
        }
        for token in lines.flat_map(|line| line.split_whitespace()) {
            v.push(token.parse::<T>().unwrap_or_else(|_| {
                alp_ioutil::abort_msg("Njn::Matrix::in : bad value for the MatrixIO::SYMMETRIC")
            }));
        }
        if m * (m + 1) / 2 != v.len() {
            alp_ioutil::abort_msg(
                "Njn::Matrix::in : rows for the MatrixIO::SYMMETRIC have incorrect lengths.",
            );
        }
        if m != self.getM() || n != self.getN() {
            self.free2();
            self.init(m, n);
        }
        let mut k = 0;
        for i in 0..m {
            for j in 0..i {
                self.d_matrix_p[i][j] = self.d_matrix_p[j][i].clone();
            }
            for j in i..n {
                self.d_matrix_p[i][j] = v[k].clone();
                k += 1;
            }
        }
    }
}

pub fn copy<S, T>(matrix_: &mut Matrix<S>, matrix0_: &Matrix<T>)
where
    S: Clone + Default + From<T>,
    T: Clone + Default,
{
    matrix_.copy_value(matrix0_.getM(), matrix0_.getN(), S::default());
    for i in 0..matrix0_.getM() {
        for j in 0..matrix0_.getN() {
            matrix_[i][j] = S::from(matrix0_[i][j].clone());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static TEST_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn test_matrix_constructors_copy_create_and_index() {
        let mut m = Matrix::from_vector(2, 3, &[1, 2, 3, 4, 5, 6]);
        assert!(m.bool_());
        assert_eq!(m.getM(), 2);
        assert_eq!(m.getN(), 3);
        assert_eq!(m[1][2], 6);
        m[1][2] = 9;
        assert_eq!(m[1][2], 9);
        m.assign_value(7);
        assert_eq!(m.getMatrix(), &[vec![7, 7, 7], vec![7, 7, 7]]);
        assert_eq!(*m.create(true), m);
        assert_eq!(m.create(false).getM(), 0);

        let mats = Matrix::matrix(2, 1, 2, 4);
        assert_eq!(mats.len(), 2);
        assert_eq!(mats[0].getMatrix(), &[vec![4, 4]]);
    }

    #[test]
    fn test_matrix_approx_and_symmetric() {
        let x = Matrix::from_vector(2, 2, &[1.0, 2.0, 2.0, 4.0]);
        let y = Matrix::from_vector(2, 2, &[1.01, 2.0, 2.0, 4.0]);
        assert!(Matrix::approx(&x, &y, 0.02));
        assert!(Matrix::relApprox(&x, &y, 0.02));
        assert!(Matrix::absRelApprox(&x, &y, 0.001, 0.02));
        assert!(Matrix::isSymmetric(&x));
        let z = Matrix::from_vector(2, 2, &[1.0, 3.0, 2.0, 4.0]);
        assert!(!Matrix::isSymmetric(&z));
    }

    #[test]
    fn test_matrix_general_human_and_machine_io() {
        let _guard = TEST_LOCK.lock().unwrap();
        let m = Matrix::from_vector(2, 2, &[1, 2, 3, 4]);
        alp_ioutil::setFormat(alp_ioutil::Format::HUMAN);
        setFormat(Format::GENERAL);
        assert_eq!(
            m.out(),
            "2\t! row dimension of matrix\n2\t! column dimension of matrix\n1\t2\t! matrix elements\n3\t4\t! matrix elements"
        );

        let mut parsed = Matrix::<i32>::new();
        alp_ioutil::setFormat(alp_ioutil::Format::HUMAN);
        setFormat(Format::GENERAL);
        parsed.in_("2 ! rows\n2 ! cols\n1 2 ! row\n3 4 ! row\n");
        assert_eq!(parsed, m);

        alp_ioutil::setFormat(alp_ioutil::Format::MACHINE);
        setFormat(Format::GENERAL);
        assert_eq!(m.out(), "1\t2\n3\t4");
        alp_ioutil::setFormat(alp_ioutil::Format::MACHINE);
        setFormat(Format::GENERAL);
        parsed.in_("5 6\n7 8\n");
        assert_eq!(parsed.getMatrix(), &[vec![5, 6], vec![7, 8]]);
    }

    #[test]
    fn test_matrix_symmetric_io() {
        let _guard = TEST_LOCK.lock().unwrap();
        let m = Matrix::from_vector(2, 2, &[1, 2, 2, 4]);
        setFormat(Format::SYMMETRIC);
        assert_eq!(m.out(), "1\t2\n4");

        let mut parsed = Matrix::<i32>::new();
        setFormat(Format::SYMMETRIC);
        parsed.in_("1 2\n4\n");
        assert_eq!(parsed, m);
    }
}
