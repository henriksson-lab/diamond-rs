#![allow(non_snake_case)]

use std::fmt::{self, Display};
use std::str::FromStr;

use crate::stats::alp_ioutil;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Vector<T> {
    d_m: usize,
    d_vector_p: Vec<T>,
}

impl<T> Default for Vector<T> {
    fn default() -> Self {
        Self {
            d_m: 0,
            d_vector_p: Vec::new(),
        }
    }
}

impl<T: Clone + Default> Vector<T> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_slice(begin_: &[T]) -> Self {
        let mut vector = Self::default();
        vector.copy(begin_);
        vector
    }

    pub fn with_value(m_: usize, a_: T) -> Self {
        let mut vector = Self::default();
        vector.copy_value(m_, a_);
        vector
    }

    pub fn create(&self, isCopy_: bool) -> Box<Self> {
        if isCopy_ {
            Box::new(self.clone())
        } else {
            Box::new(Self::default())
        }
    }

    pub fn copy(&mut self, begin_: &[T]) {
        if begin_.len() != self.getM() {
            self.free2();
            self.init(begin_.len());
        }
        for (dst, src) in self.d_vector_p.iter_mut().zip(begin_) {
            *dst = src.clone();
        }
    }

    pub fn copy_value(&mut self, m_: usize, a_: T) {
        if m_ != self.getM() {
            self.free2();
            self.init(m_);
        }
        for dst in &mut self.d_vector_p {
            *dst = a_.clone();
        }
    }

    pub fn assign_value(&mut self, a_: T) -> &mut Self {
        self.copy_value(self.getM(), a_);
        self
    }

    pub fn bool_(&self) -> bool {
        self.getM() > 0
    }

    pub fn size(&self) -> usize {
        self.getM()
    }

    pub fn getM(&self) -> usize {
        self.d_m
    }

    pub fn getVector(&self) -> &[T] {
        &self.d_vector_p
    }

    pub fn getVector_mut(&mut self) -> &mut [T] {
        &mut self.d_vector_p
    }

    fn init(&mut self, m_: usize) {
        self.d_vector_p = vec![T::default(); m_];
        self.d_m = m_;
    }

    fn free2(&mut self) {
        self.d_vector_p.clear();
        self.d_m = 0;
    }
}

impl<T> std::ops::Index<usize> for Vector<T> {
    type Output = T;

    fn index(&self, i_: usize) -> &Self::Output {
        &self.d_vector_p[i_]
    }
}

impl<T> std::ops::IndexMut<usize> for Vector<T> {
    fn index_mut(&mut self, i_: usize) -> &mut Self::Output {
        &mut self.d_vector_p[i_]
    }
}

impl<T> Vector<T>
where
    T: Clone + Default + Display + FromStr,
{
    pub fn out(&self) -> String {
        let mut out = String::new();
        match alp_ioutil::getFormat() {
            alp_ioutil::Format::MACHINE => {
                for i in 0..self.getM() {
                    if i != 0 {
                        out.push('\n');
                    }
                    out.push_str(&self.getVector()[i].to_string());
                }
            }
            alp_ioutil::Format::HUMAN | alp_ioutil::Format::FORMAT => {
                out.push_str(&format!("{}\t! dimension of vector\n", self.getM()));
                for i in 0..self.getM() {
                    if i != 0 {
                        out.push('\t');
                    }
                    out.push_str(&self.getVector()[i].to_string());
                }
                out.push_str("\t! vector elements");
            }
        }
        alp_ioutil::clearFormat();
        out
    }

    pub fn in_(&mut self, input: &str) {
        match alp_ioutil::getFormat() {
            alp_ioutil::Format::MACHINE => {
                let mut v = Vec::new();
                for token in input.split_whitespace() {
                    match token.parse::<T>() {
                        Ok(value) => v.push(value),
                        Err(_) => alp_ioutil::abort_msg("Njn::Vector::in : bad d_vector_p"),
                    }
                }
                if v.len() != self.getM() {
                    self.free2();
                    self.init(v.len());
                }
                self.d_vector_p.clone_from(&v);
            }
            alp_ioutil::Format::HUMAN | alp_ioutil::Format::FORMAT => {
                let mut offset = 0;
                let mut s = String::new();
                if !alp_ioutil::getLine(input, &mut offset, &mut s, alp_ioutil::getTerminator()) {
                    alp_ioutil::abort_msg("Njn::Vector::in : bad m");
                }
                let m: usize = s
                    .split_whitespace()
                    .next()
                    .and_then(|x| x.parse().ok())
                    .unwrap_or_else(|| alp_ioutil::abort_msg("Njn::Vector::in : bad m"));
                if m != self.getM() {
                    self.free2();
                    self.init(m);
                }

                if !alp_ioutil::getLine(input, &mut offset, &mut s, alp_ioutil::getTerminator()) {
                    alp_ioutil::abort_msg("Njn::Vector::in : bad d_vector_p [0]");
                }
                let mut tokens = s.split_whitespace();
                for i in 0..self.getM() {
                    self.d_vector_p[i] = tokens
                        .next()
                        .and_then(|x| x.parse::<T>().ok())
                        .unwrap_or_else(|| {
                            alp_ioutil::abort_msg(&format!(
                                "Njn::Vector::in : bad d_vector_p [{}]",
                                i
                            ))
                        });
                }
            }
        }
        alp_ioutil::clearFormat();
    }
}

impl<T> fmt::Display for Vector<T>
where
    T: Clone + Default + Display + FromStr,
{
    fn fmt(&self, ostr_: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(ostr_, "{}", self.out())
    }
}

impl<T> FromStr for Vector<T>
where
    T: Clone + Default + Display + FromStr,
{
    type Err = String;

    fn from_str(istr_: &str) -> Result<Self, Self::Err> {
        let mut vector_ = Self::new();
        vector_.in_(istr_);
        Ok(vector_)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

    #[test]
    fn test_vector_copy_create_assign_and_index() {
        let mut v = Vector::from_slice(&[1, 2, 3]);
        assert!(v.bool_());
        assert_eq!(v.size(), 3);
        assert_eq!(v[1], 2);
        v[1] = 9;
        assert_eq!(v.getVector(), &[1, 9, 3]);
        v.assign_value(4);
        assert_eq!(v.getVector(), &[4, 4, 4]);
        assert_eq!(*v.create(true), v);
        assert_eq!(v.create(false).size(), 0);
    }

    #[test]
    fn test_vector_human_and_machine_io() {
        let _guard = TEST_LOCK.lock().unwrap();
        let v = Vector::from_slice(&[1, 2, 3]);
        alp_ioutil::setFormat(alp_ioutil::Format::HUMAN);
        assert_eq!(
            v.out(),
            "3\t! dimension of vector\n1\t2\t3\t! vector elements"
        );

        let mut parsed = Vector::<i32>::new();
        alp_ioutil::setFormat(alp_ioutil::Format::HUMAN);
        parsed.in_("3 ! m\n1 2 3 ! values\n");
        assert_eq!(parsed, v);

        alp_ioutil::setFormat(alp_ioutil::Format::MACHINE);
        assert_eq!(v.out(), "1\n2\n3");
        alp_ioutil::setFormat(alp_ioutil::Format::MACHINE);
        parsed.in_("4\n5\n6\n");
        assert_eq!(parsed.getVector(), &[4, 5, 6]);

        alp_ioutil::setFormat(alp_ioutil::Format::MACHINE);
        assert_eq!(format!("{}", v), "1\n2\n3");
        alp_ioutil::setFormat(alp_ioutil::Format::MACHINE);
        let parsed_from_stream: Vector<i32> = "7\n8\n9\n".parse().unwrap();
        assert_eq!(parsed_from_stream.getVector(), &[7, 8, 9]);
    }
}
