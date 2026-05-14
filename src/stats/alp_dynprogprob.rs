#![allow(non_snake_case)]

pub type ValueFct = fn(i64, usize) -> i64;

#[derive(Clone)]
pub struct DynProgProb {
    d_step: usize,
    d_array_p: [Vec<f64>; 2],
    d_arrayCapacity: usize,
    d_valueBegin: i64,
    d_valueLower: i64,
    d_valueUpper: i64,
    d_valueFct: Option<ValueFct>,
    d_dimInputProb: usize,
    d_inputProb_p: Vec<f64>,
}

impl crate::stats::alp_dynprogprobproto::DynProgProbProto for DynProgProb {
    fn bool_(&self) -> bool {
        self.bool_()
    }

    fn clear_default(&mut self) {
        self.clear_default();
    }

    fn update(&mut self) {
        self.update();
    }

    fn getProb(&self, value_: i64) -> f64 {
        self.getProb(value_)
    }

    fn getStep(&self) -> usize {
        self.getStep()
    }

    fn getValueLower(&self) -> i64 {
        self.getValueLower()
    }

    fn getValueUpper(&self) -> i64 {
        self.getValueUpper()
    }
}

impl std::fmt::Debug for DynProgProb {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DynProgProb")
            .field("d_step", &self.d_step)
            .field("d_arrayCapacity", &self.d_arrayCapacity)
            .field("d_valueBegin", &self.d_valueBegin)
            .field("d_valueLower", &self.d_valueLower)
            .field("d_valueUpper", &self.d_valueUpper)
            .field("d_dimInputProb", &self.d_dimInputProb)
            .finish()
    }
}

impl DynProgProb {
    pub const ARRAY_CAPACITY: usize = 256;

    pub fn new(
        valueFct_: Option<ValueFct>,
        dimInputProb_: usize,
        inputProb_: Option<&[f64]>,
        valueLower_: i64,
        valueUpper_: i64,
        prob_: Option<&[f64]>,
    ) -> Self {
        let mut this = Self {
            d_step: 0,
            d_array_p: [Vec::new(), Vec::new()],
            d_arrayCapacity: 0,
            d_valueBegin: 0,
            d_valueLower: 0,
            d_valueUpper: 0,
            d_valueFct: None,
            d_dimInputProb: 0,
            d_inputProb_p: Vec::new(),
        };
        this.setValueFct(valueFct_);
        this.setInput(dimInputProb_, inputProb_);
        this.clear(valueLower_, valueUpper_, prob_);
        this
    }

    pub fn bool_(&self) -> bool {
        self.getArrayCapacity() != 0
            && self.d_valueFct.is_some()
            && self.d_dimInputProb != 0
            && !self.d_inputProb_p.is_empty()
    }

    pub fn copy(&mut self, dynProgProb_: &Self) {
        *self = dynProgProb_.clone();
    }

    pub fn clear(&mut self, valueLower_: i64, valueUpper_: i64, prob_: Option<&[f64]>) {
        if let Some(prob) = prob_ {
            assert!(valueLower_ < valueUpper_);
            for &p in prob {
                assert!(0.0 <= p);
            }
            self.clear_capacity(valueLower_, (valueUpper_ - valueLower_) as usize);
            self.d_valueLower = valueLower_;
            self.d_valueUpper = valueUpper_;
            let array_capacity = self.getArrayCapacity();
            self.d_array_p[0].copy_from_slice(&prob[..array_capacity]);
            return;
        }

        assert!(valueLower_ <= 0 && 0 <= valueUpper_);
        if valueLower_ == 0 && valueUpper_ == 0 {
            self.clear_capacity(-(Self::ARRAY_CAPACITY as i64 / 2) + 1, Self::ARRAY_CAPACITY);
        } else {
            self.clear_capacity(valueLower_, (valueUpper_ - valueLower_) as usize);
        }
        self.d_valueLower = 0;
        self.d_valueUpper = 1;
        let pos = self.getArrayPos(0) as usize;
        self.d_array_p[0][pos] = 1.0;
    }

    pub fn clear_default(&mut self) {
        self.clear(0, 0, None);
    }

    pub fn setValueFct(&mut self, valueFct_: Option<ValueFct>) {
        self.d_valueFct = valueFct_;
    }

    pub fn setInput(&mut self, dimInputProb_: usize, inputProb_: Option<&[f64]>) {
        if dimInputProb_ != self.getDimInputProb() {
            self.d_inputProb_p = vec![0.0; dimInputProb_];
            self.d_dimInputProb = dimInputProb_;
        }
        if self.getDimInputProb() > 0 {
            let dim_input_prob = self.getDimInputProb();
            self.d_inputProb_p
                .copy_from_slice(&inputProb_.expect("inputProb_")[..dim_input_prob]);
        }
    }

    pub fn update(&mut self) {
        assert!(self.getValueFct().is_some());
        assert!(self.getDimInputProb() != 0);
        assert!(!self.getInputProb().is_empty());

        const ARRAY_FAC: usize = 2;
        let value_fct = self.getValueFct().unwrap();
        let old_idx = self.d_step % 2;
        let new_idx = (self.d_step + 1) % 2;
        self.d_array_p[new_idx].fill(0.0);

        let mut valueLower = i64::MAX;
        let mut valueUpper = i64::MIN;

        for i in self.getValueLower()..self.getValueUpper() {
            let old_prob = self.d_array_p[old_idx][self.getArrayPos(i) as usize];
            if old_prob == 0.0 {
                continue;
            }

            for j in 0..self.getDimInputProb() {
                let input_prob = self.d_inputProb_p[j];
                if input_prob == 0.0 {
                    continue;
                }

                let value = value_fct(i, j);
                while value < self.getValueBegin() || self.getValueEnd() <= value {
                    let mut valueBegin = self.getValueBegin();
                    if value < self.getValueBegin() {
                        valueBegin -= (ARRAY_FAC as i64 - 1) * self.getArrayCapacity() as i64;
                    }
                    self.reserve(ARRAY_FAC * self.getArrayCapacity());
                    self.setValueBegin(valueBegin);
                }

                if value < valueLower {
                    valueLower = value;
                }
                if valueUpper < value {
                    valueUpper = value;
                }
                let pos = self.getArrayPos(value) as usize;
                self.d_array_p[new_idx][pos] += old_prob * input_prob;
            }
        }

        self.d_valueLower = valueLower;
        self.d_valueUpper = valueUpper + 1;
        self.d_step += 1;
    }

    pub fn getProb(&self, value_: i64) -> f64 {
        if value_ < self.getValueBegin() {
            return 0.0;
        }
        if self.getValueEnd() <= value_ {
            return 0.0;
        }
        self.d_array_p[self.getStep() % 2][self.getArrayPos(value_) as usize]
    }

    pub fn getStep(&self) -> usize {
        self.d_step
    }

    pub fn getArray(&self) -> &[Vec<f64>; 2] {
        &self.d_array_p
    }

    pub fn getArrayCapacity(&self) -> usize {
        self.d_arrayCapacity
    }

    pub fn getValueBegin(&self) -> i64 {
        self.d_valueBegin
    }

    pub fn getValueLower(&self) -> i64 {
        self.d_valueLower
    }

    pub fn getValueUpper(&self) -> i64 {
        self.d_valueUpper
    }

    pub fn getValueFct(&self) -> Option<ValueFct> {
        self.d_valueFct
    }

    pub fn getDimInputProb(&self) -> usize {
        self.d_dimInputProb
    }

    pub fn getInputProb(&self) -> &[f64] {
        &self.d_inputProb_p
    }

    fn getValue(&self, arrayPos_: usize) -> i64 {
        arrayPos_ as i64 + self.getValueBegin()
    }

    fn clear_capacity(&mut self, valueBegin_: i64, arrayCapacity_: usize) {
        self.init(arrayCapacity_);
        self.d_valueBegin = valueBegin_;
        self.d_step = 0;
    }

    fn init(&mut self, arrayCapacity_: usize) {
        self.d_array_p = [vec![0.0; arrayCapacity_], vec![0.0; arrayCapacity_]];
        self.d_arrayCapacity = arrayCapacity_;
    }

    pub(crate) fn getArrayPos(&self, value_: i64) -> i64 {
        value_ - self.getValueBegin()
    }

    pub(crate) fn getValueEnd(&self) -> i64 {
        self.getValue(self.getArrayCapacity())
    }

    pub(crate) fn reserve(&mut self, arrayCapacity_: usize) {
        assert!(self.getArrayCapacity() < arrayCapacity_);
        for i in 0..2 {
            self.d_array_p[i].resize(arrayCapacity_, 0.0);
        }
        self.d_arrayCapacity = arrayCapacity_;
    }

    pub(crate) fn setValueBegin(&mut self, valueBegin_: i64) {
        assert!(valueBegin_ <= self.getValueBegin());
        let offSet = (self.getValueBegin() - valueBegin_) as usize;
        if offSet == 0 {
            return;
        }
        assert!(offSet < self.getArrayCapacity());
        for i in 0..2 {
            let old = self.d_array_p[i].clone();
            self.d_array_p[i].fill(0.0);
            let len = self.getArrayCapacity() - offSet;
            self.d_array_p[i][offSet..offSet + len].copy_from_slice(&old[..len]);
        }
        self.d_valueBegin = valueBegin_;
    }

    pub(crate) fn lgetStep(&mut self) -> &mut usize {
        &mut self.d_step
    }

    pub(crate) fn lgetArray(&mut self) -> &mut [Vec<f64>; 2] {
        &mut self.d_array_p
    }

    pub(crate) fn lgetArrayCapacity(&mut self) -> &mut usize {
        &mut self.d_arrayCapacity
    }

    pub(crate) fn lgetValueBegin(&mut self) -> &mut i64 {
        &mut self.d_valueBegin
    }

    pub(crate) fn lgetValueLower(&mut self) -> &mut i64 {
        &mut self.d_valueLower
    }

    pub(crate) fn lgetValueUpper(&mut self) -> &mut i64 {
        &mut self.d_valueUpper
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn add_input(old: i64, state: usize) -> i64 {
        old + state as i64
    }

    #[test]
    fn test_dyn_prog_prob_default_and_update() {
        let mut dp = DynProgProb::new(Some(add_input), 2, Some(&[0.5, 0.5]), 0, 0, None);
        assert!(dp.bool_());
        assert_eq!(dp.getProb(0), 1.0);
        dp.update();
        assert_eq!(dp.getStep(), 1);
        assert_eq!(dp.getProb(0), 0.5);
        assert_eq!(dp.getProb(1), 0.5);
        dp.update();
        assert!((dp.getProb(0) - 0.25).abs() < 1e-15);
        assert!((dp.getProb(1) - 0.5).abs() < 1e-15);
        assert!((dp.getProb(2) - 0.25).abs() < 1e-15);
    }

    #[test]
    fn test_dyn_prog_prob_initial_probabilities_and_reserve_left() {
        fn minus_input(old: i64, state: usize) -> i64 {
            old - state as i64 - 200
        }

        let mut dp = DynProgProb::new(
            Some(minus_input),
            2,
            Some(&[0.0, 1.0]),
            -1,
            1,
            Some(&[0.25, 0.75]),
        );
        assert_eq!(dp.getProb(-1), 0.25);
        assert_eq!(dp.getProb(0), 0.75);
        dp.update();
        assert!(dp.getValueBegin() < -200);
        assert_eq!(dp.getProb(-202), 0.25);
        assert_eq!(dp.getProb(-201), 0.75);
    }
}
