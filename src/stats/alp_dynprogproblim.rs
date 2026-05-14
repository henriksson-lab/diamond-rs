#![allow(non_snake_case)]

use crate::stats::alp_dynprogprob::{DynProgProb, ValueFct};

#[derive(Debug, Clone)]
pub struct DynProgProbLim {
    base: DynProgProb,
    d_probLost: f64,
}

impl DynProgProbLim {
    pub fn new(
        valueFct_: Option<ValueFct>,
        dimInputProb_: usize,
        inputProb_: Option<&[f64]>,
        valueLower_: i64,
        valueUpper_: i64,
        prob_: Option<&[f64]>,
    ) -> Self {
        Self {
            base: DynProgProb::new(
                valueFct_,
                dimInputProb_,
                inputProb_,
                valueLower_,
                valueUpper_,
                prob_,
            ),
            d_probLost: 0.0,
        }
    }

    pub fn copy(&mut self, dynProgProbLim_: &Self) {
        self.base.copy(&dynProgProbLim_.base);
        self.d_probLost = dynProgProbLim_.getProbLost();
    }

    pub fn copy_base(&mut self, dynProgProb_: &DynProgProb, probLost_: f64) {
        self.base.copy(dynProgProb_);
        self.d_probLost = probLost_;
    }

    pub fn setLimits(&mut self, valueBegin_: i64, valueEnd_: i64) {
        assert!(valueBegin_ < valueEnd_);

        for value in self.getValueLower()..valueBegin_ {
            self.d_probLost += self.getProb(value);
        }
        for value in valueEnd_..self.getValueUpper() {
            self.d_probLost += self.getProb(value);
        }

        let arrayCapacity = (valueEnd_ - valueBegin_) as usize;
        if self.getArrayCapacity() <= arrayCapacity {
            self.reserve(arrayCapacity);
            self.setValueBegin(valueBegin_);
        } else {
            self.setValueBegin(valueBegin_);
            self.reserve(arrayCapacity);
        }
    }

    pub fn update(&mut self) {
        assert!(self.getValueFct().is_some());
        assert!(self.getDimInputProb() != 0);
        assert!(!self.getInputProb().is_empty());
        assert!(self.getArrayCapacity() > 0);

        let value_fct = self.getValueFct().unwrap();
        let old_idx = self.getStep() % 2;
        let new_idx = (self.getStep() + 1) % 2;
        self.base.lgetArray()[new_idx].fill(0.0);
        let old_array = self.base.getArray()[old_idx].clone();
        let input_prob = self.getInputProb().to_vec();

        let mut valueLower = i64::MAX;
        let mut valueUpper = i64::MIN;
        for i in self.getValueLower()..self.getValueUpper() {
            let old_prob = old_array[self.base.getArrayPos(i) as usize];
            if old_prob == 0.0 {
                continue;
            }
            for (j, &ip) in input_prob.iter().enumerate() {
                if ip == 0.0 {
                    continue;
                }
                let value = value_fct(i, j);
                let prob = old_prob * ip;
                if value < self.getValueBegin() || self.base.getValueEnd() <= value {
                    self.d_probLost += prob;
                } else {
                    if value < valueLower {
                        valueLower = value;
                    }
                    if valueUpper < value {
                        valueUpper = value;
                    }
                    let pos = self.base.getArrayPos(value) as usize;
                    self.base.lgetArray()[new_idx][pos] += prob;
                }
            }
        }

        *self.base.lgetValueLower() = valueLower;
        *self.base.lgetValueUpper() = valueUpper + 1;
        *self.base.lgetStep() += 1;
    }

    pub fn clear(&mut self, valueLower_: i64, valueUpper_: i64, prob_: Option<&[f64]>) {
        self.base.clear(valueLower_, valueUpper_, prob_);
        self.d_probLost = 0.0;
    }

    pub fn clear_default(&mut self) {
        self.clear(0, 0, None);
    }

    pub fn setValueFct(&mut self, valueFct_: Option<ValueFct>) {
        self.base.setValueFct(valueFct_);
    }

    pub fn setInput(&mut self, dimInputProb_: usize, inputProb_: Option<&[f64]>) {
        self.base.setInput(dimInputProb_, inputProb_);
    }

    pub fn getProb(&self, value_: i64) -> f64 {
        self.base.getProb(value_)
    }

    pub fn getStep(&self) -> usize {
        self.base.getStep()
    }

    pub fn getArray(&self) -> &[Vec<f64>; 2] {
        self.base.getArray()
    }

    pub fn getArrayCapacity(&self) -> usize {
        self.base.getArrayCapacity()
    }

    pub fn getValueBegin(&self) -> i64 {
        self.base.getValueBegin()
    }

    pub fn getValueLower(&self) -> i64 {
        self.base.getValueLower()
    }

    pub fn getValueUpper(&self) -> i64 {
        self.base.getValueUpper()
    }

    pub fn getValueFct(&self) -> Option<ValueFct> {
        self.base.getValueFct()
    }

    pub fn getDimInputProb(&self) -> usize {
        self.base.getDimInputProb()
    }

    pub fn getInputProb(&self) -> &[f64] {
        self.base.getInputProb()
    }

    pub fn getProbLost(&self) -> f64 {
        self.d_probLost
    }

    fn reserve(&mut self, arrayCapacity_: usize) {
        if arrayCapacity_ == self.getArrayCapacity() {
            return;
        }
        if self.getArrayCapacity() < arrayCapacity_ {
            self.base.reserve(arrayCapacity_);
            return;
        }
        for i in 0..2 {
            self.base.lgetArray()[i].truncate(arrayCapacity_);
        }
        *self.base.lgetArrayCapacity() = arrayCapacity_;
    }

    fn setValueBegin(&mut self, valueBegin_: i64) {
        if valueBegin_ <= self.getValueBegin() {
            self.base.setValueBegin(valueBegin_);
            return;
        }
        let offSet = (valueBegin_ - self.getValueBegin()) as usize;
        for i in 0..2 {
            let old = self.base.lgetArray()[i].clone();
            self.base.lgetArray()[i].fill(0.0);
            if offSet < self.getArrayCapacity() {
                let len = self.getArrayCapacity() - offSet;
                self.base.lgetArray()[i][..len].copy_from_slice(&old[offSet..offSet + len]);
            }
        }
        *self.base.lgetValueBegin() = valueBegin_;
    }
}

impl crate::stats::alp_dynprogprobproto::DynProgProbProto for DynProgProbLim {
    fn bool_(&self) -> bool {
        self.base.bool_()
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

#[cfg(test)]
mod tests {
    use super::*;

    fn add_input(old: i64, state: usize) -> i64 {
        old + state as i64
    }

    #[test]
    fn test_dyn_prog_prob_lim_loses_out_of_range_probability() {
        let mut dp = DynProgProbLim::new(Some(add_input), 2, Some(&[0.5, 0.5]), -1, 2, None);
        dp.setLimits(0, 2);
        dp.update();
        assert_eq!(dp.getProb(0), 0.5);
        assert_eq!(dp.getProb(1), 0.5);
        dp.update();
        assert_eq!(dp.getProb(0), 0.25);
        assert_eq!(dp.getProb(1), 0.5);
        assert!((dp.getProbLost() - 0.25).abs() < 1e-15);
    }
}
