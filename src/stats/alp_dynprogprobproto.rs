#![allow(non_snake_case)]

pub trait DynProgProbProto {
    fn bool_(&self) -> bool;
    fn clear_default(&mut self);
    fn update(&mut self);
    fn getProb(&self, value_: i64) -> f64;
    fn getStep(&self) -> usize;
    fn getValueLower(&self) -> i64;
    fn getValueUpper(&self) -> i64;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::alp_dynprogprob::DynProgProb;

    fn add_input(old: i64, state: usize) -> i64 {
        old + state as i64
    }

    #[test]
    fn test_dyn_prog_prob_proto_trait_object() {
        let mut dp = DynProgProb::new(Some(add_input), 2, Some(&[0.25, 0.75]), 0, 0, None);
        let proto: &mut dyn DynProgProbProto = &mut dp;
        assert!(proto.bool_());
        assert_eq!(proto.getStep(), 0);
        assert_eq!(proto.getValueLower(), 0);
        assert_eq!(proto.getValueUpper(), 1);
        proto.update();
        assert_eq!(proto.getStep(), 1);
        assert_eq!(proto.getProb(0), 0.25);
        assert_eq!(proto.getProb(1), 0.75);
        proto.clear_default();
        assert_eq!(proto.getStep(), 0);
        assert_eq!(proto.getProb(0), 1.0);
    }
}
