//! Hierarchy-preserving facade for
//! `diamond/src/lib/alp/njn_dynprogprobproto.cpp`.
//!
//! The implementation file contains only the abstract base class's empty
//! virtual destructor. Rust trait objects already drop their concrete value
//! through the vtable, so no explicit destructor function is required.

pub use crate::stats::alp_dynprogprob::ValueFct;
pub use crate::stats::alp_dynprogprobproto::DynProgProbProto;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::alp_dynprogprob::DynProgProb;
    use std::cell::Cell;
    use std::rc::Rc;

    struct DropProbe {
        dropped: Rc<Cell<bool>>,
    }

    impl DynProgProbProto for DropProbe {
        fn bool_(&self) -> bool {
            true
        }
        fn clear_default(&mut self) {}
        fn update(&mut self) {}
        fn getProb(&self, _value: i64) -> f64 {
            0.0
        }
        fn getStep(&self) -> usize {
            0
        }
        fn getValueLower(&self) -> i64 {
            0
        }
        fn getValueUpper(&self) -> i64 {
            0
        }
    }

    impl Drop for DropProbe {
        fn drop(&mut self) {
            self.dropped.set(true);
        }
    }

    #[test]
    fn trait_object_dispatches_to_existing_translation() {
        fn add_input(old: i64, state: usize) -> i64 {
            old + state as i64
        }
        let mut value = DynProgProb::new(Some(add_input), 2, Some(&[0.25, 0.75]), 0, 0, None);
        let object: &mut dyn DynProgProbProto = &mut value;
        object.update();
        assert_eq!(object.getProb(0), 0.25);
        assert_eq!(object.getProb(1), 0.75);
    }

    #[test]
    fn boxed_trait_object_runs_concrete_destructor() {
        let dropped = Rc::new(Cell::new(false));
        {
            let _object: Box<dyn DynProgProbProto> = Box::new(DropProbe {
                dropped: Rc::clone(&dropped),
            });
        }
        assert!(dropped.get());
    }

    #[test]
    fn upstream_translation_unit_remains_destructor_only() {
        let source = include_str!("../../../diamond/src/lib/alp/njn_dynprogprobproto.cpp");
        let executable = source
            .split("using namespace Njn;")
            .nth(1)
            .expect("namespace declaration")
            .trim();
        assert_eq!(
            executable, "DynProgProbProto::~DynProgProbProto () {}",
            "audit this facade if executable code is added upstream"
        );
    }
}
