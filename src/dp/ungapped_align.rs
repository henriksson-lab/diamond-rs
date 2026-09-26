//! Compatibility facade mirroring the original `dp/ungapped_align.cpp` file.

pub use super::ungapped::{
    make_clipped_anchor, make_null_anchor, score_range, score_range_s, self_score,
    self_score_with_cbs, trivial, trivial_at, ungapped_window, xdrop_anchored_left,
    xdrop_anchored_right, xdrop_ungapped, xdrop_ungapped_from_anchor, xdrop_ungapped_right,
    xdrop_ungapped_with_cbs, xdrop_ungapped_with_float_cbs, xdrop_ungapped_with_identities,
    DiagonalSegment,
};
