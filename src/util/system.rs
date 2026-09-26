//! Platform system helpers.

pub mod get_rss;
pub use get_rss::{get_current_rss, get_peak_rss};

pub mod system;
pub use system::*;
