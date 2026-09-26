//! Runtime search configuration.

pub mod config;
pub mod double_indexed;
pub mod main;
pub mod tools;

pub use config::{Algo, Config as SearchConfig, ConfigInput as SearchConfigInput, Round};
