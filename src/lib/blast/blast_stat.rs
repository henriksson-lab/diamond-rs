//! Inventory marker for `diamond/src/lib/blast/blast_stat.cpp`.
//!
//! The audited upstream translation unit contains only a three-line creation
//! comment: it defines no functions, tables, constants, or static
//! initialization. The large API declared by `blast_stat.h` is implemented in
//! other translation units—principally `blastn_score.cpp`, with a separate
//! Karlin implementation in `stats/comp_based_stats.cpp`—and belongs to those
//! files' audits. Keeping this physical module prevents the empty vendor file
//! from being mistaken for an unaudited numerical implementation.
