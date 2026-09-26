//! Inventory marker for `diamond/src/lib/blast/blast_stat.cpp`.
//!
//! The audited upstream translation unit contains only a three-line creation
//! comment: it defines no functions, tables, constants, or static
//! initialization. The large API declared by `blast_stat.h` is implemented in
//! other translation units—principally `blastn_score.cpp`, with a separate
//! Karlin implementation in `stats/comp_based_stats.cpp`—and belongs to those
//! files' audits. Keeping this physical module prevents the empty vendor file
//! from being mistaken for an unaudited numerical implementation.

#[cfg(test)]
mod tests {
    #[test]
    fn upstream_translation_unit_remains_comment_only() {
        let upstream = include_str!("../../../diamond/src/lib/blast/blast_stat.cpp");
        let executable_lines = upstream
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty() && !line.starts_with("//"))
            .collect::<Vec<_>>();
        assert!(
            executable_lines.is_empty(),
            "blast_stat.cpp gained executable content: {executable_lines:?}"
        );
    }
}
