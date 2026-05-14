pub use crate::stats::score_matrix::{CutoffTable, CutoffTable2D};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::score_matrix::ScoreMatrix;

    #[test]
    fn test_scores_cutoff_table_reexports() {
        let sm = ScoreMatrix::new("blosum62", 11, 1, 0, 1, 0).unwrap();
        let table = CutoffTable::new(&sm, 1e-3);
        let table2 = CutoffTable2D::new(&sm, 1e-3);
        assert_eq!(table.call(64), table.get(64));
        assert_eq!(table2.call(64, 256), table2.get(64, 256));
    }
}
