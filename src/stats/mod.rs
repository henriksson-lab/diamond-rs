pub mod alp;
pub mod alp_approx;
pub mod alp_doubletype;
pub mod alp_dynprogprob;
pub mod alp_dynprogproblim;
pub mod alp_dynprogprobproto;
pub mod alp_function;
pub mod alp_integer;
pub mod alp_ioutil;
pub mod alp_localmaxstat;
pub mod alp_localmaxstat_matrix;
pub mod alp_localmaxstat_util;
pub mod alp_matrix;
pub mod alp_memutil;
pub mod alp_random;
pub mod alp_root;
pub mod alp_uniform;
pub mod alp_vector;
pub mod cbs;
pub mod linear_algebra;
pub mod matrices;
pub mod pvalues;
pub mod score_matrix;
pub mod sls_alignment_evaluer;
pub mod sls_alp;
pub mod sls_alp_data;
pub mod sls_alp_regression;
pub mod sls_alp_sim;
pub mod sls_basic;
pub mod standard_matrix;
pub mod target_freq;

pub fn approx_id(raw_score: i32, range1: i32, range2: i32) -> f64 {
    let m = range1.max(range2);
    if m == 0 {
        100.0
    } else {
        (raw_score as f64 / m as f64 * 16.56 + 11.41).clamp(0.0, 100.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_approx_id_bounds() {
        assert_eq!(approx_id(10, 0, 0), 100.0);
        assert_eq!(approx_id(-100, 10, 10), 0.0);
        assert_eq!(approx_id(1000, 10, 10), 100.0);
        assert!((approx_id(30, 10, 20) - 36.25).abs() < 1e-9);
    }
}
