//! Frozen GNU11/Linux x86_64 Eigen vectors; no C++ runtime dependency.
//! Three libm/NaN-sensitive cases additionally have independent mathematical
//! checks. Their captured bits are asserted on the captured target profile;
//! the other 51 cases retain exact comparisons on every target.
use diamond::mcl::mcl::{DenseMatrix, Mcl, SparseMatrix, Triplet};

#[test]
fn source_mcl_numerical_stages() {
    let mut mismatches = Vec::new();
    for line in include_str!("data/mcl_numeric/source.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let f: Vec<_> = line.split('\t').collect();
        let name = f[0];
        let mode = f[1];
        let op = f[2];
        let r: f32 = f[3].parse().unwrap();
        let expansion = f[4].parse().unwrap();
        let iterations = f[5].parse().unwrap();
        let threads = f[6].parse().unwrap();
        let n: usize = f[7].parse().unwrap();
        let input: Vec<f32> = f[8].split(',').map(|v| v.parse().unwrap()).collect();
        let mut triplets = if op == "build" {
            input[1..]
                .chunks_exact(3)
                .map(|t| Triplet::new(t[0] as usize, t[1] as usize, t[2]))
                .collect()
        } else {
            input
                .iter()
                .enumerate()
                .map(|(i, &v)| Triplet::new(i % n, i / n, v))
                .collect()
        };
        let mut mcl = Mcl::default();
        let result: Result<Vec<u32>, String> = if mode == "dense" {
            let matrix = if op == "build" {
                Mcl::get_dense_matrix_and_clear(
                    &(0..n).collect::<Vec<_>>(),
                    &mut triplets,
                    r != 0.0,
                )
            } else {
                DenseMatrix::from_triplets(n, &triplets, false)
            };
            let out = match op {
                "build" => matrix,
                "gamma" => Mcl::get_gamma_dense(&matrix, r),
                "exp" => Mcl::get_exp_dense(&matrix, r),
                "exp-retain" => {
                    let mut out = DenseMatrix::from_triplets(
                        n,
                        &(0..n * n)
                            .map(|i| Triplet::new(i % n, i / n, 42.0))
                            .collect::<Vec<_>>(),
                        false,
                    );
                    Mcl::get_exp_dense_into(&matrix, &mut out, r);
                    out
                }
                "markov" => {
                    let mut out = matrix;
                    mcl.markov_process_dense(&mut out, r, expansion, iterations);
                    out
                }
                _ => panic!("unknown fixture operation"),
            };
            Ok((0..n * n).map(|i| out[(i % n, i / n)].to_bits()).collect())
        } else {
            let matrix = if op == "build" {
                Mcl::get_sparse_matrix_and_clear(
                    &(0..n).collect::<Vec<_>>(),
                    &mut triplets,
                    r != 0.0,
                )
            } else {
                // Eigen sparseView omits numerical zeros in the oracle wrapper.
                triplets.retain(|t| t.value != 0.0);
                SparseMatrix::from_triplets(n, &triplets, false)
            };
            let out = match op {
                "build" => Ok(matrix),
                "gamma" => Ok(Mcl::get_gamma_sparse(&matrix, r, threads)),
                "exp" => Mcl::get_exp_sparse(&matrix, r, threads),
                "markov" => {
                    let mut out = matrix;
                    mcl.markov_process_sparse(&mut out, r, expansion, iterations, threads)
                        .map(|()| out)
                }
                _ => panic!("unknown fixture operation"),
            };
            out.map(|out| {
                let mut values = vec![0u32; n * n];
                for t in out.triplets() {
                    values[t.col * n + t.row] = t.value.to_bits();
                }
                values
            })
        };
        if f[9] == "0" {
            let expected: Vec<u32> = f[10]
                .split(',')
                .map(|v| u32::from_str_radix(v, 16).unwrap())
                .collect();
            match result {
                Ok(actual) => {
                    assert_eq!(actual.len(), expected.len(), "{name} matrix dimensions");
                    let portable_case = matches!(
                        name,
                        "dense-markov-converged"
                            | "dense-gamma-empty-column"
                            | "dense-fractional-complex-cycle-8"
                    );
                    if portable_case {
                        check_libm_sensitive_case(name, &actual, r, n);
                    }
                    // The fixture was captured with GNU11/libm on x86_64 Linux
                    // using baseline SSE. Darwin libm differs even on x86_64;
                    // ARM also has a different default sign for 0/0 NaNs.
                    let captured_profile = cfg!(all(
                        target_os = "linux",
                        target_env = "gnu",
                        target_arch = "x86_64"
                    ));
                    for (i, (&got, &want)) in actual.iter().zip(&expected).enumerate() {
                        if got != want && (captured_profile || !portable_case) {
                            mismatches.push(format!(
                                "{name} row={} col={}: {got:08x} != Source {want:08x}",
                                i % n,
                                i / n
                            ));
                        }
                    }
                    if mcl.failed_to_converge != f[11].parse::<usize>().unwrap() {
                        mismatches.push(format!("{name} convergence counter differs"));
                    }
                }
                Err(e) => mismatches.push(format!("{name} unexpected error {e}")),
            }
        } else if result.as_ref().err().map(String::as_str) != Some(f[12]) {
            mismatches.push(format!(
                "{name} expected Source error {:?}, got {:?}",
                f[12], result
            ));
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} exact Source mismatches:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
}

fn check_libm_sensitive_case(name: &str, bits: &[u32], power: f32, n: usize) {
    let values: Vec<f64> = bits.iter().map(|&v| f32::from_bits(v) as f64).collect();
    match name {
        "dense-gamma-empty-column" => {
            // Gamma computes 0/0 in the empty column. IEEE NaN classification
            // is defined here, while its sign/payload is target arithmetic.
            assert!(values[0].is_nan() && values[1].is_nan());
            assert_eq!(bits[2], (1.0f32 / 5.0).to_bits());
            assert_eq!(bits[3], (4.0f32 / 5.0).to_bits());
        }
        "dense-markov-converged" => {
            // Both columns converge to the first absorbing state. The positive
            // residual depends on powf rounding, and is below the algorithm's
            // own f32::EPSILON stopping threshold on both captured platforms.
            assert_eq!(n, 2);
            assert_eq!(bits[0], 1.0f32.to_bits());
            assert_eq!(bits[2], 1.0f32.to_bits());
            assert_eq!(bits[1], bits[3]);
            assert!(values[1] > 0.0 && values[1] < f32::EPSILON as f64);
            assert_eq!(values[0] + values[1], 1.0);
            assert_eq!(values[2] + values[3], 1.0);
        }
        "dense-fractional-complex-cycle-8" => {
            assert_eq!(n, 8);
            // Independent real binomial series for (I + P/2)^power, where P
            // shifts a column down one position and P^8=I. Source dense.h:46
            // uses V*d.real()*V.inverse(), so its output here is the symmetric
            // half-sum of this power and its transpose: imaginary eigenvalue
            // powers are discarded before reconstruction. No eigensolver,
            // complex logarithm/power, or platform libm is used in this oracle.
            // |P/2|<1 gives absolute convergence; 200 terms leave <2^-190
            // truncation error. Compare the complete matrix, including small
            // negative entries, with an explicit n*epsilon*L1-series rounding
            // budget for the f32 matrix computation (not Source-bit parity).
            let mut column = vec![0.0f64; n];
            let mut coefficient = 1.0f64;
            let mut absolute_sum = 0.0f64;
            for k in 0..200 {
                column[k % n] += coefficient;
                absolute_sum += coefficient.abs();
                coefficient *= (power as f64 - k as f64) * 0.5 / (k + 1) as f64;
            }
            let rounding_budget = n as f64 * f32::EPSILON as f64 * absolute_sum;
            for col in 0..n {
                for row in 0..n {
                    let actual = values[col * n + row];
                    let offset = (row + n - col) % n;
                    let expected = 0.5 * (column[offset] + column[(n - offset) % n]);
                    assert!(actual.is_finite());
                    assert!(
                        (actual - expected).abs() <= rounding_budget,
                        "{name} ({row},{col}) actual={actual} binomial={expected} budget={rounding_budget}"
                    );
                }
            }
        }
        _ => unreachable!("not a profile-sensitive fixture"),
    }
}
