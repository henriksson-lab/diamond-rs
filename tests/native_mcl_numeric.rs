//! Frozen exact GNU/Eigen Source numerical vectors; no C++ runtime dependency.
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
                    for (i, (&got, &want)) in actual.iter().zip(&expected).enumerate() {
                        if got != want {
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
