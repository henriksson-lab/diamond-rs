//! Original-object TANTAN scale/probability fixtures, without runtime C++ FFI.
//!
//! The original GNU/Linux object uses FMA in LU decomposition and triangular
//! solves. Full bit patterns and RNG continuation are retained deliberately:
//! platform math differences must be investigated rather than hidden by a
//! loose lambda tolerance.
use diamond::tantan::LambdaCalculator;
use diamond::util::compat_rng;

fn bits_list(value: &str) -> Option<Vec<u64>> {
    (value != "-").then(|| {
        value
            .split(',')
            .map(|v| u64::from_str_radix(v, 16).unwrap())
            .collect()
    })
}

#[test]
fn scale_probabilities_state_and_rng_match_original_object() {
    // Keep every RNG-consuming fixture in this one test. The compatibility RNG
    // is process-global; separate parallel tests must not reseed it mid-call.
    let inputs = include_str!("data/tantan_lambda/matrices.tsv");
    let expected = include_str!("data/tantan_lambda/states.tsv");
    let mut calculator = LambdaCalculator::new();
    assert!(calculator.is_bad());
    assert_eq!(calculator.lambda().to_bits(), (-1.0f64).to_bits());
    assert!(calculator.letter_probs1().is_none());
    assert!(calculator.letter_probs2().is_none());
    let mut count = 0;
    for (input, state) in inputs.lines().zip(expected.lines()) {
        let columns: Vec<_> = input.split('\t').collect();
        let result: Vec<_> = state.split('\t').collect();
        assert_eq!(columns[0], result[0]);
        let label = columns[0];
        let seed: i64 = columns[1].parse().unwrap();
        let size: usize = columns[2].parse().unwrap();
        if seed >= 0 {
            compat_rng::seed(seed as u32);
        }
        let scores: Vec<i32> = columns[3]
            .split(',')
            .filter(|v| !v.is_empty())
            .map(|v| v.parse().unwrap())
            .collect();
        assert_eq!(scores.len(), size * size);
        let matrix: Vec<Vec<i32>> = if size == 0 {
            Vec::new()
        } else {
            scores.chunks(size).map(<[i32]>::to_vec).collect()
        };
        calculator.calculate(&matrix);
        assert_eq!(calculator.is_bad(), result[1] == "1", "{label}: state");
        assert_eq!(
            calculator.lambda().to_bits(),
            u64::from_str_radix(result[2], 16).unwrap(),
            "{label}: lambda bits"
        );
        for (index, probabilities) in [calculator.letter_probs1(), calculator.letter_probs2()]
            .into_iter()
            .enumerate()
        {
            assert_eq!(
                probabilities.map(|p| p.iter().map(|v| v.to_bits()).collect::<Vec<_>>()),
                bits_list(result[index + 3]),
                "{label}: probabilities {}",
                index + 1
            );
        }
        assert_eq!(
            compat_rng::next_i31(),
            result[5].parse::<i32>().unwrap(),
            "{label}: RNG continuation"
        );
        // Exercise reset after each completed calculation, then allow the next
        // calculate call to reuse the same object as the source driver does.
        calculator.set_bad();
        assert!(calculator.is_bad());
        assert_eq!(calculator.lambda().to_bits(), (-1.0f64).to_bits());
        assert!(calculator.letter_probs1().is_none());
        assert!(calculator.letter_probs2().is_none());
        count += 1;
    }
    assert_eq!(inputs.lines().count(), expected.lines().count());
    assert_eq!(count, 36);
}
