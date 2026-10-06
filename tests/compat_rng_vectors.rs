//! Platform reference vectors for DIAMOND's C-runtime RNG compatibility path.
//!
//! The actual assertions run in a spawned copy of this test binary so their
//! calls to process-global `srand`/`rand` cannot perturb another test. Vectors
//! are hardcoded only for GNU libc. Linux musl is verified live against its
//! target C++ runtime by the separate oracle; macOS libc and the Windows CRT
//! remain deliberately unverified rather than guessed.

#[path = "../src/util/compat_rng.rs"]
#[cfg(all(target_os = "linux", target_env = "gnu"))]
mod compat_rng;

use std::process::Command;

#[cfg(all(target_os = "linux", target_env = "gnu"))]
use diamond::basic::value::AMINO_ACID_COUNT;
#[cfg(all(target_os = "linux", target_env = "gnu"))]
use diamond::data::fasta::read_fasta_amino_acid;
#[cfg(all(target_os = "linux", target_env = "gnu"))]
use diamond::masking::tantan::TantanMasker;
#[cfg(all(target_os = "linux", target_env = "gnu"))]
use diamond::stats::matrices::BLOSUM62;
#[cfg(all(target_os = "linux", target_env = "gnu"))]
use diamond::stats::pvalues::pvalues;
#[cfg(all(target_os = "linux", target_env = "gnu"))]
use diamond::stats::sls_alignment_evaluer::AlignmentEvaluer;
#[cfg(all(target_os = "linux", target_env = "gnu"))]
use diamond::stats::sls_basic::AlignmentEvaluerParametersWithErrors;
#[cfg(all(target_os = "linux", target_env = "gnu"))]
use diamond::tantan::LambdaCalculator;

const CHILD_ENV: &str = "DIAMOND_COMPAT_RNG_VECTOR_CHILD";
const CHILD_TEST: &str = "compat_rng_reference_vector_child";

#[test]
#[cfg_attr(
    not(all(target_os = "linux", target_env = "gnu")),
    ignore = "no authoritative C rand vectors captured for this target"
)]
fn compat_rng_reference_vectors_in_isolated_subprocess() {
    for repetition in 0..3 {
        let output = Command::new(std::env::current_exe().unwrap())
            .args(["--ignored", "--exact", CHILD_TEST, "--test-threads=1"])
            .env(CHILD_ENV, "1")
            .output()
            .unwrap();

        assert!(
            output.status.success(),
            "isolated RNG vector check failed on repetition {repetition}\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr),
        );
    }
}

#[test]
#[ignore = "subprocess-only helper invoked by the reference-vector test"]
fn compat_rng_reference_vector_child() {
    assert_eq!(std::env::var_os(CHILD_ENV).as_deref(), Some("1".as_ref()));

    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    {
        // Verified against the GNU libc implementation on the Linux CI host.
        assert_vector(
            1,
            [
                1_804_289_383,
                846_930_886,
                1_681_692_777,
                1_714_636_915,
                1_957_747_793,
                424_238_335,
                719_885_386,
                1_649_760_492,
            ],
        );
        assert_vector(
            12_345,
            [
                383_100_999,
                858_300_821,
                357_768_173,
                455_528_251,
                133_005_921,
                116_285_904,
                591_987_137,
                102_557_902,
            ],
        );

        // Keep downstream snapshots in this same subprocess: all three paths
        // intentionally share the process-global C runtime RNG state. The
        // standard-normal and TANTAN values were captured independently from
        // the pinned C++ sources on this GNU/Linux target.
        compat_rng::c_srand(12_345);
        let normals = std::array::from_fn::<_, 8, _>(|_| pvalues::standard_normal());
        assert_f64_snapshot(
            &normals,
            &[
                13_832_805_915_482_062_431,
                4_601_712_150_747_084_250,
                4_612_189_209_984_272_505,
                4_609_585_598_460_568_136,
                13_820_271_598_230_441_230,
                4_597_971_959_140_252_177,
                4_602_858_202_703_976_877,
                13_832_311_949_733_793_287,
            ],
            2.0e-14,
            "standard_normal",
        );

        compat_rng::c_srand(1);
        let mut lambda = LambdaCalculator::new();
        lambda.calculate_flat_i8(&BLOSUM62.scores, AMINO_ACID_COUNT, 20);
        assert_f64_snapshot(
            &[lambda.lambda()],
            &[4_599_508_864_510_404_738],
            2.0e-14,
            "TANTAN lambda",
        );
        let probability_bits = [
            4_590_311_920_184_368_722,
            4_589_536_674_147_392_867,
            4_585_777_739_314_088_530,
            4_587_890_597_674_949_845,
            4_582_692_031_376_121_149,
            4_586_294_421_086_378_088,
            4_586_993_740_036_516_179,
            4_589_712_300_121_341_509,
        ];
        assert_f64_snapshot(
            &lambda.letter_probs1().unwrap()[..8],
            &probability_bits,
            2.0e-14,
            "TANTAN row probabilities",
        );
        assert_f64_snapshot(
            &lambda.letter_probs2().unwrap()[..8],
            &probability_bits,
            2.0e-14,
            "TANTAN column probabilities",
        );

        let mut evaluator = AlignmentEvaluer::new();
        evaluator
            .init_parameters_with_errors(&parameters_with_errors())
            .unwrap();
        assert_f64_snapshot(
            &evaluator.parameters().m_LambdaSbs[..8],
            &[
                4_598_360_626_142_876_594,
                4_598_589_354_046_589_136,
                4_598_536_399_833_527_353,
                4_598_434_270_926_483_113,
                4_598_417_383_537_481_528,
                4_598_406_022_183_420_933,
                4_598_545_512_816_164_735,
                4_598_554_019_085_063_951,
            ],
            2.0e-14,
            "initParametersWithErrors lambda samples",
        );
        assert_f64_snapshot(
            &evaluator.parameters().m_TauSbs[..8],
            &[
                4_620_718_581_650_682_260,
                4_620_615_237_306_002_225,
                4_620_601_757_903_572_318,
                4_620_681_998_564_254_315,
                4_620_564_055_266_977_452,
                4_620_720_593_140_562_828,
                4_620_615_240_279_182_955,
                4_620_723_563_896_496_957,
            ],
            2.0e-13,
            "initParametersWithErrors tau samples",
        );

        // The production masking shape computes the RNG-dependent likelihood
        // matrix once, before dispatching immutable masker state to workers.
        // Exercise that boundary rather than racing C's process-global RNG.
        compat_rng::c_srand(1);
        let masker = std::sync::Arc::new(TantanMasker::new(&BLOSUM62, 0.9));
        let input = read_fasta_amino_acid(
            b">repeat\nPPTPPTPPTPPTPPTPPTPPTPPTPPTPPTPPTPPTPPTPPT\n".as_slice(),
        )
        .unwrap()
        .remove(0)
        .sequence;
        let mut expected = input.clone();
        masker.mask(&mut expected);
        let handles = (0..4)
            .map(|_| {
                let masker = masker.clone();
                let input = input.clone();
                std::thread::spawn(move || {
                    let mut actual = input;
                    masker.mask(&mut actual);
                    actual
                })
            })
            .collect::<Vec<_>>();
        for handle in handles {
            assert_eq!(handle.join().unwrap(), expected);
        }
    }

    #[cfg(not(all(target_os = "linux", target_env = "gnu")))]
    panic!("no authoritative C rand vectors captured for this target");
}

#[cfg(all(target_os = "linux", target_env = "gnu"))]
fn assert_f64_snapshot(actual: &[f64], expected_bits: &[u64], tolerance: f64, label: &str) {
    assert_eq!(actual.len(), expected_bits.len(), "{label} length");
    for (index, (&actual, &bits)) in actual.iter().zip(expected_bits).enumerate() {
        let expected = f64::from_bits(bits);
        assert!(
            (actual - expected).abs() <= tolerance,
            "{label}[{index}] changed: actual={actual:?}, expected={expected:?}"
        );
    }
}

#[cfg(all(target_os = "linux", target_env = "gnu"))]
fn parameters_with_errors() -> AlignmentEvaluerParametersWithErrors {
    AlignmentEvaluerParametersWithErrors {
        d_lambda: 0.267,
        d_lambda_error: 0.001,
        d_k: 0.041,
        d_k_error: 0.001,
        d_a1: 1.9,
        d_a1_error: 0.01,
        d_b1: 4.0,
        d_b1_error: 0.01,
        d_a2: 2.1,
        d_a2_error: 0.01,
        d_b2: 5.0,
        d_b2_error: 0.01,
        d_alpha1: 1.7,
        d_alpha1_error: 0.01,
        d_beta1: 6.0,
        d_beta1_error: 0.01,
        d_alpha2: 1.8,
        d_alpha2_error: 0.01,
        d_beta2: 7.0,
        d_beta2_error: 0.01,
        d_sigma: 43.0,
        d_sigma_error: 0.01,
        d_tau: 8.0,
        d_tau_error: 0.01,
    }
}

#[cfg(all(target_os = "linux", target_env = "gnu"))]
fn assert_vector<const N: usize>(seed: u32, expected: [i32; N]) {
    compat_rng::c_srand(seed);
    let actual = std::array::from_fn(|_| compat_rng::c_rand());
    assert_eq!(actual, expected, "C rand vector changed for seed {seed}");
}
