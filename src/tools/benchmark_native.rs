//! Concrete execution of the active standard-build benchmark workload.
//! Fixed iteration counts and even the source's misleading labels are retained.
use std::hint::black_box;
use std::io::{self, Write};
use std::sync::Arc;
use std::time::Instant;

use crate::basic::value::Letter;
use crate::dp::anchored::{self, Target};
use crate::dp::score_profile::{make_profile16, make_profile8};
use crate::dp::swipe::anchored_wrapper::{anchored_swipe, AnchoredSwipeConfig};
use crate::dp::swipe::{self, Anchor, CarryOver, DpTarget, Flags, HspValues, Params};
use crate::stats::cbs::{CbsMode, MatrixAdjustRule, TargetMatrix, TargetMatrixAdjustment};
use crate::stats::hauser_correction::HauserCorrection;
use crate::stats::score_matrix::ScoreMatrix;
use crate::util::simd::Arch;

#[derive(Clone, Copy)]
pub(crate) struct Iterations {
    pub matrix: usize,
    pub anchored: usize,
    pub swipe: usize,
    pub diagonal: usize,
    pub banded: usize,
    pub evalue: usize,
    pub ungapped: usize,
    pub transpose: usize,
}
impl Default for Iterations {
    fn default() -> Self {
        Self {
            matrix: 10_000,
            anchored: 10_000,
            swipe: 1_000,
            diagonal: 100_000,
            banded: 10_000,
            evalue: 1_000_000,
            ungapped: 10_000_000,
            transpose: 10_000_000,
        }
    }
}

fn measure(
    out: &mut impl Write,
    label: &str,
    n: usize,
    cells: usize,
    scale: f64,
    unit: &str,
    mut f: impl FnMut(usize),
) -> io::Result<()> {
    let start = Instant::now();
    for i in 0..n {
        f(i);
    }
    let elapsed = start.elapsed().as_nanos() as f64 / (n * cells) as f64 * scale;
    writeln!(
        out,
        "{label}{} {unit}",
        crate::output::format::format_default_float(elapsed)
    )
}
fn encode(bytes: &[u8]) -> io::Result<Vec<Letter>> {
    let mut fasta = b">benchmark\n".to_vec();
    fasta.extend_from_slice(bytes);
    fasta.push(b'\n');
    Ok(crate::data::fasta::read_fasta_amino_acid(&fasta[..])?
        .remove(0)
        .sequence)
}
fn params<'a>(query: &'a [Letter], bias: &'a [i8], matrix: &'a ScoreMatrix) -> Params<'a> {
    let mut p = Params::new(query, matrix);
    p.composition_bias = Some(bias);
    p.v = HspValues::NONE;
    p.max_evalue = 0.001;
    p.cutoff_score_8bit = 240;
    p.band_bin = 24;
    p.col_bin = 400;
    p
}

/// Run the real kernels. Reduced counts exist only for internal conformance
/// tests; the command always supplies the source's full counts.
pub(crate) fn run(
    out: &mut impl Write,
    err: &mut impl Write,
    arch: Arch,
    counts: Iterations,
) -> io::Result<()> {
    let matrix = ScoreMatrix::new("blosum62", -1, -1, 0, 1, 0).map_err(io::Error::other)?;
    run_with_matrix(out, err, arch, counts, &matrix)
}

pub(crate) fn run_with_matrix(
    out: &mut impl Write,
    err: &mut impl Write,
    arch: Arch,
    counts: Iterations,
    matrix: &ScoreMatrix,
) -> io::Result<()> {
    run_with_matrix_and_controls(out, err, arch, counts, matrix, 2000, 1.0, 40)
}

pub(crate) fn run_with_matrix_and_controls(
    out: &mut impl Write,
    err: &mut impl Write,
    arch: Arch,
    counts: Iterations,
    matrix: &ScoreMatrix,
    cbs_it_limit: i32,
    transcript_len_estimate: f64,
    cbs_window: i32,
) -> io::Result<()> {
    use super::benchmark::{S1, S2, S3, S4};
    let (s1, s2, s3, s4) = (encode(S1)?, encode(S2)?, encode(S3)?, encode(S4)?);
    let avx2 = matches!(arch, Arch::Avx2 | Arch::Avx512);
    let simd = avx2 || matches!(arch, Arch::Sse4_1 | Arch::Neon);
    let mut adjusted_mode = false;

    // Both optimizer calls are commented out in the pinned source.
    let mut zeros = [0.0_f64; 400];
    let start = Instant::now();
    for _ in 0..counts.matrix {
        black_box(&mut zeros).fill(0.0);
    }
    let elapsed = start.elapsed().as_micros() as f64 / counts.matrix as f64;
    for row in zeros.chunks_exact(20) {
        for value in row {
            write!(out, "{value:.6} ")?;
        }
        writeln!(out)?;
    }
    writeln!(
        err,
        "Matrix adjust (openblas):\t\t\t{} ms",
        crate::output::format::format_default_float(elapsed)
    )?;

    if avx2 {
        let query = &s1[..128];
        let target = &s2[..128];
        let profile = make_profile16(query, None, 0, &matrix);
        // Source also constructs the unused int8 profiles/target vector.
        let profile8 = make_profile8(query, None, 0, &matrix);
        black_box((profile8.pointers(0), vec![profile8.clone(); 32]));
        let mut targets = (0..16)
            .map(|_| Target::new(target.to_vec(), -32, 32, 0, 128, 0, false))
            .collect::<Vec<_>>();
        measure(
            err,
            "Anchored Swipe (int16_t):\t",
            counts.anchored,
            128 * 64 * 16,
            1000.0,
            "ps/Cell",
            |_| {
                black_box(anchored::smith_waterman_simd_profiled(
                    query,
                    &mut targets,
                    &matrix,
                    Some(&profile),
                    Some(&profile),
                ));
                black_box(targets[0].score);
                #[cfg(test)]
                if counts.anchored == 1 {
                    assert_eq!(targets[0].score, 80, "source AVX2 anchored probe");
                }
            },
        )?;
        let mut targets = swipe::targets();
        for _ in 0..16 {
            targets[0].push_back(DpTarget::new(
                target,
                128,
                -32,
                32,
                0,
                128,
                CarryOver::default(),
                Anchor::default(),
            ));
        }
        let mut cfg = AnchoredSwipeConfig::new(query, &matrix);
        cfg.transcript_len_estimate = transcript_len_estimate;
        cfg.max_evalue = 0.001;
        measure(
            err,
            "Anchored Swipe2 (int16_t):\t",
            counts.anchored,
            128 * 64 * 16,
            1000.0,
            "ps/Cell",
            |_| {
                let result = anchored_swipe(&mut targets, &cfg);
                #[cfg(test)]
                if counts.anchored == 1 {
                    assert!(
                        result.is_empty(),
                        "source anchored pipeline probe: {result:?}"
                    );
                }
                black_box(result);
            },
        )?;
    }
    if simd {
        let query = &s3[..255];
        let bias = HauserCorrection::new(&s3, &matrix, cbs_window as u32 as usize);
        let mut statistics = crate::basic::statistics::Statistics::default();
        // Construct the actual adjusted matrix even though the source later
        // assigns it to destroyed entries in the now-empty int8 target bin.
        black_box(
            TargetMatrix::from_composition_adjustment(
                &crate::stats::cbs::compute_composition(&s3),
                s3.len() as i32,
                CbsMode::MatrixAdjust,
                &s4,
                &mut statistics,
                &matrix,
                MatrixAdjustRule::UserSpecifiedRelEntropy,
                TargetMatrixAdjustment {
                    matrix_scale: 1,
                    joint_probs: matrix.joint_probs().unwrap(),
                    background_freqs: matrix.background_freqs(),
                    freq_ratios: matrix.freq_ratios(),
                    tolerance: 0.0001,
                    max_iterations: cbs_it_limit,
                },
            )
            .map_err(io::Error::other)?,
        );
        let mut p = params(query, &bias.int8, &matrix);
        p.transcript_len_estimate = transcript_len_estimate;
        p.flags = Flags::FULL_MATRIX;
        let mut targets = swipe::targets();
        let target: Arc<[Letter]> = s4.clone().into();
        for _ in 0..32 {
            targets[0].push_back(DpTarget::new(
                target.clone(),
                s4.len() as i32,
                0,
                0,
                0,
                0,
                CarryOver::default(),
                Anchor::default(),
            ));
        }
        let cells = query.len() * s4.len() * if avx2 { 32 } else { 16 };
        for (phase, label) in [
            "SWIPE (int8_t):\t\t\t",
            "SWIPE (int16_t):\t\t",
            "SWIPE (int8_t, Stats):\t\t",
            "SWIPE (int8_t, MatrixAdjust):\t",
            "SWIPE (int8_t, CBS):\t\t",
            "SWIPE (int8_t, TB):\t\t",
        ]
        .iter()
        .enumerate()
        {
            if phase == 1 {
                targets[1] = targets[0].clone();
                targets[0].clear();
            }
            // Keep the active bin intact. Writing matrix pointers into the
            // destroyed bin-zero objects has no defined observable result.
            measure(err, label, counts.swipe, cells, 1000.0, "ps/Cell", |_| {
                let result = swipe::swipe(&targets, &mut p);
                #[cfg(test)]
                if counts.swipe == 1 {
                    assert_eq!(result.len(), 32, "source SWIPE probe");
                    assert!(
                        result.iter().all(|hsp| hsp.score == 124),
                        "source SWIPE score probe: {result:?}"
                    );
                }
                black_box(result);
            })?;
        }
        adjusted_mode = true;
        let bias = HauserCorrection::new(&s1, &matrix, cbs_window as u32 as usize);
        let profile = make_profile8(&s1, Some(&bias.int8), 0, &matrix);
        let mut scores = [0; 128];
        measure(
            err,
            "Diagonal scores:\t\t",
            counts.diagonal,
            s2.len() * 128,
            1000.0,
            "ps/Cell",
            |i| {
                crate::dp::scan_diags::scan_diags128(
                    &profile,
                    &s2,
                    -32,
                    0,
                    s2.len() as i32,
                    &mut scores,
                );
                // Source reads scores[i & 128], out of bounds for half the loop.
                // The volatile value is never output. Preserve the real scan and
                // omit that undefined read rather than fabricate a score.
                #[cfg(test)]
                if counts.diagonal == 1 {
                    assert_eq!(
                        (scores[0], scores[127]),
                        (16, 11),
                        "source AVX2 diagonal probe"
                    );
                }
                black_box(scores.get(i & 128));
            },
        )?;
    }
    let mut targets = swipe::targets();
    for _ in 0..8 {
        targets[1].push_back(DpTarget::new(
            s2.as_slice(),
            s2.len() as i32,
            -32,
            32,
            0,
            0,
            CarryOver::default(),
            Anchor::default(),
        ));
    }
    let bias = HauserCorrection::new(&s1, &matrix, cbs_window as u32 as usize);
    let mut p = params(&s1, &bias.int8, &matrix);
    p.transcript_len_estimate = transcript_len_estimate;
    // TargetIterator's constructor requires a matrix in MATRIX_ADJUST mode.
    // Source swipe() leaves that mode set; banded_swipe() has no matrices.
    if adjusted_mode
        && targets[1]
            .as_slice()
            .iter()
            .any(|target| target.matrix.is_none())
    {
        return Err(io::Error::other(
            "TargetIterator: No matrix provided for adjusted matrix.",
        ));
    }
    for (phase, label) in [
        "Banded SWIPE (int16_t, CBS):\t",
        "Banded SWIPE (int16_t):\t\t",
        "Banded SWIPE (int16_t, CBS, TB):",
    ]
    .iter()
    .enumerate()
    {
        if phase == 2 {
            p.v = HspValues::TRANSCRIPT;
        }
        measure(
            err,
            label,
            counts.banded,
            s1.len() * 65 * 16,
            1000.0,
            "ps/Cell",
            |_| {
                black_box(swipe::swipe(&targets, &mut p));
            },
        )?;
    }
    let mut x = 0.0;
    measure(err, "Evalue:\t\t\t\t", counts.evalue, 1, 1.0, "ns", |i| {
        x += matrix.evalue_norm_query_len(i as i32, 300);
        black_box(x);
    })?;
    measure(
        err,
        "Evalue (ALP):\t\t\t",
        counts.evalue,
        1,
        1.0,
        "ns",
        |_| {
            x += matrix.evalue(300, 300, 300);
            black_box(x);
        },
    )?;
    measure(
        err,
        "Scalar ungapped extension:\t",
        counts.ungapped,
        64,
        1000.0,
        "ps/Cell",
        |_| {
            black_box(crate::dp::ungapped_window::ungapped_window_score(
                &s1[34..],
                &s2[33..],
                0,
                0,
                64,
                &matrix,
            ));
        },
    )?;
    for vectorized in [false, true] {
        let input = [0_i8; 256];
        let rows = input.chunks_exact(16).collect::<Vec<_>>();
        let mut output = [0; 256];
        let label = if vectorized {
            "Transpose (16x16, vectorized):\t"
        } else {
            "Transpose (16x16, scalar):\t"
        };
        measure(
            err,
            label,
            counts.transpose,
            256,
            1000.0,
            "ps/Letter",
            |_| {
                if vectorized {
                    crate::util::simd::transpose(black_box(&rows), 16, &mut output, 16);
                } else {
                    super::benchmark::transpose_scalar::<16>(black_box(&rows), 16, &mut output);
                }
                black_box(output[0]);
            },
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn once() -> Iterations {
        Iterations {
            matrix: 1,
            anchored: 1,
            swipe: 1,
            diagonal: 1,
            banded: 1,
            evalue: 1,
            ungapped: 1,
            transpose: 1,
        }
    }
    #[test]
    fn real_simd_workloads_reach_the_source_missing_matrix_failure() {
        if !matches!(
            crate::util::simd::arch(),
            Arch::Avx2 | Arch::Avx512 | Arch::Sse4_1 | Arch::Neon
        ) {
            return;
        }
        let (mut out, mut err) = (Vec::new(), Vec::new());
        let error = run(&mut out, &mut err, crate::util::simd::arch(), once()).unwrap_err();
        assert_eq!(
            error.to_string(),
            "TargetIterator: No matrix provided for adjusted matrix."
        );
        assert_eq!(
            String::from_utf8(out)
                .unwrap()
                .split_whitespace()
                .filter(|value| *value == "0.000000")
                .count(),
            400
        );
        let err = String::from_utf8(err).unwrap();
        assert!(err.contains("SWIPE (int8_t, MatrixAdjust):"));
        assert!(err.contains("Diagonal scores:"));
        assert!(!err.contains("Banded SWIPE (int16_t, CBS):"));
    }
    #[test]
    fn generic_workloads_execute_banded_evalue_ungapped_and_transpose() {
        let (mut out, mut err) = (Vec::new(), Vec::new());
        run(&mut out, &mut err, Arch::Generic, once()).unwrap();
        let err = String::from_utf8(err).unwrap();
        for label in [
            "Banded SWIPE (int16_t, CBS, TB):",
            "Evalue (ALP):",
            "Scalar ungapped extension:",
            "Transpose (16x16, vectorized):",
        ] {
            assert!(err.contains(label), "{err}");
        }
        assert!(!err.contains("SWIPE (int8_t):"));
    }
}
