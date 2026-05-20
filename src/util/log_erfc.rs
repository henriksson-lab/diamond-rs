//! High-precision `erfc(x)` and `log(erfc(x))`, faithfully ported from GSL's
//! `gsl_sf_erfc_e` / `gsl_sf_log_erfc_e` in `diamond/src/lib/gsl/gsl.h:296-407`.
//!
//! The Chebyshev coefficient tables (`erfc_xlt1_data`, `erfc_x15_data`,
//! `erfc_x510_data`) and `erfc8_sum` rational approximation are copied verbatim
//! from the C++ source. These are public-domain numerical recipes (GSL is
//! GPL-licensed software, but the polynomial-fit data is not copyrightable);
//! the original GSL implementation by G. Jungman is the reference.
//!
//! Used by `stats::pvalues::log_area` for the log-space E-value computation,
//! where direct `erfc(x)` evaluated via Abramowitz & Stegun's ~1e-7
//! approximation produces a 1-ULP drift at the `%.2e` rounding boundary.

const GSL_DBL_EPSILON: f64 = 2.220_446_049_250_313_1e-16;
const GSL_ROOT6_DBL_EPSILON: f64 = 2.460_783_300_575_925_1e-3;
const M_SQRTPI: f64 = 1.772_453_850_905_516;

struct ChebSeries {
    c: &'static [f64],
    order: usize,
    a: f64,
    b: f64,
}

/// Evaluate a Chebyshev series at `x`, matching GSL `cheb_eval_e`.
/// Returns the value only; the C++ `result.err` is not used by callers in this codebase.
fn cheb_eval(cs: &ChebSeries, x: f64) -> f64 {
    let mut d = 0.0f64;
    let mut dd = 0.0f64;
    let y = (2.0 * x - cs.a - cs.b) / (cs.b - cs.a);
    let y2 = 2.0 * y;
    for j in (1..=cs.order).rev() {
        let temp = d;
        d = y2 * d - dd + cs.c[j];
        dd = temp;
    }
    y * d - dd + 0.5 * cs.c[0]
}

// Chebyshev fit for erfc((t+1)/2), -1 < t < 1
static ERFC_XLT1_DATA: [f64; 20] = [
    1.060_734_164_217_699_8,
    -0.425_824_458_043_810_4,
    0.049_552_626_796_204_34,
    0.004_492_934_887_683_827_5,
    -0.001_291_941_046_584_969_5,
    -0.000_018_363_892_921_493_963,
    0.000_022_111_147_040_995_263,
    -5.233_374_852_342_572e-7,
    -2.781_847_888_335_379e-7,
    1.411_580_927_488_131e-8,
    2.725_712_963_305_617e-9,
    -2.063_439_048_720_706_4e-10,
    -2.142_739_919_967_854e-11,
    2.229_902_555_393_582e-12,
    1.362_500_746_506_983e-13,
    -1.951_440_109_222_931e-14,
    -6.856_271_692_317_046e-16,
    1.445_064_928_696_999_3e-16,
    2.459_353_064_605_365e-18,
    -9.295_995_612_205_234e-19,
];
const ERFC_XLT1_CS: ChebSeries = ChebSeries {
    c: &ERFC_XLT1_DATA,
    order: 19,
    a: -1.0,
    b: 1.0,
};

// Chebyshev fit for erfc(x) * exp(x^2), 1 < x < 5, x = 2t + 3
static ERFC_X15_DATA: [f64; 25] = [
    0.440_458_320_243_381_1,
    -0.143_958_836_762_168_34,
    0.044_786_499_817_939_265,
    -0.013_343_124_200_271_212,
    0.003_824_682_739_750_47,
    -0.001_058_699_227_195_126_5,
    0.000_283_859_419_210_073_75,
    -0.000_073_906_170_662_206_76,
    0.000_018_725_312_521_489_18,
    -4.625_309_811_649_194e-6,
    1.115_586_572_444_328_6e-6,
    -2.630_986_626_508_341_3e-7,
    6.074_621_227_245_518e-8,
    -1.374_608_655_398_654_5e-8,
    3.051_570_519_054_751_5e-9,
    -6.651_747_897_203_108e-10,
    1.424_833_462_732_078e-10,
    -3.001_411_273_953_239e-11,
    6.221_717_926_453_481e-12,
    -1.269_946_392_256_685e-12,
    2.553_858_830_332_576e-13,
    -5.062_582_375_070_387e-14,
    9.897_054_094_783_273e-15,
    -1.906_859_787_891_921_8e-15,
    3.508_266_480_327_378_5e-16,
];
const ERFC_X15_CS: ChebSeries = ChebSeries {
    c: &ERFC_X15_DATA,
    order: 24,
    a: -1.0,
    b: 1.0,
};

// Chebyshev fit for erfc(x) * x * exp(x^2), 5 < x < 10, x = (5t + 15)/2
static ERFC_X510_DATA: [f64; 20] = [
    1.116_849_901_235_457,
    0.003_736_240_359_381_998_4,
    -0.000_916_623_948_045_470_2,
    0.000_199_094_325_044_940_83,
    -0.000_040_276_384_918_650_07,
    7.765_152_646_970_61e-6,
    -1.444_647_942_066_890_7e-6,
    2.613_119_303_434_64e-7,
    -4.618_330_266_348_441_5e-8,
    8.002_531_115_129_437e-9,
    -1.362_911_148_627_930_3e-9,
    2.285_704_830_901_608_7e-10,
    -3.780_225_215_632_518e-11,
    6.172_536_838_745_283e-12,
    -9.960_192_909_553_169e-13,
    1.589_531_437_069_807_8e-13,
    -2.510_459_710_471_625e-14,
    3.926_078_289_891_258e-15,
    -6.079_706_193_841_603_5e-16,
    9.126_006_072_647_948e-17,
];
const ERFC_X510_CS: ChebSeries = ChebSeries {
    c: &ERFC_X510_DATA,
    order: 19,
    a: -1.0,
    b: 1.0,
};

/// Estimate `erfc(x)` * `exp(x*x)` for 8 < x < 100 (Hart et al rational approx).
fn erfc8_sum(x: f64) -> f64 {
    let p: [f64; 6] = [
        2.978_865_626_393_992_8,
        7.409_740_605_964_742,
        6.160_209_853_109_630,
        5.019_049_726_784_267,
        1.275_366_644_729_966,
        0.564_189_583_547_755_1,
    ];
    let q: [f64; 7] = [
        3.369_075_206_982_753,
        9.608_965_327_192_788,
        17.081_440_747_466_006,
        12.048_951_927_855_129,
        9.396_034_016_235_054,
        2.260_528_520_767_327,
        1.0,
    ];
    let mut num = p[5];
    for i in (0..=4).rev() {
        num = x * num + p[i];
    }
    let mut den = q[6];
    for i in (0..=5).rev() {
        den = x * den + q[i];
    }
    num / den
}

fn log_erfc8(x: f64) -> f64 {
    erfc8_sum(x).ln() - x * x
}

/// `erfc(x)` matching GSL `gsl_sf_erfc_e`.
pub fn erfc(x: f64) -> f64 {
    let ax = x.abs();
    let e_val = if ax <= 1.0 {
        let t = 2.0 * ax - 1.0;
        cheb_eval(&ERFC_XLT1_CS, t)
    } else if ax <= 5.0 {
        let ex2 = (-x * x).exp();
        let t = 0.5 * (ax - 3.0);
        ex2 * cheb_eval(&ERFC_X15_CS, t)
    } else if ax < 10.0 {
        let exterm = (-x * x).exp() / ax;
        let t = (2.0 * ax - 15.0) / 5.0;
        exterm * cheb_eval(&ERFC_X510_CS, t)
    } else {
        erfc8_sum(ax) * (-ax * ax).exp()
    };
    if x < 0.0 {
        2.0 - e_val
    } else {
        e_val
    }
}

/// `log(erfc(x))` matching GSL `gsl_sf_log_erfc_e`. Numerically stable for large x.
pub fn log_erfc(x: f64) -> f64 {
    if x * x < 10.0 * GSL_ROOT6_DBL_EPSILON {
        // Taylor series of -1/2 log(erfc(sqrt(pi) * y)) around y = 0
        let y = x / M_SQRTPI;
        use std::f64::consts::PI;
        let c3 = (4.0 - PI) / 3.0;
        let c4 = 2.0 * (1.0 - PI / 3.0);
        let c5 = -0.001_829_764_677_455_021;
        let c6 = 0.026_296_515_210_574_65;
        let c7 = -0.016_215_753_788_354_04;
        let c8 = 0.001_259_939_617_621_16;
        let c9 = 0.005_569_646_491_38;
        let c10 = -0.004_556_333_980_2;
        let c11 = 0.000_946_158_903_2;
        let c12 = 0.001_320_024_317_4;
        let c13 = -0.001_429_06;
        let c14 = 0.000_482_04;
        let series = c8 + y * (c9 + y * (c10 + y * (c11 + y * (c12 + y * (c13 + c14 * y)))));
        let series = y
            * (1.0 + y * (1.0 + y * (c3 + y * (c4 + y * (c5 + y * (c6 + y * (c7 + y * series)))))));
        -2.0 * series
    } else if x > 8.0 {
        log_erfc8(x)
    } else {
        erfc(x).ln()
    }
}

// `GSL_DBL_EPSILON` is part of the GSL `result.err` tracking we don't carry over.
// Keep the constant declared so the port is auditable against the C++ source.
#[allow(dead_code)]
const _UNUSED: f64 = GSL_DBL_EPSILON;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_erfc_basic() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-15);
        assert!((erfc(1.0) - 0.157_299_207_050_285_13).abs() < 1e-13);
        // erfc(2) ≈ 0.0046777349810472667
        assert!((erfc(2.0) - 0.004_677_734_981_047_267).abs() < 1e-14);
        // erfc(-1) = 2 - erfc(1)
        assert!((erfc(-1.0) - (2.0 - 0.157_299_207_050_285_13)).abs() < 1e-13);
    }

    #[test]
    fn test_log_erfc_large() {
        // log(erfc(10)) ≈ -102.879889 — direct ln(erfc(10)) underflows since
        // erfc(10) ≈ 2e-45 is still representable, but for larger x it goes to 0.
        let v = log_erfc(10.0);
        assert!((v - (-102.879_889_024_844_89)).abs() < 1e-10);
        // log(erfc(20)) ≈ -403.569343 — direct ln would be -inf here.
        let v = log_erfc(20.0);
        assert!((v - (-403.569_343_334_104_25)).abs() < 1e-8);
    }

    #[test]
    fn test_log_erfc_small() {
        // For small |x|, uses the series. log(erfc(0)) = log(1) = 0
        assert!(log_erfc(0.0).abs() < 1e-14);
        // log(erfc(0.01)) ≈ -0.01128...
        assert!((log_erfc(0.01) - erfc(0.01).ln()).abs() < 1e-13);
    }

    #[test]
    fn test_log_erfc_negative() {
        // log(erfc(-1)) = log(2 - erfc(1))
        let v = log_erfc(-1.0);
        let expected = (2.0 - 0.157_299_207_050_285_13_f64).ln();
        assert!((v - expected).abs() < 1e-13);
    }
}
