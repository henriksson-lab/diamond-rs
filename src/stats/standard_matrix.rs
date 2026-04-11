use crate::basic::value::AMINO_ACID_COUNT;

/// Statistical parameters for a scoring matrix at given gap penalties.
#[derive(Debug, Clone, Copy)]
pub struct Parameters {
    pub gap_exist: f64,
    pub gap_extend: f64,
    pub reserved: f64,
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
    pub alpha: f64,
    pub beta: f64,
    pub c: f64,
    pub alpha_v: f64,
    pub sigma: f64,
}

/// A standard amino acid scoring matrix (BLOSUM, PAM).
pub struct StandardMatrix {
    pub default_gap_open: i32,
    pub default_gap_extend: i32,
    /// Statistical parameters for various gap penalty combinations.
    /// First entry is ungapped (gap_exist and gap_extend are f64::MAX).
    pub parameters: &'static [Parameters],
    /// 26x26 scoring matrix (AMINO_ACID_COUNT x AMINO_ACID_COUNT).
    pub scores: [i8; AMINO_ACID_COUNT * AMINO_ACID_COUNT],
}

impl StandardMatrix {
    /// Look up the score for aligning two amino acid letters.
    #[inline]
    pub fn score(&self, a: i8, b: i8) -> i8 {
        self.scores[a as usize * AMINO_ACID_COUNT + b as usize]
    }

    /// Get the statistical parameters for the given gap penalties.
    /// Returns the ungapped parameters if no matching gap penalties are found.
    pub fn constants(&self, gap_exist: i32, gap_extend: i32) -> &Parameters {
        for p in self.parameters.iter().skip(1) {
            if p.gap_exist as i32 == gap_exist && p.gap_extend as i32 == gap_extend {
                return p;
            }
        }
        &self.parameters[0]
    }

    /// Get the ungapped statistical parameters (first entry).
    pub fn ungapped_constants(&self) -> &Parameters {
        &self.parameters[0]
    }
}
