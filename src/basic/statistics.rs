//! Search statistics translated from `diamond/src/basic/statistics.h` and the
//! `Statistics::print` implementation in `diamond/src/basic/basic.cpp`.

use std::fmt::Write as _;
use std::io::{self, Write};
use std::sync::{LazyLock, Mutex};

pub type StatType = i64;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(usize)]
pub enum StatValue {
    SeedHits,
    TentativeMatches0,
    TentativeMatches1,
    TentativeMatches2,
    TentativeMatches3,
    TentativeMatches4,
    TentativeMatchesX,
    Matches,
    Aligned,
    Gapped,
    Duplicates,
    GappedHits,
    QuerySeeds,
    QuerySeedsHit,
    RefSeeds,
    RefSeedsHit,
    QuerySize,
    RefSize,
    OutHits,
    OutMatches,
    CollisionLookups,
    Qcov,
    BiasErrors,
    ScoreTotal,
    AlignedQlen,
    Pairwise,
    HighSim,
    SearchTempSpace,
    SecondaryHits,
    ErasedHits,
    SquaredError,
    Cells,
    TargetHits0,
    TargetHits1,
    TargetHits2,
    TargetHits3,
    TargetHits3Cbs,
    TargetHits4,
    TargetHits5,
    TargetHits6,
    TimeGreedyExt,
    LowComplexitySeeds,
    SwipeRealign,
    Ext8,
    Ext16,
    Ext32,
    GappedFilterTargets,
    GappedFilterHits1,
    GappedFilterHits2,
    GrossDpCells,
    NetDpCells,
    TimeTargetSort,
    TimeSw,
    TimeExt,
    TimeGappedFilter,
    TimeLoadHitTargets,
    TimeChaining,
    TimeLoadSeedHits,
    TimeSortSeedHits,
    TimeSortTargetsByScore,
    TimeTargetParallel,
    TimeTracebackSw,
    TimeTraceback,
    HardQueries,
    TimeMatrixAdjust,
    MatrixAdjustCount,
    CompBasedStatsCount,
    FailedCompBasedStats,
    MaskedLazy,
    SwipeTasksTotal,
    SwipeTasksAsync,
    TrivialAln,
    TimeExt32,
    ExtOverflow8,
    ExtWasted16,
    DpCells8,
    DpCells16,
    DpCells32,
    TimeProfile,
    TimeAnchoredSwipe,
    TimeAnchoredSwipeAlloc,
    TimeAnchoredSwipeSort,
    TimeAnchoredSwipeAdd,
    TimeAnchoredSwipeOutput,
    TimeProfileGeneration,
    ExtensionsRecompute,
    TimeSearch,
    SeedsHit,
    Count,
}

impl StatValue {
    pub const COUNT: usize = StatValue::Count as usize;
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Statistics {
    pub data: [StatType; StatValue::COUNT],
}

/// The three destinations used by C++ `Statistics::print` (`log_stream`,
/// `verbose_stream`, and `message_stream`), represented as owned text so the
/// caller controls where and whether each stream is emitted.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StatisticsReport {
    pub diagnostics: String,
    pub verbose: String,
    pub message: String,
}

impl Statistics {
    /// Matches C++ `Statistics::Statistics()`.
    pub fn new() -> Self {
        let mut value = Self {
            data: [0; StatValue::COUNT],
        };
        value.reset();
        value
    }

    /// Matches C++ `Statistics::reset()`.
    pub fn reset(&mut self) {
        self.data.fill(0);
    }

    /// Matches C++ `Statistics::operator+=(rhs)`.
    pub fn add_assign_statistics(&mut self, rhs: &Statistics) -> &mut Self {
        for i in 0..StatValue::COUNT {
            self.data[i] += rhs.data[i];
        }
        self
    }

    /// Matches C++ `Statistics::inc(v, n)`.
    pub fn inc(&mut self, v: StatValue, n: StatType) {
        self.data[v as usize] += n;
    }

    /// Matches C++ `Statistics::max(v, n)`.
    pub fn max(&mut self, v: StatValue, n: StatType) {
        self.data[v as usize] = self.data[v as usize].max(n);
    }

    /// Matches C++ `Statistics::get(v)`.
    pub fn get(&self, v: StatValue) -> StatType {
        self.data[v as usize]
    }

    /// Format the report produced by C++ `Statistics::print` without relying
    /// on its process-wide logging streams.
    pub fn report(&self) -> StatisticsReport {
        let d = |v| self.get(v);
        let pct = |numerator: StatType, denominator: StatType| {
            numerator as f64 * 100.0 / denominator as f64
        };
        let seconds = |v| d(v) as f64 / 1e6;
        let mut diagnostics = String::new();
        writeln!(
            diagnostics,
            "Seeds hit             = {}",
            d(StatValue::SeedsHit)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Hits (filter stage 0) = {}",
            d(StatValue::SeedHits)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Hits (filter stage 1) = {} ({} %)",
            d(StatValue::TentativeMatches1),
            pct(d(StatValue::TentativeMatches1), d(StatValue::SeedHits))
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Hits (filter stage 2) = {} ({} %)",
            d(StatValue::TentativeMatches2),
            pct(
                d(StatValue::TentativeMatches2),
                d(StatValue::TentativeMatches1)
            )
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Hits (filter stage 3) = {} ({} %)",
            d(StatValue::TentativeMatches3),
            pct(
                d(StatValue::TentativeMatches3),
                d(StatValue::TentativeMatches2)
            )
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Target hits (stage 0) = {}",
            d(StatValue::TargetHits0)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Target hits (stage 1) = {}",
            d(StatValue::TargetHits1)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Target hits (stage 2) = {}",
            d(StatValue::TargetHits2)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Target hits (stage 3) = {} ({} ({}%) with CBS)",
            d(StatValue::TargetHits3),
            d(StatValue::TargetHits3Cbs),
            pct(d(StatValue::TargetHits3Cbs), d(StatValue::TargetHits3))
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Target hits (stage 4) = {}",
            d(StatValue::TargetHits4)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Target hits (stage 5) = {}",
            d(StatValue::TargetHits5)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Target hits (stage 6) = {}",
            d(StatValue::TargetHits6)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Swipe realignments    = {}",
            d(StatValue::SwipeRealign)
        )
        .unwrap();
        if d(StatValue::MaskedLazy) != 0 {
            writeln!(
                diagnostics,
                "Lazy maskings         = {}",
                d(StatValue::MaskedLazy)
            )
            .unwrap();
        }
        writeln!(
            diagnostics,
            "Matrix adjusts        = {}",
            d(StatValue::MatrixAdjustCount)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Comp. based stats     = {}",
            d(StatValue::CompBasedStatsCount)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Failed cbs            = {}",
            d(StatValue::FailedCompBasedStats)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Extensions (8 bit)    = {}",
            d(StatValue::Ext8)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Extensions (16 bit)   = {}",
            d(StatValue::Ext16)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Extensions (32 bit)   = {}",
            d(StatValue::Ext32)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Extensions (Recompute)= {}",
            d(StatValue::ExtensionsRecompute)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Overflows (8 bit)     = {}",
            d(StatValue::ExtOverflow8)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Wasted (16 bit)       = {}",
            d(StatValue::ExtWasted16)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Effort (Extension)    = {}",
            2 * d(StatValue::Ext16) + d(StatValue::Ext8)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Effort (Cells)        = {}",
            2 * d(StatValue::DpCells16) + d(StatValue::DpCells8)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Cells (8 bit)         = {}",
            d(StatValue::DpCells8)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Cells (16 bit)        = {}",
            d(StatValue::DpCells16)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "SWIPE tasks           = {}",
            d(StatValue::SwipeTasksTotal)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "SWIPE tasks (async)   = {}",
            d(StatValue::SwipeTasksAsync)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Trivial aln           = {}",
            d(StatValue::TrivialAln)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Hard queries          = {}",
            d(StatValue::HardQueries)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Gapped filter (targets) = {}",
            d(StatValue::GappedFilterTargets)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Gapped filter (hits) stage 1 = {}",
            d(StatValue::GappedFilterHits1)
        )
        .unwrap();
        writeln!(
            diagnostics,
            "Gapped filter (hits) stage 2 = {}",
            d(StatValue::GappedFilterHits2)
        )
        .unwrap();

        let timings = [
            (
                "Time (search)                ",
                StatValue::TimeSearch,
                "wall",
            ),
            (
                "Time (Load seed hit targets) ",
                StatValue::TimeLoadHitTargets,
                "CPU",
            ),
            (
                "Time (Sort targets by score) ",
                StatValue::TimeSortTargetsByScore,
                "CPU",
            ),
            (
                "Time (Gapped filter)         ",
                StatValue::TimeGappedFilter,
                "CPU",
            ),
            (
                "Time (Matrix adjust)         ",
                StatValue::TimeMatrixAdjust,
                "CPU",
            ),
            (
                "Time (Profile generation)    ",
                StatValue::TimeProfileGeneration,
                "CPU",
            ),
            (
                "Time (Chaining)              ",
                StatValue::TimeChaining,
                "CPU",
            ),
            (
                "Time (DP target sorting)     ",
                StatValue::TimeTargetSort,
                "CPU",
            ),
            (
                "Time (Query profiles)        ",
                StatValue::TimeProfile,
                "CPU",
            ),
            ("Time (Smith Waterman)        ", StatValue::TimeSw, "CPU"),
            (
                "Time (Anchored SWIPE Alloc)  ",
                StatValue::TimeAnchoredSwipeAlloc,
                "CPU",
            ),
            (
                "Time (Anchored SWIPE Sort)   ",
                StatValue::TimeAnchoredSwipeSort,
                "CPU",
            ),
            (
                "Time (Anchored SWIPE Add)    ",
                StatValue::TimeAnchoredSwipeAdd,
                "CPU",
            ),
            (
                "Time (Anchored SWIPE Output) ",
                StatValue::TimeAnchoredSwipeOutput,
                "CPU",
            ),
            (
                "Time (Anchored SWIPE)        ",
                StatValue::TimeAnchoredSwipe,
                "CPU",
            ),
            (
                "Time (Smith Waterman TB)     ",
                StatValue::TimeTracebackSw,
                "CPU",
            ),
            ("Time (Smith Waterman-32)     ", StatValue::TimeExt32, "CPU"),
            (
                "Time (Traceback)             ",
                StatValue::TimeTraceback,
                "CPU",
            ),
            (
                "Time (Target parallel)       ",
                StatValue::TimeTargetParallel,
                "wall",
            ),
            (
                "Time (Load seed hits)        ",
                StatValue::TimeLoadSeedHits,
                "wall",
            ),
            (
                "Time (Sort seed hits)        ",
                StatValue::TimeSortSeedHits,
                "wall",
            ),
            ("Time (Extension)             ", StatValue::TimeExt, "wall"),
        ];
        for (label, value, clock) in timings {
            writeln!(diagnostics, "{}= {}s ({})", label, seconds(value), clock).unwrap();
        }

        let verbose = format!(
            "Temporary disk space used (search): {} GB\n",
            d(StatValue::SearchTempSpace) as f64 / (1u64 << 30) as f64
        );
        let message = format!(
            "Reported {} pairwise alignments, {} HSPs.\n{} queries aligned.\n",
            d(StatValue::Pairwise),
            d(StatValue::Matches),
            d(StatValue::Aligned)
        );
        StatisticsReport {
            diagnostics,
            verbose,
            message,
        }
    }

    /// Write the three report channels to caller-provided destinations.
    pub fn print<D: Write, V: Write, M: Write>(
        &self,
        diagnostics: &mut D,
        verbose: &mut V,
        message: &mut M,
    ) -> io::Result<()> {
        let report = self.report();
        diagnostics.write_all(report.diagnostics.as_bytes())?;
        verbose.write_all(report.verbose.as_bytes())?;
        message.write_all(report.message.as_bytes())
    }
}

impl Default for Statistics {
    fn default() -> Self {
        Self::new()
    }
}

impl std::ops::AddAssign<&Statistics> for Statistics {
    fn add_assign(&mut self, rhs: &Statistics) {
        self.add_assign_statistics(rhs);
    }
}

pub static STATISTICS: LazyLock<Mutex<Statistics>> =
    LazyLock::new(|| Mutex::new(Statistics::new()));

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_statistics_init_reset_inc_max_get() {
        let mut stats = Statistics::new();
        assert_eq!(stats.get(StatValue::SeedHits), 0);
        stats.inc(StatValue::SeedHits, 1);
        stats.inc(StatValue::SeedHits, 4);
        assert_eq!(stats.get(StatValue::SeedHits), 5);
        stats.max(StatValue::SeedHits, 3);
        assert_eq!(stats.get(StatValue::SeedHits), 5);
        stats.max(StatValue::SeedHits, 8);
        assert_eq!(stats.get(StatValue::SeedHits), 8);
        stats.reset();
        assert_eq!(stats.get(StatValue::SeedHits), 0);
    }

    #[test]
    fn test_statistics_add_assign() {
        let mut a = Statistics::new();
        let mut b = Statistics::new();
        a.inc(StatValue::Matches, 2);
        b.inc(StatValue::Matches, 3);
        b.inc(StatValue::Aligned, 7);
        a += &b;
        assert_eq!(a.get(StatValue::Matches), 5);
        assert_eq!(a.get(StatValue::Aligned), 7);
    }

    #[test]
    fn test_statistics_global() {
        let mut stats = STATISTICS.lock().unwrap();
        stats.reset();
        stats.inc(StatValue::TimeSearch, 11);
        assert_eq!(stats.get(StatValue::TimeSearch), 11);
        stats.reset();
    }

    #[test]
    fn test_statistics_indices_match_cpp_order() {
        assert_eq!(StatValue::SeedHits as usize, 0);
        assert_eq!(StatValue::Matches as usize, 7);
        assert_eq!(StatValue::Ext8 as usize, 43);
        assert_eq!(StatValue::Ext16 as usize, 44);
        assert_eq!(StatValue::Ext32 as usize, 45);
        assert_eq!(StatValue::SeedsHit as usize, StatValue::COUNT - 1);
        assert_eq!(StatValue::COUNT, 88);
    }

    #[test]
    fn test_statistics_report_preserves_stream_routing_and_derived_values() {
        let mut stats = Statistics::new();
        stats.inc(StatValue::SeedsHit, 12);
        stats.inc(StatValue::SeedHits, 10);
        stats.inc(StatValue::TentativeMatches1, 5);
        stats.inc(StatValue::TentativeMatches2, 4);
        stats.inc(StatValue::TentativeMatches3, 2);
        stats.inc(StatValue::TargetHits3, 8);
        stats.inc(StatValue::TargetHits3Cbs, 2);
        stats.inc(StatValue::Ext8, 3);
        stats.inc(StatValue::Ext16, 4);
        stats.inc(StatValue::MaskedLazy, 7);
        stats.inc(StatValue::TimeSearch, 1_500_000);
        stats.inc(StatValue::SearchTempSpace, 1 << 30);
        stats.inc(StatValue::Pairwise, 6);
        stats.inc(StatValue::Matches, 9);
        stats.inc(StatValue::Aligned, 3);

        let report = stats.report();
        assert!(report.diagnostics.contains("Seeds hit             = 12\n"));
        assert!(report
            .diagnostics
            .contains("Hits (filter stage 1) = 5 (50 %)\n"));
        assert!(report
            .diagnostics
            .contains("Target hits (stage 3) = 8 (2 (25%) with CBS)\n"));
        assert!(report.diagnostics.contains("Lazy maskings         = 7\n"));
        assert!(report.diagnostics.contains("Effort (Extension)    = 11\n"));
        assert!(report
            .diagnostics
            .contains("Time (search)                = 1.5s (wall)\n"));
        assert_eq!(report.verbose, "Temporary disk space used (search): 1 GB\n");
        assert_eq!(
            report.message,
            "Reported 6 pairwise alignments, 9 HSPs.\n3 queries aligned.\n"
        );
    }

    #[test]
    fn test_statistics_print_writes_separate_destinations() {
        let stats = Statistics::new();
        let (mut diagnostics, mut verbose, mut message) = (Vec::new(), Vec::new(), Vec::new());
        stats
            .print(&mut diagnostics, &mut verbose, &mut message)
            .unwrap();
        assert!(String::from_utf8(diagnostics)
            .unwrap()
            .starts_with("Seeds hit"));
        assert!(String::from_utf8(verbose)
            .unwrap()
            .starts_with("Temporary disk space"));
        assert!(String::from_utf8(message)
            .unwrap()
            .starts_with("Reported 0"));
    }
}
