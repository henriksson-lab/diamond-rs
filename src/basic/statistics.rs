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

impl Statistics {
    /// Matches C++ `Statistics::Statistics`.
    pub fn new() -> Self {
        let mut value = Self {
            data: [0; StatValue::COUNT],
        };
        value.reset();
        value
    }

    /// Matches C++ `Statistics::reset`.
    pub fn reset(&mut self) {
        self.data.fill(0);
    }

    /// Matches C++ `Statistics::operator+=`.
    pub fn add_assign_statistics(&mut self, rhs: &Statistics) -> &mut Self {
        for i in 0..StatValue::COUNT {
            self.data[i] += rhs.data[i];
        }
        self
    }

    /// Matches C++ `Statistics::inc`.
    pub fn inc(&mut self, v: StatValue, n: StatType) {
        self.data[v as usize] += n;
    }

    /// Matches C++ `Statistics::max`.
    pub fn max(&mut self, v: StatValue, n: StatType) {
        self.data[v as usize] = self.data[v as usize].max(n);
    }

    /// Matches C++ `Statistics::get`.
    pub fn get(&self, v: StatValue) -> StatType {
        self.data[v as usize]
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
}
