use std::collections::{BTreeSet, HashMap};

use crate::align::hsp::Hsp;
use crate::output::intermediate::IntermediateRecord;
use crate::util::interval::IntervalPartition;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CullingResult {
    Finished = 0,
    Next = 1,
    Include = 2,
}

#[derive(Debug, Clone)]
pub struct GlobalCulling {
    pub max_target_seqs: i64,
    pub n: i64,
    pub top_score: f64,
    pub taxon_count: HashMap<u32, u32>,
}

impl GlobalCulling {
    pub fn new(max_target_seqs: i64) -> Self {
        Self {
            max_target_seqs,
            n: 0,
            top_score: 0.0,
            taxon_count: HashMap::new(),
        }
    }

    /// Matches C++ `GlobalCulling::cull(const Target&)`.
    pub fn cull_target<F>(
        &self,
        filter_score: i32,
        taxon_rank_ids: &[u32],
        taxon_k: u32,
        toppercent: Option<f64>,
        mut bitscore: F,
    ) -> (CullingResult, f64)
    where
        F: FnMut(i32) -> f64,
    {
        if self.top_score == 0.0 {
            return (CullingResult::Include, 0.0);
        }
        if taxon_k != 0 {
            let mut taxons_exceeded = 0;
            for &taxon in taxon_rank_ids {
                if self.taxon_count.get(&taxon).copied().unwrap_or(0) >= taxon_k {
                    taxons_exceeded += 1;
                }
            }
            if taxons_exceeded == taxon_rank_ids.len() {
                return (CullingResult::Next, 0.0);
            }
        }
        if let Some(toppercent) = toppercent {
            (
                if (1.0 - bitscore(filter_score) / self.top_score) * 100.0 <= toppercent {
                    CullingResult::Include
                } else {
                    CullingResult::Finished
                },
                0.0,
            )
        } else {
            (
                if self.n < self.max_target_seqs {
                    CullingResult::Include
                } else {
                    CullingResult::Finished
                },
                0.0,
            )
        }
    }

    /// Matches C++ `GlobalCulling::cull(const vector<IntermediateRecord>&, const set<TaxId>&)`.
    pub fn cull_records<F>(
        &self,
        target_hsp: &[IntermediateRecord],
        taxon_ids: &BTreeSet<u32>,
        taxon_k: u32,
        global_ranking_targets: i64,
        toppercent: Option<f64>,
        mut bitscore: F,
    ) -> CullingResult
    where
        F: FnMut(u32) -> f64,
    {
        if self.top_score == 0.0 {
            return CullingResult::Include;
        }
        if taxon_k != 0 {
            let mut taxons_exceeded = 0;
            for &taxon in taxon_ids {
                if self.taxon_count.get(&taxon).copied().unwrap_or(0) >= taxon_k {
                    taxons_exceeded += 1;
                }
            }
            if taxons_exceeded == taxon_ids.len() {
                return CullingResult::Next;
            }
        }
        if global_ranking_targets != 0 {
            if self.n < global_ranking_targets {
                CullingResult::Include
            } else {
                CullingResult::Finished
            }
        } else if let Some(toppercent) = toppercent {
            if (1.0 - bitscore(target_hsp[0].score) / self.top_score) * 100.0 <= toppercent {
                CullingResult::Include
            } else {
                CullingResult::Finished
            }
        } else if self.n < self.max_target_seqs {
            CullingResult::Include
        } else {
            CullingResult::Finished
        }
    }

    /// Matches C++ `GlobalCulling::add(const Target&)`.
    pub fn add_target<F>(
        &mut self,
        filter_score: i32,
        taxon_rank_ids: &[u32],
        taxon_k: u32,
        mut bitscore: F,
    ) where
        F: FnMut(i32) -> f64,
    {
        if self.top_score == 0.0 {
            self.top_score = bitscore(filter_score);
        }
        self.n += 1;
        if taxon_k != 0 {
            for &taxon in taxon_rank_ids {
                *self.taxon_count.entry(taxon).or_insert(0) += 1;
            }
        }
    }

    /// Matches C++ `GlobalCulling::add(const vector<IntermediateRecord>&, const set<TaxId>&)`.
    pub fn add_records<F>(
        &mut self,
        target_hsp: &[IntermediateRecord],
        taxon_ids: &BTreeSet<u32>,
        taxon_k: u32,
        mut bitscore: F,
    ) where
        F: FnMut(u32) -> f64,
    {
        if self.top_score == 0.0 {
            self.top_score = bitscore(target_hsp[0].score);
        }
        self.n += 1;
        if taxon_k != 0 {
            for &taxon in taxon_ids {
                *self.taxon_count.entry(taxon).or_insert(0) += 1;
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct RangeCulling {
    pub p: IntervalPartition,
}

impl RangeCulling {
    pub fn new(max_target_seqs: i64) -> Self {
        Self {
            p: IntervalPartition::new(max_target_seqs),
        }
    }

    /// Matches C++ `RangeCulling::cull(const Target&)`.
    pub fn cull_hsps(
        &self,
        hsps: &[Hsp],
        toppercent: Option<f64>,
        query_range_cover: f64,
    ) -> (CullingResult, f64) {
        let mut c = 0;
        let mut l = 0;
        for hsp in hsps {
            if let Some(toppercent) = toppercent {
                let cutoff = (hsp.score as f64 / (1.0 - toppercent / 100.0)) as i32;
                c += self.p.covered_max_score(hsp.query_source_range, cutoff);
            } else {
                c += self.p.covered(hsp.query_source_range);
            }
            l += hsp.query_source_range.length();
        }
        let cov = c as f64 / l as f64;
        (
            if cov * 100.0 < query_range_cover {
                CullingResult::Include
            } else {
                CullingResult::Next
            },
            cov,
        )
    }

    /// Matches C++ `RangeCulling::cull(const vector<IntermediateRecord>&, const set<TaxId>&)`.
    pub fn cull_records(
        &self,
        target_hsp: &[IntermediateRecord],
        toppercent: Option<f64>,
        query_range_cover: f64,
    ) -> CullingResult {
        let mut c = 0;
        let mut l = 0;
        for hsp in target_hsp {
            let range = hsp.absolute_query_range();
            if let Some(toppercent) = toppercent {
                let cutoff = (hsp.score as f64 / (1.0 - toppercent / 100.0)) as i32;
                c += self.p.covered_max_score(range, cutoff);
            } else {
                c += self.p.covered(range);
            }
            l += range.length();
        }
        if c as f64 / l as f64 * 100.0 < query_range_cover {
            CullingResult::Include
        } else {
            CullingResult::Next
        }
    }

    /// Matches C++ `RangeCulling::add(const Target&)`.
    pub fn add_hsps(&mut self, hsps: &[Hsp]) {
        for hsp in hsps {
            self.p.insert(hsp.query_source_range, hsp.score);
        }
    }

    /// Matches C++ `RangeCulling::add(const vector<IntermediateRecord>&, const set<TaxId>&)`.
    pub fn add_records(&mut self, target_hsp: &[IntermediateRecord]) {
        for hsp in target_hsp {
            self.p.insert(hsp.absolute_query_range(), hsp.score as i32);
        }
    }
}

#[derive(Debug, Clone)]
pub enum TargetCulling {
    Global(GlobalCulling),
    Range(RangeCulling),
}

impl TargetCulling {
    /// Matches C++ `TargetCulling::get(const int64_t max_target_seqs)`.
    pub fn get(max_target_seqs: i64, query_range_culling: bool) -> Self {
        if query_range_culling {
            Self::Range(RangeCulling::new(max_target_seqs))
        } else {
            Self::Global(GlobalCulling::new(max_target_seqs))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::interval::Interval;

    #[test]
    fn test_global_culling_target_toppercent_and_limit() {
        let mut c = GlobalCulling::new(2);
        assert_eq!(
            c.cull_target(100, &[], 0, None, |score| score as f64).0,
            CullingResult::Include
        );
        c.add_target(100, &[], 0, |score| score as f64);
        c.add_target(95, &[], 0, |score| score as f64);
        assert_eq!(
            c.cull_target(90, &[], 0, None, |score| score as f64).0,
            CullingResult::Finished
        );
        assert_eq!(
            c.cull_target(91, &[], 0, Some(10.0), |score| score as f64)
                .0,
            CullingResult::Include
        );
        assert_eq!(
            c.cull_target(89, &[], 0, Some(10.0), |score| score as f64)
                .0,
            CullingResult::Finished
        );
    }

    #[test]
    fn test_global_culling_taxon_k() {
        let mut c = GlobalCulling::new(10);
        c.add_target(100, &[7, 8], 1, |score| score as f64);
        assert_eq!(
            c.cull_target(90, &[7], 1, None, |score| score as f64).0,
            CullingResult::Next
        );
        assert_eq!(
            c.cull_target(90, &[7, 9], 1, None, |score| score as f64).0,
            CullingResult::Include
        );
    }

    #[test]
    fn test_range_culling_hsps_and_records() {
        let mut range = RangeCulling::new(1);
        let mut h1 = Hsp::new();
        h1.score = 100;
        h1.query_source_range = Interval::new(0, 100);
        range.add_hsps(&[h1.clone()]);

        let mut h2 = Hsp::new();
        h2.score = 90;
        h2.query_source_range = Interval::new(50, 150);
        assert_eq!(
            range.cull_hsps(&[h2.clone()], None, 90.0).0,
            CullingResult::Include
        );
        assert_eq!(range.cull_hsps(&[h2], None, 50.0).0, CullingResult::Next);

        let record = IntermediateRecord {
            score: 90,
            query_begin: 25,
            query_end: 74,
            ..IntermediateRecord::default()
        };
        assert_eq!(
            range.cull_records(&[record], None, 50.0),
            CullingResult::Next
        );
    }

    #[test]
    fn test_target_culling_get() {
        assert!(matches!(
            TargetCulling::get(5, false),
            TargetCulling::Global(_)
        ));
        assert!(matches!(
            TargetCulling::get(5, true),
            TargetCulling::Range(_)
        ));
    }
}
