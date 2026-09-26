use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

type Fold = (char, i32);
type Family = (char, i32, i32, i32);

#[derive(Debug, Clone)]
pub struct RocConfig {
    pub cut_bar: bool,
    pub log_evalue_scale: f64,
    pub query_count: usize,
    pub output_hits: bool,
    pub output_fp: bool,
    pub no_forward_fp: bool,
    pub check_multi_target: bool,
    pub family_cap: usize,
    pub get_roc: bool,
}

impl Default for RocConfig {
    fn default() -> Self {
        Self {
            cut_bar: false,
            log_evalue_scale: 10.0,
            query_count: 1,
            output_hits: false,
            output_fp: false,
            no_forward_fp: false,
            check_multi_target: false,
            family_cap: 0,
            get_roc: false,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct FamilyMapping {
    values: HashMap<String, Vec<usize>>,
}

impl FamilyMapping {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn parse(
        text: &str,
        cut_bar: bool,
        family_ids: &mut BTreeMap<Family, usize>,
        family_folds: &mut BTreeMap<usize, Fold>,
    ) -> Result<Self, String> {
        let mut mapping = Self::new();
        for (line_no, line) in text.lines().enumerate() {
            if line.is_empty() {
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() < 7 {
                return Err(format!("Format error on mapping line {}.", line_no + 1));
            }
            let mut accession = fields[1].to_owned();
            let class = fields[3];
            if accession.is_empty() || class.chars().count() != 1 {
                return Err(format!("Format error on mapping line {}.", line_no + 1));
            }
            let parse = |field: &str| {
                field
                    .parse::<i32>()
                    .map_err(|_| format!("Format error on mapping line {}.", line_no + 1))
            };
            let family = (
                class.chars().next().unwrap(),
                parse(fields[4])?,
                parse(fields[5])?,
                parse(fields[6])?,
            );
            let next = family_ids.len();
            let id = *family_ids.entry(family).or_insert(next);
            if cut_bar {
                accession = accession
                    .rsplit('|')
                    .next()
                    .unwrap_or(&accession)
                    .to_owned();
            }
            mapping.values.entry(accession).or_default().push(id);
            family_folds.insert(id, (family.0, family.1));
        }
        Ok(mapping)
    }

    pub fn get(&self, accession: &str) -> &[usize] {
        self.values.get(accession).map(Vec::as_slice).unwrap_or(&[])
    }

    pub fn len(&self) -> usize {
        self.values.values().map(Vec::len).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
}

fn coverage(count: usize, family: usize, family_count: &[usize]) -> f64 {
    let total = family_count[family];
    if total == 0 {
        1.0
    } else {
        count as f64 / total as f64
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Histogram {
    pub bin_offset: i32,
    pub false_positives: Vec<usize>,
    pub coverage: Vec<f64>,
    scale: f64,
}

impl Histogram {
    pub const MAX_EV: f64 = 10_000.0;

    pub fn new(scale: f64) -> Result<Self, String> {
        if !scale.is_finite() || scale <= 0.0 {
            return Err("log_evalue_scale must be positive and finite".to_owned());
        }
        let bin_offset = (-((f64::MIN_EXP as f64) * std::f64::consts::LN_2 * scale).floor()) as i32;
        let bin_count = bin_offset + (Self::MAX_EV.ln() * scale).round() as i32 + 1;
        let bin_count =
            usize::try_from(bin_count).map_err(|_| "Invalid histogram size".to_owned())?;
        Ok(Self {
            bin_offset,
            false_positives: vec![0; bin_count],
            coverage: vec![0.0; bin_count],
            scale,
        })
    }

    pub fn bin(&self, evalue: f64) -> Result<usize, String> {
        if !evalue.is_finite() || evalue < 0.0 {
            return Err("Invalid E-value".to_owned());
        }
        if evalue == 0.0 {
            return Ok(0);
        }
        let bin = ((evalue.ln() * self.scale).round() as i32 + self.bin_offset).max(0);
        let bin = usize::try_from(bin).map_err(|_| "Evalue exceeds binning range.".to_owned())?;
        if bin >= self.coverage.len() {
            return Err("Evalue exceeds binning range.".to_owned());
        }
        Ok(bin)
    }

    pub fn add_assign(&mut self, other: &Self) -> Result<(), String> {
        if self.coverage.len() != other.coverage.len() {
            return Err("Histogram sizes differ".to_owned());
        }
        for i in 0..self.coverage.len() {
            self.false_positives[i] += other.false_positives[i];
            self.coverage[i] += other.coverage[i];
        }
        Ok(())
    }

    pub fn format(&self, query_count: usize) -> Result<String, String> {
        if query_count == 0 {
            return Err("query_count must be nonzero".to_owned());
        }
        let mut out = String::new();
        for i in 0..self.coverage.len() {
            writeln!(
                out,
                "{}\t{}",
                self.coverage[i] / query_count as f64,
                self.false_positives[i] as f64 / query_count as f64
            )
            .unwrap();
        }
        Ok(out)
    }
}

#[derive(Debug)]
pub struct QueryStats {
    query: String,
    last_subject: String,
    count: Vec<usize>,
    query_family: Vec<bool>,
    query_fold: BTreeSet<Fold>,
    family_idx: BTreeMap<usize, usize>,
    false_positives: Vec<usize>,
    previous_targets: HashSet<String>,
    true_positives: Vec<Vec<usize>>,
    pub have_rev_hit: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HitClass {
    Unknown,
    TruePositive,
    FalsePositive,
}

impl QueryStats {
    pub fn new(
        query: &str,
        families: usize,
        query_mapping: &FamilyMapping,
        family_folds: &BTreeMap<usize, Fold>,
        histogram_bins: usize,
        config: &RocConfig,
    ) -> Self {
        let mapped = query_mapping.get(query);
        let mut query_family = vec![false; families];
        let mut query_fold = BTreeSet::new();
        let mut family_idx = BTreeMap::new();
        for &family in mapped {
            if config.output_hits {
                query_family[family] = true;
            }
            if !config.no_forward_fp {
                query_fold.insert(family_folds[&family]);
            }
            if config.get_roc {
                let next = family_idx.len();
                family_idx.entry(family).or_insert(next);
            }
        }
        let n = family_idx.len();
        Self {
            query: query.to_owned(),
            last_subject: String::new(),
            count: vec![0; families],
            query_family,
            query_fold,
            family_idx,
            false_positives: vec![0; histogram_bins],
            previous_targets: HashSet::new(),
            true_positives: vec![vec![0; histogram_bins]; n],
            have_rev_hit: false,
        }
    }

    pub fn add_family_hit(&mut self, family: usize, bin: usize) {
        if let Some(&index) = self.family_idx.get(&family) {
            self.true_positives[index][bin] += 1;
        }
    }

    pub fn add(
        &mut self,
        subject: &str,
        evalue: f64,
        mapping: &FamilyMapping,
        family_folds: &BTreeMap<usize, Fold>,
        histogram: &Histogram,
        config: &RocConfig,
    ) -> Result<HitClass, String> {
        if (self.have_rev_hit && !config.get_roc) || subject == self.last_subject {
            return Ok(HitClass::Unknown);
        }
        if config.check_multi_target && !self.previous_targets.insert(subject.to_owned()) {
            return Ok(HitClass::Unknown);
        }
        self.last_subject = subject.to_owned();
        let bin = config.get_roc.then(|| histogram.bin(evalue)).transpose()?;
        if subject.starts_with('\\') {
            self.have_rev_hit = true;
            if let Some(bin) = bin {
                self.false_positives[bin] += 1;
            }
            return Ok(HitClass::FalsePositive);
        }
        let families = mapping.get(subject);
        if families.is_empty() {
            return Err("Accession not mapped.".to_owned());
        }
        let mut match_query = false;
        let mut same_fold = false;
        for &family in families {
            if !self.have_rev_hit {
                self.count[family] += 1;
            }
            if let Some(bin) = bin {
                self.add_family_hit(family, bin);
            }
            match_query |= config.output_hits && self.query_family[family];
            same_fold |= !config.no_forward_fp && self.query_fold.contains(&family_folds[&family]);
        }
        if !config.no_forward_fp && !same_fold {
            self.have_rev_hit = true;
            if let Some(bin) = bin {
                self.false_positives[bin] += 1;
            }
            Ok(HitClass::FalsePositive)
        } else if match_query {
            Ok(HitClass::TruePositive)
        } else {
            Ok(HitClass::Unknown)
        }
    }

    pub fn auc1(
        &self,
        family_count: &[usize],
        query_mapping: &FamilyMapping,
    ) -> Result<f64, String> {
        let families = query_mapping.get(&self.query);
        if families.is_empty() {
            return Err("Query accession not mapped.".to_owned());
        }
        Ok(families
            .iter()
            .map(|&family| coverage(self.count[family], family, family_count))
            .sum::<f64>()
            / families.len() as f64)
    }

    pub fn family_count(&self) -> usize {
        self.family_idx.len()
    }

    pub fn update_hist(&self, histogram: &mut Histogram, family_count: &[usize]) {
        let mut false_positives = 0;
        let mut true_positives = vec![0; self.family_count()];
        for bin in 0..histogram.coverage.len() {
            false_positives += self.false_positives[bin];
            histogram.false_positives[bin] += false_positives;
            let mut cov = 0.0;
            for (&family, &index) in &self.family_idx {
                true_positives[index] += self.true_positives[index][bin];
                cov += coverage(true_positives[index], family, family_count);
            }
            if !true_positives.is_empty() {
                histogram.coverage[bin] += cov / true_positives.len() as f64;
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct RocResult {
    pub output: String,
    pub histogram: Histogram,
    pub records: usize,
    pub queries: usize,
    pub queries_with_fp: usize,
    pub mappings: usize,
    pub query_mappings: usize,
    pub families: usize,
}

struct RocState<'a> {
    mapping: &'a FamilyMapping,
    query_mapping: &'a FamilyMapping,
    family_folds: &'a BTreeMap<usize, Fold>,
    family_count: &'a [usize],
    config: &'a RocConfig,
}

fn query_roc(
    buf: &str,
    histogram: &mut Histogram,
    state: &RocState<'_>,
) -> Result<(f64, String, bool), String> {
    let query = buf
        .lines()
        .next()
        .and_then(|line| line.split('\t').next())
        .unwrap_or("");
    let mut stats = QueryStats::new(
        query,
        state.family_count.len(),
        state.query_mapping,
        state.family_folds,
        histogram.coverage.len(),
        state.config,
    );
    let mut output = String::new();
    for line in buf.lines() {
        if line.is_empty() || (stats.have_rev_hit && !state.config.get_roc) {
            break;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 2 {
            return Err("Format error.".to_owned());
        }
        let evalue = if state.config.get_roc {
            fields
                .get(2)
                .ok_or_else(|| "Format error.".to_owned())?
                .parse()
                .map_err(|_| "Format error.".to_owned())?
        } else {
            0.0
        };
        let class = stats.add(
            fields[1],
            evalue,
            state.mapping,
            state.family_folds,
            histogram,
            state.config,
        )?;
        if (class == HitClass::TruePositive && state.config.output_hits)
            || (class == HitClass::FalsePositive && state.config.output_fp)
        {
            writeln!(output, "{line}").unwrap();
        }
    }
    let auc = stats.auc1(state.family_count, state.query_mapping)?;
    if state.config.get_roc {
        stats.update_hist(histogram, state.family_count);
    }
    if !state.config.output_hits && !state.config.output_fp {
        writeln!(output, "{}\t{}", stats.query, auc).unwrap();
    }
    Ok((auc, output, stats.have_rev_hit))
}

fn worker(
    groups: &[String],
    histogram: &mut Histogram,
    state: &RocState<'_>,
) -> Result<(String, usize), String> {
    let mut output = String::new();
    let mut queries_with_fp = 0;
    for group in groups {
        let (_, text, have_fp) = query_roc(group, histogram, state)?;
        output.push_str(&text);
        queries_with_fp += have_fp as usize;
    }
    Ok((output, queries_with_fp))
}

pub fn roc_from_text(
    mapping_text: &str,
    query_mapping_text: &str,
    alignments: &str,
    config: &RocConfig,
) -> Result<RocResult, String> {
    let mut family_ids = BTreeMap::new();
    let mut family_folds = BTreeMap::new();
    let mapping = FamilyMapping::parse(
        mapping_text,
        config.cut_bar,
        &mut family_ids,
        &mut family_folds,
    )?;
    let query_mapping = FamilyMapping::parse(
        query_mapping_text,
        config.cut_bar,
        &mut family_ids,
        &mut family_folds,
    )?;
    let mut family_count = vec![0; family_ids.len()];
    for families in mapping.values.values() {
        for &family in families {
            if config.family_cap == 0 {
                family_count[family] += 1;
            } else {
                family_count[family] = config.family_cap;
            }
        }
    }
    let mut histogram = Histogram::new(config.log_evalue_scale)?;
    let state = RocState {
        mapping: &mapping,
        query_mapping: &query_mapping,
        family_folds: &family_folds,
        family_count: &family_count,
        config,
    };
    let mut groups: Vec<String> = Vec::new();
    for line in alignments.lines().filter(|line| !line.is_empty()) {
        let query = line.split('\t').next().unwrap_or("");
        if groups
            .last()
            .and_then(|group| group.lines().next())
            .and_then(|first| first.split('\t').next())
            != Some(query)
        {
            groups.push(String::new());
        }
        writeln!(groups.last_mut().unwrap(), "{line}").unwrap();
    }
    // The C++ command feeds the same independent query groups to up to six
    // workers. This deterministic adapter owns one shard; callers can shard
    // the groups externally without reintroducing process-global queues.
    let (output, queries_with_fp) = worker(&groups, &mut histogram, &state)?;
    Ok(RocResult {
        output,
        histogram,
        records: alignments.lines().filter(|line| !line.is_empty()).count(),
        queries: groups.len(),
        queries_with_fp,
        mappings: mapping.len(),
        query_mappings: query_mapping.len(),
        families: family_ids.len(),
    })
}

pub fn roc(
    family_map: impl AsRef<Path>,
    family_map_query: impl AsRef<Path>,
    alignment_file: impl AsRef<Path>,
    config: &RocConfig,
) -> Result<RocResult, String> {
    let mapping = fs::read_to_string(family_map).map_err(|error| error.to_string())?;
    let queries = fs::read_to_string(family_map_query).map_err(|error| error.to_string())?;
    let alignments = fs::read_to_string(alignment_file).map_err(|error| error.to_string())?;
    roc_from_text(&mapping, &queries, &alignments, config)
}

#[cfg(test)]
mod tests {
    use super::*;

    const MAP: &str = "x\ta\tx\tA\t1\t1\t1\nx\tb\tx\tA\t1\t1\t1\nx\tc\tx\tB\t2\t1\t1\n";
    const QMAP: &str = "x\tq\tx\tA\t1\t1\t1\n";

    #[test]
    fn mapping_cut_bar_and_registry_are_stable() {
        let mut ids = BTreeMap::new();
        let mut folds = BTreeMap::new();
        let parsed =
            FamilyMapping::parse("x\tsp|a\tx\tA\t1\t2\t3\n", true, &mut ids, &mut folds).unwrap();
        assert_eq!(parsed.get("a"), [0]);
        assert_eq!(folds[&0], ('A', 1));
    }

    #[test]
    fn histogram_bins_zero_and_accumulates() {
        let mut a = Histogram::new(10.0).unwrap();
        let mut b = Histogram::new(10.0).unwrap();
        assert_eq!(a.bin(0.0).unwrap(), 0);
        let bin = a.bin(1.0).unwrap();
        b.false_positives[bin] = 2;
        a.add_assign(&b).unwrap();
        assert_eq!(a.false_positives[bin], 2);
    }

    #[test]
    fn roc_classifies_same_family_and_first_false_positive() {
        let config = RocConfig {
            output_hits: true,
            output_fp: true,
            get_roc: true,
            query_count: 1,
            ..Default::default()
        };
        let result =
            roc_from_text(MAP, QMAP, "q\ta\t1e-20\nq\tc\t1e-5\nq\tb\t1e-3\n", &config).unwrap();
        assert_eq!(result.records, 3);
        assert_eq!(result.queries, 1);
        assert_eq!(result.queries_with_fp, 1);
        assert!(result.output.contains("q\ta\t1e-20"));
        assert!(result.output.contains("q\tc\t1e-5"));
        assert!(result.histogram.false_positives.iter().any(|&n| n == 1));
    }

    #[test]
    fn roc_auc_output_and_duplicate_suppression() {
        let config = RocConfig::default();
        let result = roc_from_text(MAP, QMAP, "q\ta\nq\ta\nq\tb\n", &config).unwrap();
        assert_eq!(result.output, "q\t1\n");
        assert_eq!(result.queries_with_fp, 0);
    }
}
