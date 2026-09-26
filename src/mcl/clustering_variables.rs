use crate::align::hsp::HspContext;
use crate::dp::swipe::HspValues;
use crate::output::format::OutputFlags;
use std::collections::BTreeMap;
use std::fmt;
use std::sync::LazyLock;

type Getter = fn(&HspContext, bool) -> f64;

/// One numeric HSP property available to MCL expressions.
///
/// The C++ implementation uses one subclass per property. A function pointer
/// retains the same immutable polymorphic registry without heap allocation.
#[derive(Debug, Clone, Copy)]
pub struct Variable {
    name: &'static str,
    pub hsp_values: HspValues,
    pub flags: OutputFlags,
    getter: Getter,
}

impl Variable {
    const fn new(
        name: &'static str,
        hsp_values: HspValues,
        flags: OutputFlags,
        getter: Getter,
    ) -> Self {
        Self {
            name,
            hsp_values,
            flags,
            getter,
        }
    }

    pub const fn get_name(&self) -> &'static str {
        self.name
    }

    /// Evaluate this property.
    ///
    /// `query_translated` replaces the alignment-mode process global used by
    /// C++ `HspContext::blast_query_frame()`.
    pub fn get(&self, context: &HspContext, query_translated: bool) -> f64 {
        (self.getter)(context, query_translated)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnknownVariableError {
    key: String,
}

impl fmt::Display for UnknownVariableError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Unknown variable: {}", self.key)
    }
}

impl std::error::Error for UnknownVariableError {}

/// Immutable counterpart of C++ `StaticVariableRegistry`.
pub struct StaticVariableRegistry {
    reg_map: BTreeMap<&'static str, Variable>,
}

impl StaticVariableRegistry {
    pub fn new() -> Self {
        let mut reg_map = BTreeMap::new();
        for variable in VARIABLES {
            reg_map.insert(variable.get_name(), variable);
        }
        Self { reg_map }
    }

    pub fn get(&self, key: &str) -> Result<&Variable, UnknownVariableError> {
        self.reg_map.get(key).ok_or_else(|| UnknownVariableError {
            key: key.to_owned(),
        })
    }

    pub fn has(&self, key: &str) -> bool {
        self.reg_map.contains_key(key)
    }

    /// Return keys in the lexicographic order provided by C++ `std::map`.
    pub fn get_keys(&self) -> Vec<&'static str> {
        self.reg_map.keys().copied().collect()
    }
}

impl Default for StaticVariableRegistry {
    fn default() -> Self {
        Self::new()
    }
}

/// Facade matching the source's process-wide `VariableRegistry`.
pub struct VariableRegistry;

impl VariableRegistry {
    pub fn get(key: &str) -> Result<&'static Variable, UnknownVariableError> {
        VARIABLE_REGISTRY.get(key)
    }

    pub fn has(key: &str) -> bool {
        VARIABLE_REGISTRY.has(key)
    }

    pub fn get_keys() -> Vec<&'static str> {
        VARIABLE_REGISTRY.get_keys()
    }
}

static VARIABLE_REGISTRY: LazyLock<StaticVariableRegistry> =
    LazyLock::new(StaticVariableRegistry::new);

fn query_length(r: &HspContext, _: bool) -> f64 {
    r.query_source_len as f64
}

fn subject_length(r: &HspContext, _: bool) -> f64 {
    r.subject_len as f64
}

fn query_start(r: &HspContext, _: bool) -> f64 {
    (r.oriented_query_range().begin + 1) as f64
}

fn query_end(r: &HspContext, _: bool) -> f64 {
    (r.oriented_query_range().end + 1) as f64
}

fn subject_start(r: &HspContext, _: bool) -> f64 {
    (r.subject_range().begin + 1) as f64
}

fn subject_end(r: &HspContext, _: bool) -> f64 {
    r.subject_range().end as f64
}

fn evalue(r: &HspContext, _: bool) -> f64 {
    r.evalue()
}

fn bit_score(r: &HspContext, _: bool) -> f64 {
    r.bit_score()
}

fn raw_score(r: &HspContext, _: bool) -> f64 {
    r.score() as f64
}

fn length(r: &HspContext, _: bool) -> f64 {
    r.length() as f64
}

fn percent_identical_matches(r: &HspContext, _: bool) -> f64 {
    r.identities() as f64 * 100.0 / r.length() as f64
}

fn number_identical_matches(r: &HspContext, _: bool) -> f64 {
    r.identities() as f64
}

fn number_mismatches(r: &HspContext, _: bool) -> f64 {
    r.mismatches() as f64
}

fn number_positive_matches(r: &HspContext, _: bool) -> f64 {
    r.positives() as f64
}

fn number_gap_openings(r: &HspContext, _: bool) -> f64 {
    r.gap_openings() as f64
}

fn number_gaps(r: &HspContext, _: bool) -> f64 {
    r.gaps() as f64
}

fn percentage_positive_matches(r: &HspContext, _: bool) -> f64 {
    r.positives() as f64 * 100.0 / r.length() as f64
}

fn query_frame(r: &HspContext, query_translated: bool) -> f64 {
    r.blast_query_frame(query_translated) as f64
}

fn query_coverage_per_hsp(r: &HspContext, _: bool) -> f64 {
    r.query_source_range().length() as f64 * 100.0 / r.query_source_len as f64
}

fn subject_coverage_per_hsp(r: &HspContext, _: bool) -> f64 {
    r.subject_range().length() as f64 * 100.0 / r.subject_len as f64
}

fn normalized_bit_score_global(r: &HspContext, _: bool) -> f64 {
    r.bit_score() / r.query_self_aln_score.max(r.target_self_aln_score) * 100.0
}

const VARIABLES: [Variable; 21] = [
    Variable::new("qlen", HspValues::NONE, OutputFlags::NONE, query_length),
    Variable::new("slen", HspValues::NONE, OutputFlags::NONE, subject_length),
    Variable::new(
        "qstart",
        HspValues::QUERY_START,
        OutputFlags::NONE,
        query_start,
    ),
    Variable::new("qend", HspValues::QUERY_END, OutputFlags::NONE, query_end),
    Variable::new(
        "sstart",
        HspValues::TARGET_START,
        OutputFlags::NONE,
        subject_start,
    ),
    Variable::new(
        "send",
        HspValues::TARGET_END,
        OutputFlags::NONE,
        subject_end,
    ),
    Variable::new("evalue", HspValues::NONE, OutputFlags::NONE, evalue),
    Variable::new("bitscore", HspValues::NONE, OutputFlags::NONE, bit_score),
    Variable::new("score", HspValues::NONE, OutputFlags::NONE, raw_score),
    Variable::new("length", HspValues::LENGTH, OutputFlags::NONE, length),
    Variable::new(
        "pident",
        HspValues(HspValues::LENGTH.0 | HspValues::IDENT.0),
        OutputFlags::NONE,
        percent_identical_matches,
    ),
    Variable::new(
        "nident",
        HspValues::IDENT,
        OutputFlags::NONE,
        number_identical_matches,
    ),
    Variable::new(
        "mismatch",
        HspValues::MISMATCHES,
        OutputFlags::NONE,
        number_mismatches,
    ),
    Variable::new(
        "positive",
        HspValues::TRANSCRIPT,
        OutputFlags::NONE,
        number_positive_matches,
    ),
    Variable::new(
        "gapopen",
        HspValues::GAP_OPENINGS,
        OutputFlags::NONE,
        number_gap_openings,
    ),
    Variable::new("gaps", HspValues::GAPS, OutputFlags::NONE, number_gaps),
    Variable::new(
        "ppos",
        HspValues::TRANSCRIPT,
        OutputFlags::NONE,
        percentage_positive_matches,
    ),
    Variable::new("qframe", HspValues::NONE, OutputFlags::NONE, query_frame),
    Variable::new(
        "qcovhsp",
        HspValues::QUERY_COORDS,
        OutputFlags::NONE,
        query_coverage_per_hsp,
    ),
    Variable::new(
        "scovhsp",
        HspValues::TARGET_COORDS,
        OutputFlags::NONE,
        subject_coverage_per_hsp,
    ),
    Variable::new(
        "normalized_bitscore_global",
        HspValues::NONE,
        OutputFlags::SELF_ALN_SCORES,
        normalized_bit_score_global,
    ),
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;
    use crate::util::interval::Interval;

    fn context() -> HspContext {
        let hsp = Hsp {
            score: 60,
            evalue: 2.5e-8,
            bit_score: 50.0,
            frame: 1,
            length: 12,
            identities: 9,
            mismatches: 3,
            positives: 10,
            gap_openings: 2,
            gaps: 4,
            query_source_range: Interval::new(4, 13),
            subject_range: Interval::new(7, 17),
            ..Hsp::default()
        };
        HspContext::new(
            hsp,
            0,
            0,
            vec![vec![0; 10]; 6],
            30,
            "query",
            0,
            40,
            "subject",
            1,
            1,
            vec![0; 40],
            80.0,
            125.0,
        )
    }

    fn value(name: &str, context: &HspContext) -> f64 {
        VariableRegistry::get(name).unwrap().get(context, true)
    }

    #[test]
    fn registry_contains_every_source_variable_in_map_order() {
        assert_eq!(
            VariableRegistry::get_keys(),
            vec![
                "bitscore",
                "evalue",
                "gapopen",
                "gaps",
                "length",
                "mismatch",
                "nident",
                "normalized_bitscore_global",
                "pident",
                "positive",
                "ppos",
                "qcovhsp",
                "qend",
                "qframe",
                "qlen",
                "qstart",
                "score",
                "scovhsp",
                "send",
                "slen",
                "sstart",
            ]
        );
        assert!(VariableRegistry::has("pident"));
        assert!(!VariableRegistry::has("unknown"));
        assert_eq!(
            VariableRegistry::get("unknown").unwrap_err().to_string(),
            "Unknown variable: unknown"
        );
    }

    #[test]
    fn variables_evaluate_all_source_formulas() {
        let context = context();
        let expected = [
            ("qlen", 30.0),
            ("slen", 40.0),
            ("qstart", 5.0),
            ("qend", 13.0),
            ("sstart", 8.0),
            ("send", 17.0),
            ("evalue", 2.5e-8),
            ("bitscore", 50.0),
            ("score", 60.0),
            ("length", 12.0),
            ("pident", 75.0),
            ("nident", 9.0),
            ("mismatch", 3.0),
            ("positive", 10.0),
            ("gapopen", 2.0),
            ("gaps", 4.0),
            ("ppos", 1000.0 / 12.0),
            ("qframe", 2.0),
            ("qcovhsp", 30.0),
            ("scovhsp", 25.0),
            ("normalized_bitscore_global", 40.0),
        ];
        for (name, expected) in expected {
            assert!((value(name, &context) - expected).abs() < 1e-12, "{name}");
        }
        assert_eq!(
            VariableRegistry::get("qframe")
                .unwrap()
                .get(&context, false),
            0.0
        );
    }

    #[test]
    fn metadata_matches_hsp_and_output_requirements() {
        let pident = VariableRegistry::get("pident").unwrap();
        assert_eq!(pident.hsp_values, HspValues::LENGTH | HspValues::IDENT);
        assert_eq!(
            VariableRegistry::get("gaps").unwrap().hsp_values,
            HspValues::GAPS
        );
        assert_eq!(
            VariableRegistry::get("positive").unwrap().hsp_values,
            HspValues::TRANSCRIPT
        );
        let normalized = VariableRegistry::get("normalized_bitscore_global").unwrap();
        assert_eq!(normalized.hsp_values, HspValues::NONE);
        assert_eq!(normalized.flags, OutputFlags::SELF_ALN_SCORES);
    }

    #[test]
    fn floating_point_edge_cases_follow_cpp_arithmetic() {
        let context = HspContext::default();
        assert!(value("pident", &context).is_nan());
        assert!(value("ppos", &context).is_nan());
        assert!(value("qcovhsp", &context).is_nan());
        assert!(value("scovhsp", &context).is_nan());
        assert!(value("normalized_bitscore_global", &context).is_nan());
    }
}
