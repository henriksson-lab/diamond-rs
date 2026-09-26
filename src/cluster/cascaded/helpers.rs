//! Cascaded-clustering round policy and edge callbacks.
//!
//! This mirrors `diamond/src/cluster/cascaded/helpers.cpp` and the inline
//! callback surface in `cascaded.h`. Upstream reads global configuration;
//! Rust exposes that state explicitly through `CascadedHelpersConfig`.

use crate::basic::value::OId;
use crate::output::edge::EdgeData;
use crate::util::algo::Edge;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CascadedHelpersConfig {
    pub cluster_steps: Vec<String>,
    pub connected_component_depth: Vec<String>,
}

/// Upstream `Cascaded::get_key()`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Cascaded;

impl Cascaded {
    pub const fn get_key() -> &'static str {
        "cascaded"
    }
}

/// Select default clustering rounds when no explicit override is configured.
pub fn cluster_steps(approx_id: f64, linear: bool) -> Vec<String> {
    let mut steps = vec!["faster_lin".to_string()];
    if approx_id < 90.0 {
        steps.push("fast_lin".to_string());
    }
    if approx_id < 40.0 {
        steps.push("linclust-20_lin".to_string());
    } else if approx_id < 80.0 {
        steps.push("linclust-40_lin".to_string());
    }
    if linear {
        return steps;
    }
    if approx_id < 80.0 {
        steps.push("default".to_string());
    } else {
        steps.push("fast".to_string());
    }
    if approx_id < 50.0 {
        steps.push("more-sensitive".to_string());
    }
    steps
}

pub fn cluster_steps_with_config(
    approx_id: f64,
    linear: bool,
    config: &CascadedHelpersConfig,
) -> Vec<String> {
    if !config.cluster_steps.is_empty() {
        return config.cluster_steps.clone();
    }
    cluster_steps(approx_id, linear)
}

pub fn is_linclust(steps: &[String]) -> bool {
    steps.iter().all(|step| step.ends_with("_lin"))
}

pub fn default_round_approx_id(_steps: i32) -> Vec<String> {
    Vec::new()
}

pub fn default_round_cov(_steps: i32) -> Vec<String> {
    Vec::new()
}

struct RoundCcdParser;

impl RoundCcdParser {
    /// Private C++ `round_ccd(const string&)` overload.
    fn round_ccd(depth: &str) -> Result<i32, String> {
        // `stringstream >> int` accepts leading whitespace, but the following
        // `!eof()` rejects any trailing character, including whitespace.
        let value = depth.trim_start();
        if value.is_empty() || value.trim_end().len() != value.len() {
            return Err("Invalid number format for --connected-component-depth".to_string());
        }
        value
            .parse::<i32>()
            .map_err(|_| "Invalid number format for --connected-component-depth".to_string())
    }
}

pub fn round_ccd(round: i32, round_count: i32, linear: bool) -> Result<i32, String> {
    round_ccd_with_config(
        round,
        round_count,
        linear,
        &CascadedHelpersConfig::default(),
    )
}

pub fn round_ccd_with_config(
    round: i32,
    round_count: i32,
    linear: bool,
    config: &CascadedHelpersConfig,
) -> Result<i32, String> {
    let depths = &config.connected_component_depth;
    if depths.len() > 1 && depths.len() != round_count as usize {
        return Err("Parameter count for --connected-component-depth has to be 1 or the number of cascaded clustering rounds.".to_string());
    }
    if depths.is_empty() {
        return Ok(0);
    }
    if depths.len() > 1 {
        return RoundCcdParser::round_ccd(&depths[round as usize]);
    }
    if (round == round_count - 1) ^ linear {
        Ok(1)
    } else {
        RoundCcdParser::round_ccd(&depths[0])
    }
}

pub trait EdgeCallback {
    fn consume(&mut self, bytes: &[u8]) -> Result<(), String>;
    fn count(&self) -> i64;
    fn edges(&self) -> &[Edge<OId>];
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct CallbackUnidirectional {
    pub member_cover: f64,
    pub edge_file: Vec<Edge<OId>>,
    pub count: i64,
}

impl CallbackUnidirectional {
    pub fn new(member_cover: f64) -> Self {
        Self {
            member_cover,
            edge_file: Vec::new(),
            count: 0,
        }
    }
}

impl EdgeCallback for CallbackUnidirectional {
    fn consume(&mut self, bytes: &[u8]) -> Result<(), String> {
        for edge in decode_edges(bytes)? {
            if f64::from(edge.qcovhsp) >= self.member_cover {
                self.edge_file
                    .push(Edge::new(edge.target, edge.query, edge.evalue));
                self.count += 1;
            }
            if f64::from(edge.scovhsp) >= self.member_cover {
                self.edge_file
                    .push(Edge::new(edge.query, edge.target, edge.evalue));
                self.count += 1;
            }
        }
        Ok(())
    }

    fn count(&self) -> i64 {
        self.count
    }

    fn edges(&self) -> &[Edge<OId>] {
        &self.edge_file
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct CallbackBidirectional {
    pub edge_file: Vec<Edge<OId>>,
    pub count: i64,
}

impl EdgeCallback for CallbackBidirectional {
    fn consume(&mut self, bytes: &[u8]) -> Result<(), String> {
        for edge in decode_edges(bytes)? {
            if edge.query != edge.target {
                self.edge_file
                    .push(Edge::new(edge.target, edge.query, edge.evalue));
                self.edge_file
                    .push(Edge::new(edge.query, edge.target, edge.evalue));
                self.count += 2;
            }
        }
        Ok(())
    }

    fn count(&self) -> i64 {
        self.count
    }

    fn edges(&self) -> &[Edge<OId>] {
        &self.edge_file
    }
}

fn decode_edges(bytes: &[u8]) -> Result<Vec<EdgeData>, String> {
    if bytes.len() % EdgeData::SIZE != 0 {
        return Err("Invalid edge buffer size".to_string());
    }
    let mut edges = Vec::with_capacity(bytes.len() / EdgeData::SIZE);
    for record in bytes.chunks_exact(EdgeData::SIZE) {
        edges.push(EdgeData {
            query: u64::from_ne_bytes(record[0..8].try_into().unwrap()),
            target: u64::from_ne_bytes(record[8..16].try_into().unwrap()),
            qcovhsp: f32::from_ne_bytes(record[16..20].try_into().unwrap()),
            scovhsp: f32::from_ne_bytes(record[20..24].try_into().unwrap()),
            evalue: f64::from_ne_bytes(record[24..32].try_into().unwrap()),
        });
    }
    Ok(edges)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_string()).collect()
    }

    fn edge_bytes(edges: &[EdgeData]) -> Vec<u8> {
        let mut bytes = Vec::new();
        for edge in edges {
            edge.write(&mut bytes).unwrap();
        }
        bytes
    }

    #[test]
    fn cluster_step_thresholds_and_override_are_exact() {
        assert_eq!(cluster_steps(90.0, false), strings(&["faster_lin", "fast"]));
        assert_eq!(
            cluster_steps(80.0, false),
            strings(&["faster_lin", "fast_lin", "fast"])
        );
        assert_eq!(
            cluster_steps(79.999, false),
            strings(&["faster_lin", "fast_lin", "linclust-40_lin", "default"])
        );
        assert_eq!(
            cluster_steps(39.999, false),
            strings(&[
                "faster_lin",
                "fast_lin",
                "linclust-20_lin",
                "default",
                "more-sensitive"
            ])
        );
        assert_eq!(
            cluster_steps(39.999, true),
            strings(&["faster_lin", "fast_lin", "linclust-20_lin"])
        );
        let config = CascadedHelpersConfig {
            cluster_steps: strings(&["custom", "steps"]),
            ..Default::default()
        };
        assert_eq!(
            cluster_steps_with_config(0.0, false, &config),
            config.cluster_steps
        );
        assert!(is_linclust(&[]));
        assert!(is_linclust(&strings(&["a_lin", "b_lin"])));
        assert!(!is_linclust(&strings(&["a_lin", "default"])));
    }

    #[test]
    fn connected_component_depth_count_selection_and_parse_match() {
        let empty = CascadedHelpersConfig::default();
        assert_eq!(round_ccd_with_config(0, 3, false, &empty).unwrap(), 0);

        let one = CascadedHelpersConfig {
            connected_component_depth: strings(&[" 7"]),
            ..Default::default()
        };
        assert_eq!(round_ccd_with_config(0, 3, false, &one).unwrap(), 7);
        assert_eq!(round_ccd_with_config(2, 3, false, &one).unwrap(), 1);
        assert_eq!(round_ccd_with_config(0, 3, true, &one).unwrap(), 1);
        assert_eq!(round_ccd_with_config(2, 3, true, &one).unwrap(), 7);

        let per_round = CascadedHelpersConfig {
            connected_component_depth: strings(&["2", "3", "4"]),
            ..Default::default()
        };
        assert_eq!(round_ccd_with_config(1, 3, false, &per_round).unwrap(), 3);
        assert!(round_ccd_with_config(0, 2, false, &per_round).is_err());
        let invalid = CascadedHelpersConfig {
            connected_component_depth: strings(&["7 "]),
            ..Default::default()
        };
        assert_eq!(
            round_ccd_with_config(0, 2, false, &invalid).unwrap_err(),
            "Invalid number format for --connected-component-depth"
        );
        assert!(default_round_approx_id(3).is_empty());
        assert!(default_round_cov(3).is_empty());
    }

    #[test]
    fn unidirectional_callback_applies_cover_independently_and_in_order() {
        let bytes = edge_bytes(&[
            EdgeData {
                query: 1,
                target: 2,
                qcovhsp: 80.0,
                scovhsp: 79.0,
                evalue: 1e-5,
            },
            EdgeData {
                query: 3,
                target: 4,
                qcovhsp: 90.0,
                scovhsp: 95.0,
                evalue: 2e-5,
            },
        ]);
        let mut callback = CallbackUnidirectional::new(80.0);
        callback.consume(&bytes).unwrap();
        assert_eq!(callback.count(), 3);
        assert_eq!(
            callback.edges(),
            &[
                Edge::new(2, 1, 1e-5),
                Edge::new(4, 3, 2e-5),
                Edge::new(3, 4, 2e-5)
            ]
        );
    }

    #[test]
    fn bidirectional_callback_suppresses_self_edges_and_rejects_partial_records() {
        let bytes = edge_bytes(&[
            EdgeData {
                query: 5,
                target: 5,
                qcovhsp: 100.0,
                scovhsp: 100.0,
                evalue: 0.0,
            },
            EdgeData {
                query: 6,
                target: 7,
                qcovhsp: 1.0,
                scovhsp: 2.0,
                evalue: 3.0,
            },
        ]);
        let mut callback = CallbackBidirectional::default();
        callback.consume(&bytes).unwrap();
        assert_eq!(callback.count(), 2);
        assert_eq!(
            callback.edges(),
            &[Edge::new(7, 6, 3.0), Edge::new(6, 7, 3.0)]
        );
        assert_eq!(
            callback.consume(&[0; 31]).unwrap_err(),
            "Invalid edge buffer size"
        );
    }

    #[test]
    fn cascaded_key_is_stable() {
        assert_eq!(Cascaded::get_key(), "cascaded");
    }
}
