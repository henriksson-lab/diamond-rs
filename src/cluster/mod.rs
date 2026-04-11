use crate::basic::value::OId;

/// A protein sequence cluster — a centroid with its member sequences.
#[derive(Debug, Clone)]
pub struct Cluster {
    /// The representative (centroid) sequence ID.
    pub centroid: OId,
    /// Member sequence IDs (including the centroid).
    pub members: Vec<OId>,
}

/// Result of a clustering run.
#[derive(Debug, Clone)]
pub struct ClusterResult {
    /// All clusters.
    pub clusters: Vec<Cluster>,
    /// Mapping from sequence OId to its centroid OId.
    pub mapping: Vec<OId>,
}

impl ClusterResult {
    pub fn new() -> Self {
        ClusterResult {
            clusters: Vec::new(),
            mapping: Vec::new(),
        }
    }

    /// Number of clusters.
    pub fn num_clusters(&self) -> usize {
        self.clusters.len()
    }

    /// Build the mapping vector from the clusters.
    pub fn build_mapping(&mut self, num_sequences: usize) {
        self.mapping = vec![0; num_sequences];
        for cluster in &self.clusters {
            for &member in &cluster.members {
                self.mapping[member as usize] = cluster.centroid;
            }
        }
    }

    /// Get the centroid for a given sequence.
    pub fn get_centroid(&self, seq_id: OId) -> OId {
        self.mapping[seq_id as usize]
    }
}

impl Default for ClusterResult {
    fn default() -> Self {
        Self::new()
    }
}

/// Clustering algorithm types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClusterAlgo {
    /// Cascaded multi-round clustering (default).
    Cascaded,
    /// Linear-time clustering.
    Linclust,
    /// Deep clustering.
    DeepClust,
}

impl ClusterAlgo {
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "cascaded" => Some(ClusterAlgo::Cascaded),
            "linclust" => Some(ClusterAlgo::Linclust),
            "deepclust" => Some(ClusterAlgo::DeepClust),
            _ => None,
        }
    }
}

/// Write cluster results in TSV format (centroid_id \t member_id).
pub fn write_clusters<W: std::io::Write>(
    writer: &mut W,
    result: &ClusterResult,
    ids: &[String],
) -> std::io::Result<()> {
    for cluster in &result.clusters {
        let centroid_name = &ids[cluster.centroid as usize];
        for &member in &cluster.members {
            let member_name = &ids[member as usize];
            writeln!(writer, "{}\t{}", centroid_name, member_name)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cluster_result() {
        let mut result = ClusterResult::new();
        result.clusters.push(Cluster {
            centroid: 0,
            members: vec![0, 1, 2],
        });
        result.clusters.push(Cluster {
            centroid: 3,
            members: vec![3, 4],
        });
        result.build_mapping(5);

        assert_eq!(result.num_clusters(), 2);
        assert_eq!(result.get_centroid(0), 0);
        assert_eq!(result.get_centroid(1), 0);
        assert_eq!(result.get_centroid(4), 3);
    }

    #[test]
    fn test_write_clusters() {
        let result = ClusterResult {
            clusters: vec![Cluster {
                centroid: 0,
                members: vec![0, 1],
            }],
            mapping: vec![],
        };
        let ids = vec!["seq_A".to_string(), "seq_B".to_string()];
        let mut buf = Vec::new();
        write_clusters(&mut buf, &result, &ids).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert_eq!(output, "seq_A\tseq_A\nseq_A\tseq_B\n");
    }

    #[test]
    fn test_cluster_algo_parse() {
        assert_eq!(ClusterAlgo::parse("cascaded"), Some(ClusterAlgo::Cascaded));
        assert_eq!(ClusterAlgo::parse("linclust"), Some(ClusterAlgo::Linclust));
        assert_eq!(ClusterAlgo::parse("unknown"), None);
    }
}
