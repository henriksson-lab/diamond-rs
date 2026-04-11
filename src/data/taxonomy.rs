use crate::basic::value::TaxId;

/// A taxonomy node in the NCBI taxonomy tree.
#[derive(Debug, Clone)]
pub struct TaxonomyNode {
    pub taxid: TaxId,
    pub parent: TaxId,
    pub rank: String,
    pub name: String,
}

/// Taxonomy tree for taxonomic classification.
#[derive(Debug, Clone, Default)]
pub struct TaxonomyTree {
    nodes: Vec<TaxonomyNode>,
    parent_map: Vec<TaxId>, // index by taxid, value is parent taxid
}

impl TaxonomyTree {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a node to the tree.
    pub fn add_node(&mut self, node: TaxonomyNode) {
        let taxid = node.taxid as usize;
        if taxid >= self.parent_map.len() {
            self.parent_map.resize(taxid + 1, 0);
        }
        self.parent_map[taxid] = node.parent;
        self.nodes.push(node);
    }

    /// Get the parent taxid.
    pub fn parent(&self, taxid: TaxId) -> TaxId {
        let idx = taxid as usize;
        if idx < self.parent_map.len() {
            self.parent_map[idx]
        } else {
            0
        }
    }

    /// Get the lineage (path to root) for a taxid.
    pub fn lineage(&self, taxid: TaxId) -> Vec<TaxId> {
        let mut path = Vec::new();
        let mut current = taxid;
        while current > 1 {
            path.push(current);
            let parent = self.parent(current);
            if parent == current {
                break;
            }
            current = parent;
        }
        if current == 1 {
            path.push(1); // root
        }
        path
    }

    /// Find the lowest common ancestor of two taxids.
    pub fn lca(&self, a: TaxId, b: TaxId) -> TaxId {
        let lineage_a: std::collections::HashSet<TaxId> =
            self.lineage(a).into_iter().collect();
        for taxid in self.lineage(b) {
            if lineage_a.contains(&taxid) {
                return taxid;
            }
        }
        1 // root
    }

    /// Number of nodes.
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }
}

/// Mapping from sequence accessions to taxids.
#[derive(Debug, Clone, Default)]
pub struct TaxonMapping {
    mapping: std::collections::HashMap<String, TaxId>,
}

impl TaxonMapping {
    pub fn new() -> Self {
        Self::default()
    }

    /// Load mapping from a TSV file (accession\ttaxid).
    pub fn load_tsv<R: std::io::BufRead>(reader: R) -> std::io::Result<Self> {
        let mut mapping = std::collections::HashMap::new();
        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                if let Ok(taxid) = parts[1].parse::<TaxId>() {
                    mapping.insert(parts[0].to_string(), taxid);
                }
            }
        }
        Ok(TaxonMapping { mapping })
    }

    /// Look up the taxid for an accession.
    pub fn get(&self, accession: &str) -> Option<TaxId> {
        self.mapping.get(accession).copied()
    }

    pub fn len(&self) -> usize {
        self.mapping.len()
    }

    pub fn is_empty(&self) -> bool {
        self.mapping.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn test_taxonomy_tree() {
        let mut tree = TaxonomyTree::new();
        tree.add_node(TaxonomyNode {
            taxid: 1,
            parent: 1,
            rank: "root".to_string(),
            name: "root".to_string(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 2,
            parent: 1,
            rank: "superkingdom".to_string(),
            name: "Bacteria".to_string(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 10,
            parent: 2,
            rank: "phylum".to_string(),
            name: "Proteobacteria".to_string(),
        });

        assert_eq!(tree.parent(10), 2);
        assert_eq!(tree.parent(2), 1);

        let lineage = tree.lineage(10);
        assert_eq!(lineage, vec![10, 2, 1]);
    }

    #[test]
    fn test_lca() {
        let mut tree = TaxonomyTree::new();
        tree.add_node(TaxonomyNode {
            taxid: 1, parent: 1, rank: "root".into(), name: "root".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 2, parent: 1, rank: "".into(), name: "A".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 3, parent: 2, rank: "".into(), name: "B".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 4, parent: 2, rank: "".into(), name: "C".into(),
        });

        assert_eq!(tree.lca(3, 4), 2);
        assert_eq!(tree.lca(3, 2), 2);
    }

    #[test]
    fn test_taxon_mapping() {
        let tsv = "acc1\t9606\nacc2\t10090\n";
        let mapping =
            TaxonMapping::load_tsv(BufReader::new(tsv.as_bytes())).unwrap();
        assert_eq!(mapping.get("acc1"), Some(9606));
        assert_eq!(mapping.get("acc2"), Some(10090));
        assert_eq!(mapping.get("acc3"), None);
    }

    #[test]
    fn test_load_real_mapping() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/diamond/src/test/acc2taxid.tsv"
        );
        if let Ok(file) = std::fs::File::open(path) {
            let mapping =
                TaxonMapping::load_tsv(BufReader::new(file)).unwrap();
            assert!(!mapping.is_empty());
        }
    }
}
