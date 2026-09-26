use crate::basic::value::TaxId;
use crate::data::blastdb::taxdmp::read_names_dmp;

// Preserve the established `data::taxonomy::{Rank, TaxonomyNodes}` API while
// mirroring the upstream `data/taxonomy_nodes.{h,cpp}` source boundary.
pub use crate::data::taxonomy_nodes::{Rank, TaxonomyNodes};

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Taxonomy {
    pub name: Vec<String>,
}

impl Taxonomy {
    pub fn init(&mut self, namesdmp: &str) -> Result<usize, String> {
        if namesdmp.is_empty() {
            return Ok(0);
        }
        self.load_names(namesdmp)
    }

    pub fn load_names(&mut self, file_name: &str) -> Result<usize, String> {
        let mut n = 0usize;
        read_names_dmp(file_name, |id, name| {
            if id >= 0 {
                let i = id as usize;
                self.name.resize(i + 1, String::new());
                self.name[i] = name.to_string();
                n += 1;
            }
        })
        .map_err(|e| e.to_string())?;
        Ok(n)
    }
}

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
        if a == b || b <= 0 {
            return a;
        }
        if a <= 0 {
            return b;
        }
        let lineage_a: std::collections::HashSet<TaxId> = self.lineage(a).into_iter().collect();
        for taxid in self.lineage(b) {
            if lineage_a.contains(&taxid) {
                return taxid;
            }
        }
        1 // root
    }

    /// Matches C++ `SequenceFile::lineage(taxid)`.
    pub fn lineage_root_to_taxid(&self, taxid: TaxId) -> Vec<TaxId> {
        let mut out = Vec::new();
        let mut current = taxid;
        while current > 1 {
            out.push(current);
            let parent = self.parent(current);
            if parent == current || parent <= 0 {
                return Vec::new();
            }
            current = parent;
        }
        if current <= 0 {
            return Vec::new();
        }
        out.reverse();
        out
    }

    /// Matches C++ `SequenceFile::taxon_scientific_name(taxid)`.
    pub fn taxon_scientific_name(&self, taxid: TaxId) -> String {
        self.nodes
            .iter()
            .find(|node| node.taxid == taxid && !node.name.is_empty())
            .map(|node| node.name.clone())
            .unwrap_or_else(|| taxid.to_string())
    }

    pub fn rank(&self, taxid: TaxId) -> Option<&str> {
        self.nodes
            .iter()
            .find(|node| node.taxid == taxid)
            .map(|node| node.rank.as_str())
    }

    /// Matches C++ `SequenceFile::rank_taxid(taxid, rank)`.
    pub fn rank_taxid(&self, mut taxid: TaxId, rank: &str) -> TaxId {
        const MAX_LINEAGE: usize = 64;
        let mut n = 0usize;
        loop {
            if self.rank(taxid) == Some(rank) {
                return taxid;
            }
            if taxid <= 1 {
                return 0;
            }
            n += 1;
            if n > MAX_LINEAGE {
                panic!("Path in taxonomy too long (rank_taxid).");
            }
            taxid = self.parent(taxid);
        }
    }

    /// Matches C++ `SequenceFile::rank_taxid(taxids, rank)`.
    pub fn rank_taxids(&self, taxids: &[TaxId], rank: &str) -> std::collections::BTreeSet<TaxId> {
        let mut out = std::collections::BTreeSet::new();
        for &taxid in taxids {
            let ranked = self.rank_taxid(taxid, rank);
            if ranked > 0 {
                out.insert(ranked);
            }
        }
        out
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
    use crate::util::io::{Deserializer, Serializer, VecStream};
    use std::collections::BTreeSet;
    use std::io::BufReader;
    use std::io::Write;

    #[test]
    fn test_rank_cpp_order_names_and_parse() {
        assert_eq!(Rank::COUNT, 49);
        assert_eq!(Rank::None as usize, 0);
        assert_eq!(Rank::Superkingdom as usize, 1);
        assert_eq!(Rank::Phylum as usize, 9);
        assert_eq!(Rank::Species as usize, 34);
        assert_eq!(Rank::Subvariety as usize, 48);
        assert_eq!(Rank::NAMES[0], "no rank");
        assert_eq!(Rank::NAMES[48], "subvariety");
        assert_eq!(Rank::predefined("forma specialis"), Some(41));
        assert_eq!(Rank::predefined("missing"), None);
        assert_eq!(Rank::parse("genus").unwrap(), Rank::Genus);
        assert_eq!(Rank::parse("genus").unwrap().to_string(), "genus");
        assert!(Rank::parse("bad rank").is_err());
        assert_eq!(Rank::from_index(12), Some(Rank::ClassRank));
        assert_eq!(Rank::from_index(49), None);
    }

    #[test]
    fn test_taxonomy_load_names_and_init() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-taxonomy-names-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "1\t|\troot\t|\t\t|\tscientific name\t|").unwrap();
            writeln!(f, "2\t|\tBacteria\t|\t\t|\tscientific name\t|").unwrap();
            writeln!(f, "2\t|\tbacteria\t|\t\t|\tcommon name\t|").unwrap();
            writeln!(f, "1224\t|\tPseudomonadota\t|\t\t|\tscientific name\t|").unwrap();
        }

        let mut taxonomy = Taxonomy::default();
        assert_eq!(taxonomy.init(""), Ok(0));
        assert!(taxonomy.name.is_empty());
        assert_eq!(taxonomy.load_names(path.to_str().unwrap()).unwrap(), 3);
        assert_eq!(taxonomy.name[1], "root");
        assert_eq!(taxonomy.name[2], "Bacteria");
        assert_eq!(taxonomy.name[1224], "Pseudomonadota");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_taxonomy_nodes_from_nodes_dmp() {
        let path = std::env::temp_dir().join(format!(
            "diamond-rs-taxonomy-nodes-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "1\t|\t1\t|\tno rank\t|\t").unwrap();
            writeln!(f, "2\t|\t1\t|\tsuperkingdom\t|\t").unwrap();
            writeln!(f, "561\t|\t1224\t|\tspecies\t|\t").unwrap();
            writeln!(f, "1224\t|\t2\t|\tphylum\t|\t").unwrap();
        }

        let nodes = TaxonomyNodes::from_nodes_dmp(path.to_str().unwrap()).unwrap();
        assert_eq!(nodes.len(), 1225);
        assert_eq!(nodes.max(), 1224);
        assert_eq!(nodes.get_parent(561).unwrap(), 1224);
        assert_eq!(nodes.get_parent(1224).unwrap(), 2);
        assert_eq!(nodes.rank(561), Rank::Species as i32);
        assert_eq!(nodes.rank(1224), Rank::Phylum as i32);
        assert_eq!(nodes.rank(-1), -1);
        assert!(nodes.get_parent(5000).is_err());
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_taxonomy_nodes_contained() {
        let nodes = TaxonomyNodes {
            parent: vec![0, 1, 1, 2, 2, 1],
            rank: vec![Rank::None; 6],
        };
        let mut filter = BTreeSet::new();
        filter.insert(2);
        assert!(nodes.contained(4, &filter).unwrap());
        assert!(nodes.contained(3, &filter).unwrap());
        assert!(!nodes.contained(5, &filter).unwrap());
        assert!(nodes.contained_any(&[5, 4], &filter).unwrap());
        assert!(!nodes.contained_any(&[5], &filter).unwrap());
        filter.insert(1);
        assert!(nodes.contained(5, &filter).unwrap());
        assert!(nodes.contained_any(&[], &filter).unwrap());
    }

    #[test]
    fn test_taxonomy_nodes_save_and_load() {
        let nodes = TaxonomyNodes {
            parent: vec![0, 1, 1, 2, 1224],
            rank: vec![
                Rank::None,
                Rank::None,
                Rank::Superkingdom,
                Rank::Phylum,
                Rank::Species,
            ],
        };
        let mut out = Serializer::new(VecStream::new());
        nodes.save(&mut out).unwrap();
        let stream = out.into_inner().unwrap();
        let mut input = Deserializer::new(VecStream::from_vec(stream.data().to_vec()));
        let roundtrip = TaxonomyNodes::from_deserializer(&mut input, 131).unwrap();
        assert_eq!(roundtrip, nodes);

        let mut legacy_input = Deserializer::new(VecStream::from_vec(stream.data().to_vec()));
        let legacy = TaxonomyNodes::from_deserializer(&mut legacy_input, 130).unwrap();
        assert_eq!(legacy.parents(), nodes.parents());
        assert!(legacy.ranks().is_empty());
    }

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
            taxid: 1,
            parent: 1,
            rank: "root".into(),
            name: "root".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 2,
            parent: 1,
            rank: "".into(),
            name: "A".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 3,
            parent: 2,
            rank: "".into(),
            name: "B".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 4,
            parent: 2,
            rank: "".into(),
            name: "C".into(),
        });

        assert_eq!(tree.lca(3, 4), 2);
        assert_eq!(tree.lca(3, 2), 2);
        assert_eq!(tree.lca(0, 4), 4);
        assert_eq!(tree.lca(3, 0), 3);
    }

    #[test]
    fn test_taxonomy_tree_cpp_lineage_and_names() {
        let mut tree = TaxonomyTree::new();
        tree.add_node(TaxonomyNode {
            taxid: 1,
            parent: 1,
            rank: "root".into(),
            name: "root".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 2,
            parent: 1,
            rank: "superkingdom".into(),
            name: "Bacteria".into(),
        });
        tree.add_node(TaxonomyNode {
            taxid: 10,
            parent: 2,
            rank: "phylum".into(),
            name: "Proteobacteria".into(),
        });
        assert_eq!(tree.lineage_root_to_taxid(10), vec![2, 10]);
        assert_eq!(tree.taxon_scientific_name(10), "Proteobacteria");
        assert_eq!(tree.taxon_scientific_name(999), "999");
        assert_eq!(tree.rank(10), Some("phylum"));
        assert_eq!(tree.rank_taxid(10, "superkingdom"), 2);
        assert_eq!(tree.rank_taxid(2, "phylum"), 0);
        assert_eq!(
            tree.rank_taxids(&[10, 999], "superkingdom")
                .into_iter()
                .collect::<Vec<_>>(),
            vec![2]
        );
    }

    #[test]
    fn test_taxon_mapping() {
        let tsv = "acc1\t9606\nacc2\t10090\n";
        let mapping = TaxonMapping::load_tsv(BufReader::new(tsv.as_bytes())).unwrap();
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
            let mapping = TaxonMapping::load_tsv(BufReader::new(file)).unwrap();
            assert!(!mapping.is_empty());
        }
    }
}
