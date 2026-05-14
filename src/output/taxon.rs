use std::io::{self, Write};

use crate::align::hsp::HspContext;
use crate::basic::value::{OId, TaxId};
use crate::data::taxonomy::TaxonomyTree;
use crate::output::format::format_evalue;
use crate::util::sequence::ID_DELIMITERS;

#[derive(Debug, Clone)]
pub struct TaxonFormat {
    pub taxid: TaxId,
    pub evalue: f64,
    pub needs_taxon_id_lists: bool,
    pub needs_taxon_nodes: bool,
    pub needs_taxon_scientific_names: bool,
    pub include_lineage: bool,
}

impl TaxonFormat {
    pub fn new(include_lineage: bool) -> Self {
        Self {
            taxid: 0,
            evalue: f64::MAX,
            needs_taxon_id_lists: true,
            needs_taxon_nodes: true,
            needs_taxon_scientific_names: include_lineage,
            include_lineage,
        }
    }

    pub fn print_match(&mut self, r: &HspContext, subject_taxids: &[TaxId], tree: &TaxonomyTree) {
        if subject_taxids.is_empty() {
            return;
        }
        self.evalue = self.evalue.min(r.evalue());
        for &taxid in subject_taxids {
            self.taxid = tree.lca(self.taxid, taxid);
        }
    }

    pub fn print_query_epilog<W: Write>(
        &self,
        writer: &mut W,
        query_title: &str,
        tree: &TaxonomyTree,
    ) -> io::Result<()> {
        write!(writer, "{}\t{}\t", query_id(query_title), self.taxid)?;
        if self.taxid > 0 {
            write!(writer, "{}", format_evalue(self.evalue))?;
        } else {
            write!(writer, "0")?;
        }
        if self.include_lineage {
            write!(
                writer,
                "\t{}",
                if self.taxid > 0 {
                    taxon_lineage(self.taxid, tree)
                } else {
                    "N/A".to_string()
                }
            )?;
        }
        writeln!(writer)
    }

    pub fn reset(&mut self) {
        self.taxid = 0;
        self.evalue = f64::MAX;
    }
}

pub fn taxon_lineage(taxid: TaxId, tree: &TaxonomyTree) -> String {
    let lin = tree.lineage_root_to_taxid(taxid);
    if lin.is_empty() {
        return "N/A".to_string();
    }
    let mut out = tree.taxon_scientific_name(lin[0]);
    for &taxid in &lin[1..] {
        out.push_str("; ");
        out.push_str(&tree.taxon_scientific_name(taxid));
    }
    out
}

pub fn query_id(title: &str) -> &str {
    title
        .split(|c: char| ID_DELIMITERS.contains(c))
        .next()
        .unwrap_or("")
}

pub fn print_match(
    format: &mut TaxonFormat,
    r: &HspContext,
    subject_taxids: &[TaxId],
    tree: &TaxonomyTree,
) {
    format.print_match(r, subject_taxids, tree);
}

pub fn print_query_epilog<W: Write>(
    format: &TaxonFormat,
    writer: &mut W,
    query_title: &str,
    tree: &TaxonomyTree,
) -> io::Result<()> {
    format.print_query_epilog(writer, query_title, tree)
}

pub fn subject_taxids<'a>(subject_oid: OId, taxids_by_oid: &'a [Vec<TaxId>]) -> &'a [TaxId] {
    taxids_by_oid
        .get(subject_oid as usize)
        .map(Vec::as_slice)
        .unwrap_or(&[])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::Hsp;
    use crate::data::taxonomy::TaxonomyNode;

    fn test_tree() -> TaxonomyTree {
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
        tree.add_node(TaxonomyNode {
            taxid: 20,
            parent: 2,
            rank: "phylum".into(),
            name: "Firmicutes".into(),
        });
        tree
    }

    #[test]
    fn test_taxon_lineage() {
        let tree = test_tree();
        assert_eq!(taxon_lineage(10, &tree), "Bacteria; Proteobacteria");
        assert_eq!(taxon_lineage(0, &tree), "N/A");
    }

    #[test]
    fn test_query_id_delimiters() {
        assert_eq!(query_id("q1 description"), "q1");
        assert_eq!(query_id("q1\x01other"), "q1");
        assert_eq!(query_id("q1\tother"), "q1");
    }

    #[test]
    fn test_taxon_format_match_and_epilog_with_lineage() {
        let tree = test_tree();
        let mut hsp = Hsp::new();
        hsp.evalue = 1.0e-20;
        let r = HspContext::new(
            hsp,
            0,
            0,
            Vec::new(),
            0,
            "query one",
            0,
            0,
            "",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        );
        let mut format = TaxonFormat::new(true);
        format.print_match(&r, &[10, 20], &tree);
        assert_eq!(format.taxid, 2);
        assert_eq!(format.evalue, 1.0e-20);
        let mut out = Vec::new();
        format
            .print_query_epilog(&mut out, &r.query_title, &tree)
            .unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "query\t2\t1.00e-20\tBacteria\n"
        );
    }

    #[test]
    fn test_taxon_format_empty_and_subject_lookup() {
        let tree = test_tree();
        let r = HspContext::default();
        let mut format = TaxonFormat::new(false);
        format.print_match(&r, &[], &tree);
        let mut out = Vec::new();
        format
            .print_query_epilog(&mut out, "query two", &tree)
            .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "query\t0\t0\n");
        let taxids = vec![vec![10], vec![20, 10]];
        assert_eq!(subject_taxids(1, &taxids), &[20, 10]);
        assert!(subject_taxids(7, &taxids).is_empty());
    }
}
