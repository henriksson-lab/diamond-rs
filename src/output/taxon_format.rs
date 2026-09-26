//! Mirrored facade for `diamond/src/output/taxon_format.cpp` and the
//! `TaxonFormat` declaration in `output_format.h`.

pub use super::taxon::{
    print_match, print_query_epilog, query_id, subject_taxids, taxon_lineage, TaxonFormat,
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::hsp::{Hsp, HspContext};
    use crate::data::taxonomy::{TaxonomyNode, TaxonomyTree};
    use crate::dp::swipe::HspValues;
    use crate::output::format::{OutputFlags, OutputFormatKind};

    fn tree() -> TaxonomyTree {
        let mut tree = TaxonomyTree::new();
        for (taxid, parent, name) in [
            (1, 1, "root"),
            (2, 1, "Bacteria"),
            (10, 2, "Proteobacteria"),
            (11, 10, "Species A"),
            (20, 2, "Firmicutes"),
        ] {
            tree.add_node(TaxonomyNode {
                taxid,
                parent,
                rank: String::new(),
                name: name.into(),
            });
        }
        tree
    }

    fn context(evalue: f64) -> HspContext {
        let mut hsp = Hsp::new();
        hsp.evalue = evalue;
        HspContext::new(
            hsp,
            0,
            0,
            Vec::new(),
            0,
            "query one\tcomment",
            0,
            0,
            "",
            0,
            0,
            Vec::new(),
            0.0,
            0.0,
        )
    }

    #[test]
    fn constructor_state_matches_cpp_header_and_config() {
        let without_lineage = TaxonFormat::new(false).output_format_spec();
        assert_eq!(without_lineage.code, OutputFormatKind::Taxon);
        assert_eq!(without_lineage.hsp_values, HspValues::NONE);
        assert_eq!(without_lineage.flags, OutputFlags::DEFAULT_REPORT_UNALIGNED);
        assert!(without_lineage.needs_taxon_id_lists);
        assert!(without_lineage.needs_taxon_nodes);
        assert!(!without_lineage.needs_taxon_scientific_names);

        let with_lineage = TaxonFormat::new(true).output_format_spec();
        assert!(with_lineage.needs_taxon_scientific_names);
    }

    #[test]
    fn facade_aggregation_and_epilog_are_byte_exact() {
        let tree = tree();
        let first = context(1.0e-5);
        let second = context(1.0e-20);
        let mut format = TaxonFormat::new(true);
        print_match(&mut format, &first, &[11], &tree);
        print_match(&mut format, &second, &[20], &tree);

        assert_eq!(format.taxid, 2);
        assert_eq!(format.evalue, 1.0e-20);
        let mut output = Vec::new();
        print_query_epilog(&format, &mut output, &first.query_title, &tree).unwrap();
        assert_eq!(output, b"query\t2\t1.00e-20\tBacteria\n");
    }

    #[test]
    fn clone_copies_accumulator_state_and_reset_restores_constructor_values() {
        let tree = tree();
        let mut format = TaxonFormat::new(false);
        print_match(&mut format, &context(0.25), &[10], &tree);
        let cloned = format.clone();
        assert_eq!(cloned.taxid, 10);
        assert_eq!(cloned.evalue, 0.25);

        format.reset();
        assert_eq!(format.taxid, 0);
        assert_eq!(format.evalue, f64::MAX);
    }
}
