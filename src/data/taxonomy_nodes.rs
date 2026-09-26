//! Taxonomy node storage translated from
//! `diamond/src/data/taxonomy_nodes.{h,cpp}`.

use crate::basic::value::TaxId;
use crate::data::blastdb::taxdmp::read_nodes_dmp;
use crate::util::io::{Deserializer, IoResult, Serializer, StreamEntity};
use std::collections::BTreeSet;

/// Stable on-disk taxonomic rank values used by DIAMOND databases.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Rank {
    None = 0,
    Superkingdom = 1,
    CellularRoot = 2,
    AcellularRoot = 3,
    Domain = 4,
    Realm = 5,
    Kingdom = 6,
    Subkingdom = 7,
    Superphylum = 8,
    Phylum = 9,
    Subphylum = 10,
    Superclass = 11,
    ClassRank = 12,
    Subclass = 13,
    Infraclass = 14,
    Cohort = 15,
    Subcohort = 16,
    Superorder = 17,
    Order = 18,
    Suborder = 19,
    Infraorder = 20,
    Parvorder = 21,
    Superfamily = 22,
    Family = 23,
    Subfamily = 24,
    Tribe = 25,
    Subtribe = 26,
    Genus = 27,
    Subgenus = 28,
    Section = 29,
    Subsection = 30,
    Series = 31,
    SpeciesGroup = 32,
    SpeciesSubgroup = 33,
    Species = 34,
    Subspecies = 35,
    Varietas = 36,
    Forma = 37,
    Strain = 38,
    Biotype = 39,
    Clade = 40,
    FormaSpecialis = 41,
    Genotype = 42,
    Isolate = 43,
    Morph = 44,
    Pathogroup = 45,
    Serogroup = 46,
    Serotype = 47,
    Subvariety = 48,
}

impl Default for Rank {
    fn default() -> Self {
        Self::None
    }
}

impl Rank {
    pub const COUNT: usize = 49;
    pub const NAMES: [&'static str; Self::COUNT] = [
        "no rank",
        "superkingdom",
        "cellular root",
        "acellular root",
        "domain",
        "realm",
        "kingdom",
        "subkingdom",
        "superphylum",
        "phylum",
        "subphylum",
        "superclass",
        "class",
        "subclass",
        "infraclass",
        "cohort",
        "subcohort",
        "superorder",
        "order",
        "suborder",
        "infraorder",
        "parvorder",
        "superfamily",
        "family",
        "subfamily",
        "tribe",
        "subtribe",
        "genus",
        "subgenus",
        "section",
        "subsection",
        "series",
        "species group",
        "species subgroup",
        "species",
        "subspecies",
        "varietas",
        "forma",
        "strain",
        "biotype",
        "clade",
        "forma specialis",
        "genotype",
        "isolate",
        "morph",
        "pathogroup",
        "serogroup",
        "serotype",
        "subvariety",
    ];

    /// Matches the integer constructor in `taxonomy_nodes.h`.
    pub fn new(index: usize) -> Self {
        Self::from_index(index).expect("Invalid taxonomic rank index")
    }

    /// Matches `Rank::Rank(const char*)` with a recoverable Rust error.
    pub fn parse(name: &str) -> Result<Self, String> {
        Self::predefined(name)
            .and_then(Self::from_index)
            .ok_or_else(|| format!("Invalid taxonomic rank: {name}"))
    }

    /// Matches `Rank::predefined`; `None` corresponds to C++'s `-1`.
    pub fn predefined(name: &str) -> Option<usize> {
        Self::NAMES.iter().position(|&candidate| candidate == name)
    }

    pub fn name(self) -> &'static str {
        Self::NAMES[self as usize]
    }

    pub fn from_index(index: usize) -> Option<Self> {
        Some(match index {
            0 => Self::None,
            1 => Self::Superkingdom,
            2 => Self::CellularRoot,
            3 => Self::AcellularRoot,
            4 => Self::Domain,
            5 => Self::Realm,
            6 => Self::Kingdom,
            7 => Self::Subkingdom,
            8 => Self::Superphylum,
            9 => Self::Phylum,
            10 => Self::Subphylum,
            11 => Self::Superclass,
            12 => Self::ClassRank,
            13 => Self::Subclass,
            14 => Self::Infraclass,
            15 => Self::Cohort,
            16 => Self::Subcohort,
            17 => Self::Superorder,
            18 => Self::Order,
            19 => Self::Suborder,
            20 => Self::Infraorder,
            21 => Self::Parvorder,
            22 => Self::Superfamily,
            23 => Self::Family,
            24 => Self::Subfamily,
            25 => Self::Tribe,
            26 => Self::Subtribe,
            27 => Self::Genus,
            28 => Self::Subgenus,
            29 => Self::Section,
            30 => Self::Subsection,
            31 => Self::Series,
            32 => Self::SpeciesGroup,
            33 => Self::SpeciesSubgroup,
            34 => Self::Species,
            35 => Self::Subspecies,
            36 => Self::Varietas,
            37 => Self::Forma,
            38 => Self::Strain,
            39 => Self::Biotype,
            40 => Self::Clade,
            41 => Self::FormaSpecialis,
            42 => Self::Genotype,
            43 => Self::Isolate,
            44 => Self::Morph,
            45 => Self::Pathogroup,
            46 => Self::Serogroup,
            47 => Self::Serotype,
            48 => Self::Subvariety,
            _ => return None,
        })
    }
}

impl std::fmt::Display for Rank {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

/// Parent and rank arrays indexed directly by NCBI taxonomy id.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TaxonomyNodes {
    pub(crate) parent: Vec<TaxId>,
    pub(crate) rank: Vec<Rank>,
}

impl TaxonomyNodes {
    /// Matches `TaxonomyNodes(const string&)`.
    pub fn from_nodes_dmp(file_name: &str) -> Result<Self, String> {
        let mut nodes = Self::default();
        let mut rank_error = None;
        read_nodes_dmp(file_name, |taxid, parent, rank_name| {
            if rank_error.is_some() {
                return;
            }
            if taxid < 0 {
                rank_error = Some(format!("Invalid negative taxon id: {taxid}"));
                return;
            }
            let rank = match Rank::parse(rank_name) {
                Ok(rank) => rank,
                Err(error) => {
                    rank_error = Some(error);
                    return;
                }
            };
            let index = taxid as usize;
            nodes.parent.resize(index + 1, 0);
            nodes.parent[index] = parent;
            nodes.rank.resize(index + 1, Rank::None);
            nodes.rank[index] = rank;
        })
        .map_err(|error| error.to_string())?;
        if let Some(error) = rank_error {
            return Err(error);
        }
        Ok(nodes)
    }

    /// Matches the database constructor. Builds before 131 contain no rank
    /// array; newer builds store one rank byte per parent entry.
    pub fn from_deserializer<S: StreamEntity>(
        input: &mut Deserializer<S>,
        db_build: u32,
    ) -> IoResult<Self> {
        let count = input.read_value::<u32>()? as usize;
        let mut parent = Vec::with_capacity(count);
        for _ in 0..count {
            parent.push(input.read_i32()?);
        }

        let mut rank = Vec::new();
        if db_build >= 131 {
            let mut raw = vec![0u8; parent.len()];
            input.read_exact(&mut raw)?;
            rank.reserve(raw.len());
            for value in raw {
                rank.push(Rank::from_index(value as usize).unwrap_or(Rank::None));
            }
        }
        Ok(Self { parent, rank })
    }

    /// Matches `TaxonomyNodes::save`'s binary representation.
    pub fn save<S: StreamEntity>(&self, out: &mut Serializer<S>) -> IoResult<()> {
        out.write_value(self.parent.len() as u32)?;
        for &parent in &self.parent {
            out.write_i32(parent)?;
        }
        for &rank in &self.rank {
            out.write_raw(&[rank as u8])?;
        }
        Ok(())
    }

    pub fn get_parent(&self, taxid: TaxId) -> Result<TaxId, String> {
        if taxid < 0 || taxid as usize >= self.parent.len() {
            return Err(format!("No taxonomy node found for taxon id {taxid}"));
        }
        Ok(self.parent[taxid as usize])
    }

    pub fn rank(&self, taxid: TaxId) -> i32 {
        if taxid < 0 || taxid as usize >= self.rank.len() {
            -1
        } else {
            self.rank[taxid as usize] as i32
        }
    }

    /// Test whether `query` or one of its ancestors is in `filter`.
    pub fn contained(&self, query: TaxId, filter: &BTreeSet<TaxId>) -> Result<bool, String> {
        const MAX_LINEAGE: usize = 64;
        if filter.contains(&1) {
            return Ok(true);
        }
        let mut parent = query;
        let mut depth = 0usize;
        while parent > 1 && !filter.contains(&parent) {
            parent = self.get_parent(parent)?;
            if parent <= 0 {
                return Ok(false);
            }
            depth += 1;
            if depth > MAX_LINEAGE {
                return Err("Path in taxonomy too long (contained).".to_string());
            }
        }
        Ok(parent > 1)
    }

    /// Vector overload of C++ `contained`.
    pub fn contained_any(&self, query: &[TaxId], filter: &BTreeSet<TaxId>) -> Result<bool, String> {
        if filter.contains(&1) {
            return Ok(true);
        }
        for &taxid in query {
            if self.contained(taxid, filter)? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    pub fn max(&self) -> TaxId {
        self.parent.len().saturating_sub(1) as TaxId
    }

    pub fn len(&self) -> usize {
        self.parent.len()
    }

    pub fn is_empty(&self) -> bool {
        self.parent.is_empty()
    }

    pub fn parents(&self) -> &[TaxId] {
        &self.parent
    }

    pub fn ranks(&self) -> &[Rank] {
        &self.rank
    }

    /// Counts computed by the original `save` implementation for its status
    /// report, exposed without coupling storage to global logging streams.
    pub fn rank_counts(&self) -> [usize; Rank::COUNT] {
        let mut counts = [0; Rank::COUNT];
        for &rank in &self.rank {
            counts[rank as usize] += 1;
        }
        counts
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::io::VecStream;
    use std::io::Write;

    fn temp_file(label: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-taxonomy-nodes-{label}-{}-{}.dmp",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    #[test]
    fn rank_names_indices_and_parsing_match_cpp() {
        assert_eq!(Rank::COUNT, 49);
        for (index, name) in Rank::NAMES.iter().enumerate() {
            let rank = Rank::parse(name).unwrap();
            assert_eq!(rank as usize, index);
            assert_eq!(rank.name(), *name);
            assert_eq!(Rank::predefined(name), Some(index));
        }
        assert_eq!(Rank::from_index(Rank::COUNT), None);
        assert_eq!(Rank::predefined("root"), None);
        assert_eq!(
            Rank::parse("root").unwrap_err(),
            "Invalid taxonomic rank: root"
        );
    }

    #[test]
    fn nodes_dmp_constructor_returns_invalid_rank_error() {
        let path = temp_file("invalid-rank");
        {
            let mut file = std::fs::File::create(&path).unwrap();
            writeln!(file, "1\t|\t1\t|\tno rank\t|\t").unwrap();
            writeln!(file, "2\t|\t1\t|\tnot a rank\t|\t").unwrap();
        }
        let error = TaxonomyNodes::from_nodes_dmp(path.to_str().unwrap()).unwrap_err();
        assert_eq!(error, "Invalid taxonomic rank: not a rank");
        std::fs::remove_file(path).unwrap();
    }

    #[test]
    fn save_load_layout_and_legacy_rank_gate_match_cpp() {
        let nodes = TaxonomyNodes {
            parent: vec![0, 1, 1, 2],
            rank: vec![Rank::None, Rank::None, Rank::Superkingdom, Rank::Species],
        };
        let mut output = Serializer::new(VecStream::new());
        nodes.save(&mut output).unwrap();
        let bytes = output.into_inner().unwrap();
        assert_eq!(&bytes.data()[..4], &4u32.to_ne_bytes());

        let mut input = Deserializer::new(VecStream::from_vec(bytes.data().to_vec()));
        assert_eq!(
            TaxonomyNodes::from_deserializer(&mut input, 131).unwrap(),
            nodes
        );

        let mut legacy = Deserializer::new(VecStream::from_vec(bytes.data().to_vec()));
        let legacy = TaxonomyNodes::from_deserializer(&mut legacy, 130).unwrap();
        assert_eq!(legacy.parents(), nodes.parents());
        assert!(legacy.ranks().is_empty());
    }

    #[test]
    fn parent_rank_containment_and_counts_cover_header_api() {
        let nodes = TaxonomyNodes {
            parent: vec![0, 1, 1, 2, 3, 0],
            rank: vec![
                Rank::None,
                Rank::None,
                Rank::Superkingdom,
                Rank::Phylum,
                Rank::Species,
                Rank::None,
            ],
        };
        assert_eq!(nodes.get_parent(4), Ok(3));
        assert_eq!(
            nodes.get_parent(-1).unwrap_err(),
            "No taxonomy node found for taxon id -1"
        );
        assert_eq!(nodes.rank(4), Rank::Species as i32);
        assert_eq!(nodes.rank(99), -1);
        assert_eq!(nodes.max(), 5);

        let filter = BTreeSet::from([2]);
        assert!(nodes.contained(4, &filter).unwrap());
        assert!(!nodes.contained(5, &filter).unwrap());
        assert!(nodes.contained_any(&[5, 4], &filter).unwrap());
        assert_eq!(nodes.rank_counts()[Rank::None as usize], 3);
        assert_eq!(nodes.rank_counts()[Rank::Species as usize], 1);
    }
}
