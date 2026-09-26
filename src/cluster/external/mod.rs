//! External-memory clustering compatibility hierarchy.

pub mod align;
pub mod cluster;
pub mod external;
pub mod output;
pub mod pair_table;
pub mod seed_table;

pub use align::{align, align_rep, AlignConfig, ChunkSequences};
pub use cluster::{
    cluster, cluster_bidirectional, compute_closure, compute_closure_from_assignment_file, get_reps,
};
pub use external::*;
pub use output::{
    merge, output, output_accs, output_accs_round1, output_oids, read_clustering, AccMapping,
    ExternalOutputConfig, ExternalOutputSummary,
};
pub use pair_table::{
    build_pair_table, get_pairs_mutual_cov, get_pairs_uni_cov, Bucket, PairFileArray,
    PairTableConfig, RadixedTable, SeedEntry, RADIX_BITS, RADIX_COUNT,
};
pub use seed_table::{
    build_seed_table, DefaultSeedMasker, FastaVolumeReader, SeedFileArray, SeedMasker,
    SeedSequence, SeedTableConfig, SeedTableStats, SequenceVolumeReader, Volume, VolumedFile,
};
