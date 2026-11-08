//! Source registrations with no effect in the pinned search workflow.
//!
//! These options still pass the canonical argument parser's type arity and
//! workflow checks. Their primitive scalar readers do not reject numeric text;
//! preserving their no-effect behavior must not select a legacy backend.
//! Constructor initialization and active consumers remain outside these lists.
//! The second group names fields consumed only by getseq/legacy tools, plus
//! seed-freq whose sole Source access requires the absent FREQUENCY_MASKING macro.

pub(crate) fn contains(name: &str) -> bool {
    globally_unused(name) || unused_by_search(name) || other_workflows_only(name)
}

fn globally_unused(name: &str) -> bool {
    matches!(
        name,
        "--space-penalty" // only unrelated Aligner function parameters use this name
        | "--score-ratio" // only a local DNA mapping-quality variable uses this name
        | "--raw" // unrelated stream buffer locals; Config::raw is never read
        | "--file1" // Config::file1 is never read
        | "--file2" // Config::file2 is never read
        | "--neighbors" // unrelated graph locals; Config::neighbors is never read
        | "--match1" // match_file1
        | "--match2" // match_file2
        | "--neighborhood-score" // neighborhood_score
        | "--superblock" // superblock
        | "--log-subject" // log_subject
        | "--fetch-size" // fetch_size
        | "--target-fetch-size" // target_fetch_size
        | "--family-counts" // family_counts_file
        | "--sort-join" // sort_join
        | "--simple-freq" // simple_freq
        | "--freq-treshold" // freq_treshold
        | "--use-dataset-field" // use_dataset_field
        | "--no-dict" // no_dict
        | "--upgma-edge-limit" // upgma_edge_limit
        | "--tree" // tree_file
        | "--upgma-dist" // upgma_dist
        | "--upgma-input" // upgma_input
        | "--tantan-maxRepeatOffset" // tantan_maxRepeatOffset
        | "--tantan-ungapped" // tantan_ungapped
        | "--no-swipe-realign" // no_swipe_realign
        | "--min-realign-overhang" // min_realign_overhang
        | "--ext-yield" // ext_min_yield
        | "--full-sw-len" // full_sw_len
        | "--relaxed-evalue-factor" // relaxed_evalue_factor
        | "--fast-tsv" // fast_tsv
        | "--query-memory" // query_memory
        | "--memory-intervals" // memory_intervals
        | "--seed-hit-density" // seedhit_density
        | "--chunk-size-multiplier" // chunk_size_multiplier
        | "--no-ref-masking" // no_ref_masking
        | "--target-bias" // target_bias
        | "--deque_bucket_size" // deque_bucket_size
        | "--no-reorder" // no_reorder
        | "--key2" // key2
        | "--motif-mask-file" // motif_mask_file
        | "--bootstrap-block" // bootstrap_block
        | "--centroid-factor" // centroid_factor
        | "--timeout" // timeout
        | "--resume" // resume
        | "--reassign-overlap" // reassign_overlap
        | "--reassign-ratio" // reassign_ratio
        | "--reassign-max" // reassign_max
        | "--classic-band" // classic_band
        | "--no_8bit_extension" // no_8bit_extension
        | "--max-indirection" // max_indirection
        | "--promiscuous-seed-ratio" // promiscuous_seed_ratio
    )
}

fn unused_by_search(name: &str) -> bool {
    matches!(
        name,
        "--hard-masked" | "--seq" | "--reverse" | "--oid-list"
        // BasicConfig still validates both primitives before dispatch.
        // dbtype changes makedb only; extension selects the DNA extension
        // engine, while this routing predicate is for blastp/blastx.
        | "--dbtype"
        | "--extension"
        // Only the obsolete Seed class uses this global; modern search and
        // listseeds construct PackedSeed through Shape instead.
        | "--seed-weight"
        | "--seed-freq"
        // These globals are consumed only by the nucleotide branch. Registry
        // validation still controls their availability in each build profile.
        | "--zdrop" | "--repetition-cutoff" | "--chaining-out"
        | "--align-long-reads" | "--best-hsp-only" | "--chain-pen-gap-scale"
        | "--chain-pen-skip-scale" | "--penalty" | "--reward"
        | "--chain-align-cutoff" | "--min-chain-score" | "--max-overlap-extension"
        | "--zdrop-extension" | "--zdrop-global" | "--band-extension" | "--band-global"
    )
}

// These fields have real consumers in makedb, clustering, view, ROC or tool
// command branches. The ordinary protein/translated search call graph never
// enters those branches; Config performs no post-store validation of them.
// In particular auto_create returns a FastaFile, not the makedb command.
fn other_workflows_only(name: &str) -> bool {
    matches!(
        name,
        "--in"
            | "--taxonmap"
            | "--taxonnodes"
            | "--taxonnames"
            | "--cluster-steps"
            | "--round-coverage"
            | "--round-approx-id"
            | "--member-cover"
            | "--connected-component-depth"
            | "--no-reassign"
            | "--centroid-out"
            | "--edges"
            | "--edge-format"
            | "--symmetric"
            | "--clusters"
            | "--weighted-gvc"
            | "--strict-gvc"
            | "--forwardonly"
            | "--type"
            | "--query-count"
            | "--cluster-similarity"
            | "--cluster-graph-file"
            | "--cluster-restart"
            | "--mcl-expansion"
            | "--mcl-inflation"
            | "--mcl-chunk-size"
            | "--mcl-max-iterations"
            | "--mcl-sparsity-switch"
            | "--mcl-nonsymmetric"
            | "--mcl-stats"
            | "--cluster-algo"
            | "--graph-algo"
            | "--tsv-read-size"
            | "--output-hits"
            | "--no-forward-fp"
            | "--output-fp"
            | "--family-cap"
    )
}
