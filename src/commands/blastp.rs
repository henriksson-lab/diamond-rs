use std::io::{self, BufWriter, Write};
use std::ops::Range;
use std::path::Path;
use std::time::Instant;

use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::align::hsp::Match;
use crate::align::target::{extend as extend_targets, GappedScoreConfig};
use crate::align::ungapped::UngappedStageConfig;
use crate::basic::reduction::Reduction;
use crate::basic::seed::{seed_partition, seedp_count, seedp_mask};
use crate::basic::shape::Shape;
use crate::basic::statistics::Statistics;
use crate::basic::value::{Letter, SequenceType};
use crate::config::Sensitivity;
use crate::data::block::Block;
use crate::data::fasta;
use crate::data::seed_histogram::SeedPartitionRange;
use crate::dp::swipe::{Flags, HspValues};
use crate::masking::{MaskingAlgo, MaskingMode};
use crate::output::format::{self, FieldId, Hsp as OutputHsp};
use crate::search::hit::Hit;
use crate::search::left_most::{left_most_filter_with_range, Context as LeftMostContext};
use crate::search::{parallel, sensitivity};
use crate::stats::cbs::CbsMode;
use crate::stats::score_matrix::{CutoffTable2D, ScoreMatrix};
use crate::util::algo::PatternMatcher;

#[cfg(all(target_os = "linux", target_env = "gnu"))]
unsafe extern "C" {
    fn malloc_trim(pad: usize) -> i32;
}

fn trim_freed_heap_pages() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    unsafe {
        let _ = malloc_trim(0);
    }
}

/// Configuration for a blastp run.
pub struct BlastpConfig {
    pub query_files: Vec<String>,
    pub database: String,
    pub output: Option<String>,
    pub matrix: String,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub max_evalue: f64,
    pub max_target_seqs: i64,
    pub ext_chunk_size: i64,
    pub toppercent: Option<f64>,
    pub global_ranking_targets: i64,
    pub min_id: f64,
    pub threads: i32,
    pub outfmt: Vec<String>,
    pub sensitivity: Sensitivity,
    pub masking: MaskingMode,
    pub motif_masking: String,
    pub min_query_len: usize,
    pub query_cover: f64,
    pub subject_cover: f64,
    pub comp_based_stats: CbsMode,
    pub no_self_hits: bool,
    /// Ungapped extension x-drop threshold in bits (`--xdrop`).
    /// C++ default is 12.3 (`diamond/src/basic/config.cpp:427`) and converts
    /// it to a raw score with the active score matrix at startup.
    pub ungapped_xdrop_bits: f64,
}

/// Run a simplified blastp search.
///
/// This implements the basic pipeline:
/// 1. Load query and database sequences
/// 2. Extract seeds from both using shapes
/// 3. Find seed matches (hash join)
/// 4. Extend hits with ungapped x-drop
/// 5. Perform gapped Smith-Waterman alignment
/// 6. Filter by e-value and output
pub fn run(config: &BlastpConfig) -> io::Result<()> {
    let start = Instant::now();

    // Honor `--threads N`. Previously the value was parsed but ignored, so
    // rayon fell back to its global pool (= `RAYON_NUM_THREADS` env or core
    // count). `try_build_global` is a no-op once the global pool exists, so
    // first call wins — subsequent invocations (e.g. multiple blastp runs in
    // the same process) silently keep the original count. For a one-shot CLI
    // run this is correct; for library use callers should install their own
    // pool. `config.threads <= 0` keeps the default.
    if config.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads as usize)
            .build_global();
    }

    // Load database sequences - try DMND first, then FASTA.
    // C++ `auto_append_extension_if_exists` (`config.cpp:770`) APPENDS `.dmnd`
    // only when the file as-given doesn't exist — it never STRIPS an existing
    // extension. So `--db nr.fasta` resolves to `nr.fasta` in C++ (then parsed
    // as FASTA), but `--db nr.fasta` was resolving to `nr.dmnd` in Rust because
    // `with_extension("dmnd")` REPLACES rather than appends. Fix: only consult
    // a `<db>.dmnd` companion when the given path doesn't exist.
    let (mut db_records, _db_from_dmnd) = {
        let db_path = Path::new(&config.database);
        let dmnd_path = if db_path.extension().is_some_and(|e| e == "dmnd") {
            db_path.to_path_buf()
        } else if db_path.exists() {
            // File-as-given exists (typically a FASTA). Don't try to substitute
            // a `.dmnd` sibling — leave the path alone so the FASTA branch
            // below runs. We set `dmnd_path` to the given path; `exists()` is
            // true but it's not a DMND file, so the DMND read will fail format
            // detection and we fall through to FASTA. To avoid that wasted
            // open, mark with a non-existent suffix instead.
            let mut p = db_path.as_os_str().to_owned();
            p.push(".dmnd");
            std::path::PathBuf::from(p)
        } else {
            // No file as-given — append `.dmnd` (don't replace).
            let mut p = db_path.as_os_str().to_owned();
            p.push(".dmnd");
            std::path::PathBuf::from(p)
        };

        if dmnd_path.exists() {
            let (header, records) = crate::data::dmnd_reader::read_dmnd(&dmnd_path)?;
            eprintln!(
                "Database: {} sequences, {} letters (DMND)",
                header.sequences, header.letters
            );
            (records, true)
        } else {
            let fasta_path = if db_path.extension().is_none() {
                db_path.with_extension("faa")
            } else {
                db_path.to_path_buf()
            };
            let records = fasta::read_fasta_file(&fasta_path, SequenceType::AminoAcid)?;
            eprintln!("Database: {} sequences (FASTA)", records.len());
            (records, false)
        }
    };

    // Compute total database letters for E-value normalization
    let db_letters: u64 = db_records.iter().map(|r| r.sequence.len() as u64).sum();

    // Load scoring matrix with database size
    let score_matrix = ScoreMatrix::new(
        &config.matrix,
        config.gap_open,
        config.gap_extend,
        0,
        1,
        db_letters,
    )
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    // Load query sequences
    let mut query_records = Vec::new();
    for qf in &config.query_files {
        let records = fasta::read_fasta_file(Path::new(qf), SequenceType::AminoAcid)?;
        query_records.extend(records);
    }
    eprintln!("Queries: {} sequences", query_records.len());

    use rayon::iter::IntoParallelRefMutIterator;
    match config.masking {
        MaskingMode::None => {}
        MaskingMode::Tantan => {
            // C++ constructs the tantan masker from the active ScoreMatrix and
            // applies it at search time to both targets and queries
            // (`Masking::Masking(const ScoreMatrix&)`, then `mask_seqs` in
            // `run/double_indexed.cpp`). This matters for non-default matrices
            // such as BLOSUM45: reusing makedb's default BLOSUM62 mask leaves
            // low-complexity regions under-masked and inflates self scores.
            let tantan_masker =
                crate::masking::tantan::TantanMasker::from_score_matrix(&score_matrix, 0.9);
            db_records.par_iter_mut().for_each(|r| {
                crate::masking::remove_bit_mask(&mut r.sequence);
                tantan_masker.mask(&mut r.sequence);
            });
            query_records.par_iter_mut().for_each(|r| {
                crate::masking::remove_bit_mask(&mut r.sequence);
                tantan_masker.mask(&mut r.sequence);
            });
        }
        MaskingMode::BlastSeg => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "native blastp does not implement --masking seg; use --legacy",
            ));
        }
    }

    // Convert tantan soft masks (SEED_MASK bit) to hard masks (MASK_LETTER = X).
    // C++ blastp calls `mask_seqs(..., hard_mask=true)` on both queries and
    // targets (run/double_indexed.cpp:127 and :719), which `Masking::operator()`
    // dispatches to `tantan::mask` in mode 1 — REPLACING masked letters with
    // `value_traits.mask_char` (= MASK_LETTER = 23) rather than OR-ing in the
    // high bit. Downstream the SIMD score profile then scores those positions
    // as `BLOSUM62(X, X) = -1`, not 0.
    //
    // Keeping the soft mask here makes our DP zero those positions instead of
    // applying the -1 BLOSUM penalty — for a self-self of a heavily-masked
    // protein (Q8QZQ8: 120+ masked residues), that's 120+ extra score relative
    // to C++ and shifts which targets fit inside `-k 25`.
    use crate::basic::value::{MASK_LETTER, SEED_MASK};
    let hard_mask = |seq: &mut [Letter]| {
        for l in seq.iter_mut() {
            if *l & SEED_MASK != 0 {
                *l = MASK_LETTER;
            }
        }
    };
    if config.masking == MaskingMode::Tantan {
        db_records
            .par_iter_mut()
            .for_each(|r| hard_mask(&mut r.sequence));
        query_records
            .par_iter_mut()
            .for_each(|r| hard_mask(&mut r.sequence));
    }

    // Motif masking — ports C++ `Block::soft_mask(MOTIF)` invoked from
    // `enum_seeds` (enum_seeds.h:202). At default sensitivity DIAMOND
    // hard-masks any 8-letter window matching a curated motif before seed
    // enumeration; the seed iterator then drops seeds overlapping a motif.
    // After enumeration C++ restores the original letters via
    // `Block::remove_soft_masking` so alignment runs against the real
    // sequence — we mirror that with `restore_motifs` further down.
    //
    // Run in parallel across sequences with rayon. C++ also parallelises
    // motif masking via `mask_seqs` (`masking/masking.cpp:172`).
    let soft_masking = sensitivity::soft_masking_algo(
        &sensitivity::get_traits(config.sensitivity),
        &config.motif_masking,
        false,
        false,
    )
    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let use_motif_masking = soft_masking == MaskingAlgo::Motif;
    let db_motif_saves: Vec<Vec<crate::masking::motifs::MotifMaskEntry>> = if use_motif_masking {
        db_records
            .par_iter_mut()
            .map(|r| crate::masking::motifs::mask_motifs(&mut r.sequence))
            .collect()
    } else {
        vec![Vec::new(); db_records.len()]
    };
    let query_motif_saves: Vec<Vec<crate::masking::motifs::MotifMaskEntry>> = if use_motif_masking {
        query_records
            .par_iter_mut()
            .map(|r| crate::masking::motifs::mask_motifs(&mut r.sequence))
            .collect()
    } else {
        vec![Vec::new(); query_records.len()]
    };

    // Set up seed extraction using sensitivity-appropriate shapes
    let reduction = Reduction::default_reduction();
    let shape_codes = sensitivity::get_shape_codes(config.sensitivity);
    let shapes: Vec<Shape> = shape_codes
        .iter()
        .map(|code| Shape::from_code(code, &reduction))
        .collect();
    let shape_patterns: Vec<u32> = shapes.iter().map(|shape| shape.mask).collect();
    let left_most_contexts: Vec<LeftMostContext<'_>> = (0..shapes.len())
        .map(|sid| {
            let previous = if sid == 0 {
                &shape_patterns[0..0]
            } else {
                &shape_patterns[0..sid]
            };
            LeftMostContext {
                previous_matcher: PatternMatcher::new(previous),
                current_matcher: PatternMatcher::new(&shape_patterns[0..=sid]),
                short_query_ungapped_cutoff: score_matrix.rawscore_int(25.0),
                seedp_mask: seedp_mask(10),
                reduction: &reduction,
            }
        })
        .collect();
    eprintln!(
        "Sensitivity: {:?}, shapes: {} (weights: {})",
        config.sensitivity,
        shapes.len(),
        shapes
            .iter()
            .map(|s| s.weight.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );

    // Build partitioned seed arrays and join for each shape (parallel)
    let db_seqs: Vec<&[Letter]> = db_records.iter().map(|r| r.sequence.as_slice()).collect();
    let query_seqs: Vec<&[Letter]> = query_records
        .iter()
        .map(|r| r.sequence.as_slice())
        .collect();
    // Per-sensitivity seed filters matching C++ DIAMOND:
    //   - `seed_cut * ln(2) * shape.weight` is the entropy floor for
    //     `seed_is_complex` (search/setup.cpp).
    //   - Frequent-seed masking is guarded by C++ `config.freq_masking`; the
    //     native default path leaves it disabled.
    let traits = sensitivity::get_traits(config.sensitivity);
    // Process shapes in order and let each shape's seed-array build/join use
    // the Rayon pool internally. Running the outer shape loop in parallel
    // competes with the inner reference seed-array work and scales worse on
    // small real query blocks; collecting in shape order preserves the previous
    // output order.
    let mut all_seed_matches = Vec::new();
    for (shape_id, shape) in shapes.iter().enumerate() {
        let complexity_cut = traits.seed_cut * std::f64::consts::LN_2 * shape.weight as f64;
        let mut matches = parallel::find_seed_matches_partitioned_filtered_min_query_len(
            &query_seqs,
            &db_seqs,
            shape,
            &reduction,
            complexity_cut,
            0.0,
            config.min_query_len,
        );
        for m in &mut matches {
            m.shape_id = shape_id as u32;
        }
        all_seed_matches.append(&mut matches);
    }
    let raw_seed_matches = all_seed_matches.len();
    // C++ stage0/stage2 preserves shape + seed-partition join emission order
    // inside a query. Only group by query so range building is cheap without
    // imposing a Rust-only target/position tie-breaker on ranking boundaries.
    all_seed_matches.sort_by_key(|m| m.query_id);

    // Stage-1 hamming filter: drop seed hits whose 48-letter window matches
    // fewer than `min_identities` letters. Ports C++ search/hamming/kernel.h.
    all_seed_matches = crate::search::hamming_filter::apply_hamming_filter(
        all_seed_matches,
        &query_seqs,
        &db_seqs,
        traits.min_identities,
    );
    eprintln!(
        "Seed matches: {} (raw) -> {} (hamming)",
        raw_seed_matches,
        all_seed_matches.len(),
    );

    // Restore motif-masked positions before alignment so the gapped extension
    // scores against the real letters (C++ `Block::remove_soft_masking`).
    // The earlier `db_seqs`/`query_seqs` slices are dropped here so we can
    // re-borrow the underlying records mutably.
    {
        let _ = &db_seqs;
        let _ = &query_seqs;
    }
    drop(db_seqs);
    drop(query_seqs);
    trim_freed_heap_pages();
    // C++ `enum_seeds.h:223` calls `seqs.remove_soft_masking(template_len, mask_seeds)`
    // with `mask_seeds=true` only on the query side, propagating SEED_MASK over
    // `[motif_begin - template_len + 1, motif_begin + motif_len)`. The DB-side
    // restore (stage0.cpp:142-144 with `mask_seeds=false`) just rewrites
    // letters. `template_len` is `max(shape.length_)` — derive from active shapes.
    let template_len = shapes.iter().map(|s| s.length).max().unwrap_or(0);
    for (record, saved) in db_records.iter_mut().zip(db_motif_saves.iter()) {
        crate::masking::motifs::restore_motifs(&mut record.sequence, saved, 0);
    }
    for (record, saved) in query_records.iter_mut().zip(query_motif_saves.iter()) {
        crate::masking::motifs::restore_motifs(&mut record.sequence, saved, template_len);
    }

    let mut db_block = Block::new();
    for (idx, record) in db_records.iter().enumerate() {
        db_block
            .push_back(
                &record.sequence,
                Some(&record.id),
                None,
                idx as u64,
                SequenceType::AminoAcid,
                0,
                false,
            )
            .map_err(io::Error::other)?;
    }

    // Parse output format. Only tabular is implemented in the native pipeline.
    // For PAF/SAM/XML/pairwise/DAA the user must use --legacy.
    let fields = if config.outfmt.is_empty() || config.outfmt[0] == "6" || config.outfmt[0] == "tab"
    {
        if config.outfmt.len() > 1 {
            config.outfmt[1..]
                .iter()
                .filter_map(|f| FieldId::from_name(f))
                .collect()
        } else {
            format::DEFAULT_TABULAR_FIELDS.to_vec()
        }
    } else {
        return Err(io::Error::other(format!(
            "Output format '{}' is not implemented in the native blastp pipeline. \
             Supported: 6/tab. Use --legacy to route to the C++ engine for other formats.",
            config.outfmt[0]
        )));
    };

    // Set up output writer
    let output: Box<dyn Write> = match &config.output {
        Some(path) => Box::new(BufWriter::new(std::fs::File::create(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut writer = output;
    const GAPPED_FILTER_EVALUE1: f64 = 2000.0;
    const GAPPED_FILTER_DIAG_BITS: f64 = 12.0;
    const GAPPED_FILTER_WINDOW: i32 = 200;

    let gapped_filter_evalue = traits.gapped_filter_evalue as f64;
    let cutoff_gapped1 = if gapped_filter_evalue != 0.0 {
        Some(CutoffTable2D::new(&score_matrix, GAPPED_FILTER_EVALUE1))
    } else {
        None
    };
    let cutoff_gapped2 = if gapped_filter_evalue != 0.0 {
        Some(CutoffTable2D::new(&score_matrix, gapped_filter_evalue))
    } else {
        None
    };
    let gapped_filter_diag_score = score_matrix.rawscore_int(GAPPED_FILTER_DIAG_BITS);

    // `all_seed_matches` is sorted by query id before the hamming filter, and
    // the filter preserves order. Store dense ranges instead of hashing each
    // query into a small Vec of references.
    let mut query_match_ranges: Vec<Range<usize>> = vec![0..0; query_records.len()];
    let mut i = 0usize;
    while i < all_seed_matches.len() {
        let query_id = all_seed_matches[i].query_id as usize;
        let begin = i;
        while i < all_seed_matches.len() && all_seed_matches[i].query_id as usize == query_id {
            i += 1;
        }
        if query_id < query_match_ranges.len() {
            query_match_ranges[query_id] = begin..i;
        }
    }

    // Process each query in input order. Native blastp keeps lazy target
    // masking disabled, so the C++-style extension path can read the shared
    // target block without serializing all queries.
    let process_query =
        |(query_idx, query_rec): (usize, &fasta::FastaRecord)| -> io::Result<Vec<u8>> {
            let query = &query_rec.sequence;

            // Compute ungapped score cutoff for this query length.
            // C++ `CutoffTable` (`util/scores/cutoff_table.h`) uses the NORMALIZED
            // 1e9-letter database size via `bitscore_norm`/`rawscore`, NOT the
            // actual DB size — the cutoff is a property of the query length only.
            //
            // The evalue threshold comes from sensitivity traits:
            //   - Faster/Fast/Shapes6x10/Shapes30x10/Linclust*: 0 (no filter)
            //   - Default/MidSensitive/Sensitive/MoreSensitive: 10000
            //   - VerySensitive: 100000
            //   - UltraSensitive: 300000
            // The previous hardcoded `10000.0` matched only Default and lost
            // recall on VerySensitive/UltraSensitive runs and over-filtered
            // Faster/Fast/Shapes paths.
            let ungapped_evalue = traits.ungapped_evalue as f64;
            // C++ `Search::ungapped_cutoff` (`diamond/src/search/stage2.h:39-54`)
            // dispatches on query length: for `query_len <= short_query_max_len`
            // (default 60) it returns a fixed cutoff derived from
            // `short_query_ungapped_bitscore` (default 25.0 bits); only longer
            // queries go through the 1e9-letter `CutoffTable`. Previously this
            // call site went straight to `cutoff_table`, over-filtering very
            // short queries. The `query_translated && qlen <= 85` branch is
            // blastx-only.
            let short_query_ungapped_cutoff = if ungapped_evalue > 0.0 {
                score_matrix.rawscore_int(25.0)
            } else {
                0
            };
            let ungapped_cutoff = crate::search::stage2::ungapped_cutoff(
                query.len() as i32,
                ungapped_evalue,
                60,
                short_query_ungapped_cutoff,
                false,
                |q| score_matrix.ungapped_cutoff(q as usize, ungapped_evalue),
                |q| score_matrix.ungapped_cutoff(q as usize, ungapped_evalue),
            );

            // CBS (composition-based statistics) correction per query position.
            // C++ default is comp-based-stats=1 (Hauser correction, window=40).
            // `Sequence::operator[]` strips the soft-mask bit while preserving
            // the residue, so compute Hauser from the query as stored rather
            // than converting masked positions to X.
            let query_cbs = if config.comp_based_stats.hauser() {
                crate::stats::cbs::hauser_correction(query, &score_matrix)
            } else {
                Vec::new()
            };
            let query_comp = crate::stats::cbs::compute_composition(query);
            let ungapped_cfg = UngappedStageConfig {
                comp_based_stats: config.comp_based_stats,
                xdrop: score_matrix.rawscore_int(config.ungapped_xdrop_bits),
                ..UngappedStageConfig::default()
            };
            let ext_mode = sensitivity::default_ext_mode(config.sensitivity);
            let gapped_cfg = GappedScoreConfig {
                comp_based_stats_hauser: config.comp_based_stats.hauser(),
                comp_based_stats_matrix_adjust: config.comp_based_stats.matrix_adjust(),
                query_cover: config.query_cover,
                subject_cover: config.subject_cover,
                no_self_hits: config.no_self_hits,
                max_evalue: config.max_evalue,
                min_id: config.min_id,
                max_target_seqs: config.max_target_seqs,
                ext_chunk_size: config.ext_chunk_size,
                toppercent: config.toppercent,
                global_ranking_targets: config.global_ranking_targets,
                gapped_filter_evalue,
                sensitivity: config.sensitivity,
                ..GappedScoreConfig::default()
            };

            // Get seed matches for this query.
            let query_seed_matches = &all_seed_matches[query_match_ranges[query_idx].clone()];

            // Stage-2 filter: 96-letter ungapped window walk (C++
            // `dp/ungapped_align.cpp:ungapped_window`). Surviving hits are
            // passed to the same batched extension/ranking path used by the
            // lower-level Rust translation of C++ `Extension::extend`.
            let mut hits = Vec::new();
            hits.reserve(query_seed_matches.len());
            let ref_seq_data = db_block.seqs().data();
            let mut group_begin = 0usize;
            while group_begin < query_seed_matches.len() {
                let sid = query_seed_matches[group_begin].shape_id as usize;
                let q_pos = query_seed_matches[group_begin].query_pos as usize;
                let mut group_end = group_begin + 1;
                while group_end < query_seed_matches.len()
                    && query_seed_matches[group_end].shape_id as usize == sid
                    && query_seed_matches[group_end].query_pos as usize == q_pos
                {
                    group_end += 1;
                }

                let q_start = q_pos.saturating_sub(crate::dp::ungapped_window::UNGAPPED_WINDOW);
                let window_left = q_pos - q_start;
                let q_end =
                    (q_start + crate::dp::ungapped_window::UNGAPPED_WINDOW * 2).min(query.len());
                let window_clipped = q_end - q_start;
                let shape = &shapes[sid];
                let interval_mod = (q_pos % 32) as i32;
                let interval_overhang = (window_left as i32 - interval_mod).max(0) as usize;
                let left_q_start = q_start + interval_overhang;
                let left_seed_offset = window_left.saturating_sub(interval_overhang);
                if left_q_start >= q_end {
                    group_begin = group_end;
                    continue;
                }
                let left_len = q_end - left_q_start;
                let query_window = &query[q_start..q_end];

                let mut chunk_begin = group_begin;
                while chunk_begin < group_end {
                    let chunk_end = (chunk_begin + 32).min(group_end);
                    let chunk = &query_seed_matches[chunk_begin..chunk_end];
                    let subjects = chunk
                        .iter()
                        .map(|hit| {
                            db_block
                                .seqs()
                                .position(hit.ref_id as usize, hit.ref_pos as usize)
                                as isize
                                - window_left as isize
                        })
                        .collect::<Vec<_>>();

                    let subject_storage;
                    let subject_windows = if subjects.iter().all(|&start| start >= 0) {
                        subjects
                            .iter()
                            .map(|&start| &ref_seq_data[start as usize..])
                            .collect::<Vec<_>>()
                    } else {
                        subject_storage = subjects
                            .iter()
                            .map(|&start| {
                                (0..window_clipped)
                                    .map(|n| {
                                        let si = start + n as isize;
                                        if si < 0 {
                                            crate::basic::value::DELIMITER_LETTER
                                        } else {
                                            ref_seq_data
                                                .get(si as usize)
                                                .copied()
                                                .unwrap_or(crate::basic::value::DELIMITER_LETTER)
                                        }
                                    })
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>();
                        subject_storage
                            .iter()
                            .map(Vec::as_slice)
                            .collect::<Vec<_>>()
                    };
                    let scores = if ungapped_cutoff == 0 {
                        vec![i32::MAX; chunk.len()]
                    } else {
                        crate::dp::simd_ungapped::window_ungapped_best(
                            query_window,
                            &subject_windows,
                            window_clipped,
                            &score_matrix,
                        )
                    };

                    for (chunk_idx, hit) in chunk.iter().enumerate() {
                        let score = scores[chunk_idx];
                        if score <= ungapped_cutoff {
                            continue;
                        }
                        let subject = db_block
                            .seqs()
                            .position(hit.ref_id as usize, hit.ref_pos as usize);
                        let left_subject_start = subjects[chunk_idx] + interval_overhang as isize;
                        let subject_storage;
                        let subject_window = if left_subject_start >= 0
                            && (left_subject_start as usize + left_len) <= ref_seq_data.len()
                        {
                            &ref_seq_data[left_subject_start as usize
                                ..left_subject_start as usize + left_len]
                        } else {
                            subject_storage = (0..left_len)
                                .map(|n| {
                                    let si = left_subject_start + n as isize;
                                    if si < 0 {
                                        crate::basic::value::DELIMITER_LETTER
                                    } else {
                                        ref_seq_data
                                            .get(si as usize)
                                            .copied()
                                            .unwrap_or(crate::basic::value::DELIMITER_LETTER)
                                    }
                                })
                                .collect::<Vec<_>>();
                            &subject_storage
                        };
                        let chunked = traits.index_chunks > 1;
                        let use_left_most_range = chunked
                            && (config.ext_chunk_size == 0 || config.ext_chunk_size <= 128)
                            && config.max_target_seqs <= 25;
                        let current_range = if use_left_most_range {
                            let partitions = seedp_count(10) as usize;
                            let chunk_size = partitions.div_ceil(traits.index_chunks as usize);
                            let partition = seed_partition(hit.seed, seedp_mask(10)) as usize;
                            let begin = (partition / chunk_size) * chunk_size;
                            let end = (begin + chunk_size).min(partitions);
                            Some(SeedPartitionRange::with_bounds(begin as u32, end as u32))
                        } else {
                            None
                        };
                        let skip_left_most = config.sensitivity >= Sensitivity::VerySensitive;
                        if !skip_left_most
                            && !left_most_filter_with_range(
                                &query[left_q_start..q_end],
                                subject_window,
                                left_seed_offset as i32,
                                shape.length,
                                &left_most_contexts[sid],
                                sid == 0,
                                shape,
                                ungapped_cutoff,
                                chunked,
                                traits.min_identities,
                                current_range,
                            )
                        {
                            continue;
                        }
                        hits.push(Hit::with_score(
                            0,
                            subject as u64,
                            hit.query_pos as u32,
                            if score == i32::MAX {
                                u16::MAX
                            } else {
                                score.min(u16::MAX as i32) as u16
                            },
                        ));
                    }
                    chunk_begin = chunk_end;
                }
                group_begin = group_end;
            }
            let mut stat = Statistics::new();
            let output_hsp_values = HspValues::COORDS
                | HspValues::IDENT
                | HspValues::LENGTH
                | HspValues::MISMATCHES
                | HspValues::GAP_OPENINGS;
            let query_cbs_arg: &[Vec<i8>] = if query_cbs.is_empty() {
                &[]
            } else {
                std::slice::from_ref(&query_cbs)
            };
            let matches = extend_targets(
                query_idx as u32,
                &mut hits,
                std::slice::from_ref(query),
                &query_rec.id,
                query.len() as i32,
                query_cbs_arg,
                &query_comp,
                &db_block,
                &mut stat,
                Flags::NONE,
                ext_mode,
                &gapped_cfg,
                &ungapped_cfg,
                &score_matrix,
                output_hsp_values,
                |query_len, target_len| {
                    cutoff_gapped1
                        .as_ref()
                        .map_or(-1, |table| table.call(query_len, target_len))
                },
                |query_len, target_len| {
                    cutoff_gapped2
                        .as_ref()
                        .map_or(-1, |table| table.call(query_len, target_len))
                },
                score_matrix.gap_open(),
                score_matrix.gap_extend(),
                gapped_filter_diag_score,
                GAPPED_FILTER_WINDOW,
                Option::<fn(u32, &crate::align::gapped_filter::SeedHitList) -> Vec<Match>>::None,
            );
            let mut target_results: Vec<(u32, f64, OutputHsp)> = Vec::new();
            for m in matches {
                let Some(best_hsp) = m.hsps.first() else {
                    continue;
                };
                let target_id = m.target_block_id;
                let evalue = best_hsp.evalue;

                let hsp = OutputHsp {
                    score: best_hsp.score,
                    evalue,
                    bit_score: best_hsp.bit_score,
                    query_range: (best_hsp.query_range.begin, best_hsp.query_range.end),
                    subject_range: (best_hsp.subject_range.begin, best_hsp.subject_range.end),
                    query_source_range: (best_hsp.query_range.begin, best_hsp.query_range.end),
                    subject_source_range: (
                        best_hsp.subject_range.begin,
                        best_hsp.subject_range.end,
                    ),
                    frame: best_hsp.frame,
                    length: best_hsp.length,
                    identities: best_hsp.identities,
                    mismatches: best_hsp.mismatches,
                    positives: best_hsp.positives,
                    gap_openings: best_hsp.gap_openings,
                    gaps: best_hsp.gaps,
                };

                target_results.push((target_id, evalue, hsp));
            }

            // Output results into a per-query buffer.
            let mut buf: Vec<u8> = Vec::new();
            for (target_id, _, hsp) in &target_results {
                format::write_tabular_row(
                    &mut buf,
                    &query_rec.id,
                    &db_records[*target_id as usize].id,
                    hsp,
                    &fields,
                    query.len() as i32,
                    db_records[*target_id as usize].sequence.len() as i32,
                )?;
            }
            Ok(buf)
        };
    let per_query_output: Vec<Vec<u8>> = if rayon::current_num_threads() == 1 {
        query_records
            .iter()
            .enumerate()
            .map(process_query)
            .collect::<io::Result<Vec<_>>>()?
    } else {
        query_records
            .par_iter()
            .enumerate()
            .map(process_query)
            .collect::<io::Result<Vec<_>>>()?
    };

    let mut total_alignments = 0u64;
    for buf in &per_query_output {
        // Each line in buf corresponds to one alignment; count newlines.
        let buf: &Vec<u8> = buf;
        total_alignments += buf.iter().filter(|&&b: &&u8| b == b'\n').count() as u64;
        writer.write_all(buf)?;
    }
    writer.flush()?;

    let elapsed = start.elapsed();
    eprintln!(
        "Reported {} alignments in {:.1}s",
        total_alignments,
        elapsed.as_secs_f64()
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blastp_with_dmnd() {
        let query = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/5.faa");
        let db = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/data.dmnd");
        let output_path = std::env::temp_dir().join("test_blastp_dmnd.out");

        let config = BlastpConfig {
            query_files: vec![query.to_string()],
            database: db.to_string(),
            output: Some(output_path.to_string_lossy().to_string()),
            matrix: "blosum62".to_string(),
            gap_open: 11,
            gap_extend: 1,
            max_evalue: 0.001,
            max_target_seqs: 25,
            ext_chunk_size: 0,
            toppercent: None,
            global_ranking_targets: 0,
            min_id: 0.0,
            threads: 1,
            outfmt: vec![],
            sensitivity: Sensitivity::Default,
            masking: MaskingMode::Tantan,
            motif_masking: String::new(),
            min_query_len: 0,
            query_cover: 0.0,
            subject_cover: 0.0,
            comp_based_stats: CbsMode::Hauser,
            no_self_hits: false,
            ungapped_xdrop_bits: 12.3,
        };

        let result = run(&config);
        assert!(
            result.is_ok(),
            "blastp with DMND failed: {:?}",
            result.err()
        );

        let output = std::fs::read_to_string(&output_path).unwrap();
        assert!(!output.is_empty(), "blastp produced no output");
        let _ = std::fs::remove_file(&output_path);
    }

    #[test]
    fn test_blastp_self() {
        let fasta_path = concat!(env!("CARGO_MANIFEST_DIR"), "/diamond/src/test/1.faa");
        let output_path = std::env::temp_dir().join("test_blastp_self.out");

        let config = BlastpConfig {
            query_files: vec![fasta_path.to_string()],
            database: fasta_path.to_string(),
            output: Some(output_path.to_string_lossy().to_string()),
            matrix: "blosum62".to_string(),
            gap_open: 11,
            gap_extend: 1,
            max_evalue: 0.001,
            max_target_seqs: 25,
            ext_chunk_size: 0,
            toppercent: None,
            global_ranking_targets: 0,
            min_id: 0.0,
            threads: 1,
            outfmt: vec![],
            sensitivity: Sensitivity::Default,
            masking: MaskingMode::Tantan,
            motif_masking: String::new(),
            min_query_len: 0,
            query_cover: 0.0,
            subject_cover: 0.0,
            comp_based_stats: CbsMode::Hauser,
            no_self_hits: false,
            ungapped_xdrop_bits: 12.3,
        };

        let result = run(&config);
        assert!(result.is_ok(), "blastp failed: {:?}", result.err());

        // Check output file exists and has content
        let output = std::fs::read_to_string(&output_path).unwrap();
        assert!(!output.is_empty(), "blastp produced no output");
        // Self-alignment should find at least one hit
        assert!(output.lines().count() >= 1);

        let _ = std::fs::remove_file(&output_path);
    }
}
