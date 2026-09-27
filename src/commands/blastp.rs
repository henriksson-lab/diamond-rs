use std::io::{self, BufWriter, Write};
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
use crate::search::seed_match::SeedMatch;
use crate::search::{parallel, sensitivity};
use crate::stats::cbs::CbsMode;
use crate::stats::score_matrix::{CutoffTable2D, ScoreMatrix};
use crate::util::algo::PatternMatcher;

#[inline]
fn stage2_query_bounds(query_len: usize, seed_pos: usize) -> (usize, usize) {
    let window = crate::dp::ungapped_window::UNGAPPED_WINDOW;
    (
        seed_pos.saturating_sub(window),
        seed_pos.saturating_add(window).min(query_len),
    )
}

#[derive(Clone, Copy)]
struct StoredHit {
    subject: u64,
    query_id: u32,
    seed_offset: u32,
    score: u16,
}

#[derive(Clone, Copy)]
struct CompactHit {
    subject: u64,
    seed_offset: u32,
    score: u16,
}

/// Consume one stage-1 partition and retain only candidates that pass the
/// ungapped window filter.  This is deliberately partition-local: retaining
/// all Hamming survivors was the multi-gigabyte RSS bottleneck on repetitive
/// real databases.
fn filter_partition_to_hits(
    matches: &[SeedMatch],
    queries: &[Vec<Letter>],
    db_block: &Block,
    cutoffs: &[i32],
    score_matrix: &ScoreMatrix,
    shape: &Shape,
    left_most_context: &LeftMostContext<'_>,
    first_shape: bool,
    chunked: bool,
    use_left_most_range: bool,
    skip_left_most: bool,
    index_chunks: usize,
    min_identities: u32,
) -> Vec<StoredHit> {
    let mut out = Vec::new();
    let ref_seq_data = db_block.seqs().data();
    let mut group_begin = 0usize;
    while group_begin < matches.len() {
        let query_id = matches[group_begin].query_id as usize;
        let q_pos = matches[group_begin].query_pos as usize;
        let mut group_end = group_begin + 1;
        while group_end < matches.len()
            && matches[group_end].query_id as usize == query_id
            && matches[group_end].query_pos as usize == q_pos
        {
            group_end += 1;
        }
        let Some(query) = queries.get(query_id) else {
            group_begin = group_end;
            continue;
        };
        let (q_start, q_end) = stage2_query_bounds(query.len(), q_pos);
        let window_left = q_pos - q_start;
        let window_clipped = q_end - q_start;
        let query_window = &query[q_start..q_end];
        let cutoff = cutoffs[query_id];
        let interval_mod = (q_pos % 32) as i32;
        let interval_overhang = (window_left as i32 - interval_mod).max(0) as usize;
        let left_q_start = q_start + interval_overhang;
        let left_seed_offset = window_left.saturating_sub(interval_overhang);
        if left_q_start >= q_end {
            group_begin = group_end;
            continue;
        }
        let left_len = q_end - left_q_start;

        for chunk in matches[group_begin..group_end].chunks(32) {
            let mut subject_starts = [0isize; 32];
            for (start, hit) in subject_starts.iter_mut().zip(chunk) {
                *start = db_block
                    .seqs()
                    .position(hit.ref_id as usize, hit.ref_pos as usize)
                    as isize
                    - window_left as isize;
            }
            let subject_storage;
            let mut subject_windows = [&[][..]; 32];
            if subject_starts[..chunk.len()]
                .iter()
                .all(|&start| start >= 0)
            {
                for (window, &start) in subject_windows.iter_mut().zip(&subject_starts) {
                    *window = &ref_seq_data[start as usize..];
                }
            } else {
                subject_storage = subject_starts[..chunk.len()]
                    .iter()
                    .map(|&start| {
                        (0..window_clipped)
                            .map(|n| {
                                let pos = start + n as isize;
                                if pos < 0 {
                                    crate::basic::value::DELIMITER_LETTER
                                } else {
                                    ref_seq_data
                                        .get(pos as usize)
                                        .copied()
                                        .unwrap_or(crate::basic::value::DELIMITER_LETTER)
                                }
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
                for (window, storage) in subject_windows.iter_mut().zip(&subject_storage) {
                    *window = storage;
                }
            }
            let mut scores = [i32::MAX; 32];
            if cutoff != 0 {
                crate::dp::simd_ungapped::window_ungapped_best_into(
                    query_window,
                    &subject_windows[..chunk.len()],
                    window_clipped,
                    score_matrix,
                    &mut scores[..chunk.len()],
                );
            }
            for (hit, &score) in chunk.iter().zip(&scores) {
                if score <= cutoff {
                    continue;
                }
                let subject = db_block
                    .seqs()
                    .position(hit.ref_id as usize, hit.ref_pos as usize);
                let left_subject_start =
                    subject as isize - window_left as isize + interval_overhang as isize;
                let subject_storage;
                let subject_window = if left_subject_start >= 0
                    && left_subject_start as usize + left_len <= ref_seq_data.len()
                {
                    &ref_seq_data
                        [left_subject_start as usize..left_subject_start as usize + left_len]
                } else {
                    subject_storage = (0..left_len)
                        .map(|n| {
                            let pos = left_subject_start + n as isize;
                            if pos < 0 {
                                crate::basic::value::DELIMITER_LETTER
                            } else {
                                ref_seq_data
                                    .get(pos as usize)
                                    .copied()
                                    .unwrap_or(crate::basic::value::DELIMITER_LETTER)
                            }
                        })
                        .collect::<Vec<_>>();
                    &subject_storage
                };
                let current_range = if use_left_most_range {
                    let partitions = seedp_count(10) as usize;
                    let chunk_size = partitions.div_ceil(index_chunks);
                    let partition = seed_partition(hit.seed, seedp_mask(10)) as usize;
                    let begin = (partition / chunk_size) * chunk_size;
                    let end = (begin + chunk_size).min(partitions);
                    Some(SeedPartitionRange::with_bounds(begin as u32, end as u32))
                } else {
                    None
                };
                if !skip_left_most
                    && !left_most_filter_with_range(
                        &query[left_q_start..q_end],
                        subject_window,
                        left_seed_offset as i32,
                        shape.length,
                        left_most_context,
                        first_shape,
                        shape,
                        cutoff,
                        chunked,
                        min_identities,
                        current_range,
                    )
                {
                    continue;
                }
                out.push(StoredHit {
                    subject: subject as u64,
                    query_id: hit.query_id,
                    seed_offset: hit.query_pos,
                    score: if score == i32::MAX {
                        u16::MAX
                    } else {
                        score.min(u16::MAX as i32) as u16
                    },
                });
            }
        }
        group_begin = group_end;
    }
    out
}

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

    use crate::basic::value::{MASK_LETTER, SEED_MASK};
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

    // Keep compact, alignment-ready sequence storage before motif masking.
    // This duplicates only the residue bytes (a few MB in the benchmark), and
    // lets partition-local stage 2 run while the record copies remain masked
    // for seed enumeration.  It replaces the former multi-GB pair retention.
    let mut stage2_queries: Vec<Vec<Letter>> = query_records
        .iter()
        .map(|record| record.sequence.clone())
        .collect();
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
    let hamming_restores = |motifs: &[Vec<crate::masking::motifs::MotifMaskEntry>]| {
        motifs
            .iter()
            .map(|saved| {
                saved
                    .iter()
                    .map(|entry| (entry.pos, entry.original))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    let db_hamming_restores = hamming_restores(&db_motif_saves);
    let query_hamming_restores = hamming_restores(&query_motif_saves);

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
    let ungapped_evalue = traits.ungapped_evalue as f64;
    let short_query_ungapped_cutoff = if ungapped_evalue > 0.0 {
        score_matrix.rawscore_int(25.0)
    } else {
        0
    };
    let ungapped_cutoffs: Vec<i32> = stage2_queries
        .iter()
        .map(|query| {
            crate::search::stage2::ungapped_cutoff(
                query.len() as i32,
                ungapped_evalue,
                60,
                short_query_ungapped_cutoff,
                false,
                |q| score_matrix.ungapped_cutoff(q as usize, ungapped_evalue),
                |q| score_matrix.ungapped_cutoff(q as usize, ungapped_evalue),
            )
        })
        .collect();
    // Left-most filtering consumes SEED_MASK bits.  Gather them first without
    // retaining joined pairs, then the search pass can discard rejected pairs
    // inside each partition instead of keeping them until query extension.
    let mut low_complexity_query_positions = Vec::new();
    for shape in &shapes {
        let complexity_cut = traits.seed_cut * std::f64::consts::LN_2 * shape.weight as f64;
        low_complexity_query_positions.extend(
            parallel::collect_low_complexity_positions_partitioned_min_query_len(
                &query_seqs,
                &db_seqs,
                shape,
                &reduction,
                complexity_cut,
                config.min_query_len,
            ),
        );
    }
    let template_len = shapes.iter().map(|shape| shape.length).max().unwrap_or(0);
    for (query, saved) in stage2_queries.iter_mut().zip(&query_motif_saves) {
        crate::masking::motifs::restore_motifs(query, saved, template_len);
    }
    for &(query_id, query_pos) in &low_complexity_query_positions {
        if let Some(letter) = stage2_queries
            .get_mut(query_id as usize)
            .and_then(|query| query.get_mut(query_pos as usize))
        {
            *letter |= SEED_MASK;
        }
    }
    // The masking prepass builds and releases full seed arrays. Return those
    // pages before the search pass so they do not inflate its RSS high-water
    // mark on large databases.
    trim_freed_heap_pages();

    let chunked = traits.index_chunks > 1;
    let use_left_most_range = chunked
        && (config.ext_chunk_size == 0 || config.ext_chunk_size <= 128)
        && config.max_target_seqs <= 25;
    let skip_left_most = config.sensitivity >= Sensitivity::VerySensitive;
    let mut hits_by_query: Vec<Vec<CompactHit>> =
        (0..stage2_queries.len()).map(|_| Vec::new()).collect();
    let mut retained_hit_count = 0usize;
    let mut raw_seed_matches = 0usize;
    for (shape_id, shape) in shapes.iter().enumerate() {
        let complexity_cut = traits.seed_cut * std::f64::consts::LN_2 * shape.weight as f64;
        let (shape_partitions, shape_raw_seed_matches) =
            parallel::map_seed_matches_partitioned_streaming_hamming_min_query_len(
                &query_seqs,
                &query_hamming_restores,
                &db_seqs,
                &db_hamming_restores,
                shape,
                &reduction,
                complexity_cut,
                config.min_query_len,
                traits.min_identities,
                |matches| {
                    filter_partition_to_hits(
                        matches,
                        &stage2_queries,
                        &db_block,
                        &ungapped_cutoffs,
                        &score_matrix,
                        shape,
                        &left_most_contexts[shape_id],
                        shape_id == 0,
                        chunked,
                        use_left_most_range,
                        skip_left_most,
                        traits.index_chunks as usize,
                        traits.min_identities,
                    )
                },
            );
        raw_seed_matches += shape_raw_seed_matches;
        for shape_hits in shape_partitions {
            retained_hit_count += shape_hits.len();
            for hit in shape_hits {
                if let Some(query_hits) = hits_by_query.get_mut(hit.query_id as usize) {
                    query_hits.push(CompactHit {
                        subject: hit.subject,
                        seed_offset: hit.seed_offset,
                        score: hit.score,
                    });
                }
            }
        }
    }
    // C++ stage0/stage2 preserves shape + seed-partition join emission order
    // inside a query. Only group by query so range building is cheap without
    // imposing a Rust-only target/position tie-breaker on ranking boundaries.
    eprintln!(
        "Seed matches: {} (raw) -> {} (left-most)",
        raw_seed_matches, retained_hit_count,
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
    for (record, saved) in db_records.iter_mut().zip(db_motif_saves.iter()) {
        crate::masking::motifs::restore_motifs(&mut record.sequence, saved, 0);
    }
    for (record, saved) in query_records.iter_mut().zip(query_motif_saves.iter()) {
        crate::masking::motifs::restore_motifs(&mut record.sequence, saved, template_len);
    }
    // C++ `Search::mask_seeds` sets the high bit on every query occurrence in
    // a joined seed group rejected by the low-complexity test. Stage 2's
    // left-most filter uses those bits to suppress alternative seed starts.
    // The join workers return the compact (query id, position) list so we can
    // apply it after releasing the immutable sequence views.
    for (query_id, query_pos) in low_complexity_query_positions {
        if let Some(letter) = query_records
            .get_mut(query_id as usize)
            .and_then(|record| record.sequence.get_mut(query_pos as usize))
        {
            *letter |= SEED_MASK;
        }
    }

    // Stage 2 used the pre-motif copies above; release both seed-enumeration
    // record storage and the temporary query copy before extension.
    drop(db_records);
    drop(stage2_queries);
    drop(db_motif_saves);
    drop(query_motif_saves);
    trim_freed_heap_pages();
    let db_ids = db_block.ids().map_err(io::Error::other)?;

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

    // Process each query in input order. Native blastp keeps lazy target
    // masking disabled, so the C++-style extension path can read the shared
    // target block without serializing all queries.
    let process_query =
        |(query_idx, query_rec): (usize, &fasta::FastaRecord)| -> io::Result<Vec<u8>> {
            let query = &query_rec.sequence;

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

            let query_hits = &hits_by_query[query_idx];
            let mut hits = Vec::with_capacity(query_hits.len());
            hits.extend(
                query_hits
                    .iter()
                    .map(|hit| Hit::with_score(0, hit.subject, hit.seed_offset, hit.score)),
            );
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
            // Format each completed match directly into its query buffer. The
            // match list is already in output order, so retaining a second
            // vector of copied HSP summaries only adds allocation and traffic.
            let mut buf: Vec<u8> = Vec::new();
            for m in matches {
                let Some(best_hsp) = m.hsps.first() else {
                    continue;
                };
                let target_id = m.target_block_id as usize;

                let hsp = OutputHsp {
                    score: best_hsp.score,
                    evalue: best_hsp.evalue,
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
                let target_title = std::str::from_utf8(db_ids.get(target_id)).map_err(|error| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid UTF-8 database title: {error}"),
                    )
                })?;
                format::write_tabular_row(
                    &mut buf,
                    &query_rec.id,
                    target_title,
                    &hsp,
                    &fields,
                    query.len() as i32,
                    db_block.seqs().length(target_id) as i32,
                )?;
            }
            Ok(buf)
        };
    let mut total_alignments = 0u64;
    // Preserve input order without retaining every query's formatted output.
    // A moderately sized batch gives Rayon enough work to balance variable
    // protein lengths while keeping buffered output proportional to threads.
    let output_batch_size = rayon::current_num_threads().max(1) * 32;
    for (batch_idx, query_batch) in query_records.chunks(output_batch_size).enumerate() {
        let query_begin = batch_idx * output_batch_size;
        let batch_output: Vec<Vec<u8>> = if rayon::current_num_threads() == 1 {
            query_batch
                .iter()
                .enumerate()
                .map(|(offset, record)| process_query((query_begin + offset, record)))
                .collect::<io::Result<Vec<_>>>()?
        } else {
            query_batch
                .par_iter()
                .enumerate()
                .map(|(offset, record)| process_query((query_begin + offset, record)))
                .collect::<io::Result<Vec<_>>>()?
        };
        for buf in batch_output {
            // Each line in buf corresponds to one alignment; count newlines.
            total_alignments += buf.iter().filter(|&&b| b == b'\n').count() as u64;
            writer.write_all(&buf)?;
        }
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
    fn stage2_window_keeps_anchor_relative_right_edge_when_left_clipped() {
        assert_eq!(stage2_query_bounds(232, 3), (0, 51));
        assert_eq!(stage2_query_bounds(232, 100), (52, 148));
        assert_eq!(stage2_query_bounds(120, 100), (52, 120));
    }

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
