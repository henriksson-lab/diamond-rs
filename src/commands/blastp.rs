use std::collections::HashMap;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::align::gapped_filter::SeedHit;
use crate::align::target::{align_work_targets, GappedScoreConfig};
use crate::align::ungapped::{ungapped_stage_target, UngappedStageConfig};
use crate::basic::reduction::Reduction;
use crate::basic::shape::Shape;
use crate::basic::statistics::Statistics;
use crate::basic::value::{Letter, SequenceType};
use crate::config::Sensitivity;
use crate::data::block::Block;
use crate::data::fasta;
use crate::dp::swipe::{Flags, HspValues};
use crate::masking::{MaskingAlgo, MaskingMode};
use crate::output::format::{self, FieldId, Hsp as OutputHsp};
use crate::search::{parallel, seed_match, sensitivity};
use crate::stats::cbs::CbsMode;
use crate::stats::score_matrix::ScoreMatrix;

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
    let mut db_records = {
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
            records
        } else {
            let fasta_path = if db_path.extension().is_none() {
                db_path.with_extension("faa")
            } else {
                db_path.to_path_buf()
            };
            let records = fasta::read_fasta_file(&fasta_path, SequenceType::AminoAcid)?;
            eprintln!("Database: {} sequences (FASTA)", records.len());
            records
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
            // Apply tantan soft masking for seed filtering (high bit set on masked positions).
            // The DMND reader returns sequences that already have tantan marks (set during
            // makedb), so this is a no-op for the canonical workflow — but we keep the
            // call so the FASTA->blastp path (no precomputed marks) still masks.
            db_records.par_iter_mut().for_each(|r| {
                crate::masking::mask_sequence(&mut r.sequence, MaskingAlgo::Tantan);
            });
            query_records.par_iter_mut().for_each(|r| {
                crate::masking::mask_sequence(&mut r.sequence, MaskingAlgo::Tantan);
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
    //   - `freq_sd` (default 50.0 for blastp) is the multiplier for the
    //     mean+SD frequent-seed mask (data/frequent_seeds.cpp).
    let traits = sensitivity::get_traits(config.sensitivity);
    // Process shapes in parallel: each shape's seed extraction is independent
    // (no shared state). C++ also processes shapes concurrently (interleaved
    // with chunks via the index_chunks loop).
    let shape_results: Vec<Vec<seed_match::SeedMatch>> = shapes
        .par_iter()
        .map(|shape| {
            let complexity_cut = traits.seed_cut * std::f64::consts::LN_2 * shape.weight as f64;
            parallel::find_seed_matches_partitioned_filtered_min_query_len(
                &query_seqs,
                &db_seqs,
                shape,
                &reduction,
                complexity_cut,
                traits.freq_sd,
                config.min_query_len,
            )
        })
        .collect();
    let mut all_seed_matches = Vec::new();
    for matches in shape_results {
        all_seed_matches.extend(matches);
    }
    let raw_seed_matches = all_seed_matches.len();

    // Dedupe (q_id, q_pos, r_id, r_pos) across shapes. C++ doesn't dedupe at
    // this stage, but the downstream pipeline largely collapses duplicates
    // anyway. Keeping this dedup avoids a large runtime hit in the native
    // benchmark path and does not affect the current residual divergent rows.
    all_seed_matches.sort_by_key(|m| (m.query_id, m.query_pos, m.ref_id, m.ref_pos));
    all_seed_matches.dedup_by_key(|m| (m.query_id, m.query_pos, m.ref_id, m.ref_pos));
    let after_dedup = all_seed_matches.len();

    // Stage-1 hamming filter: drop seed hits whose 48-letter window matches
    // fewer than `min_identities` letters. Ports C++ search/hamming/kernel.h.
    all_seed_matches = crate::search::hamming_filter::apply_hamming_filter(
        all_seed_matches,
        &query_seqs,
        &db_seqs,
        traits.min_identities,
    );
    eprintln!(
        "Seed matches: {} (raw) -> {} (dedup) -> {} (hamming)",
        raw_seed_matches,
        after_dedup,
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

    // Group all seed matches by query
    let mut query_matches: HashMap<u32, Vec<&seed_match::SeedMatch>> = HashMap::new();
    for m in &all_seed_matches {
        query_matches.entry(m.query_id).or_default().push(m);
    }

    // Process each query — parallel over queries (each query is independent).
    // Each worker formats its output into a per-query buffer; we write the
    // buffers to the final writer in input order to keep output deterministic.
    let per_query_output: Vec<Vec<u8>> = query_records
        .par_iter()
        .enumerate()
        .map(|(query_idx, query_rec)| -> io::Result<Vec<u8>> {
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

            // Get seed matches for this query
            let empty_matches = Vec::new();
            let query_seed_matches = query_matches
                .get(&(query_idx as u32))
                .unwrap_or(&empty_matches);

            // Group by target
            let mut target_hits: HashMap<u32, Vec<&&seed_match::SeedMatch>> = HashMap::new();
            for m in query_seed_matches {
                target_hits.entry(m.ref_id).or_default().push(m);
            }

            // Process each target
            let mut target_results: Vec<(u32, f64, OutputHsp)> = Vec::new();

            for (&target_id, hits) in &target_hits {
                let target = &db_records[target_id as usize].sequence;

                // Stage-2 filter: 96-letter ungapped window walk (C++
                // `dp/ungapped_align.cpp:ungapped_window`). This is the score C++
                // uses for the seed-hit cutoff check — deliberately less sensitive
                // than x-drop so marginal targets get rejected here rather than
                // sailing through to full alignment.
                let mut any_seed_passes_window = false;
                for hit in hits {
                    let s = crate::dp::ungapped_window::filter_ungapped_score(
                        query,
                        target,
                        hit.query_pos as usize,
                        hit.ref_pos as usize,
                        crate::dp::ungapped_window::UNGAPPED_WINDOW,
                        &score_matrix,
                    );
                    if s > ungapped_cutoff {
                        any_seed_passes_window = true;
                        break;
                    }
                }
                if !any_seed_passes_window {
                    continue;
                }

                let mut seed_hits = Vec::with_capacity(hits.len());
                for hit in hits {
                    let score = crate::dp::ungapped_window::filter_ungapped_score(
                        query,
                        target,
                        hit.query_pos as usize,
                        hit.ref_pos as usize,
                        crate::dp::ungapped_window::UNGAPPED_WINDOW,
                        &score_matrix,
                    );
                    seed_hits.push(SeedHit::new(
                        hit.query_pos as i32,
                        hit.ref_pos as i32,
                        score,
                        0,
                    ));
                }

                let query_comp = crate::stats::cbs::compute_composition(query);
                let mut stat = Statistics::new();
                // C++ `raw_ungapped_xdrop` is computed at startup as
                // `score_matrix.rawscore(config.ungapped_xdrop)` where
                // `ungapped_xdrop` is expressed in bits (`--xdrop`, default
                // 12.3). Compute it from the active matrix so non-BLOSUM62
                // runs and user overrides match C++.
                let ungapped_cfg = UngappedStageConfig {
                    comp_based_stats: config.comp_based_stats,
                    xdrop: score_matrix.rawscore_int(config.ungapped_xdrop_bits),
                    ..UngappedStageConfig::default()
                };
                // C++ uses `BandedSlow` for `MoreSensitive`/`VerySensitive`/
                // `UltraSensitive` and `BandedFast` for the others (per
                // `Extension::default_ext_mode` in `extend.cpp:73-86`).
                // Hardcoding `BandedFast` here over-truncates the band at high
                // sensitivities, dropping marginal alignments C++ would find.
                let ext_mode = sensitivity::default_ext_mode(config.sensitivity);
                let work_target = ungapped_stage_target(
                    &mut seed_hits,
                    std::slice::from_ref(query),
                    std::slice::from_ref(&query_cbs),
                    &query_comp,
                    target_id,
                    target.len() as i32,
                    &mut stat,
                    &db_block,
                    ext_mode,
                    &ungapped_cfg,
                    &score_matrix,
                );
                if work_target.hsp[0].is_empty() {
                    continue;
                }

                let gapped_cfg = GappedScoreConfig {
                    comp_based_stats_hauser: config.comp_based_stats.hauser(),
                    comp_based_stats_matrix_adjust: config.comp_based_stats.matrix_adjust(),
                    query_cover: config.query_cover,
                    subject_cover: config.subject_cover,
                    no_self_hits: config.no_self_hits,
                    max_evalue: config.max_evalue,
                    min_id: config.min_id,
                    max_target_seqs: config.max_target_seqs,
                    sensitivity: config.sensitivity,
                    ..GappedScoreConfig::default()
                };
                let mut work_targets = vec![work_target];
                let aligned = align_work_targets(
                    &mut work_targets,
                    std::slice::from_ref(query),
                    Some(&query_rec.id),
                    std::slice::from_ref(&query_cbs),
                    query.len() as i32,
                    Flags::NONE,
                    HspValues::COORDS
                        | HspValues::IDENT
                        | HspValues::LENGTH
                        | HspValues::MISMATCHES
                        | HspValues::GAP_OPENINGS,
                    ext_mode,
                    &mut stat,
                    &gapped_cfg,
                    &score_matrix,
                );
                let Some(best_target) = aligned.first() else {
                    continue;
                };
                let Some(best_hsp) = best_target.hsp[0].first() else {
                    continue;
                };

                let evalue = best_hsp.evalue;
                if evalue > config.max_evalue {
                    continue;
                }

                let pident = if best_hsp.length > 0 {
                    100.0 * best_hsp.identities as f64 / best_hsp.length as f64
                } else {
                    0.0
                };
                if pident < config.min_id {
                    continue;
                }
                let qcov = if query.is_empty() {
                    0.0
                } else {
                    100.0 * (best_hsp.query_range.end - best_hsp.query_range.begin).max(0) as f64
                        / query.len() as f64
                };
                if qcov < config.query_cover {
                    continue;
                }
                let scov = if target.is_empty() {
                    0.0
                } else {
                    100.0
                        * (best_hsp.subject_range.end - best_hsp.subject_range.begin).max(0) as f64
                        / target.len() as f64
                };
                if scov < config.subject_cover {
                    continue;
                }
                if config.no_self_hits
                    && query == target
                    && query_rec.id == db_records[target_id as usize].id
                {
                    continue;
                }

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

            // Final top-N sort. C++ `Match::cmp_evalue` (`diamond/src/align/extend.h:58-63`)
            // orders by `(evalue ASC, score DESC, target_block_id ASC)`. The score
            // tie-break is critical at the `-k 25` boundary: when two targets share
            // the same evalue, C++ keeps the higher-scoring one. Without it Rust
            // falls straight to block_id and can swap a higher-score target out
            // for a lower-score one — produces the residual cpp-only / rust-only
            // flips around the cutoff.
            target_results.sort_by(|a, b| {
                a.1.partial_cmp(&b.1)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| b.2.score.cmp(&a.2.score))
                    .then(a.0.cmp(&b.0))
            });

            // Apply max target seqs
            let max_targets = config.max_target_seqs as usize;
            if target_results.len() > max_targets {
                target_results.truncate(max_targets);
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
        })
        .collect::<io::Result<Vec<_>>>()?;

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
