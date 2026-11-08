//! Concrete adapters for Source's shared external linclust job pipeline.
use std::cell::Cell;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::basic::reduction::Reduction;
use crate::basic::shape::Shape;
use crate::basic::value::{Letter, SEED_MASK};
use crate::cluster::cascaded::wrapper::CascadedWrapperConfig;
use crate::cluster::external::{
    self, ExternalConfig, ExternalPipeline, ExternalRecordKind, Job, PairFileArray,
    PairTableConfig, RadixedTable, SeedMasker, SeedSequence, SeedTableConfig, SequenceVolumeReader,
    VolumedFile,
};
use crate::config::Sensitivity;
use crate::masking::{MaskingAlgo, MaskingStat};
use crate::search::sensitivity;
use crate::stats::score_matrix::ScoreMatrix;

struct Reader {
    errors: crate::util::io::IoErrorOwner,
}
impl SequenceVolumeReader for Reader {
    fn visit_volume_bytes_while_with_cleanup(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        on_error: &mut dyn FnMut(),
        visitor: &mut dyn FnMut(external::seed_table::ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        external::reader::visit_volume_bytes_while_owned_with_cleanup(
            path,
            should_read,
            visitor,
            on_open,
            on_error,
            &self.errors,
        )
    }
    fn visit_volume_bytes_while_with_open(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(external::seed_table::ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        external::reader::visit_volume_bytes_while_owned(
            path,
            should_read,
            visitor,
            on_open,
            &self.errors,
        )
    }
    fn fork_reader(&self) -> Option<Box<dyn SequenceVolumeReader + Send>> {
        Some(Box::new(Self {
            errors: self.errors.clone(),
        }))
    }
    fn visit_volume_bytes(
        &mut self,
        path: &Path,
        visitor: &mut dyn FnMut(external::seed_table::ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        external::reader::visit_volume_bytes_owned(path, visitor, &mut || Ok(()), &self.errors)
    }
    fn visit_volume_bytes_while(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        visitor: &mut dyn FnMut(external::seed_table::ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        external::reader::visit_volume_bytes_while_owned(
            path,
            should_read,
            visitor,
            &mut || Ok(()),
            &self.errors,
        )
    }
    fn visit_volume_bytes_with_open(
        &mut self,
        path: &Path,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(external::seed_table::ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        external::reader::visit_volume_bytes_owned(path, visitor, on_open, &self.errors)
    }
    fn read_volume(&mut self, path: &Path) -> Result<Vec<SeedSequence>, String> {
        let mut records = Vec::new();
        self.visit_volume(path, &mut |record| {
            records.push(record);
            Ok(())
        })?;
        Ok(records)
    }
    fn visit_volume(
        &mut self,
        path: &Path,
        visit: &mut dyn FnMut(SeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        external::reader::visit_volume_owned(path, visit, &mut || Ok(()), &self.errors)
    }
    fn visit_volume_with_open(
        &mut self,
        path: &Path,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(SeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        external::reader::visit_volume_owned(path, visitor, on_open, &self.errors)
    }
}

#[derive(Clone)]
struct Masker(Arc<crate::masking::tantan::TantanMasker>, i32);
impl SeedMasker for Masker {
    fn fork_masker(&self) -> Option<Box<dyn SeedMasker + Send>> {
        Some(Box::new(self.clone()))
    }
    fn mask(&mut self, sequence: &mut [Letter], tantan: bool, motif: bool) -> MaskingStat {
        let mut stats = MaskingStat::new();
        if tantan {
            self.0.mask_bit(sequence);
            let count = sequence
                .iter()
                .filter(|&&letter| letter & SEED_MASK != 0)
                .count();
            crate::masking::bit_to_hard_mask(sequence);
            stats.add(MaskingAlgo::Tantan, count as u64);
        }
        if motif {
            let masked =
                crate::masking::motifs::mask_motifs_with_limit(sequence, self.1.max(0) as usize);
            stats.add(MaskingAlgo::Motif, masked.len() as u64);
        }
        stats
    }
}

struct Pipeline {
    reader: Reader,
    masker: Masker,
    matrix: Arc<ScoreMatrix>,
    sensitivity: Cell<Sensitivity>,
    threads: usize,
    member_cover: f64,
    mutual_cover: Option<f64>,
    approx_min_id: f64,
    max_evalue: f64,
    soft_masking: Vec<u8>,
    motif_masking: Vec<u8>,
    alignment_output: Option<PathBuf>,
    pair_files: Option<PairFileArray>,
    tuning: super::blastp::SearchTuning,
}
impl SequenceVolumeReader for Pipeline {
    fn visit_volume_bytes_while_with_cleanup(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        on_error: &mut dyn FnMut(),
        visitor: &mut dyn FnMut(external::seed_table::ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        self.reader.visit_volume_bytes_while_with_cleanup(
            path,
            should_read,
            on_open,
            on_error,
            visitor,
        )
    }
    fn visit_volume_bytes_while_with_open(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(external::seed_table::ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        self.reader
            .visit_volume_bytes_while_with_open(path, should_read, on_open, visitor)
    }
    fn fork_reader(&self) -> Option<Box<dyn SequenceVolumeReader + Send>> {
        self.reader.fork_reader()
    }
    fn visit_volume_bytes(
        &mut self,
        path: &Path,
        visitor: &mut dyn FnMut(external::seed_table::ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        self.reader.visit_volume_bytes(path, visitor)
    }
    fn visit_volume_bytes_while(
        &mut self,
        path: &Path,
        should_read: &mut dyn FnMut() -> bool,
        visitor: &mut dyn FnMut(external::seed_table::ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        self.reader
            .visit_volume_bytes_while(path, should_read, visitor)
    }
    fn visit_volume_bytes_with_open(
        &mut self,
        path: &Path,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(external::seed_table::ByteSeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        self.reader
            .visit_volume_bytes_with_open(path, on_open, visitor)
    }
    fn read_volume(&mut self, path: &Path) -> Result<Vec<SeedSequence>, String> {
        self.reader.read_volume(path)
    }
    fn visit_volume(
        &mut self,
        path: &Path,
        visit: &mut dyn FnMut(SeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        self.reader.visit_volume(path, visit)
    }
    fn visit_volume_with_open(
        &mut self,
        path: &Path,
        on_open: &mut dyn FnMut() -> Result<(), String>,
        visitor: &mut dyn FnMut(SeedSequence) -> Result<(), String>,
    ) -> Result<(), String> {
        self.reader.visit_volume_with_open(path, on_open, visitor)
    }
}
impl ExternalPipeline for Pipeline {
    fn io_errors(&self) -> crate::util::io::IoErrorOwner {
        self.reader.errors.clone()
    }
    fn shape_count(&self, sensitivity: Sensitivity) -> usize {
        self.sensitivity.set(sensitivity);
        sensitivity::get_shape_codes(sensitivity).len()
    }
    fn build_seed_table(
        &mut self,
        job: &Job,
        volumes: &VolumedFile,
        shape: usize,
    ) -> Result<RadixedTable, String> {
        let sensitivity = self.sensitivity.get();
        let reduction = Reduction::default_reduction();
        let shapes: Vec<_> = sensitivity::get_shape_codes(sensitivity)
            .iter()
            .map(|code| Shape::from_code(code, &reduction))
            .collect();
        let traits = sensitivity::get_traits(sensitivity);
        let config = SeedTableConfig {
            threads: self.threads,
            soft_masking: &self.soft_masking,
            motif_masking: &self.motif_masking,
            sketch_size: self.tuning.setup.sketch_size.unwrap_or(0),
            sensitivity_sketch_size: traits.sketch_size,
            seed_cut: self.tuning.setup.seed_cut,
            sensitivity_seed_cut: traits.seed_cut,
            shapes: &shapes,
            reduction: &reduction,
        };
        external::seed_table::build_seed_table(
            job,
            volumes,
            shape,
            &config,
            &mut self.reader,
            &mut self.masker,
        )
        .map(|(table, _)| table)
    }
    fn radix_sort(
        &mut self,
        job: &Job,
        table: RadixedTable,
        kind: ExternalRecordKind,
        bits_unsorted: u32,
    ) -> Result<RadixedTable, String> {
        external::storage::radix_sort(job, table, kind, bits_unsorted)
    }
    fn begin_pair_table(&mut self, job: &Job) -> Result<(), String> {
        self.pair_files = Some(PairFileArray::new_bounded_owned(
            job.base_dir(None).join("pair_table"),
            job.worker_id(),
            job.io_errors(),
        )?);
        Ok(())
    }
    fn build_pair_table(
        &mut self,
        job: &Job,
        seeds: &RadixedTable,
        shape: usize,
        max_oid: u64,
    ) -> Result<RadixedTable, String> {
        let min_length_ratio = self.mutual_cover.map_or(0.0, |coverage| {
            if self.sensitivity.get() < Sensitivity::Linclust40 {
                (coverage / 100.0 + 0.05).min(1.0)
            } else {
                coverage / 100.0 - 0.05
            }
        });
        external::pair_table::build_pair_table(
            job,
            seeds,
            shape as i32,
            i64::try_from(max_oid).map_err(|_| "External OID exceeds signed64 range".to_owned())?,
            self.pair_files
                .as_mut()
                .ok_or("Pair table owner is absent")?,
            PairTableConfig {
                threads: self.threads,
                mutual_cover_present: self.mutual_cover.is_some(),
                min_length_ratio,
            },
        )
    }
    fn finish_pair_table(&mut self) -> Result<(), String> {
        self.pair_files
            .take()
            .ok_or("Pair table owner is absent")?
            .finish()?;
        Ok(())
    }
    fn align(
        &mut self,
        job: &Job,
        chunk_count: usize,
        max_oid: u64,
    ) -> Result<RadixedTable, String> {
        Arc::get_mut(&mut self.matrix)
            .expect(
                "External alignment workers must release the scoring matrix before the next round",
            )
            .set_db_letters(1_000_000_000);
        external::align::align_owned(
            job,
            chunk_count,
            i64::try_from(max_oid).map_err(|_| "External OID exceeds signed64 range".to_owned())?,
            Arc::clone(&self.matrix),
            &external::AlignConfig {
                threads: self.threads,
                approx_min_identity: self.approx_min_id,
                member_cover: self.member_cover,
                mutual_cover: self.mutual_cover,
                alignment_output: self.alignment_output.clone(),
                max_evalue: self.max_evalue,
                cutoff_score_8bit: self.tuning.cutoff_score_8bit,
                max_swipe_dp: self.tuning.max_swipe_dp,
                swipe_task_size: self.tuning.swipe_task_size,
                band_bin: self.tuning.band_bin,
                col_bin: self.tuning.col_bin,
                approx_backtrace: self.tuning.approx_backtrace,
                ..Default::default()
            },
            &mut self.reader,
        )
        .map(|(table, _)| table)
    }
}

pub(crate) fn run(
    config: &CascadedWrapperConfig,
    output: &Path,
    chunk_size: u64,
    alignment_output: Option<&Path>,
    matrix: ScoreMatrix,
    tuning: &super::blastp::SearchTuning,
    errors: &crate::util::io::IoErrorOwner,
) -> Result<(), String> {
    if output.as_os_str().is_empty() {
        return Err("Option missing: output file (--out/-o)".to_owned());
    }
    let database = config
        .database
        .as_deref()
        .ok_or("Option missing: database file (--db/-d)")?;
    let volumes = external::read_volumed_file_manifest_owned(database, errors)?;
    // Actual Source FileArray/BufferArray cleanup is noexcept. Keep this
    // invocation policy out of the reusable library's fallible Job APIs.
    let errors = errors.with_source_destructors();
    let mut setup = ExternalConfig {
        output_file: output.to_path_buf(),
        oid_output: config.oid_output,
        // Source External::Job places its shared queues under this versioned
        // directory; the library Job API already accepts the final root.
        root_dir: config.parallel_tmpdir.join(format!(
            "diamond-tmp-{}",
            crate::basic::consts::VERSION_STRING
        )),
        memory_limit: config.memory_limit,
        threads: config.threads.max(1) as usize,
        mutual_cover: config.thresholds.mutual_cover,
        member_cover: config.thresholds.member_cover.unwrap_or(80.0),
        approx_min_id: config.core.approx_min_id,
        linclust_chunk_size: chunk_size,
        cluster_steps: config.core.helpers.cluster_steps.clone(),
        ..Default::default()
    };
    let mut backend = Pipeline {
        reader: Reader {
            errors: errors.clone(),
        },
        masker: Masker(
            Arc::new(crate::masking::tantan::TantanMasker::from_score_matrix(
                &matrix,
                tuning.tantan_min_mask_prob as f32,
            )),
            tuning.max_motif_len,
        ),
        matrix: Arc::new(matrix),
        sensitivity: Cell::new(Sensitivity::Default),
        threads: setup.threads,
        member_cover: setup.member_cover,
        mutual_cover: setup.mutual_cover,
        approx_min_id: setup.approx_min_id,
        max_evalue: config.core.max_evalue,
        soft_masking: tuning.setup.soft_masking.clone(),
        motif_masking: config.core.motif_masking.clone(),
        alignment_output: alignment_output.map(Path::to_path_buf),
        pair_files: None,
        tuning: tuning.clone(),
    };
    external::external(&mut setup, &volumes, &mut backend)?;
    Ok(())
}
