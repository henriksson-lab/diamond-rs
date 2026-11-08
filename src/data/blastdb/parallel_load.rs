//! Concrete `SequenceFile::load_parallel` for BLAST raw chunks.
use super::{
    blastdb::BlastDB,
    volume::{DbFilter, DecodedPackage, RawChunk},
};

pub(crate) fn load(
    file: &mut BlastDB,
    max_letters: u64,
    load_threads: i32,
    threads: i32,
    minichunk: usize,
    filter: Option<&DbFilter>,
) -> Result<(DecodedPackage, usize), String> {
    let t = load_threads.min(threads);
    let workers = if t >= 3 { (t - 2) as usize } else { 1 };
    let flags = file.flags();
    let mut loaded_letters = 0u64;
    let mut processed = 0usize;
    let mut combined = DecodedPackage::default();
    crate::util::parallel::ordered_pipeline::run(
        workers,
        || {
            if loaded_letters >= max_letters {
                return Ok(None);
            }
            let letters = (max_letters - loaded_letters).min(minichunk as u64) as usize;
            let raw = file.raw_chunk(letters, flags)?;
            if raw.empty() {
                return Ok(None);
            }
            loaded_letters = loaded_letters.saturating_add(raw.letters as u64);
            processed = processed
                .checked_add((raw.end - raw.begin) as usize)
                .ok_or("Reference sequence count overflow.")?;
            Ok(Some(raw))
        },
        |raw: RawChunk| raw.decode(flags, filter, None),
        |package| {
            combined.ids.append(&package.ids);
            combined.seqs.append(&package.seqs);
            combined.oids.extend(package.oids);
            combined.taxids.extend(package.taxids);
            Ok(())
        },
    )?;
    Ok((combined, processed))
}
