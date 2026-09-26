//! Final multinode clustering merge from
//! `diamond/src/cluster/multinode/output.cpp`.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::basic::value::OId;

const NIL: OId = OId::MAX;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Volume {
    pub path: PathBuf,
    pub oid_begin: OId,
    pub oid_end: OId,
    pub record_count: OId,
}

impl Volume {
    pub fn new(path: impl Into<PathBuf>, oid_begin: OId, oid_end: OId, record_count: OId) -> Self {
        Self {
            path: path.into(),
            oid_begin,
            oid_end,
            record_count,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MergeConfig {
    pub output_file: PathBuf,
    pub oid_output: bool,
}

/// Narrow ownership adapter for the `Job` surface used by this source file.
pub trait MultinodeOutputJob {
    fn max_oid(&self) -> OId;
    fn root_dir(&self) -> &Path;
    fn base_dir(&self, round: i32) -> PathBuf;
    fn round(&self) -> i32;
    fn log(&mut self, message: &str) -> Result<(), String>;
}

/// C++ file-local `read_clusters`.
pub fn read_clusters(path: &Path, max_oid: OId) -> Result<Vec<OId>, String> {
    let bytes =
        fs::read(path).map_err(|_| format!("Error opening clustering file: {}", path.display()))?;
    let len = usize::try_from(max_oid)
        .ok()
        .and_then(|value| value.checked_add(1))
        .ok_or_else(|| "Clustering OID range does not fit in memory.".to_owned())?;
    let mut mapping = vec![NIL; len];
    let mut fields = bytes
        .split(u8::is_ascii_whitespace)
        .filter(|field| !field.is_empty());
    loop {
        let Some(rep_text) = fields.next() else {
            break;
        };
        let Some(member_text) = fields.next() else {
            break;
        };
        let Ok(rep) = std::str::from_utf8(rep_text)
            .ok()
            .unwrap_or("")
            .parse::<OId>()
        else {
            break;
        };
        let Ok(member) = std::str::from_utf8(member_text)
            .ok()
            .unwrap_or("")
            .parse::<OId>()
        else {
            break;
        };
        let member = usize::try_from(member)
            .ok()
            .filter(|member| *member < mapping.len())
            .ok_or_else(|| format!("Cluster member OID out of bounds: {member}"))?;
        mapping[member] = rep;
    }
    Ok(mapping)
}

/// C++ file-local `chain_round`.
pub fn chain_round(mapping: &mut [OId], path: &Path, max_oid: OId) -> Result<(), String> {
    let next = read_clusters(path, max_oid)?;
    for centroid in mapping {
        if *centroid != NIL {
            let index = usize::try_from(*centroid)
                .ok()
                .filter(|index| *index < next.len())
                .ok_or_else(|| format!("Cluster representative OID out of bounds: {centroid}"))?;
            *centroid = next[index];
        }
    }
    Ok(())
}

/// C++ file-local `build_merged`.
pub fn build_merged<J: MultinodeOutputJob>(job: &J) -> Result<Vec<OId>, String> {
    let mut mapping = read_clusters(&job.base_dir(0).join("clusters.tsv"), job.max_oid())?;
    for round in 1..=job.round() {
        chain_round(
            &mut mapping,
            &job.base_dir(round).join("clusters.tsv"),
            job.max_oid(),
        )?;
    }
    Ok(mapping)
}

/// C++ file-local `output_oids`.
pub fn output_oids<J: MultinodeOutputJob>(
    job: &J,
    merged: &[OId],
    output_file: &Path,
) -> Result<OId, String> {
    let file = File::create(output_file)
        .map_err(|_| format!("Error opening output file: {}", output_file.display()))?;
    let mut output = BufWriter::new(file);
    let mut cluster_count = 0;
    for oid in 0..=job.max_oid() {
        let centroid = merged[oid as usize];
        if centroid == oid {
            cluster_count += 1;
        }
        writeln!(output, "{}\t{}", centroid as i64, oid as i64)
            .map_err(|error| error.to_string())?;
    }
    output.flush().map_err(|error| error.to_string())?;
    Ok(cluster_count)
}

/// C++ file-local `output_accs`.
pub fn output_accs<J: MultinodeOutputJob>(
    job: &J,
    merged: &[OId],
    volumes: &[Volume],
    output_file: &Path,
) -> Result<OId, String> {
    let len = usize::try_from(job.max_oid())
        .ok()
        .and_then(|value| value.checked_add(1))
        .ok_or_else(|| "Accession OID range does not fit in memory.".to_owned())?;
    let mut accessions = vec![Vec::<u8>::new(); len];
    for (volume_index, volume) in volumes.iter().enumerate() {
        let path = job
            .root_dir()
            .join("accessions")
            .join(format!("{volume_index}.txt"));
        let bytes = fs::read(&path)
            .map_err(|_| format!("Error opening accessions file: {}", path.display()))?;
        let mut oid = volume.oid_begin;
        let line_count = if bytes.is_empty() {
            0
        } else {
            bytes
                .split(|byte| *byte == b'\n')
                .count()
                .saturating_sub(usize::from(bytes.last() == Some(&b'\n')))
        };
        for line in bytes.split(|byte| *byte == b'\n').take(line_count) {
            if oid >= volume.oid_end {
                break;
            }
            accessions[oid as usize] = line.to_vec();
            oid += 1;
        }
    }

    let file = File::create(output_file)
        .map_err(|_| format!("Error opening output file: {}", output_file.display()))?;
    let mut output = BufWriter::new(file);
    let mut cluster_count = 0;
    for oid in 0..=job.max_oid() {
        let centroid = merged[oid as usize];
        let centroid_index = usize::try_from(centroid)
            .ok()
            .filter(|index| *index < accessions.len())
            .ok_or_else(|| format!("Cluster representative OID out of bounds: {centroid}"))?;
        output
            .write_all(&accessions[centroid_index])
            .and_then(|_| output.write_all(b"\t"))
            .and_then(|_| output.write_all(&accessions[oid as usize]))
            .and_then(|_| output.write_all(b"\n"))
            .map_err(|error| error.to_string())?;
        if centroid == oid {
            cluster_count += 1;
        }
    }
    output.flush().map_err(|error| error.to_string())?;
    Ok(cluster_count)
}

/// C++ `merge(Job&, const VolumedFile&)` with global output settings made
/// explicit.
pub fn merge<J: MultinodeOutputJob>(
    job: &mut J,
    volumes: &[Volume],
    config: &MergeConfig,
) -> Result<OId, String> {
    job.log("Merging clusterings")?;
    let merged = build_merged(job)?;
    let cluster_count = if config.oid_output {
        output_oids(job, &merged, &config.output_file)?
    } else {
        output_accs(job, &merged, volumes, &config.output_file)?
    };
    job.log(&format!("Total clusters: {cluster_count}"))?;
    Ok(cluster_count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

    fn temp_dir(name: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "diamond_multinode_output_{name}_{}_{}",
            std::process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    struct TestJob {
        max_oid: OId,
        root: PathBuf,
        round: i32,
        logs: Vec<String>,
    }

    impl MultinodeOutputJob for TestJob {
        fn max_oid(&self) -> OId {
            self.max_oid
        }

        fn root_dir(&self) -> &Path {
            &self.root
        }

        fn base_dir(&self, round: i32) -> PathBuf {
            self.root.join(format!("round{round}"))
        }

        fn round(&self) -> i32 {
            self.round
        }

        fn log(&mut self, message: &str) -> Result<(), String> {
            self.logs.push(message.to_owned());
            Ok(())
        }
    }

    fn make_job(root: PathBuf, round: i32) -> TestJob {
        TestJob {
            max_oid: 4,
            root,
            round,
            logs: Vec::new(),
        }
    }

    #[test]
    fn read_and_chain_preserve_nil_and_follow_representatives() {
        let root = temp_dir("chain");
        let first = root.join("first.tsv");
        let second = root.join("second.tsv");
        fs::write(&first, "0 0\n0 1\n2 2\n2 3\n").unwrap();
        fs::write(&second, "4 0\n2 2\n").unwrap();
        let mut mapping = read_clusters(&first, 4).unwrap();
        chain_round(&mut mapping, &second, 4).unwrap();
        assert_eq!(mapping, [4, 4, 2, 2, NIL]);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn oid_output_is_complete_ordered_and_prints_nil_as_minus_one() {
        let root = temp_dir("oids");
        let job = make_job(root.clone(), 0);
        let path = root.join("out.tsv");
        let merged = [0, 0, 2, 2, NIL];
        assert_eq!(output_oids(&job, &merged, &path).unwrap(), 2);
        assert_eq!(
            fs::read_to_string(&path).unwrap(),
            "0\t0\n0\t1\n2\t2\n2\t3\n-1\t4\n"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn accession_output_respects_volume_ranges_and_empty_missing_lines() {
        let root = temp_dir("accs");
        fs::create_dir_all(root.join("accessions")).unwrap();
        fs::write(root.join("accessions/0.txt"), "a0\na1\nextra\n").unwrap();
        fs::write(root.join("accessions/1.txt"), "a2\na3\n").unwrap();
        let job = make_job(root.clone(), 0);
        let volumes = [
            Volume::new("unused0", 0, 2, 2),
            Volume::new("unused1", 2, 5, 3),
        ];
        let path = root.join("out.tsv");
        let merged = [0, 0, 2, 2, 2];
        assert_eq!(output_accs(&job, &merged, &volumes, &path).unwrap(), 2);
        assert_eq!(
            fs::read_to_string(&path).unwrap(),
            "a0\ta0\na0\ta1\na2\ta2\na2\ta3\na2\t\n"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn merge_chains_all_rounds_selects_format_and_logs_count() {
        let root = temp_dir("merge");
        fs::create_dir_all(root.join("round0")).unwrap();
        fs::create_dir_all(root.join("round1")).unwrap();
        fs::write(
            root.join("round0/clusters.tsv"),
            "0 0\n0 1\n2 2\n2 3\n4 4\n",
        )
        .unwrap();
        fs::write(root.join("round1/clusters.tsv"), "0 0\n0 2\n4 4\n").unwrap();
        let output = root.join("merged.tsv");
        let mut job = make_job(root.clone(), 1);
        let count = merge(
            &mut job,
            &[],
            &MergeConfig {
                output_file: output.clone(),
                oid_output: true,
            },
        )
        .unwrap();
        assert_eq!(count, 2);
        assert_eq!(
            fs::read_to_string(output).unwrap(),
            "0\t0\n0\t1\n0\t2\n0\t3\n4\t4\n"
        );
        assert_eq!(job.logs, ["Merging clusterings", "Total clusters: 2"]);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn opening_errors_match_upstream_context() {
        let root = temp_dir("errors");
        assert_eq!(
            read_clusters(&root.join("missing.tsv"), 1).unwrap_err(),
            format!(
                "Error opening clustering file: {}",
                root.join("missing.tsv").display()
            )
        );
        let job = make_job(root.clone(), 0);
        let error = output_accs(
            &job,
            &[0, 0, 0, 0, 0],
            &[Volume::new("unused", 0, 5, 5)],
            &root.join("out"),
        )
        .unwrap_err();
        assert!(error.contains("Error opening accessions file:"));
        fs::remove_dir_all(root).unwrap();
    }
}
