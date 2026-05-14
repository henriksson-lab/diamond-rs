use clap::{Args, Parser, Subcommand};

use crate::util::io::Compressor;

/// DIAMOND protein aligner — Rust port
#[derive(Parser, Debug)]
#[command(name = "diamond", version, about = "DIAMOND protein aligner")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Build DIAMOND database from FASTA file
    Makedb(MakedbArgs),
    /// Align protein sequences against a protein database
    Blastp(AlignArgs),
    /// Align DNA sequences against a protein database
    Blastx(AlignArgs),
    /// View DIAMOND alignment archive (DAA) file
    View(ViewArgs),
    /// Print database information
    Dbinfo(DbinfoArgs),
    /// Cluster protein sequences
    Cluster(ClusterArgs),
    /// Linear-time clustering
    Linclust(ClusterArgs),
    /// Deep clustering
    Deepclust(ClusterArgs),
    /// Retrieve sequences from database
    Getseq(GetseqArgs),
    /// Display version
    Version,
    /// Run built-in regression tests
    Test(TestArgs),
}

/// Arguments for building a database.
#[derive(Args, Debug)]
pub struct MakedbArgs {
    /// Input reference file(s) in FASTA format
    #[arg(long = "in", required = true)]
    pub input: Vec<String>,

    /// Database output file
    #[arg(short = 'd', long = "db", required = true)]
    pub database: String,

    /// Number of threads
    #[arg(short = 'p', long = "threads", default_value = "0")]
    pub threads: i32,
}

/// Common alignment arguments (blastp/blastx).
#[derive(Args, Debug)]
pub struct AlignArgs {
    /// Query input file(s)
    #[arg(short = 'q', long = "query")]
    pub query: Vec<String>,

    /// Database file
    #[arg(short = 'd', long = "db")]
    pub database: Option<String>,

    /// Output file
    #[arg(short = 'o', long = "out")]
    pub output: Option<String>,

    /// Number of threads
    #[arg(short = 'p', long = "threads", default_value = "0")]
    pub threads: i32,

    // --- Sensitivity ---
    /// Enable faster mode
    #[arg(long)]
    pub faster: bool,

    /// Enable fast mode
    #[arg(long)]
    pub fast: bool,

    /// Enable mid-sensitive mode
    #[arg(long = "mid-sensitive")]
    pub mid_sensitive: bool,

    /// Enable sensitive mode
    #[arg(long)]
    pub sensitive: bool,

    /// Enable more-sensitive mode
    #[arg(long = "more-sensitive")]
    pub more_sensitive: bool,

    /// Enable very-sensitive mode
    #[arg(long = "very-sensitive")]
    pub very_sensitive: bool,

    /// Enable ultra-sensitive mode
    #[arg(long = "ultra-sensitive")]
    pub ultra_sensitive: bool,

    // --- Scoring ---
    /// Scoring matrix
    #[arg(long = "matrix", default_value = "blosum62")]
    pub matrix: String,

    /// Gap open penalty
    #[arg(long = "gapopen", default_value = "-1")]
    pub gap_open: i32,

    /// Gap extension penalty
    #[arg(long = "gapextend", default_value = "-1")]
    pub gap_extend: i32,

    /// Frame shift penalty (blastx)
    #[arg(short = 'F', long = "frameshift", default_value = "0")]
    pub frame_shift: i32,

    // --- Filtering ---
    /// Maximum e-value to report
    #[arg(short = 'e', long = "evalue", default_value = "0.001")]
    pub max_evalue: f64,

    /// Minimum identity percentage
    #[arg(long = "id", default_value = "0")]
    pub min_id: f64,

    /// Minimum query cover percentage
    #[arg(long = "query-cover", default_value = "0")]
    pub query_cover: f64,

    /// Minimum subject cover percentage
    #[arg(long = "subject-cover", default_value = "0")]
    pub subject_cover: f64,

    /// Maximum number of target sequences
    #[arg(short = 'k', long = "max-target-seqs", default_value = "25")]
    pub max_target_seqs: i64,

    /// Report alignments within this percentage of top score
    #[arg(long = "top")]
    pub top: Option<f64>,

    /// Minimum bit score
    #[arg(long = "min-score")]
    pub min_score: Option<f64>,

    /// Maximum HSPs per target
    #[arg(long = "max-hsps", default_value = "1")]
    pub max_hsps: u32,

    // --- Output ---
    /// Output format (0/5/6/100/101/102/103/104 + field names)
    #[arg(short = 'f', long = "outfmt")]
    pub outfmt: Vec<String>,

    /// Compression (0=none, 1=gzip)
    #[arg(long = "compress", default_value = "0")]
    pub compress: String,

    /// Print header lines
    #[arg(long = "header")]
    pub header: Option<u32>,

    // --- Performance ---
    /// Sequence block size in billions of letters
    #[arg(short = 'b', long = "block-size", default_value = "2.0")]
    pub block_size: f64,

    /// Number of index chunks
    #[arg(short = 'c', long = "index-chunks", default_value = "4")]
    pub index_chunks: u32,

    /// Composition-based statistics mode
    #[arg(long = "comp-based-stats", default_value = "1")]
    pub comp_based_stats: String,

    /// Masking algorithm
    #[arg(long = "masking", default_value = "tantan")]
    pub masking: String,

    /// Motif masking
    #[arg(long = "motif-masking")]
    pub motif_masking: Option<u32>,

    /// Soft masking
    #[arg(long = "soft-masking")]
    pub soft_masking: Option<u32>,

    /// Algorithm (0=double-indexed, 1=query-indexed)
    #[arg(long = "algo")]
    pub algo: Option<u32>,

    /// Enable swipe mode
    #[arg(long)]
    pub swipe: bool,

    /// Global ranking
    #[arg(long = "global-ranking")]
    pub global_ranking: Option<u32>,

    /// Approximate identity threshold
    #[arg(long = "approx-id")]
    pub approx_id: Option<f64>,

    /// File buffer size
    #[arg(long = "file-buffer-size")]
    pub file_buffer_size: Option<u64>,

    /// Temporary directory
    #[arg(short = 't', long = "tmpdir")]
    pub tmpdir: Option<String>,

    /// Query parallel limit
    #[arg(long = "query-parallel-limit")]
    pub query_parallel_limit: Option<u32>,

    /// Ungapped x-drop
    #[arg(short = 'x', long = "xdrop")]
    pub ungapped_xdrop: Option<f64>,

    /// DNA strand (both/plus/minus)
    #[arg(long = "strand")]
    pub strand: Option<String>,

    /// Query genetic code
    #[arg(long = "query-gencode", default_value = "1")]
    pub query_gencode: u32,

    /// Minimum ORF length
    #[arg(long = "min-orf")]
    pub min_orf: Option<u32>,
}

/// Arguments for viewing DAA files.
#[derive(Args, Debug)]
pub struct ViewArgs {
    /// DAA input file
    #[arg(short = 'a', long = "daa")]
    pub daa: String,

    /// Output file
    #[arg(short = 'o', long = "out")]
    pub output: Option<String>,

    /// Output format
    #[arg(short = 'f', long = "outfmt")]
    pub outfmt: Vec<String>,
}

/// Arguments for database info.
#[derive(Args, Debug)]
pub struct DbinfoArgs {
    /// Database file
    #[arg(short = 'd', long = "db", required = true)]
    pub database: String,
}

/// Arguments for clustering.
#[derive(Args, Debug)]
pub struct ClusterArgs {
    /// Database file
    #[arg(short = 'd', long = "db", required = true)]
    pub database: String,

    /// Output file
    #[arg(short = 'o', long = "out", required = true)]
    pub output: String,

    /// Number of threads
    #[arg(short = 'p', long = "threads", default_value = "0")]
    pub threads: i32,

    /// Minimum member coverage
    #[arg(long = "member-cover", default_value = "80")]
    pub member_cover: f64,

    /// Approximate identity threshold
    #[arg(long = "approx-id", default_value = "0")]
    pub approx_id: f64,
}

/// Arguments for sequence retrieval.
#[derive(Args, Debug)]
pub struct GetseqArgs {
    /// Database file
    #[arg(short = 'd', long = "db", required = true)]
    pub database: String,

    /// Sequence IDs (comma-separated)
    #[arg(long = "seq")]
    pub seq: Option<String>,
}

/// Arguments for testing.
#[derive(Args, Debug)]
pub struct TestArgs {
    /// Enable bootstrap mode (print hashes)
    #[arg(long)]
    pub bootstrap: bool,

    /// Enable debug logging
    #[arg(long = "log")]
    pub log: bool,
}

/// Sensitivity level.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Sensitivity {
    Faster = 0,
    Fast = 1,
    Default = 2,
    Linclust40 = 3,
    Linclust20 = 4,
    Shapes6x10 = 5,
    Shapes30x10 = 6,
    MidSensitive = 7,
    Sensitive = 8,
    MoreSensitive = 9,
    VerySensitive = 10,
    UltraSensitive = 11,
}

impl AlignArgs {
    /// Determine the sensitivity level from the flags.
    pub fn sensitivity(&self) -> Sensitivity {
        if self.ultra_sensitive {
            Sensitivity::UltraSensitive
        } else if self.very_sensitive {
            Sensitivity::VerySensitive
        } else if self.more_sensitive {
            Sensitivity::MoreSensitive
        } else if self.sensitive {
            Sensitivity::Sensitive
        } else if self.mid_sensitive {
            Sensitivity::MidSensitive
        } else if self.fast {
            Sensitivity::Fast
        } else if self.faster {
            Sensitivity::Faster
        } else {
            Sensitivity::Default
        }
    }
}

pub fn block_size(
    memory_limit: i64,
    db_letters: i64,
    sensitivity: Sensitivity,
    lin: bool,
    thread_count: i32,
) -> (f64, i32) {
    const AVG_SEQ_LENGTH_EST: f64 = 200.0;
    let reduction = crate::basic::reduction::Reduction::default_reduction();
    let traits = crate::search::sensitivity::get_traits(sensitivity);
    let m = memory_limit as f64 / 1e9;
    let min = traits.minimizer_window;
    let sketch_size = traits.sketch_size;
    let max_c = if min > 0 || sketch_size > 0 {
        1
    } else if lin {
        16
    } else {
        4
    };
    let shape_weight = crate::basic::shape::Shape::from_code(
        crate::search::sensitivity::get_shape_codes(sensitivity)[0],
        &reduction,
    )
    .weight;
    let max_b = if lin {
        32768.0
    } else if sensitivity <= Sensitivity::Default {
        12.0
    } else if sensitivity <= Sensitivity::MoreSensitive {
        6.0
    } else {
        1.6
    };

    assert!(min == 0 || sketch_size == 0);
    let mut c = 0;
    let mut b;
    loop {
        c += 1;
        let mut seeds_per_letter = if sketch_size > 0 {
            sketch_size as f64 / AVG_SEQ_LENGTH_EST
        } else {
            1.0
        } / c as f64;
        if min > 0 {
            seeds_per_letter /= min as f64 / 2.0;
        }
        let seedp_bits =
            crate::search::sensitivity::seedp_bits(shape_weight, thread_count, c, &reduction);
        let seedp_per_chunk = crate::basic::seed::seedp_count(seedp_bits) / c as u64;
        let hash_join_factor = 1.0 + thread_count as f64 / seedp_per_chunk as f64;
        let seed_array_entry_size = 18.0 * hash_join_factor;
        b = m / (seed_array_entry_size * seeds_per_letter + 2.0);
        if !((b as i64) * 1_000_000_000 < db_letters && b < max_b && c < max_c) {
            break;
        }
    }
    b = b.min(max_b);
    (b.max(0.001), c)
}

pub fn set_string_option<T>(s: &str, name: &str, values: &[(&str, T)]) -> Result<T, String>
where
    T: Copy + Default,
{
    if s.is_empty() {
        return Ok(T::default());
    }
    for (key, value) in values {
        if s == *key {
            return Ok(*value);
        }
    }
    let allowed = values.iter().fold(String::new(), |mut acc, (key, _)| {
        acc.push(' ');
        acc.push_str(key);
        acc
    });
    Err(format!(
        "Invalid argument for option {name}. Allowed values are:{allowed}"
    ))
}

pub fn compressor(compression: &str) -> Result<Compressor, String> {
    if compression.is_empty() || compression == "0" {
        Ok(Compressor::None)
    } else if compression == "1" {
        Ok(Compressor::Zlib)
    } else if compression == "zstd" {
        Ok(Compressor::Zstd)
    } else {
        Err(format!("Invalid compression algorithm: {compression}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cli_parse_version() {
        let cli = Cli::try_parse_from(["diamond", "version"]);
        assert!(cli.is_ok());
    }

    #[test]
    fn test_cli_parse_blastp() {
        let cli = Cli::try_parse_from([
            "diamond",
            "blastp",
            "-q",
            "query.faa",
            "-d",
            "db",
            "-o",
            "out.txt",
        ]);
        assert!(cli.is_ok());
        if let Command::Blastp(args) = cli.unwrap().command {
            assert_eq!(args.query, vec!["query.faa"]);
            assert_eq!(args.database, Some("db".to_string()));
            assert_eq!(args.output, Some("out.txt".to_string()));
            assert_eq!(args.sensitivity(), Sensitivity::Default);
        }
    }

    #[test]
    fn test_cli_parse_sensitivity() {
        let cli = Cli::try_parse_from([
            "diamond",
            "blastp",
            "-q",
            "q.faa",
            "-d",
            "db",
            "--more-sensitive",
        ]);
        if let Command::Blastp(args) = cli.unwrap().command {
            assert_eq!(args.sensitivity(), Sensitivity::MoreSensitive);
        }
    }

    #[test]
    fn test_cli_parse_makedb() {
        let cli = Cli::try_parse_from(["diamond", "makedb", "--in", "ref.faa", "-d", "db"]);
        assert!(cli.is_ok());
    }

    #[test]
    fn test_block_size_calculation() {
        assert_eq!(block_size(0, 0, Sensitivity::Default, false, 1), (0.001, 1));
        let (b, c) = block_size(8_000_000_000, 5_000_000_000, Sensitivity::Default, false, 4);
        assert!(b > 0.0);
        assert!((1..=4).contains(&c));
        let (_, lin_c) = block_size(8_000_000_000, 5_000_000_000, Sensitivity::Default, true, 4);
        assert!((1..=16).contains(&lin_c));
    }

    #[test]
    fn test_set_string_option() {
        let values = [("one", 1_i32), ("two", 2_i32)];
        assert_eq!(set_string_option("", "mode", &values).unwrap(), 0);
        assert_eq!(set_string_option("two", "mode", &values).unwrap(), 2);
        assert_eq!(
            set_string_option("three", "mode", &values).unwrap_err(),
            "Invalid argument for option mode. Allowed values are: one two"
        );
    }

    #[test]
    fn test_compressor() {
        assert_eq!(compressor("").unwrap(), Compressor::None);
        assert_eq!(compressor("0").unwrap(), Compressor::None);
        assert_eq!(compressor("1").unwrap(), Compressor::Zlib);
        assert_eq!(compressor("zstd").unwrap(), Compressor::Zstd);
        assert_eq!(
            compressor("gzip").unwrap_err(),
            "Invalid compression algorithm: gzip"
        );
    }
}
