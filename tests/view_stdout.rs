use std::fs::File;
use std::io::Write;
use std::process::{Command, Stdio};

use diamond::basic::packed_transcript::{EditOperation, PackedOperation};
use diamond::basic::value::SequenceType;
use diamond::data::daa::{
    compute_flag, finish_daa_query_record, write_daa_query_record, BlockType, DaaFile, DaaHeader1,
    DaaHeader2, DaaQueryRecord,
};

fn write_fixture(path: &std::path::Path) {
    let mut query = Vec::new();
    let seek = write_daa_query_record(&mut query, "query0", &[0, 1, 2], SequenceType::AminoAcid);
    query.extend_from_slice(&0u32.to_ne_bytes());
    query.push(compute_flag(30, 0, 0, false));
    query.push(30);
    query.push(0);
    query.push(0);
    query.push(PackedOperation::from_op_count(EditOperation::Match, 3).code);
    query.push(PackedOperation::terminator().code);
    finish_daa_query_record(&mut query, seek);
    query.extend_from_slice(&0u32.to_ne_bytes());

    let mut header = DaaHeader2::with_params(1, 100, 11, 1, 2, -3, 0.7, 1.2, 0.001, "BLOSUM62", 2);
    header.db_seqs_used = 1;
    header.query_records = 1;
    header.block_type[0] = BlockType::Alignments as u8;
    header.block_type[1] = BlockType::RefNames as u8;
    header.block_type[2] = BlockType::RefLengths as u8;
    header.block_size[0] = query.len() as u64;
    header.block_size[1] = 5;
    header.block_size[2] = 4;

    let mut file = File::create(path).unwrap();
    DaaHeader1::new().write_to(&mut file).unwrap();
    header.write_to(&mut file).unwrap();
    file.write_all(&query).unwrap();
    file.write_all(b"ref0\0").unwrap();
    file.write_all(&100u32.to_ne_bytes()).unwrap();
}

#[test]
fn native_view_writes_text_and_seekable_daa_to_stdout() {
    let directory =
        std::env::temp_dir().join(format!("diamond-view-stdout-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let input = directory.join("input.daa");
    let daa_stdout = directory.join("stdout.daa");
    write_fixture(&input);

    let text = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args([
            "view",
            "-a",
            input.to_str().unwrap(),
            "-f",
            "6",
            "qseqid",
            "sseqid",
            "score",
            "-p1",
        ])
        .output()
        .unwrap();
    assert!(
        text.status.success(),
        "{}",
        String::from_utf8_lossy(&text.stderr)
    );
    assert_eq!(text.stdout, b"query0\tref0\t30\n");

    let redirected = File::create(&daa_stdout).unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args(["view", "-a", input.to_str().unwrap(), "-f100", "-p1"])
        .stdout(Stdio::from(redirected))
        .status()
        .unwrap();
    assert!(status.success());

    let mut daa = DaaFile::open(&daa_stdout).unwrap();
    assert_eq!(daa.query_records(), 1);
    assert_eq!(daa.ref_name(0), "ref0");
    let (buffer, query_number) = daa.read_query_buffer().unwrap().unwrap();
    let record = DaaQueryRecord::from_buffer(&daa, &buffer, query_number).unwrap();
    assert_eq!(record.query_name, "query0");

    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn native_view_matches_upstream_help_and_console_verbosity_policy() {
    let directory =
        std::env::temp_dir().join(format!("diamond-view-presentation-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let input = directory.join("input.daa");
    write_fixture(&input);

    let help = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .arg("view")
        .output()
        .unwrap();
    assert!(!help.status.success());
    let help_stdout = String::from_utf8(help.stdout).unwrap();
    assert!(help_stdout.starts_with("Options:\n--threads"));
    assert!(help_stdout.contains("--outfmt                 output format\n"));
    assert!(help_stdout.contains(
        "\tnormalized_bitscore means Bitscore normalized by maximum self alignment score\n"
    ));
    assert!(help_stdout.contains(
        "\t        ssubvariety means Unique subject subvariety(s), separated by a ';'\n"
    ));
    assert!(help_stdout
        .ends_with("--forwardonly            only show alignments of forward strand\n\n"));
    assert!(String::from_utf8(help.stderr)
        .unwrap()
        .ends_with("Error: The view command requires a DAA (option -a) input file.\n"));

    let base = [
        "view",
        "-a",
        input.to_str().unwrap(),
        "-f",
        "6",
        "qseqid",
        "-p1",
    ];
    let stdout_default = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args(base)
        .output()
        .unwrap();
    assert!(stdout_default.status.success());
    assert_eq!(stdout_default.stdout, b"query0\n");
    assert!(stdout_default.stderr.is_empty());

    let stdout_verbose = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args(base)
        .arg("--verbose")
        .output()
        .unwrap();
    assert!(stdout_verbose.status.success());
    let verbose_stderr = String::from_utf8(stdout_verbose.stderr).unwrap();
    assert!(verbose_stderr.starts_with("diamond v2.1.24.178 "));
    assert!(verbose_stderr.contains("#CPU threads: 1\n"));
    assert!(verbose_stderr
        .contains("Scoring parameters: (Matrix=BLOSUM62 Lambda=0.267 K=0.041 Penalties=11/1)\n"));
    assert!(verbose_stderr.contains("Build version = 178\n"));
    assert!(verbose_stderr.contains("#Target sequences to report alignments for: 25\n"));
    assert!(!verbose_stderr.contains("DP fields:"));

    let output = directory.join("view.tsv");
    let file_default = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args(base)
        .args(["-o", output.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(file_default.status.success());
    let default_stderr = String::from_utf8(file_default.stderr).unwrap();
    assert!(default_stderr.starts_with("diamond v2.1.24.178 "));
    assert!(!default_stderr.contains("Build version ="));

    let file_quiet = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args(base)
        .args(["-o", output.to_str().unwrap(), "--quiet"])
        .output()
        .unwrap();
    assert!(file_quiet.status.success());
    assert!(file_quiet.stderr.is_empty());

    let file_log = Command::new(env!("CARGO_BIN_EXE_diamond"))
        .args(base)
        .args(["-o", output.to_str().unwrap(), "--quiet", "--log"])
        .output()
        .unwrap();
    assert!(file_log.status.success());
    let log_stderr = String::from_utf8(file_log.stderr).unwrap();
    assert!(log_stderr.contains(" --quiet --log\n"));
    assert!(log_stderr.contains("Build version = 178\n"));
    assert!(log_stderr.contains("TRACE_PT_MEMBUF=OFF\n"));
    assert!(log_stderr.contains("MAX_SHAPE_LEN=19 SEQ_MASK STRICT_BAND\n"));
    assert!(log_stderr.contains("Min length ratio: 0\n"));
    assert!(log_stderr.contains("#Target sequences to report alignments for: 25\n"));
    assert!(log_stderr.contains("DP fields: 0\n"));

    std::fs::remove_dir_all(directory).unwrap();
}
