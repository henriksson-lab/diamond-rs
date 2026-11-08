//! Native search headers, gzip and configured I/O against frozen source output.
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

fn database(root: &Path, directory: &Path, tax: bool) -> PathBuf {
    let path = directory.join(if tax { "tax.dmnd" } else { "plain.dmnd" });
    let input = root.join(if tax {
        "tests/data/search_blocks/taxonomy/reference.faa"
    } else {
        "tests/data/tabular_context/reference.faa"
    });
    let mut command = Command::new(env!("CARGO_BIN_EXE_diamond"));
    command
        .args(["makedb", "--in"])
        .arg(input)
        .arg("-d")
        .arg(&path)
        .arg("-p1");
    if tax {
        let data = root.join("tests/data/search_blocks/taxonomy");
        command
            .arg("--taxonmap")
            .arg(data.join("mapping.tsv"))
            .arg("--taxonnodes")
            .arg(data.join("nodes.dmp"))
            .arg("--taxonnames")
            .arg(data.join("names.dmp"));
    }
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    path
}
fn decode(bytes: &[u8], compressed: bool) -> Vec<u8> {
    if !compressed {
        return bytes.to_vec();
    }
    let mut result = Vec::new();
    flate2::read::MultiGzDecoder::new(bytes)
        .read_to_end(&mut result)
        .unwrap();
    result
}
fn expected_metadata(bytes: Vec<u8>, format: &str, invocation: &str, database: &str) -> Vec<u8> {
    if format == "edge" || format == "100" {
        return bytes;
    }
    let mut text = String::from_utf8(bytes).unwrap();
    for marker in ["# Invocation: ", "\tCL:"] {
        if let Some(start) = text.find(marker) {
            let begin = start + marker.len();
            let end = begin + text[begin..].find('\n').unwrap();
            text.replace_range(begin..end, invocation);
        }
    }
    if format == "xml" {
        let marker = "<BlastOutput_db>";
        let begin = text.find(marker).unwrap() + marker.len();
        let end = begin + text[begin..].find("</BlastOutput_db>").unwrap();
        text.replace_range(begin..end, database);
    }
    text.into_bytes()
}
fn wire(actual: &[u8], expected: &[u8], label: &str) {
    let offset = actual
        .iter()
        .zip(expected)
        .position(|(a, b)| a != b)
        .unwrap_or(actual.len().min(expected.len()));
    assert!(
        actual == expected,
        "{label}: first differing byte {offset}; lengths {} and {}",
        actual.len(),
        expected.len()
    );
}
#[test]
fn native_search_output_options_match_source_bytes_and_metadata() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let fixtures = root.join("tests/data/output_options");
    let dir = std::env::temp_dir().join(format!("diamond-output-options-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let plain = database(root, &dir, false);
    let taxonomy = database(root, &dir, true);
    let mut checked = 0;
    let mut encoded = 0;
    for row in std::fs::read_to_string(fixtures.join("cases.tsv"))
        .unwrap()
        .lines()
    {
        let c: Vec<_> = row.split('\t').collect();
        assert_eq!(c.len(), 10);
        let source = std::fs::read(fixtures.join(format!("{}.wire", c[0]))).unwrap();
        let source_decoded = decode(&source, c[6] == "1");
        for (threads, memory) in [("1", "16G"), ("4", "0G")] {
            let output = dir.join(if ["gzip_existing", "plain_named_gzip"].contains(&c[9]) {
                "output.gz"
            } else {
                "output"
            });
            let db = if c[3] == "tax" { &taxonomy } else { &plain };
            let mut command = Command::new(env!("CARGO_BIN_EXE_diamond"));
            command
                .arg(c[1])
                .arg("-q")
                .arg(root.join("tests/data/tabular_context").join(c[4]))
                .arg("-d")
                .arg(db)
                .args([
                    "--masking",
                    "0",
                    "--comp-based-stats",
                    "1",
                    "--max-hsps",
                    "0",
                    "-p",
                    threads,
                    "--memory-limit",
                    memory,
                    "-b",
                    c[5],
                    "-f",
                    c[2],
                ]);
            if ["6", "104"].contains(&c[2]) {
                if c[3] == "tax" {
                    command.args(["qseqid", "staxids", "sscinames", "sfamily", "slineages"]);
                } else {
                    command.args([
                        "qseqid",
                        "sseqid",
                        "qseq",
                        "sseq",
                        "cigar",
                        "full_qseq",
                        "full_qqual",
                        "hspnum",
                        "bitscore",
                    ]);
                }
            }
            if c[2] == "102" {
                command.arg("--include-lineage");
            }
            match c[8] {
                "absent" => (),
                "bare" => {
                    command.arg("--header");
                }
                "repeat" => {
                    command.args(["--header", "simple", "--header"]);
                }
                value => {
                    command.args(["--header", value]);
                }
            }
            if c[6] != "none" {
                command.args(["--compress", c[6]]);
            }
            if c[7] != "default" {
                command.args(["--file-buffer-size", c[7]]);
            }
            if ["gzip_no_append", "daa_no_append"].contains(&c[9]) {
                command.arg("--no-auto-append");
            }
            command.arg("-o").arg(&output);
            let invocation = std::iter::once(command.get_program())
                .chain(command.get_args())
                .map(|s| s.to_string_lossy())
                .collect::<Vec<_>>()
                .join(" ");
            let result = command.output().unwrap();
            assert!(
                result.status.success(),
                "{}: {}",
                c[0],
                String::from_utf8_lossy(&result.stderr)
            );
            let actual_path = if c[2] == "100" && c[9] != "daa_no_append" {
                PathBuf::from(format!("{}.daa", output.display()))
            } else if c[6] == "1"
                && c[9] != "gzip_no_append"
                && !output.to_string_lossy().ends_with(".gz")
            {
                PathBuf::from(format!("{}.gz", output.display()))
            } else {
                output.clone()
            };
            let actual = std::fs::read(actual_path).unwrap();
            let decoded = decode(&actual, c[6] == "1");
            let expected = expected_metadata(
                source_decoded.clone(),
                c[2],
                &invocation,
                &db.to_string_lossy(),
            );
            wire(&decoded, &expected, c[0]);
            if expected == source_decoded {
                wire(&actual, &source, c[0]);
                encoded += 1;
            }
            if threads == "1" && ["case-001", "case-081", "case-154"].contains(&c[0]) {
                let args = std::iter::once(command.get_program())
                    .chain(command.get_args())
                    .map(|s| s.to_string_lossy().into_owned())
                    .collect::<Vec<_>>();
                assert_eq!(
                    diamond::run(&args.iter().map(String::as_str).collect::<Vec<_>>()),
                    0
                );
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 348);
    assert!(encoded > 100);
    std::fs::remove_dir_all(dir).unwrap();
}

#[test]
fn large_gzip_stream_matches_source_across_codec_scratch_boundaries() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let expected = std::fs::read(root.join("tests/data/output_options/large.gzip")).unwrap();
    assert!(expected.len() > 1 << 20);
    let mut state = 0x6a09e667f3bcc909u64;
    let alphabet = b"ARNDCQEGHILKMFPSTWYV";
    let mut raw = Vec::with_capacity(2_401_200);
    for _ in 0..1200 {
        for _ in 0..2000 {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            raw.push(alphabet[((state >> 32) % 20) as usize]);
        }
        raw.push(b'\n');
    }
    wire(&decode(&expected, true), &raw, "source large generator");
    let directory =
        std::env::temp_dir().join(format!("diamond-output-large-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    for buffer in [65536, 67_108_864] {
        let output = directory.join("output.gz");
        let mut sink = diamond::util::io::OutputFile::with_buffer_size(
            output.to_str().unwrap(),
            diamond::util::io::Compressor::Zlib,
            "wb",
            buffer,
        )
        .unwrap();
        for chunk in raw.chunks(2001) {
            sink.write_raw(chunk).unwrap();
        }
        sink.close().unwrap();
        wire(
            &std::fs::read(output).unwrap(),
            &expected,
            &format!("large/buffer{buffer}"),
        );
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn output_option_protocol_and_file_creation_match_source() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let dir = std::env::temp_dir().join(format!("diamond-output-errors-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let db = database(root, &dir, false);
    for (index, (format, flags, error, created, suffix)) in [
        (
            "6",
            vec!["--header", "bad"],
            "Invalid value for string field: bad. Permitted values: 0, simple, verbose",
            true,
            "",
        ),
        (
            "6",
            vec!["--header", "simple", "verbose"],
            "Invalid header format: simple verbose",
            true,
            "",
        ),
        (
            "6",
            vec!["--header=simple"],
            "Invalid option: header=simple",
            false,
            "",
        ),
        (
            "6",
            vec!["--compress"],
            "Invalid parameter count for option '--compress'",
            false,
            "",
        ),
        (
            "6",
            vec!["--compress", "1", "0"],
            "Invalid parameter count for option '--compress'",
            false,
            "",
        ),
        (
            "6",
            vec!["--compress=1"],
            "Invalid option: compress=1",
            false,
            "",
        ),
        (
            "6",
            vec!["--compress", "bad"],
            "Invalid compression algorithm: bad",
            false,
            "",
        ),
        (
            "6",
            vec!["--compress", "zstd"],
            "Executable was not compiled with ZStd support.",
            true,
            ".zst",
        ),
        (
            "6",
            vec!["--file-buffer-size"],
            "Invalid parameter count for option '--file-buffer-size'",
            false,
            "",
        ),
        (
            "6",
            vec!["--file-buffer-size", "7", "8"],
            "Invalid parameter count for option '--file-buffer-size'",
            false,
            "",
        ),
        (
            "6",
            vec!["--file-buffer-size=7"],
            "Invalid option: file-buffer-size=7",
            false,
            "",
        ),
        (
            "6",
            vec!["--no-auto-append", "1"],
            "Invalid parameter count for option '--no-auto-append'",
            false,
            "",
        ),
        (
            "6",
            vec!["--no-auto-append=1"],
            "Invalid option: no-auto-append=1",
            false,
            "",
        ),
        (
            "100",
            vec!["--compress", "0"],
            "Compression is not supported for DAA format.",
            false,
            ".daa",
        ),
        (
            "100",
            vec!["--compress", "1"],
            "Compression is not supported for DAA format.",
            false,
            ".daa.gz",
        ),
        (
            "100",
            vec!["--compress", "zstd"],
            "Compression is not supported for DAA format.",
            false,
            ".daa.zst",
        ),
        (
            "6",
            vec!["--file-buffer-size", "0"],
            "Invalid peek",
            false,
            "",
        ),
        (
            "6",
            vec!["--file-buffer-size", "1"],
            "Invalid peek",
            false,
            "",
        ),
        (
            "6",
            vec!["--file-buffer-size", "2"],
            "Invalid peek",
            false,
            "",
        ),
        (
            "6",
            vec!["--file-buffer-size", "3"],
            "Invalid peek",
            false,
            "",
        ),
    ]
    .into_iter()
    .enumerate()
    {
        let base = dir.join(format!("error-{index}"));
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("blastp")
            .arg("-q")
            .arg(root.join("tests/data/tabular_context/query.faa"))
            .arg("-d")
            .arg(&db)
            .args(["-f", format])
            .args(flags)
            .arg("-o")
            .arg(&base)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(
            String::from_utf8_lossy(&result.stderr).ends_with(&format!("Error: {error}\n")),
            "{index}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let file = PathBuf::from(format!("{}{suffix}", base.display()));
        assert_eq!(file.exists(), created, "{index}");
        if created {
            assert!(std::fs::read(file).unwrap().is_empty());
        }
    }
    // The source catches bad_alloc for impossible positive capacities as well
    // as negative ones. The third case checks the signed parser's saturation.
    // On 32-bit targets also reject a positive value that would wrap to zero.
    for value in ["-1", "9223372036854775807", "9223372036854775808"]
        .into_iter()
        .chain((usize::BITS == 32).then_some("4294967296"))
    {
        let output = dir.join(format!("allocation-{value}"));
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("blastp")
            .arg("-q")
            .arg(root.join("tests/data/tabular_context/query.faa"))
            .arg("-d")
            .arg(&db)
            .args(["--file-buffer-size", value])
            .arg("-o")
            .arg(&output)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1), "buffer {value}");
        assert!(
            String::from_utf8_lossy(&result.stderr).ends_with("Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage.\n"),
            "buffer {value}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(!output.exists(), "buffer {value} created output");
    }
    std::fs::remove_dir_all(dir).unwrap();
}
