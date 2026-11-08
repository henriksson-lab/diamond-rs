//! Snapshots produced by pinned upstream DIAMOND 2.1.24, not the Rust graph
//! implementation. The executable must reach native GVC without an FFI build.

#[path = "support/source_startup.rs"]
mod source_startup;

use std::path::Path;
use std::process::Command;

#[test]
fn native_gvc_matches_upstream_graph_modes_and_ordering() {
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gvc");
    let root = std::env::temp_dir().join(format!("diamond-native-gvc-{}", std::process::id()));
    std::fs::create_dir_all(&root).unwrap();
    for (mode, options) in [
        ("default", vec![]),
        ("cover50", vec!["--member-cover", "50"]),
        ("cover95", vec!["--member-cover", "95"]),
        ("no-reassign", vec!["--no-reassign"]),
        ("depth1", vec!["--connected-component-depth", "1"]),
        ("depth2", vec!["--connected-component-depth", "2"]),
        ("triplet", vec!["--edge-format", "triplet"]),
        ("symmetric", vec!["--edge-format", "triplet", "--symmetric"]),
    ] {
        let triplet = mode == "triplet" || mode == "symmetric";
        let output = root.join(format!("{mode}.tsv"));
        let centroids = root.join(format!("{mode}.centroids"));
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("greedy-vertex-cover")
            .arg("-d")
            .arg(fixtures.join("mapping.tsv"))
            .arg("--edges")
            .arg(fixtures.join(if triplet {
                "triplet-edges.tsv"
            } else {
                "edges.tsv"
            }))
            .arg("-o")
            .arg(&output)
            .arg("--centroid-out")
            .arg(&centroids)
            .args(["--quiet", "-p1"])
            .args(options)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{mode}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(source_startup::stdout_body(&result), b"", "{mode}: stdout");
        assert_eq!(
            result.stderr,
            source_startup::age_warning_stderr(),
            "{mode}: quiet diagnostics"
        );
        for (actual, extension) in [(&output, "tsv"), (&centroids, "centroids")] {
            assert_eq!(
                std::fs::read(actual).unwrap(),
                std::fs::read(fixtures.join(format!("{mode}.{extension}"))).unwrap(),
                "{mode}.{extension}"
            );
        }
    }
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn native_gvc_uses_upstream_option_cardinality_and_symmetry_checks() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gvc");
    for (options, expected) in [
        (
            vec!["--member-cover"],
            "Invalid parameter count for option '--member-cover'",
        ),
        (
            vec!["--symmetric", "1"],
            "Invalid parameter count for option '--symmetric'",
        ),
        (vec!["--member-cover=90"], "Invalid option: member-cover=90"),
        (
            vec!["--symmetric"],
            "--symmetric requires triplet edge format",
        ),
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_diamond"))
            .arg("greedy-vertex-cover")
            .arg("-d")
            .arg(fixture.join("mapping.tsv"))
            .arg("--edges")
            .arg(fixture.join("edges.tsv"))
            .args(options)
            .output()
            .unwrap();
        assert!(!result.status.success());
        let stderr = String::from_utf8_lossy(&result.stderr);
        assert!(stderr.contains(expected), "{stderr}");
        assert!(!stderr.contains("FFI"), "{stderr}");
    }
}
