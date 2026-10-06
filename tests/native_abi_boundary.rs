use std::path::{Path, PathBuf};

const ALLOWED_NATIVE_ABI_FILES: &[&str] = &["src/ffi.rs"];

fn rust_files(root: &Path, files: &mut Vec<PathBuf>) {
    for entry in std::fs::read_dir(root).expect("read source directory") {
        let entry = entry.expect("read source entry");
        let path = entry.path();
        if path.is_dir() {
            rust_files(&path, files);
        } else if path.extension().is_some_and(|extension| extension == "rs") {
            files.push(path);
        }
    }
}

#[test]
fn handwritten_native_abi_is_confined_to_reviewed_compatibility_modules() {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let source = manifest.join("src");
    let mut files = Vec::new();
    rust_files(&source, &mut files);

    let mut violations = Vec::new();
    for path in files {
        let relative = path.strip_prefix(manifest).unwrap();
        let relative = relative.to_string_lossy().replace('\\', "/");
        let text = std::fs::read_to_string(&path).expect("read Rust source");
        for (line_index, line) in text.lines().enumerate() {
            let code = line.trim_start();
            if code.starts_with("//") {
                continue;
            }
            // Match visibility-qualified declarations and exported functions
            // too; checking only line prefixes would let `pub extern` bypass
            // the guard accidentally.
            let declares_native_abi = code.contains("extern \"")
                || code.contains("#[link(")
                || code.contains("#[link_name")
                || code.contains("#[link_ordinal")
                || code.contains("#[no_mangle")
                || code.contains("#[unsafe(no_mangle")
                || code.contains("#[export_name")
                || code.contains("#[unsafe(export_name");
            if declares_native_abi && !ALLOWED_NATIVE_ABI_FILES.contains(&relative.as_str()) {
                violations.push(format!("{relative}:{}: {code}", line_index + 1));
            }
        }
    }

    assert!(
        violations.is_empty(),
        "new handwritten native ABI declarations must be placed in an explicitly reviewed \
         compatibility module:\n{}",
        violations.join("\n")
    );
}
