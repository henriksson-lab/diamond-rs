#!/usr/bin/env python3
"""Check source inventory claims; this is a developer audit, not a parity test.

Requires the adjacent upstream checkout. Does not affect Cargo/package builds.
Reports missing translations and unmapped implementation files, including .cc
files that the historical .cpp/.c-only count missed.
"""

import argparse
import csv
import json
from pathlib import Path
import re


def audit(repo):
    source = repo / "diamond/src"
    if not source.is_dir():
        raise SystemExit("Inventory auditing requires the adjacent diamond/ source checkout.")
    rows = list(csv.DictReader((repo / "translation/file_map.csv").open()))
    originals = {
        p.relative_to(repo).as_posix()
        for p in source.rglob("*")
        if p.suffix in {".c", ".cc", ".cpp", ".cxx"}
    }
    mapped = {row["original_file"] for row in rows}
    missing = []
    for row in rows:
        files = re.findall(r"\bsrc/[A-Za-z0-9_./-]+\.rs", row["rust_file"])
        absent = [p for p in files if not (repo / p).is_file()]
        if absent:
            missing.append({"original": row["original_file"], "status": row["status"], "missing_rust_files": absent})
    return {
        "original_implementation_files": len(originals),
        "mapped_original_files": len(mapped),
        "unmapped_originals": sorted(originals - mapped),
        "missing_originals": sorted(mapped - originals),
        "missing_rust_modules": missing,
        "pending_audits": [row["original_file"] for row in rows if row["status"] not in {"translated", "parity-tested", "replaced"}],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    args = parser.parse_args()
    report = audit(args.repo.resolve())
    print(json.dumps(report, indent=2))
    return int(any(report[key] for key in ["unmapped_originals", "missing_originals", "missing_rust_modules", "pending_audits"]))


if __name__ == "__main__":
    raise SystemExit(main())
