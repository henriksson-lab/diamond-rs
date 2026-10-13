#!/usr/bin/env python3
"""Prepare a contiguous natural UniProt FASTA/taxid benchmark subset.

The input records are copied verbatim and never duplicated.  The mapping is
derived from the accession in the UniProt identifier and the record's OX
field, so the paired files exercise DIAMOND's real accession parser and
taxonomy join without inventing sequences or taxonomic assignments.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import TextIO


OX = re.compile(r"(?:^|\s)OX=(\d+)(?:\s|$)")


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def uniprot_accession(header: str) -> str:
    identifier = header[1:].split(None, 1)[0]
    fields = identifier.split("|")
    if len(fields) >= 3 and fields[0] in {"sp", "tr"}:
        return fields[1]
    return identifier


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_fasta", type=Path)
    parser.add_argument("output_mapping", type=Path)
    parser.add_argument("--records", type=int, default=150_000)
    args = parser.parse_args()
    if args.records <= 0:
        parser.error("--records must be positive")

    count = 0
    with (
        open_text(args.input) as source,
        args.output_fasta.open("wt", encoding="utf-8", newline="") as fasta,
        args.output_mapping.open("wt", encoding="utf-8", newline="") as mapping,
    ):
        mapping.write("accession.version\ttaxid\n")
        for line in source:
            if line.startswith(">"):
                if count == args.records:
                    break
                match = OX.search(line)
                if match is None:
                    raise ValueError(f"record {count + 1} has no OX taxid: {line.rstrip()}")
                mapping.write(f"{uniprot_accession(line)}\t{match.group(1)}\n")
                count += 1
            fasta.write(line)

    if count != args.records:
        raise ValueError(f"requested {args.records} records but input contained {count}")


if __name__ == "__main__":
    main()
