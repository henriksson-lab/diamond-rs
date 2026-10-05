#!/usr/bin/env python3
"""Small, dependency-free semantic comparators for VALIDATION_PLAN.md."""

from __future__ import annotations

import argparse
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


def fail(message: str) -> None:
    print(message, file=sys.stderr)
    raise SystemExit(1)


def bytes_of(path: str) -> bytes:
    return Path(path).read_bytes()


def compare_exact(left: str, right: str) -> None:
    a, b = bytes_of(left), bytes_of(right)
    if a != b:
        offset = next((i for i, pair in enumerate(zip(a, b)) if pair[0] != pair[1]), min(len(a), len(b)))
        fail(f"byte mismatch at offset {offset}: lengths {len(a)} != {len(b)}")


def normalized_lines(path: str) -> list[bytes]:
    return sorted(line.rstrip(b"\r") for line in bytes_of(path).splitlines())


def compare_sorted_lines(left: str, right: str) -> None:
    a, b = normalized_lines(left), normalized_lines(right)
    if a != b:
        fail(f"sorted records differ: {len(a)} != {len(b)} lines")


def load_json_stream(path: str):
    text = Path(path).read_text(encoding="utf-8")
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        values = []
        for number, line in enumerate(text.splitlines(), 1):
            if line.strip():
                try:
                    values.append(json.loads(line))
                except json.JSONDecodeError as error:
                    fail(f"{path}:{number}: invalid JSON: {error}")
        return values


def compare_json(left: str, right: str) -> None:
    if load_json_stream(left) != load_json_stream(right):
        fail("parsed JSON values differ")


def normalize_sam(path: str) -> list[str]:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    # Program version and full invocation necessarily name different binaries.
    return [line for line in lines if not line.startswith("@PG\t")]


def compare_sam(left: str, right: str) -> None:
    if normalize_sam(left) != normalize_sam(right):
        fail("SAM differs after removing @PG program metadata")


VOLATILE_XML_TAGS = {"BlastOutput_program", "BlastOutput_version", "BlastOutput_reference"}


def xml_value(node: ET.Element):
    tag = node.tag.rsplit("}", 1)[-1]
    if tag in VOLATILE_XML_TAGS:
        return None
    children = [value for child in node if (value := xml_value(child)) is not None]
    return tag, (node.text or "").strip(), tuple(sorted(node.attrib.items())), tuple(children)


def compare_xml(left: str, right: str) -> None:
    try:
        a = xml_value(ET.parse(left).getroot())
        b = xml_value(ET.parse(right).getroot())
    except ET.ParseError as error:
        fail(f"invalid XML: {error}")
    if a != b:
        fail("XML structure differs after removing version/reference metadata")


def fasta_ids(path: str) -> list[str]:
    result = []
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith(">"):
                result.append(line[1:].split()[0])
    return result


def cluster_pairs(path: str) -> list[tuple[str, str]]:
    result = []
    with Path(path).open(encoding="utf-8") as stream:
        for number, line in enumerate(stream, 1):
            fields = line.rstrip("\r\n").split("\t")
            if not line.strip():
                continue
            if len(fields) != 2 or not all(fields):
                fail(f"{path}:{number}: expected two nonempty tab-separated fields")
            result.append((fields[0], fields[1]))
    return result


def cluster_stats(path: str, fasta: str | None) -> dict[str, int | str]:
    pairs = cluster_pairs(path)
    members = [member for _, member in pairs]
    counts = Counter(members)
    duplicates = sorted(member for member, count in counts.items() if count != 1)
    if duplicates:
        fail(f"{path}: members occurring other than once: {duplicates[:5]}")
    if fasta:
        expected = Counter(fasta_ids(fasta))
        if counts != expected:
            fail(f"{path}: cluster members do not equal FASTA identifiers")
    sizes = Counter(centroid for centroid, _ in pairs)
    import hashlib
    digest = hashlib.sha256(
        b"".join(f"{centroid}\t{member}\n".encode() for centroid, member in sorted(pairs))
    ).hexdigest()
    return {
        "pairs": len(pairs),
        "clusters": len(sizes),
        "singletons": sum(size == 1 for size in sizes.values()),
        "largest": max(sizes.values(), default=0),
        "sorted_sha256": digest,
    }


def compare_cluster(left: str, right: str, fasta: str | None) -> None:
    a = cluster_stats(left, fasta)
    b = cluster_stats(right, fasta)
    print(json.dumps({"left": a, "right": b}, sort_keys=True))
    if sorted(cluster_pairs(left)) != sorted(cluster_pairs(right)):
        fail("sorted centroid/member pairs differ")


def essential_error(path: str) -> str:
    text = Path(path).read_text(encoding="utf-8", errors="replace").lower()
    text = re.sub(r"diamond v?\d[^\n]*", "", text)
    text = re.sub(r"(?:[a-z]:)?[/\\][^\s:]+", "<path>", text)
    return " ".join(text.split())


def compare_errors(left: str, right: str) -> None:
    a, b = essential_error(left), essential_error(right)
    gap_terms = ("gap", "penalt", "support", "matrix")
    if not a or not b or not any(term in a and term in b for term in gap_terms):
        fail(f"essential diagnostics differ: {a!r} vs {b!r}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=["exact", "sorted-lines", "json", "sam", "xml", "cluster", "errors"])
    parser.add_argument("left")
    parser.add_argument("right")
    parser.add_argument("--fasta")
    args = parser.parse_args()
    if args.kind == "exact":
        compare_exact(args.left, args.right)
    elif args.kind == "sorted-lines":
        compare_sorted_lines(args.left, args.right)
    elif args.kind == "json":
        compare_json(args.left, args.right)
    elif args.kind == "sam":
        compare_sam(args.left, args.right)
    elif args.kind == "xml":
        compare_xml(args.left, args.right)
    elif args.kind == "cluster":
        compare_cluster(args.left, args.right, args.fasta)
    else:
        compare_errors(args.left, args.right)


if __name__ == "__main__":
    main()
