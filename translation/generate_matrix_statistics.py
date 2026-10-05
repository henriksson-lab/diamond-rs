#!/usr/bin/env python3
"""Extract StandardMatrix statistical tables from the vendored C++ initializers."""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "diamond/src/stats/matrices"
OUTPUT = ROOT / "src/stats/matrix_data.rs"
MATRICES = ["blosum45", "blosum50", "blosum62", "blosum80", "blosum90", "pam250", "pam30", "pam70"]


def balanced_body(text: str, start: int) -> str:
    depth = 0
    for index in range(start, len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1:index]
    raise ValueError("unbalanced initializer")


def fields(body: str) -> list[str]:
    result, start, depth = [], 0, 0
    for index, char in enumerate(body):
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
        elif char == "," and depth == 0:
            result.append(body[start:index].strip())
            start = index + 1
    tail = body[start:].strip()
    if tail:
        result.append(tail)
    return result


def numbers(field: str) -> list[str]:
    return re.findall(r"(?<![A-Za-z_])[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?", field)


def rust_array(name: str, values: list[str]) -> str:
    rows = []
    for offset in range(0, len(values), 6):
        rows.append("    " + ", ".join(values[offset:offset + 6]) + ",")
    return f"pub const {name}: [f64; {len(values)}] = [\n" + "\n".join(rows) + "\n];\n"


def rust_matrix(name: str, values: list[str], width: int) -> str:
    rows = ["    [" + ", ".join(values[offset:offset + width]) + "]," for offset in range(0, len(values), width)]
    return f"pub const {name}: [[f64; {width}]; {len(values) // width}] = [\n" + "\n".join(rows) + "\n];\n"


parts = [
    "//! Generated from `diamond/src/stats/matrices/*.h`; do not edit by hand.\n",
    "use super::matrices;\nuse super::standard_matrix::StandardMatrix;\n",
    "pub struct MatrixStatistics {\n"
    "    pub joint_probs: &'static [f64; 20 * 20],\n"
    "    pub background_freqs: &'static [f64; 20],\n"
    "    pub freq_ratios: &'static [[f64; 28]; 28],\n"
    "}\n",
]

for matrix in MATRICES:
    text = (SOURCE / f"{matrix}.h").read_text()
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    marker = re.search(rf"const\s+StandardMatrix\s+{matrix}\s*\{{", text)
    if marker is None:
        raise ValueError(f"initializer not found: {matrix}")
    body = balanced_body(text, marker.end() - 1)
    initializers = fields(body)
    if len(initializers) != 7:
        raise ValueError(f"{matrix}: expected 7 fields, got {len(initializers)}")
    joint, background, ratios = map(numbers, initializers[4:7])
    if (len(joint), len(background), len(ratios)) != (400, 20, 784):
        raise ValueError(f"{matrix}: bad table sizes {(len(joint), len(background), len(ratios))}")
    prefix = matrix.upper()
    parts.extend([
        rust_array(f"{prefix}_JOINT_PROBS", joint),
        rust_array(f"{prefix}_BACKGROUND_FREQS", background),
        rust_matrix(f"{prefix}_FREQ_RATIOS", ratios, 28),
    ])

parts.append("pub fn get(matrix: &'static StandardMatrix) -> MatrixStatistics {\n")
for matrix in MATRICES:
    prefix = matrix.upper()
    parts.append(
        f"    if matrix.scores == matrices::{prefix}.scores {{\n"
        f"        return MatrixStatistics {{ joint_probs: &{prefix}_JOINT_PROBS, "
        f"background_freqs: &{prefix}_BACKGROUND_FREQS, freq_ratios: &{prefix}_FREQ_RATIOS }};\n"
        "    }\n"
    )
parts.append('    unreachable!("unknown static StandardMatrix")\n}\n')

OUTPUT.write_text("\n".join(parts))
