These graph regression snapshots were generated with the adjacent upstream
DIAMOND 2.1.24 executable (`diamond/build/diamond`, baseline commit
`1d162b4fefb5e6b4d868c24e6ad551cad2dcf246`) on 2026-11-06.

Each snapshot uses `greedy-vertex-cover -d mapping.tsv --edges edges.tsv
-p1 --quiet -o MODE.tsv --centroid-out MODE.centroids`, with the additional
options recorded in `tests/native_gvc.rs`. Triplet and symmetric cases use
`triplet-edges.tsv` instead. Both output files are compared byte-for-byte.

The deterministic synthetic graph covers directional coverage, cutoff
boundaries, self edges, equal weights, isolated vertices, reassignment,
connected components, and symmetric/directed triplet interpretation. These
are algorithm tests, not biological performance benchmarks.
