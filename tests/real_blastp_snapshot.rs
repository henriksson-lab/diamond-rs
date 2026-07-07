use std::collections::BTreeMap;
use std::fs;

const CPP_MOUSE_500Q: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/cpp_blastp_mouse_500q.tsv"
);

#[derive(Debug, Clone, PartialEq, Eq)]
struct BlastRow {
    qseqid: String,
    sseqid: String,
    pident: String,
    length: u32,
    mismatch: u32,
    gapopen: u32,
    qstart: u32,
    qend: u32,
    sstart: u32,
    send: u32,
    evalue: String,
    bitscore: String,
}

impl BlastRow {
    fn from_line(line: &str, line_no: usize) -> Self {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            cols.len(),
            12,
            "line {} should have BLAST tabular outfmt 6 columns: {}",
            line_no,
            line
        );
        Self {
            qseqid: cols[0].to_string(),
            sseqid: cols[1].to_string(),
            pident: cols[2].to_string(),
            length: cols[3].parse().expect("length"),
            mismatch: cols[4].parse().expect("mismatch"),
            gapopen: cols[5].parse().expect("gapopen"),
            qstart: cols[6].parse().expect("qstart"),
            qend: cols[7].parse().expect("qend"),
            sstart: cols[8].parse().expect("sstart"),
            send: cols[9].parse().expect("send"),
            evalue: cols[10].to_string(),
            bitscore: cols[11].to_string(),
        }
    }

    fn bitscore_value(&self) -> f64 {
        self.bitscore.parse().expect("bitscore")
    }

    fn evalue_value(&self) -> f64 {
        self.evalue.parse().expect("evalue")
    }
}

fn rows() -> Vec<BlastRow> {
    fs::read_to_string(CPP_MOUSE_500Q)
        .expect("read C++ mouse blastp snapshot")
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.is_empty())
        .map(|(i, line)| BlastRow::from_line(line, i + 1))
        .collect()
}

fn find_row<'a>(rows: &'a [BlastRow], qseqid: &str, sseqid: &str) -> &'a BlastRow {
    rows.iter()
        .find(|row| row.qseqid == qseqid && row.sseqid == sseqid)
        .unwrap_or_else(|| panic!("missing row {qseqid}/{sseqid}"))
}

fn assert_row_exact(row: &BlastRow, expected: [&str; 12]) {
    assert_eq!(row.qseqid, expected[0]);
    assert_eq!(row.sseqid, expected[1]);
    assert_eq!(row.pident, expected[2]);
    assert_eq!(row.length.to_string(), expected[3]);
    assert_eq!(row.mismatch.to_string(), expected[4]);
    assert_eq!(row.gapopen.to_string(), expected[5]);
    assert_eq!(row.qstart.to_string(), expected[6]);
    assert_eq!(row.qend.to_string(), expected[7]);
    assert_eq!(row.sstart.to_string(), expected[8]);
    assert_eq!(row.send.to_string(), expected[9]);
    assert_eq!(row.evalue, expected[10]);
    assert_eq!(row.bitscore, expected[11]);
}

#[test]
fn cpp_mouse_snapshot_has_expected_shape() {
    let rows = rows();
    assert_eq!(rows.len(), 6242);
    assert_row_exact(
        &rows[0],
        [
            "Gm10778", "Gm10778", "100", "113", "0", "0", "1", "113", "1", "113", "1.19e-76", "248",
        ],
    );
}

#[test]
fn cpp_mouse_snapshot_pins_known_scoring_divergence_rows() {
    let rows = rows();
    for expected in [
        [
            "Chd5", "Chd5", "100", "1946", "0", "0", "1", "1946", "1", "1946", "0.0", "3624",
        ],
        [
            "Chd5", "Chd3", "70.3", "1857", "398", "23", "133", "1908", "183", "1967", "0.0",
            "2444",
        ],
        [
            "Chd5",
            "Chd1",
            "42.5",
            "717",
            "335",
            "14",
            "512",
            "1205",
            "315",
            "977",
            "2.36e-155",
            "523",
        ],
        [
            "Chd5",
            "Chd1l",
            "43.3",
            "563",
            "269",
            "12",
            "689",
            "1232",
            "29",
            "560",
            "2.65e-121",
            "406",
        ],
        [
            "Chd5", "Btaf1", "28.6", "594", "317", "19", "702", "1226", "1265", "1820", "8.52e-50",
            "195",
        ],
        [
            "Olfr1490",
            "Olfr1490",
            "100",
            "321",
            "0",
            "0",
            "1",
            "321",
            "1",
            "321",
            "8.95e-232",
            "630",
        ],
        [
            "Olfr1024",
            "Olfr1024",
            "100",
            "327",
            "0",
            "0",
            "1",
            "327",
            "1",
            "327",
            "6.88e-239",
            "649",
        ],
        [
            "Obscn", "Obscn", "100", "8886", "0", "0", "1", "8886", "1", "8886", "0.0", "17357",
        ],
        [
            "Ttn", "Ttn", "100", "35213", "0", "0", "1", "35213", "1", "35213", "0.0", "64489",
        ],
    ] {
        assert_row_exact(find_row(&rows, expected[0], expected[1]), expected);
    }
}

#[test]
fn cpp_mouse_snapshot_perfect_self_hits_are_well_formed() {
    let rows = rows();
    let self_hits: Vec<&BlastRow> = rows
        .iter()
        .filter(|row| {
            row.qseqid == row.sseqid && row.pident == "100" && row.mismatch == 0 && row.gapopen == 0
        })
        .collect();
    assert!(
        self_hits.len() >= 400,
        "expected many real perfect self hits, saw {}",
        self_hits.len()
    );

    for row in self_hits {
        assert_eq!(
            row.qstart, row.sstart,
            "{}/{} qstart/sstart",
            row.qseqid, row.sseqid
        );
        assert_eq!(
            row.qend, row.send,
            "{}/{} qend/send",
            row.qseqid, row.sseqid
        );
        assert_eq!(
            row.qend - row.qstart + 1,
            row.length,
            "{}/{} aligned span",
            row.qseqid,
            row.sseqid
        );
        assert!(
            row.evalue_value().is_finite(),
            "{}/{} evalue",
            row.qseqid,
            row.sseqid
        );
        assert!(
            row.bitscore_value() > 0.0,
            "{}/{} bitscore",
            row.qseqid,
            row.sseqid
        );
    }
}

#[test]
fn cpp_mouse_snapshot_query_groups_are_ranked_and_capped() {
    let rows = rows();
    let mut groups: BTreeMap<&str, (usize, f64)> = BTreeMap::new();

    for row in &rows {
        let evalue = row.evalue_value();
        let entry = groups.entry(&row.qseqid).or_insert((0usize, f64::INFINITY));
        assert!(
            entry.0 < 25,
            "query {} has more than 25 reported rows",
            row.qseqid
        );
        if entry.1.is_finite() {
            assert!(
                evalue >= entry.1,
                "query {} is not sorted by evalue: {} after {} at subject {}",
                row.qseqid,
                evalue,
                entry.1,
                row.sseqid
            );
        }
        *entry = (entry.0 + 1, evalue);
    }

    assert!(groups.len() >= 200, "expected many real query groups");
}
