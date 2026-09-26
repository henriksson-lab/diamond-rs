//! Regression case inventory from `diamond/src/test/test_cases.cpp`.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TestCase {
    pub description: &'static str,
    pub command_line: &'static str,
}

pub const TEST_CASES: &[TestCase] = &[
    TestCase {
        description: "blastp (default)",
        command_line: "blastp -p1",
    },
    TestCase {
        description: "blastp (multithreaded)",
        command_line: "blastp -p4",
    },
    TestCase {
        description: "blastp (blocked)",
        command_line: "blastp -c1 -b0.00002 -p4",
    },
    TestCase {
        description: "blastp (more-sensitive)",
        command_line: "blastp --more-sensitive -c1 -p4",
    },
    TestCase {
        description: "blastp (very-sensitive)",
        command_line: "blastp --very-sensitive -c1 -p4",
    },
    TestCase {
        description: "blastp (ultra-sensitive)",
        command_line: "blastp --ultra-sensitive -c1 -p4",
    },
    TestCase {
        description: "blastp (max-hsps)",
        command_line: "blastp --more-sensitive -c1 -p4 --max-hsps 0",
    },
    TestCase {
        description: "blastp (target-parallel)",
        command_line: "blastp --more-sensitive -c1 -p4 --query-parallel-limit 1",
    },
    TestCase {
        description: "blastp (query-indexed)",
        command_line: "blastp --more-sensitive -c1 -p4 --algo 1",
    },
    TestCase {
        description: "blastp (comp-based-stats 0)",
        command_line: "blastp --more-sensitive -c1 -p4 --comp-based-stats 0",
    },
    TestCase {
        description: "blastp (comp-based-stats 2)",
        command_line: "blastp --more-sensitive -c1 -p4 --comp-based-stats 2",
    },
    TestCase {
        description: "blastp (comp-based-stats 3)",
        command_line: "blastp --more-sensitive -c1 -p4 --comp-based-stats 3",
    },
    TestCase {
        description: "blastp (comp-based-stats 4)",
        command_line: "blastp --more-sensitive -c1 -p4 --comp-based-stats 4",
    },
    TestCase {
        description: "blastp (target seqs)",
        command_line: "blastp -k3 -c1 -p4",
    },
    TestCase {
        description: "blastp (top)",
        command_line: "blastp --top 10 -p4",
    },
    TestCase {
        description: "blastp (evalue)",
        command_line: "blastp -e10000 --more-sensitive -c1 -p4",
    },
    TestCase {
        description: "blastp (blosum50)",
        command_line: "blastp --matrix blosum50 -p4",
    },
    TestCase {
        description: "blastp (pairwise format)",
        command_line: "blastp -c1 -f0 -p4",
    },
    TestCase {
        description: "blastp (XML format)",
        command_line: "blastp -c1 -f xml -p4",
    },
    TestCase {
        description: "blastp (PAF format)",
        command_line: "blastp -c1 -f paf -p1",
    },
];

pub const REF_HASHES: &[u64] = &[
    0x36bf16afef49c7ad,
    0x36bf16afef49c7ad,
    0x36bf16afef49c7ad,
    0x7ed13391c638dc2e,
    0x61ac7ee1bb73d36d,
    0xd62b1c97fb27608f,
    0x2dd4b2985c1bebd2,
    0x7ed13391c638dc2e,
    0x7ed13391c638dc2e,
    0x9a20976998759371,
    0xa67de9d0530d5968,
    0xa67de9d0530d5968,
    0x3d593e440ca8eb97,
    0x487a213a131d4958,
    0x201a627d0d128fd5,
    0xe787dcb23cc5b120,
    0x5aa4baf48a888be9,
    0xa2519e06e3bfa2fd,
    0xae983ea5eb1cd6f4,
    0x67b3a14cdd541dc3,
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inventory_is_complete_and_paired() {
        assert_eq!(TEST_CASES.len(), 20);
        assert_eq!(REF_HASHES.len(), TEST_CASES.len());
        assert_eq!(TEST_CASES[0].command_line, "blastp -p1");
        assert_eq!(REF_HASHES[0], 0x36bf16afef49c7ad);
        assert_eq!(TEST_CASES[19].description, "blastp (PAF format)");
        assert_eq!(REF_HASHES[19], 0x67b3a14cdd541dc3);
    }

    #[test]
    fn source_and_rust_have_the_same_literal_inventory() {
        let source =
            String::from_utf8_lossy(include_bytes!("../../diamond/src/test/test_cases.cpp"));
        for case in TEST_CASES {
            assert!(source.contains(&format!("\"{}\"", case.description)));
            assert!(source.contains(&format!("\"{}\"", case.command_line)));
        }
        for hash in REF_HASHES {
            assert!(source.contains(&format!("0x{hash:016x}")));
        }
    }
}
