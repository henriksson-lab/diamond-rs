//! NCBI BLAST nucleotide scoring/statistics support from
//! `diamond/src/lib/blast/blastn_score.cpp`.

use super::blast_message::{
    blast_message_write, BlastMessage, BlastSeverity, BLAST_MESSAGE_NO_CONTEXT,
};
use crate::stats::cbs::{self, BlastScoreFreq as BorrowedScoreFreq};

pub const BLAST_SCORE_MIN: i32 = i16::MIN as i32;
pub const BLAST_SCORE_MAX: i32 = i16::MAX as i32;
pub const BLASTNA_SIZE: usize = 16;
pub const BLASTAA_SIZE: usize = 28;
pub const BLASTNA_SEQ_CODE: u8 = 99;
pub const BLASTAA_SEQ_CODE: u8 = 11;
pub const NCBI4NA_SEQ_CODE: u8 = 4;
pub const BLASTERR_INVALIDPARAM: i16 = 75;

pub const NUCLEOTIDE_QUERY_MASK: u32 = 1 << 2;
pub const NUCLEOTIDE_SUBJECT_MASK: u32 = 1 << 3;
pub const BLAST_TYPE_BLASTN: u32 = NUCLEOTIDE_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK;
pub const BLAST_TYPE_TBLASTX: u32 = NUCLEOTIDE_QUERY_MASK | NUCLEOTIDE_SUBJECT_MASK | (1 << 5);

pub const NCBI4NA_TO_BLASTNA: [u8; BLASTNA_SIZE] =
    [15, 0, 1, 6, 2, 4, 9, 13, 3, 8, 5, 12, 7, 11, 10, 14];
pub const BLASTNA_TO_NCBI4NA: [u8; BLASTNA_SIZE] =
    [1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0];
pub const BLASTNA_TO_IUPACNA: &[u8; BLASTNA_SIZE] = b"ACGTRYMKWSBDHVN-";
pub const NCBI4NA_TO_IUPACNA: &[u8; BLASTNA_SIZE] = b"-ACMGRSVTWYHKDBN";
const NCBISTDAA: &[u8] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
const IDENTITY_SYMBOLS: &[u8] = b"ARNDCQEGHILKMFPSTWYVBJZX*";

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BlastKarlinBlock {
    pub lambda: f64,
    pub k: f64,
    pub log_k: f64,
    pub h: f64,
    pub param_c: f64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BlastGumbelBlock {
    pub lambda: f64,
    pub c: f64,
    pub g: f64,
    pub a: f64,
    pub alpha: f64,
    pub sigma: f64,
    pub a_un: f64,
    pub alpha_un: f64,
    pub b: f64,
    pub beta: f64,
    pub tau: f64,
    pub db_length: i64,
    pub filled: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BlastScoreMatrix {
    pub data: Vec<Vec<i32>>,
    pub ncols: usize,
    pub nrows: usize,
    pub freqs: Vec<f64>,
    pub lambda: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BlastScoreFreq {
    pub score_min: i32,
    pub score_max: i32,
    pub obs_min: i32,
    pub obs_max: i32,
    pub score_avg: f64,
    pub probabilities: Vec<f64>,
}

impl BlastScoreFreq {
    fn get(&self, score: i32) -> f64 {
        self.probabilities[(score - self.score_min) as usize]
    }

    fn set(&mut self, score: i32, value: f64) {
        self.probabilities[(score - self.score_min) as usize] = value;
    }

    fn borrowed(&self) -> BorrowedScoreFreq<'_> {
        BorrowedScoreFreq {
            score_min: self.score_min,
            score_max: self.score_max,
            obs_min: self.obs_min,
            obs_max: self.obs_max,
            score_avg: self.score_avg,
            score_probs: &self.probabilities,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BlastResFreq {
    pub alphabet_code: u8,
    pub probabilities: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BlastScoreBlock {
    pub protein_alphabet: bool,
    pub alphabet_code: u8,
    pub alphabet_size: usize,
    pub alphabet_start: usize,
    pub name: Option<String>,
    pub comments: Vec<String>,
    pub matrix: BlastScoreMatrix,
    pub matrix_only_scoring: bool,
    pub complexity_adjusted_scoring: bool,
    pub loscore: i32,
    pub hiscore: i32,
    pub penalty: i32,
    pub reward: i32,
    pub scale_factor: f64,
    pub read_in_matrix: bool,
    pub score_freqs: Vec<Option<BlastScoreFreq>>,
    pub karlin_std: Vec<Option<BlastKarlinBlock>>,
    pub karlin_gap_std: Vec<Option<BlastKarlinBlock>>,
    pub karlin_psi: Vec<Option<BlastKarlinBlock>>,
    pub karlin_gap_psi: Vec<Option<BlastKarlinBlock>>,
    pub karlin_ideal: Option<BlastKarlinBlock>,
    pub gumbel: Option<BlastGumbelBlock>,
    pub ambiguous_residues: Vec<u8>,
    pub round_down: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlastScoringOptions {
    pub matrix: Option<String>,
    pub matrix_path: Option<String>,
    pub reward: i16,
    pub penalty: i16,
    pub gapped_calculation: bool,
    pub complexity_adjusted_scoring: bool,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub is_ooframe: bool,
    pub shift_pen: i32,
    pub program_number: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedScoreMatrix {
    pub symbols: &'static [u8],
    pub scores: Vec<i32>,
    pub default_score: i32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PsiBlastScoreMatrix {
    pub pssm: BlastScoreMatrix,
    pub frequency_ratios: Vec<Vec<f64>>,
    pub karlin: BlastKarlinBlock,
}

pub fn psi_deallocate_matrix<T>(
    _matrix: Option<Vec<Vec<T>>>,
    _ncols: usize,
) -> Option<Vec<Vec<T>>> {
    None
}

pub fn psi_allocate_matrix<T: Default + Clone>(
    ncols: usize,
    nrows: usize,
    _data_type_size: usize,
) -> Option<Vec<Vec<T>>> {
    let mut matrix = Vec::new();
    matrix.try_reserve_exact(ncols).ok()?;
    for _ in 0..ncols {
        let mut row = Vec::new();
        row.try_reserve_exact(nrows).ok()?;
        row.resize(nrows, T::default());
        matrix.push(row);
    }
    Some(matrix)
}

pub fn blast_karlin_e_to_s_simple(e: f64, kbp: &BlastKarlinBlock, search_space: i64) -> i32 {
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return BLAST_SCORE_MIN;
    }
    ((kbp.k * search_space as f64 / e.max(1.0e-297)).ln() / kbp.lambda).ceil() as i32
}

pub fn blast_karlin_s_to_e_simple(score: i32, kbp: &BlastKarlinBlock, search_space: i64) -> f64 {
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return -1.0;
    }
    search_space as f64 * (-kbp.lambda * score as f64 + kbp.log_k).exp()
}

pub fn s_blast_gumbel_blk_new() -> Option<BlastGumbelBlock> {
    Some(BlastGumbelBlock::default())
}

pub fn s_blast_gumbel_blk_free(_block: Option<BlastGumbelBlock>) -> Option<BlastGumbelBlock> {
    None
}

pub fn s_blast_score_matrix_new(ncols: usize, nrows: usize) -> Option<BlastScoreMatrix> {
    Some(BlastScoreMatrix {
        data: psi_allocate_matrix(ncols, nrows, std::mem::size_of::<i32>())?,
        ncols,
        nrows,
        freqs: vec![0.0; ncols],
        lambda: 0.0,
    })
}

pub fn s_blast_score_matrix_free(mut matrix: Option<BlastScoreMatrix>) -> Option<BlastScoreMatrix> {
    if let Some(matrix) = matrix.as_mut() {
        for row in &mut matrix.data {
            row.clear();
        }
        matrix.data.clear();
        matrix.freqs.clear();
    }
    drop(matrix);
    None
}

pub fn blast_scoring_options_free(
    _options: Option<BlastScoringOptions>,
) -> Option<BlastScoringOptions> {
    None
}

pub fn blast_gcd(mut a: i32, mut b: i32) -> i32 {
    b = b.abs();
    if b > a {
        std::mem::swap(&mut a, &mut b);
    }
    while b != 0 {
        let c = a % b;
        a = b;
        b = c;
    }
    a
}

pub fn blast_nint(x: f64) -> i64 {
    (x + if x >= 0.0 { 0.5 } else { -0.5 }) as i64
}

pub fn blast_karlin_blk_copy(
    destination: Option<&mut BlastKarlinBlock>,
    source: Option<&BlastKarlinBlock>,
) -> i16 {
    match (destination, source) {
        (Some(destination), Some(source)) => {
            *destination = *source;
            0
        }
        _ => -1,
    }
}

pub fn s_ncbism_starts_with(value: &str, prefix: &str) -> i32 {
    i32::from(
        value.len() >= prefix.len()
            && value
                .bytes()
                .zip(prefix.bytes())
                .all(|(value, prefix)| value.to_ascii_lowercase() == prefix),
    )
}

pub fn ncbism_get_index(matrix: &PackedScoreMatrix, mut amino_acid: i32) -> i32 {
    if amino_acid >= 0 && (amino_acid as usize) < NCBISTDAA.len() + 1 {
        amino_acid = NCBISTDAA.get(amino_acid as usize).copied().unwrap_or(0) as i32;
    } else if (0..=u8::MAX as i32).contains(&amino_acid) {
        amino_acid = (amino_acid as u8).to_ascii_uppercase() as i32;
    }
    matrix
        .symbols
        .iter()
        .position(|&symbol| symbol as i32 == amino_acid)
        .map_or(-1, |index| index as i32)
}

pub fn ncbism_get_score(matrix: &PackedScoreMatrix, aa1: i32, aa2: i32) -> i32 {
    let i1 = ncbism_get_index(matrix, aa1);
    let i2 = ncbism_get_index(matrix, aa2);
    if i1 >= 0 && i2 >= 0 {
        matrix.scores[i1 as usize * matrix.symbols.len() + i2 as usize]
    } else {
        matrix.default_score
    }
}

pub fn ncbism_get_standard_matrix(name: &str) -> Option<PackedScoreMatrix> {
    // Preserve the vendored switch fall-through: BLOSUM falls into PAM and
    // PAM falls into IDENTITY, so only IDENTITY reaches the return value.
    if !name.to_ascii_lowercase().starts_with("identity") {
        return None;
    }
    let width = IDENTITY_SYMBOLS.len();
    let mut scores = vec![-5; width * width];
    for index in 0..width {
        scores[index * width + index] = 9;
    }
    Some(PackedScoreMatrix {
        symbols: IDENTITY_SYMBOLS,
        scores,
        default_score: -5,
    })
}

pub fn blast_score_blk_new(alphabet: u8, contexts: i32) -> Option<BlastScoreBlock> {
    let alphabet_size = if alphabet == BLASTNA_SEQ_CODE {
        BLASTNA_SIZE
    } else {
        BLASTAA_SIZE
    };
    let context_count = usize::try_from(contexts).ok()?;
    Some(BlastScoreBlock {
        protein_alphabet: alphabet == BLASTAA_SEQ_CODE,
        alphabet_code: alphabet,
        alphabet_size,
        alphabet_start: 0,
        name: None,
        comments: Vec::new(),
        matrix: s_blast_score_matrix_new(alphabet_size, alphabet_size)?,
        matrix_only_scoring: false,
        complexity_adjusted_scoring: false,
        loscore: 0,
        hiscore: 0,
        penalty: 0,
        reward: 0,
        scale_factor: 1.0,
        read_in_matrix: false,
        score_freqs: vec![None; context_count],
        karlin_std: vec![None; context_count],
        karlin_gap_std: vec![None; context_count],
        karlin_psi: vec![None; context_count],
        karlin_gap_psi: vec![None; context_count],
        karlin_ideal: None,
        gumbel: std::env::var_os("OLD_FSC")
            .is_none()
            .then(BlastGumbelBlock::default),
        ambiguous_residues: Vec::new(),
        round_down: false,
    })
}

pub fn blast_score_blk_free(mut block: Option<BlastScoreBlock>) -> Option<BlastScoreBlock> {
    if let Some(block) = block.as_mut() {
        for context in &mut block.score_freqs {
            let _ = blast_score_freq_free(context.take());
        }
        for contexts in [
            &mut block.karlin_std,
            &mut block.karlin_gap_std,
            &mut block.karlin_psi,
            &mut block.karlin_gap_psi,
        ] {
            for context in contexts {
                let _ = blast_karlin_blk_free(context.take());
            }
        }
        let _ = blast_karlin_blk_free(block.karlin_ideal.take());
        let _ = s_blast_gumbel_blk_free(block.gumbel.take());
        block.ambiguous_residues.clear();
        block.comments.clear();
    }
    drop(block);
    None
}

pub fn blast_karlin_blk_new() -> Option<BlastKarlinBlock> {
    Some(BlastKarlinBlock::default())
}

pub fn blast_karlin_blk_free(_block: Option<BlastKarlinBlock>) -> Option<BlastKarlinBlock> {
    None
}

pub fn blast_score_freq_new(score_min: i32, score_max: i32) -> Option<BlastScoreFreq> {
    if blast_score_chk(score_min, score_max) != 0 {
        return None;
    }
    Some(BlastScoreFreq {
        score_min,
        score_max,
        obs_min: 0,
        obs_max: 0,
        score_avg: 0.0,
        probabilities: vec![0.0; (score_max - score_min + 1) as usize],
    })
}

pub fn blast_score_freq_free(_freq: Option<BlastScoreFreq>) -> Option<BlastScoreFreq> {
    None
}

pub fn s_psi_blast_score_matrix_free(
    _matrix: Option<PsiBlastScoreMatrix>,
) -> Option<PsiBlastScoreMatrix> {
    None
}

pub fn blast_res_freq_new(block: Option<&BlastScoreBlock>) -> Option<BlastResFreq> {
    let block = block?;
    Some(BlastResFreq {
        alphabet_code: block.alphabet_code,
        probabilities: vec![0.0; block.alphabet_size],
    })
}

pub fn blast_res_freq_free(_freq: Option<BlastResFreq>) -> Option<BlastResFreq> {
    None
}

fn iupac_to_blastna(character: u8) -> usize {
    BLASTNA_TO_IUPACNA
        .iter()
        .position(|&value| value == character.to_ascii_uppercase())
        .unwrap_or(15)
}

fn amino_to_ncbistdaa(character: u8) -> usize {
    NCBISTDAA
        .iter()
        .position(|&value| value == character.to_ascii_uppercase())
        .unwrap_or(0)
}

pub fn blast_score_blk_protein_matrix_read(block: &mut BlastScoreBlock, text: &str) -> i16 {
    if block.alphabet_size != BLASTAA_SIZE {
        return 2;
    }
    for row in &mut block.matrix.data {
        row.fill(BLAST_SCORE_MIN);
    }
    let mut lines = text.lines();
    let columns = loop {
        let Some(line) = lines.next() else { return 2 };
        let trimmed = line.trim();
        if let Some(comment) = trimmed.strip_prefix('#') {
            block.comments.push(comment.to_owned());
            continue;
        }
        let columns = trimmed
            .split('#')
            .next()
            .unwrap_or("")
            .split_whitespace()
            .map(|token| amino_to_ncbistdaa(token.as_bytes()[0]))
            .collect::<Vec<_>>();
        if !columns.is_empty() {
            break columns;
        }
    };
    if columns.len() <= 1 || columns.len() > BLASTAA_SIZE {
        return 2;
    }
    let mut rows = 0usize;
    for line in lines {
        let data = line.split('#').next().unwrap_or("").trim();
        if data.is_empty() {
            continue;
        }
        let mut fields = data.split_whitespace();
        let Some(symbol) = fields.next() else {
            continue;
        };
        let row = amino_to_ncbistdaa(symbol.as_bytes()[0]);
        let scores = fields.collect::<Vec<_>>();
        if scores.len() != columns.len() || rows >= BLASTAA_SIZE {
            return 2;
        }
        for (&column, token) in columns.iter().zip(scores) {
            let score = if token.eq_ignore_ascii_case("na") {
                BLAST_SCORE_MIN
            } else {
                let Ok(value) = token.parse::<f64>() else {
                    return 2;
                };
                if !(BLAST_SCORE_MIN as f64..=BLAST_SCORE_MAX as f64).contains(&value) {
                    return 2;
                }
                blast_nint(value) as i32
            };
            block.matrix.data[row][column] = score;
        }
        rows += 1;
    }
    if rows <= 1 {
        return 2;
    }
    let (x, u, o, c) = (
        amino_to_ncbistdaa(b'X'),
        amino_to_ncbistdaa(b'U'),
        amino_to_ncbistdaa(b'O'),
        amino_to_ncbistdaa(b'C'),
    );
    for index in 0..block.alphabet_size {
        block.matrix.data[u][index] = block.matrix.data[c][index];
        block.matrix.data[index][u] = block.matrix.data[index][c];
        block.matrix.data[o][index] = block.matrix.data[x][index];
        block.matrix.data[index][o] = block.matrix.data[index][x];
    }
    0
}

pub fn blast_score_blk_max_score_set(block: &mut BlastScoreBlock) -> i16 {
    block.loscore = BLAST_SCORE_MAX;
    block.hiscore = BLAST_SCORE_MIN;
    for &score in block.matrix.data.iter().flatten() {
        if score <= BLAST_SCORE_MIN || score >= BLAST_SCORE_MAX {
            continue;
        }
        block.loscore = block.loscore.min(score);
        block.hiscore = block.hiscore.max(score);
    }
    block.loscore = block.loscore.max(BLAST_SCORE_MIN);
    block.hiscore = block.hiscore.min(BLAST_SCORE_MAX);
    0
}

pub fn blast_score_blk_nucleotide_matrix_read(block: &mut BlastScoreBlock, text: &str) -> i16 {
    if block.alphabet_size != BLASTNA_SIZE {
        return 2;
    }
    for row in &mut block.matrix.data {
        row.fill(BLAST_SCORE_MIN);
    }
    block.matrix.freqs.fill(0.0);
    let mut alphabet = Vec::new();
    let mut matrix_rows = Vec::<Vec<i32>>::new();
    let mut frequency_count = 0;
    for raw in text.lines() {
        let line = raw.trim();
        if line.starts_with('#') {
            if let Some(rest) = line.find("FREQS").map(|index| &line[index + 5..]) {
                let fields = rest.split_whitespace().collect::<Vec<_>>();
                if fields.len() % 2 != 0 {
                    return 2;
                }
                for pair in fields.chunks_exact(2) {
                    let Ok(value) = pair[1].parse::<f64>() else {
                        return 2;
                    };
                    block.matrix.freqs[iupac_to_blastna(pair[0].as_bytes()[0])] = value;
                    frequency_count += 1;
                }
            } else {
                block.comments.push(line.to_owned());
            }
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if alphabet.is_empty()
            && fields
                .iter()
                .all(|field| field.as_bytes()[0].is_ascii_alphabetic())
        {
            alphabet = fields
                .iter()
                .map(|field| field.as_bytes()[0].to_ascii_uppercase())
                .collect();
            continue;
        }
        let score_fields = if fields
            .first()
            .is_some_and(|field| field.len() == 1 && field.as_bytes()[0].is_ascii_alphabetic())
        {
            &fields[1..]
        } else {
            &fields[..]
        };
        let scores = score_fields
            .iter()
            .map(|value| value.parse::<i32>())
            .collect::<Result<Vec<_>, _>>();
        let Ok(scores) = scores else { return 2 };
        if scores.len() != alphabet.len() {
            return 2;
        }
        matrix_rows.push(scores);
    }
    if frequency_count != 4 || matrix_rows.len() != alphabet.len() {
        return 2;
    }
    for (row_symbol, scores) in alphabet.iter().zip(&matrix_rows) {
        let row = iupac_to_blastna(*row_symbol);
        for (&column_symbol, &score) in alphabet.iter().zip(scores) {
            block.matrix.data[row][iupac_to_blastna(column_symbol)] = score;
        }
    }
    let mut lower = 0.0;
    let mut lambda = 0.5;
    loop {
        let sum = matrix_lambda_sum(block, lambda);
        if sum >= 1.0 {
            break;
        }
        lower = lambda;
        lambda *= 2.0;
    }
    let mut upper = lambda;
    while upper - lower > 0.00001 {
        lambda = (lower + upper) / 2.0;
        if matrix_lambda_sum(block, lambda) >= 1.0 {
            upper = lambda;
        } else {
            lower = lambda;
        }
    }
    block.matrix.lambda = lambda;
    for index in 0..BLASTNA_SIZE {
        block.matrix.data[BLASTNA_SIZE - 1][index] = i32::MIN / 2;
        block.matrix.data[index][BLASTNA_SIZE - 1] = i32::MIN / 2;
    }
    0
}

fn matrix_lambda_sum(block: &BlastScoreBlock, lambda: f64) -> f64 {
    let mut sum = 0.0;
    for row in 0..block.alphabet_size {
        for column in 0..block.alphabet_size {
            if block.matrix.freqs[row] != 0.0 && block.matrix.freqs[column] != 0.0 {
                sum += block.matrix.freqs[row]
                    * block.matrix.freqs[column]
                    * (lambda * block.matrix.data[row][column] as f64).exp();
            }
        }
    }
    sum
}

pub fn blast_score_blk_nucl_matrix_create(block: &mut BlastScoreBlock) -> i16 {
    if block.alphabet_size != BLASTNA_SIZE {
        return 1;
    }
    let mut degeneracy = [0_i32; BLASTNA_SIZE];
    degeneracy[..4].fill(1);
    for index in 4..BLASTNA_SIZE {
        degeneracy[index] = (0..4)
            .filter(|&base| BLASTNA_TO_NCBI4NA[index] & BLASTNA_TO_NCBI4NA[base] != 0)
            .count() as i32;
    }
    for first in 0..BLASTNA_SIZE {
        for second in first..BLASTNA_SIZE {
            let score = if BLASTNA_TO_NCBI4NA[first] & BLASTNA_TO_NCBI4NA[second] != 0 {
                blast_nint(
                    (((degeneracy[second] - 1) * block.penalty + block.reward) as f64)
                        / degeneracy[second] as f64,
                ) as i32
            } else {
                block.penalty
            };
            block.matrix.data[first][second] = score;
            block.matrix.data[second][first] = score;
        }
    }
    for index in 0..BLASTNA_SIZE {
        block.matrix.data[BLASTNA_SIZE - 1][index] = i32::MIN / 2;
        block.matrix.data[index][BLASTNA_SIZE - 1] = i32::MIN / 2;
    }
    0
}

pub fn blast_score_blk_protein_matrix_load(block: &mut BlastScoreBlock) -> i16 {
    let Some(packed) = block.name.as_deref().and_then(ncbism_get_standard_matrix) else {
        return 1;
    };
    for row in &mut block.matrix.data {
        row.fill(BLAST_SCORE_MIN);
    }
    let (x, u, o, c, gap) = (
        amino_to_ncbistdaa(b'X'),
        amino_to_ncbistdaa(b'U'),
        amino_to_ncbistdaa(b'O'),
        amino_to_ncbistdaa(b'C'),
        amino_to_ncbistdaa(b'-'),
    );
    for row in 0..block.alphabet_size {
        for column in 0..block.alphabet_size {
            if [u, o, gap].contains(&row) || [u, o, gap].contains(&column) {
                continue;
            }
            block.matrix.data[row][column] = ncbism_get_score(&packed, row as i32, column as i32);
        }
    }
    for index in 0..block.alphabet_size {
        block.matrix.data[u][index] = block.matrix.data[c][index];
        block.matrix.data[index][u] = block.matrix.data[index][c];
        block.matrix.data[o][index] = block.matrix.data[x][index];
        block.matrix.data[index][o] = block.matrix.data[index][x];
    }
    0
}

pub fn blast_score_blk_matrix_fill(block: &mut BlastScoreBlock, has_path_callback: bool) -> i16 {
    let status = if block.alphabet_code == BLASTNA_SEQ_CODE {
        if block.read_in_matrix && has_path_callback {
            -1
        } else {
            blast_score_blk_nucl_matrix_create(block)
        }
    } else {
        blast_score_blk_protein_matrix_load(block)
    };
    if status != 0 {
        return status;
    }
    blast_score_blk_max_score_set(block)
}

pub fn blast_score_freq_calc(
    block: Option<&BlastScoreBlock>,
    score_freq: Option<&mut BlastScoreFreq>,
    first: &BlastResFreq,
    second: &BlastResFreq,
) -> i16 {
    let (Some(block), Some(score_freq)) = (block, score_freq) else {
        return 1;
    };
    if block.loscore < score_freq.score_min || block.hiscore > score_freq.score_max {
        return 1;
    }
    score_freq.probabilities.fill(0.0);
    for row in block.alphabet_start..block.alphabet_start + block.alphabet_size {
        for column in block.alphabet_start..block.alphabet_start + block.alphabet_size {
            let score = block.matrix.data[row][column];
            if score >= block.loscore {
                let value =
                    score_freq.get(score) + first.probabilities[row] * second.probabilities[column];
                score_freq.set(score, value);
            }
        }
    }
    let mut sum = 0.0;
    let mut observed_min = BLAST_SCORE_MIN;
    let mut observed_max = BLAST_SCORE_MIN;
    for score in score_freq.score_min..=score_freq.score_max {
        if score_freq.get(score) > 0.0 {
            sum += score_freq.get(score);
            observed_max = score;
            if observed_min == BLAST_SCORE_MIN {
                observed_min = score;
            }
        }
    }
    score_freq.obs_min = observed_min;
    score_freq.obs_max = observed_max;
    score_freq.score_avg = 0.0;
    if sum.abs() > 0.0001 {
        for score in observed_min..=observed_max {
            let probability = score_freq.get(score) / sum;
            score_freq.set(score, probability);
            score_freq.score_avg += score as f64 * probability;
        }
    }
    0
}

pub fn blast_score_blk_kbp_ideal_calc(block: Option<&mut BlastScoreBlock>) -> i16 {
    let Some(block) = block else { return 1 };
    let Some(mut residue_frequency) = blast_res_freq_new(Some(block)) else {
        return 1;
    };
    blast_res_freq_std_comp(block, &mut residue_frequency);
    let Some(mut score_frequency) = blast_score_freq_new(block.loscore, block.hiscore) else {
        return 1;
    };
    if blast_score_freq_calc(
        Some(block),
        Some(&mut score_frequency),
        &residue_frequency,
        &residue_frequency,
    ) != 0
    {
        return 1;
    }
    let mut ideal = BlastKarlinBlock::default();
    let _ = blast_karlin_blk_ungapped_calc(Some(&mut ideal), Some(&score_frequency));
    block.karlin_ideal = Some(ideal);
    0
}

pub fn blast_score_chk(low: i32, high: i32) -> i16 {
    if low >= 0
        || high <= 0
        || low < BLAST_SCORE_MIN
        || high > BLAST_SCORE_MAX
        || high - low > BLAST_SCORE_MAX - BLAST_SCORE_MIN
    {
        1
    } else {
        0
    }
}

pub fn blast_expm1(x: f64) -> f64 {
    let absolute = x.abs();
    if absolute > 0.33 {
        return x.exp() - 1.0;
    }
    if absolute < 1.0e-16 {
        return x;
    }
    x * (1.0
        + x * (1.0 / 2.0
            + x * (1.0 / 6.0
                + x * (1.0 / 24.0
                    + x * (1.0 / 120.0
                        + x * (1.0 / 720.0
                            + x * (1.0 / 5040.0
                                + x * (1.0 / 40320.0
                                    + x * (1.0 / 362880.0
                                        + x * (1.0 / 3628800.0
                                            + x * (1.0 / 39916800.0
                                                + x * (1.0 / 479001600.0
                                                    + x / 6227020800.0))))))))))))
}

pub fn blast_powi(mut value: f64, mut exponent: i32) -> f64 {
    if exponent == 0 {
        return 1.0;
    }
    if value == 0.0 {
        return if exponent < 0 { f64::INFINITY } else { 0.0 };
    }
    if exponent < 0 {
        value = 1.0 / value;
        exponent = -exponent;
    }
    let mut result = 1.0;
    while exponent > 0 {
        if exponent & 1 != 0 {
            result *= value;
        }
        exponent /= 2;
        value *= value;
    }
    result
}

pub fn blast_karlin_lambda_nr(score_freq: &BlastScoreFreq, initial_guess: f64) -> f64 {
    if blast_score_chk(score_freq.obs_min, score_freq.obs_max) != 0 {
        return -1.0;
    }
    cbs::blast_karlin_lambda_nr(&score_freq.borrowed(), initial_guess)
}

pub fn nlm_karlin_lambda_nr(
    probabilities: &[f64],
    probability_min: i32,
    low: i32,
    high: i32,
    lambda_zero: f64,
) -> f64 {
    let average = (low..=high)
        .map(|score| score as f64 * probabilities[(score - probability_min) as usize])
        .sum();
    cbs::blast_karlin_lambda_nr(
        &BorrowedScoreFreq {
            score_min: probability_min,
            score_max: probability_min + probabilities.len() as i32 - 1,
            obs_min: low,
            obs_max: high,
            score_avg: average,
            score_probs: probabilities,
        },
        lambda_zero,
    )
}

pub fn blast_karlin_l_to_h(score_freq: &BlastScoreFreq, lambda: f64) -> f64 {
    if lambda < 0.0 || blast_score_chk(score_freq.obs_min, score_freq.obs_max) != 0 {
        return -1.0;
    }
    let exponential = (-lambda).exp();
    let mut sum = score_freq.obs_min as f64 * score_freq.get(score_freq.obs_min);
    for score in score_freq.obs_min + 1..=score_freq.obs_max {
        sum = score as f64 * score_freq.get(score) + exponential * sum;
    }
    let scale = blast_powi(exponential, score_freq.obs_max);
    if scale > 0.0 {
        lambda * sum / scale
    } else {
        lambda * (lambda * score_freq.obs_max as f64 + sum.ln()).exp()
    }
}

pub fn blast_karlin_lh_to_k(score_freq: &BlastScoreFreq, mut lambda: f64, h: f64) -> f64 {
    if lambda <= 0.0 || h <= 0.0 || score_freq.score_avg >= 0.0 {
        return -1.0;
    }
    let mut low = score_freq.obs_min;
    let mut high = score_freq.obs_max;
    let mut divisor = -low;
    for offset in 1..=high - low {
        if score_freq.get(low + offset) != 0.0 && divisor > 1 {
            divisor = blast_gcd(divisor, offset);
        }
    }
    high /= divisor;
    low /= divisor;
    lambda *= divisor as f64;
    let mut first_term = h / lambda;
    let exp_minus_lambda = (-lambda).exp();
    if low == -1 && high == 1 {
        let p_low = score_freq.get(low * divisor);
        let p_high = score_freq.get(high * divisor);
        return (p_low - p_high).powi(2) / p_low;
    }
    if low == -1 || high == 1 {
        if high != 1 {
            let average = score_freq.score_avg / divisor as f64;
            first_term = average * average / first_term;
        }
        return first_term * (1.0 - exp_minus_lambda);
    }
    // General Karlin-Altschul convolution, matching the 100-iteration source.
    let range = high - low;
    let capacity = (100 * range + 1) as usize;
    let mut probabilities = vec![0.0; capacity];
    probabilities[0] = 1.0;
    let mut outer_sum = 0.0;
    let mut low_alignment = 0;
    let mut high_alignment = 0;
    let mut inner_sum = 1.0;
    let mut iteration = 0;
    while iteration < 100 && inner_sum > 0.0001 {
        low_alignment += low;
        high_alignment += high;
        let width = high_alignment - low_alignment;
        let previous = probabilities.clone();
        for score in low_alignment..=high_alignment {
            let mut value = 0.0;
            for step in low..=high {
                let prior_score = score - step;
                if prior_score >= low_alignment - low && prior_score <= high_alignment - high {
                    let prior_index = (prior_score - (low_alignment - low)) as usize;
                    value += previous[prior_index] * score_freq.get(step * divisor);
                }
            }
            probabilities[(score - low_alignment) as usize] = value;
        }
        inner_sum = 0.0;
        for score in low_alignment..0 {
            inner_sum +=
                probabilities[(score - low_alignment) as usize] * (lambda * score as f64).exp();
        }
        for score in 0..=high_alignment {
            inner_sum += probabilities[(score - low_alignment) as usize];
        }
        iteration += 1;
        outer_sum += inner_sum / iteration as f64;
        let _ = width;
    }
    -(-2.0 * outer_sum).exp() / (first_term * blast_expm1(-lambda))
}

pub fn blast_karlin_blk_ungapped_calc(
    block: Option<&mut BlastKarlinBlock>,
    score_freq: Option<&BlastScoreFreq>,
) -> i16 {
    let (Some(block), Some(score_freq)) = (block, score_freq) else {
        return 1;
    };
    block.lambda = blast_karlin_lambda_nr(score_freq, 0.5);
    if block.lambda >= 0.0 {
        block.h = blast_karlin_l_to_h(score_freq, block.lambda);
    }
    if block.lambda >= 0.0 && block.h >= 0.0 {
        block.k = blast_karlin_lh_to_k(score_freq, block.lambda, block.h);
    }
    if block.lambda < 0.0 || block.h < 0.0 || block.k < 0.0 {
        block.lambda = -1.0;
        block.h = -1.0;
        block.k = -1.0;
        block.log_k = f64::INFINITY;
        1
    } else {
        block.log_k = block.k.ln();
        0
    }
}

pub fn blast_res_freq_normalize(
    block: &BlastScoreBlock,
    frequency: &mut BlastResFreq,
    norm: f64,
) -> i16 {
    if norm == 0.0 {
        return 1;
    }
    let range = block.alphabet_start..block.alphabet_start + block.alphabet_size;
    let mut sum = 0.0;
    for index in range.clone() {
        if frequency.probabilities[index] < 0.0 {
            return 1;
        }
        sum += frequency.probabilities[index];
    }
    if sum <= 0.0 {
        return 0;
    }
    for index in range {
        frequency.probabilities[index] = frequency.probabilities[index] / sum * norm;
    }
    0
}

pub fn blast_res_freq_std_comp(block: &BlastScoreBlock, frequency: &mut BlastResFreq) -> i16 {
    frequency.probabilities[..4].fill(25.0);
    blast_res_freq_normalize(block, frequency, 1.0);
    0
}

pub const BLAST_NUM_STAT_VALUES: usize = 11;
pub type StatRow = [f64; BLAST_NUM_STAT_VALUES];

const V1_5: &[StatRow] = &[
    [0., 0., 1.39, 0.747, 1.38, 1., 0., 100., 0., 0., 0.],
    [3., 3., 1.39, 0.747, 1.38, 1., 0., 100., 0., 0., 0.],
];
const V1_4: &[StatRow] = &[
    [0., 0., 1.383, 0.738, 1.36, 1.02, 0., 100., 0., 0., 0.],
    [1., 2., 1.36, 0.67, 1.2, 1.1, 0., 98., 0., 0., 0.],
    [0., 2., 1.26, 0.43, 0.90, 1.4, -1., 91., 0., 0., 0.],
    [2., 1., 1.35, 0.61, 1.1, 1.2, -1., 98., 0., 0., 0.],
    [1., 1., 1.22, 0.35, 0.72, 1.7, -3., 88., 0., 0., 0.],
];
const V2_7: &[StatRow] = &[
    [0., 0., 0.69, 0.73, 1.34, 0.515, 0., 100., 0., 0., 0.],
    [2., 4., 0.68, 0.67, 1.2, 0.55, 0., 99., 0., 0., 0.],
    [0., 4., 0.63, 0.43, 0.90, 0.7, -1., 91., 0., 0., 0.],
    [4., 2., 0.675, 0.62, 1.1, 0.6, -1., 98., 0., 0., 0.],
    [2., 2., 0.61, 0.35, 0.72, 1.7, -3., 88., 0., 0., 0.],
];
const V1_3: &[StatRow] = &[
    [0., 0., 1.374, 0.711, 1.31, 1.05, 0., 100., 0., 0., 0.],
    [2., 2., 1.37, 0.70, 1.2, 1.1, 0., 99., 0., 0., 0.],
    [1., 2., 1.35, 0.64, 1.1, 1.2, -1., 98., 0., 0., 0.],
    [0., 2., 1.25, 0.42, 0.83, 1.5, -2., 91., 0., 0., 0.],
    [2., 1., 1.34, 0.60, 1.1, 1.2, -1., 97., 0., 0., 0.],
    [1., 1., 1.21, 0.34, 0.71, 1.7, -2., 88., 0., 0., 0.],
];
const V2_5: &[StatRow] = &[
    [0., 0., 0.675, 0.65, 1.1, 0.6, -1., 99., 0., 0., 0.],
    [2., 4., 0.67, 0.59, 1.1, 0.6, -1., 98., 0., 0., 0.],
    [0., 4., 0.62, 0.39, 0.78, 0.8, -2., 91., 0., 0., 0.],
    [4., 2., 0.67, 0.61, 1., 0.65, -2., 98., 0., 0., 0.],
    [2., 2., 0.56, 0.32, 0.59, 0.95, -4., 82., 0., 0., 0.],
];
const V1_2: &[StatRow] = &[
    [0., 0., 1.28, 0.46, 0.85, 1.5, -2., 96., 0., 0., 0.],
    [2., 2., 1.33, 0.62, 1.1, 1.2, 0., 99., 0., 0., 0.],
    [1., 2., 1.30, 0.52, 0.93, 1.4, -2., 97., 0., 0., 0.],
    [0., 2., 1.19, 0.34, 0.66, 1.8, -3., 89., 0., 0., 0.],
    [3., 1., 1.32, 0.57, 1., 1.3, -1., 99., 0., 0., 0.],
    [2., 1., 1.29, 0.49, 0.92, 1.4, -1., 96., 0., 0., 0.],
    [1., 1., 1.14, 0.26, 0.52, 2.2, -5., 85., 0., 0., 0.],
];
const V2_3: &[StatRow] = &[
    [0., 0., 0.55, 0.21, 0.46, 1.2, -5., 87., 0., 0., 0.],
    [4., 4., 0.63, 0.42, 0.84, 0.75, -2., 99., 0., 0., 0.],
    [2., 4., 0.615, 0.37, 0.72, 0.85, -3., 97., 0., 0., 0.],
    [0., 4., 0.55, 0.21, 0.46, 1.2, -5., 87., 0., 0., 0.],
    [3., 3., 0.615, 0.37, 0.68, 0.9, -3., 97., 0., 0., 0.],
    [6., 2., 0.63, 0.42, 0.84, 0.75, -2., 99., 0., 0., 0.],
    [5., 2., 0.625, 0.41, 0.78, 0.8, -2., 99., 0., 0., 0.],
    [4., 2., 0.61, 0.35, 0.68, 0.9, -3., 96., 0., 0., 0.],
    [2., 2., 0.515, 0.14, 0.33, 1.55, -9., 81., 0., 0., 0.],
];
const V3_4: &[StatRow] = &[
    [6., 3., 0.389, 0.25, 0.56, 0.7, -5., 95., 0., 0., 0.],
    [5., 3., 0.375, 0.21, 0.47, 0.8, -6., 92., 0., 0., 0.],
    [4., 3., 0.351, 0.14, 0.35, 1., -9., 86., 0., 0., 0.],
    [6., 2., 0.362, 0.16, 0.45, 0.8, -4., 88., 0., 0., 0.],
    [5., 2., 0.330, 0.092, 0.28, 1.2, -13., 81., 0., 0., 0.],
    [4., 2., 0.281, 0.046, 0.16, 1.8, -23., 69., 0., 0., 0.],
];
const V4_5: &[StatRow] = &[
    [0., 0., 0.22, 0.061, 0.22, 1., -15., 74., 0., 0., 0.],
    [6., 5., 0.28, 0.21, 0.47, 0.6, -7., 93., 0., 0., 0.],
    [5., 5., 0.27, 0.17, 0.39, 0.7, -9., 90., 0., 0., 0.],
    [4., 5., 0.25, 0.10, 0.31, 0.8, -10., 83., 0., 0., 0.],
    [3., 5., 0.23, 0.065, 0.25, 0.9, -11., 76., 0., 0., 0.],
];
const V1_1: &[StatRow] = &[
    [3., 2., 1.09, 0.31, 0.55, 2., -2., 99., 0., 0., 0.],
    [2., 2., 1.07, 0.27, 0.49, 2.2, -3., 97., 0., 0., 0.],
    [1., 2., 1.02, 0.21, 0.36, 2.8, -6., 92., 0., 0., 0.],
    [0., 2., 0.80, 0.064, 0.17, 4.8, -16., 72., 0., 0., 0.],
    [4., 1., 1.08, 0.28, 0.54, 2., -2., 98., 0., 0., 0.],
    [3., 1., 1.06, 0.25, 0.46, 2.3, -4., 96., 0., 0., 0.],
    [2., 1., 0.99, 0.17, 0.30, 3.3, -10., 90., 0., 0., 0.],
];
const V3_2: &[StatRow] = &[[5., 5., 0.208, 0.030, 0.072, 2.9, -47., 77., 0., 0., 0.]];
const V5_4: &[StatRow] = &[
    [10., 6., 0.163, 0.068, 0.16, 1., -19., 85., 0., 0., 0.],
    [8., 6., 0.146, 0.039, 0.11, 1.3, -29., 76., 0., 0., 0.],
];

pub fn s_adjust_gap_parameters_by_gcd(
    normal: &mut [StatRow],
    linear: Option<&mut StatRow>,
    gap_open_max: &mut i32,
    gap_extend_max: &mut i32,
    divisor: i32,
) -> i16 {
    if divisor == 1 {
        return 0;
    }
    if normal.is_empty() {
        return 1;
    }
    *gap_open_max *= divisor;
    *gap_extend_max *= divisor;
    for row in normal {
        row[0] *= divisor as f64;
        row[1] *= divisor as f64;
        row[2] /= divisor as f64;
        row[5] /= divisor as f64;
    }
    if let Some(row) = linear {
        row[0] *= divisor as f64;
        row[1] *= divisor as f64;
        row[2] /= divisor as f64;
        row[5] /= divisor as f64;
    }
    0
}

pub fn s_split_array_of_8(
    input: Option<&[StatRow]>,
) -> Result<(Vec<StatRow>, Option<StatRow>, bool), i16> {
    let input = input.ok_or(-1_i16)?;
    if input.is_empty() {
        return Err(-1);
    }
    if input[0][0] == 0.0 && input[0][1] == 0.0 {
        Ok((input[1..].to_vec(), Some(input[0]), true))
    } else {
        Ok((input.to_vec(), None, false))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct NuclValues {
    pub normal: Vec<StatRow>,
    pub non_affine: Option<StatRow>,
    pub gap_open_max: i32,
    pub gap_extend_max: i32,
    pub round_down: bool,
}

pub fn s_get_nucl_values_array(
    mut reward: i32,
    mut penalty: i32,
    errors: Option<&mut Option<Box<BlastMessage>>>,
) -> Result<NuclValues, i16> {
    let divisor = blast_gcd(reward, penalty);
    if divisor != 1 {
        reward /= divisor;
        penalty /= divisor;
    }
    let (table, open, extend, round_down) = match (reward, penalty) {
        (1, -5) => (V1_5, 3, 3, false),
        (1, -4) => (V1_4, 2, 2, false),
        (2, -7) => (V2_7, 4, 4, true),
        (1, -3) => (V1_3, 2, 2, false),
        (2, -5) => (V2_5, 4, 4, true),
        (1, -2) => (V1_2, 2, 2, false),
        (2, -3) => (V2_3, 6, 4, true),
        (3, -4) => (V3_4, 6, 3, true),
        (1, -1) => (V1_1, 4, 2, false),
        (3, -2) => (V3_2, 5, 5, false),
        (4, -5) => (V4_5, 12, 8, false),
        (5, -4) => (V5_4, 25, 10, false),
        _ => {
            if let Some(errors) = errors {
                let text = format!("Substitution scores {reward} and {penalty} are not supported");
                blast_message_write(
                    Some(errors),
                    BlastSeverity::Error,
                    BLAST_MESSAGE_NO_CONTEXT,
                    &text,
                );
            }
            return Err(-1);
        }
    };
    let (mut normal, mut non_affine, _) = s_split_array_of_8(Some(table))?;
    let mut gap_open_max = open;
    let mut gap_extend_max = extend;
    s_adjust_gap_parameters_by_gcd(
        &mut normal,
        non_affine.as_mut(),
        &mut gap_open_max,
        &mut gap_extend_max,
        divisor,
    );
    Ok(NuclValues {
        normal,
        non_affine,
        gap_open_max,
        gap_extend_max,
        round_down,
    })
}

pub fn blast_karlin_blk_nucl_gapped_calc(
    block: &mut BlastKarlinBlock,
    gap_open: i32,
    gap_extend: i32,
    reward: i32,
    penalty: i32,
    ungapped: &BlastKarlinBlock,
    round_down: &mut bool,
    errors: Option<&mut Option<Box<BlastMessage>>>,
) -> i16 {
    let mut errors = errors;
    let values = match s_get_nucl_values_array(reward, penalty, errors.as_deref_mut()) {
        Ok(values) => values,
        Err(status) => return status,
    };
    *round_down = values.round_down;
    if gap_open == 0 && gap_extend == 0 {
        if let Some(row) = values.non_affine {
            block.lambda = row[2];
            block.k = row[3];
            block.log_k = block.k.ln();
            block.h = row[4];
            return 0;
        }
    }
    if let Some(row) = values
        .normal
        .iter()
        .find(|row| row[0] == gap_open as f64 && row[1] == gap_extend as f64)
    {
        block.lambda = row[2];
        block.k = row[3];
        block.log_k = block.k.ln();
        block.h = row[4];
        return 0;
    }
    if gap_open >= values.gap_open_max && gap_extend >= values.gap_extend_max {
        *block = *ungapped;
        return 0;
    }
    if let Some(errors) = errors {
        let mut text = format!("Gap existence and extension values {gap_open} and {gap_extend} are not supported for substitution scores {reward} and {penalty}\n");
        for row in &values.normal {
            text.push_str(&format!(
                "{} and {} are supported existence and extension values\n",
                row[0] as i32, row[1] as i32
            ));
        }
        text.push_str(&format!("{} and {} are supported existence and extension values\nAny values more stringent than {} and {} are supported\n", values.gap_open_max, values.gap_extend_max, values.gap_open_max, values.gap_extend_max));
        blast_message_write(
            Some(errors),
            BlastSeverity::Error,
            BLAST_MESSAGE_NO_CONTEXT,
            &text,
        );
        return 1;
    }
    // The C implementation returns success without modifying kbp when no
    // error-return pointer was supplied.
    0
}

pub fn blast_str_to_upper(value: Option<&str>) -> Option<String> {
    value.map(|value| {
        value
            .bytes()
            .map(|byte| byte.to_ascii_uppercase() as char)
            .collect()
    })
}

pub const fn blast_subject_is_nucleotide(program: u32) -> bool {
    program & NUCLEOTIDE_SUBJECT_MASK != 0
}

pub const fn blast_query_is_nucleotide(program: u32) -> bool {
    program & NUCLEOTIDE_QUERY_MASK != 0
}

pub const fn blast_program_is_nucleotide(program: u32) -> bool {
    blast_query_is_nucleotide(program) && blast_subject_is_nucleotide(program)
}

pub fn blast_scoring_options_set_matrix(
    options: &mut BlastScoringOptions,
    matrix_name: Option<&str>,
) -> i16 {
    if let Some(matrix_name) = matrix_name {
        options.matrix = blast_str_to_upper(Some(matrix_name));
    }
    0
}

pub fn blast_scoring_options_new(program: u32) -> Result<BlastScoringOptions, i16> {
    let nucleotide = blast_program_is_nucleotide(program);
    Ok(BlastScoringOptions {
        matrix: (!nucleotide).then(|| "BLOSUM62".to_owned()),
        matrix_path: None,
        reward: if nucleotide { 1 } else { 0 },
        penalty: if nucleotide { -3 } else { 0 },
        gapped_calculation: program != BLAST_TYPE_TBLASTX,
        complexity_adjusted_scoring: false,
        gap_open: if nucleotide { 5 } else { 11 },
        gap_extend: if nucleotide { 2 } else { 1 },
        is_ooframe: false,
        shift_pen: if nucleotide { 0 } else { i16::MAX as i32 },
        program_number: program,
    })
}

#[allow(clippy::too_many_arguments)]
pub fn blast_fill_scoring_options(
    options: Option<&mut BlastScoringOptions>,
    program: u32,
    greedy_extension: bool,
    penalty: i32,
    reward: i32,
    matrix: Option<&str>,
    gap_open: i32,
    gap_extend: i32,
) -> i16 {
    let Some(options) = options else {
        return BLASTERR_INVALIDPARAM;
    };
    if !blast_program_is_nucleotide(program) {
        if matrix.is_some() {
            blast_scoring_options_set_matrix(options, matrix);
        }
    } else {
        if penalty != 0 {
            options.penalty = penalty as i16;
        }
        if reward != 0 {
            options.reward = reward as i16;
        }
        if greedy_extension {
            options.gap_open = 0;
            options.gap_extend = 0;
        } else {
            options.gap_open = 5;
            options.gap_extend = 2;
        }
    }
    if gap_open >= 0 {
        options.gap_open = gap_open;
    }
    if gap_extend >= 0 {
        options.gap_extend = gap_extend;
    }
    options.program_number = program;
    0
}

pub fn blast_score_set_ambig_res(block: Option<&mut BlastScoreBlock>, residue: char) -> i16 {
    let Some(block) = block else { return 1 };
    let character = residue as u8;
    let encoded = if block.alphabet_code == BLASTAA_SEQ_CODE {
        amino_to_ncbistdaa(character) as u8
    } else if block.alphabet_code == BLASTNA_SEQ_CODE {
        iupac_to_blastna(character) as u8
    } else if block.alphabet_code == NCBI4NA_SEQ_CODE {
        BLASTNA_TO_NCBI4NA[iupac_to_blastna(character)]
    } else {
        0
    };
    block.ambiguous_residues.push(encoded);
    0
}

pub fn blast_score_blk_matrix_init(
    program: u32,
    options: Option<&BlastScoringOptions>,
    block: Option<&mut BlastScoreBlock>,
    has_path_callback: bool,
) -> i16 {
    let (Some(options), Some(block)) = (options, block) else {
        return 1;
    };
    block.matrix_only_scoring = false;
    if program == BLAST_TYPE_BLASTN {
        blast_score_set_ambig_res(Some(block), 'N');
        blast_score_set_ambig_res(Some(block), '-');
        if options.penalty == 0 && options.reward == 0 {
            block.matrix_only_scoring = true;
            block.penalty = -3;
            block.reward = 1;
        } else {
            block.penalty = options.penalty as i32;
            block.reward = options.reward as i32;
        }
        if let Some(matrix) = options
            .matrix
            .as_deref()
            .filter(|matrix| !matrix.is_empty())
        {
            block.read_in_matrix = true;
            block.name = Some(matrix.to_owned());
        } else {
            block.read_in_matrix = false;
            block.name = Some(format!("blastn matrix:{} {}", block.reward, block.penalty));
        }
    } else {
        block.read_in_matrix = true;
        blast_score_set_ambig_res(Some(block), 'X');
        block.name = blast_str_to_upper(options.matrix.as_deref());
    }
    blast_score_blk_matrix_fill(block, has_path_callback)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected}, got {actual}"
        );
    }

    #[test]
    fn scalar_helpers_preserve_blast_rounding_and_statistics_conversions() {
        assert_eq!(blast_gcd(18, -12), 6);
        assert_eq!(blast_nint(2.5), 3);
        assert_eq!(blast_nint(-2.5), -3);
        for value in [-0.4_f64, -1.0e-17, 0.0, 0.2, 0.5] {
            close(blast_expm1(value), value.exp_m1(), 2.0e-14);
        }
        assert_eq!(blast_powi(2.0, 10), 1024.0);
        assert_eq!(blast_powi(2.0, -3), 0.125);
        assert!(blast_powi(0.0, -1).is_infinite());

        let kbp = BlastKarlinBlock {
            lambda: 0.5,
            k: 0.2,
            log_k: 0.2_f64.ln(),
            h: 1.0,
            param_c: 0.0,
        };
        let e = blast_karlin_s_to_e_simple(11, &kbp, 1_000_000);
        assert_eq!(blast_karlin_e_to_s_simple(e, &kbp, 1_000_000), 11);
        assert_eq!(blast_karlin_e_to_s_simple(-1.0, &kbp, 100), 1374);
        let invalid = BlastKarlinBlock {
            lambda: -1.0,
            ..kbp
        };
        assert_eq!(
            blast_karlin_e_to_s_simple(1.0, &invalid, 100),
            BLAST_SCORE_MIN
        );
    }

    #[test]
    fn nucleotide_matrix_creation_handles_ambiguous_bases_and_sentinel() {
        let mut block = blast_score_blk_new(BLASTNA_SEQ_CODE, 2).unwrap();
        block.reward = 1;
        block.penalty = -3;
        assert_eq!(blast_score_blk_nucl_matrix_create(&mut block), 0);
        assert_eq!(block.matrix.data[0][0], 1);
        assert_eq!(block.matrix.data[0][1], -3);
        // M = A/C: one match averaged with one mismatch.
        assert_eq!(block.matrix.data[0][6], -1);
        assert_eq!(block.matrix.data[6][0], -1);
        assert_eq!(block.matrix.data[15][0], i32::MIN / 2);
        assert_eq!(block.matrix.data[0][15], i32::MIN / 2);
        assert_eq!(blast_score_blk_max_score_set(&mut block), 0);
        assert_eq!((block.loscore, block.hiscore), (-3, 1));
        assert_eq!(block.score_freqs.len(), 2);
    }

    #[test]
    fn custom_nucleotide_matrix_reader_parses_frequencies_and_finds_lambda() {
        let text = "# FREQS A .25 C .25 G .25 T .25\n\
                    A C G T\n\
                    A 1 -3 -3 -3\n\
                    C -3 1 -3 -3\n\
                    G -3 -3 1 -3\n\
                    T -3 -3 -3 1\n";
        let mut block = blast_score_blk_new(BLASTNA_SEQ_CODE, 1).unwrap();
        assert_eq!(blast_score_blk_nucleotide_matrix_read(&mut block, text), 0);
        assert_eq!(block.matrix.freqs[..4], [0.25; 4]);
        assert_eq!(block.matrix.data[0][0], 1);
        assert_eq!(block.matrix.data[0][1], -3);
        close(matrix_lambda_sum(&block, block.matrix.lambda), 1.0, 2.0e-5);
        assert_eq!(block.matrix.data[15][3], i32::MIN / 2);
    }

    #[test]
    fn score_frequency_and_ungapped_karlin_values_match_binary_distribution() {
        let mut frequency = blast_score_freq_new(-1, 1).unwrap();
        frequency.set(-1, 0.75);
        frequency.set(1, 0.25);
        frequency.obs_min = -1;
        frequency.obs_max = 1;
        frequency.score_avg = -0.5;
        let mut kbp = BlastKarlinBlock::default();
        assert_eq!(
            blast_karlin_blk_ungapped_calc(Some(&mut kbp), Some(&frequency)),
            0
        );
        close(kbp.lambda, 3.0_f64.ln(), 1.0e-5);
        close(kbp.k, 1.0 / 3.0, 1.0e-5);
        assert!(kbp.h > 0.0);
        close(kbp.log_k, kbp.k.ln(), 1.0e-12);
    }

    #[test]
    fn nucleotide_gapped_tables_scale_and_report_unsupported_parameters() {
        let values = s_get_nucl_values_array(1, -3, None).unwrap();
        assert_eq!(values.normal.len(), 5);
        assert_eq!(values.non_affine.unwrap()[2], 1.374);
        assert_eq!((values.gap_open_max, values.gap_extend_max), (2, 2));
        assert!(!values.round_down);

        let scaled = s_get_nucl_values_array(2, -6, None).unwrap();
        assert_eq!((scaled.gap_open_max, scaled.gap_extend_max), (4, 4));
        close(scaled.normal[0][2], values.normal[0][2] / 2.0, 1.0e-12);

        let ungapped = BlastKarlinBlock {
            lambda: 2.0,
            k: 3.0,
            log_k: 4.0,
            h: 5.0,
            param_c: 6.0,
        };
        let mut result = BlastKarlinBlock::default();
        let mut round_down = true;
        assert_eq!(
            blast_karlin_blk_nucl_gapped_calc(
                &mut result,
                2,
                2,
                1,
                -3,
                &ungapped,
                &mut round_down,
                None
            ),
            0
        );
        close(result.lambda, 1.37, 1.0e-12);
        close(result.k, 0.70, 1.0e-12);
        close(result.h, 1.2, 1.0e-12);

        let mut errors = None;
        assert_eq!(
            blast_karlin_blk_nucl_gapped_calc(
                &mut result,
                0,
                1,
                1,
                -3,
                &ungapped,
                &mut round_down,
                Some(&mut errors)
            ),
            1
        );
        assert!(errors.unwrap().message.contains("are not supported"));
        assert_eq!(
            blast_karlin_blk_nucl_gapped_calc(
                &mut result,
                9,
                9,
                1,
                -3,
                &ungapped,
                &mut round_down,
                None
            ),
            0
        );
        assert_eq!(result, ungapped);
    }

    #[test]
    fn standard_matrix_and_scoring_options_preserve_active_vendor_behavior() {
        assert!(ncbism_get_standard_matrix("BLOSUM62").is_none());
        assert!(ncbism_get_standard_matrix("PAM30").is_none());
        let identity = ncbism_get_standard_matrix("iDeNtItY").unwrap();
        assert_eq!(ncbism_get_score(&identity, b'A' as i32, b'A' as i32), 9);
        assert_eq!(ncbism_get_score(&identity, b'A' as i32, b'R' as i32), -5);
        assert_eq!(s_ncbism_starts_with("id", "identity"), 0);

        let mut options = blast_scoring_options_new(BLAST_TYPE_BLASTN).unwrap();
        assert_eq!(
            (
                options.reward,
                options.penalty,
                options.gap_open,
                options.gap_extend
            ),
            (1, -3, 5, 2)
        );
        assert_eq!(
            blast_fill_scoring_options(
                Some(&mut options),
                BLAST_TYPE_BLASTN,
                true,
                -5,
                2,
                None,
                -1,
                -1
            ),
            0
        );
        assert_eq!(
            (
                options.reward,
                options.penalty,
                options.gap_open,
                options.gap_extend
            ),
            (2, -5, 0, 0)
        );
        assert_eq!(
            blast_fill_scoring_options(None, BLAST_TYPE_BLASTN, false, 0, 0, None, -1, -1),
            BLASTERR_INVALIDPARAM
        );
    }
}
