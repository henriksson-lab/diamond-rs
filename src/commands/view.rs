use std::fs::File;
use std::io::{self, Write};

use crate::data::daa::{
    finish_daa_from_file, init_daa, view_query_daa, view_query_edge, view_query_null,
    view_query_paf, view_query_pairwise, view_query_sam, view_query_tabular, view_query_xml,
    DaaFile, DaaQueryRecord, DaaViewConfig,
};
use crate::output::format::{print_tabular_footer, print_tabular_header, Header, TabularFormat};
use crate::output::{pairwise, sam, xml};
use crate::stats::score_matrix::ScoreMatrix;

#[derive(Debug, Clone)]
pub struct ViewConfig {
    pub daa_file: String,
    pub output: String,
    pub outfmt: Vec<String>,
    pub max_target_seqs: i64,
    pub toppercent: Option<f64>,
    pub forward_only: bool,
    pub report_unaligned: bool,
    pub sam_qlen_field: bool,
    pub invocation: String,
}

/// Run the view command — show DAA file metadata.
///
/// Non-DAA record rendering still uses FFI via the binary fallback.
pub fn run(daa_file: &str) -> io::Result<()> {
    let daa = DaaFile::open(daa_file)?;

    let mode_str = match daa.mode() {
        2 => "blastp",
        3 => "blastx",
        4 => "blastn",
        _ => "unknown",
    };

    eprintln!("DAA file: {}", daa_file);
    eprintln!("Mode: {}", mode_str);
    eprintln!("Database sequences: {}", daa.db_seqs());
    eprintln!("Database sequences used: {}", daa.db_seqs_used());
    eprintln!("Database letters: {}", daa.db_letters());
    eprintln!("Query records: {}", daa.query_records());
    eprintln!("Matrix: {}", daa.score_matrix());
    eprintln!(
        "Gap penalties: {}/{}",
        daa.gap_open_penalty(),
        daa.gap_extension_penalty()
    );

    if daa.db_seqs_used() > 0 {
        eprintln!(
            "First reference: {} ({})",
            daa.ref_name(0),
            daa.ref_len_at(0)
        );
    }

    Ok(())
}

/// Matches the DAA-output branch of C++ `view_daa`.
pub fn run_daa(config: &ViewConfig) -> io::Result<()> {
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
    };

    eprintln!(
        "Scoring parameters: matrix={} gap_open={} gap_extend={}",
        daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty()
    );
    eprintln!("Build version = {}", daa.diamond_build());
    eprintln!("DB sequences = {}", daa.db_seqs());
    eprintln!("DB sequences used = {}", daa.db_seqs_used());
    eprintln!("DB letters = {}", daa.db_letters());

    let mut out = File::create(&config.output)?;
    init_daa(&mut out)?;
    while let Some((buf, query_num)) = daa.read_query_buffer()? {
        let r = DaaQueryRecord::from_buffer(&daa, &buf, query_num)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut query_out = Vec::new();
        view_query_daa(&r, &daa, &mut query_out, &score_matrix, &view_config)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        std::io::Write::write_all(&mut out, &query_out)?;
    }
    finish_daa_from_file(&mut out, &daa)
}

/// Matches the tabular-output branch of C++ `view_daa`.
pub fn run_tabular(config: &ViewConfig) -> io::Result<()> {
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let outfmt = if config.outfmt.is_empty() {
        vec!["tab".to_string()]
    } else {
        config.outfmt.clone()
    };
    let is_json = outfmt
        .first()
        .is_some_and(|f| f == "json-flat" || f == "104");
    let outfmt_refs = outfmt.iter().map(String::as_str).collect::<Vec<_>>();
    let format = TabularFormat::new(
        &outfmt_refs,
        is_json,
        0,
        daa.mode() == 3,
        &daa.score_matrix(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
    };

    eprintln!(
        "Scoring parameters: matrix={} gap_open={} gap_extend={}",
        daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty()
    );
    eprintln!("Build version = {}", daa.diamond_build());
    eprintln!("DB sequences = {}", daa.db_seqs());
    eprintln!("DB sequences used = {}", daa.db_seqs_used());
    eprintln!("DB letters = {}", daa.db_letters());

    let mut out = File::create(&config.output)?;
    let header = print_tabular_header(&format.fields, Header::None, format.is_json, "", "")
        .map_err(io::Error::other)?;
    out.write_all(header.as_bytes())?;
    while let Some((buf, query_num)) = daa.read_query_buffer()? {
        let r = DaaQueryRecord::from_buffer(&daa, &buf, query_num)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut query_out = Vec::new();
        view_query_tabular(
            &r,
            &daa,
            &mut query_out,
            &score_matrix,
            &format,
            &view_config,
            config.report_unaligned,
        )
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        out.write_all(&query_out)?;
    }
    out.write_all(print_tabular_footer(format.is_json).as_bytes())?;
    Ok(())
}

/// Matches the PAF-output branch of C++ `view_daa`.
pub fn run_paf(config: &ViewConfig) -> io::Result<()> {
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
    };

    eprintln!(
        "Scoring parameters: matrix={} gap_open={} gap_extend={}",
        daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty()
    );
    eprintln!("Build version = {}", daa.diamond_build());
    eprintln!("DB sequences = {}", daa.db_seqs());
    eprintln!("DB sequences used = {}", daa.db_seqs_used());
    eprintln!("DB letters = {}", daa.db_letters());

    let mut out = File::create(&config.output)?;
    while let Some((buf, query_num)) = daa.read_query_buffer()? {
        let r = DaaQueryRecord::from_buffer(&daa, &buf, query_num)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut query_out = Vec::new();
        view_query_paf(&r, &daa, &mut query_out, &score_matrix, &view_config)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        out.write_all(&query_out)?;
    }
    Ok(())
}

/// Matches the SAM-output branch of C++ `view_daa`.
pub fn run_sam(config: &ViewConfig) -> io::Result<()> {
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
    };

    eprintln!(
        "Scoring parameters: matrix={} gap_open={} gap_extend={}",
        daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty()
    );
    eprintln!("Build version = {}", daa.diamond_build());
    eprintln!("DB sequences = {}", daa.db_seqs());
    eprintln!("DB sequences used = {}", daa.db_seqs_used());
    eprintln!("DB letters = {}", daa.db_letters());

    let mut out = File::create(&config.output)?;
    sam::print_header(&mut out, daa.mode(), &config.invocation)?;
    while let Some((buf, query_num)) = daa.read_query_buffer()? {
        let r = DaaQueryRecord::from_buffer(&daa, &buf, query_num)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut query_out = Vec::new();
        view_query_sam(
            &r,
            &daa,
            &mut query_out,
            &score_matrix,
            &view_config,
            config.sam_qlen_field,
        )
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        out.write_all(&query_out)?;
    }
    Ok(())
}

/// Matches the pairwise-output branch of C++ `view_daa`.
pub fn run_pairwise(config: &ViewConfig) -> io::Result<()> {
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
    };

    eprintln!(
        "Scoring parameters: matrix={} gap_open={} gap_extend={}",
        daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty()
    );
    eprintln!("Build version = {}", daa.diamond_build());
    eprintln!("DB sequences = {}", daa.db_seqs());
    eprintln!("DB sequences used = {}", daa.db_seqs_used());
    eprintln!("DB letters = {}", daa.db_letters());

    let mut out = File::create(&config.output)?;
    pairwise::write_header(&mut out)?;
    while let Some((buf, query_num)) = daa.read_query_buffer()? {
        let r = DaaQueryRecord::from_buffer(&daa, &buf, query_num)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut query_out = Vec::new();
        view_query_pairwise(&r, &daa, &mut query_out, &score_matrix, &view_config)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        out.write_all(&query_out)?;
    }
    pairwise::print_footer(&mut out)?;
    Ok(())
}

/// Matches the XML-output branch of C++ `view_daa`.
pub fn run_xml(config: &ViewConfig) -> io::Result<()> {
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
    };

    eprintln!(
        "Scoring parameters: matrix={} gap_open={} gap_extend={}",
        daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty()
    );
    eprintln!("Build version = {}", daa.diamond_build());
    eprintln!("DB sequences = {}", daa.db_seqs());
    eprintln!("DB sequences used = {}", daa.db_seqs_used());
    eprintln!("DB letters = {}", daa.db_letters());

    let mut out = File::create(&config.output)?;
    if let Some((buf, query_num)) = daa.read_query_buffer()? {
        let r = DaaQueryRecord::from_buffer(&daa, &buf, query_num)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        xml::print_header(
            &mut out,
            daa.mode(),
            &daa.score_matrix(),
            daa.gap_open_penalty(),
            daa.gap_extension_penalty(),
            daa.evalue(),
            &r.query_name,
            r.query_len() as u32,
            "",
        )?;
        let mut query_out = Vec::new();
        view_query_xml(&r, &daa, &mut query_out, &score_matrix, &view_config)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        out.write_all(&query_out)?;
        while let Some((buf, query_num)) = daa.read_query_buffer()? {
            let r = DaaQueryRecord::from_buffer(&daa, &buf, query_num)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let mut query_out = Vec::new();
            view_query_xml(&r, &daa, &mut query_out, &score_matrix, &view_config)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            out.write_all(&query_out)?;
        }
    } else {
        xml::print_header(
            &mut out,
            daa.mode(),
            &daa.score_matrix(),
            daa.gap_open_penalty(),
            daa.gap_extension_penalty(),
            daa.evalue(),
            "",
            0,
            "",
        )?;
    }
    xml::print_footer(&mut out)?;
    Ok(())
}

/// Matches the null-output branch of C++ `view_daa`.
pub fn run_null(config: &ViewConfig) -> io::Result<()> {
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
    };

    eprintln!(
        "Scoring parameters: matrix={} gap_open={} gap_extend={}",
        daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty()
    );
    eprintln!("Build version = {}", daa.diamond_build());
    eprintln!("DB sequences = {}", daa.db_seqs());
    eprintln!("DB sequences used = {}", daa.db_seqs_used());
    eprintln!("DB letters = {}", daa.db_letters());

    let mut out = File::create(&config.output)?;
    while let Some((buf, query_num)) = daa.read_query_buffer()? {
        let r = DaaQueryRecord::from_buffer(&daa, &buf, query_num)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        view_query_null(&r, &daa, &score_matrix, &view_config)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }
    out.flush()
}

/// Matches the edge-output branch of C++ `view_daa`.
pub fn run_edge(config: &ViewConfig) -> io::Result<()> {
    let mut daa = DaaFile::open(&config.daa_file)?;
    let score_matrix = ScoreMatrix::new(
        &daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty(),
        0,
        1,
        daa.db_letters(),
    )
    .map_err(io::Error::other)?;
    let view_config = DaaViewConfig {
        max_target_seqs: config.max_target_seqs,
        toppercent: config.toppercent,
        forward_only: config.forward_only,
    };

    eprintln!(
        "Scoring parameters: matrix={} gap_open={} gap_extend={}",
        daa.score_matrix(),
        daa.gap_open_penalty(),
        daa.gap_extension_penalty()
    );
    eprintln!("Build version = {}", daa.diamond_build());
    eprintln!("DB sequences = {}", daa.db_seqs());
    eprintln!("DB sequences used = {}", daa.db_seqs_used());
    eprintln!("DB letters = {}", daa.db_letters());

    let mut out = File::create(&config.output)?;
    while let Some((buf, query_num)) = daa.read_query_buffer()? {
        let r = DaaQueryRecord::from_buffer(&daa, &buf, query_num)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut query_out = Vec::new();
        view_query_edge(&r, &daa, &mut query_out, &score_matrix, &view_config)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        out.write_all(&query_out)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::packed_transcript::{EditOperation, PackedOperation};
    use crate::basic::value::SequenceType;
    use crate::data::daa::{
        compute_flag, write_daa_query_record, BlockType, DaaHeader1, DaaHeader2,
    };
    use std::io::Write;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static VIEW_TEST_FILE_ID: AtomicUsize = AtomicUsize::new(0);

    #[test]
    fn test_view_stub() {
        // View command shows DAA metadata — full record parsing uses FFI
    }

    #[test]
    fn test_run_daa_rewrites_archive() {
        let id = VIEW_TEST_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let input = std::env::temp_dir().join(format!(
            "diamond_rs_view_input_{}_{}.daa",
            std::process::id(),
            id
        ));
        let output = std::env::temp_dir().join(format!(
            "diamond_rs_view_output_{}_{}.daa",
            std::process::id(),
            id
        ));

        let mut query = Vec::new();
        let seek_pos =
            write_daa_query_record(&mut query, "query0", &[0, 1, 2], SequenceType::AminoAcid);
        crate::data::daa::finish_daa_query_record(&mut query, seek_pos);
        query.extend_from_slice(&0u32.to_ne_bytes());

        let mut h2 = DaaHeader2::with_params(1, 100, 11, 1, 2, -3, 0.7, 1.2, 0.001, "BLOSUM62", 2);
        h2.db_seqs_used = 1;
        h2.query_records = 1;
        h2.block_type[0] = BlockType::Alignments as u8;
        h2.block_type[1] = BlockType::RefNames as u8;
        h2.block_type[2] = BlockType::RefLengths as u8;
        h2.block_size[0] = query.len() as u64;
        h2.block_size[1] = 5;
        h2.block_size[2] = 4;

        let mut f = File::create(&input).unwrap();
        DaaHeader1::new().write_to(&mut f).unwrap();
        h2.write_to(&mut f).unwrap();
        f.write_all(&query).unwrap();
        f.write_all(b"ref0\0").unwrap();
        f.write_all(&100u32.to_ne_bytes()).unwrap();

        run_daa(&ViewConfig {
            daa_file: input.to_string_lossy().to_string(),
            output: output.to_string_lossy().to_string(),
            outfmt: vec!["daa".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            invocation: String::new(),
        })
        .unwrap();

        let mut out_daa = DaaFile::open(&output).unwrap();
        assert_eq!(out_daa.query_records(), 1);
        assert_eq!(out_daa.db_seqs_used(), 1);
        assert_eq!(out_daa.ref_name(0), "ref0");
        let (buf, query_num) = out_daa.read_query_buffer().unwrap().unwrap();
        let r = DaaQueryRecord::from_buffer(&out_daa, &buf, query_num).unwrap();
        assert_eq!(r.query_name, "query0");
        assert!(out_daa.read_query_buffer().unwrap().is_none());

        let _ = std::fs::remove_file(input);
        let _ = std::fs::remove_file(output);
    }

    #[test]
    fn test_run_tabular_writes_alignment_rows() {
        let id = VIEW_TEST_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let input = std::env::temp_dir().join(format!(
            "diamond_rs_view_tab_input_{}_{}.daa",
            std::process::id(),
            id
        ));
        let output = std::env::temp_dir().join(format!(
            "diamond_rs_view_tab_output_{}_{}.txt",
            std::process::id(),
            id
        ));
        let output_json = std::env::temp_dir().join(format!(
            "diamond_rs_view_json_output_{}_{}.json",
            std::process::id(),
            id
        ));
        let output_paf = std::env::temp_dir().join(format!(
            "diamond_rs_view_paf_output_{}_{}.paf",
            std::process::id(),
            id
        ));
        let output_sam = std::env::temp_dir().join(format!(
            "diamond_rs_view_sam_output_{}_{}.sam",
            std::process::id(),
            id
        ));
        let output_pairwise = std::env::temp_dir().join(format!(
            "diamond_rs_view_pairwise_output_{}_{}.txt",
            std::process::id(),
            id
        ));
        let output_xml = std::env::temp_dir().join(format!(
            "diamond_rs_view_xml_output_{}_{}.xml",
            std::process::id(),
            id
        ));
        let output_null = std::env::temp_dir().join(format!(
            "diamond_rs_view_null_output_{}_{}.txt",
            std::process::id(),
            id
        ));
        let output_edge = std::env::temp_dir().join(format!(
            "diamond_rs_view_edge_output_{}_{}.bin",
            std::process::id(),
            id
        ));

        let mut query = Vec::new();
        let seek_pos =
            write_daa_query_record(&mut query, "query0", &[0, 1, 2], SequenceType::AminoAcid);
        let flag = compute_flag(30, 0, 0, false);
        query.extend_from_slice(&0u32.to_ne_bytes());
        query.push(flag);
        query.push(30);
        query.push(0);
        query.push(0);
        query.push(PackedOperation::from_op_count(EditOperation::Match, 3).code);
        query.push(PackedOperation::terminator().code);
        crate::data::daa::finish_daa_query_record(&mut query, seek_pos);
        query.extend_from_slice(&0u32.to_ne_bytes());

        let mut h2 = DaaHeader2::with_params(1, 100, 11, 1, 2, -3, 0.7, 1.2, 0.001, "BLOSUM62", 2);
        h2.db_seqs_used = 1;
        h2.query_records = 1;
        h2.block_type[0] = BlockType::Alignments as u8;
        h2.block_type[1] = BlockType::RefNames as u8;
        h2.block_type[2] = BlockType::RefLengths as u8;
        h2.block_size[0] = query.len() as u64;
        h2.block_size[1] = 5;
        h2.block_size[2] = 4;

        let mut f = File::create(&input).unwrap();
        DaaHeader1::new().write_to(&mut f).unwrap();
        h2.write_to(&mut f).unwrap();
        f.write_all(&query).unwrap();
        f.write_all(b"ref0\0").unwrap();
        f.write_all(&100u32.to_ne_bytes()).unwrap();

        run_tabular(&ViewConfig {
            daa_file: input.to_string_lossy().to_string(),
            output: output.to_string_lossy().to_string(),
            outfmt: vec![
                "tab".to_string(),
                "qseqid".to_string(),
                "sseqid".to_string(),
                "score".to_string(),
                "length".to_string(),
                "pident".to_string(),
                "hspnum".to_string(),
            ],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            invocation: String::new(),
        })
        .unwrap();

        let out = std::fs::read_to_string(&output).unwrap();
        assert_eq!(out, "query0\tref0\t30\t3\t100\t0\n");

        run_tabular(&ViewConfig {
            daa_file: input.to_string_lossy().to_string(),
            output: output_json.to_string_lossy().to_string(),
            outfmt: vec![
                "104".to_string(),
                "qseqid".to_string(),
                "sseqid".to_string(),
                "score".to_string(),
            ],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            invocation: String::new(),
        })
        .unwrap();

        let out = std::fs::read_to_string(&output_json).unwrap();
        assert_eq!(
            out,
            "[\n\t{\n\t\"qseqid\":\"query0\",\n\t\"sseqid\":\"ref0\",\n\t\"score\":30\n\t}\n]"
        );

        run_paf(&ViewConfig {
            daa_file: input.to_string_lossy().to_string(),
            output: output_paf.to_string_lossy().to_string(),
            outfmt: vec!["103".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            invocation: String::new(),
        })
        .unwrap();

        let out = std::fs::read_to_string(&output_paf).unwrap();
        assert!(out.starts_with("query0\t3\t0\t2\t+\tref0\t100\t0\t2\t3\t3\t255\tAS:i:"));
        assert!(out.contains("\tZR:i:30\tZE:f:"));

        run_sam(&ViewConfig {
            daa_file: input.to_string_lossy().to_string(),
            output: output_sam.to_string_lossy().to_string(),
            outfmt: vec!["101".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: true,
            invocation: "diamond view -f sam".to_string(),
        })
        .unwrap();

        let out = std::fs::read_to_string(&output_sam).unwrap();
        assert!(out.starts_with("@HD\tVN:1.5\tSO:query\n"));
        assert!(out.contains("@PG\tPN:DIAMOND\tVN:2.1.24\tCL:diamond view -f sam\n"));
        assert!(out.contains("\nquery0\t0\tref0\t1\t255\t3M\t*\t0\t0\tARN\t*"));
        assert!(out.contains("\tNM:i:0\tZL:i:100\tZR:i:30\t"));
        assert!(out.contains("\tZI:i:100\tZF:i:1\tZS:i:1\tMD:Z:3\tZQ:i:3\n"));

        run_pairwise(&ViewConfig {
            daa_file: input.to_string_lossy().to_string(),
            output: output_pairwise.to_string_lossy().to_string(),
            outfmt: vec!["0".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            invocation: String::new(),
        })
        .unwrap();

        let out = std::fs::read_to_string(&output_pairwise).unwrap();
        assert!(out.starts_with("BLASTP 2.3.0+\n\n\nQuery= query0\n\nLength=3\n\n"));
        assert!(out.contains(">ref0\nLength=100\n\n"));
        assert!(out.contains("Score = "));
        assert!(out.contains("Identities = 3/3 (100%), Positives = 3/3 (100%), Gaps = 0/3 (0%)"));
        assert!(out.contains("Query  1  ARN 3\n"));
        assert!(out.contains("Sbjct  1  ARN 3\n"));

        run_xml(&ViewConfig {
            daa_file: input.to_string_lossy().to_string(),
            output: output_xml.to_string_lossy().to_string(),
            outfmt: vec!["5".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            invocation: String::new(),
        })
        .unwrap();

        let out = std::fs::read_to_string(&output_xml).unwrap();
        assert!(out.starts_with("<?xml version=\"1.0\"?>\n"));
        assert!(out.contains("<BlastOutput_program>blastp</BlastOutput_program>"));
        assert!(out.contains("<BlastOutput_query-def>query0</BlastOutput_query-def>"));
        assert!(out.contains("<Iteration_query-ID>Query_1</Iteration_query-ID>"));
        assert!(out.contains("<Hit_id>ref0</Hit_id>"));
        assert!(out.contains("<Hsp_qseq>ARN</Hsp_qseq>"));
        assert!(out.contains("<Hsp_hseq>ARN</Hsp_hseq>"));
        assert!(out.ends_with("</BlastOutput_iterations>\n</BlastOutput>"));

        run_null(&ViewConfig {
            daa_file: input.to_string_lossy().to_string(),
            output: output_null.to_string_lossy().to_string(),
            outfmt: vec!["null".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            invocation: String::new(),
        })
        .unwrap();

        assert!(std::fs::read(&output_null).unwrap().is_empty());

        run_edge(&ViewConfig {
            daa_file: input.to_string_lossy().to_string(),
            output: output_edge.to_string_lossy().to_string(),
            outfmt: vec!["edge".to_string()],
            max_target_seqs: 25,
            toppercent: None,
            forward_only: false,
            report_unaligned: false,
            sam_qlen_field: false,
            invocation: String::new(),
        })
        .unwrap();

        let out = std::fs::read(&output_edge).unwrap();
        assert_eq!(out.len(), crate::output::edge::EdgeData::SIZE);
        assert_eq!(u64::from_ne_bytes(out[0..8].try_into().unwrap()), 0);
        assert_eq!(u64::from_ne_bytes(out[8..16].try_into().unwrap()), 0);
        assert_eq!(f32::from_ne_bytes(out[16..20].try_into().unwrap()), 100.0);
        assert_eq!(f32::from_ne_bytes(out[20..24].try_into().unwrap()), 3.0);

        let _ = std::fs::remove_file(input);
        let _ = std::fs::remove_file(output);
        let _ = std::fs::remove_file(output_json);
        let _ = std::fs::remove_file(output_paf);
        let _ = std::fs::remove_file(output_sam);
        let _ = std::fs::remove_file(output_pairwise);
        let _ = std::fs::remove_file(output_xml);
        let _ = std::fs::remove_file(output_null);
        let _ = std::fs::remove_file(output_edge);
    }
}
