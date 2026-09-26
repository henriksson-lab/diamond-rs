//! Translation of `diamond/src/data/block/block_wrapper.cpp`.

use super::Block;
use crate::basic::sequence::Sequence;
use crate::basic::value::{Letter, Loc, OId, TaxId};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BlockWrapperSeqInfo {
    pub pos: u64,
    pub seq_len: u32,
}

impl BlockWrapperSeqInfo {
    pub fn new(pos: u64, len: usize) -> Self {
        Self {
            pos,
            seq_len: len as u32,
        }
    }
}

pub struct BlockWrapper<'a> {
    block: &'a Block,
    oid: OId,
}

impl<'a> BlockWrapper<'a> {
    pub fn new(block: &'a Block) -> Self {
        Self { block, oid: 0 }
    }

    pub fn file_count(&self) -> i64 {
        1
    }

    pub fn files_synced(&self) -> bool {
        true
    }

    pub fn read_seqinfo(&mut self) -> Result<BlockWrapperSeqInfo, String> {
        if self.oid >= self.block.seqs().len() as OId {
            self.oid += 1;
            return Ok(BlockWrapperSeqInfo::new(0, 0));
        }
        let len = self.block.seqs().length(self.oid as usize);
        if len == 0 {
            return Err("Database with sequence length 0 is not supported".to_owned());
        }
        let out = BlockWrapperSeqInfo::new(self.oid, len as usize);
        self.oid += 1;
        Ok(out)
    }

    pub fn putback_seqinfo(&mut self) {
        self.oid -= 1;
    }

    pub fn close(&mut self) {}

    pub fn set_seqinfo_ptr(&mut self, oid: OId) {
        self.oid = oid;
    }

    pub fn tell_seq(&self) -> OId {
        self.oid
    }

    pub fn eof(&self) -> bool {
        self.oid >= self.block.seqs().len() as OId
    }

    pub fn init_seq_access(&mut self) {
        self.set_seqinfo_ptr(0);
    }

    pub fn read_seq(&mut self) -> Result<(Vec<Letter>, String, Option<Vec<u8>>), String> {
        Err("Operation not supported".to_owned())
    }

    pub fn create_partition_balanced(&mut self, _max_letters: i64) -> Result<(), String> {
        Err("Operation not supported".to_owned())
    }

    pub fn save_partition(
        &mut self,
        _partition_file_name: &str,
        _annotation: &str,
    ) -> Result<(), String> {
        Err("Operation not supported".to_owned())
    }

    pub fn get_n_partition_chunks(&mut self) -> Result<i32, String> {
        Err("Operation not supported".to_owned())
    }

    pub fn init_seqinfo_access(&mut self) {}

    pub fn seek_chunk(
        &mut self,
        _chunk_i: i32,
        _offset: usize,
        _n_seqs: i64,
    ) -> Result<(), String> {
        Err("Operation not supported".to_owned())
    }

    pub fn seqid(&self, oid: OId, _all: bool, _full_titles: bool) -> Result<String, String> {
        Ok(std::str::from_utf8(self.block.ids()?.get(oid as usize))
            .unwrap_or("")
            .to_owned())
    }

    pub fn id_len(
        &self,
        seq_info: &BlockWrapperSeqInfo,
        _seq_info_next: &BlockWrapperSeqInfo,
    ) -> Result<usize, String> {
        Ok(self.block.ids()?.length(seq_info.pos as usize) as usize)
    }

    pub fn seek_offset(&mut self, _p: usize) {}

    pub fn read_seq_data(&self, dst: &mut Vec<Letter>, len: usize, pos: &mut usize, _seek: bool) {
        dst.clear();
        dst.reserve(len + 2);
        dst.push(Sequence::DELIMITER);
        dst.extend_from_slice(self.block.seqs().get(*pos));
        dst.push(Sequence::DELIMITER);
        *pos += 1;
    }

    pub fn read_id_data(&self, oid: i64, _all: bool, _full_titles: bool) -> Result<String, String> {
        Ok(std::str::from_utf8(self.block.ids()?.get(oid as usize))
            .unwrap_or("")
            .to_owned())
    }

    pub fn skip_id_data(&mut self) {}

    pub fn sequence_count(&self) -> u64 {
        self.block.seqs().len() as u64
    }

    pub fn letters(&self) -> u64 {
        self.block.seqs().letters()
    }

    pub fn db_version(&self) -> Result<i32, String> {
        Err("Operation not supported".to_owned())
    }

    pub fn program_build_version(&self) -> Result<i32, String> {
        Err("Operation not supported".to_owned())
    }

    pub fn build_version(&mut self) -> Result<i32, String> {
        Err("Operation not supported".to_owned())
    }

    pub fn filter_by_accession(&mut self, _file_name: &str) -> Result<(), String> {
        Err("Operation not supported".to_owned())
    }

    pub fn file_name(&self) -> String {
        String::new()
    }

    pub fn taxids(&self, _oid: usize) -> Result<Vec<TaxId>, String> {
        Err("Operation not supported".to_owned())
    }

    pub fn seq_data(&self, _oid: usize, _dst: &mut Vec<Letter>) -> Result<(), String> {
        Err("Operation not supported".to_owned())
    }

    pub fn seq_length(&self, _oid: usize) -> Result<Loc, String> {
        Err("Operation not supported".to_owned())
    }

    pub fn end_random_access(&mut self, _dictionary: bool) {}
}
