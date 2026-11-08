macro_rules! search_message {
    ($($argument:tt)*) => {{
        $crate::util::log_stream::message_stream().lock().unwrap()
            .write(format!("{}\n", format_args!($($argument)*)))?;
    }};
}
pub mod argument_registry;
pub mod benchmark;
#[cfg(feature = "dna")]
pub mod blastn;
pub mod blastp;
pub mod blastx;
pub mod cluster_cmd;
mod cluster_database;
pub mod dbinfo;
pub mod deprecated;
pub mod dispatch;
mod external_cluster;
pub mod extra;
#[cfg(feature = "extra")]
mod extra_tools;
pub mod getseq;
pub mod greedy_vertex_cover;
pub mod makedb;
pub mod makeidx;
#[cfg(feature = "mcl")]
pub mod mcl_cmd;
pub mod merge_daa;
pub mod query_output;
pub mod realign;
pub mod reassign;
pub mod recluster;
pub mod search_algorithm;
mod search_blocks;
mod search_logging;
mod search_output;
pub mod test_cmd;
pub mod view;

pub mod source_config;
