//! Production owner-based services for Source's double-indexed control loops.
//!
//! Query/reference associated types are real loaded owners, not descriptive
//! stand-ins. Backends provide concrete I/O and search phases; this module owns
//! master/query/iteration/reference ordering, including lazy reference claims.
use std::io;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReferenceFlow {
    Continue,
    /// A concrete single-reference output path completed the round without
    /// needing the ordinary final EOF load.
    FinishRound,
}

#[derive(Clone, Copy, Debug)]
pub struct QueryCompletion {
    /// Actual Source counter after reference/global/join processing. For MPP
    /// this comes from claimed jobs/join ownership, not the number of polls.
    pub reference_counter: usize,
}

pub trait DoubleIndexedBackend: Sized {
    type Query;
    type Round;
    type Reference;
    type Output;

    fn database_sequence_count(&self) -> u64;
    fn sensitivity_round_count(&self) -> usize;
    /// Ordinary Source's first read is indeterminate; native initializes it
    /// to zero. Provided owners may carry an explicit, already established
    /// compatibility state for that first master-query ordering decision.
    fn initial_reference_counter(&self) -> usize {
        0
    }
    /// MPP samples completed claimed jobs instead of the ordinary scan end.
    fn log_after_reference_scan(&self) -> bool {
        true
    }
    /// Record the primary failure while actual backend owners are still alive.
    fn capture_source_error(&self, error: io::Error) -> io::Error {
        error
    }
    fn begin_master(&mut self) -> io::Result<()>;
    fn load_query(&mut self, ordinal: usize) -> io::Result<Option<Self::Query>>;
    fn query_is_selected(&self, ordinal: usize, query: &Self::Query) -> bool;
    fn prepare_query(
        &mut self,
        query: &mut Self::Query,
        ordinal: usize,
        previous_reference: usize,
    ) -> io::Result<()>;
    fn query_has_pending_alignments(&self, query: &Self::Query) -> bool;
    fn begin_round(&mut self, query: &mut Self::Query, ordinal: usize) -> io::Result<Self::Round>;
    /// For ordinary search `cursor` identifies the next sequential block.
    /// MPP backends claim one opaque queue job here and retain its token in
    /// Reference; they must not preclaim/materialize all job indices.
    fn load_reference(
        &mut self,
        query: &mut Self::Query,
        round: &mut Self::Round,
        cursor: usize,
    ) -> io::Result<Option<Self::Reference>>;
    fn run_reference(
        &mut self,
        query: &mut Self::Query,
        round: &mut Self::Round,
        reference: Self::Reference,
    ) -> io::Result<ReferenceFlow>;
    fn finish_round(
        &mut self,
        query: &mut Self::Query,
        round: Self::Round,
        loaded_references: usize,
    ) -> io::Result<()>;
    fn finish_query(&mut self, query: Self::Query) -> io::Result<QueryCompletion>;
    fn stop_exists(&self) -> bool;
    fn finish(self) -> io::Result<Self::Output>;
}

pub fn run<B: DoubleIndexedBackend>(mut backend: B) -> io::Result<B::Output> {
    let result = run_body(&mut backend);
    if let Err(error) = result {
        return Err(backend.capture_source_error(error));
    }
    backend.finish()
}

fn run_body<B: DoubleIndexedBackend>(backend: &mut B) -> io::Result<()> {
    let _ = crate::util::system::log_rss();
    backend.begin_master()?;
    let mut query_ordinal = 0usize;
    // Initialize the undefined Source first-read boundary from its owner's
    // explicit compatibility state, then retain actual traversal counters.
    let mut previous_reference = backend.initial_reference_counter();
    while backend.database_sequence_count() != 0 {
        let _ = crate::util::system::log_rss();
        let Some(mut query) = backend.load_query(query_ordinal)? else {
            break;
        };
        if backend.query_is_selected(query_ordinal, &query) {
            backend.prepare_query(&mut query, query_ordinal, previous_reference)?;
            let _ = crate::util::system::log_rss();
            for round_ordinal in 0..backend.sensitivity_round_count() {
                if !backend.query_has_pending_alignments(&query) {
                    break;
                }
                let mut round = backend.begin_round(&mut query, round_ordinal)?;
                let _ = crate::util::system::log_rss();
                let mut reference_cursor = 0usize;
                while let Some(reference) =
                    backend.load_reference(&mut query, &mut round, reference_cursor)?
                {
                    let flow = backend.run_reference(&mut query, &mut round, reference)?;
                    reference_cursor = reference_cursor
                        .checked_add(1)
                        .ok_or_else(|| io::Error::other("Reference block index overflow"))?;
                    if flow == ReferenceFlow::FinishRound {
                        break;
                    }
                }
                if backend.log_after_reference_scan() {
                    let _ = crate::util::system::log_rss();
                }
                backend.finish_round(&mut query, round, reference_cursor)?;
            }
            let _ = crate::util::system::log_rss();
            previous_reference = backend.finish_query(query)?.reference_counter;
            if backend.stop_exists() {
                if let Ok(mut stream) = crate::util::log_stream::message_stream().lock() {
                    let _ = stream.write("Encountered 'stop' file, shutting down run\n");
                }
                break;
            }
        }
        query_ordinal = query_ordinal
            .checked_add(1)
            .ok_or_else(|| io::Error::other("Query block index overflow"))?;
    }
    Ok(())
}
