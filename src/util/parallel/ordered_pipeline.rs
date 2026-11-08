//! Bounded producer/parallel decoder/ordered consumer, matching Source loaders.
use std::collections::BTreeMap;
use std::sync::{
    atomic::{AtomicBool, Ordering},
    mpsc, Mutex,
};

pub fn run<I, O, P, D, E>(
    workers: usize,
    mut produce: P,
    decode: D,
    mut emit: E,
) -> Result<(), String>
where
    I: Send,
    O: Send,
    P: FnMut() -> Result<Option<I>, String>,
    D: Fn(I) -> Result<O, String> + Sync,
    E: FnMut(O) -> Result<(), String> + Send,
{
    let workers = workers.max(1);
    let depth = workers.saturating_mul(4).max(1);
    let (input_tx, input_rx) = mpsc::sync_channel::<(usize, I)>(depth);
    let input_rx = Mutex::new(input_rx);
    let (output_tx, output_rx) = mpsc::sync_channel::<(usize, Result<O, String>)>(depth);
    let failed = AtomicBool::new(false);
    std::thread::scope(|scope| {
        let failed_ref = &failed;
        let writer = scope.spawn(move || {
            let mut next = 0;
            let mut pending = BTreeMap::new();
            let mut error = None;
            while let Ok((id, value)) = output_rx.recv() {
                pending.insert(id, value);
                while let Some(value) = pending.remove(&next) {
                    if error.is_none() {
                        if let Err(message) = value.and_then(&mut emit) {
                            error = Some(message);
                            failed_ref.store(true, Ordering::Release);
                        }
                    }
                    next += 1;
                }
            }
            if let Some(error) = error {
                Err(error)
            } else if !pending.is_empty() {
                Err("Ordered decoder lost a package.".into())
            } else {
                Ok(())
            }
        });
        let mut decoders = Vec::with_capacity(workers);
        for _ in 0..workers {
            let output_tx = output_tx.clone();
            let input_rx = &input_rx;
            let decode = &decode;
            decoders.push(scope.spawn(move || {
                loop {
                    let next = input_rx
                        .lock()
                        .map_err(|_| "Raw decoder input lock poisoned.")?
                        .recv();
                    let Ok((id, value)) = next else {
                        break;
                    };
                    output_tx
                        .send((id, decode(value)))
                        .map_err(|_| "Raw decoder output closed.")?;
                }
                Ok::<(), String>(())
            }));
        }
        drop(output_tx);
        let produced = (|| {
            let mut id = 0usize;
            while !failed.load(Ordering::Acquire) {
                let Some(value) = produce()? else {
                    break;
                };
                input_tx
                    .send((id, value))
                    .map_err(|_| "Raw decoder input closed.")?;
                id = id.checked_add(1).ok_or("Raw package ordinal overflow.")?;
            }
            Ok::<(), String>(())
        })();
        drop(input_tx);
        let mut worker_error = None;
        for decoder in decoders {
            if let Err(error) = decoder
                .join()
                .map_err(|_| "Raw decoder panicked.".to_string())?
            {
                worker_error = Some(error);
            }
        }
        let emitted = writer
            .join()
            .map_err(|_| "Raw package consumer panicked.".to_string())?;
        produced?;
        if let Some(error) = worker_error {
            return Err(error);
        }
        emitted
    })
}
