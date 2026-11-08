//! Source search verbosity and command-scoped log stream lifecycle.
use crate::util::log_stream::{log_stream, message_stream, verbose_stream, MessageStream};
use std::io;

pub(super) struct SearchLogging {
    previous: [MessageStream; 3],
}

impl SearchLogging {
    pub(super) fn establish_source(verbosity: u8, logfile: bool) -> io::Result<Self> {
        let guard = Self {
            previous: [
                message_stream().lock().unwrap().clone(),
                verbose_stream().lock().unwrap().clone(),
                log_stream().lock().unwrap().clone(),
            ],
        };
        Self::reconfigure_source(verbosity, logfile);
        Ok(guard)
    }

    pub(super) fn reconfigure_source(verbosity: u8, logfile: bool) {
        // The C++ constructor intentionally leaves streams not mentioned by
        // its switch unchanged. Nested regression Configs share this scope.
        match verbosity {
            0 => *message_stream().lock().unwrap() = MessageStream::new(false, false),
            3 => {
                *log_stream().lock().unwrap() = MessageStream::new(true, logfile);
                *verbose_stream().lock().unwrap() = MessageStream::new(true, logfile);
                *message_stream().lock().unwrap() = MessageStream::new(true, logfile);
            }
            2 => *verbose_stream().lock().unwrap() = MessageStream::new(true, false),
            _ => (),
        }
    }
    pub(super) fn establish(verbosity: u8, logfile: bool) -> io::Result<Self> {
        let mut messages = message_stream().lock().unwrap();
        let mut verbose = verbose_stream().lock().unwrap();
        let mut log = log_stream().lock().unwrap();
        let guard = Self {
            previous: [messages.clone(), verbose.clone(), log.clone()],
        };
        let file = verbosity == 3 && logfile;
        *messages = MessageStream::new(verbosity > 0, file);
        *verbose = MessageStream::new(verbosity >= 2, file);
        *log = MessageStream::new(verbosity == 3, file);
        Ok(guard)
    }
    pub(super) fn begin(
        verbosity: u8,
        logfile: bool,
        invocation: &str,
        threads: usize,
    ) -> io::Result<Self> {
        let guard = Self::establish(verbosity, logfile)?;
        if verbosity > 0 {
            super::dispatch::print_banner();
        }
        log_stream()
            .lock()
            .unwrap()
            .write(format!("{invocation}\n"))?;
        message_stream()
            .lock()
            .unwrap()
            .write(format!("#CPU threads: {threads}\n"))?;
        if cfg!(debug_assertions) {
            verbose_stream()
                .lock()
                .unwrap()
                .write("Assertions enabled.\n")?;
        }
        Ok(guard)
    }
}

impl Drop for SearchLogging {
    fn drop(&mut self) {
        *message_stream().lock().unwrap() = self.previous[0].clone();
        *verbose_stream().lock().unwrap() = self.previous[1].clone();
        *log_stream().lock().unwrap() = self.previous[2].clone();
    }
}
