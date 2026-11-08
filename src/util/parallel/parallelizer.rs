//! Translation of `diamond/src/util/parallel/parallelizer.{h,cpp}`.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::Duration;

use super::filestack::FileStack;

#[derive(Debug)]
pub struct Parallelizer {
    work_directory: PathBuf,
    barrier_file: PathBuf,
    rank: i32,
    id: Vec<u8>,
    n_registered: i32,
    master_flag: bool,
    i_barrier: i32,
    initialized: bool,
    continuous_cleanup_list: Vec<PathBuf>,
    final_cleanup_list: Vec<PathBuf>,
    fs_map: HashMap<String, Arc<FileStack>>,
    errors: crate::util::io::IoErrorOwner,
}

static PARALLELIZER_INSTANCE: OnceLock<Arc<Mutex<Parallelizer>>> = OnceLock::new();

impl Default for Parallelizer {
    fn default() -> Self {
        Self::new()
    }
}

impl Parallelizer {
    pub const LOG: &'static str = "log";
    pub const COMMAND: &'static str = "command";
    pub const WORKERS: &'static str = "workers";
    pub const REGISTER: &'static str = "register";
    const BARRIER: &'static str = "barrier";
    const POLL_SLEEP_SECONDS: f64 = 0.5;
    const POLL_MAX_ITERATIONS: usize = 7_200;

    /// C++ `Parallelizer::get()` singleton.
    pub fn get() -> Arc<Mutex<Parallelizer>> {
        Arc::clone(PARALLELIZER_INSTANCE.get_or_init(|| Arc::new(Mutex::new(Self::new()))))
    }

    /// C++ `Parallelizer::Parallelizer()`.
    pub fn new() -> Self {
        Self::new_owned(Default::default())
    }

    pub fn new_owned(errors: crate::util::io::IoErrorOwner) -> Self {
        Self {
            work_directory: PathBuf::from("parallelizer"),
            barrier_file: PathBuf::new(),
            // The C++ field is not explicitly initialized; zero is the only
            // useful deterministic representation for its getter in Rust.
            rank: 0,
            id: Vec::new(),
            n_registered: 0,
            master_flag: true,
            i_barrier: 0,
            initialized: false,
            continuous_cleanup_list: Vec::new(),
            final_cleanup_list: Vec::new(),
            fs_map: HashMap::new(),
            errors,
        }
    }

    /// C++ `Parallelizer::init(tempdir)`.
    pub fn init(&mut self, tempdir: impl AsRef<Path>) -> Result<(), String> {
        let tempdir = tempdir.as_ref();
        if !tempdir.as_os_str().is_empty() {
            self.work_directory = tempdir.join(&self.work_directory);
        }
        self.errors
            .capture(crate::util::system::mkdir_path(&self.work_directory))?;

        #[cfg(unix)]
        {
            // Source uses gethostname, not an overridable environment value.
            let mut hostname = [0u8; 1024];
            let result = unsafe { libc::gethostname(hostname.as_mut_ptr().cast(), 1023) };
            if result != 0 {
                // Source ignores this failure and reads uninitialized bytes.
                // Preserve a defined diagnostic rather than emulate that UB.
                return self.errors.capture(Err(std::io::Error::other(
                    "Source gethostname failed before worker identity construction.",
                )));
            }
            self.id = hostname[..hostname.iter().position(|b| *b == 0).unwrap_or(1023)].to_vec();
            self.id
                .extend_from_slice(format!("_{}", std::process::id()).as_bytes());
        }
        #[cfg(not(unix))]
        {
            // Pinned Source Windows does not initialize hostname at all.
            // Keep the existing safe native identity; this remains qualified.
            let hostname = std::env::var("HOSTNAME")
                .or_else(|_| std::env::var("COMPUTERNAME"))
                .unwrap_or_else(|_| "localhost".to_owned());
            self.id = format!("{}_{}", hostname, std::process::id()).into_bytes();
        }
        self.create_stack_bytes(Self::LOG, &self.id.clone())?;
        self.create_stack(Self::COMMAND, "")?;
        self.create_stack(Self::WORKERS, "")?;
        self.create_stack(Self::REGISTER, "")?;

        self.barrier_file = self.work_directory.join(Self::BARRIER);
        self.log("PARALLELIZER BEGIN")?;
        self.initialized = true;
        Ok(())
    }

    /// C++ intentionally implements `clear()` as a no-op.
    pub fn clear(&mut self) {}

    pub fn get_rank(&self) -> i32 {
        self.rank
    }

    pub fn get_id(&self) -> &str {
        std::str::from_utf8(&self.id).expect("Worker identity is not UTF-8; use get_id_bytes")
    }

    pub fn get_id_bytes(&self) -> &[u8] {
        &self.id
    }

    pub fn io_errors(&self) -> crate::util::io::IoErrorOwner {
        self.errors.clone()
    }

    pub fn get_work_directory(&self) -> &Path {
        &self.work_directory
    }

    pub fn get_n_registered(&self) -> i32 {
        self.n_registered
    }

    pub fn is_master(&self) -> bool {
        self.master_flag
    }

    pub fn register_workers(&mut self, sleep_s: f64) -> Result<bool, String> {
        self.get_stack(Self::REGISTER)?
            .push_bytes(&self.id)
            .map(|_| ())?;
        Self::sleep(sleep_s);
        if self.is_master() {
            while let Some(line) = self.get_stack(Self::REGISTER)?.pop_bytes()? {
                self.get_stack(Self::WORKERS)?.push_bytes(&line)?;
                self.n_registered += 1;
            }
        }
        Ok(true)
    }

    /// C++ file-stack barrier, including its default polling duration.
    pub fn barrier(&mut self, tag: &str) -> Result<bool, String> {
        if !self.initialized {
            return Ok(false);
        }

        let cmd_file_name = self.get_barrier_file_name("cmd", tag, self.i_barrier);
        let cmd_fs = FileStack::try_new_owned(&cmd_file_name, self.errors.clone())?;
        let ack_file_name = self.get_barrier_file_name("ack", tag, self.i_barrier);
        let ack_fs = FileStack::try_new_owned(&ack_file_name, self.errors.clone())?;

        let msg = "WAIT";
        if self.is_master() {
            ack_fs.clear()?;
            cmd_fs.push_string(msg)?;
        }
        cmd_fs.poll_query(msg, Self::POLL_SLEEP_SECONDS, Self::POLL_MAX_ITERATIONS)?;
        ack_fs.push_bytes(&self.id)?;

        let msg_ok = "GOON";
        if self.is_master() {
            let n_workers = self.get_stack(Self::WORKERS)?.size()?;
            ack_fs.poll_size(
                n_workers,
                Self::POLL_SLEEP_SECONDS,
                Self::POLL_MAX_ITERATIONS,
            )?;
            cmd_fs.push_string(msg_ok)?;
        }
        cmd_fs.poll_query(msg_ok, Self::POLL_SLEEP_SECONDS, Self::POLL_MAX_ITERATIONS)?;

        if self.is_master() {
            Self::clean(&mut self.continuous_cleanup_list);
            self.continuous_cleanup_list.push(cmd_file_name);
            self.continuous_cleanup_list.push(ack_file_name);
        }

        self.i_barrier += 1;
        Ok(true)
    }

    pub fn create_stack(&mut self, tag: &str, sfx: &str) -> Result<bool, String> {
        self.create_stack_bytes(tag, sfx.as_bytes())
    }

    fn create_stack_bytes(&mut self, tag: &str, sfx: &[u8]) -> Result<bool, String> {
        if self.fs_map.contains_key(tag) {
            return Ok(false);
        }
        let mut name = tag.as_bytes().to_vec();
        if !sfx.is_empty() {
            name.push(b'_');
            name.extend_from_slice(sfx);
        }
        let name = self
            .errors
            .capture(crate::util::io::source_bytes_os(&name).map_err(std::io::Error::other))?;
        self.create_stack_from_file(tag, self.work_directory.join(name))
    }

    pub fn create_stack_from_file(
        &mut self,
        tag: &str,
        file_name: impl AsRef<Path>,
    ) -> Result<bool, String> {
        self.delete_stack(tag);
        self.fs_map.insert(
            tag.to_string(),
            Arc::new(FileStack::try_new_owned(
                file_name.as_ref(),
                self.errors.clone(),
            )?),
        );
        Ok(true)
    }

    pub fn get_stack(&self, tag: &str) -> Result<Arc<FileStack>, String> {
        self.fs_map
            .get(tag)
            .cloned()
            .ok_or_else(|| format!("FileStack tag not found: {tag}"))
    }

    pub fn delete_stack(&mut self, tag: &str) -> bool {
        self.fs_map.remove(tag).is_some()
    }

    pub fn sleep(sleep_s: f64) {
        std::thread::sleep(Duration::from_secs_f64(sleep_s));
    }

    fn clean(file_list: &mut Vec<PathBuf>) -> bool {
        for path in file_list.iter() {
            let _ = std::fs::remove_file(path);
        }
        file_list.clear();
        true
    }

    /// Testable equivalent of C++ `list_filestacks()`'s diagnostic listing.
    pub fn list_filestacks(&self) -> Vec<(String, PathBuf)> {
        let mut out = self
            .fs_map
            .iter()
            .map(|(tag, stack)| (tag.clone(), stack.file_name().to_path_buf()))
            .collect::<Vec<_>>();
        out.sort();
        out
    }

    pub fn log(&self, buf: &str) -> Result<(), String> {
        let log_stack = self.get_stack(Self::LOG)?;
        let ms = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_err(|e| e.to_string())?
            .as_millis();
        log_stack.push_string(&format!("{ms} {buf}\n")).map(|_| ())
    }

    fn get_barrier_file_name(&self, step: &str, tag: &str, i: i32) -> PathBuf {
        let mut path = self.barrier_file.as_os_str().to_os_string();
        path.push(format!("_{step}_{tag}_{i}"));
        path.into()
    }
}

impl Drop for Parallelizer {
    fn drop(&mut self) {
        if self.initialized {
            let _ = self.log("PARALLELIZER END");
            Self::clean(&mut self.continuous_cleanup_list);
            Self::clean(&mut self.final_cleanup_list);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn temp_root(test: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "diamond-rs-parallelizer-{test}-{}-{}",
            std::process::id(),
            crate::util::test_path_component()
        ))
    }

    #[cfg(unix)]
    #[test]
    fn init_does_not_create_missing_parent_directories() {
        let root = temp_root("missing-parent");
        let errors = crate::util::io::IoErrorOwner::default();
        let mut p = Parallelizer::new_owned(errors.clone());
        let error = errors.finish(p.init(root.join("nested"))).unwrap_err();
        let mut expected = b"could not create temporary directory ".to_vec();
        expected.extend_from_slice(
            &crate::util::io::source_os_bytes(root.join("nested/parallelizer").as_os_str())
                .unwrap(),
        );
        assert_eq!(
            crate::util::io::RawDiagnosticError::find(&error),
            Some(expected.as_slice())
        );
        assert!(!root.exists());
    }

    #[test]
    fn stack_constructor_failure_returns_the_captured_raw_diagnostic() {
        let root = temp_root("stack-open-error");
        std::fs::create_dir(&root).unwrap();
        let blocked = root.join("directory-instead-of-stack");
        std::fs::create_dir(&blocked).unwrap();
        let errors = crate::util::io::IoErrorOwner::default();
        let mut p = Parallelizer::new_owned(errors.clone());
        let error = errors
            .finish(p.create_stack_from_file("blocked", &blocked))
            .unwrap_err();
        let mut expected = if cfg!(windows) {
            b"Error opening file ".to_vec()
        } else {
            b"could not open file ".to_vec()
        };
        expected.extend_from_slice(&crate::util::io::source_os_bytes(blocked.as_os_str()).unwrap());
        assert_eq!(
            crate::util::io::RawDiagnosticError::find(&error),
            Some(expected.as_slice())
        );
        assert!(!p.fs_map.contains_key("blocked"));
        drop(p);
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn init_register_barrier_log_and_cleanup() {
        let root = temp_root("lifecycle");
        let _ = std::fs::remove_dir_all(&root);
        let mut p = Parallelizer::new();
        std::fs::create_dir(&root).unwrap();
        p.init(&root).unwrap();
        assert!(p.is_master());
        assert_eq!(p.get_rank(), 0);
        assert!(p.get_id_bytes().contains(&b'_'));
        assert!(p.get_work_directory().ends_with("parallelizer"));
        assert!(p.register_workers(0.0).unwrap());
        assert_eq!(p.get_n_registered(), 1);
        assert!(p.barrier("unit").unwrap());
        assert_eq!(p.i_barrier, 1);
        p.log("custom message").unwrap();
        assert!(p
            .get_stack(Parallelizer::LOG)
            .unwrap()
            .top_string()
            .unwrap()
            .unwrap()
            .contains("custom message"));
        drop(p);
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn uninitialized_and_stack_error_paths() {
        let mut p = Parallelizer::new();
        assert!(!p.barrier("before-init").unwrap());
        assert_eq!(
            p.get_stack("missing").unwrap_err(),
            "FileStack tag not found: missing"
        );
        assert!(!p.delete_stack("missing"));
        p.clear();
    }

    #[test]
    fn barrier_propagates_filestack_stop_error() {
        let root = temp_root("stop");
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(root.join("parallelizer")).unwrap();
        let mut p = Parallelizer::new();
        p.work_directory = root.join("parallelizer");
        p.barrier_file = p.work_directory.join(Parallelizer::BARRIER);
        p.id = b"worker".to_vec();
        p.master_flag = false;
        p.initialized = true;
        let command = FileStack::new(p.get_barrier_file_name("cmd", "failed", 0));
        command.push_string("STOP").unwrap();
        let error = p.barrier("failed").unwrap_err();
        assert!(error.contains("STOP on FileStack"));
        p.initialized = false;
        drop(command);
        drop(p);
        let _ = std::fs::remove_dir_all(&root);
    }
}
