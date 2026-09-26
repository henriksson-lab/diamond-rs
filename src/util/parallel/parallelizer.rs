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
    id: String,
    n_registered: i32,
    master_flag: bool,
    i_barrier: i32,
    initialized: bool,
    continuous_cleanup_list: Vec<PathBuf>,
    final_cleanup_list: Vec<PathBuf>,
    fs_map: HashMap<String, Arc<FileStack>>,
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
        Self {
            work_directory: PathBuf::from("parallelizer"),
            barrier_file: PathBuf::new(),
            // The C++ field is not explicitly initialized; zero is the only
            // useful deterministic representation for its getter in Rust.
            rank: 0,
            id: String::new(),
            n_registered: 0,
            master_flag: true,
            i_barrier: 0,
            initialized: false,
            continuous_cleanup_list: Vec::new(),
            final_cleanup_list: Vec::new(),
            fs_map: HashMap::new(),
        }
    }

    /// C++ `Parallelizer::init(tempdir)`.
    pub fn init(&mut self, tempdir: impl AsRef<Path>) -> Result<(), String> {
        let tempdir = tempdir.as_ref();
        if !tempdir.as_os_str().is_empty() {
            self.work_directory = tempdir.join(&self.work_directory);
        }
        std::fs::create_dir_all(&self.work_directory).map_err(|e| e.to_string())?;

        let hostname = std::env::var("HOSTNAME")
            .or_else(|_| std::env::var("COMPUTERNAME"))
            .unwrap_or_else(|_| "localhost".to_string());
        self.id = format!("{}_{}", hostname, std::process::id());

        self.create_stack(Self::LOG, &self.id.clone())?;
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
        &self.id
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
            .push_string(&self.id)
            .map(|_| ())?;
        Self::sleep(sleep_s);
        if self.is_master() {
            while let Some(line) = self.get_stack(Self::REGISTER)?.pop_string()? {
                self.get_stack(Self::WORKERS)?.push_string(&line)?;
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
        let cmd_fs = FileStack::new(&cmd_file_name);
        let ack_file_name = self.get_barrier_file_name("ack", tag, self.i_barrier);
        let ack_fs = FileStack::new(&ack_file_name);

        let msg = "WAIT";
        if self.is_master() {
            ack_fs.clear()?;
            cmd_fs.push_string(msg)?;
        }
        cmd_fs.poll_query(msg, Self::POLL_SLEEP_SECONDS, Self::POLL_MAX_ITERATIONS)?;
        ack_fs.push_string(&self.id)?;

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
        if self.fs_map.contains_key(tag) {
            return Ok(false);
        }
        let suffix = if sfx.is_empty() {
            String::new()
        } else {
            format!("_{sfx}")
        };
        let file_name = self.work_directory.join(format!("{tag}{suffix}"));
        self.create_stack_from_file(tag, &file_name)
    }

    pub fn create_stack_from_file(
        &mut self,
        tag: &str,
        file_name: impl AsRef<Path>,
    ) -> Result<bool, String> {
        self.delete_stack(tag);
        self.fs_map.insert(
            tag.to_string(),
            Arc::new(FileStack::new(file_name.as_ref())),
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
        PathBuf::from(format!(
            "{}_{}_{}_{}",
            self.barrier_file.display(),
            step,
            tag,
            i
        ))
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
            std::thread::current().name().unwrap_or("thread")
        ))
    }

    #[test]
    fn init_register_barrier_log_and_cleanup() {
        let root = temp_root("lifecycle");
        let _ = std::fs::remove_dir_all(&root);
        let mut p = Parallelizer::new();
        p.init(&root).unwrap();
        assert!(p.is_master());
        assert_eq!(p.get_rank(), 0);
        assert!(p.get_id().contains('_'));
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
        p.id = "worker".to_string();
        p.master_flag = false;
        p.initialized = true;
        let command = FileStack::new(p.get_barrier_file_name("cmd", "failed", 0));
        command.push_string("STOP").unwrap();
        let error = p.barrier("failed").unwrap_err();
        assert!(error.contains("STOP on FileStack"));
        p.initialized = false;
        let _ = std::fs::remove_dir_all(&root);
    }
}
