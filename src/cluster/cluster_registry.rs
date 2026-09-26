//! Clustering-algorithm registry.
//!
//! `cluster_registry.cpp` only instantiates the process-wide registry; the
//! behavior mirrored here is inline in `cluster_registry.h`.

use crate::cluster::cascaded::Cascaded;
use crate::commands::cluster_cmd::{self, ClusterConfig};
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::sync::{Arc, OnceLock};

/// Explicit-config counterpart of the C++ global-config algorithm interface.
pub trait ClusteringAlgorithm: Debug + Send + Sync {
    fn key(&self) -> &'static str;
    fn get_description(&self) -> &'static str;
    fn run(&self, config: &ClusterConfig) -> Result<(), String>;
}

impl ClusteringAlgorithm for Cascaded {
    fn key(&self) -> &'static str {
        Cascaded::get_key()
    }

    fn get_description(&self) -> &'static str {
        "Cascaded greedy vertex cover algorithm"
    }

    fn run(&self, config: &ClusterConfig) -> Result<(), String> {
        cluster_cmd::run(config).map_err(|error| error.to_string())
    }
}

/// Owning registry corresponding to C++ `ClusterRegistryStatic`.
#[derive(Debug, Clone)]
pub struct ClusterRegistryStatic {
    registry: BTreeMap<String, Arc<dyn ClusteringAlgorithm>>,
}

impl ClusterRegistryStatic {
    pub fn new() -> Self {
        let mut registry = Self {
            registry: BTreeMap::new(),
        };
        registry.register(Arc::new(Cascaded));
        registry
    }

    /// Construct a registry from algorithms in insertion order.
    ///
    /// Like C++ `regMap[key] = ptr`, a repeated key replaces the previous
    /// registration. Rust drops the replaced `Arc` safely.
    pub fn from_algorithms<I>(algorithms: I) -> Self
    where
        I: IntoIterator<Item = Arc<dyn ClusteringAlgorithm>>,
    {
        let mut registry = Self {
            registry: BTreeMap::new(),
        };
        for algorithm in algorithms {
            registry.register(algorithm);
        }
        registry
    }

    pub fn register(
        &mut self,
        algorithm: Arc<dyn ClusteringAlgorithm>,
    ) -> Option<Arc<dyn ClusteringAlgorithm>> {
        self.registry.insert(algorithm.key().to_string(), algorithm)
    }

    pub fn get(&self, key: &str) -> Result<Arc<dyn ClusteringAlgorithm>, String> {
        self.registry
            .get(key)
            .cloned()
            .ok_or_else(|| "Clustering algorithm not found.".to_string())
    }

    pub fn has(&self, key: &str) -> bool {
        self.registry.contains_key(key)
    }

    pub fn get_keys(&self) -> Vec<String> {
        self.registry.keys().cloned().collect()
    }
}

impl Default for ClusterRegistryStatic {
    fn default() -> Self {
        Self::new()
    }
}

/// Static facade corresponding to C++ `ClusterRegistry::reg` and methods.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ClusterRegistry;

impl ClusterRegistry {
    fn registry() -> &'static ClusterRegistryStatic {
        static REGISTRY: OnceLock<ClusterRegistryStatic> = OnceLock::new();
        REGISTRY.get_or_init(ClusterRegistryStatic::new)
    }

    pub fn get(key: &str) -> Result<Arc<dyn ClusteringAlgorithm>, String> {
        Self::registry().get(key)
    }

    pub fn has(key: &str) -> bool {
        Self::registry().has(key)
    }

    pub fn get_keys() -> Vec<String> {
        Self::registry().get_keys()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    #[derive(Debug)]
    struct TestAlgorithm {
        key: &'static str,
        description: &'static str,
        drops: Option<Arc<AtomicUsize>>,
    }

    impl Drop for TestAlgorithm {
        fn drop(&mut self) {
            if let Some(drops) = &self.drops {
                drops.fetch_add(1, Ordering::SeqCst);
            }
        }
    }

    impl ClusteringAlgorithm for TestAlgorithm {
        fn key(&self) -> &'static str {
            self.key
        }

        fn get_description(&self) -> &'static str {
            self.description
        }

        fn run(&self, _config: &ClusterConfig) -> Result<(), String> {
            Ok(())
        }
    }

    fn test_algorithm(
        key: &'static str,
        description: &'static str,
    ) -> Arc<dyn ClusteringAlgorithm> {
        Arc::new(TestAlgorithm {
            key,
            description,
            drops: None,
        })
    }

    #[test]
    fn default_registry_contains_cascaded_with_exact_metadata() {
        let registry = ClusterRegistryStatic::new();
        assert!(registry.has("cascaded"));
        assert!(!registry.has("Cascaded"));
        assert_eq!(registry.get_keys(), vec!["cascaded"]);
        let algorithm = registry.get("cascaded").unwrap();
        assert_eq!(algorithm.key(), "cascaded");
        assert_eq!(
            algorithm.get_description(),
            "Cascaded greedy vertex cover algorithm"
        );
    }

    #[test]
    fn unknown_key_has_exact_upstream_error() {
        assert_eq!(
            ClusterRegistryStatic::new().get("missing").unwrap_err(),
            "Clustering algorithm not found."
        );
    }

    #[test]
    fn keys_are_map_sorted_and_lookup_clones_same_registration() {
        let registry = ClusterRegistryStatic::from_algorithms([
            test_algorithm("zeta", "z"),
            test_algorithm("alpha", "a"),
            test_algorithm("middle", "m"),
        ]);
        assert_eq!(registry.get_keys(), ["alpha", "middle", "zeta"]);
        let first = registry.get("middle").unwrap();
        let second = registry.get("middle").unwrap();
        assert!(Arc::ptr_eq(&first, &second));
    }

    #[test]
    fn duplicate_key_replaces_and_drops_previous_owner() {
        let drops = Arc::new(AtomicUsize::new(0));
        let old: Arc<dyn ClusteringAlgorithm> = Arc::new(TestAlgorithm {
            key: "same",
            description: "old",
            drops: Some(drops.clone()),
        });
        let new = test_algorithm("same", "new");
        let mut registry = ClusterRegistryStatic::from_algorithms([old]);

        let replaced = registry.register(new).unwrap();
        assert_eq!(registry.get("same").unwrap().get_description(), "new");
        assert_eq!(drops.load(Ordering::SeqCst), 0);
        drop(replaced);
        assert_eq!(drops.load(Ordering::SeqCst), 1);
    }

    #[test]
    fn static_facade_is_initialized_once_and_reuses_algorithm() {
        assert!(ClusterRegistry::has("cascaded"));
        assert_eq!(ClusterRegistry::get_keys(), ["cascaded"]);
        let first = ClusterRegistry::get("cascaded").unwrap();
        let second = ClusterRegistry::get("cascaded").unwrap();
        assert!(Arc::ptr_eq(&first, &second));
    }
}
