use std::time::Duration;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ExtensionTimer {
    pub total_time: Duration,
    pub preprocessing_time: Duration,
    pub postprocessing_time: Duration,
    pub extension: Duration,
    pub next: Duration,
    pub chaining: Duration,
    pub extension_filter: Duration,
}

impl ExtensionTimer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_assign_from(&mut self, other: &mut ExtensionTimer) -> &mut Self {
        self.total_time += other.total_time;
        self.extension += other.extension;
        self.preprocessing_time += other.preprocessing_time;
        self.postprocessing_time += other.postprocessing_time;
        self.next += other.next;
        self.chaining += other.chaining;
        self.extension_filter += other.extension_filter;

        other.total_time = Duration::ZERO;
        other.extension = Duration::ZERO;
        other.preprocessing_time = Duration::ZERO;
        other.postprocessing_time = Duration::ZERO;
        other.next = Duration::ZERO;
        other.chaining = Duration::ZERO;
        other.extension_filter = Duration::ZERO;

        self
    }

    pub fn update(&mut self, operation: i32, duration: Duration) {
        match operation {
            0 => self.total_time += duration,
            1 => self.preprocessing_time += duration,
            2 => self.postprocessing_time += duration,
            3 => self.extension_filter += duration,
            4 => self.extension += duration,
            5 => self.next += duration,
            6 => self.chaining += duration,
            _ => {}
        }
    }
}

impl std::ops::AddAssign<&mut ExtensionTimer> for ExtensionTimer {
    fn add_assign(&mut self, other: &mut ExtensionTimer) {
        self.add_assign_from(other);
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TotalTime {
    pub timer: ExtensionTimer,
}

impl TotalTime {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn print(&self) {
        eprintln!(
            "Chaining: {} seconds",
            self.timer.chaining.as_nanos() as f64 / 1_000_000_000.0
        );
        eprintln!(
            "Extension-Time: {} seconds",
            self.timer.extension.as_nanos() as f64 / 1_000_000_000.0
        );
        eprintln!(
            "seed-lookup: {} seconds",
            self.timer.next.as_nanos() as f64 / 1_000_000_000.0
        );
        eprintln!(
            "build Hsp from cigar (part of extension): {} seconds",
            self.timer.postprocessing_time.as_nanos() as f64 / 1_000_000_000.0
        );
        eprintln!(
            "ungapped: {} seconds",
            self.timer.preprocessing_time.as_nanos() as f64 / 1_000_000_000.0
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extension_timer_update_and_add_assign_resets_other() {
        let ns = |n| Duration::from_nanos(n);
        let mut a = ExtensionTimer::new();
        a.update(0, ns(10));
        a.update(1, ns(20));
        a.update(2, ns(30));
        a.update(3, ns(40));
        a.update(4, ns(50));
        a.update(5, ns(60));
        a.update(6, ns(70));
        a.update(99, ns(999));

        assert_eq!(a.total_time, ns(10));
        assert_eq!(a.preprocessing_time, ns(20));
        assert_eq!(a.postprocessing_time, ns(30));
        assert_eq!(a.extension_filter, ns(40));
        assert_eq!(a.extension, ns(50));
        assert_eq!(a.next, ns(60));
        assert_eq!(a.chaining, ns(70));

        let mut b = ExtensionTimer::new();
        b.update(4, ns(7));
        a += &mut b;
        assert_eq!(a.extension, ns(57));
        assert_eq!(b.extension, Duration::ZERO);
    }
}
