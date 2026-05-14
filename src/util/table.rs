#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Table {
    data: Vec<(String, String)>,
    max_len: usize,
}

impl Table {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn call_string(&mut self, s: &str, t: &str) -> &mut Self {
        self.add_string(s, t)
    }

    pub fn call_i64(&mut self, s: &str, n: i64, unit: &str) -> &mut Self {
        self.add_i64(s, n, unit)
    }

    pub fn call_u64(&mut self, s: &str, n: u64, unit: &str) -> &mut Self {
        self.add_u64(s, n, unit)
    }

    pub fn call_i32(&mut self, s: &str, n: i32, unit: &str) -> &mut Self {
        self.add_i32(s, n, unit)
    }

    pub fn call_u32(&mut self, s: &str, n: u32, unit: &str) -> &mut Self {
        self.add_u32(s, n, unit)
    }

    pub fn call_f64(&mut self, s: &str, n: f64, unit: &str) -> &mut Self {
        self.add_f64(s, n, unit)
    }

    pub fn add_string(&mut self, s: &str, t: &str) -> &mut Self {
        self.data.push((s.to_string(), t.to_string()));
        self.max_len = self.max_len.max(s.len());
        self
    }

    pub fn add_i64(&mut self, s: &str, n: i64, unit: &str) -> &mut Self {
        self.add_string(s, &format!("{}{}", n, unit))
    }

    pub fn add_u64(&mut self, s: &str, n: u64, unit: &str) -> &mut Self {
        self.add_string(s, &format!("{}{}", n, unit))
    }

    pub fn add_i32(&mut self, s: &str, n: i32, unit: &str) -> &mut Self {
        self.add_string(s, &format!("{}{}", n, unit))
    }

    pub fn add_u32(&mut self, s: &str, n: u32, unit: &str) -> &mut Self {
        self.add_string(s, &format!("{}{}", n, unit))
    }

    pub fn add_f64(&mut self, s: &str, n: f64, unit: &str) -> &mut Self {
        let value = if n >= 100.0 {
            format!("{}{}", n.round() as i64, unit)
        } else {
            format!("{}{}", format!("{:.6}", n), unit)
        };
        self.add_string(s, &value)
    }

    pub fn data(&self) -> &[(String, String)] {
        &self.data
    }
}

impl std::fmt::Display for Table {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (key, value) in &self.data {
            writeln!(f, "{:>width$}  {}", key, value, width = self.max_len)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_table_add_and_display() {
        let mut table = Table::new();
        table
            .add_string("a", "x")
            .add_i32("long", 7, "u")
            .add_f64("dbl", 123.4, "")
            .add_f64("small", 1.25, "x");
        table
            .call_string("s", "t")
            .call_i64("i64", -1, "")
            .call_u64("u64", 2, "")
            .call_i32("i32", -3, "")
            .call_u32("u32", 4, "")
            .call_f64("f64", 5.5, "");
        assert_eq!(table.data()[0], ("a".to_string(), "x".to_string()));
        assert!(format!("{}", table).contains(" f64  5.500000\n"));
    }
}
