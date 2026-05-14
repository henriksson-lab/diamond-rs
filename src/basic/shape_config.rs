use super::consts::MAX_SHAPES;
use super::reduction::Reduction;
use super::shape::Shape;

#[derive(Clone, Debug)]
pub struct ShapeConfig {
    shapes: [Shape; MAX_SHAPES],
    n: usize,
}

impl Default for ShapeConfig {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeConfig {
    pub fn new() -> Self {
        Self {
            shapes: [Shape::default(); MAX_SHAPES],
            n: 0,
        }
    }

    pub fn from_codes(codes: &[String], count: u32, reduction: &Reduction) -> Result<Self, String> {
        let mut cfg = Self::new();
        let max_shapes = if count == 0 {
            codes.len()
        } else {
            (count as usize).min(codes.len())
        };

        for (i, code) in codes.iter().take(max_shapes).enumerate() {
            cfg.shapes[cfg.n] = Shape::from_code(code, reduction);
            if cfg.shapes[cfg.n].weight != cfg.shapes[0].weight {
                return Err("Seed shape weight has to be uniform.".to_string());
            }
            cfg.n += 1;
            assert_eq!(cfg.n, i + 1);
        }

        Ok(cfg)
    }

    pub fn count(&self) -> i32 {
        self.n as i32
    }

    pub fn get(&self, i: usize) -> &Shape {
        &self.shapes[i]
    }

    pub fn patterns(&self, begin: u32, end: u32) -> Vec<u32> {
        let mut v = Vec::new();
        for i in begin..end {
            v.push(self.shapes[i as usize].mask);
        }
        v
    }
}

impl std::ops::Index<usize> for ShapeConfig {
    type Output = Shape;

    fn index(&self, index: usize) -> &Self::Output {
        self.get(index)
    }
}

impl std::fmt::Display for ShapeConfig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.n {
            write!(f, "{}", self.shapes[i])?;
            if i < self.n - 1 {
                write!(f, ",")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shape_config_default() {
        let cfg = ShapeConfig::new();
        assert_eq!(cfg.count(), 0);
        assert_eq!(format!("{}", cfg), "");
    }

    #[test]
    fn test_shape_config_from_codes_count_display_and_patterns() {
        let reduction = Reduction::default_reduction();
        let codes = vec!["111".to_string(), "1011".to_string(), "1101".to_string()];
        let cfg = ShapeConfig::from_codes(&codes, 2, &reduction).unwrap();
        assert_eq!(cfg.count(), 2);
        assert_eq!(format!("{}", cfg), "111,1011");
        assert_eq!(cfg[0].weight, 3);
        assert_eq!(cfg.patterns(0, 2), vec![0b111, 0b1101]);
    }

    #[test]
    fn test_shape_config_rejects_nonuniform_weight() {
        let reduction = Reduction::default_reduction();
        let codes = vec!["111".to_string(), "101".to_string()];
        let err = ShapeConfig::from_codes(&codes, 0, &reduction).unwrap_err();
        assert_eq!(err, "Seed shape weight has to be uniform.");
    }
}
