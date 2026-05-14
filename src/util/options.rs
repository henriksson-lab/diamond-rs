#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OptionBase {
    pub id: String,
    pub desc: String,
    pub short_id: char,
    pub disabled: bool,
    pub group: Option<String>,
}

impl OptionBase {
    pub fn new(id: &str, short_id: char, desc: &str, disabled: bool, group: Option<&str>) -> Self {
        Self {
            id: id.to_string(),
            desc: desc.to_string(),
            short_id,
            disabled,
            group: group.map(|s| s.to_string()),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct OptionValue<T> {
    value: Option<T>,
    base: Option<OptionBase>,
}

impl<T> Default for OptionValue<T> {
    fn default() -> Self {
        Self {
            value: None,
            base: None,
        }
    }
}

impl<T> OptionValue<T> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_base(base: OptionBase) -> Self {
        Self {
            value: None,
            base: Some(base),
        }
    }

    pub fn present(&self) -> bool {
        self.value.is_some()
    }

    pub fn blank(&self) -> bool {
        !self.present()
    }

    pub fn set(&mut self, value: T) -> &mut Self {
        self.value = Some(value);
        self
    }

    pub fn set_if_blank(&mut self, value: T) -> &mut Self {
        if self.value.is_none() {
            self.value = Some(value);
        }
        self
    }

    pub fn require(&self) -> Result<(), String> {
        if self.present() {
            Ok(())
        } else if let Some(base) = &self.base {
            Err(format!(
                "Missing parameter: --{}/-{}",
                base.id, base.short_id
            ))
        } else {
            Err("Missing parameter".to_string())
        }
    }

    pub fn get(&self, default_value: T) -> T
    where
        T: Clone,
    {
        self.value.clone().unwrap_or(default_value)
    }

    pub fn get_present(&self) -> Result<T, String>
    where
        T: Clone,
    {
        self.value
            .clone()
            .ok_or_else(|| "Option::get_present".to_string())
    }

    pub fn unset(&mut self) {
        self.value = None;
    }

    pub fn as_ref(&self) -> Option<&T> {
        self.value.as_ref()
    }

    pub fn set_base_ptr(&mut self, base: OptionBase) {
        self.base = Some(base);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_option_value_basic() {
        let mut opt =
            OptionValue::with_base(OptionBase::new("threads", 'p', "threads", false, None));
        assert!(opt.blank());
        assert!(opt.require().unwrap_err().contains("--threads/-p"));
        assert_eq!(opt.get(4), 4);
        opt.set(8);
        assert!(opt.present());
        assert_eq!(opt.get(4), 8);
        assert_eq!(opt.get_present().unwrap(), 8);
        opt.set_if_blank(16);
        assert_eq!(opt.get_present().unwrap(), 8);
        opt.unset();
        opt.set_if_blank(16);
        assert_eq!(opt.get_present().unwrap(), 16);
    }

    #[test]
    fn test_option_value_string_and_base() {
        let mut opt = OptionValue::<String>::new();
        assert!(opt.get_present().is_err());
        opt.set_base_ptr(OptionBase::new("db", 'd', "database", false, Some("io")));
        opt.set("nr".to_string());
        assert_eq!(opt.as_ref().unwrap(), "nr");
        assert_eq!(opt.get("fallback".to_string()), "nr");
    }
}
