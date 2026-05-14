use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FieldValue<T> {
    pub v: T,
    pub primary: bool,
}

impl<T> FieldValue<T> {
    pub fn new(v: T, primary: bool) -> Self {
        Self { v, primary }
    }

    pub fn primary(v: T) -> Self {
        Self { v, primary: true }
    }
}

impl<T> From<T> for FieldValue<T> {
    fn from(v: T) -> Self {
        Self::primary(v)
    }
}

pub type EMap<T> = BTreeMap<T, String>;
pub type SEMap<T> = BTreeMap<String, FieldValue<T>>;

pub fn enum_to_string<T: Copy + Eq>(
    v: T,
    map: &[(T, &'static str)],
) -> Result<&'static str, String> {
    for &(key, value) in map {
        if key == v {
            return Ok(value);
        }
    }
    Err("Invalid conversion from enum to string.".to_string())
}

pub fn enum_from_string<T: Copy>(
    s: &str,
    map: &[(&'static str, FieldValue<T>)],
) -> Result<T, String> {
    for &(key, value) in map {
        if key == s {
            return Ok(value.v);
        }
    }

    let mut permitted = String::new();
    let mut n = 0usize;
    for &(key, value) in map {
        if !value.primary {
            continue;
        }
        if n > 0 {
            permitted.push_str(", ");
        }
        permitted.push_str(key);
        n += 1;
    }
    Err(format!(
        "Invalid value for string field: {}. Permitted values: {}",
        s, permitted
    ))
}

pub fn to_string<T>(v: T, map: &EMap<T>) -> Result<String, String>
where
    T: Ord,
{
    map.get(&v)
        .cloned()
        .ok_or_else(|| "Invalid conversion from enum to string.".to_string())
}

pub fn from_string<T>(s: &str, map: &SEMap<T>) -> Result<T, String>
where
    T: Copy,
{
    if let Some(value) = map.get(s) {
        return Ok(value.v);
    }

    let mut permitted = String::new();
    let mut n = 0usize;
    for (key, value) in map {
        if !value.primary {
            continue;
        }
        if n > 0 {
            permitted.push_str(", ");
        }
        permitted.push_str(key);
        n += 1;
    }
    Err(format!(
        "Invalid value for string field: {}. Permitted values: {}",
        s, permitted
    ))
}

pub trait FlagBits {
    type Bits: Copy
        + From<u8>
        + PartialEq
        + std::ops::BitAnd<Output = Self::Bits>
        + std::ops::Not<Output = Self::Bits>;

    fn bits(self) -> Self::Bits;
}

pub fn flag_all<T: FlagBits>(a: T, b: T) -> bool {
    let c = b.bits();
    (a.bits() & c) == c
}

pub fn flag_any<T: FlagBits>(a: T, b: T) -> bool {
    (a.bits() & b.bits()) != T::Bits::from(0)
}

pub fn flag_only<T: FlagBits>(a: T, b: T) -> bool {
    (a.bits() & !b.bits()) == T::Bits::from(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
    enum Mode {
        A,
        B,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    struct Flags(u32);

    impl FlagBits for Flags {
        type Bits = u32;

        fn bits(self) -> Self::Bits {
            self.0
        }
    }

    #[test]
    fn test_enum_string_helpers() {
        let to = [(Mode::A, "a"), (Mode::B, "b")];
        let from = [
            ("a", FieldValue::primary(Mode::A)),
            ("alias-a", FieldValue::new(Mode::A, false)),
            ("b", FieldValue::primary(Mode::B)),
        ];
        assert_eq!(enum_to_string(Mode::B, &to).unwrap(), "b");
        assert_eq!(enum_from_string("alias-a", &from).unwrap(), Mode::A);
        assert_eq!(
            enum_from_string("x", &from).unwrap_err(),
            "Invalid value for string field: x. Permitted values: a, b"
        );

        let mut to_map = EMap::new();
        to_map.insert(Mode::A, "a".to_string());
        to_map.insert(Mode::B, "b".to_string());
        let mut from_map = SEMap::new();
        from_map.insert("a".to_string(), FieldValue::from(Mode::A));
        from_map.insert("alias-a".to_string(), FieldValue::new(Mode::A, false));
        from_map.insert("b".to_string(), FieldValue::from(Mode::B));
        assert_eq!(to_string(Mode::A, &to_map).unwrap(), "a");
        assert_eq!(from_string("b", &from_map).unwrap(), Mode::B);
        assert_eq!(
            from_string("x", &from_map).unwrap_err(),
            "Invalid value for string field: x. Permitted values: a, b"
        );
    }

    #[test]
    fn test_flag_helpers() {
        assert!(flag_all(Flags(0b111), Flags(0b011)));
        assert!(!flag_all(Flags(0b101), Flags(0b011)));
        assert!(flag_any(Flags(0b100), Flags(0b110)));
        assert!(!flag_any(Flags(0b001), Flags(0b110)));
        assert!(flag_only(Flags(0b010), Flags(0b110)));
        assert!(!flag_only(Flags(0b001), Flags(0b110)));
    }
}
