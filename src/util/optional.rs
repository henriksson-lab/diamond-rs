#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct InPlace;

pub const IN_PLACE: InPlace = InPlace;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NullOpt;

pub const NULLOPT: NullOpt = NullOpt;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BadOptionalAccess;

impl std::fmt::Display for BadOptionalAccess {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("bad optional access")
    }
}

impl std::error::Error for BadOptionalAccess {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Optional<T> {
    value: Option<T>,
}

impl<T> Default for Optional<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Optional<T> {
    pub fn new() -> Self {
        Self { value: None }
    }

    pub fn nullopt(_nullopt: NullOpt) -> Self {
        Self::new()
    }

    pub fn some(value: T) -> Self {
        Self { value: Some(value) }
    }

    pub fn in_place_with(_in_place: InPlace, value: T) -> Self {
        Self::some(value)
    }

    pub fn in_place(_in_place: InPlace) -> Self
    where
        T: Default,
    {
        Self::some(T::default())
    }

    pub fn assign_nullopt(&mut self, _nullopt: NullOpt) -> &mut Self {
        self.reset();
        self
    }

    pub fn assign(&mut self, value: T) -> &mut Self {
        self.value = Some(value);
        self
    }

    pub fn has_value(&self) -> bool {
        self.value.is_some()
    }

    pub fn value(&self) -> Result<&T, BadOptionalAccess> {
        self.value.as_ref().ok_or(BadOptionalAccess)
    }

    pub fn value_mut(&mut self) -> Result<&mut T, BadOptionalAccess> {
        self.value.as_mut().ok_or(BadOptionalAccess)
    }

    pub fn as_ref(&self) -> Option<&T> {
        self.value.as_ref()
    }

    pub fn as_mut(&mut self) -> Option<&mut T> {
        self.value.as_mut()
    }

    pub fn deref(&self) -> Result<&T, BadOptionalAccess> {
        self.value()
    }

    pub fn deref_mut(&mut self) -> Result<&mut T, BadOptionalAccess> {
        self.value_mut()
    }

    pub fn value_or(&self, default_value: T) -> T
    where
        T: Clone,
    {
        self.value.clone().unwrap_or(default_value)
    }

    pub fn reset(&mut self) {
        self.value = None;
    }

    pub fn take(&mut self) -> Self {
        Self {
            value: self.value.take(),
        }
    }

    pub fn emplace(&mut self, value: T) -> &mut T {
        self.value = Some(value);
        self.value.as_mut().unwrap()
    }

    pub fn swap(&mut self, other: &mut Self) {
        std::mem::swap(&mut self.value, &mut other.value);
    }
}

impl<T> PartialEq<NullOpt> for Optional<T> {
    fn eq(&self, _other: &NullOpt) -> bool {
        !self.has_value()
    }
}

impl<T> PartialEq<Optional<T>> for NullOpt {
    fn eq(&self, other: &Optional<T>) -> bool {
        !other.has_value()
    }
}

impl<T> From<T> for Optional<T> {
    fn from(value: T) -> Self {
        Self::some(value)
    }
}

pub fn swap<T>(a: &mut Optional<T>, b: &mut Optional<T>) {
    a.swap(b);
}

pub fn make_optional<T>(value: T) -> Optional<T> {
    Optional::some(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_optional_value_reset_and_access() {
        let mut x = Optional::<i32>::new();
        assert!(!x.has_value());
        assert_eq!(x.value().unwrap_err().to_string(), "bad optional access");
        assert_eq!(x.value_or(7), 7);

        x.assign(3);
        assert!(x.has_value());
        assert_eq!(*x.value().unwrap(), 3);
        *x.value_mut().unwrap() = 4;
        assert_eq!(*x.deref().unwrap(), 4);
        assert_eq!(x.as_ref(), Some(&4));
        *x.as_mut().unwrap() = 5;
        assert_eq!(*x.deref_mut().unwrap(), 5);
        let y = x.take();
        assert_eq!(y.value_or(0), 5);
        assert_eq!(x, NULLOPT);
        assert_eq!(NULLOPT, x);
        x.assign(4);
        x.reset();
        assert_eq!(x, Optional::nullopt(NULLOPT));
    }

    #[test]
    fn test_optional_emplace_swap_and_make() {
        let defaulted = Optional::<String>::in_place(IN_PLACE);
        assert_eq!(defaulted.value().unwrap(), "");

        let mut a = Optional::in_place_with(IN_PLACE, "a".to_string());
        let mut b = make_optional("b".to_string());
        assert_eq!(a.value().unwrap(), "a");
        assert_eq!(b.value().unwrap(), "b");

        swap(&mut a, &mut b);
        assert_eq!(a.value().unwrap(), "b");
        assert_eq!(b.value().unwrap(), "a");

        b.assign_nullopt(NULLOPT);
        assert!(!b.has_value());
        b.emplace("c".to_string());
        assert_eq!(b.value_or("x".to_string()), "c");
    }
}
