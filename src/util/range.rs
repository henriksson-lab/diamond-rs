use std::ops::Sub;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Range<It> {
    begin: It,
    end: It,
}

impl<It> Range<It>
where
    It: Copy,
{
    pub fn new(begin: It, end: It) -> Self {
        Self { begin, end }
    }

    pub fn begin(&self) -> It {
        self.begin
    }

    pub fn end(&self) -> It {
        self.end
    }
}

impl<It> Range<It>
where
    It: Copy + Sub<Output = usize>,
{
    pub fn size(&self) -> usize {
        self.end - self.begin
    }
}

impl<It> Default for Range<It>
where
    It: Default,
{
    fn default() -> Self {
        Self {
            begin: It::default(),
            end: It::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_range_default_and_accessors() {
        let r = Range::<usize>::default();
        assert_eq!(r.begin(), 0);
        assert_eq!(r.end(), 0);
        assert_eq!(r.size(), 0);

        let r = Range::new(3usize, 9usize);
        assert_eq!(r.begin(), 3);
        assert_eq!(r.end(), 9);
        assert_eq!(r.size(), 6);
    }
}
