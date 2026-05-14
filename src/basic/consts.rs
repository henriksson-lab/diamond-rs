pub const BUILD_VERSION: u32 = 178;
pub const MAX_SEED_WEIGHT: usize = 32;
pub const MAX_SHAPES: usize = 64;
pub const VERSION_STRING: &str = "2.1.24";
pub const PROGRAM_NAME: &str = "diamond";
pub const MAX_CONTEXT: i32 = 6;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_consts_match_cpp() {
        assert_eq!(BUILD_VERSION, 178);
        assert_eq!(MAX_SEED_WEIGHT, 32);
        assert_eq!(MAX_SHAPES, 64);
        assert_eq!(VERSION_STRING, "2.1.24");
        assert_eq!(PROGRAM_NAME, "diamond");
        assert_eq!(MAX_CONTEXT, 6);
    }
}
